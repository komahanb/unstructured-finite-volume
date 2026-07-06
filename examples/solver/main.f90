!=====================================================================!
! Generic finite volume driver - the problem lives entirely in a config
! file, nothing is hardcoded here.
!
!   ./run <config-file>
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program solver

  use iso_fortran_env          , only : dp => real64
  use class_config             , only : config
  use class_gmsh_loader        , only : gmsh_loader
  use class_mesh               , only : mesh
  use class_assembler          , only : assembler, CONVECTION_UPWIND
  use class_diffusion_flux  , only : diffusion_flux, constant_source
  use class_advection_flux     , only : advection_diffusion_flux
  use interface_linear_solver  , only : linear_solver
  use class_conjugate_gradient , only : conjugate_gradient
  use class_gmres              , only : gmres_solver
  use class_normal_cg          , only : normal_cg, CGNR_METHOD, CGNE_METHOD
  use class_partitioned_assembler, only : partitioned_assembler
  use class_sor                , only : sor
  use class_gauss_seidel       , only : gauss_seidel
  use class_gauss_jacobi       , only : gauss_jacobi
  use class_csr                , only : csr_matrix
  use class_algebraic_multigrid, only : algebraic_multigrid
  use class_bdf                , only : bdf
  use class_paraview_writer    , only : paraview_writer
  use class_gmsh_writer        , only : gmsh_writer
  use class_string             , only : string
  use module_verbosity         , only : set_verbosity

  implicit none

  type(config)                       :: cfg
  character(len=256)                 :: cfgfile
  class(gmsh_loader)  , allocatable  :: gmsh
  class(mesh)         , allocatable  :: grid
  class(assembler)    , allocatable  :: fvm
  class(linear_solver), allocatable  :: lsolver
  type(csr_matrix)                   :: amg_op       ! assembled operator (for pcg-amg)
  type(algebraic_multigrid)                          :: amg_precond  ! SA-AMG preconditioner
  type(bdf)                          :: ti
  class(paraview_writer), allocatable :: pw
  real(dp), allocatable              :: x(:)
  type(string)                       :: labels(1)
  integer                            :: i

  ! Read the run
  call get_command_argument(1, cfgfile)
  if (len_trim(cfgfile) .eq. 0) then
     print *, "usage: run <config-file>"
     error stop
  end if
  cfg = config(trim(cfgfile))
  call cfg % print()
  call set_verbosity(cfg % verbosity)

  ! Geometry
  allocate(gmsh, source = gmsh_loader(trim(cfg % meshfile % str)))
  allocate(grid, source = mesh(gmsh))

  ! Assembler + physics (set the equation before the boundaries). The
  ! distributed solve uses the same plain cg on a PARTITIONED system.
  if (trim(cfg % solver % str) .eq. "distributed_cg") then
     allocate(fvm, source = partitioned_assembler(grid))
  else
     allocate(fvm, source = assembler(grid))
  end if
  select case (trim(cfg % equation % str))
  case ("diffusion")
     call fvm % set_equation(diffusion_flux(cfg % kappa), constant_source(cfg % source))
     ! the diffusion operator is symmetric: declare it on the configured
     ! instance so transpose products (cgnr/cgne, adjoints) are an
     ! explicit, verifiable claim
     fvm % operator_is_symmetric = .true.
  case ("advection_diffusion")
     call fvm % set_equation(advection_diffusion_flux(cfg % velocity, cfg % kappa), &
          &                  constant_source(cfg % source))
     if (trim(cfg % convection % str) .eq. "upwind") &
          & call fvm % set_convection_scheme(CONVECTION_UPWIND)
  case default
     print *, "unknown equation '", trim(cfg % equation % str), "'"
     error stop
  end select

  ! Boundary conditions, by physical group name
  do i = 1, cfg % nbc
     select case (trim(cfg % bc_kind(i) % str))
     case ("dirichlet")
        call fvm % set_dirichlet(trim(cfg % bc_name(i) % str), cfg % bc_c(i))
     case ("neumann")
        call fvm % set_neumann(trim(cfg % bc_name(i) % str), cfg % bc_c(i))
     case ("robin")
        call fvm % set_robin(trim(cfg % bc_name(i) % str), cfg % bc_a(i), cfg % bc_b(i), cfg % bc_c(i))
     end select
  end do

  ! partition the system across the images (after boundary conditions)
  select type (fvm)
  type is (partitioned_assembler)
     call fvm % setup_partition()
  end select

  if (cfg % dt .gt. 0.0_dp) then

     ! Transient: backward-euler (bdf order 1) march from a zero initial state
     ti = bdf(fvm, cfg % tinit, cfg % tfinal, cfg % dt, max_order = 1)
     call ti % integrate()
     x = real(ti % U(ti % num_steps, :, 1), dp)

  else

     ! Steady: pick the linear solver
     select case (trim(cfg % solver % str))
     case ("cg")
        allocate(lsolver, source = conjugate_gradient(&
             & max_tol=cfg % max_tol, max_it=cfg % max_it, print_level=1))
     case ("sor")
        allocate(lsolver, source = sor(omega=cfg % omega, &
             & max_tol=cfg % max_tol, max_it=cfg % max_it, print_level=1))
     case ("gs")
        allocate(lsolver, source = gauss_seidel(&
             & max_tol=cfg % max_tol, max_it=cfg % max_it, print_level=1))
     case ("gj")
        allocate(lsolver, source = gauss_jacobi(&
             & max_tol=cfg % max_tol, max_it=cfg % max_it, print_level=1))
     case ("pcg-amg")
        ! smoothed-aggregation AMG-preconditioned CG: assemble the operator
        ! once, build the multigrid hierarchy, hand it to CG as the precond
        call fvm % get_operator_csr(amg_op)
        call amg_precond % setup(amg_op)
        allocate(lsolver, source = conjugate_gradient(&
             & max_tol=cfg % max_tol, max_it=cfg % max_it, print_level=1, &
             & precond=amg_precond))
     case ("gmres")
        ! restarted GMRES on the assembled operator (for nonsymmetric /
        ! advection problems where CG does not apply)
        allocate(lsolver, source = gmres_solver(&
             & max_tol=cfg % max_tol, max_it=cfg % max_it, print_level=1))
     case ("cgnr")
        ! CG on the normal equations A^T A (robust baseline; kappa^2 convergence)
        allocate(lsolver, source = normal_cg(&
             & max_tol=cfg % max_tol, max_it=cfg % max_it, method=CGNR_METHOD, print_level=1))
     case ("cgne")
        ! CG on the normal equations A A^T (Craig)
        allocate(lsolver, source = normal_cg(&
             & max_tol=cfg % max_tol, max_it=cfg % max_it, method=CGNE_METHOD, print_level=1))
     case ("distributed_cg")
        ! the system is partitioned across the images (above); the solver
        ! is the ordinary cg - parallelism lives on the system side.
        ! (serial build / one image reduces to plain CG); SPD only
        allocate(lsolver, source = conjugate_gradient(&
             & max_tol=cfg % max_tol, max_it=cfg % max_it, print_level=1))
     case default
        print *, "unknown solver '", trim(cfg % solver % str), "'"
        error stop
     end select

     call lsolver % solve(fvm, x)

  end if

  print *, "solution: min/max/mean =", minval(x), maxval(x), sum(x)/size(x)

  ! Output - paraview (.vtu) and gmsh (.msh) views of the same field
  labels(1) = string("phi")
  allocate(pw, source = paraview_writer(fvm % grid))
  call pw % write(trim(cfg % output % str), reshape(x, [size(x), 1]), labels)

  ! gmsh post-processing file (input mesh + cell field); open: gmsh <out>.msh
  gmsh_output: block

    type(gmsh_writer)             :: gw
    character(len=:), allocatable :: gout
    integer                       :: idot

    if (size(x) .eq. fvm % grid % num_cells) then

       gout = trim(cfg % output % str)
       idot = index(gout, ".vtu", back = .true.)
       if (idot .gt. 0) gout = gout(1:idot-1)//".msh"

       gw = gmsh_writer(trim(cfg % meshfile % str))
       call gw % write(gout, fvm % grid % cell_numbers, x, "phi")

    end if

  end block gmsh_output

end program solver
