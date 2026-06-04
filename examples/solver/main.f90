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
  use class_assembler          , only : assembler
  use class_diffusion_flux  , only : diffusion_flux, constant_source
  use interface_linear_solver  , only : linear_solver
  use class_conjugate_gradient , only : conjugate_gradient
  use class_sor                , only : sor
  use class_gauss_seidel       , only : gauss_seidel
  use class_gauss_jacobi       , only : gauss_jacobi
  use class_time_integrator    , only : time_integrator
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
  type(time_integrator)              :: ti
  class(paraview_writer), allocatable :: pw
  real(dp), allocatable              :: x(:), phi0(:)
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

  ! Assembler + physics (set the equation before the boundaries)
  allocate(fvm, source = assembler(grid))
  select case (trim(cfg % equation % str))
  case ("diffusion")
     call fvm % set_equation(diffusion_flux(cfg % kappa), constant_source(cfg % source))
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

  if (cfg % dt .gt. 0.0_dp) then

     ! Transient: backward-euler march from a zero initial state
     ti = time_integrator(fvm, cfg % tinit, cfg % tfinal, cfg % dt, cfg % max_it, cfg % max_tol)
     allocate(phi0(fvm % num_state_vars)); phi0 = 0.0_dp
     call ti % integrate(phi0, x)

  else

     ! Steady: pick the linear solver
     select case (trim(cfg % solver % str))
     case ("cg")
        allocate(lsolver, source = conjugate_gradient(FVAssembler=fvm, &
             & max_tol=cfg % max_tol, max_it=cfg % max_it, print_level=1))
     case ("sor")
        allocate(lsolver, source = sor(FVAssembler=fvm, omega=cfg % omega, &
             & max_tol=cfg % max_tol, max_it=cfg % max_it, print_level=1))
     case ("gs")
        allocate(lsolver, source = gauss_seidel(FVAssembler=fvm, &
             & max_tol=cfg % max_tol, max_it=cfg % max_it, print_level=1))
     case ("gj")
        allocate(lsolver, source = gauss_jacobi(FVAssembler=fvm, &
             & max_tol=cfg % max_tol, max_it=cfg % max_it, print_level=1))
     case default
        print *, "unknown solver '", trim(cfg % solver % str), "'"
        error stop
     end select

     call lsolver % solve(x)

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
