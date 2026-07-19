!=====================================================================!
! Test loading of mesh and mesh pre-processing
!=====================================================================!

program test_mesh

  use iso_fortran_env          , only : dp => REAL64
  use class_gmsh_loader        , only : gmsh_loader
  use class_mesh               , only : mesh
  use class_assembler          , only : assembler
  use interface_linear_solver  , only : linear_solver, AUTO
  use class_conjugate_gradient , only : conjugate_gradient
  use class_gauss_jacobi       , only : gauss_jacobi
  use class_gauss_seidel       , only : gauss_seidel
  use class_sor                , only : sor
  use class_paraview_writer    , only : paraview_writer
  use class_string             , only : string

  implicit none

  character(len=*)      , parameter   :: filename = "../box-3.msh"
  class(gmsh_loader)    , allocatable :: gmsh
  class(mesh)           , allocatable :: grid
  class(linear_solver)  , allocatable :: solver
  class(assembler)      , allocatable :: FVMAssembler
  class(paraview_writer), allocatable :: paraview

  meshing : block

    ! Create a mesh object
    allocate(gmsh, source = gmsh_loader(filename))
    allocate(grid, source = mesh(gmsh))
    deallocate(gmsh)

  end block meshing

  assembly : block

    ! Create an assembler object for assembling the linear system
    ! Geometry and meshing
    allocate(FVMAssembler, source = assembler(grid))

    ! Boundary conditions, by physical group name (box-3.msh)
    call FVMAssembler % set_dirichlet("front" , 5.0d0)
    call FVMAssembler % set_dirichlet("bottom", 10.0d0)
    call FVMAssembler % set_dirichlet("right" , 15.0d0)
    call FVMAssembler % set_dirichlet("top"   , 0.0d0)
    call FVMAssembler % set_dirichlet("left"  , 0.0d0)
    call FVMAssembler % set_dirichlet("back"  , 0.0d0)

    ! Also supply
    ! allocate(FVMAssembler, source = assembler(grid,physics_list))
    ! physics with tags Assembler combines Geometry and Physics ( EQNS
    ! + BC) to provide linear/nonlinear systems

  end block assembly

  cg_solver : block

    real(dp) , parameter   :: max_tol     = 100.0d0*epsilon(1.0d0)
    integer  , parameter   :: max_it      = 100
    integer  , parameter   :: print_level = 1
    real(dp) , allocatable :: x(:)
    real(dp) , allocatable :: phi(:,:)
    type(string), allocatable :: field_names(:)
    integer :: i

    allocate(solver, &
         & source      = conjugate_gradient( &
         & max_tol     = max_tol, &
         & max_it      = max_it, &
         & print_level = print_level))

    ! Solve using solver method
    call solver % solve(FVMassembler, x)
    print *, 'cg solution = '
    do i = 1, min(10, size(x))
       print *, i,  x(i)
    end do

    ! Writes the mesh for tecplot
    ! call FVMassembler % write_solution("mesh-cg.dat", x)

    allocate(field_names(1))
    field_names(1) = string("phi")
    allocate(paraview, source = paraview_writer(FVMAssembler % grid))
    allocate(phi(size(x),1))
    phi(:,1) = x

    call paraview % write("output.vtu", phi, field_names)

    deallocate(phi)
    deallocate(paraview)
    deallocate(x)
    deallocate(solver)

  end block cg_solver

  ! auto-tuned sor: the solver works its own parameter out - measures the
  ! sweep's convergence factor at entry and sets the optimal omega, with
  ! the rollback gate live during the outer iteration. must agree with cg.
  auto_sor : block

    real(dp) , parameter   :: max_tol = 1.0d-10
    integer  , parameter   :: max_it  = 500
    real(dp) , allocatable :: x(:), xref(:)
    type(sor), allocatable :: tuned
    real(dp)               :: omega0, diff

    allocate(tuned, source = sor(omega = 1.0d0, max_it = max_it, &
         & max_tol = max_tol, print_level = 1))
    tuned % tuning = AUTO
    omega0 = tuned % omega

    call tuned % solve(FVMassembler, x)

    allocate(solver, source = conjugate_gradient(max_tol = max_tol, &
         & max_it = 100, print_level = 0))
    call solver % solve(FVMassembler, xref)

    ! the iteration diagnostic lives on the object (no module-scope
    ! side channels): a completed march must have populated it
    if (solver % last_inner_iters .le. 0) then
       write(*,'(1x,a)') "FAIL : last_inner_iters not populated by the march"
       error stop
    else
       write(*,'(1x,a,i0,a)') "PASS : last_inner_iters populated (", &
            & solver % last_inner_iters, " inner iterations)"
    end if

    diff = norm2(x - xref)/norm2(xref)
    if (abs(tuned % omega - omega0) .gt. 0.0d0 .and. diff .lt. 1.0d-8) then
       write(*,'(1x,a,f8.5,a,es10.3)') &
            & "PASS : auto-tuned sor (omega ", tuned % omega, ") matches cg, diff ", diff
    else
       write(*,'(1x,a,f8.5,a,es10.3)') &
            & "FAIL : auto-tuned sor (omega ", tuned % omega, "), diff vs cg ", diff
       error stop
    end if

    deallocate(x, xref, tuned, solver)

  end block auto_sor

  deallocate(grid)
  deallocate(FVMAssembler)

  contains

end program test_mesh
