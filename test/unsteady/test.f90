!=====================================================================!
! Drive the assembled system through the linear solver and write the
! solution - the unsteady-solve driver.
!=====================================================================!

program test_mesh

  use iso_fortran_env            , only : dp => REAL64
  use class_gmsh_loader          , only : gmsh_loader
  use class_mesh                 , only : mesh
  use class_assembler            , only : assembler
  use interface_linear_solver    , only : linear_solver
  !use interface_nonlinear_solver , only : nonlinear_solver
  use class_conjugate_gradient   , only : conjugate_gradient

  implicit none

  character(len=*)        , parameter   :: filename = "../rectangle.msh"
  !class(gmsh_loader)      , allocatable :: gmsh
  !class(physics)          , allocatable :: heat
  class(mesh)             , allocatable :: grid
  class(assembler)        , allocatable :: FVMAssembler
  class(linear_solver)    , allocatable :: linear
  !class(nonlinear_solver) , allocatable :: nonlinear

  real(dp) , parameter   :: max_tol     = 100.0d0*epsilon(1.0d0)
  integer  , parameter   :: max_it      = 100
  integer  , parameter   :: print_level = 1
  real(dp) , allocatable :: x(:)
  integer :: i

  ! Create a mesh object.
  allocate(grid, source = mesh(gmsh_loader(filename)))

  ! Create the assembler object coordinating the geometry and the physics.
  allocate(FVMAssembler, source = assembler(grid))!, physics))

  ! Set the boundary conditions by physical group name (rectangle.msh, 2D).
  call FVMAssembler % set_dirichlet("BoundaryLeft"  , 1.0d0)
  call FVMAssembler % set_dirichlet("BoundaryRight" , 0.0d0)
  call FVMAssembler % set_neumann  ("BoundaryTop"   , 0.0d0)
  call FVMAssembler % set_neumann  ("BoundaryBottom", 0.0d0)

  ! Set up the linear solver.
  allocate(linear, &
       & source      = conjugate_gradient( &
       & max_tol     = max_tol, &
       & max_it      = max_it, &
       & print_level = print_level))

  ! Solve using the solver method.
  call linear % solve(FVMassembler, x)
  print *, 'cg solution = '
  do i = 1, min(10, size(x))
     print *, i,  x(i)
  end do

  ! Write the solution in Tecplot format.
  call FVMassembler % write_solution("mesh-cg.dat", x)

  deallocate(x)
  !deallocate(linear,nonlinear, FVMAssembler, physics, grid)

end program test_mesh
