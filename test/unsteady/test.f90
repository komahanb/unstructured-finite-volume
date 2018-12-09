!=====================================================================!
! Test loading of mesh and mesh pre-processing
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
  
  ! Create a mesh object
  allocate(grid, source = mesh(gmsh_loader(filename)))
  
  ! Assembler Object coordinating geometry and physics
  allocate(FVMAssembler, source = assembler(grid))!, physics))

  ! Nonlinear Solution
!!$  allocate(nonlinear, &
!!$       & source = newton( &
!!$       & FVMAssembler, max_tol, max_it, print_level) &
!!$       & )
  
  ! Linear Solution
  allocate(linear, &
       & source      = conjugate_gradient( &
       & FVAssembler = FVMassembler, &
       & max_tol     = max_tol, &
       & max_it      = max_it, &
       & print_level = print_level))
  
  ! Solve using solver method
  call linear % solve(x)
  print *, 'cg solution = '
  do i = 1, min(10, size(x))
     print *, i,  x(i)
  end do
  
  ! Writes the mesh for tecplot
  call FVMassembler % write_solution("mesh-cg.dat", x)
  
  deallocate(x)
  !deallocate(linear,nonlinear, FVMAssembler, physics, grid)
  
  contains

end program test_mesh
