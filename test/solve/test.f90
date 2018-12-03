!=====================================================================!
! Test loading of mesh and mesh pre-processing
!=====================================================================!

program test_mesh

  use iso_fortran_env          , only : dp => REAL64
  use class_gmsh_loader        , only : gmsh_loader
  use class_mesh               , only : mesh
  use class_assembler          , only : assembler
  use interface_linear_solver  , only : linear_solver
  use class_conjugate_gradient , only : conjugate_gradient
  use class_gauss_jacobi       , only : gauss_jacobi
  use class_gauss_seidel       , only : gauss_seidel

  implicit none
  
  character(len=*)     , parameter   :: filename = "../square-10.msh"
  class(gmsh_loader)   , allocatable :: gmsh
  class(mesh)          , allocatable :: grid
  class(linear_solver) , allocatable :: CG
  class(assembler)     , allocatable :: FVMAssembler

  meshing : block

    ! Create a mesh object
    allocate(gmsh, source =  gmsh_loader(filename))
    allocate(grid, source = mesh(gmsh))
    deallocate(gmsh)
    
  end block meshing

  assembly : block

    ! Create an assembler object for assembling the linear system
    ! Geometry and meshing
    allocate(FVMAssembler, source = assembler(grid))

    ! Also supply
    ! allocate(FVMAssembler, source = assembler(grid,physics_list)) 
    ! physics with tags Assembler combines Geometry and Physics ( EQNS
    ! + BC) to provide linear/nonlinear systems

  end block assembly
  
  seidel_solver : block

    real(dp) , parameter   :: max_tol     = 100.0d0*epsilon(1.0d0)
    integer  , parameter   :: max_it      = 100
    integer  , parameter   :: print_level = 1
    real(dp) , allocatable :: x(:)
    integer :: i

    allocate(CG, &
         & source      = gauss_seidel( &
         & FVAssembler = FVMassembler, &
         & max_tol     = max_tol, &
         & max_it      = max_it, &
         & print_level = print_level))

    ! Solve using Seidel method
    call CG % solve(x)
    print *, 'seidel solution = '
    do i = 1, min(10, size(x))
       print *, i,  x(i)
    end do

    ! Writes the mesh for tecplot
    call FVMassembler % write_solution("mesh-seidel.dat", x)

    deallocate(x)   
    deallocate(CG)
    
  end block seidel_solver

  jacobi_solver : block

    real(dp) , parameter   :: max_tol     = 100.0d0*epsilon(1.0d0)
    integer  , parameter   :: max_it      = 100
    integer  , parameter   :: print_level = 1
    real(dp) , allocatable :: x(:)
    integer :: i

    allocate(CG, &
         & source      = gauss_jacobi( &
         & FVAssembler = FVMassembler, &
         & max_tol     = max_tol, &
         & max_it      = max_it, &
         & print_level = print_level))

    ! Solve using Jacobi method
    call CG % solve(x)
    print *, 'jacobi solution = '
    do i = 1, min(10, size(x))
       print *, i,  x(i)
    end do

    ! Writes the mesh for tecplot
    call FVMassembler % write_solution("mesh-jacobi.dat", x)

    deallocate(x)   
    deallocate(CG)
    
  end block jacobi_solver

  cg_solver : block

    real(dp) , parameter   :: max_tol     = 100.0d0*epsilon(1.0d0)
    integer  , parameter   :: max_it      = 100
    integer  , parameter   :: print_level = 1
    real(dp) , allocatable :: x(:)
    integer :: i

    allocate(CG, &
         & source      = conjugate_gradient( &
         & FVAssembler = FVMassembler, &
         & max_tol     = max_tol, &
         & max_it      = max_it, &
         & print_level = print_level))

    ! Solve using CG method
    call CG % solve(x)
    print *, 'cg solution = '
    do i = 1, min(10, size(x))
       print *, i,  x(i)
    end do

    ! Writes the mesh for tecplot
    call FVMassembler % write_solution("mesh-cg.dat", x)

    deallocate(x)   
    deallocate(CG)
    
  end block cg_solver

  deallocate(grid)
  deallocate(FVMAssembler)

  contains

end program test_mesh
