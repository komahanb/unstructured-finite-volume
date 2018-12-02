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

  implicit none
  
  character(len=*)     , parameter   :: filename = "rectangle.msh"
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

    real(dp), allocatable :: A(:,:)
    real(dp), allocatable :: AT(:,:)

    ! Create an assembler object for assembling the linear system
    ! Geometry and meshing
    allocate(FVMAssembler, source = assembler(grid))

    call FVMAssembler % get_jacobian(A, filter=FVMAssembler % UPPER_TRIANGLE)
    allocate(AT, mold = A); AT = transpose(A)   
    print *, "asymmetry", maxval(abs(A-AT))

    call print (A)

    stop

    ! Also supply
    ! allocate(FVMAssembler, source = assembler(grid,physics_list)) 
    ! physics with tags Assembler combines Geometry and Physics ( EQNS
    ! + BC) to provide linear/nonlinear systems

  end block assembly

!!$  solver : block
!!$
!!$    real(dp) , parameter   :: max_tol     = 100.0d0*epsilon(1.0d0)
!!$    integer  , parameter   :: max_it      = 10000
!!$    integer  , parameter   :: print_level = 1
!!$    real(dp) , allocatable :: x(:)
!!$    integer :: i
!!$
!!$    allocate(CG, &
!!$         & source      = conjugate_gradient( &
!!$         & FVAssembler = FVMassembler, &
!!$         & max_tol     = max_tol, &
!!$         & max_it      = max_it, &
!!$         & print_level = print_level))
!!$
!!$    ! Solve using CG method
!!$    call CG % solve(x)
!!$    print *, 'solution = '
!!$    do i = 1, min(10, size(x))
!!$       print *, i,  x(i)
!!$    end do
!!$
!!$    ! Writes the mesh for tecplot
!!$    call FVMassembler % write_solution("mesh.dat", x)
!!$
!!$    deallocate(x)   
!!$    deallocate(CG)
!!$
!!$  end block solver

  deallocate(grid)
  deallocate(FVMAssembler)

  contains
  
  subroutine print(matrix)

    real(dp), intent(in) :: matrix(:,:)
    integer :: i, j, m, n

    m = size(matrix, 1)
    n = size(matrix, 2)

    do i = 1, m
       write(*,'(100g15.1)') (matrix(i,j),j=1,n)
    end do

!!$    do i = 1, m
!!$       print *, sum(matrix(:,i))
!!$    end do

  end subroutine print

end program test_mesh
