!=====================================================================!
! Test loading of mesh and mesh pre-processing
!=====================================================================!

program test_mesh

  use iso_fortran_env          , only : dp => REAL64
  use class_gmsh_loader        , only : gmsh_loader
  use class_mesh               , only : mesh
  use class_assembler          , only : assembler, DIAGONAL, LOWER_TRIANGLE, UPPER_TRIANGLE

  implicit none

  character(len=*)     , parameter   :: filename = "../box-36.msh"
  class(gmsh_loader)   , allocatable :: gmsh
  class(mesh)          , allocatable :: grid
  class(assembler)     , allocatable :: FVMAssembler

  meshing : block

    ! Create a mesh object
    allocate(gmsh, source =  gmsh_loader(filename))
    allocate(grid, source = mesh(gmsh))
    deallocate(gmsh)

  end block meshing

  assembly : block

    real(dp), allocatable :: A(:,:)
    real(dp), allocatable :: D(:,:)
    real(dp), allocatable :: U(:,:)
    real(dp), allocatable :: L(:,:)

    real(dp), allocatable :: AT(:,:)

    ! Create an assembler object for assembling the linear system
    ! Geometry and meshing
    allocate(FVMAssembler, source = assembler(grid))

    !-----------------------------------------------------------------!
    ! Assemble parts of the jacobian and test
    !-----------------------------------------------------------------!

    print *, 'getting upper triangle'
    call FVMAssembler % get_jacobian(U, filter = UPPER_TRIANGLE)
    call print (U)

    print *, 'getting lower triangle'
    call FVMAssembler % get_jacobian(L, filter = LOWER_TRIANGLE)
    call print (L)

    print *, 'getting diagonal matrix'
    call FVMAssembler % get_jacobian(D, filter = DIAGONAL)
    call print (D)

    print *, 'getting full jacobian'
    call FVMAssembler % get_jacobian(A)
    call print (A)

    ! Check consistency of matrix assembly
    if (maxval(abs(A-L-U-D)) .gt. tiny(1.0d0)) then
       error stop "error in assembly"
    else
       print *, 'passed assembly test'
    end if

    print *, 'performing symmetry test'
    call FVMAssembler % get_transpose_jacobian(AT)
    call print (AT)
    print *, "asymmetry", (maxval(abs(A-AT)) .gt. tiny(1.0d0))

    deallocate(A, D, U, L, AT)

  end block assembly

  ! analytic consistency checks of the jacobian-vector product, exact
  ! to machine precision
  consistency_checks: block

    real(dp) :: defect

    defect = FVMAssembler % verify_parts_consistency()
    if (defect .lt. 1.0d-13) then
       write(*,'(1x,a,es10.3)') "PASS : parts consistency (D+L+U)v = Av, defect ", defect
    else
       write(*,'(1x,a,es10.3)') "FAIL : parts consistency, defect ", defect
       error stop
    end if

    defect = FVMAssembler % verify_transpose_consistency()
    if (defect .lt. 1.0d-13) then
       write(*,'(1x,a,es10.3)') "PASS : transpose consistency <w,Av> = <A^T w,v>, defect ", defect
    else
       write(*,'(1x,a,es10.3)') "FAIL : transpose consistency, defect ", defect
       error stop
    end if

  end block consistency_checks

  deallocate(grid)
  deallocate(FVMAssembler)

  contains

  subroutine print(matrix)

    real(dp), intent(in) :: matrix(:,:)
    integer :: i, j, m, n

    m = size(matrix, 1)
    n = size(matrix, 2)

    do i = 1, min(10,m)
       write(*,'(100g15.1)') (matrix(i,j),j=1,min(10,n))
    end do

!!$    do i = 1, m
!!$       print *, sum(matrix(:,i))
!!$    end do

  end subroutine print

end program test_mesh
