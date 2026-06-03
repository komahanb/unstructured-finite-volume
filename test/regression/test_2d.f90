!=====================================================================!
! Two-dimensional verification on the unit-square mesh square-10. Two
! analytic oracles:
!
!   1. constant dirichlet on all sides  -> u = 5 everywhere
!   2. left=0, right=1, top/bottom insulated (neumann 0) -> u = x,
!      which the finite volume laplacian reproduces exactly on the
!      orthogonal quad mesh
!
! A nonzero exit (error stop) means a 2d check failed.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_2d

  use iso_fortran_env          , only : dp => REAL64
  use class_gmsh_loader        , only : gmsh_loader
  use class_mesh               , only : mesh
  use class_assembler          , only : assembler
  use class_conjugate_gradient , only : conjugate_gradient

  implicit none

  integer  :: nfail
  real(dp) :: tol

  nfail = 0
  tol   = 100.0_dp*epsilon(1.0_dp)

  call check_constant(nfail)
  call check_linear(nfail)

  write(*,*) "============================================="
  if (nfail .eq. 0) then
     write(*,*) "all 2d checks passed"
  else
     write(*,*) nfail, " 2d checks FAILED"
     error stop
  end if

contains

  !===================================================================!
  ! Build a fresh assembler on a mesh file
  !===================================================================!

  subroutine make(meshfile, fvm)

    character(len=*)             , intent(in)  :: meshfile
    class(assembler), allocatable, intent(out) :: fvm

    class(gmsh_loader), allocatable :: gl
    class(mesh)       , allocatable :: grid

    allocate(gl, source = gmsh_loader(meshfile))
    allocate(grid, source = mesh(gl))
    allocate(fvm, source = assembler(grid))

  end subroutine make

  !===================================================================!
  ! Report the outcome of one check
  !===================================================================!

  subroutine report(name, ok, nfail)

    character(len=*), intent(in)    :: name
    logical         , intent(in)    :: ok
    integer         , intent(inout) :: nfail

    if (ok) then
       write(*,'(1x,a,a)') "PASS : ", name
    else
       write(*,'(1x,a,a)') "FAIL : ", name
       nfail = nfail + 1
    end if

  end subroutine report

  !===================================================================!
  ! Constant dirichlet on every side must give a uniform field
  !===================================================================!

  subroutine check_constant(nfail)

    integer, intent(inout) :: nfail

    class(assembler)         , allocatable :: fvm
    class(conjugate_gradient), allocatable :: cg
    real(dp)                 , allocatable :: x(:)

    call make("../square-10.msh", fvm)

    call fvm % set_dirichlet("BoundaryLeft"  , 5.0_dp)
    call fvm % set_dirichlet("BoundaryRight" , 5.0_dp)
    call fvm % set_dirichlet("BoundaryTop"   , 5.0_dp)
    call fvm % set_dirichlet("BoundaryBottom", 5.0_dp)

    allocate(cg, source = conjugate_gradient(fvm, 500, tol, 0))
    call cg % solve(x)

    call report("2d constant dirichlet -> u = 5", &
         & maxval(abs(x - 5.0_dp)) .lt. 1.0e-10_dp, nfail)

  end subroutine check_constant

  !===================================================================!
  ! Linear conduction across x with insulated top/bottom must give the
  ! exact field u = x (checked at the cell centres).
  !===================================================================!

  subroutine check_linear(nfail)

    integer, intent(inout) :: nfail

    class(assembler)         , allocatable :: fvm
    class(conjugate_gradient), allocatable :: cg
    real(dp)                 , allocatable :: x(:)
    real(dp)                               :: err
    integer                                :: icell

    call make("../square-10.msh", fvm)

    call fvm % set_dirichlet("BoundaryLeft" , 0.0_dp)
    call fvm % set_dirichlet("BoundaryRight", 1.0_dp)
    call fvm % set_neumann  ("BoundaryTop"   , 0.0_dp)
    call fvm % set_neumann  ("BoundaryBottom", 0.0_dp)

    allocate(cg, source = conjugate_gradient(fvm, 500, tol, 0))
    call cg % solve(x)

    err = 0.0_dp
    do icell = 1, fvm % grid % num_cells
       err = max(err, abs(x(icell) - fvm % grid % cell_centers(1, icell)))
    end do

    ! u = x is the exact discrete field; the residual is the iterative
    ! solve tolerance, so check it reproduces the linear field to ~1e-6
    call report("2d linear conduction -> u = x", err .lt. 1.0e-6_dp, nfail)

  end subroutine check_linear

end program test_2d
