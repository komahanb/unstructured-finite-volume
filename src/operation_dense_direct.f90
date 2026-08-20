!=====================================================================!
! Dense direct solve: Gaussian elimination with partial pivoting
! as a concrete minimizer. The matrix is not passed in; it is
! assembled by applying the attached operation's matvec to each
! basis vector,
!
!      A(:,j) = matvec(e_j),
!
! and A x = rhs is then eliminated. matvec and the norm come from
! the attached operation, and the achieved residual is measured
! through them; no matrix representation is owned here, since
! operation_stencil already provides one - a stencil compiled from
! an operation, or the transpose of one, is attached like any other
! operation.
!
! A direct solve is a single pass, so the tolerance and iteration
! budget inherited from the minimizer family are unused. The one
! numerical check is on the pivot: a pivot at or below
! singular_tolerance stops the program, because continuing would
! divide by a value indistinguishable from zero.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_dense_direct

  use iso_fortran_env    , only : dp => REAL64
  use operation_minimization , only : minimizer

  implicit none

  private
  public :: dense_direct

  type, extends(minimizer) :: dense_direct

     real(dp) :: singular_tolerance = 1.0e-14_dp

   contains

     procedure :: name  => dense_direct_name
     procedure :: solve => dense_direct_solve

  end type dense_direct

contains

  pure function dense_direct_name(this) result(name)

    class(dense_direct), intent(in) :: this
    character(len=:), allocatable :: name

    associate (u1 => this); end associate

    name = 'dense direct'

  end function dense_direct_name

  !===================================================================!
  ! Solve A x = rhs where A(:, j) = matvec(e_j). Checks, each
  ! stopping the program: singular_tolerance must be positive,
  ! size(x) must equal size(rhs) because x is written in place,
  ! and every pivot must exceed singular_tolerance.
  !===================================================================!

  subroutine dense_direct_solve(this, rhs, x, achieved)

    class(dense_direct), intent(inout) :: this
    real(dp), intent(in)    :: rhs(:)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out)   :: achieved

    real(dp), allocatable :: a(:,:), basis(:), y(:), r(:), row(:)
    real(dp) :: swap_value, factor
    integer :: n, j, k, p, i

    if (this % singular_tolerance <= 0.0_dp) then
       error stop 'dense_direct: singular tolerance is positive'
    end if

    if (size(x) /= size(rhs)) then
       error stop 'dense_direct: solution size matches rhs'
    end if

    n = size(rhs)

    !----------------------------------------------------------------!
    ! Assemble the matrix: one matvec per basis vector, one column
    ! each.
    !----------------------------------------------------------------!

    allocate(a(n, n), basis(n))

    do j = 1, n
       basis    = 0.0_dp
       basis(j) = 1.0_dp
       call this % matvec(basis, y)
       a(:, j) = y
    end do

    !----------------------------------------------------------------!
    ! Gaussian elimination with partial pivoting, on the assembled
    ! matrix and a copy of the right-hand side.
    !----------------------------------------------------------------!

    r = rhs
    allocate(row(n))

    do k = 1, n

       p = k - 1 + maxloc(abs(a(k:n, k)), dim=1)

       if (abs(a(p, k)) <= this % singular_tolerance) then
          error stop 'dense_direct: pivot is nonsingular'
       end if

       if (p /= k) then
          row        = a(k, :)
          a(k, :)    = a(p, :)
          a(p, :)    = row
          swap_value = r(k)
          r(k)       = r(p)
          r(p)       = swap_value
       end if

       do i = k + 1, n
          factor    = a(i, k) / a(k, k)
          a(i, k:n) = a(i, k:n) - factor * a(k, k:n)
          r(i)      = r(i) - factor * r(k)
       end do

    end do

    do k = n, 1, -1
       x(k) = (r(k) - dot_product(a(k, k+1:n), x(k+1:n))) / a(k, k)
    end do

    !----------------------------------------------------------------!
    ! Measure the residual of the computed solution through the
    ! attached operation's matvec and norm.
    !----------------------------------------------------------------!

    call this % matvec(x, y)
    achieved = this % norm(rhs - y)

  end subroutine dense_direct_solve

end module operation_dense_direct
