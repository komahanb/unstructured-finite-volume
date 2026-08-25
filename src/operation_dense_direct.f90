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
! singular_tolerance times the largest entry of the matrix cannot be
! divided by, being indistinguishable from zero at the matrix's own
! scale; the default is the spacing of the build's kind. By default
! that stops the program, the matrix being singular where the caller
! expected it not to be.
!
! A caller whose matrix is a tangent frozen at an intermediate
! iterate expects no such thing, since singularity there is a fact
! about the iterate and not a fault. Such a caller sets
! singular_reported, and a singular pivot then leaves the unknown
! unchanged and reports huge(1.0_dp) as the achieved residual - a
! value no completed elimination produces - for the outer iteration
! to act on.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_dense_direct

  use util_precision  , only : dp, spacing_at_one
  use operation_minimization , only : minimizer
  use util_factorisation, only : dense_factorisation

  implicit none

  private
  public :: dense_direct

  type, extends(minimizer) :: dense_direct

     real(dp) :: singular_tolerance = spacing_at_one

     ! Whether a singular pivot is reported through the achieved
     ! residual instead of stopping the program. False by default,
     ! a singular matrix being a fault wherever the caller has not
     ! said otherwise.
     logical  :: singular_reported = .false.

     ! The factors kept, and the stamp of the statement they belong
     ! to. A statement stamped the same is not formed or factorised
     ! again; one stamped zero always is. A stamp of the opposite sign
     ! is the transpose of the statement kept, and is substituted
     ! against the same factors the other way round.
     type(dense_factorisation), private :: factor
     integer                  , private :: kept_stamp = 0

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
  ! and every pivot must exceed singular_tolerance times the largest
  ! entry of the matrix.
  !===================================================================!

  subroutine dense_direct_solve(this, rhs, x, achieved)

    class(dense_direct), intent(inout) :: this
    real(dp), intent(in)    :: rhs(:)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out)   :: achieved

    real(dp), allocatable :: a(:,:), basis(:), y(:), solution(:)
    integer :: n, j
    logical :: kept

    if (this % singular_tolerance <= 0.0_dp) then
       error stop 'dense_direct: singular tolerance is positive'
    end if

    if (size(x) /= size(rhs)) then
       error stop 'dense_direct: solution size matches rhs'
    end if

    n = size(rhs)

    kept = this % action % stamp() /= 0 .and. &
         & abs(this % action % stamp()) == abs(this % kept_stamp) .and. &
         & this % factor % order() == n

    !----------------------------------------------------------------!
    ! Assemble the matrix - one matvec per basis vector, one column
    ! each - and factorise, unless the factors kept are this
    ! statement's already.
    !----------------------------------------------------------------!

    if (.not. kept) then
       allocate(a(n, n), basis(n))
       do j = 1, n
          basis    = 0.0_dp
          basis(j) = 1.0_dp
          call this % matvec(basis, y)
          a(:, j) = y
       end do
       call this % factor % factorise(a, this % singular_tolerance * maxval(abs(a)))
       this % kept_stamp = this % action % stamp()
    end if

    !----------------------------------------------------------------!
    ! Factorise, and substitute. A singular pivot is reported through
    ! the achieved residual where the caller asked for that, and stops
    ! the program otherwise.
    !----------------------------------------------------------------!

    if (this % factor % singular()) then
       if (this % singular_reported) then
          achieved = huge(1.0_dp)
          return
       end if
       error stop 'dense_direct: the pivot is singular'
    end if

    call this % factor % substitute(rhs, solution, &
         & transposed = (this % action % stamp() < 0) .neqv. (this % kept_stamp < 0))
    x = solution

    !----------------------------------------------------------------!
    ! Measure the residual of the computed solution through the
    ! attached operation's matvec and norm.
    !----------------------------------------------------------------!

    call this % matvec(x, y)
    achieved = this % norm(rhs - y)

  end subroutine dense_direct_solve

end module operation_dense_direct
