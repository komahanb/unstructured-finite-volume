!=====================================================================!
! Dense direct solve: Gaussian elimination with partial pivoting
! as a concrete minimizer. The matrix is not passed in; it is
! assembled by applying the stated operation's matvec to each
! basis vector,
!
!      A(:,j) = matvec(e_j),
!
! and A x = rhs is then eliminated. matvec and the norm come from
! the stated operation, and the achieved residual is measured
! through them; no matrix representation is owned here, since
! operation_stencil already provides one - a stencil built explicitly
! from an operation, or the transpose of one, is stated like any
! other operation.
!
! A direct solve is a single pass, so the tolerance and iteration
! limit inherited from the minimizer family are unused. The one
! numerical check is on the pivot: a pivot at or below
! singular_tolerance times the largest entry of the matrix cannot be
! divided by, being indistinguishable from zero at the matrix's own
! scale; the default is the spacing of the build's kind. By default
! that stops the program, the matrix being singular where the caller
! expected it not to be.
!
! A caller whose matrix is a tangent frozen at an intermediate
! iterate does not expect that, since singularity there is a property
! of the iterate and not a fault. Such a caller sets
! singular_reported, and a singular pivot then leaves the unknown
! unchanged and reports huge(1.0_dp) as the achieved residual - a
! value no completed elimination produces - for the outer iteration
! to act on.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_dense_direct

  use util_precision  , only : dp, spacing_at_one
  use operation_stencil, only : compile_matrix_from_action
  use operation_minimization , only : minimizer
  use util_factorisation, only : dense_factorisation
  use util_tally        , only : tally_record, linear_solves

  implicit none

  private
  public :: dense_direct

  type, extends(minimizer) :: dense_direct

     real(dp) :: singular_tolerance = spacing_at_one

     ! Whether a singular pivot is reported through the achieved
     ! residual instead of stopping the program. False by default,
     ! a singular matrix being a fault wherever the caller has not
     ! specified otherwise.
     logical  :: singular_reported = .false.

     ! The factors retained, and the version of the statement they
     ! belong to. A statement with the same version is not formed or
     ! factorised again; one with version zero always is. A statement
     ! with the same version but the opposite orientation is the
     ! transpose of the one retained, and is substituted against the
     ! same factors transposed.
     type(dense_factorisation), private :: factor
     integer                  , private :: retained_version  = 0
     logical                  , private :: retained_transposed = .false.

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

    real(dp), allocatable :: a(:,:), constant(:), r(:), solution(:)
    integer :: n
    logical :: kept

    if (this % singular_tolerance <= 0.0_dp) then
       error stop 'dense_direct: singular tolerance is positive'
    end if

    if (size(x) /= size(rhs)) then
       error stop 'dense_direct: solution size matches rhs'
    end if

    call tally_record(linear_solves)

    n = size(rhs)

    kept = this % action % version() /= 0 .and. &
         & this % action % version() == this % retained_version .and. &
         & this % factor % order() == n

    !----------------------------------------------------------------!
    ! Assemble the matrix - one matvec per basis vector, one column
    ! each - and factorise, unless the factors retained are already
    ! this statement's.
    !----------------------------------------------------------------!

    if (.not. kept) then
       call compile_matrix_from_action(this % action, this % graph, this % unknown_domain, &
            & this % num_unknowns, n, this % num_components, a, constant, stored=this % stored)
       call this % factor % factorise(a, this % singular_tolerance * maxval(abs(a)))
       this % retained_version  = this % action % version()
       this % retained_transposed = this % action % transpose_version()
    end if

    !----------------------------------------------------------------!
    ! Factorise, and substitute. A singular pivot is reported through
    ! the achieved residual where the caller requested that, and stops
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
         & transposed = this % action % transpose_version() .neqv. this % retained_transposed)
    x = solution

    !----------------------------------------------------------------!
    ! Measure the residual of the computed solution through the
    ! stated operation's matvec and norm.
    !----------------------------------------------------------------!

    call this % imbalance(rhs, x, r)
    achieved = this % norm(r)

  end subroutine dense_direct_solve

end module operation_dense_direct
