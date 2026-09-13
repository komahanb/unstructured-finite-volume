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
  use operation_minimization , only : minimizer, restrict, solve_result, SOLVE_SINGULAR, SOLVE_CONTINUE, SOLVE_EXHAUSTED
  use iso_fortran_env   , only : int64
  use util_factorisation, only : dense_factorisation
  use util_tally        , only : linear_solves, factorisations

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
     procedure :: restrict => dense_direct_restrict
     procedure :: storage_entries => dense_direct_storage_entries
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
  ! The factors retained belong to the whole statement and are not
  ! those of any restriction of it: the retained version is cleared so
  ! the selected statement is factorised when first solved.
  !===================================================================!

  subroutine dense_direct_restrict(this, selected)

    class(dense_direct), intent(inout) :: this
    integer            , intent(in)    :: selected(:)

    call restrict(this, selected)
    this % retained_version    = 0
    this % retained_transposed = .false.

  end subroutine dense_direct_restrict

  !===================================================================!
  ! The assembled matrix and its factorisation, n^2 entries each,
  ! both stored while the matrix is factorised.
  !===================================================================!

  pure integer(int64) function dense_direct_storage_entries(this, num_unknowns) result(entries)

    class(dense_direct), intent(in) :: this
    integer            , intent(in) :: num_unknowns

    associate (u1 => this); end associate
    entries = int(num_unknowns, int64) * int(num_unknowns, int64)
    if (entries <= huge(entries) - entries) then
       entries = 2_int64 * entries
    else
       entries = huge(entries)
    end if

  end function dense_direct_storage_entries

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
    logical :: factors_current
    type(solve_result) :: outcome
    character(len=250) :: message

    if (this % singular_tolerance <= 0.0_dp) then
       write(message,'(a,es12.5)') 'dense_direct: singular_tolerance must be positive; actual = ', &
            & this % singular_tolerance
       error stop trim(message)
    end if

    if (size(x) /= size(rhs)) then
       write(message,'(a,i0,a,i0)') 'dense_direct: size(x) must equal size(rhs); size(x) = ', &
            & size(x), ', size(rhs) = ', size(rhs)
       error stop trim(message)
    end if

    call this % record_event(linear_solves)

    call this % initialize_residual_history()
    call this % imbalance(rhs, x, r)
    achieved = this % norm(r)
    call this % record_residual_norm(achieved)
    if (achieved == 0.0_dp) then
       call this % record_result(achieved, 0)
       return
    end if

    n = size(rhs)

    factors_current = this % action % version() /= 0 .and. &
         & this % action % version() == this % retained_version .and. &
         & this % factor % order() == n

    !----------------------------------------------------------------!
    ! Assemble the matrix - one matvec per basis vector, one column
    ! each - and factorise, unless the factors retained are already
    ! this statement's.
    !----------------------------------------------------------------!

    if (.not. factors_current) then
       call compile_matrix_from_action(this % action, this % graph, this % unknown_domain, &
            & this % num_unknowns, n, this % num_components, a, constant, stored=this % stored)
       call this % record_event(factorisations)
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
          call this % record_result(this % norm(r), 0, SOLVE_SINGULAR)
          return
       end if
       error stop 'dense_direct: factorisation found a singular pivot, and singular_reported is &
            &false so the failure is not returned through achieved'
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
    call this % record_result(achieved, 1)
    outcome = this % result()
    if (outcome % reason == SOLVE_CONTINUE) call this % record_result(achieved, 1, SOLVE_EXHAUSTED)

  end subroutine dense_direct_solve

end module operation_dense_direct
