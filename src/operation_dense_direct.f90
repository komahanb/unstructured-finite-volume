!=====================================================================!
! Module that states the dense direct solve: Gaussian elimination with
! partial pivoting as a concrete minimizer. The matrix is not passed
! in; it is assembled by applying the stated operation's matvec to
! each basis vector,
!
!      A(:,j) = matvec(e_j),
!
!      A = [ matvec(e_1) | matvec(e_2) | ... | matvec(e_n) ]
!
! and A x = rhs is then eliminated. The matvec and the norm come from
! the stated operation, and the achieved residual is measured through
! them. No matrix representation is owned here, since
! operation_stencil already provides one: a stencil built explicitly
! from an operation, or the transpose of one, is stated like any other
! operation.
!
! A direct solve is a single pass, so the tolerance and the iteration
! limit inherited from the minimizer family are unused. The one
! numerical check is on the pivot. A pivot at or below
! singular_tolerance times the largest entry of its column cannot be
! divided by, being indistinguishable from zero at its column's own
! scale; the default is the spacing of the build's kind:
!
!   |pivot| <= singular_tolerance * max |A(:, j)|      singular
!
! A matrix whose columns differ in scale, such as the reactions of an
! index-three system, whose pivots fall as the square of the step, is
! singular only where a column is. By default a singular pivot stops
! the program, the matrix being singular where the caller expected it
! not to be.
!
! A caller whose matrix is a jacobian frozen at an intermediate
! iterate does not expect that, since singularity there is a property
! of the iterate and not a fault. Such a caller sets
! singular_reported, and a singular pivot then leaves the unknown
! unchanged and reports huge(1.0_dp) as the achieved residual, a value
! no completed elimination produces, for the outer iteration to act
! on.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_dense_direct

  use util_precision     , only : dp, spacing_at_one
  use operation_stencil  , only : compile_matrix_from_action, compile_entries_from_action
  use operation_minimizer, only : minimizer, restrict, solve_result, SOLVE_SINGULAR, SOLVE_CONTINUE, SOLVE_EXHAUSTED
  use iso_fortran_env    , only : int64
  use util_factorisation , only : dense_factorisation
  use util_tally         , only : linear_solves, factorisations

  implicit none

  private

  public :: dense_direct

  !-------------------------------------------------------------------!
  ! Type that models the dense direct solve, a member of the minimizer
  ! hierarchy
  !-------------------------------------------------------------------!

  type, extends(minimizer) :: dense_direct

     real(dp) :: singular_tolerance = spacing_at_one

     ! Whether a singular pivot is reported through the achieved
     ! residual instead of stopping the program. False by default, a
     ! singular matrix being a fault wherever the caller has not
     ! specified otherwise.
     logical  :: singular_reported = .false.

     ! Whether the matrix is equilibrated before it is factorised (see
     ! util_factorisation): for rows of different units
     logical  :: equilibrated = .false.

     ! Whether the matrix is permuted to its block triangular form and
     ! its diagonal blocks alone are factorised (see util_factorisation):
     ! for the step of a family of several stages, each reading the
     ! stages before it
     logical  :: block_triangular = .false.

     ! The factors retained, and the version of the statement they
     ! belong to. A statement with the same version is not formed or
     ! factorised again; one with version zero always is. A statement
     ! with the same version but the opposite orientation is the
     ! transpose of the one retained, and is substituted against the
     ! same factors transposed.
     type(dense_factorisation), private :: factor
     integer                  , private :: retained_version    = 0
     logical                  , private :: retained_transposed = .false.

   contains

     procedure :: name            => dense_direct_name
     procedure :: restrict        => dense_direct_restrict
     procedure :: storage_entries => dense_direct_storage_entries
     procedure :: solve           => dense_direct_solve

  end type dense_direct

  ! A dense direct solve stated by its options
  interface dense_direct
     module procedure dense_direct_of
  end interface dense_direct

contains

  !===================================================================!
  ! Returns the name of the solve
  !===================================================================!

  pure function dense_direct_name(this) result(name)

    class(dense_direct), intent(in)  :: this
    character(len=:)   , allocatable :: name

    ! Names this, which the name does not read, to keep the compiler
    ! from warning that the argument is unused
    associate (u1 => this); end associate

    name = 'dense direct'

  end function dense_direct_name

  !===================================================================!
  ! Restricts the solve to the selected unknowns. The factors retained
  ! belong to the whole statement and are not those of any restriction
  ! of it: the retained version is cleared so that the selected
  ! statement is factorised when first solved.
  !===================================================================!

  subroutine dense_direct_restrict(this, selected)

    class(dense_direct), intent(inout) :: this
    integer            , intent(in)    :: selected(:)

    call restrict(this, selected)
    this % retained_version    = 0
    this % retained_transposed = .false.

  end subroutine dense_direct_restrict

  !===================================================================!
  ! Returns the entries stored while the matrix is factorised: the
  ! assembled matrix and its factorisation, n^2 entries each, the count
  ! set to huge where 2 n^2 would overflow
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
  ! Solves A x = rhs where A(:, j) = matvec(e_j), and returns the norm
  ! of the residual achieved:
  !
  !   r = rhs - A x            zero: done
  !        |
  !   factors current? ------ yes ------------------------+
  !        | no                                           |
  !        v                                              |
  !   block triangular, its entries stated:               |
  !        factorise the entries                          |
  !   otherwise:                                          |
  !        assemble A by columns, factorise               |
  !        |                                              |
  !        v                                              v
  !   singular? -- reported: achieved = huge; else stop
  !        |
  !   substitute, transposed where the orientation differs;
  !   measure the residual through the stated operation
  !
  ! Checks, each stopping the program: singular_tolerance must be
  ! positive, size(x) must equal size(rhs) because x is written in
  ! place, and every pivot must exceed singular_tolerance times the
  ! largest entry of its column.
  !===================================================================!

  subroutine dense_direct_solve(this, rhs, x, achieved)

    class(dense_direct), intent(inout) :: this
    real(dp)           , intent(in)    :: rhs(:)
    real(dp)           , intent(inout) :: x(:)
    real(dp)           , intent(out)   :: achieved

    real(dp), allocatable :: a(:,:), constant(:), r(:), solution(:), weights(:)
    integer , allocatable :: rows(:), columns(:)
    integer               :: n
    logical               :: factors_current, entries_stated
    type(solve_result)    :: outcome
    character(len=250)    :: message

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

    !-----------------------------------------------------------------!
    ! Assemble the matrix, one matvec per basis vector and one column
    ! each, and factorise, unless the factors retained are already this
    ! statement's.
    !-----------------------------------------------------------------!

    ! Block triangular, a stencil is factorised from its entries, with
    ! no matrix of the order of the system formed
    entries_stated = .false.
    if (.not. factors_current .and. this % block_triangular) then
       call compile_entries_from_action(this % action, n, rows, columns, weights, entries_stated)
    end if
    if (entries_stated) then
       call this % record_event(factorisations)
       call this % factor % factorise_entries(n, rows, columns, weights, this % singular_tolerance, &
            & this % equilibrated)
       this % retained_version    = this % action % version()
       this % retained_transposed = this % action % transpose_version()
    else if (.not. factors_current) then
       call compile_matrix_from_action(this % action, this % graph, this % unknown_domain, &
            & this % num_unknowns, n, this % num_components, a, constant, stored=this % stored)
       call this % record_event(factorisations)
       ! The scale of each column bounds its pivots; equilibrated, every
       ! column's is one
       if (this % equilibrated) then
          call this % factor % factorise(a, this % singular_tolerance, equilibrated=.true., &
               & block_triangular=this % block_triangular)
       else
          call this % factor % factorise(a, this % singular_tolerance, maxval(abs(a), dim=1), &
               & equilibrated=.false., block_triangular=this % block_triangular)
       end if
       this % retained_version    = this % action % version()
       this % retained_transposed = this % action % transpose_version()
    end if

    !-----------------------------------------------------------------!
    ! Substitute. A singular pivot is reported through the achieved
    ! residual where the caller requested that, and stops the program
    ! otherwise.
    !-----------------------------------------------------------------!

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

    !-----------------------------------------------------------------!
    ! Measure the residual of the computed solution through the stated
    ! operation's matvec and norm.
    !-----------------------------------------------------------------!

    call this % imbalance(rhs, x, r)
    achieved = this % norm(r)
    call this % record_result(achieved, 1)
    outcome = this % result()
    if (outcome % reason == SOLVE_CONTINUE) call this % record_result(achieved, 1, SOLVE_EXHAUSTED)

  end subroutine dense_direct_solve

  !===================================================================!
  ! Creates a dense direct solve stated by its options: equilibrated,
  ! its rows and columns scaled to a unit largest entry before the
  ! factorisation; block_triangular, its diagonal blocks alone
  ! factorised
  !===================================================================!

  function dense_direct_of(equilibrated, block_triangular) result(this)
    logical, intent(in), optional :: equilibrated, block_triangular
    type(dense_direct)            :: this
    if (present(equilibrated)) this % equilibrated = equilibrated
    if (present(block_triangular)) this % block_triangular = block_triangular
  end function dense_direct_of

end module operation_dense_direct
