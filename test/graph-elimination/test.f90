!=====================================================================!
! Sparse elimination stores coefficients, not dependency paths, and
! every construction is admitted within a declared storage limit.
!
! The overlapping system: unknowns 1 and 2 retained, every later row
! eliminated, each reading the two rows before it and both retained
! columns twice over, so the dependency paths to column 1 grow as a
! Fibonacci recurrence while M has at most two columns per row.
!=====================================================================!

module elimination_fixture
  use iso_fortran_env, only : int64
  use util_precision, only : dp
  use operation_elimination, only : elimination, entry_bytes
  use operation_dense_direct, only : dense_direct
  use operation_gmres, only : gmres
  use operation_minimization, only : solve_result, SOLVE_STORAGE_EXCEEDED, absolute
  use operation_stencil, only : stencil
  implicit none
contains

  subroutine report(satisfied, description, failures)
    logical, intent(in) :: satisfied
    character(len=*), intent(in) :: description
    integer, intent(inout) :: failures
    if (satisfied) then
       write(*,'(a,a)') ' PASS : ', description
    else
       failures = failures + 1
       write(*,'(a,a)') ' FAIL : ', description
    end if
  end subroutine report

  subroutine overlapping_system(n, rows, columns, weights)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: rows(:), columns(:)
    real(dp), allocatable, intent(out) :: weights(:)
    integer :: count, i
    allocate(rows(8 * n), columns(8 * n), weights(8 * n))
    count = 0
    call entry(1, 1, 3.0_dp)
    call entry(1, n, -0.4_dp)
    call entry(2, 2, 2.0_dp)
    call entry(2, n - 1, 0.2_dp)
    do i = 3, n
       call entry(i, i, 2.0_dp)
       ! Two routes to retained column 1, including direct duplicates.
       call entry(i, 1, 0.1_dp)
       call entry(i, 1, 0.2_dp)
       call entry(i, 2, 0.2_dp)
       if (mod(i, 4) == 0) call entry(i, 2, -0.2_dp)
       if (i > 3) call entry(i, i - 1, -0.4_dp)
       if (i > 4) call entry(i, i - 2, -0.6_dp)
    end do
    rows = rows(1:count)
    columns = columns(1:count)
    weights = weights(1:count)
  contains
    subroutine entry(row, column, weight)
      integer, intent(in) :: row, column
      real(dp), intent(in) :: weight
      count = count + 1
      rows(count) = row
      columns(count) = column
      weights(count) = weight
    end subroutine entry
  end subroutine overlapping_system

  ! the right-hand side of the exact solution under the triples
  subroutine right_hand_side(rows, columns, weights, exact, rhs)
    integer, intent(in) :: rows(:), columns(:)
    real(dp), intent(in) :: weights(:), exact(:)
    real(dp), allocatable, intent(out) :: rhs(:)
    integer :: j
    allocate(rhs(size(exact)), source=0.0_dp)
    do j = 1, size(rows)
       rhs(rows(j)) = rhs(rows(j)) + weights(j) * exact(columns(j))
    end do
  end subroutine right_hand_side

  ! the counts of the coupling blocks, N and the eliminated rows: the
  ! substitution account less nnz(M), and the retained block's count
  subroutine block_counts(rows, columns, weights, eliminated, coupling, nkk)
    integer, intent(in) :: rows(:), columns(:)
    real(dp), intent(in) :: weights(:)
    logical, intent(in) :: eliminated(:)
    integer(int64), intent(out) :: coupling, nkk
    integer :: j
    coupling = int(count(eliminated), int64)
    nkk = 0
    do j = 1, size(rows)
       if (eliminated(rows(j)) .and. eliminated(columns(j))) then
          if (rows(j) /= columns(j) .and. weights(j) /= 0.0_dp) coupling = coupling + 1
       else if (eliminated(rows(j)) .neqv. eliminated(columns(j))) then
          coupling = coupling + 1
       else
          nkk = nkk + 1
       end if
    end do
  end subroutine block_counts

  subroutine check_substitution(n, failures)
    integer, intent(in) :: n
    integer, intent(inout) :: failures
    type(stencil) :: matrix
    type(elimination) :: solver
    type(dense_direct) :: reference
    integer, allocatable :: rows(:), columns(:)
    real(dp), allocatable :: weights(:), rhs(:)
    real(dp) :: exact(n), x(n), full(n), achieved, direct_achieved
    integer(int64) :: coupling, nkk
    integer :: i, repetition
    character(len=64) :: described

    call overlapping_system(n, rows, columns, weights)
    matrix = stencil(rows, columns, weights, spread(0.0_dp, 1, n), 'overlapping eliminated dependencies')
    allocate(solver % inner, source=dense_direct())
    solver % eliminated = [(i > 2, i = 1, n)]
    do repetition = 1, 2
       exact = [(1.0_dp + 0.1_dp * sin(real(i + repetition, dp)), i = 1, n)]
       call right_hand_side(rows, columns, weights, exact, rhs)
       call solver % state(matrix, matrix % pattern, matrix % pattern % vertex_set(), n)
       call reference % state(matrix, matrix % pattern, matrix % pattern % vertex_set(), n)
       x = 0.0_dp
       full = 0.0_dp
       call solver % solve(rhs, x, achieved)
       call reference % solve(rhs, full, direct_achieved)
       write(described, '(a,i0,a,i0)') 'sparse elimination, order ', n, ', statement ', repetition
       call report(maxval(abs(x - exact)) <= 1.0e-12_dp .and. maxval(abs(x - full)) <= 1.0e-12_dp &
            & .and. achieved <= 1.0e-12_dp, trim(described), failures)
    end do
    ! the substitution account counts at most nk coefficients of M
    ! per eliminated row, whatever the number of paths
    call block_counts(rows, columns, weights, solver % eliminated, coupling, nkk)
    write(*, '(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') '        accounts: input ', solver % storage % input, &
         & ' substitution ', solver % storage % substitution, ' temporary ', solver % storage % temporary, &
         & ' schur ', solver % storage % schur, ' factorisation ', solver % storage % factorisation, &
         & ' bytes ', solver % storage % bytes()
    call report(solver % storage % substitution <= coupling + int(n - 2, int64) * 2_int64 .and. &
         & solver % storage % substitution > coupling .and. solver % storage % schur <= 4_int64 .and. &
         & solver % storage % input == int(size(rows), int64) + nkk .and. &
         & solver % storage % temporary == int(n + 2, int64) + 2_int64 * solver % storage % schur .and. &
         & solver % storage % factorisation == 8_int64, &
         & 'the accounts count coefficients of M and of the complement, not dependency paths', failures)
  end subroutine check_substitution

  !===================================================================!
  ! Cancellation: row 4 reads column 1 directly and through row 3
  ! with opposite signs, so its coefficient of M vanishes exactly, is
  ! not stored, and the product J_KE M contributes no fill through it.
  !===================================================================!

  subroutine check_cancellation(failures)
    integer, intent(inout) :: failures
    integer, parameter :: n = 4
    integer :: rows(11), columns(11)
    real(dp) :: weights(11), exact(n), x(n), full(n), achieved, direct_achieved
    real(dp), allocatable :: rhs(:)
    type(stencil) :: matrix
    type(elimination) :: solver
    type(dense_direct) :: reference
    integer(int64) :: coupling, nkk
    integer :: i

    rows    = [1, 1, 2, 2, 1, 2, 3, 3, 4, 4, 4]
    columns = [1, 2, 1, 2, 4, 3, 3, 1, 4, 1, 3]
    weights = [3.0_dp, -0.5_dp, -0.25_dp, 2.0_dp, 0.7_dp, 0.3_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    matrix = stencil(rows, columns, weights, spread(0.0_dp, 1, n), 'cancelling dependencies')
    allocate(solver % inner, source=dense_direct())
    solver % eliminated = [.false., .false., .true., .true.]
    exact = [(1.0_dp + 0.1_dp * sin(real(i, dp)), i = 1, n)]
    call right_hand_side(rows, columns, weights, exact, rhs)
    call solver % state(matrix, matrix % pattern, matrix % pattern % vertex_set(), n)
    call reference % state(matrix, matrix % pattern, matrix % pattern % vertex_set(), n)
    x = 0.0_dp
    full = 0.0_dp
    call solver % solve(rhs, x, achieved)
    call reference % solve(rhs, full, direct_achieved)
    call block_counts(rows, columns, weights, solver % eliminated, coupling, nkk)
    call report(maxval(abs(x - exact)) <= 1.0e-13_dp .and. maxval(abs(x - full)) <= 1.0e-13_dp, &
         & 'a cancelling substitution row solves the whole system', failures)
    call report(solver % storage % substitution == coupling + 1_int64 .and. solver % storage % schur == nkk, &
         & 'an exactly cancelled coefficient of M is not stored and adds no fill to the complement', failures)
  end subroutine check_cancellation

  !===================================================================!
  ! Large true fill: three eliminated rows each reading every retained
  ! column, every retained row reading one of them: the complement is
  ! dense, nk^2 entries, and the account equals that count.
  !===================================================================!

  subroutine check_fill(failures)
    integer, intent(inout) :: failures
    integer, parameter :: nk = 12, ne = 3, n = nk + ne
    integer, allocatable :: rows(:), columns(:)
    real(dp), allocatable :: weights(:), rhs(:)
    real(dp) :: exact(n), x(n), full(n), achieved, direct_achieved
    type(stencil) :: matrix
    type(elimination) :: solver
    type(dense_direct) :: reference
    integer(int64) :: coupling, nkk
    integer :: i, e, count

    allocate(rows(4 * nk + ne * (nk + 1)), columns(4 * nk + ne * (nk + 1)), weights(4 * nk + ne * (nk + 1)))
    count = 0
    do i = 1, nk
       call entry(i, i, 4.0_dp)
       if (i > 1)  call entry(i, i - 1, -1.0_dp)
       if (i < nk) call entry(i, i + 1, -1.0_dp)
       call entry(i, nk + mod(i, ne) + 1, 0.05_dp)
    end do
    do e = 1, ne
       call entry(nk + e, nk + e, 2.0_dp)
       do i = 1, nk
          call entry(nk + e, i, 0.1_dp / real(e, dp))
       end do
    end do
    rows = rows(1:count)
    columns = columns(1:count)
    weights = weights(1:count)
    matrix = stencil(rows, columns, weights, spread(0.0_dp, 1, n), 'dense complement')
    allocate(solver % inner, source=dense_direct())
    solver % eliminated = [(i > nk, i = 1, n)]
    exact = [(1.0_dp + 0.1_dp * sin(real(i, dp)), i = 1, n)]
    call right_hand_side(rows, columns, weights, exact, rhs)
    call solver % state(matrix, matrix % pattern, matrix % pattern % vertex_set(), n)
    call reference % state(matrix, matrix % pattern, matrix % pattern % vertex_set(), n)
    x = 0.0_dp
    full = 0.0_dp
    call solver % solve(rhs, x, achieved)
    call reference % solve(rhs, full, direct_achieved)
    call block_counts(rows, columns, weights, solver % eliminated, coupling, nkk)
    call report(maxval(abs(x - exact)) <= 1.0e-12_dp .and. maxval(abs(x - full)) <= 1.0e-12_dp, &
         & 'a dense complement solves the whole system', failures)
    call report(solver % storage % schur == int(nk, int64) ** 2 .and. &
         & solver % storage % substitution == coupling + int(ne * nk, int64) .and. &
         & solver % storage % factorisation == 2_int64 * int(nk, int64) ** 2, &
         & 'the complement account is nk^2 for true fill and M has ne nk coefficients', failures)
  contains
    subroutine entry(row, column, weight)
      integer, intent(in) :: row, column
      real(dp), intent(in) :: weight
      count = count + 1
      rows(count) = row
      columns(count) = column
      weights(count) = weight
    end subroutine entry
  end subroutine check_fill

  !===================================================================!
  ! The exact boundary: the limit equal to the accounts' total admits
  ! the statement and changes no digit; one entry less refuses it
  ! before any allocation, reports SOLVE_STORAGE_EXCEEDED with the
  ! unknowns unchanged, and a sufficient limit afterwards solves again.
  !===================================================================!

  subroutine check_limit(failures)
    integer, intent(inout) :: failures
    integer, parameter :: n = 40
    integer, allocatable :: rows(:), columns(:)
    real(dp), allocatable :: weights(:), rhs(:)
    real(dp) :: exact(n), x(n), unbounded(n), seed(n), achieved, initial
    real(dp), allocatable :: seed_image(:)
    type(stencil) :: matrix
    type(elimination) :: solver
    type(solve_result) :: outcome
    integer(int64) :: total
    integer :: i

    call overlapping_system(n, rows, columns, weights)
    matrix = stencil(rows, columns, weights, spread(0.0_dp, 1, n), 'overlapping eliminated dependencies')
    allocate(solver % inner, source=dense_direct())
    solver % eliminated = [(i > 2, i = 1, n)]
    exact = [(1.0_dp + 0.1_dp * sin(real(i, dp)), i = 1, n)]
    call right_hand_side(rows, columns, weights, exact, rhs)
    seed = [(0.5_dp + 0.01_dp * real(i, dp), i = 1, n)]

    call solver % state(matrix, matrix % pattern, matrix % pattern % vertex_set(), n)
    total = solver % storage % total()
    unbounded = seed
    call solver % solve(rhs, unbounded, achieved)
    outcome = solver % result()
    call report(outcome % converged() .and. maxval(abs(unbounded - exact)) <= 1.0e-12_dp .and. &
         & solver % storage % bytes() == total * int(entry_bytes, int64), &
         & 'the default limit admits the statement and the bytes are the entries at their size', failures)

    solver % max_entries = int(total)
    call solver % state(matrix, matrix % pattern, matrix % pattern % vertex_set(), n)
    x = seed
    call solver % solve(rhs, x, achieved)
    outcome = solver % result()
    call report(outcome % converged() .and. all(x == unbounded) .and. solver % storage % total() == total, &
         & 'a limit equal to the total admits the statement with the same digits', failures)

    solver % max_entries = int(total) - 1
    call solver % state(matrix, matrix % pattern, matrix % pattern % vertex_set(), n)
    outcome = solver % result()
    call report(outcome % reason == SOLVE_STORAGE_EXCEEDED .and. outcome % failed() .and. &
         & solver % storage % total() == total .and. solver % storage % limit == total - 1_int64, &
         & 'one entry below the total refuses the statement, the accounts naming the requirement', failures)
    x = seed
    call solver % solve(rhs, x, achieved)
    outcome = solver % result()
    call right_hand_side(rows, columns, weights, seed, seed_image)
    initial = norm2(rhs - seed_image)
    call report(outcome % reason == SOLVE_STORAGE_EXCEEDED .and. outcome % failed() .and. all(x == seed) .and. &
         & abs(achieved - initial) <= 1.0e-14_dp * initial .and. outcome % iterations == 0 .and. &
         & outcome % description() == 'storage limit exceeded', &
         & 'the refused solve returns the initial residual with the unknowns unchanged', failures)

    solver % max_entries = int(total)
    call solver % state(matrix, matrix % pattern, matrix % pattern % vertex_set(), n)
    x = seed
    call solver % solve(rhs, x, achieved)
    outcome = solver % result()
    unbounded = seed
    solver % max_entries = huge(1)
    call solver % state(matrix, matrix % pattern, matrix % pattern % vertex_set(), n)
    call solver % solve(rhs, unbounded, achieved)
    call report(outcome % converged() .and. all(x == unbounded), &
         & 'a sufficient limit after a refusal solves again with the same digits', failures)
  end subroutine check_limit

  !===================================================================!
  ! Arithmetic overflow: 50 000 retained unknowns with a dense direct
  ! inner declare 2 nk^2 = 5 x 10^9 entries, beyond the largest count
  ! an index array addresses, and are refused before the 60 GB
  ! allocation whatever the limit; GMRES over the same retained set
  ! declares 31 nk and solves.
  !===================================================================!

  subroutine check_overflow(failures)
    integer, intent(inout) :: failures
    integer, parameter :: nk = 50000, n = nk + 1
    integer, allocatable :: rows(:), columns(:)
    real(dp), allocatable :: weights(:), exact(:), rhs(:), x(:)
    real(dp) :: achieved
    type(stencil) :: matrix
    type(elimination) :: solver
    type(gmres) :: krylov
    type(solve_result) :: outcome
    integer :: i

    rows    = [(i, i = 1, n), n]
    columns = [(i, i = 1, n), 1]
    weights = [(2.0_dp + 0.5_dp * sin(real(i, dp)), i = 1, n), 0.5_dp]
    matrix = stencil(rows, columns, weights, spread(0.0_dp, 1, n), 'one eliminated row')
    allocate(solver % inner, source=dense_direct())
    solver % eliminated = [(i > nk, i = 1, n)]
    call solver % state(matrix, matrix % pattern, matrix % pattern % vertex_set(), n)
    outcome = solver % result()
    call report(outcome % reason == SOLVE_STORAGE_EXCEEDED .and. &
         & solver % storage % factorisation == 2_int64 * int(nk, int64) ** 2 .and. &
         & solver % storage % total() > int(huge(1), int64) .and. solver % storage % limit == int(huge(1), int64), &
         & 'a factorisation requirement beyond the largest index count is refused before its allocation', failures)

    deallocate(solver % inner)
    krylov % tolerance = 1.0e-14_dp
    krylov % criterion = absolute
    krylov % max_iterations = 10
    allocate(solver % inner, source=krylov)
    allocate(exact(n), x(n))
    exact = [(1.0_dp + 0.1_dp * sin(real(i, dp)), i = 1, n)]
    call right_hand_side(rows, columns, weights, exact, rhs)
    call solver % state(matrix, matrix % pattern, matrix % pattern % vertex_set(), n)
    x = 0.0_dp
    call solver % solve(rhs, x, achieved)
    outcome = solver % result()
    call report(outcome % converged() .and. maxval(abs(x - exact)) <= 1.0e-12_dp .and. &
         & solver % storage % factorisation == 31_int64 * int(nk, int64) .and. &
         & solver % storage % schur == int(nk, int64) .and. solver % storage % substitution == 3_int64, &
         & 'the same system under a Krylov inner is admitted and solved', failures)
  end subroutine check_overflow

end module elimination_fixture

program test_graph_elimination
  use elimination_fixture, only : check_substitution, check_cancellation, check_fill, check_limit, check_overflow
  implicit none
  integer :: failures
  failures = 0
  call check_substitution(8, failures)
  call check_substitution(40, failures)
  call check_cancellation(failures)
  call check_fill(failures)
  call check_limit(failures)
  call check_overflow(failures)
  if (failures /= 0) error stop 'elimination: a check failed'
  write(*,'(a)') ' PASS : eliminated rows sum repeated columns before their dependants read them, within the limit'
end program test_graph_elimination
