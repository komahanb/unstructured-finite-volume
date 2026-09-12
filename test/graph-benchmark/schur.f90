!=====================================================================!
! Elimination: substitution depth, repeated dependency paths, retained
! size and the actual nonzero counts, measured separately.
!
! Unknowns 1..nk are retained, nk+1..nk+ne eliminated. The retained
! block is tridiagonal (4, -1); retained row r reads ke eliminated
! columns; eliminated row i has diagonal 2, reads the p eliminated rows
! before it (weight -1/(2p), so I + N is well conditioned) and c
! retained columns, cycling through K. The number of dependency paths
! from row i to K grows as a p-step recurrence, exponential in ne for
! p >= 2, while the distinct retained columns are at most nk.
!
! The symbolic oracle recomputes the patterns of M = (I+N)^-1 J_EK,
! the uncombined product triples J_KE M and the combined Schur fill
! from the same input. The numerical oracle is the residual of the
! whole system from the input triples and, below dense_limit unknowns,
! the dense direct solution.
!
! usage: schur ne nk p c ke tolerance dense_limit
!=====================================================================!

program elimination_scaling

  use util_precision        , only : dp
  use operation_elimination , only : elimination
  use operation_dense_direct, only : dense_direct
  use operation_stencil     , only : stencil
  use benchmark_measurement , only : phase, begin_phase, end_phase, record, peak_rss_kilobytes, &
       &                             argument_integer, argument_real

  implicit none

  character(len=:), allocatable :: tokens
  character(len=128) :: written
  integer :: ne, nk, p, c, ke, n, count, i, j, q, r, e, num_failures, dense_limit
  integer :: nnz_input, nnz_n, nnz_ek, nnz_ke, nnz_kk, nnz_m, uncombined, schur_nnz, num_active
  integer, allocatable :: rows(:), columns(:), pattern_first(:), pattern_column(:), column_row(:), active(:)
  real(dp), allocatable :: weights(:), exact(:), rhs(:), x(:), residual(:), full(:), paths(:)
  real(dp) :: tolerance, achieved, dense_achieved, relative_residual, relative_error, dense_departure
  real(dp) :: paths_maximum, paths_total
  logical :: paths_exceeded
  type(stencil) :: matrix
  type(elimination) :: solver
  type(dense_direct) :: reference
  type(phase) :: measured

  ne = argument_integer(1, 160)
  nk = argument_integer(2, 8)
  p  = argument_integer(3, 2)
  c  = argument_integer(4, 2)
  ke = argument_integer(5, 2)
  tolerance = argument_real(6, 1.0e-9_dp)
  dense_limit = argument_integer(7, 400)
  if (ne < 1 .or. nk < 1 .or. p < 1 .or. c < 1 .or. ke < 1) error stop 'schur: every count is positive'
  n = nk + ne
  write(written, '(a,i0,a,i0,a,i0,a,i0,a,i0)') 'suite=schur ne=', ne, ' nk=', nk, ' p=', p, ' c=', c, ' ke=', ke
  tokens = trim(written)
  num_failures = 0

  ! ---- the input triples ----
  nnz_kk = nk + 2 * (nk - 1)
  nnz_ke = nk * ke
  nnz_ek = ne * c
  nnz_n = 0
  do i = 1, ne
     nnz_n = nnz_n + min(p, i - 1)
  end do
  nnz_input = nnz_kk + nnz_ke + ne + nnz_n + nnz_ek
  allocate(rows(nnz_input), columns(nnz_input), weights(nnz_input))
  count = 0
  do r = 1, nk
     call entry(r, r, 4.0_dp)
     if (r > 1)  call entry(r, r - 1, -1.0_dp)
     if (r < nk) call entry(r, r + 1, -1.0_dp)
     do j = 1, ke
        call entry(r, nk + mod((r - 1) * ke + j - 1, ne) + 1, 0.1_dp)
     end do
  end do
  do i = 1, ne
     e = nk + i
     call entry(e, e, 2.0_dp)
     do q = 1, min(p, i - 1)
        call entry(e, e - q, -0.5_dp / real(p, dp))
     end do
     do j = 1, c
        call entry(e, mod((i - 1) * c + j - 1, nk) + 1, 0.3_dp / real(c, dp))
     end do
  end do
  if (count /= nnz_input) error stop 'schur: the triple count is stated in advance'

  ! ---- the symbolic oracle ----
  ! patterns of M by eliminated row, in the natural substitution order
  allocate(pattern_first(ne + 1), pattern_column(ne * nk))
  allocate(column_row(nk), source=0)
  allocate(active(nk))
  pattern_first(1) = 1
  do i = 1, ne
     num_active = 0
     do j = 1, c
        call activate(mod((i - 1) * c + j - 1, nk) + 1, i)
     end do
     do q = 1, min(p, i - 1)
        do j = pattern_first(i - q), pattern_first(i - q + 1) - 1
           call activate(pattern_column(j), i)
        end do
     end do
     pattern_column(pattern_first(i):pattern_first(i) + num_active - 1) = active(1:num_active)
     pattern_first(i + 1) = pattern_first(i) + num_active
  end do
  nnz_m = pattern_first(ne + 1) - 1
  ! the uncombined product J_KE M and the combined Schur fill
  uncombined = nnz_kk
  schur_nnz = 0
  column_row = 0
  do r = 1, nk
     num_active = 0
     call activate(r, ne + r)
     if (r > 1)  call activate(r - 1, ne + r)
     if (r < nk) call activate(r + 1, ne + r)
     do j = 1, ke
        i = mod((r - 1) * ke + j - 1, ne) + 1
        uncombined = uncombined + pattern_first(i + 1) - pattern_first(i)
        do q = pattern_first(i), pattern_first(i + 1) - 1
           call activate(pattern_column(q), ne + r)
        end do
     end do
     schur_nnz = schur_nnz + num_active
  end do
  ! dependency paths from every eliminated row to K
  allocate(paths(ne))
  paths_exceeded = .false.
  do i = 1, ne
     paths(i) = 1.0_dp
     do q = 1, min(p, i - 1)
        paths(i) = paths(i) + paths(i - q)
     end do
     if (paths(i) > 1.0e300_dp) then
        paths(i) = 1.0e300_dp
        paths_exceeded = .true.
     end if
  end do
  paths_maximum = maxval(paths)
  paths_total = 0.0_dp
  do r = 1, nk
     do j = 1, ke
        paths_total = min(paths_total + paths(mod((r - 1) * ke + j - 1, ne) + 1), 1.0e300_dp)
     end do
  end do

  ! ---- the numerical statement ----
  allocate(exact(n), rhs(n), x(n), full(n), residual(n))
  exact = [(1.0_dp + 0.1_dp * sin(real(i, dp)), i = 1, n)]
  rhs = 0.0_dp
  do j = 1, count
     rhs(rows(j)) = rhs(rows(j)) + weights(j) * exact(columns(j))
  end do

  call begin_phase(measured)
  matrix = stencil(rows, columns, weights, spread(0.0_dp, 1, n), 'eliminated dependencies')
  call end_phase(measured)
  call record(tokens, 'input_stencil', measured)

  allocate(solver % inner, source=dense_direct())
  solver % eliminated = [(i > nk, i = 1, n)]
  call begin_phase(measured)
  call solver % state(matrix, matrix % pattern, matrix % pattern % vertex_set(), n)
  call end_phase(measured)
  call record(tokens, 'state', measured)

  x = 0.0_dp
  call begin_phase(measured)
  call solver % solve(rhs, x, achieved)
  call end_phase(measured)
  call record(tokens, 'solve', measured)

  residual = -rhs
  do j = 1, count
     residual(rows(j)) = residual(rows(j)) + weights(j) * x(columns(j))
  end do
  relative_residual = maxval(abs(residual)) / maxval(abs(rhs))
  relative_error = maxval(abs(x - exact)) / maxval(abs(exact))
  call verified(relative_residual <= tolerance, 'whole_system_residual_within_tolerance')
  call verified(relative_error <= tolerance, 'solution_within_tolerance_of_exact')

  dense_departure = -1.0_dp
  if (n <= dense_limit) then
     call begin_phase(measured)
     call reference % state(matrix, matrix % pattern, matrix % pattern % vertex_set(), n)
     full = 0.0_dp
     call reference % solve(rhs, full, dense_achieved)
     call end_phase(measured)
     call record(tokens, 'dense_reference', measured)
     dense_departure = maxval(abs(x - full)) / maxval(abs(full))
     call verified(dense_departure <= tolerance, 'agrees_with_dense_direct')
  end if

  write(*, '(a,1x,a,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0)') &
       & 'summary', tokens, 'n=', n, 'nnz_input=', nnz_input, 'nnz_kk=', nnz_kk, 'nnz_ke=', nnz_ke, &
       & 'nnz_ek=', nnz_ek, 'nnz_n=', nnz_n, 'nnz_m=', nnz_m, 'schur_uncombined=', uncombined, &
       & 'schur_nnz=', schur_nnz, 'stored_entries=', nnz_ke + nnz_ek + nnz_n + ne + schur_nnz
  write(*, '(a,1x,a,1x,a,es12.4e3,1x,a,es12.4e3,1x,a,l1,1x,a,es12.4e3,1x,a,es12.4e3,1x,a,es12.4e3,1x,a,es12.4e3,1x,a,i0)') &
       & 'numbers', tokens, 'paths_maximum=', paths_maximum, 'paths_total=', paths_total, &
       & 'paths_exceeded=', paths_exceeded, 'relative_residual=', relative_residual, &
       & 'relative_error=', relative_error, 'reported_residual=', achieved, 'dense_departure=', dense_departure, &
       & 'peak_rss_kilobytes=', peak_rss_kilobytes()
  if (num_failures > 0) error stop 'schur: a verification failed'

contains

  subroutine entry(row, column, weight)
    integer, intent(in) :: row, column
    real(dp), intent(in) :: weight
    count = count + 1
    rows(count) = row
    columns(count) = column
    weights(count) = weight
  end subroutine entry

  ! one marker per retained column: a column is active for a row once
  subroutine activate(column, row_version)
    integer, intent(in) :: column, row_version
    if (column_row(column) == row_version) return
    column_row(column) = row_version
    num_active = num_active + 1
    active(num_active) = column
  end subroutine activate

  subroutine verified(condition, name)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: name
    if (.not. condition) num_failures = num_failures + 1
    write(*, '(a,1x,a,1x,a,a,1x,a,l1)') 'verification', tokens, 'name=', name, 'satisfied=', condition
  end subroutine verified

end program elimination_scaling
