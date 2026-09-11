module elimination_fixture
  use util_precision, only : dp
  use operation_elimination, only : elimination
  use operation_dense_direct, only : dense_direct
  use operation_stencil, only : stencil
  implicit none
contains

  subroutine check_substitution(n, failures)
    integer, intent(in) :: n
    integer, intent(inout) :: failures
    type(stencil) :: matrix
    type(elimination) :: solver
    type(dense_direct) :: reference
    integer :: rows(8*n), columns(8*n), count, i, j, repetition
    real(dp) :: weights(8*n), exact(n), rhs(n), x(n), full(n), achieved, direct_achieved

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
    matrix = stencil(rows(1:count), columns(1:count), weights(1:count), &
         & spread(0.0_dp, 1, n), 'overlapping eliminated dependencies')
    allocate(solver % inner, source=dense_direct())
    solver % eliminated = [(i > 2, i = 1, n)]
    do repetition = 1, 2
       exact = [(1.0_dp + 0.1_dp * sin(real(i + repetition, dp)), i = 1, n)]
       rhs = 0.0_dp
       do j = 1, count
          rhs(rows(j)) = rhs(rows(j)) + weights(j) * exact(columns(j))
       end do
       call solver % state(matrix, matrix % pattern, matrix % pattern % vertex_set(), n)
       call reference % state(matrix, matrix % pattern, matrix % pattern % vertex_set(), n)
       x = 0.0_dp
       full = 0.0_dp
       call solver % solve(rhs, x, achieved)
       call reference % solve(rhs, full, direct_achieved)
       if (maxval(abs(x - exact)) > 1.0e-12_dp .or. &
            & maxval(abs(x - full)) > 1.0e-12_dp .or. achieved > 1.0e-12_dp) then
          failures = failures + 1
          write(*,'(a,i0,3es12.3)') ' FAIL : sparse elimination, order ', n, &
               & maxval(abs(x - exact)), maxval(abs(x - full)), achieved
       else
          write(*,'(a,i0,a,i0)') ' PASS : sparse elimination, order ', n, ', statement ', repetition
       end if
    end do

  contains

    subroutine entry(row, column, weight)
      integer, intent(in) :: row, column
      real(dp), intent(in) :: weight
      count = count + 1
      rows(count) = row
      columns(count) = column
      weights(count) = weight
    end subroutine entry

  end subroutine check_substitution
end module elimination_fixture

program test_graph_elimination
  use elimination_fixture, only : check_substitution
  implicit none
  integer :: failures
  failures = 0
  call check_substitution(8, failures)
  call check_substitution(40, failures)
  if (failures /= 0) error stop 'elimination: a sparse substitution changed the solution'
  write(*,'(a)') ' PASS : eliminated rows sum repeated columns before their dependants read them'
end program test_graph_elimination
