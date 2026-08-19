!=====================================================================!
! The dense direct minimizer suite: one more concretion of the
! graph minimization tower, proven over stencil statements - a
! matrix is a graph with weights on its edges, and dense_direct
! eliminates whatever the attached operation's matvec answers,
! probed one basis direction per column. Exact solutions on 1x1
! and 2x2 systems, survival of a zero leading pivot by partial
! pivoting, repeated solves on one attached value, an achieved
! norm the tower itself certifies, an unmutated stencil, and the
! dense-array adapter door the GTI drivers walk through.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_graph_dense_direct

  use iso_fortran_env     , only : dp => REAL64
  use class_graph_stencil , only : stencil_operator
  use class_graph_dense_direct, only : dense_direct, &
       & solve_dense_matrix_with_dense_direct

  implicit none

  type(stencil_operator) :: one_by_one, two_by_two, pivoting
  type(dense_direct)     :: solver

  real(dp), allocatable :: x(:), xa(:), w_before(:), w_after(:)
  real(dp) :: achieved
  real(dp) :: amat(2,2)
  integer  :: nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "graph dense direct minimizer suite"
  write(*,'(1x,a)') "============================================="

  !-------------------------------------------------------------------!
  ! 1x1: the statement [2] x = [6] on a one-vertex, one-edge graph.
  !-------------------------------------------------------------------!

  one_by_one = stencil_operator([1], [1], [2.0_dp], [0.0_dp], 'one by one')

  call solver % attach(one_by_one, one_by_one % pattern, &
       & one_by_one % pattern % vertex_set(), &
       & one_by_one % pattern % num_vertices())

  allocate(x(1)); x = 0.0_dp
  call solver % solve([6.0_dp], x, achieved)

  call report(matches(x, [3.0_dp], 1.0e-12_dp) .and. achieved <= 1.0e-12_dp, &
       & "1x1 stencil: [2] x = [6] gives x = [3]", nfail)
  deallocate(x)

  !-------------------------------------------------------------------!
  ! 2x2: A = [2 1; 1 3] as four weighted edges, x = [1, 2],
  ! b = [4, 7].
  !-------------------------------------------------------------------!

  two_by_two = stencil_operator([1, 1, 2, 2], [1, 2, 1, 2], &
       & [2.0_dp, 1.0_dp, 1.0_dp, 3.0_dp], [0.0_dp, 0.0_dp], 'two by two')

  call solver % attach(two_by_two, two_by_two % pattern, &
       & two_by_two % pattern % vertex_set(), &
       & two_by_two % pattern % num_vertices())

  allocate(x(2)); x = 0.0_dp
  call solver % solve([4.0_dp, 7.0_dp], x, achieved)

  call report(matches(x, [1.0_dp, 2.0_dp], 1.0e-12_dp), &
       & "2x2 stencil: A x = b recovers x = [1, 2] exactly", nfail)

  call report(achieved <= 1.0e-12_dp, &
       & "achieved is the tower's own norm of rhs - matvec(x), small", nfail)

  !-------------------------------------------------------------------!
  ! Repeated solves on one attached value: no memory of the first
  ! carries into the second.
  !-------------------------------------------------------------------!

  x = 0.0_dp
  call solver % solve([2.0_dp, 1.0_dp], x, achieved)
  call report(matches(x, [1.0_dp, 0.0_dp], 1.0e-12_dp), &
       & "repeated solves on one attached dense_direct are lawful", nfail)
  deallocate(x)

  !-------------------------------------------------------------------!
  ! The stencil is read, never written: its weights survive the
  ! solves untouched.
  !-------------------------------------------------------------------!

  call two_by_two % weights % get_real_vector(w_before)
  allocate(x(2)); x = 0.0_dp
  call solver % solve([4.0_dp, 7.0_dp], x, achieved)
  call two_by_two % weights % get_real_vector(w_after)
  call report(matches(w_after, w_before, 0.0_dp + tiny(1.0_dp)), &
       & "the stencil's weights survive the solve unmutated", nfail)
  deallocate(x)

  !-------------------------------------------------------------------!
  ! Partial pivoting: A(1,1) = 0 and the system is nonsingular -
  ! A = [0 1; 1 1], x = [3, 5], b = [5, 8].
  !-------------------------------------------------------------------!

  pivoting = stencil_operator([1, 2, 2], [2, 1, 2], &
       & [1.0_dp, 1.0_dp, 1.0_dp], [0.0_dp, 0.0_dp], 'zero leading pivot')

  call solver % attach(pivoting, pivoting % pattern, &
       & pivoting % pattern % vertex_set(), &
       & pivoting % pattern % num_vertices())

  allocate(x(2)); x = 0.0_dp
  call solver % solve([5.0_dp, 8.0_dp], x, achieved)

  call report(matches(x, [3.0_dp, 5.0_dp], 1.0e-12_dp), &
       & "a zero leading pivot still solves exactly, by partial pivoting", nfail)
  deallocate(x)

  !-------------------------------------------------------------------!
  ! The adapter door: a plain dense array laid on a stencil and
  ! solved through the same minimizer face - the seat the GTI
  ! drivers call.
  !-------------------------------------------------------------------!

  amat = reshape([2.0_dp, 1.0_dp, 1.0_dp, 3.0_dp], [2, 2])
  call solve_dense_matrix_with_dense_direct(amat, [4.0_dp, 7.0_dp], &
       & 1.0e-14_dp, xa, achieved)

  call report(matches(xa, [1.0_dp, 2.0_dp], 1.0e-12_dp) .and. &
       & achieved <= 1.0e-12_dp, &
       & "the dense-array adapter solves through the same tower", nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all dense direct minimizer checks passed"
  else
     error stop
  end if

contains

  pure function matches(values, expected, tolerance) result(ok)

    real(dp), intent(in) :: values(:), expected(:)
    real(dp), intent(in) :: tolerance
    logical :: ok

    ok = size(values) == size(expected)
    if (ok) ok = all(abs(values - expected) <= tolerance)

  end function matches

  subroutine report(ok, label, nfail)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: label
    integer, intent(inout) :: nfail
    if (ok) then
       write(*,'(1x,a,a)') "PASS : ", label
    else
       write(*,'(1x,a,a)') "FAIL : ", label
       nfail = nfail + 1
    end if
  end subroutine report

end program test_graph_dense_direct
