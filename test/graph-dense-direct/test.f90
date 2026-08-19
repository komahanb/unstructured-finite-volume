!=====================================================================!
! Tests for class_graph_dense_direct. dense_direct assembles its
! matrix by applying the attached operation's matvec to each basis
! vector, then eliminates with partial pivoting. Covered here:
! exact solutions of 1x1 and 2x2 stencil systems; the achieved
! residual, which must be the attached operation's norm of
! rhs - matvec(x); repeated solves on one attached operation; the
! stencil weights being left unmodified by a solve; a zero leading
! pivot handled by partial pivoting;
! solve_dense_matrix_with_dense_direct, which lays a plain dense
! array on a stencil and solves it through the same interface; and
! dense_matrix_of, which rebuilds an operation's matrix column by
! column.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_graph_dense_direct

  use iso_fortran_env     , only : dp => REAL64
  use class_graph_stencil , only : stencil_operator
  use class_graph_dense_direct, only : dense_direct, &
       & solve_dense_matrix_with_dense_direct, dense_matrix_of

  implicit none

  type(stencil_operator) :: one_by_one, two_by_two, pivoting
  type(dense_direct)     :: solver

  real(dp), allocatable :: x(:), xa(:), w_before(:), w_after(:)
  real(dp), allocatable :: arebuilt(:,:)
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
       & "achieved equals the attached operation's norm of rhs - matvec(x)", &
       & nfail)

  !-------------------------------------------------------------------!
  ! Repeated solves on one attached operation: the second solve
  ! must not be affected by the first.
  !-------------------------------------------------------------------!

  x = 0.0_dp
  call solver % solve([2.0_dp, 1.0_dp], x, achieved)
  call report(matches(x, [1.0_dp, 0.0_dp], 1.0e-12_dp), &
       & "a second solve on one attached dense_direct is exact", nfail)
  deallocate(x)

  !-------------------------------------------------------------------!
  ! A solve must not modify the stencil's weights.
  !-------------------------------------------------------------------!

  call two_by_two % weights % get_real_vector(w_before)
  allocate(x(2)); x = 0.0_dp
  call solver % solve([4.0_dp, 7.0_dp], x, achieved)
  call two_by_two % weights % get_real_vector(w_after)
  call report(matches(w_after, w_before, 0.0_dp + tiny(1.0_dp)), &
       & "the solve leaves the stencil's weights unmodified", nfail)
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
  ! solve_dense_matrix_with_dense_direct: a plain dense array is
  ! laid on a stencil and solved through the same minimizer
  ! interface.
  !-------------------------------------------------------------------!

  amat = reshape([2.0_dp, 1.0_dp, 1.0_dp, 3.0_dp], [2, 2])
  call solve_dense_matrix_with_dense_direct(amat, [4.0_dp, 7.0_dp], &
       & 1.0e-14_dp, xa, achieved)

  call report(matches(xa, [1.0_dp, 2.0_dp], 1.0e-12_dp) .and. &
       & achieved <= 1.0e-12_dp, &
       & "the dense-array adapter solves through the same minimizer", nfail)

  !-------------------------------------------------------------------!
  ! dense_matrix_of: one apply per basis column. Applied to the
  ! stencil that carries A it must reproduce A exactly.
  !-------------------------------------------------------------------!

  call dense_matrix_of(two_by_two, two_by_two % pattern, 2, arebuilt)

  call report(size(arebuilt, 1) == 2 .and. size(arebuilt, 2) == 2 .and. &
       & matches([arebuilt(1,1), arebuilt(2,1), arebuilt(1,2), &
       &          arebuilt(2,2)], &
       &         [amat(1,1), amat(2,1), amat(1,2), amat(2,2)], &
       &         0.0_dp + tiny(1.0_dp)), &
       & "dense_matrix_of reassembles the stencil's matrix exactly", nfail)

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
