!=====================================================================!
! Tests for operation_dense_direct. dense_direct assembles its
! matrix by applying the attached operation's matvec to each basis
! vector, then eliminates with partial pivoting. Covered here:
! exact solutions of 1x1 and 2x2 stencil systems; the achieved
! residual, which must be the attached operation's norm of
! rhs - matvec(x); repeated solves on one attached operation; the
! stencil weights being left unmodified by a solve; a zero leading
! pivot handled by partial pivoting; a stencil compiled from an
! operation by evaluation on the standard basis, constants split
! from weights; and the transpose of a stencil, solved through the
! same minimizer.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_graph_dense_direct

  use iso_fortran_env     , only : dp => REAL64
  use field_calculus, only : field
  use field_stored  , only : stored_field
  use operation_stencil , only : stencil
  use operation_dense_direct, only : dense_direct

  implicit none

  type(stencil) :: one_by_one, two_by_two, pivoting
  type(stencil) :: compiled, affine, upper, transposed, twice
  type(stored_field) :: basis
  class(field), allocatable :: image
  type(dense_direct)     :: solver

  real(dp), allocatable :: x(:), w_before(:), w_after(:)
  real(dp), allocatable :: col1(:), col2(:), c(:)
  real(dp) :: achieved
  real(dp) :: amat(2,2), umat(2,2)
  integer  :: nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "graph dense direct minimizer suite"
  write(*,'(1x,a)') "============================================="

  !-------------------------------------------------------------------!
  ! 1x1: the statement [2] x = [6] on a one-vertex, one-edge graph.
  !-------------------------------------------------------------------!

  one_by_one = stencil([1], [1], [2.0_dp], [0.0_dp], 'one by one')

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

  two_by_two = stencil([1, 1, 2, 2], [1, 2, 1, 2], &
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

  call two_by_two % weights % real_vector(w_before)
  allocate(x(2)); x = 0.0_dp
  call solver % solve([4.0_dp, 7.0_dp], x, achieved)
  call two_by_two % weights % real_vector(w_after)
  call report(matches(w_after, w_before, 0.0_dp + tiny(1.0_dp)), &
       & "the solve leaves the stencil's weights unmodified", nfail)
  deallocate(x)

  !-------------------------------------------------------------------!
  ! Partial pivoting: A(1,1) = 0 and the system is nonsingular -
  ! A = [0 1; 1 1], x = [3, 5], b = [5, 8].
  !-------------------------------------------------------------------!

  pivoting = stencil([1, 2, 2], [2, 1, 2], &
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
  ! A stencil compiled from an operation by evaluation on the
  ! standard basis: compiled from two_by_two, its columns are the
  ! columns of A exactly and its constants are zero; compiled from
  ! an affine stencil, the constants are split from the weights.
  !-------------------------------------------------------------------!

  amat = reshape([2.0_dp, 1.0_dp, 1.0_dp, 3.0_dp], [2, 2])
  compiled = stencil(two_by_two, two_by_two % pattern, 2)
  call column_of(compiled, 1, col1)
  call column_of(compiled, 2, col2)
  call compiled % constants % real_vector(c)
  call report(matches([col1, col2], [amat(:, 1), amat(:, 2)], &
       &         0.0_dp + tiny(1.0_dp)) .and. all(c == 0.0_dp) .and. &
       & compiled % name() == 'two by two', &
       & "a stencil compiled from an operation carries its matrix exactly", nfail)

  affine = stencil([1, 1, 2, 2], [1, 2, 1, 2], &
       & [2.0_dp, 1.0_dp, 1.0_dp, 3.0_dp], [5.0_dp, -1.0_dp], 'affine')
  compiled = stencil(affine, affine % pattern, 2, 'compiled affine')
  call compiled % constants % real_vector(c)
  call column_of(compiled, 1, col1)
  call column_of(compiled, 2, col2)
  call report(matches(c, [5.0_dp, -1.0_dp], 0.0_dp + tiny(1.0_dp)) .and. &
       & matches([col1, col2], [7.0_dp, 0.0_dp, 6.0_dp, 2.0_dp], &
       &         0.0_dp + tiny(1.0_dp)), &
       & "compiling an affine operation splits the constant from the weights", &
       & nfail)

  !-------------------------------------------------------------------!
  ! The transpose of a stencil: for U = [2 1; 0 3] the transposed
  ! stencil applied to e_j returns row j of U, its constants are
  ! zero, transposing twice returns U, and U^T x = [2, 7] solved
  ! through dense_direct gives x = [1, 2].
  !-------------------------------------------------------------------!

  umat = reshape([2.0_dp, 0.0_dp, 1.0_dp, 3.0_dp], [2, 2])
  upper = stencil(umat, 'upper')
  transposed = upper % transpose()
  call column_of(transposed, 1, col1)
  call column_of(transposed, 2, col2)
  call transposed % constants % real_vector(c)
  call report(matches([col1, col2], [umat(1, :), umat(2, :)], &
       &         0.0_dp + tiny(1.0_dp)) .and. all(c == 0.0_dp) .and. &
       & transposed % name() == 'transpose of upper', &
       & "the transpose applied to e_j returns row j, constants dropped", nfail)

  twice = transposed % transpose()
  call column_of(twice, 1, col1)
  call column_of(twice, 2, col2)
  call report(matches([col1, col2], [umat(:, 1), umat(:, 2)], &
       &         0.0_dp + tiny(1.0_dp)), &
       & "transposing twice returns the stencil", nfail)

  call solver % attach(transposed, transposed % pattern, &
       & transposed % pattern % vertex_set(), &
       & transposed % pattern % num_vertices())
  allocate(x(2)); x = 0.0_dp
  call solver % solve([2.0_dp, 7.0_dp], x, achieved)
  call report(matches(x, [1.0_dp, 2.0_dp], 1.0e-12_dp) .and. &
       & achieved <= 1.0e-12_dp, &
       & "U^T x = [2, 7] through dense_direct gives x = [1, 2]", nfail)
  deallocate(x)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all dense direct minimizer checks passed"
  else
     error stop
  end if

contains

  !-------------------------------------------------------------------!
  ! Column j of a two-vertex stencil: the stencil applied to e_j on
  ! its own pattern.
  !-------------------------------------------------------------------!

  subroutine column_of(statement, j, column)

    type(stencil), intent(in) :: statement
    integer, intent(in) :: j
    real(dp), allocatable, intent(out) :: column(:)

    real(dp) :: e(2)

    e    = 0.0_dp
    e(j) = 1.0_dp
    basis = stored_field('e', statement % pattern % vertex_set(), 2)
    call basis % set_real_vector(e)
    call statement % apply(statement % pattern, [basis], image)
    call image % real_vector(column)

  end subroutine column_of

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
