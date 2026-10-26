!=====================================================================!
! Module that states the conjugate gradient iteration on the
! minimizer hierarchy.
!
! For a symmetric positive operator A, each step makes the residual
! orthogonal to every direction searched so far, at the cost of one
! product with A and two inner products:
!
!      alpha = (r, r) / (p, A p)         the step length
!      x <- x + alpha p                  along this direction
!      r <- r - alpha A p                the residual
!      p <- r + beta p                   the next direction,
!                                        conjugate to all before it
!
!   x0 ---- alpha0 p0 ----> x1 ---- alpha1 p1 ----> x2 ---- ...
!   p0 = r0                 p1 = r1 + beta0 p0      (p1, A p0) = 0
!
! Everything else is inherited: the product with A from the attached
! operation, the inner product from the sum reduction with its
! measure, and the norm from the norm reduction. This module states
! the iteration and nothing more.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_conjugate_gradient

  use util_precision     , only : dp
  use operation_minimizer, only : minimizer, solve_result, SOLVE_BREAKDOWN, SOLVE_CONTINUE
  use util_tally         , only : linear_solves

  implicit none

  private

  public :: conjugate_gradient

  !-------------------------------------------------------------------!
  ! Type that models the conjugate gradient iteration, a member of the
  ! minimizer hierarchy
  !-------------------------------------------------------------------!

  type, extends(minimizer) :: conjugate_gradient

   contains

     procedure :: name => conjugate_gradient_name
     procedure :: solve

  end type conjugate_gradient

contains

  !===================================================================!
  ! Returns the name of the iteration
  !===================================================================!

  pure function conjugate_gradient_name(this) result(name)

    class(conjugate_gradient), intent(in)  :: this
    character(len=:)         , allocatable :: name

    ! Names this, which the name does not read, to keep the compiler
    ! from warning that the argument is unused
    associate (u1 => this); end associate

    name = 'conjugate gradient'

  end function conjugate_gradient_name

  !===================================================================!
  ! Solves A x = rhs from the given x, and returns the norm of the
  ! residual achieved.
  !
  ! When the recurrence's residual shows the iteration has stopped, the
  ! residual is measured again through the operator itself, since a
  ! recurrence can lose its last digits:
  !
  !   r = r - alpha A p
  !       |
  !   stopped? ---- no ----> the next direction
  !       | yes
  !       v
  !   r = rhs - A x, measured
  !       |
  !   stopped? ---- yes ---> return
  !       | no: the recurrence had drifted
  !       v
  !   restart from r: p = r
  !
  ! The iteration also returns, as a breakdown, when (p, A p) or
  ! (r, r) is at or below tiny(1.0_dp), the smallest positive normal
  ! number.
  !===================================================================!

  subroutine solve(this, rhs, x, achieved)

    class(conjugate_gradient), intent(inout) :: this
    real(dp)                 , intent(in)    :: rhs(:)
    real(dp)                 , intent(inout) :: x(:)
    real(dp)                 , intent(out)   :: achieved

    real(dp), allocatable :: r(:), p(:), ap(:)
    real(dp)              :: rr, rr_next, alpha, beta, pap
    integer               :: it
    type(solve_result)    :: outcome

    call this % record_event(linear_solves)

    ! The first residual is the first direction
    call this % imbalance(rhs, x, r)
    p = r

    rr   = this % inner_product(r, r)
    call this % initialize_residual_history()
    achieved = this % norm(r)
    if (this % terminated(achieved, 0)) return

    do it = 1, this % max_iterations

       call this % matvec(p, ap)
       pap = this % inner_product(p, ap)
       if (abs(pap) <= tiny(1.0_dp) .or. rr <= tiny(1.0_dp)) then
          call this % imbalance(rhs, x, r)
          achieved = this % norm(r)
          call this % record_result(achieved, it - 1, SOLVE_BREAKDOWN)
          return
       end if

       ! The step along p
       alpha = rr / pap
       x = x + alpha * p
       r = r - alpha * ap

       achieved = this % norm(r)
       if (this % terminated(achieved, it)) then
          ! Every returned residual is measured through the stated
          ! operator itself
          call this % imbalance(rhs, x, r)
          achieved = this % norm(r)
          call this % record_result(achieved, it)
          outcome = this % result()
          if (outcome % reason /= SOLVE_CONTINUE) return
          rr = this % inner_product(r, r)
          p = r
          cycle
       end if

       ! The next direction, conjugate to the ones before it
       rr_next = this % inner_product(r, r)
       beta    = rr_next / rr
       rr      = rr_next

       p = r + beta * p

    end do

  end subroutine solve

end module operation_conjugate_gradient
