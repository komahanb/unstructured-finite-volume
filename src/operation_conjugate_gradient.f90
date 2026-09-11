!=====================================================================!
! The conjugate gradient iteration on the minimizer hierarchy.
!
! For a symmetric positive operator, the optimal correction in every
! direction already searched, at the cost of one matvec and two
! inner products a step:
!
!      alpha = (r, r) / (p, A p)         the step length
!      x <- x + alpha p                  along this direction
!      r <- r - alpha A p                the residual
!      p <- r + beta p                   the next direction,
!                                        conjugate to all before it
!
! Everything the iteration needs is inherited: matvec from the attached
! operation, the inner product from the sum reduction with its
! measure, the norm from the norm reduction. This file states the
! iteration and nothing else.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_conjugate_gradient

  use util_precision  , only : dp
  use operation_minimization, only : minimizer, solve_result, SOLVE_BREAKDOWN, SOLVE_CONTINUE
  use util_tally, only : tally_record, linear_solves

  implicit none

  private
  public :: conjugate_gradient

  type, extends(minimizer) :: conjugate_gradient

   contains

     procedure :: name => conjugate_gradient_name
     procedure :: solve

  end type conjugate_gradient

contains

  pure function conjugate_gradient_name(this) result(name)

    class(conjugate_gradient), intent(in) :: this
    character(len=:), allocatable :: name

    associate (u1 => this); end associate

    name = 'conjugate gradient'

  end function conjugate_gradient_name


  subroutine solve(this, rhs, x, achieved)

    class(conjugate_gradient), intent(inout) :: this
    real(dp), intent(in)    :: rhs(:)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out)   :: achieved

    real(dp), allocatable :: r(:), p(:), ap(:)
    real(dp) :: rr, rr_next, alpha, beta, pap
    integer :: it
    type(solve_result) :: outcome

    call tally_record(linear_solves)

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

       alpha = rr / pap
       x = x + alpha * p
       r = r - alpha * ap

       achieved = this % norm(r)
       if (this % terminated(achieved, it)) then
          ! A recurrence can lose its last digits. Every returned
          ! residual is measured through the stated operator itself.
          call this % imbalance(rhs, x, r)
          achieved = this % norm(r)
          call this % record_result(achieved, it)
          outcome = this % result()
          if (outcome % reason /= SOLVE_CONTINUE) return
          rr = this % inner_product(r, r)
          p = r
          cycle
       end if

       rr_next = this % inner_product(r, r)
       beta    = rr_next / rr
       rr      = rr_next

       p = r + beta * p

    end do

  end subroutine solve

end module operation_conjugate_gradient
