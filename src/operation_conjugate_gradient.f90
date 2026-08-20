!=====================================================================!
! The conjugate gradient iteration on the tower.
!
! For a symmetric positive operator, the best correction in every
! direction already searched, at the cost of one matvec and two
! inner products a step:
!
!      alpha = (r, r) / (p, A p)         walk this far
!      x <- x + alpha p                  along this direction
!      r <- r - alpha A p                what remains
!      p <- r + beta p                   the next direction, kept
!                                        conjugate to all before it
!
! Everything it needs is inherited: matvec from the attached
! operation, the inner product from the sum reduction with its
! measure, the norm from the norm reduction. This file states the
! iteration and nothing else.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_conjugate_gradient

  use iso_fortran_env  , only : dp => REAL64
  use operation_minimization, only : minimizer

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

    real(dp), allocatable :: r(:), p(:), ap(:), y(:)
    real(dp) :: rr, rr_next, alpha, beta, pap, goal
    integer :: it

    call this % matvec(x, y)
    r = rhs - y
    p = r

    rr   = this % inner_product(r, r)
    goal = this % tolerance * (1.0_dp + this % norm(rhs))

    do it = 1, this % max_iterations

       achieved = this % norm(r)
       if (achieved < goal) return

       call this % matvec(p, ap)
       pap = this % inner_product(p, ap)
       if (abs(pap) < tiny(1.0_dp)) return

       alpha = rr / pap
       x = x + alpha * p
       r = r - alpha * ap

       rr_next = this % inner_product(r, r)
       beta    = rr_next / rr
       rr      = rr_next

       p = r + beta * p

    end do

    achieved = this % norm(r)

  end subroutine solve

end module operation_conjugate_gradient
