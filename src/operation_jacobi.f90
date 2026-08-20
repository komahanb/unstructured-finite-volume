!=====================================================================!
! The jacobi iteration on the tower.
!
! The oldest idea in iterative solving: correct every cell by its
! own residual over its own diagonal, all cells at once,
!
!      x  <-  x + omega * (rhs - A x) / diag
!
! Everything it needs is inherited: matvec from the attached
! operation, the diagonal probed by colour, the norm from the
! reduction. This file states the iteration and nothing else.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_jacobi

  use iso_fortran_env  , only : dp => REAL64
  use operation_minimization, only : minimizer

  implicit none

  private
  public :: jacobi

  type, extends(minimizer) :: jacobi

     real(dp) :: omega = 1.0_dp

   contains

     procedure :: name => jacobi_name
     procedure :: solve

  end type jacobi

contains

  pure function jacobi_name(this) result(name)

    class(jacobi), intent(in) :: this
    character(len=:), allocatable :: name

    associate (u1 => this); end associate

    name = 'jacobi'

  end function jacobi_name


  subroutine solve(this, rhs, x, achieved)

    class(jacobi), intent(inout) :: this
    real(dp), intent(in)    :: rhs(:)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out)   :: achieved

    real(dp), allocatable :: d(:), y(:), r(:)
    real(dp) :: goal
    integer :: it, v

    call this % diagonal(d)

    ! A zero diagonal cannot correct its cell; leave that cell alone
    ! rather than divide by nothing.
    do v = 1, size(d)
       if (abs(d(v)) < tiny(1.0_dp)) d(v) = huge(1.0_dp)
    end do

    goal = this % tolerance * (1.0_dp + this % norm(rhs))

    do it = 1, this % max_iterations

       call this % matvec(x, y)
       r = rhs - y

       achieved = this % norm(r)
       if (achieved < goal) return

       x = x + this % omega * r / d

    end do

  end subroutine solve

end module operation_jacobi
