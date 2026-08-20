!=====================================================================!
! The gauss-seidel iteration on the tower, swept by colour.
!
! Jacobi corrects every cell from the old state; gauss-seidel lets
! each correction see the ones already made. On a graph the safe
! order is the colouring the sweep_order delegation already answers:
! all cells of one colour share no face, so a whole colour updates
! at once, each colour seeing every colour before it,
!
!      for each colour:  r = rhs - A x       (x already partly new)
!                        x <- x + omega * r / diag   on that colour
!
! SOR is not another solver. It is this one at omega away from one -
! a parameter, absorbed, exactly as the admission law orders.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_gauss_seidel

  use iso_fortran_env  , only : dp => REAL64
  use operation_minimization, only : minimizer

  implicit none

  private
  public :: gauss_seidel

  type, extends(minimizer) :: gauss_seidel

     real(dp) :: omega = 1.0_dp

   contains

     procedure :: name => gauss_seidel_name
     procedure :: solve

  end type gauss_seidel

contains

  pure function gauss_seidel_name(this) result(name)

    class(gauss_seidel), intent(in) :: this
    character(len=:), allocatable :: name

    associate (u1 => this); end associate

    name = 'gauss-seidel'

  end function gauss_seidel_name


  subroutine solve(this, rhs, x, achieved)

    class(gauss_seidel), intent(inout) :: this
    real(dp), intent(in)    :: rhs(:)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out)   :: achieved

    real(dp), allocatable :: d(:), y(:), r(:)
    integer , allocatable :: colours(:)
    real(dp) :: goal
    integer :: it, col, v

    call this % diagonal(d)
    do v = 1, size(d)
       if (abs(d(v)) < tiny(1.0_dp)) d(v) = huge(1.0_dp)
    end do

    call this % sweep_order(colours)

    goal = this % tolerance * (1.0_dp + this % norm(rhs))

    do it = 1, this % max_iterations

       do col = 1, maxval(colours)
          call this % matvec(x, y)
          r = rhs - y
          do v = 1, size(x)
             if (colours(v) == col) x(v) = x(v) + this % omega * r(v) / d(v)
          end do
       end do

       call this % matvec(x, y)
       achieved = this % norm(rhs - y)
       if (achieved < goal) return

    end do

  end subroutine solve

end module operation_gauss_seidel
