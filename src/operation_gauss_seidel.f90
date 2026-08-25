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

  use util_precision  , only : dp
  use operation_minimization, only : minimizer
  use util_tally, only : tally_record, linear_solves

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
        integer :: it, col, v

    call tally_record(linear_solves)

    call this % diagonal(d)
    do v = 1, size(d)
       if (abs(d(v)) < tiny(1.0_dp)) d(v) = huge(1.0_dp)
    end do

    call this % sweep_order(colours)

    call this % begin_imbalance()

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
       call this % note_imbalance(achieved)
       if (this % converged(achieved)) return
       if (this % exhausted(it)) return

    end do

  end subroutine solve

end module operation_gauss_seidel
