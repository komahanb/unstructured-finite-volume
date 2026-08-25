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
! The residual measured at the top of an iteration is the first
! colour's, so a sweep of c colours costs c products with the
! operator. Which colouring is swept is a question the type answers
! for itself - jacobi is this iteration over one colour class, and
! says so by extension.
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
     procedure :: colouring
     procedure :: solve

  end type gauss_seidel

contains

  pure function gauss_seidel_name(this) result(name)

    class(gauss_seidel), intent(in) :: this
    character(len=:), allocatable :: name

    associate (u1 => this); end associate

    name = 'gauss-seidel'

  end function gauss_seidel_name

  !===================================================================!
  ! The colour classes the sweep runs over: one entry per unknown,
  ! from the coupling's own colouring.
  !===================================================================!

  subroutine colouring(this, n, colours)

    class(gauss_seidel), intent(in)  :: this
    integer            , intent(in)  :: n
    integer, allocatable, intent(out) :: colours(:)

    call this % sweep_order(colours)

    if (size(colours) /= n) then
       error stop 'gauss_seidel: one colour per unknown'
    end if

  end subroutine colouring


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

    call this % colouring(size(x), colours)

    call this % begin_imbalance()

    do it = 1, this % max_iterations

       call this % matvec(x, y)
       r = rhs - y

       achieved = this % norm(r)
       if (this % halted(achieved, it)) return

       do col = 1, maxval(colours)
          if (col > 1) then
             call this % matvec(x, y)
             r = rhs - y
          end if
          do v = 1, size(x)
             if (colours(v) == col) x(v) = x(v) + this % omega * r(v) / d(v)
          end do
       end do

    end do

    call this % matvec(x, y)
    achieved = this % norm(rhs - y)

  end subroutine solve

end module operation_gauss_seidel
