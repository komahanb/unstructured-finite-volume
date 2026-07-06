#include "scalar.fpp"

!=====================================================================!
! The concrete state: the solution and its time derivatives up to the
! differential order, U(nvars, order+1), with the linearization
! coefficients coeff(order+1) that define how a correction moves it:
!
!     U(:,k) <- U(:,k) + coeff(k) * correction
!
! the condensed newton update, owned by the state. Order 0 with
! coefficient [1] is the plain solution vector u <- u + du - the linear
! solver's state is the degenerate instance of this one class.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_differential_state

  use iso_fortran_env, only : dp => REAL64
  use interface_state, only : state

  implicit none

  private
  public :: differential_state

  type, extends(state) :: differential_state

     type(scalar), allocatable :: U(:,:)      ! (nvars, order+1)
     type(scalar), allocatable :: coeff(:)    ! (order+1) update weights

   contains

     procedure :: update
     procedure :: values

  end type differential_state

  interface differential_state
     module procedure create
  end interface differential_state

contains

  !===================================================================!
  ! A zero state of the given size and differential order. The default
  ! coefficients are [1, 0, ...] (a correction moves the solution
  ! only); a linearization supplies its own, one weight per order.
  !===================================================================!

  pure type(differential_state) function create(nvars, order, coeff) result(this)

    integer     , intent(in)           :: nvars
    integer     , intent(in)           :: order
    type(scalar), intent(in), optional :: coeff(:)

    allocate(this % U(nvars, order + 1))
    this % U = 0.0_dp

    allocate(this % coeff(order + 1))
    this % coeff    = 0.0_dp
    this % coeff(1) = 1.0_dp

    if (present(coeff)) then
       if (size(coeff) .ne. order + 1) then
          error stop "differential_state: coeff must have one weight per order"
       end if
       this % coeff = coeff
    end if

  end function create

  !===================================================================!
  ! Apply one correction to every derivative order, weighted by the
  ! linearization coefficients. At order 0 this is u <- u + du.
  !===================================================================!

  pure subroutine update(this, correction)

    class(differential_state), intent(inout) :: this
    real(dp)                 , intent(in)    :: correction(:)

    integer :: k

    do k = 1, size(this % coeff)
       this % U(:,k) = this % U(:,k) + this % coeff(k)*correction
    end do

  end subroutine update

  !===================================================================!
  ! The dof array the system evaluates
  !===================================================================!

  pure function values(this) result(v)

    class(differential_state), intent(in) :: this
    type(scalar), allocatable :: v(:,:)

    v = this % U

  end function values

end module class_differential_state
