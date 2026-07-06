#include "scalar.fpp"

!=====================================================================!
! Abstract state: the object a marcher advances. interface_field draws
! the picture - discretising the field turns it into its state, the
! finite vector of coefficients U - and this interface is that state,
! as a type.
!
! The contract has two parts. update is where the families genuinely
! differ: each concrete state knows how a correction moves it (a plain
! vector adds the correction; a state carrying time derivatives applies
! one correction to every derivative order at once, weighted by the
! linearization coefficients). values presents the state to the system
! for residual evaluation, as the (nvars, order+1) array matching the
! assembler's state seat S.
!
! The caller owns the state object: a march receives it, advances it,
! and the answer is returned where it was given - no side channels.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module interface_state

  use iso_fortran_env, only : dp => REAL64

  implicit none

  private
  public :: state

  type, abstract :: state

   contains

     ! this <- this + correction; the concrete knows its rule
     procedure(update_interface), deferred :: update

     ! the dof array (nvars, order+1) the system evaluates
     procedure(values_interface), deferred :: values

  end type state

  abstract interface

     pure subroutine update_interface(this, correction)
       import :: state, dp
       class(state), intent(inout) :: this
       real(dp)    , intent(in)    :: correction(:)
     end subroutine update_interface

     pure function values_interface(this) result(v)
       import :: state
       class(state), intent(in) :: this
       type(scalar), allocatable :: v(:,:)
     end function values_interface

  end interface

contains

end module interface_state
