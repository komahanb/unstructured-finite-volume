#include "scalar.fpp"

!=====================================================================!
! Class that implements two dimensional unsteady heat equation's
! residual and jacobians
!
! Author: Komahan Boopathy
!=====================================================================!

module class_physics_unsteady_heat

  use interface_physics, only : physics

  implicit none

  private
  public :: unsteady_heat

  !-------------------------------------------------------------------!
  ! Type that implements first order transport equations
  !-------------------------------------------------------------------!
  
  type, extends(physics) :: unsteady_heat

     type(scalar) :: diffusivity 
     
   contains
     
     ! Implement deferred procedures from superclasses
     procedure :: add_residual
     procedure :: add_jacobian_vector_product => add_jacobian_vector_product
     procedure :: add_initial_condition
     
     ! Destructor
     final :: destruct
     
  end type unsteady_heat
  
  ! Interface to construct the physical system
  interface unsteady_heat
     procedure construct
  end interface unsteady_heat
  
contains
  
  !===================================================================!
  ! Constructor for unsteady heat physics
  !===================================================================!
  
  type(unsteady_heat) function construct(diffusivity) result (this)

    type(scalar), intent(in) :: diffusivity

    call this % set_description('Unsteady Heat')

    ! Set time order of physical system
    call this % set_differential_order(1)

    ! System parameters
    this % diffusivity = diffusivity

    ! Use FD jacobian or supplied jacobian
    this % approximate_jacobian = .false.

    ! Set the number of state variables per node. Temperature is the
    ! only unknown
    call this % set_num_state_vars(1)

  end function construct
  
  !===================================================================!
  ! Destructor for unsteady heat physics
  !===================================================================!
  
  pure subroutine destruct(this)

    type(unsteady_heat), intent(inout) :: this

    deallocate(this % description)

  end subroutine destruct
    
  !================================================================!
  ! Implements Interface for residual assembly R(U,xi)
  !================================================================!
  
  impure subroutine add_residual(this, residual, U, xi)

    class(unsteady_heat), intent(in)    :: this
    type(scalar)        , intent(inout) :: residual(:)
    type(scalar)        , intent(in)    :: U(:,:)
    type(scalar)        , intent(in)    :: xi(:)

  end subroutine add_residual

  !================================================================!
  ! Routine to return the product of jacobian matrix with a compati-
  ! ble vector pdt <---- [scalar(i)*dR(U,X)/dU(i)]*vec
  !================================================================!

  impure subroutine add_jacobian_vector_product(this, pdt, vec, &
       & scalars, U, xi)

    class(unsteady_heat) , intent(in)    :: this
    type(scalar)         , intent(inout) :: pdt(:)
    type(scalar)         , intent(in)    :: vec(:)
    type(scalar)         , intent(in)    :: scalars(:)
    type(scalar)         , intent(in)    :: U(:,:)
    type(scalar)         , intent(in)    :: xi(:)

  end subroutine add_jacobian_vector_product
  
  !================================================================!
  ! Supplying the initial condition to march in time
  !================================================================!
  
  impure subroutine add_initial_condition(this, U, xi)

    class(unsteady_heat) , intent(in)    :: this
    type(scalar)         , intent(inout) :: U(:,:)
    type(scalar)         , intent(in)    :: xi(:)

  end subroutine add_initial_condition

end module class_physics_unsteady_heat
