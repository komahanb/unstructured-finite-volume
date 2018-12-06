#include "scalar.fpp"

!=====================================================================!
! Module that contains common procedures for any physical system
! subject to governing equations
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module physics_interface

  implicit none
  
  private

  public :: physics
 
  !-------------------------------------------------------------------!
  ! Type that models any physical system
  !-------------------------------------------------------------------!
  
  type, abstract :: physics

     type(character(len=:)), allocatable :: description
     type(integer)                       :: num_state_vars
     type(logical)                       :: approximate_jacobian
     logical                             :: sparse ! temp hack
     type(scalar)          , allocatable :: X(:,:) ! temp maybe? contain ref to a class?

   contains  

     ! Provided procedures
     procedure :: get_num_state_vars      , set_num_state_vars
     procedure :: get_description         , set_description

     ! Dummy implementation. Interface it later
     procedure :: add_jacobian_vector_product

     ! Deferred procedures
     procedure(add_residual_interface), deferred :: add_residual
     procedure(add_jacobian_interface), deferred :: add_jacobian
     !procedure(add_jacobian_vector_product_interface), deferred :: add_jacobian_vector_product

  end type physics

  ! Interfaces to deferred procedures
  abstract interface

     !----------------------------------------------------------------!
     ! Interface for residual assembly
     !----------------------------------------------------------------!

     pure subroutine add_residual_interface(this, residual, U, X)
       
       import :: physics
       
       class(physics) , intent(inout) :: this
       type(scalar)   , intent(inout) :: residual(:)
       type(scalar)   , intent(in)    :: U(:,:)
       type(scalar)   , intent(in)    :: X(:,:)

     end subroutine add_residual_interface

     !----------------------------------------------------------------!
     ! Interface for jacobian assembly
     !----------------------------------------------------------------!
     
     pure subroutine add_jacobian_interface(this, jacobian, coeff, U, X)

       import :: physics

       class(physics) , intent(inout) :: this
       type(scalar)   , intent(inout) :: jacobian(:,:)
       type(scalar)   , intent(in)    :: coeff(:)
       type(scalar)   , intent(in)    :: U(:,:)
       type(scalar)   , intent(in)    :: X(:,:)
       
     end subroutine add_jacobian_interface

  end interface

contains
  
  !=====================================================================!
  ! Routine to return the product of jacobian matrix with a compatible
  ! vector
  ! =====================================================================!

  pure subroutine add_jacobian_vector_product(this, pdt, vec, scalars, U, X)

    class(physics) , intent(inout) :: this
    type(scalar)   , intent(inout) :: pdt(:)
    !type(scalar)   , intent(inout) :: jacobian(:,:)

    type(scalar)   , intent(in)    :: vec(:)
    type(scalar)   , intent(in)    :: scalars(:)
    type(scalar)   , intent(in)    :: U(:,:)
    type(scalar)   , intent(in)    :: X(:,:)

  end subroutine add_jacobian_vector_product
 
  !===================================================================!
  ! Returns the number of state variables in the physical system
  !===================================================================!
  
  pure type(integer) function get_num_state_vars(this)

    class(physics), intent(in) :: this

    get_num_state_vars = this % num_state_vars

  end function get_num_state_vars

  !===================================================================!
  ! Sets the number of state variables in the physical system
  !===================================================================!
  
  pure subroutine set_num_state_vars(this, num_state_vars)

    class(physics), intent(inout) :: this
    type(integer)  , intent(in)    :: num_state_vars

    this % num_state_vars  = num_state_vars

  end subroutine set_num_state_vars
  
  !===================================================================!
  ! Returns the description set for the physical system
  !===================================================================!
  
  pure type(character) function get_description(this)

    class(physics), intent(in) :: this

    get_description = this % description

  end function get_description

  !===================================================================!
  ! Sets the description for physical system
  !===================================================================!

  pure subroutine set_description(this, description)

    class(physics), intent(inout) :: this
    type(character(len=*)), intent(in)    :: description
    
    allocate(this % description, source = trim(description))

  end subroutine set_description

end module physics_interface
