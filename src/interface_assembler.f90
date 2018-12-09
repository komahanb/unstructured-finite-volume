#include "scalar.fpp"

!=====================================================================!
! Module that contains class that assembles residual and performs
! matrix vector products on all elements in mesh
!
! Author: Komahan Boopathy
!=====================================================================!

module interface_assembler

  use iso_fortran_env  , only : dp => REAL64
  implicit none
  
  private
  public :: assembler
 
  !===================================================================!
  ! Assembler for the physical system
  !===================================================================!
  
  type, abstract :: assembler

     type(character(len=:)), allocatable :: description

     type(integer) :: num_state_vars
     type(logical) :: approximate_jacobian
     type(integer) :: differential_order

     real(dp), allocatable :: S(:,:)

   contains  

     ! Deferred procedures
     procedure(add_residual_interface)               , deferred :: add_residual
     procedure(add_jacobian_vector_product_interface), deferred :: add_jacobian_vector_product
     procedure(add_initial_condition_interface)      , deferred :: add_initial_condition
     
     ! Assembler knows the size of state array
     procedure :: create_vector
     procedure :: create_state

     ! Provided procedures
     procedure :: get_num_state_vars, set_num_state_vars
     procedure :: get_description   , set_description
  
     ! Defined procedures
     procedure :: get_differential_order
     procedure :: set_differential_order

  end type assembler

  !===================================================================!
  ! Interfaces for deferred procedures
  !===================================================================!

  abstract interface

     !================================================================!
     ! Interface for residual assembly R(U,xi)
     !================================================================!

     impure subroutine add_residual_interface(this, residual, filter)

       import :: assembler

       class(assembler), intent(in)           :: this
       type(scalar)    , intent(inout)        :: residual(:)
       type(integer)   , intent(in), optional :: filter

     end subroutine add_residual_interface

     !================================================================!
     ! Routine to return the product of jacobian matrix with a compati-
     ! ble vector pdt <---- [scalar(i)*dR(U,X)/dU(i)]*vec
     !================================================================!
     
     impure subroutine add_jacobian_vector_product_interface(this, pdt, vec, scalars, filter)

       import :: assembler

       class(assembler) , intent(in)    :: this
       type(scalar)     , intent(inout) :: pdt(:)
       type(scalar)     , intent(in)    :: vec(:)
       type(scalar)     , intent(in)    :: scalars(:)
       type(integer)    , intent(in), optional :: filter

     end subroutine add_jacobian_vector_product_interface

     !================================================================!
     ! Supplying the initial condition to march in time
     !================================================================!

     impure subroutine add_initial_condition_interface(this, U)

       import :: assembler

       class(assembler), intent(in)    :: this
       type(scalar)    , intent(inout) :: U(:,:)

     end subroutine add_initial_condition_interface

  end interface

contains
  
  !===================================================================!
  ! Returns the number of state variables in the physical system
  !===================================================================!
  
  pure type(integer) function get_num_state_vars(this)

    class(assembler), intent(in) :: this

    get_num_state_vars = this % num_state_vars

  end function get_num_state_vars

  !===================================================================!
  ! Sets the number of state variables in the physical system
  !===================================================================!
  
  pure subroutine set_num_state_vars(this, num_state_vars)

    class(assembler), intent(inout) :: this
    type(integer)  , intent(in)   :: num_state_vars

    this % num_state_vars  = num_state_vars

  end subroutine set_num_state_vars
  
  !===================================================================!
  ! Returns the description set for the physical system
  !===================================================================!
  
  pure type(character) function get_description(this)

    class(assembler), intent(in) :: this

    get_description = this % description

  end function get_description

  !===================================================================!
  ! Sets the description for physical system
  !===================================================================!

  pure subroutine set_description(this, description)

    class(assembler), intent(inout) :: this
    type(character(len=*)), intent(in) :: description
    
    allocate(this % description, source = trim(description))

  end subroutine set_description

  !===================================================================!
  ! Returns the highest order of time derivative in the assembler
  !===================================================================!
  
  pure type(integer) function get_differential_order(this)

    class(assembler), intent(in) :: this

    get_differential_order = this % differential_order

  end function get_differential_order

  !===================================================================!
  ! Sets the highest order of time derivative in the assembler
  !===================================================================!

  pure subroutine set_differential_order(this, order)

    class(assembler), intent(inout) :: this
    type(integer) , intent(in)    :: order

    this % differential_order = order

  end subroutine set_differential_order

  !===================================================================!
  ! Create a state vector and sets values if a scalar is supplied
  !===================================================================!

  subroutine create_vector(this, x, val)
    
    class(assembler), intent(in)               :: this
    real(dp)        , intent(out), allocatable :: x(:)
    real(dp)        , intent(in) , optional    :: val
    
    if (allocated(x)) error stop "vector already allocated"
    allocate(x(this % num_state_vars))
    if (present(val))  x = val
    
  end subroutine create_vector

  !===================================================================!
  ! Create a state vector and sets values if a scalar is supplied
  !===================================================================!

  subroutine create_state(this, S, val)
    
    class(assembler), intent(in)               :: this
    real(dp)        , intent(out), allocatable :: S(:,:)
    real(dp)        , intent(in) , optional    :: val
    
    if (allocated(S)) error stop "vector already allocated"
    allocate( &
         & S( &
         & this % num_state_vars, &
         & this % get_differential_order() + 1 &
         & ))
    if (present(val))  S = val
    
  end subroutine create_state
  
end module interface_assembler
