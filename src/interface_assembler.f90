#include "scalar.fpp"

!=====================================================================!
! Module that contains class that assembles residual and performs
! matrix vector products on all elements in mesh
!
! Author: Komahan Boopathy
!=====================================================================!

module interface_assembler

  use iso_fortran_env, only : dp => REAL64
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

   contains  

     ! Deferred procedures
     procedure(add_residual_interface)               , deferred :: add_residual
     procedure(add_jacobian_vector_product_interface), deferred :: add_jacobian_vector_product
     procedure(add_initial_condition_interface)      , deferred :: add_initial_condition

     ! Provided procedures
     procedure :: get_num_state_vars, set_num_state_vars
     procedure :: get_description   , set_description
  
     ! Defined procedures
     procedure :: get_differential_order
     procedure :: set_differential_order

     ! Default finite difference jacobian implementation
     procedure :: add_jacobian => add_jacobian_fd

  end type assembler

  !===================================================================!
  ! Interfaces for deferred procedures
  !===================================================================!

  abstract interface

     !================================================================!
     ! Interface for residual assembly R(U,xi)
     !================================================================!

     impure subroutine add_residual_interface(this, residual, U, xi)

       import :: assembler

       class(assembler), intent(in)    :: this
       type(scalar)  , intent(inout) :: residual(:)
       type(scalar)  , intent(in)    :: U(:,:)
       type(scalar)  , intent(in)    :: xi(:)

     end subroutine add_residual_interface

     !================================================================!
     ! Routine to return the product of jacobian matrix with a compati-
     ! ble vector pdt <---- [scalar(i)*dR(U,X)/dU(i)]*vec
     !================================================================!
     
     impure subroutine add_jacobian_vector_product_interface(this, pdt, vec, &
          & scalars, U, xi)

       import :: assembler

       class(assembler) , intent(in)    :: this
       type(scalar)   , intent(inout) :: pdt(:)
       type(scalar)   , intent(in)    :: vec(:)
       type(scalar)   , intent(in)    :: scalars(:)
       type(scalar)   , intent(in)    :: U(:,:)
       type(scalar)   , intent(in)    :: xi(:)

     end subroutine add_jacobian_vector_product_interface

     !================================================================!
     ! Supplying the initial condition to march in time
     !================================================================!

     impure subroutine add_initial_condition_interface(this, U, xi)

       import :: assembler

       class(assembler), intent(in)    :: this
       type(scalar)  , intent(inout) :: U(:,:)
       type(scalar)  , intent(in)    :: xi(:)

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
  ! Jacobian assembly at each time step.
  !===================================================================!
  
  impure subroutine add_jacobian_fd(this, jacobian, coeff, U, xi)

    class(assembler) , intent(inout) :: this
    type(scalar)   , intent(inout) :: jacobian(:,:)
    type(scalar)   , intent(in)    :: coeff(:)
    type(scalar)   , intent(in)    :: U(:,:)
    type(scalar)   , intent(in)    :: xi(:)

    type(integer)             :: nvars, dorder
    type(integer)             :: n, m
    type(scalar), allocatable :: R(:), Rtmp(:), Utmp(:,:)
    real(8)     , parameter   :: dh = 1.0d-8

    nvars = this % get_num_state_vars()
    dorder = this % get_differential_order()

    ! Make a copy of state variables
    allocate(Utmp, source=U)

    ! Allocate space for residual perturbations
    allocate(R(nvars))
    allocate(Rtmp(nvars))

    ! Make a residual call with original variables
    R = 0.0d0
    call this % add_residual(R, U, xi)

    diff_order: do n = 1, dorder + 1

       loop_vars: do m = 1, nvars

          ! Perturb the m-th variable of n-th order
          Utmp(n,m) = Utmp(n,m) + dh

          ! Make a residual call with the perturbed variable
          rtmp = 0.0d0
          call this % add_residual(Rtmp, Utmp, xi)

          ! Unperturb (restore) the m-th variable of n-th order
          Utmp(n+1,m) = U(n+1,m)

          ! Approximate the jacobian with respect to the m-th variable
          jacobian(:,m) = jacobian(:,m) + coeff(n)*(Rtmp-R)/dh

       end do loop_vars

    end do diff_order

    ! Freeup memory
    deallocate(R, Rtmp)

  end subroutine add_jacobian_fd

end module interface_assembler
