#include "scalar.fpp"

!=====================================================================!
! Module that contains common procedures for any time-dependent
! physical system subject to governing equations of n-th order
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module dynamic_physics_interface

  use physics_interface, only : physics

  implicit none
  
  private
  
  public :: dynamics
  
  type, abstract, extends(physics) :: dynamics

     type(integer) :: differential_order ! order of the differential equation
     
   contains
     
     ! Defined procedures
     procedure :: get_differential_order
     procedure :: set_differential_order

     ! Default finite difference jacobian implementation
     procedure :: add_jacobian => add_jacobian_fd

     ! Deferred procedure to subtypes
     procedure(initial_condition_interface), deferred :: get_initial_condition

  end type dynamics

  ! Interfaces to deferred procedures
  interface

     !----------------------------------------------------------------!
     ! Supplying the initial condition to march in time
     !----------------------------------------------------------------!
     
     pure subroutine initial_condition_interface(this, U, X)

       import :: dynamics

       class(dynamics), intent(in)  :: this
       type(scalar)   , intent(inout) :: U(:,:)
       type(scalar)   , intent(in)    :: X(:,:)

     end subroutine initial_condition_interface

  end interface

contains
  
  !===================================================================!
  ! Returns the highest order of time derivative in the physics
  !===================================================================!
  
  pure type(integer) function get_differential_order(this)

    class(dynamics), intent(in) :: this

    get_differential_order = this % differential_order

  end function get_differential_order

  !===================================================================!
  ! Sets the highest order of time derivative in the physics
  !===================================================================!
  
  pure subroutine set_differential_order(this, order)

    class(dynamics), intent(inout) :: this
    type(integer), intent(in)      :: order

    this % differential_order = order

  end subroutine set_differential_order
  
  !===================================================================!
  ! Jacobian assembly at each time step.
  !===================================================================!
  
  pure subroutine add_jacobian_fd(this, jacobian, coeff, U, X)

    class(dynamics) , intent(inout) :: this
    type(scalar)    , intent(inout) :: jacobian(:,:)
    type(scalar)    , intent(in)    :: coeff(:)
    type(scalar)    , intent(in)    :: U(:,:)
    type(scalar)    , intent(in)    :: X(:,:)

    type(integer)             :: nvars, dorder
    type(integer)             :: n, m    
    type(scalar), allocatable :: R(:), Rtmp(:), Utmp(:,:)
    real(8)     , parameter   ::  dh = 1.0d-8

    nvars = this % get_num_state_vars()
    dorder = this % get_differential_order()

    ! Make a copy of state variables
    allocate(Utmp, source=U)

    ! Allocate space for residual perturbations
    allocate(R(nvars))
    allocate(Rtmp(nvars))

    ! Make a residual call with original variables
    R = 0.0d0
    call this % add_residual(R, U, X)

    diff_order: do n = 1, dorder + 1

       loop_vars: do m = 1, nvars

          ! Perturb the m-th variable of n-th order
          Utmp(n,m) = Utmp(n,m) + dh
          
          ! Make a residual call with the perturbed variable
          rtmp = 0.0d0
          call this % add_residual(Rtmp, Utmp, X)

          ! Unperturb (restore) the m-th variable of n-th order
          Utmp(n+1,m) = U(n+1,m)

          ! Approximate the jacobian with respect to the m-th variable
          jacobian(:,m) = jacobian(:,m) + coeff(n)*(Rtmp-R)/dh

       end do loop_vars

    end do diff_order

    ! Freeup memory
    deallocate(R, Rtmp)

  end subroutine add_jacobian_fd
  
end module dynamic_physics_interface
