#include "scalar.fpp"

!=====================================================================!
! Abstract time integrator. Marches the semi-discrete system
! R(t, U) = 0 supplied by a class(assembler), where the per-step state
! U = [q, qdot, qddot, ...] is stored as (num_state_vars, order+1).
! Concrete schemes (bdf, ...) implement `step` (advance one step,
! solving the nonlinear system when implicit) and `get_bandwidth` (how
! many past steps the multistep scheme needs).
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module interface_integrator

  use iso_fortran_env   , only : dp => REAL64
  use interface_assembler, only : assembler

  implicit none

  private
  public :: integrator

  !-------------------------------------------------------------------!
  ! Abstract integrator type
  !-------------------------------------------------------------------!

  type, abstract :: integrator

     class(assembler), allocatable :: system     ! the semi-discrete system
     real(dp)        , allocatable :: time(:)     ! time at each step
     type(scalar)    , allocatable :: U(:,:,:)    ! primal trajectory (step, nvars, order+1)
     type(scalar)    , allocatable :: psi(:,:)    ! adjoint trajectory (step, nvars)

     real(dp) :: tinit
     real(dp) :: tfinal
     real(dp) :: h
     logical  :: implicit
     integer  :: num_steps

   contains

     procedure :: construct
     procedure :: destruct
     procedure :: integrate

     ! Deferred to concrete schemes
     procedure(step_interface)            , deferred :: step
     procedure(get_bandwidth_interface)   , deferred :: get_bandwidth
     procedure(get_stencil_coeff_interface), deferred :: get_stencil_coeff

  end type integrator

  !-------------------------------------------------------------------!
  ! Interfaces to the deferred procedures
  !-------------------------------------------------------------------!

  abstract interface

     !================================================================!
     ! Advance one step: given the window of the last p+1 times/states,
     ! fill the newest state (solving R = 0 when implicit).
     !================================================================!

     impure subroutine step_interface(this, t, U, h, p, ierr)

       import :: integrator, dp

       class(integrator), intent(inout) :: this
       real(dp)         , intent(inout) :: t(:)
       type(scalar)     , intent(inout) :: U(:,:,:)
       integer          , intent(in)    :: p
       real(dp)         , intent(in)    :: h
       integer          , intent(out)   :: ierr

     end subroutine step_interface

     !================================================================!
     ! Number of past steps the scheme uses at the given step index
     !================================================================!

     pure integer function get_bandwidth_interface(this, step_index) result(width)

       import :: integrator

       class(integrator), intent(in) :: this
       integer          , intent(in) :: step_index

     end function get_bandwidth_interface

     !================================================================!
     ! The first-derivative stencil at bandwidth p: scoeff(i+1), i=0..p,
     ! is the coefficient multiplying U(k-i) when forming udot_k (for bdf,
     ! A(p,i+1)/h). The adjoint needs the whole stencil for the backward
     ! sweep's future-step coupling, not just the leading coefficient.
     !================================================================!

     impure subroutine get_stencil_coeff_interface(this, p, h, scoeff)

       import :: integrator, dp

       class(integrator), intent(in)               :: this
       integer          , intent(in)               :: p
       real(dp)         , intent(in)               :: h
       type(scalar)     , intent(out), allocatable :: scoeff(:)

     end subroutine get_stencil_coeff_interface

  end interface

contains

  !===================================================================!
  ! Base class constructor
  !===================================================================!

  pure subroutine construct(this, system, tinit, tfinal, h, implicit)

    class(integrator), intent(inout) :: this
    class(assembler) , intent(in)    :: system
    real(dp)         , intent(in)    :: tinit, tfinal, h
    logical          , intent(in)    :: implicit

    allocate(this % system, source = system)

    this % tinit     = tinit
    this % tfinal    = tfinal
    this % h         = h
    this % implicit  = implicit
    this % num_steps = floor((tfinal - tinit)/h) + 1

  end subroutine construct

  !===================================================================!
  ! Base class destructor
  !===================================================================!

  pure subroutine destruct(this)

    class(integrator), intent(inout) :: this

    if (allocated(this % system)) deallocate(this % system)
    if (allocated(this % time))   deallocate(this % time)
    if (allocated(this % U))      deallocate(this % U)
    if (allocated(this % psi))    deallocate(this % psi)

  end subroutine destruct

  !===================================================================!
  ! Drive the integration. FORWARD (default) marches the primal system
  ! from tinit to tfinal over the trajectory U; REVERSE drives the adjoint
  ! backward sweep over psi (wired by the concrete integrator).
  !===================================================================!

  impure subroutine integrate(this, mode)

    class(integrator), intent(inout)        :: this
    integer          , intent(in), optional :: mode

    integer :: k, p, ierr

    ! Re-entrant: a fresh march each call (the adjoint re-solves at
    ! perturbed designs), so drop any previous history first.
    if (allocated(this % time)) deallocate(this % time)
    if (allocated(this % U))    deallocate(this % U)

    ! Time history
    allocate(this % time(this % num_steps))
    this % time    = 0.0_dp
    this % time(1) = this % tinit

    ! State history: (step, nvars, order+1) to match the assembler state
    allocate(this % U( &
         & this % num_steps, &
         & this % system % get_num_state_vars(), &
         & this % system % get_differential_order() + 1))
    this % U = 0.0d0

    ! Initial condition into the first step
    call this % system % add_initial_condition(this % U(1,:,:))

    ! March one step at a time
    march: do k = 2, this % num_steps

       p = this % get_bandwidth(k)

       call this % step(this % time(k-p:k), this % U(k-p:k,:,:), this % h, p, ierr)

    end do march

  end subroutine integrate

end module interface_integrator
