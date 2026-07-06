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

  use iso_fortran_env    , only : dp => REAL64
  use module_solve_mode  , only : is_valid_mode
  use interface_marcher  , only : marcher
  use interface_assembler, only : assembler
  use interface_state    , only : state
  use class_differential_state, only : differential_state

  implicit none

  private
  public :: integrator

  !-------------------------------------------------------------------!
  ! Abstract integrator type
  !-------------------------------------------------------------------!

  type, abstract, extends(marcher) :: integrator

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

     ! the marcher contract: advance the state window along the step
     ! chain, recording the trajectory. integrate is the family wrapper
     ! (it seeds the window from the initial condition).
     procedure :: march
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
  ! The family wrapper: seed the state window from the system's initial
  ! condition and delegate to the march on the constructed system.
  !===================================================================!

  impure subroutine integrate(this, mode)

    class(integrator), intent(inout)        :: this
    integer          , intent(in), optional :: mode

    type(differential_state) :: s

    ! a wrong tag dies at the door with its name
    if (present(mode)) then
       if (.not. is_valid_mode(mode)) then
          write(*,'(1x,a,i0)') "integrate: invalid mode tag ", mode
          error stop "integrate: mode must be FORWARD or REVERSE"
       end if
    end if

    s = differential_state( &
         & this % system % get_num_state_vars(), &
         & this % system % get_differential_order())

    call this % system % add_initial_condition(s % U)

    call this % march(this % system, s, mode)

  end subroutine integrate

  !===================================================================!
  ! The marcher contract: the state window enters as the initial
  ! condition, is advanced one step at a time along the step chain
  ! (each step solving through the deferred step), and exits as the
  ! final window; the trajectory this % U is the record of the states
  ! visited - the chain's payload, not the state.
  !===================================================================!

  impure subroutine march(this, system, s, mode)

    class(integrator), intent(inout)        :: this
    class(assembler) , intent(inout)        :: system
    class(state)     , intent(inout)        :: s
    integer          , intent(in), optional :: mode

    integer :: k, p, ierr

    ! a wrong tag dies at the door with its name
    if (present(mode)) then
       if (.not. is_valid_mode(mode)) then
          write(*,'(1x,a,i0)') "integrator % march: invalid mode tag ", mode
          error stop "integrator % march: mode must be FORWARD or REVERSE"
       end if
    end if

    ! the deferred step advances the system this integrator was
    ! constructed with; a different system of the same size would march
    ! the wrong equations
    if (system % get_num_state_vars() .ne. this % system % get_num_state_vars()) then
       error stop "integrator % march: the system does not match the constructed one"
    end if

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

    ! The seed window enters the record
    select type (s)
    type is (differential_state)
       this % U(1,:,:) = s % U
    class default
       error stop "integrator % march: requires a differential_state"
    end select

    ! March one step at a time
    stepping: do k = 2, this % num_steps

       p = this % get_bandwidth(k)

       call this % step(this % time(k-p:k), this % U(k-p:k,:,:), this % h, p, ierr)

    end do stepping

    ! The state exits as the final window
    select type (s)
    type is (differential_state)
       s % U = this % U(this % num_steps,:,:)
    end select

  end subroutine march

end module interface_integrator
