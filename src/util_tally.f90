!=====================================================================!
! What a run consumes, attributed to where it was consumed.
!
! An amount is recorded under three indices: the level of the
! hierarchy the computation was inside, the derivative order being
! computed, and the kind of amount. That triple is the whole of the
! structure, and a matrix over any two of its axes is derived from it
! without any further data being stored.
!
! The level is a scope, not an argument. A caller descending the
! hierarchy opens a level and closes it again, and every amount
! recorded in between - including amounts recorded far below, inside a
! minimizer that has no record of the caller's levels - is recorded
! under it. That is how a module which stores nothing of the caller's
! hierarchy records into it: the module names the amount, and the
! scope names the level. The hierarchy itself is the caller's: its
! levels are declared, outermost first, when the tally is opened.
!
! Levels nest, so a level's elapsed time includes the time of the levels
! opened inside it.
!
! A tally is a value owned by the execution whose amounts it records:
! two executions record into two tallies, and a caller sums its
! executions' tallies into its own. No amount is stored in the module.
!
! Recording is off until open, and an amount recorded while off costs
! one test of a logical. Nothing here alters a computed value; an
! amount is observed and never fed back.
!
! An amount recorded with no level open is discarded rather than
! recorded under a level it did not occur in. An order outside the
! range open was given is discarded for the same reason.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module util_tally

  use iso_fortran_env, only : int64
  use util_precision  , only : dp

  implicit none

  private

  public :: tally
  public :: tally_num_events, tally_event_name, tally_event_of

  !-------------------------------------------------------------------!
  ! The levels an amount can be recorded under are the caller's,
  ! declared outermost first at open; the kinds of amount are named
  ! here so a caller records by name and not by a number the caller
  ! has to maintain.
  !-------------------------------------------------------------------!

  integer, parameter, public :: elapsed_time      = 1
  integer, parameter, public :: primal_loops   = 2
  integer, parameter, public :: tangent_loops  = 3
  integer, parameter, public :: adjoint_loops  = 4
  integer, parameter, public :: newton_solves  = 5
  integer, parameter, public :: linear_solves  = 6
  integer, parameter, public :: factorisations = 7

  integer, parameter :: num_events = 7

  character(len=14), parameter :: event_named(num_events) = &
       & ['elapsed_time  ', 'primal_loops  ', 'tangent_loops ', &
       &  'adjoint_loops ', 'newton_solves ', 'linear_solves ', &
       &  'factorisations']

  integer, parameter :: deepest = 32

  !-------------------------------------------------------------------!
  ! One tally per execution: the declared levels, the amounts by
  ! (level, order, event), the order now being computed and the stack
  ! of open levels with the instant each was opened at.
  !-------------------------------------------------------------------!

  type :: tally
     private
     integer :: levels_declared = 0
     character(len=:), allocatable :: level_named(:)
     logical  :: on = .false.
     real(dp), allocatable :: amounts(:,:,:)
     integer  :: highest = 0
     integer  :: order_now = 0
     integer  :: level_now(deepest) = 0
     real(dp) :: entered_at(deepest) = 0.0_dp
     integer  :: depth = 0
   contains
     procedure :: open
     procedure :: close
     procedure :: recording
     procedure :: num_levels
     procedure :: num_orders
     procedure :: level_name
     procedure :: order
     procedure :: enter
     procedure :: leave
     procedure :: record
     procedure :: amount
     procedure :: restarted
     procedure :: add
  end type tally

contains

  !===================================================================!
  ! Start recording, with storage for orders zero to highest_order, over
  ! the caller's levels, outermost first. A negative highest order, no
  ! levels, or more levels than the stack contains stops the program.
  !===================================================================!

  subroutine open(this, highest_order, levels)

    class(tally)    , intent(inout) :: this
    integer         , intent(in)    :: highest_order
    character(len=*), intent(in)    :: levels(:)

    if (highest_order < 0) then
       error stop 'util_tally: the highest order is zero or above'
    end if
    if (size(levels) < 1 .or. size(levels) > deepest) then
       error stop 'util_tally: the levels are one to the stack''s depth'
    end if

    this % highest         = highest_order
    this % levels_declared = size(levels)
    this % level_named     = levels

    if (allocated(this % amounts)) deallocate(this % amounts)
    allocate(this % amounts(this % levels_declared, 0:this % highest, num_events), source=0.0_dp)

    this % order_now = 0
    this % depth     = 0
    this % on        = .true.

  end subroutine open

  !===================================================================!
  ! Stop recording. The recorded amounts remain readable, so a caller
  ! closes before printing.
  !===================================================================!

  subroutine close(this)

    class(tally), intent(inout) :: this

    this % on = .false.

  end subroutine close

  pure logical function recording(this) result(on)

    class(tally), intent(in) :: this

    on = this % on

  end function recording

  pure integer function num_levels(this) result(n)

    class(tally), intent(in) :: this

    n = this % levels_declared

  end function num_levels

  pure integer function tally_num_events() result(n)

    n = num_events

  end function tally_num_events

  pure integer function num_orders(this) result(n)

    class(tally), intent(in) :: this

    n = this % highest

  end function num_orders

  pure function level_name(this, level) result(named)

    class(tally), intent(in) :: this
    integer     , intent(in) :: level
    character(len=:), allocatable :: named

    named = trim(this % level_named(level))

  end function level_name

  pure function tally_event_name(event) result(named)

    integer, intent(in) :: event
    character(len=:), allocatable :: named

    named = trim(event_named(event))

  end function tally_event_name

  !===================================================================!
  ! The event a name denotes, or zero where no event has that name.
  !===================================================================!

  pure integer function tally_event_of(named) result(event)

    character(len=*), intent(in) :: named

    integer :: i

    event = 0
    do i = 1, num_events
       if (trim(event_named(i)) == named) event = i
    end do

  end function tally_event_of

  !===================================================================!
  ! Which derivative order the amounts that follow belong to.
  !===================================================================!

  subroutine order(this, derivative_order)

    class(tally), intent(inout) :: this
    integer     , intent(in)    :: derivative_order

    if (.not. this % on) return

    this % order_now = derivative_order

  end subroutine order

  !===================================================================!
  ! Open a level. A depth beyond the stack, or a level that is not
  ! one of those declared, stops the program: recording under the
  ! wrong level is less useful than not recording.
  !===================================================================!

  subroutine enter(this, level)

    class(tally), intent(inout) :: this
    integer     , intent(in)    :: level

    if (.not. this % on) return

    if (level < 1 .or. level > this % levels_declared) then
       error stop 'util_tally: a level is one of those declared'
    end if
    if (this % depth == deepest) then
       error stop 'util_tally: the levels opened are within the stack'
    end if

    this % depth                    = this % depth + 1
    this % level_now(this % depth)  = level
    this % entered_at(this % depth) = clock()

  end subroutine enter

  !===================================================================!
  ! Close the innermost level, recording the time it was open. Closing
  ! one that was never opened stops the program.
  !===================================================================!

  subroutine leave(this)

    class(tally), intent(inout) :: this

    real(dp) :: elapsed

    if (.not. this % on) return

    if (this % depth == 0) then
       error stop 'util_tally: a level closed was opened'
    end if

    elapsed = clock() - this % entered_at(this % depth)
    if (this % order_now >= 0 .and. this % order_now <= this % highest) then
       this % amounts(this % level_now(this % depth), this % order_now, elapsed_time) = &
            & this % amounts(this % level_now(this % depth), this % order_now, elapsed_time) + elapsed
    end if

    this % depth = this % depth - 1

  end subroutine leave

  !===================================================================!
  ! Record one event under the level now open.
  !===================================================================!

  subroutine record(this, event)

    class(tally), intent(inout) :: this
    integer     , intent(in)    :: event

    if (.not. this % on) return
    if (this % depth == 0) return
    if (this % order_now < 0 .or. this % order_now > this % highest) return

    if (event < 1 .or. event > num_events) then
       error stop 'util_tally: an event is one of the kinds named'
    end if

    this % amounts(this % level_now(this % depth), this % order_now, event) = &
         & this % amounts(this % level_now(this % depth), this % order_now, event) + 1.0_dp

  end subroutine record

  !===================================================================!
  ! The recorded amount. An index outside the extent open allocated
  ! stops the program.
  !===================================================================!

  pure real(dp) function amount(this, level, order, event) result(elapsed)

    class(tally), intent(in) :: this
    integer     , intent(in) :: level, order, event

    if (.not. allocated(this % amounts)) then
       error stop 'util_tally: an amount is read after open'
    end if
    if (level < 1 .or. level > this % levels_declared) then
       error stop 'util_tally: a level is one of those declared'
    end if
    if (order < 0 .or. order > this % highest) then
       error stop 'util_tally: an order is within the range opened'
    end if
    if (event < 1 .or. event > num_events) then
       error stop 'util_tally: an event is one of the kinds named'
    end if

    elapsed = this % amounts(level, order, event)

  end function amount

  !===================================================================!
  ! A tally over the same levels and orders, recording if this one
  ! is, with zero amounts and no level open: the tally an execution
  ! begins with, whose amounts its caller adds into this one.
  !===================================================================!

  function restarted(this) result(initial)

    class(tally), intent(in) :: this
    type(tally) :: initial

    initial = this
    initial % order_now  = 0
    initial % level_now  = 0
    initial % entered_at = 0.0_dp
    initial % depth      = 0
    if (allocated(initial % amounts)) initial % amounts = 0.0_dp

  end function restarted

  !===================================================================!
  ! Add another tally's amounts into this one, cell by cell. A tally
  ! that was never opened adds nothing; one opened over other levels
  ! or orders stops the program, since its cells name other scopes.
  !===================================================================!

  subroutine add(this, other)

    class(tally), intent(inout) :: this
    type(tally) , intent(in)    :: other

    if (.not. allocated(other % amounts)) return
    if (.not. allocated(this % amounts)) then
       error stop 'util_tally: amounts are added into an opened tally'
    end if
    if (any(shape(other % amounts) /= shape(this % amounts))) then
       error stop 'util_tally: amounts added are over the same levels and orders'
    end if

    this % amounts = this % amounts + other % amounts

  end subroutine add

  !===================================================================!
  ! Seconds, from the system clock.
  !===================================================================!

  real(dp) function clock() result(s)

    integer(int64) :: ticks, rate

    call system_clock(ticks, rate)
    s = real(ticks, dp) / real(rate, dp)

  end function clock

end module util_tally
