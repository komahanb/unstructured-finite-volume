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
! Recording is off until tally_open, and an amount recorded while off
! costs one test of a logical. Nothing here alters a computed value;
! an amount is observed and never fed back.
!
! An amount recorded with no level open is discarded rather than
! recorded under a level it did not occur in. An order outside the
! range tally_open was given is discarded for the same reason.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module util_tally

  use iso_fortran_env, only : int64
  use util_precision  , only : dp

  implicit none

  private

  public :: tally_open, tally_close, tally_recording
  public :: tally_order, tally_enter, tally_leave
  public :: tally_record, tally_amount
  public :: tally_num_levels, tally_num_events, tally_num_orders
  public :: tally_level_name, tally_event_name, tally_event_of

  !-------------------------------------------------------------------!
  ! The levels an amount can be recorded under are the caller's,
  ! declared outermost first at tally_open; the kinds of amount are
  ! named here so a caller records by name and not by a number the
  ! caller has to maintain.
  !-------------------------------------------------------------------!

  integer, save :: num_levels = 0
  character(len=:), allocatable, save :: level_named(:)

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

  !-------------------------------------------------------------------!
  ! One tally per run. The tally observes and owns nothing the
  ! computation reads, so a single tally serves every caller.
  !-------------------------------------------------------------------!

  integer, parameter :: deepest = 32

  logical , save :: recording = .false.
  real(dp), save, allocatable :: amount(:,:,:)
  integer , save :: highest = 0
  integer , save :: order_now = 0
  integer , save :: level_now(deepest) = 0
  real(dp), save :: entered_at(deepest) = 0.0_dp
  integer , save :: depth = 0

contains

  !===================================================================!
  ! Start recording, with storage for orders zero to highest_order, over
  ! the caller's levels, outermost first. A negative highest order, no
  ! levels, or more levels than the stack contains stops the program.
  !===================================================================!

  subroutine tally_open(highest_order, levels)

    integer         , intent(in) :: highest_order
    character(len=*), intent(in) :: levels(:)

    if (highest_order < 0) then
       error stop 'util_tally: the highest order is zero or above'
    end if
    if (size(levels) < 1 .or. size(levels) > deepest) then
       error stop 'util_tally: the levels are one to the stack''s depth'
    end if

    highest     = highest_order
    num_levels  = size(levels)
    level_named = levels

    if (allocated(amount)) deallocate(amount)
    allocate(amount(num_levels, 0:highest, num_events), source=0.0_dp)

    order_now = 0
    depth     = 0
    recording = .true.

  end subroutine tally_open

  !===================================================================!
  ! Stop recording. The recorded amounts remain readable, so a caller
  ! closes before printing.
  !===================================================================!

  subroutine tally_close()

    recording = .false.

  end subroutine tally_close

  pure logical function tally_recording() result(on)

    on = recording

  end function tally_recording

  pure integer function tally_num_levels() result(n)

    n = num_levels

  end function tally_num_levels

  pure integer function tally_num_events() result(n)

    n = num_events

  end function tally_num_events

  pure integer function tally_num_orders() result(n)

    n = highest

  end function tally_num_orders

  pure function tally_level_name(level) result(named)

    integer, intent(in) :: level
    character(len=:), allocatable :: named

    named = trim(level_named(level))

  end function tally_level_name

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

  subroutine tally_order(order)

    integer, intent(in) :: order

    if (.not. recording) return

    order_now = order

  end subroutine tally_order

  !===================================================================!
  ! Open a level. A depth beyond the stack, or a level that is not
  ! one of those declared, stops the program: recording under the
  ! wrong level is less useful than not recording.
  !===================================================================!

  subroutine tally_enter(level)

    integer, intent(in) :: level

    if (.not. recording) return

    if (level < 1 .or. level > num_levels) then
       error stop 'util_tally: a level is one of those declared'
    end if
    if (depth == deepest) then
       error stop 'util_tally: the levels opened are within the stack'
    end if

    depth             = depth + 1
    level_now(depth)  = level
    entered_at(depth) = clock()

  end subroutine tally_enter

  !===================================================================!
  ! Close the innermost level, recording the time it was open. Closing
  ! one that was never opened stops the program.
  !===================================================================!

  subroutine tally_leave()

    real(dp) :: elapsed

    if (.not. recording) return

    if (depth == 0) then
       error stop 'util_tally: a level closed was opened'
    end if

    elapsed = clock() - entered_at(depth)
    if (order_now >= 0 .and. order_now <= highest) then
       amount(level_now(depth), order_now, elapsed_time) = &
            & amount(level_now(depth), order_now, elapsed_time) + elapsed
    end if

    depth = depth - 1

  end subroutine tally_leave

  !===================================================================!
  ! Record one event under the level now open.
  !===================================================================!

  subroutine tally_record(event)

    integer, intent(in) :: event

    if (.not. recording) return
    if (depth == 0) return
    if (order_now < 0 .or. order_now > highest) return

    if (event < 1 .or. event > num_events) then
       error stop 'util_tally: an event is one of the kinds named'
    end if

    amount(level_now(depth), order_now, event) = &
         & amount(level_now(depth), order_now, event) + 1.0_dp

  end subroutine tally_record

  !===================================================================!
  ! The recorded amount. An index outside the extent tally_open
  ! allocated stops the program.
  !===================================================================!

  pure real(dp) function tally_amount(level, order, event) result(elapsed)

    integer, intent(in) :: level, order, event

    if (.not. allocated(amount)) then
       error stop 'util_tally: an amount is read after tally_open'
    end if
    if (level < 1 .or. level > num_levels) then
       error stop 'util_tally: a level is one of those declared'
    end if
    if (order < 0 .or. order > highest) then
       error stop 'util_tally: an order is within the range opened'
    end if
    if (event < 1 .or. event > num_events) then
       error stop 'util_tally: an event is one of the kinds named'
    end if

    elapsed = amount(level, order, event)

  end function tally_amount

  !===================================================================!
  ! Seconds, from the system clock.
  !===================================================================!

  real(dp) function clock() result(s)

    integer(int64) :: ticks, rate

    call system_clock(ticks, rate)
    s = real(ticks, dp) / real(rate, dp)

  end function clock

end module util_tally
