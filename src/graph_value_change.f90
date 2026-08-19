!=====================================================================!
! THE VALUE CHANGE
!
! The pure attached-value member of the reversible change family:
! one concrete change that updates a value map through the one
! generic change_controller -
!
!      apply    save the seat's prior state, then vouch the new
!               values (attaching first when no seat exists)
!      check    the caller-provided verdict
!      keep     leave the map as updated
!      revert   restore the seat EXACTLY: unattached seats leave,
!               unknown seats stand untrusted again, known seats
!               get their old values back
!
! This is not a controller and chooses no policy: the change owns
! the rollback details, the map owns identity, status, and value
! storage, and the protocol owns the lifecycle. It touches values
! and never structure - the mirror image of a structural change
! under the same four verbs.
!
! Bound state is copied at bind - the element with its identity,
! the new values whole - so a caller mutating its array after
! binding changes nothing; only the map is held by pointer,
! because the map is the thing being changed.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module graph_value_change

  use iso_fortran_env      , only : dp => REAL64
  use fractal_graph        , only : graph
  use graph_value_map      , only : value_map, &
       & VALUE_UNATTACHED, VALUE_UNKNOWN, VALUE_KNOWN
  use graph_change_protocol, only : reversible_change, change_result

  implicit none

  private
  public :: value_change

  !===================================================================!
  ! One bound update: WHICH map (by pointer - the map is the thing
  ! being changed), WHICH graph (copied, identity and all), WHAT
  ! values (copied at bind), and the saved prior state revert
  ! restores exactly.
  !===================================================================!

  type, extends(reversible_change) :: value_change

     type(value_map), pointer :: values => null()
     type(graph)              :: element
     real(dp), allocatable    :: new_values(:)

     logical :: check_passes = .true.

     integer :: old_status = VALUE_UNATTACHED
     real(dp), allocatable :: old_values(:)

   contains

     procedure :: bind
     procedure :: apply
     procedure :: check
     procedure :: keep
     procedure :: revert

  end type value_change

contains

  !===================================================================!
  ! Bind one update: point at the map, copy the element and the
  ! new values, take the check verdict (true unless the caller
  ! says otherwise), and clear any saved state from an earlier
  ! life.
  !===================================================================!

  subroutine bind(this, values, element, new_values, check_passes)

    class(value_change)     , intent(inout) :: this
    type(value_map) , target, intent(inout) :: values
    type(graph)             , intent(in)    :: element
    real(dp)                , intent(in)    :: new_values(:)
    logical, optional       , intent(in)    :: check_passes

    this % values     => values
    this % element    = element
    this % new_values = new_values

    this % check_passes = .true.
    if (present(check_passes)) this % check_passes = check_passes

    this % old_status = VALUE_UNATTACHED
    if (allocated(this % old_values)) deallocate(this % old_values)

  end subroutine bind

  !===================================================================!
  ! Apply: remember what the seat was, then vouch the new values -
  ! attaching first when no seat exists. The identity gate and the
  ! known-value laws stay the map's own; this change adds only the
  ! memory revert needs.
  !===================================================================!

  subroutine apply(this, result)

    class(value_change), intent(inout) :: this
    type(change_result), intent(inout) :: result

    if (.not. associated(this % values)) then
       error stop 'value_change: value map is bound'
    end if

    this % old_status = this % values % status_of(this % element)

    if (allocated(this % old_values)) deallocate(this % old_values)
    if (this % old_status == VALUE_KNOWN) then
       call this % values % value_of(this % element, this % old_values)
    end if

    if (this % old_status == VALUE_UNATTACHED) then
       call this % values % attach_unknown(this % element)
    end if

    call this % values % mark_known(this % element, this % new_values)

    result % touches_value     = .true.
    result % touches_structure = .false.
    call result % mark_applied()

  end subroutine apply

  !===================================================================!
  ! Check: the caller-provided verdict - a proof hook, not policy.
  !===================================================================!

  subroutine check(this, result)

    class(value_change), intent(inout) :: this
    type(change_result), intent(inout) :: result

    call result % mark_checked(this % check_passes)

  end subroutine check

  !===================================================================!
  ! Keep: the map already holds the update; keeping is declining
  ! to revert.
  !===================================================================!

  subroutine keep(this, result)

    class(value_change), intent(inout) :: this
    type(change_result), intent(inout) :: result

    associate(unread => this)
    end associate

    call result % mark_kept()

  end subroutine keep

  !===================================================================!
  ! Revert: restore the seat exactly as apply found it. A seat
  ! that did not exist leaves; a seat that stood untrusted stands
  ! untrusted again; a seat that held values holds its old values.
  !===================================================================!

  subroutine revert(this, result)

    class(value_change), intent(inout) :: this
    type(change_result), intent(inout) :: result

    select case (this % old_status)

    case (VALUE_UNATTACHED)

       if (this % values % attached(this % element)) then
          call this % values % detach(this % element)
       end if

    case (VALUE_UNKNOWN)

       if (.not. this % values % attached(this % element)) then
          call this % values % attach_unknown(this % element)
       end if
       call this % values % mark_unknown(this % element)

    case (VALUE_KNOWN)

       if (.not. this % values % attached(this % element)) then
          call this % values % attach_unknown(this % element)
       end if
       call this % values % mark_known(this % element, this % old_values)

    end select

    call result % mark_reverted()

  end subroutine revert

end module graph_value_change
