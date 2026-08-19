!=====================================================================!
! GTI ATTACHED-VALUE UPDATE CHANGE
!
! The first value-side client of the reversible change protocol:
! one concrete change that updates an attached value map -
!
!      apply    save the seat's prior state, then vouch the new
!               value (attaching first when no seat exists)
!      check    the caller-provided verdict, a proof hook for now
!      keep     leave the map as updated
!      revert   restore the seat EXACTLY: unattached seats leave,
!               unknown seats stand untrusted again, known seats
!               get their old value back
!
! run by the one generic gti_change_controller - never by machinery
! of its own. This is not a controller and chooses no policy: the
! change owns the rollback details for the attached-value map, the
! map owns identity lookup, status, and value storage, and the
! protocol owns the lifecycle. A pure attached-value change:
! touches_value and never touches_structure - the mirror image of
! adaptive growth's mixed change, proving the same four verbs carry
! both.
!
! Bound state is copied at bind - the element with its identity,
! the new value whole - so a caller mutating its buffer after
! binding changes nothing; only the map is held by pointer, because
! the map is the thing being changed.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_attached_value_changes

  use fractal_graph          , only : graph
  use gti_value_buffers      , only : gti_value_buffer
  use gti_attached_value_maps, only : gti_attached_value_map, &
       & GTI_VALUE_STATUS_UNATTACHED, GTI_VALUE_STATUS_UNKNOWN, &
       & GTI_VALUE_STATUS_KNOWN
  use gti_change_protocols   , only : gti_reversible_change, &
       & gti_change_result

  implicit none

  private
  public :: gti_attached_value_change

  !===================================================================!
  ! One bound update: WHICH map (by pointer - the map is the thing
  ! being changed), WHICH graph (copied, identity and all), WHAT
  ! value (copied at bind), and the saved prior state revert
  ! restores exactly.
  !===================================================================!

  type, extends(gti_reversible_change) :: gti_attached_value_change

     type(gti_attached_value_map), pointer :: values => null()
     type(graph)                           :: element
     type(gti_value_buffer)                :: new_value

     logical :: check_passes = .true.

     integer :: old_status = GTI_VALUE_STATUS_UNATTACHED
     type(gti_value_buffer) :: old_value

   contains

     procedure :: bind
     procedure :: apply
     procedure :: check
     procedure :: keep
     procedure :: revert

  end type gti_attached_value_change

contains

  !===================================================================!
  ! Bind one update: point at the map, copy the element and the new
  ! value, take the check verdict (true unless the caller says
  ! otherwise), and clear any saved state from an earlier life.
  !===================================================================!

  subroutine bind(this, values, element, new_value, check_passes)

    class(gti_attached_value_change)     , intent(inout) :: this
    type(gti_attached_value_map), target , intent(inout) :: values
    type(graph)                          , intent(in)    :: element
    type(gti_value_buffer)               , intent(in)    :: new_value
    logical, optional                    , intent(in)    :: check_passes

    this % values    => values
    this % element   = element
    this % new_value = new_value

    this % check_passes = .true.
    if (present(check_passes)) this % check_passes = check_passes

    this % old_status = GTI_VALUE_STATUS_UNATTACHED
    call this % old_value % clear()

  end subroutine bind

  !===================================================================!
  ! Apply: remember what the seat was, then vouch the new value -
  ! attaching first when no seat exists. The identity gate and the
  ! known-value laws stay the map's own; this change adds only the
  ! memory revert needs.
  !===================================================================!

  subroutine apply(this, result)

    class(gti_attached_value_change), intent(inout) :: this
    type(gti_change_result)         , intent(inout) :: result

    if (.not. associated(this % values)) then
       error stop 'gti_attached_value_change: value map is bound'
    end if

    this % old_status = this % values % status_of(this % element)

    if (this % old_status == GTI_VALUE_STATUS_KNOWN) then
       call this % values % value_of(this % element, this % old_value)
    else
       call this % old_value % clear()
    end if

    if (this % old_status == GTI_VALUE_STATUS_UNATTACHED) then
       call this % values % attach_unknown(this % element)
    end if

    call this % values % mark_known(this % element, this % new_value)

    result % touches_value     = .true.
    result % touches_structure = .false.
    call result % mark_applied()

  end subroutine apply

  !===================================================================!
  ! Check: the caller-provided verdict - a proof hook, not policy.
  !===================================================================!

  subroutine check(this, result)

    class(gti_attached_value_change), intent(inout) :: this
    type(gti_change_result)         , intent(inout) :: result

    call result % mark_checked(this % check_passes)

  end subroutine check

  !===================================================================!
  ! Keep: the map already holds the update; keeping is declining
  ! to revert.
  !===================================================================!

  subroutine keep(this, result)

    class(gti_attached_value_change), intent(inout) :: this
    type(gti_change_result)         , intent(inout) :: result

    associate(unread => this)
    end associate

    call result % mark_kept()

  end subroutine keep

  !===================================================================!
  ! Revert: restore the seat exactly as apply found it. A seat that
  ! did not exist leaves; a seat that stood untrusted stands
  ! untrusted again; a seat that held a value holds its old value.
  !===================================================================!

  subroutine revert(this, result)

    class(gti_attached_value_change), intent(inout) :: this
    type(gti_change_result)         , intent(inout) :: result

    select case (this % old_status)

    case (GTI_VALUE_STATUS_UNATTACHED)

       if (this % values % attached(this % element)) then
          call this % values % detach(this % element)
       end if

    case (GTI_VALUE_STATUS_UNKNOWN)

       if (.not. this % values % attached(this % element)) then
          call this % values % attach_unknown(this % element)
       end if
       call this % values % mark_unknown(this % element)

    case (GTI_VALUE_STATUS_KNOWN)

       if (.not. this % values % attached(this % element)) then
          call this % values % attach_unknown(this % element)
       end if
       call this % values % mark_known(this % element, this % old_value)

    end select

    call result % mark_reverted()

  end subroutine revert

end module gti_attached_value_changes
