!=====================================================================!
! VALUE MAP: graph identity -> value status x field, stored
! outside the graph. Per graph it records whether a value is
! attached, whether it is known, and the value itself, a field
! on the graph's own domain.
!
! Statuses, closed:
!
!     VALUE_UNATTACHED   no row in the map
!     VALUE_UNKNOWN      a row with no known value
!     VALUE_KNOWN        a row storing a known field
!
! Readers accept absence: an unattached graph reads as
! VALUE_UNATTACHED. Reading a value requires KNOWN, because an
! unknown number must not be read. Attach requires an assigned
! identity and no row; every other writer requires an existing row,
! because updating a missing row or attaching twice would each leave
! the map ambiguous.
!
! Rows are keyed on type(token) identity tokens copied at attach,
! never on position, so growing, reordering, or compacting the
! container cannot redirect an attachment, and the map may outlive
! every variable that built it. No argument requires TARGET.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module map_value

  use util_precision  , only : dp
  use graph_fractal    , only : graph
  use map_token_rows   , only : identity_rows
  use field_stored, only : stored_field

  implicit none

  private
  public :: value_map
  public :: VALUE_UNATTACHED, VALUE_UNKNOWN, VALUE_KNOWN

  integer, parameter :: VALUE_UNATTACHED = 0
  integer, parameter :: VALUE_UNKNOWN    = 1
  integer, parameter :: VALUE_KNOWN      = 2

  !===================================================================!
  ! The state of one attached value: its status, and the field of
  ! numbers when the status is KNOWN. The states run parallel to the
  ! table's keys.
  !===================================================================!

  type :: value_state

     integer            :: status = VALUE_UNKNOWN
     type(stored_field) :: value

  end type value_state

  type :: value_map

     type(identity_rows)           , private :: rows
     type(value_state), allocatable, private :: states(:)

   contains

     procedure :: attach_unknown
     procedure :: mark_known
     procedure :: mark_unknown
     procedure :: detach
     procedure :: attached
     procedure :: status_of
     procedure :: value_of

  end type value_map

contains

  !===================================================================!
  ! Add one row for a graph, status UNKNOWN, no value. Attaching
  ! twice stops the program, because two rows for one identity
  ! would make lookups ambiguous.
  !===================================================================!

  subroutine attach_unknown(this, element)

    class(value_map), intent(inout) :: this
    type(graph)     , intent(in)    :: element

    integer :: at

    at = this % rows % append(element % id(), &
         & 'map_value: a value map is keyed on assigned identity', &
         & 'map_value: a value row is attached once')

    if (.not. allocated(this % states)) allocate(this % states(0))
    this % states = [this % states, value_state()]

  end subroutine attach_unknown

  !===================================================================!
  ! Store a copy of the given values in a field on the element's
  ! own domain and set the status to KNOWN. Updating an existing
  ! KNOWN row is allowed. Stops the program when the element has
  ! no row, or when values is empty, because KNOWN with no numbers
  ! could not be read back.
  !===================================================================!

  subroutine mark_known(this, element, values, num_components)

    class(value_map)  , intent(inout) :: this
    type(graph)       , intent(in)    :: element
    real(dp)          , intent(in)    :: values(:)
    integer , optional, intent(in)    :: num_components

    integer :: at, width

    at = this % rows % row(element % id(), 'map_value: an update requires an attached row')

    if (size(values) == 0) then
       error stop 'map_value: a known value has values'
    end if

    width = 1
    if (present(num_components)) width = max(num_components, 1)

    this % states(at) % value = stored_field('attached value', element, &
         & size(values) / width, num_components=width)
    call this % states(at) % value % set_real_vector(values)
    this % states(at) % status = VALUE_KNOWN

  end subroutine mark_known

  !===================================================================!
  ! Clear the value and set the status back to UNKNOWN; the row
  ! remains. Stops the program when the element has no row.
  !===================================================================!

  subroutine mark_unknown(this, element)

    class(value_map), intent(inout) :: this
    type(graph)     , intent(in)    :: element

    type(stored_field) :: empty_field
    integer :: at

    at = this % rows % row(element % id(), 'map_value: an update requires an attached row')

    this % states(at) % value  = empty_field
    this % states(at) % status = VALUE_UNKNOWN

  end subroutine mark_unknown

  !===================================================================!
  ! Remove the element's row; other rows are unaffected, because
  ! lookup is by token and not by position. Stops the program when
  ! the element has no row.
  !===================================================================!

  subroutine detach(this, element)

    class(value_map), intent(inout) :: this
    type(graph)     , intent(in)    :: element

    integer :: at

    at = this % rows % row(element % id(), 'map_value: a detach removes an attached row')

    call this % rows % remove(at)
    this % states = [this % states(1:at - 1), this % states(at + 1:)]

  end subroutine detach

  !===================================================================!
  ! Readers. Absence is an accepted input: an undeclared or
  ! unattached graph reads as not present. Only value_of stops the
  ! program, when the status is not KNOWN, because an unknown
  ! number must not be read.
  !===================================================================!

  pure logical function attached(this, element)

    class(value_map), intent(in) :: this
    type(graph)     , intent(in) :: element

    attached = this % rows % position(element % id()) /= 0

  end function attached

  pure integer function status_of(this, element) result(status)

    class(value_map), intent(in) :: this
    type(graph)     , intent(in) :: element

    integer :: at

    at = this % rows % position(element % id())

    if (at == 0) then
       status = VALUE_UNATTACHED
    else
       status = this % states(at) % status
    end if

  end function status_of

  subroutine value_of(this, element, values)

    class(value_map)     , intent(in)  :: this
    type(graph)          , intent(in)  :: element
    real(dp), allocatable, intent(out) :: values(:)

    integer :: at

    at = this % rows % row(element % id(), 'map_value: a known value is read')

    if (this % states(at) % status /= VALUE_KNOWN) then
       error stop 'map_value: a known value is read'
    end if

    call this % states(at) % value % real_vector(values)

  end subroutine value_of

end module map_value
