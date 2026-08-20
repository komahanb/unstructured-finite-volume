!=====================================================================!
! VALUE MAP: graph identity -> value status x field, stored
! outside the graph. Per graph it records whether a value is
! attached, whether it is trusted, and the value itself, a field
! on the graph's own domain.
!
! Statuses, closed:
!
!     VALUE_UNATTACHED   no row in the map
!     VALUE_UNKNOWN      a row with no trusted value
!     VALUE_KNOWN        a row holding a trusted field
!
! Readers accept absence: an unattached graph reads as
! VALUE_UNATTACHED. Reading a value requires KNOWN, because an
! untrusted number must not be consumed. Writers require an
! assigned identity and (except attach) an existing row, because
! updating a missing row or attaching twice would each leave the
! map ambiguous.
!
! Rows are keyed on type(token) identity tokens copied at attach,
! never on position, so growing, reordering, or compacting the
! container cannot redirect an attachment, and the map may outlive
! every variable that built it. No argument requires TARGET.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module map_value

  use iso_fortran_env  , only : dp => REAL64
  use graph_fractal    , only : graph
  use token_identity   , only : token
  use field_stored, only : stored_field

  implicit none

  private
  public :: value_map
  public :: VALUE_UNATTACHED, VALUE_UNKNOWN, VALUE_KNOWN

  integer, parameter :: VALUE_UNATTACHED = 0
  integer, parameter :: VALUE_UNKNOWN    = 1
  integer, parameter :: VALUE_KNOWN      = 2

  !===================================================================!
  ! One row: the copied identity token, the status, and the value
  ! field.
  !===================================================================!

  type :: value_pair

     type(token) :: identity
     integer     :: status = VALUE_UNKNOWN
     type(stored_field) :: value

  end type value_pair

  type :: value_map

     type(value_pair), allocatable, private :: rows(:)

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
  ! Index of the row keyed on the given token, or zero.
  !===================================================================!

  pure integer function row_at(this, key) result(at)

    class(value_map), intent(in) :: this
    type(token)     , intent(in) :: key

    integer :: k

    at = 0
    if (.not. allocated(this % rows)) return

    do k = 1, size(this % rows)
       if (this % rows(k) % identity % matches(key)) then
          at = k
          return
       end if
    end do

  end function row_at

  !===================================================================!
  ! Check that the element has an assigned identity before any
  ! write: an undeclared token does not match itself, so a row
  ! keyed on it could never be found again. Violation stops the
  ! program.
  !===================================================================!

  pure function writer_key(element) result(key)

    type(graph), intent(in) :: element
    type(token)             :: key

    key = element % id()
    if (.not. key % matches(key)) then
       error stop 'map_value: a value map is keyed on assigned identity'
    end if

  end function writer_key

  !===================================================================!
  ! Add one row for a graph, status UNKNOWN, no value. Attaching
  ! twice stops the program, because two rows for one identity
  ! would make lookups ambiguous.
  !===================================================================!

  subroutine attach_unknown(this, element)

    class(value_map), intent(inout) :: this
    type(graph)     , intent(in)    :: element

    type(value_pair), allocatable :: grown(:)
    type(token)                  :: key
    integer                      :: n

    key = writer_key(element)

    if (row_at(this, key) /= 0) then
       error stop 'map_value: a value row is attached once'
    end if

    if (.not. allocated(this % rows)) allocate(this % rows(0))

    n = size(this % rows)
    allocate(grown(n + 1))
    grown(1:n) = this % rows
    grown(n + 1) % identity = key
    grown(n + 1) % status   = VALUE_UNKNOWN
    call move_alloc(grown, this % rows)

  end subroutine attach_unknown

  !===================================================================!
  ! Store a copy of the given values in a field on the element's
  ! own domain and set the status to KNOWN. Updating an existing
  ! KNOWN row is allowed. Stops the program when the element has
  ! no row, or when values is empty, because KNOWN with no numbers
  ! could not be read back.
  !===================================================================!

  subroutine mark_known(this, element, values, ncomp)

    class(value_map)  , intent(inout) :: this
    type(graph)       , intent(in)    :: element
    real(dp)          , intent(in)    :: values(:)
    integer , optional, intent(in)    :: ncomp

    type(token) :: key
    integer     :: at, width

    key = writer_key(element)

    at = row_at(this, key)
    if (at == 0) then
       error stop 'map_value: an update touches an attached row'
    end if

    if (size(values) == 0) then
       error stop 'map_value: a known value has values'
    end if

    width = 1
    if (present(ncomp)) width = max(ncomp, 1)

    this % rows(at) % value = stored_field('attached value', element, &
         & size(values) / width, ncomp=width)
    call this % rows(at) % value % set_real_vector(values)
    this % rows(at) % status = VALUE_KNOWN

  end subroutine mark_known

  !===================================================================!
  ! Clear the value and set the status back to UNKNOWN; the row
  ! remains. Stops the program when the element has no row.
  !===================================================================!

  subroutine mark_unknown(this, element)

    class(value_map), intent(inout) :: this
    type(graph)     , intent(in)    :: element

    type(stored_field) :: nothing
    type(token) :: key
    integer     :: at

    key = writer_key(element)

    at = row_at(this, key)
    if (at == 0) then
       error stop 'map_value: an update touches an attached row'
    end if

    this % rows(at) % value  = nothing
    this % rows(at) % status = VALUE_UNKNOWN

  end subroutine mark_unknown

  !===================================================================!
  ! Remove the element's row; other rows are unaffected, because
  ! lookup is by token and not by position. Stops the program when
  ! the element has no row.
  !===================================================================!

  subroutine detach(this, element)

    class(value_map), intent(inout) :: this
    type(graph)     , intent(in)    :: element

    type(value_pair), allocatable :: kept(:)
    type(token)                  :: key
    integer                      :: at, n, k, m

    key = writer_key(element)

    at = row_at(this, key)
    if (at == 0) then
       error stop 'map_value: a detach removes an attached row'
    end if

    n = size(this % rows)
    allocate(kept(n - 1))
    m = 0
    do k = 1, n
       if (k == at) cycle
       m = m + 1
       kept(m) = this % rows(k)
    end do
    call move_alloc(kept, this % rows)

  end subroutine detach

  !===================================================================!
  ! Readers. Absence is an accepted input: an undeclared or
  ! unattached graph reads as not present. Only value_of stops the
  ! program, when the status is not KNOWN, because an untrusted
  ! number must not be consumed.
  !===================================================================!

  pure logical function attached(this, element)

    class(value_map), intent(in) :: this
    type(graph)     , intent(in) :: element

    type(token) :: key

    key = element % id()
    attached = row_at(this, key) /= 0

  end function attached

  pure integer function status_of(this, element) result(status)

    class(value_map), intent(in) :: this
    type(graph)     , intent(in) :: element

    type(token) :: key
    integer     :: at

    key = element % id()
    at  = row_at(this, key)

    if (at == 0) then
       status = VALUE_UNATTACHED
    else
       status = this % rows(at) % status
    end if

  end function status_of

  subroutine value_of(this, element, values)

    class(value_map)     , intent(in)  :: this
    type(graph)          , intent(in)  :: element
    real(dp), allocatable, intent(out) :: values(:)

    type(token) :: key
    integer     :: at

    key = element % id()
    at  = row_at(this, key)

    if (at == 0) then
       error stop 'map_value: a known value is read'
    end if

    if (this % rows(at) % status /= VALUE_KNOWN) then
       error stop 'map_value: a known value is read'
    end if

    call this % rows(at) % value % get_real_vector(values)

  end subroutine value_of

end module map_value
