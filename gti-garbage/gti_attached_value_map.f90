!=====================================================================!
! GRAPH-ATTACHED VALUE MAP
!
! What a graph's value IS, and whether it is trusted. Keyed on
! graph identity, stored outside the graph:
!
!     value_map : graph identity -> value status x value buffer
!
! The label map named the seat this module now fills: "when a
! second kind of metadata earns its place it will be measured and
! named on its own." This is that second association - not a name,
! a number with a status - and the split it serves is the tower's:
!
!     Structural knownness lives in fractal_graph.
!     Value knownness lives in an attached value map.
!     Reversible changes update maps, but do not define the
!     map ontology.
!
!                  KEYED BY IDENTITY, NEVER BY POSITION
!
! No integer index exists in this API. A row is found by the one
! comparison the identity module owns, so reordering, growing, or
! compacting whatever container holds the graphs can never re-aim
! an attachment - and an integer that names a position in one
! array is not an identity at all.
!
!                       THE STATUS VOCABULARY
!
! Three answers, closed:
!
!     GTI_VALUE_STATUS_UNATTACHED   no seat: the map holds nothing
!     GTI_VALUE_STATUS_UNKNOWN      a seat, but no trusted number
!     GTI_VALUE_STATUS_KNOWN        a seat holding a trusted value
!
! Unattached is absence and absence is a lawful answer, exactly as
! a set nobody named is still a set. Reading a value demands
! KNOWN: a number nobody vouched for is not read, it is refused.
! Writers are the strict half - they demand assigned identity and
! an existing seat, because updating nothing and attaching twice
! both leave two answers to one question.
!
!                         THE STORAGE LAW
!
! Rows key on type(token), copied at attach - the graph_label_map
! storage law, inherited whole. This map borrows no graph object
! merely to recognize it, so it may outlive every variable that
! built it, and no verb demands TARGET.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_attached_value_maps

  use fractal_graph    , only : graph
  use graph_identity   , only : token
  use gti_value_buffers, only : gti_value_buffer

  implicit none

  private
  public :: gti_attached_value_map
  public :: GTI_VALUE_STATUS_UNATTACHED
  public :: GTI_VALUE_STATUS_UNKNOWN
  public :: GTI_VALUE_STATUS_KNOWN

  integer, parameter :: GTI_VALUE_STATUS_UNATTACHED = 0
  integer, parameter :: GTI_VALUE_STATUS_UNKNOWN    = 1
  integer, parameter :: GTI_VALUE_STATUS_KNOWN      = 2

  !===================================================================!
  ! One row: WHICH graph - by copied token - what status its value
  ! carries, and the value itself.
  !===================================================================!

  type :: value_row

     type(token)            :: identity
     integer                :: status = GTI_VALUE_STATUS_UNKNOWN
     type(gti_value_buffer) :: value

  end type value_row

  !===================================================================!
  ! The map. The type keeps its public singular name; Fortran
  ! denies a type its host module's name, so the module speaks in
  ! the plural.
  !===================================================================!

  type :: gti_attached_value_map

     type(value_row), allocatable, private :: rows(:)

   contains

     procedure :: attach_unknown
     procedure :: mark_known
     procedure :: mark_unknown
     procedure :: detach
     procedure :: attached
     procedure :: status_of
     procedure :: value_of

  end type gti_attached_value_map

contains

  !===================================================================!
  ! Where the row for an identity is, or zero. The one comparison,
  ! and it reads nothing outside this map.
  !===================================================================!

  pure integer function row_at(this, key) result(at)

    class(gti_attached_value_map), intent(in) :: this
    type(token)                  , intent(in) :: key

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
  ! The assigned-identity gate every writer passes: an undeclared
  ! token does not match itself, and no seat can key on it.
  !===================================================================!

  pure function writer_key(element) result(key)

    type(graph), intent(in) :: element
    type(token)             :: key

    key = element % id()
    if (.not. key % matches(key)) then
       error stop 'gti_attached_value_map: a value map is keyed on assigned identity'
    end if

  end function writer_key

  !===================================================================!
  ! Open one value seat for a graph: attached, UNKNOWN, empty. A
  ! seat is attached once - a second attach would leave two
  ! answers to one question.
  !===================================================================!

  subroutine attach_unknown(this, element)

    class(gti_attached_value_map), intent(inout) :: this
    type(graph)                  , intent(in)    :: element

    type(value_row), allocatable :: grown(:)
    type(token)                  :: key
    integer                      :: n

    key = writer_key(element)

    if (row_at(this, key) /= 0) then
       error stop 'gti_attached_value_map: a value seat is attached once'
    end if

    if (.not. allocated(this % rows)) allocate(this % rows(0))

    n = size(this % rows)
    allocate(grown(n + 1))
    grown(1:n) = this % rows
    grown(n + 1) % identity = key
    grown(n + 1) % status   = GTI_VALUE_STATUS_UNKNOWN
    call move_alloc(grown, this % rows)

  end subroutine attach_unknown

  !===================================================================!
  ! Vouch for a value: the seat holds a copy of the given buffer
  ! and reads KNOWN. Updating an already-KNOWN seat is lawful -
  ! values move, identity does not - but a KNOWN claim with no
  ! numbers behind it is refused.
  !===================================================================!

  subroutine mark_known(this, element, values)

    class(gti_attached_value_map), intent(inout) :: this
    type(graph)                  , intent(in)    :: element
    type(gti_value_buffer)       , intent(in)    :: values

    type(token) :: key
    integer     :: at, n

    key = writer_key(element)

    at = row_at(this, key)
    if (at == 0) then
       error stop 'gti_attached_value_map: an update touches an attached seat'
    end if

    n = 0
    if (allocated(values % rvals)) n = size(values % rvals)
    if (n == 0) then
       error stop 'gti_attached_value_map: a known value has values'
    end if

    this % rows(at) % value  = values
    this % rows(at) % status = GTI_VALUE_STATUS_KNOWN

  end subroutine mark_known

  !===================================================================!
  ! Withdraw trust: the seat stays attached, the number leaves,
  ! and the status reads UNKNOWN again.
  !===================================================================!

  subroutine mark_unknown(this, element)

    class(gti_attached_value_map), intent(inout) :: this
    type(graph)                  , intent(in)    :: element

    type(token) :: key
    integer     :: at

    key = writer_key(element)

    at = row_at(this, key)
    if (at == 0) then
       error stop 'gti_attached_value_map: an update touches an attached seat'
    end if

    call this % rows(at) % value % clear()
    this % rows(at) % status = GTI_VALUE_STATUS_UNKNOWN

  end subroutine mark_unknown

  !===================================================================!
  ! Close one seat: the row leaves whole, every other attachment
  ! stands - identity lookup owes nothing to position.
  !===================================================================!

  subroutine detach(this, element)

    class(gti_attached_value_map), intent(inout) :: this
    type(graph)                  , intent(in)    :: element

    type(value_row), allocatable :: kept(:)
    type(token)                  :: key
    integer                      :: at, n, k, m

    key = writer_key(element)

    at = row_at(this, key)
    if (at == 0) then
       error stop 'gti_attached_value_map: a detach removes an attached seat'
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
  ! The readers. Absence is a lawful answer - an undeclared or
  ! unattached graph is simply not in the map - so the readers
  ! refuse nothing except reading a number nobody vouched for.
  !===================================================================!

  pure logical function attached(this, element)

    class(gti_attached_value_map), intent(in) :: this
    type(graph)                  , intent(in) :: element

    type(token) :: key

    key = element % id()
    attached = row_at(this, key) /= 0

  end function attached

  pure integer function status_of(this, element) result(status)

    class(gti_attached_value_map), intent(in) :: this
    type(graph)                  , intent(in) :: element

    type(token) :: key
    integer     :: at

    key = element % id()
    at  = row_at(this, key)

    if (at == 0) then
       status = GTI_VALUE_STATUS_UNATTACHED
    else
       status = this % rows(at) % status
    end if

  end function status_of

  subroutine value_of(this, element, values)

    class(gti_attached_value_map), intent(in)  :: this
    type(graph)                  , intent(in)  :: element
    type(gti_value_buffer)       , intent(out) :: values

    type(token) :: key
    integer     :: at

    key = element % id()
    at  = row_at(this, key)

    if (at == 0) then
       error stop 'gti_attached_value_map: a known value is read'
    end if

    if (this % rows(at) % status /= GTI_VALUE_STATUS_KNOWN) then
       error stop 'gti_attached_value_map: a known value is read'
    end if

    values = this % rows(at) % value

  end subroutine value_of

end module gti_attached_value_maps
