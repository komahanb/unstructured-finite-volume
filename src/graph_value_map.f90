!=====================================================================!
! VALUE MAP
!
! What a graph's value IS, and whether it is trusted. Keyed on
! graph identity, stored outside the graph:
!
!     value_map : graph identity -> value status x field
!
! The second metadata association, filling the seat the label map
! reserved: a label answers WHAT IT IS CALLED, this map answers
! WHAT NUMBER IS ATTACHED AND WHETHER ANYONE VOUCHES FOR IT. The
! split it serves is the tower's:
!
!     Structural knownness lives in fractal_graph.
!     Value knownness lives in an attached value map.
!     Updates are reversible changes; the map owns the ontology,
!     never the lifecycle.
!
! The attached value is a FIELD, minted on the element's own
! domain when it is vouched - the value carries the identity of
! what it values, and no second value carrier exists.
!
!                  KEYED BY IDENTITY, NEVER BY POSITION
!
! No integer index exists in this API. A row is found by the one
! comparison the identity module owns, so reordering, growing, or
! compacting whatever container holds the graphs can never re-aim
! an attachment.
!
!                       THE STATUS VOCABULARY
!
! Three answers, closed:
!
!     VALUE_UNATTACHED   no seat: the map holds nothing
!     VALUE_UNKNOWN      a seat, but no trusted number
!     VALUE_KNOWN        a seat holding a trusted field
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
! Rows key on type(token), copied at attach - the label map's
! storage law, inherited whole. This map borrows no graph object
! merely to recognize it, so it may outlive every variable that
! built it, and no verb demands TARGET.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module graph_value_map

  use iso_fortran_env  , only : dp => REAL64
  use fractal_graph    , only : graph
  use graph_identity   , only : token
  use class_graph_field, only : field

  implicit none

  private
  public :: value_map
  public :: VALUE_UNATTACHED, VALUE_UNKNOWN, VALUE_KNOWN

  integer, parameter :: VALUE_UNATTACHED = 0
  integer, parameter :: VALUE_UNKNOWN    = 1
  integer, parameter :: VALUE_KNOWN      = 2

  !===================================================================!
  ! One row: WHICH graph - by copied token - what status its value
  ! carries, and the value itself, a field on the element's own
  ! domain.
  !===================================================================!

  type :: value_row

     type(token) :: identity
     integer     :: status = VALUE_UNKNOWN
     type(field) :: value

  end type value_row

  type :: value_map

     type(value_row), allocatable, private :: rows(:)

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
  ! Where the row for an identity is, or zero. The one comparison,
  ! and it reads nothing outside this map.
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
  ! The assigned-identity gate every writer passes: an undeclared
  ! token does not match itself, and no seat can key on it.
  !===================================================================!

  pure function writer_key(element) result(key)

    type(graph), intent(in) :: element
    type(token)             :: key

    key = element % id()
    if (.not. key % matches(key)) then
       error stop 'graph_value_map: a value map is keyed on assigned identity'
    end if

  end function writer_key

  !===================================================================!
  ! Open one value seat for a graph: attached, UNKNOWN, empty. A
  ! seat is attached once - a second attach would leave two
  ! answers to one question.
  !===================================================================!

  subroutine attach_unknown(this, element)

    class(value_map), intent(inout) :: this
    type(graph)     , intent(in)    :: element

    type(value_row), allocatable :: grown(:)
    type(token)                  :: key
    integer                      :: n

    key = writer_key(element)

    if (row_at(this, key) /= 0) then
       error stop 'graph_value_map: a value seat is attached once'
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
  ! Vouch for a value: a field is minted on the element's own
  ! domain, holding a copy of the given numbers, and the seat
  ! reads KNOWN. Updating an already-KNOWN seat is lawful - values
  ! move, identity does not - but a KNOWN claim with no numbers
  ! behind it is refused.
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
       error stop 'graph_value_map: an update touches an attached seat'
    end if

    if (size(values) == 0) then
       error stop 'graph_value_map: a known value has values'
    end if

    width = 1
    if (present(ncomp)) width = max(ncomp, 1)

    this % rows(at) % value = field('attached value', element, &
         & size(values) / width, ncomp=width)
    call this % rows(at) % value % set_real_vector(values)
    this % rows(at) % status = VALUE_KNOWN

  end subroutine mark_known

  !===================================================================!
  ! Withdraw trust: the seat stays attached, the number leaves,
  ! and the status reads UNKNOWN again.
  !===================================================================!

  subroutine mark_unknown(this, element)

    class(value_map), intent(inout) :: this
    type(graph)     , intent(in)    :: element

    type(field) :: nothing
    type(token) :: key
    integer     :: at

    key = writer_key(element)

    at = row_at(this, key)
    if (at == 0) then
       error stop 'graph_value_map: an update touches an attached seat'
    end if

    this % rows(at) % value  = nothing
    this % rows(at) % status = VALUE_UNKNOWN

  end subroutine mark_unknown

  !===================================================================!
  ! Close one seat: the row leaves whole, every other attachment
  ! stands - identity lookup owes nothing to position.
  !===================================================================!

  subroutine detach(this, element)

    class(value_map), intent(inout) :: this
    type(graph)     , intent(in)    :: element

    type(value_row), allocatable :: kept(:)
    type(token)                  :: key
    integer                      :: at, n, k, m

    key = writer_key(element)

    at = row_at(this, key)
    if (at == 0) then
       error stop 'graph_value_map: a detach removes an attached seat'
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
       error stop 'graph_value_map: a known value is read'
    end if

    if (this % rows(at) % status /= VALUE_KNOWN) then
       error stop 'graph_value_map: a known value is read'
    end if

    call this % rows(at) % value % get_real_vector(values)

  end subroutine value_of

end module graph_value_map
