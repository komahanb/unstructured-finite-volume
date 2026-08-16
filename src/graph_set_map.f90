!=====================================================================!
! SET MAP
!
! Which representation describes which set. Keyed on graph identity,
! stored outside the graph, and outside the representation:
!
!     set graph  ->  set representation
!
! A graph answers WHICH SET. A representation answers HOW ITS MEMBERS
! ARE STORED. This map is the association between them, and it is the
! only place the two ever meet.
!
!                  IT LENDS NOTHING, AND SO IT IS SIMPLE
!
! Every question is answered by VALUE. No caller receives a pointer
! into this map's storage, so:
!
!     rows may hold their representation as an ALLOCATABLE component
!     the row array may grow by move_alloc and relocate freely
!     intrinsic assignment deep-copies, so no defined assignment
!     nothing is freed twice, so no finalizer
!
! relational_binding needed all of that machinery because it lends
! pointers into its rows; this map needs none of it because it does
! not. Not lending is what removes the lifetime problem, not care.
!
! If a caller ever genuinely needs a borrowed representation, that is
! a new lifetime gate to be written and measured FIRST - not an
! accessor to be added quietly.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module graph_set_map

  use fractal_graph          , only : graph
  use graph_set_representation, only : set_representation

  implicit none

  private
  public :: set_map

  !===================================================================!
  ! One row: the set, and how its members are stored.
  !===================================================================!

  type :: set_row
     type(graph), pointer                    :: element => null()
     class(set_representation), allocatable  :: extent
  end type set_row

  type :: set_map

     type(set_row), allocatable, private :: rows(:)

   contains

     procedure :: bind

     !----------------------------------------------------------------!
     ! The dispatched questions. Each finds the row by identity and
     ! asks the representation, answering a value.
     !----------------------------------------------------------------!

     procedure :: describes
     procedure :: size_of
     procedure :: member_of
     procedure :: members_of
     procedure :: has_in
     procedure :: index_in

  end type set_map

contains

  !===================================================================!
  ! Bind a representation to a set. A set is described once: a second
  ! binding would leave two answers to one question.
  !===================================================================!

  subroutine bind(this, element, extent)

    class(set_map)           , intent(inout)       :: this
    type(graph)              , intent(in), target  :: element
    class(set_representation), intent(in)          :: extent

    type(set_row), allocatable :: grown(:)
    integer                    :: n, k

    ! An undeclared token does not match itself.
    if (.not. element % same_as(element)) then
       error stop 'graph_set_map: a set map is keyed on assigned identity'
    end if

    if (.not. allocated(this % rows)) allocate(this % rows(0))

    do k = 1, size(this % rows)
       if (this % rows(k) % element % same_as(element)) then
          error stop 'graph_set_map: a set is described once'
       end if
    end do

    n = size(this % rows)
    allocate(grown(n + 1))
    grown(1:n) = this % rows
    grown(n + 1) % element => element
    allocate(grown(n + 1) % extent, source=extent)
    call move_alloc(grown, this % rows)

  end subroutine bind

  !===================================================================!
  ! Where the row is, or zero. Private to the dispatch below.
  !===================================================================!

  pure integer function row_of(this, element) result(at)

    class(set_map), intent(in) :: this
    type(graph)   , intent(in) :: element

    integer :: k

    at = 0
    if (.not. allocated(this % rows)) return

    do k = 1, size(this % rows)
       if (this % rows(k) % element % same_as(element)) then
          at = k
          return
       end if
    end do

  end function row_of

  pure logical function describes(this, element)

    class(set_map), intent(in) :: this
    type(graph)   , intent(in) :: element

    describes = row_of(this, element) /= 0

  end function describes

  !===================================================================!
  ! The dispatch. A set with no representation is not a set this map
  ! can answer for, and it says so rather than inventing an extent.
  !===================================================================!

  integer function size_of(this, element)

    class(set_map), intent(in) :: this
    type(graph)   , intent(in) :: element

    integer :: at

    at = row_of(this, element)
    if (at == 0) error stop 'graph_set_map: no representation describes that set'

    size_of = this % rows(at) % extent % size()

  end function size_of

  integer function member_of(this, element, position)

    class(set_map), intent(in) :: this
    type(graph)   , intent(in) :: element
    integer       , intent(in) :: position

    integer :: at

    at = row_of(this, element)
    if (at == 0) error stop 'graph_set_map: no representation describes that set'

    member_of = this % rows(at) % extent % member(position)

  end function member_of

  subroutine members_of(this, element, values)

    class(set_map)      , intent(in)  :: this
    type(graph)         , intent(in)  :: element
    integer, allocatable, intent(out) :: values(:)

    integer :: at

    at = row_of(this, element)
    if (at == 0) error stop 'graph_set_map: no representation describes that set'

    call this % rows(at) % extent % members(values)

  end subroutine members_of

  logical function has_in(this, element, value)

    class(set_map), intent(in) :: this
    type(graph)   , intent(in) :: element
    integer       , intent(in) :: value

    integer :: at

    at = row_of(this, element)
    if (at == 0) error stop 'graph_set_map: no representation describes that set'

    has_in = this % rows(at) % extent % has(value)

  end function has_in

  integer function index_in(this, element, value)

    class(set_map), intent(in) :: this
    type(graph)   , intent(in) :: element
    integer       , intent(in) :: value

    integer :: at

    at = row_of(this, element)
    if (at == 0) error stop 'graph_set_map: no representation describes that set'

    index_in = this % rows(at) % extent % local_index(value)

  end function index_in

end module graph_set_map
