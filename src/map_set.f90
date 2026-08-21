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
!                AND IT BORROWS NOTHING EITHER
!
! Lending nothing OUTWARD was only half the law. This map once keyed
! its rows on
!
!     type(graph), pointer :: element
!
! which borrows the caller's graph INWARD, merely in order to
! recognize it later. Then the map outlived its own key: with the
! binder's graph destroyed, every lookup scanned freed storage, and
! the native allocator answered CORRECTLY - 170 invalid reads under
! valgrind, in row_of and so in all six questions and in bind itself.
! An answer that is right because the page has not been reused yet is
! not an answer.
!
! A row now keeps a COPY of the identity:
!
!     identity map owns its keys by value;
!     it borrows no graph object merely to recognize it.
!
! A token is the whole of what recognition needs - matches is the one
! comparison, and it reads nothing outside this map. So bind no longer
! demands TARGET, and a map may outlive every graph variable that
! built it, which is what an association stored OUTSIDE both parties
! was supposed to mean in the first place.
!
! The token is infrastructure and stays here: nothing in the set view
! takes or returns one. A caller names a set with a type(graph), as
! before - a copy of a graph carries its token, so a copy is the same
! set, and that is exactly why the key may be copied.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module map_set

  use graph_fractal          , only : graph
  use token_identity         , only : token, index_of
  use map_set_representation, only : set_representation

  implicit none

  private
  public :: set_map

  !===================================================================!
  ! One row: WHICH set - by value - and how its members are stored.
  !===================================================================!

  type :: set_pair
     type(token)                             :: identity
     class(set_representation), allocatable  :: extent
  end type set_pair

  type :: set_map

     type(set_pair), allocatable, private :: rows(:)

   contains

     procedure :: bind

     !----------------------------------------------------------------!
     ! The dispatched questions. Each finds the row by identity and
     ! asks the representation, answering a value.
     !----------------------------------------------------------------!

     procedure :: describes
     procedure :: num_members_of
     procedure :: member_of
     procedure :: members_of
     procedure :: has
     procedure :: index_in

     !----------------------------------------------------------------!
     ! The extent itself, COPIED. This is how a compiled representation
     ! - a CSR relation's row numbering - takes the coordinates it will
     ! need for life, at the one moment a map is in scope.
     !
     ! It is not the borrowed accessor the header warns about: the
     ! answer is a fresh allocatable, so the caller owns it and the map
     ! may grow, move or die without touching it. Lending is what would
     ! need a gate. Copying needs only a reason, and compiling one is
     ! the reason.
     !----------------------------------------------------------------!

     procedure :: extent_of

  end type set_map

contains

  !===================================================================!
  ! Bind a representation to a set. A set is described once: a second
  ! binding would leave two answers to one question.
  !===================================================================!

  subroutine bind(this, element, extent)

    class(set_map)           , intent(inout) :: this
    type(graph)              , intent(in)    :: element
    class(set_representation), intent(in)    :: extent

    type(set_pair), allocatable :: grown(:)
    type(token)                :: key
    integer                    :: n

    ! An undeclared token does not match itself.
    key = element % id()
    if (.not. key % matches(key)) then
       error stop 'map_set: a set map is keyed on assigned identity'
    end if

    if (row_of(this, element) /= 0) then
       error stop 'map_set: a set is described once'
    end if

    if (.not. allocated(this % rows)) allocate(this % rows(0))

    n = size(this % rows)
    allocate(grown(n + 1))
    grown(1:n) = this % rows
    grown(n + 1) % identity = key
    allocate(grown(n + 1) % extent, source=extent)
    call move_alloc(grown, this % rows)

  end subroutine bind

  !===================================================================!
  ! Where the pair is, or zero. Private to the dispatch below. The
  ! comparison loop is token_identity's index_of - written once for
  ! every map - so this guards unallocated storage and extracts the
  ! caller's key, nothing more.
  !===================================================================!

  pure integer function row_of(this, element) result(at)

    class(set_map), intent(in) :: this
    type(graph)   , intent(in) :: element

    at = 0
    if (.not. allocated(this % rows)) return

    at = index_of(this % rows % identity, element % id())

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

  integer function num_members_of(this, element)

    class(set_map), intent(in) :: this
    type(graph)   , intent(in) :: element

    integer :: at

    at = row_of(this, element)
    if (at == 0) error stop 'map_set: no representation describes that set'

    num_members_of = this % rows(at) % extent % num_members()

  end function num_members_of

  integer function member_of(this, element, position)

    class(set_map), intent(in) :: this
    type(graph)   , intent(in) :: element
    integer       , intent(in) :: position

    integer :: at

    at = row_of(this, element)
    if (at == 0) error stop 'map_set: no representation describes that set'

    member_of = this % rows(at) % extent % member(position)

  end function member_of

  subroutine members_of(this, element, values)

    class(set_map)      , intent(in)  :: this
    type(graph)         , intent(in)  :: element
    integer, allocatable, intent(out) :: values(:)

    integer :: at

    at = row_of(this, element)
    if (at == 0) error stop 'map_set: no representation describes that set'

    call this % rows(at) % extent % members(values)

  end subroutine members_of

  logical function has(this, element, value)

    class(set_map), intent(in) :: this
    type(graph)   , intent(in) :: element
    integer       , intent(in) :: value

    integer :: at

    at = row_of(this, element)
    if (at == 0) error stop 'map_set: no representation describes that set'

    has = this % rows(at) % extent % has(value)

  end function has

  subroutine extent_of(this, element, extent)

    class(set_map)                        , intent(in)  :: this
    type(graph)                           , intent(in)  :: element
    class(set_representation), allocatable, intent(out) :: extent

    integer :: at

    at = row_of(this, element)
    if (at == 0) error stop 'map_set: no representation describes that set'

    allocate(extent, source=this % rows(at) % extent)

  end subroutine extent_of

  integer function index_in(this, element, value)

    class(set_map), intent(in) :: this
    type(graph)   , intent(in) :: element
    integer       , intent(in) :: value

    integer :: at

    at = row_of(this, element)
    if (at == 0) error stop 'map_set: no representation describes that set'

    index_in = this % rows(at) % extent % local_index(value)

  end function index_in

end module map_set
