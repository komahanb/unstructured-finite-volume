!=====================================================================!
! SET MAP
!
! Which representation describes which set. Keyed on graph identity,
! stored outside the graph, and outside the representation:
!
!     set graph  ->  set representation
!
! A graph identifies WHICH SET. A representation specifies HOW ITS
! MEMBERS ARE STORED. This map is the association between them, and it
! is the only module where the two are combined.
!
!        THE MAP RETURNS NO REFERENCES, SO ITS STORAGE IS SIMPLE
!
! Every query is returned by VALUE. No caller receives a pointer
! into this map's storage, so:
!
!     rows may store their representation as an ALLOCATABLE component
!     the row array may grow by move_alloc and relocate freely
!     intrinsic assignment deep-copies, so no defined assignment
!     nothing is freed twice, so no finalizer
!
! relational_binding needs all of that because it returns pointers
! into its rows; this map needs none of it because it returns none.
! Returning no reference is what removes the lifetime problem.
!
! If a caller needs a non-owning representation, that is a new lifetime
! check to be written and measured FIRST - not an accessor to be added
! without one.
!
!                AND THE MAP REFERENCES NOTHING EITHER
!
! Returning no reference OUTWARD was only half the law. This map once
! keyed its rows on
!
!     type(graph), pointer :: element
!
! which references the caller's graph INWARD, in order to identify it
! later. Then the map outlived its own key: with the binder's graph
! deallocated, every lookup scanned freed storage, and the native
! allocator returned CORRECT values - 170 invalid reads under
! valgrind, in every lookup and so in all six queries and in bind itself.
! A value that is correct because the page has not been reused yet is
! not a valid result.
!
! A row now stores a COPY of the identity:
!
!     identity map owns its keys by value;
!     it references no graph object in order to identify it.
!
! A token is all that identification needs - matches is the one
! comparison, and it reads nothing outside this map. So bind no longer
! requires TARGET, and a map may outlive every graph variable that
! built it, which is what an association stored OUTSIDE both objects
! means.
!
! The token is infrastructure and stays here: nothing in the set view
! takes or returns one. A caller names a set with a type(graph), as
! before - a copy of a graph stores its token, so a copy is the same
! set, and that is exactly why the key may be copied.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module map_set

  use graph_fractal          , only : graph
  use map_token_rows         , only : identity_rows
  use map_set_representation, only : set_representation

  implicit none

  private
  public :: set_map

  !===================================================================!
  ! The extent of one set: how its members are stored. The extents run
  ! parallel to the table's keys, and each stores its own kind of
  ! representation.
  !===================================================================!

  type :: extent
     class(set_representation), allocatable :: representation
  end type extent

  type :: set_map

     type(identity_rows)       , private :: rows
     type(extent), allocatable , private :: extents(:)

   contains

     procedure :: bind

     !----------------------------------------------------------------!
     ! The dispatched queries. Each finds the row by identity and
     ! calls the representation, returning a value.
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
     ! need for its whole lifetime, at the one point where a map is in
     ! scope.
     !
     ! It is not the non-owning accessor the header describes: the
     ! result is a newly allocated copy, so the caller owns it and the
     ! map may grow, relocate or be deallocated without affecting it.
     ! Returning a reference is what would need a lifetime check.
     ! Copying needs only a reason, and compiling one is the reason.
     !----------------------------------------------------------------!

     procedure :: extent_of

  end type set_map

contains

  !===================================================================!
  ! Bind a representation to a set. A set is described once: a second
  ! binding would leave two representations for one key.
  !===================================================================!

  subroutine bind(this, element, representation)

    class(set_map)           , intent(inout) :: this
    type(graph)              , intent(in)    :: element
    class(set_representation), intent(in)    :: representation

    type(extent), allocatable :: extended_extents(:)
    integer :: n, at

    at = this % rows % append(element % id(), &
         & 'map_set: a set map is keyed on assigned identity', &
         & 'map_set: a set is described once')

    if (.not. allocated(this % extents)) allocate(this % extents(0))
    n = size(this % extents)
    allocate(extended_extents(n + 1))
    extended_extents(1:n) = this % extents
    allocate(extended_extents(n + 1) % representation, source=representation)
    call move_alloc(extended_extents, this % extents)

  end subroutine bind

  pure logical function describes(this, element)

    class(set_map), intent(in) :: this
    type(graph)   , intent(in) :: element

    describes = this % rows % position(element % id()) /= 0

  end function describes

  !===================================================================!
  ! The dispatch. A set with no representation is not a set this map
  ! can evaluate, and the program stops rather than fabricating an
  ! extent.
  !===================================================================!

  integer function num_members_of(this, element)

    class(set_map), intent(in) :: this
    type(graph)   , intent(in) :: element

    integer :: at

    at = this % rows % row(element % id(), 'map_set: no representation describes that set')

    num_members_of = this % extents(at) % representation % num_members()

  end function num_members_of

  integer function member_of(this, element, position)

    class(set_map), intent(in) :: this
    type(graph)   , intent(in) :: element
    integer       , intent(in) :: position

    integer :: at

    at = this % rows % row(element % id(), 'map_set: no representation describes that set')

    member_of = this % extents(at) % representation % member(position)

  end function member_of

  subroutine members_of(this, element, values)

    class(set_map)      , intent(in)  :: this
    type(graph)         , intent(in)  :: element
    integer, allocatable, intent(out) :: values(:)

    integer :: at

    at = this % rows % row(element % id(), 'map_set: no representation describes that set')

    call this % extents(at) % representation % members(values)

  end subroutine members_of

  logical function has(this, element, value)

    class(set_map), intent(in) :: this
    type(graph)   , intent(in) :: element
    integer       , intent(in) :: value

    integer :: at

    at = this % rows % row(element % id(), 'map_set: no representation describes that set')

    has = this % extents(at) % representation % has(value)

  end function has

  subroutine extent_of(this, element, extent)

    class(set_map)                        , intent(in)  :: this
    type(graph)                           , intent(in)  :: element
    class(set_representation), allocatable, intent(out) :: extent

    integer :: at

    at = this % rows % row(element % id(), 'map_set: no representation describes that set')

    allocate(extent, source=this % extents(at) % representation)

  end subroutine extent_of

  integer function index_in(this, element, value)

    class(set_map), intent(in) :: this
    type(graph)   , intent(in) :: element
    integer       , intent(in) :: value

    integer :: at

    at = this % rows % row(element % id(), 'map_set: no representation describes that set')

    index_in = this % extents(at) % representation % local_index(value)

  end function index_in

end module map_set
