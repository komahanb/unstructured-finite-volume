!=====================================================================!
! Module that states the label map: what a set is named. It is keyed
! on graph identity and stored outside the graph:
!
!     set graph  ->  character label
!
! It is the third association of the family, and the smallest. A
! graph determines which object; a representation determines how its
! members are stored; this map determines what it is named, and
! nothing here is ever read to decide the first.
!
! A label is not an identity. Two graphs may have the same label and
! remain two sets, as two graphs with equal extensions remain two
! sets. So this map is searched by identity and never by name: there
! is no lookup from a label back to a graph, by design, because such a
! lookup would have to return one graph for a string that names two.
! Naming is not addressing.
!
!     same label, different graph   ->  still different
!     copied token, same graph      ->  same label
!
! Why this map exists. The domain capability audit found one
! requirement the re-rooted domain law could not satisfy. type(graph)
! stores no name by design, its kernel being branch(2) and a private
! token, yet production code names a derived subset after the domain
! it is derived from:
!
!     sg = subset_set(dom % name(), global_carrier, members(1:n))
!
! That is metadata, and metadata belongs beside the mathematics rather
! than inside it. A label in graph_fractal would add to the core a
! component it does not need; a label in set_representation would let
! two descriptions of one set disagree about its name. So the label is
! stored here, orthogonal to both, and a consumer that names nothing
! never allocates this map.
!
! This is a label map and not an attribute system: one noun, one role.
! When a second kind of metadata is required it will be measured and
! named separately, not added through a generalization.
!
! The unnamed result. An unbound graph is not an error but a set with
! no name: label_of returns the empty string, as the carriers already
! do (name() returns '' when no name was set). Binding is the strict
! half: a set is named once, and an undeclared token is rejected,
! because both would leave two results for one query.
!
! The storage law. Rows key on type(token), copied at bind. The map
! references no graph object to recognise it, so it may outlive every
! variable that built it, and bind requires no TARGET.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module map_label

  use graph_fractal , only : graph
  use map_token_rows, only : identity_rows
  use util_string   , only : string

  implicit none

  private

  public :: label_map

  !-------------------------------------------------------------------!
  ! Type that models the label map. The labels run parallel to the
  ! table's keys: labels(at) is the label of the set whose token the
  ! table stores at row at.
  !
  !   rows            labels
  !   at 1   token -> 'domain'
  !   at 2   token -> 'boundary'
  !   ...
  !-------------------------------------------------------------------!

  type :: label_map

     type(identity_rows)      , private :: rows
     type(string), allocatable, private :: labels(:)

   contains

     procedure :: bind    => bind_label
     procedure :: labelled
     procedure :: label_of

  end type label_map

contains

  !===================================================================!
  ! Names a set. A set is named once: a second binding would leave two
  ! results for one query, and the later one would replace the first
  ! without an error.
  !===================================================================!

  subroutine bind_label(this, element, label)

    class(label_map), intent(inout) :: this
    type(graph)     , intent(in)    :: element
    character(len=*), intent(in)    :: label

    type(string), allocatable :: extended_labels(:)
    integer                   :: n, at

    at = this % rows % append(element % id(), &
         & 'map_label: a label map is keyed on assigned identity', &
         & 'map_label: a set is named once')

    ! type(string) is finalizable, so growing by an array constructor is
    ! not admitted under -std=f2023. The labels grow by move_alloc
    ! instead, doubling their capacity so that a bind costs amortised
    ! constant time:
    !
    !   capacity   8 -> 16 -> 32 -> ...
    if (.not. allocated(this % labels)) allocate(this % labels(max(at, 8)))
    if (at > size(this % labels)) then
       n = size(this % labels)
       allocate(extended_labels(2 * n))
       extended_labels(1:n) = this % labels(1:n)
       call move_alloc(extended_labels, this % labels)
    end if
    this % labels(at) = string(label)

  end subroutine bind_label

  !===================================================================!
  ! Returns whether a name was bound to the set
  !===================================================================!

  pure logical function labelled(this, element)

    class(label_map), intent(in) :: this
    type(graph)     , intent(in) :: element

    labelled = this % rows % position(element % id()) /= 0

  end function labelled

  !===================================================================!
  ! Returns the set's name, as a value, or '' when no name was bound.
  ! The unnamed result is not a rejection, because an absent name is
  ! not a contradiction: metadata may be absent.
  !===================================================================!

  pure function label_of(this, element) result(label)

    class(label_map), intent(in)  :: this
    type(graph)     , intent(in)  :: element
    character(len=:), allocatable :: label

    integer :: at

    at = this % rows % position(element % id())

    if (at == 0) then
       label = ''
    else
       label = this % labels(at) % str
    end if

  end function label_of

end module map_label
