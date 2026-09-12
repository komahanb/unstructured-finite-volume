!=====================================================================!
! LABEL MAP
!
! What a set is NAMED. Keyed on graph identity, stored outside the
! graph:
!
!     set graph  ->  character label
!
! The third association of the family, and the smallest. A graph
! determines WHICH object; a representation determines HOW ITS
! MEMBERS ARE STORED; this map determines WHAT IT IS NAMED, and
! nothing here is ever read to decide the first.
!
!                    A LABEL IS NOT AN IDENTITY
!
! Two graphs may have the same label and remain two sets, exactly as
! two graphs with equal extensions remain two sets. So this map is
! searched BY IDENTITY and never by name: there is no lookup from a
! label back to a graph, by design, because such a lookup would have
! to return one graph for a string that names two. Naming is not
! addressing.
!
!     same label, different graph   ->  still different
!     copied token, same graph      ->  same label
!
!                    WHY THIS MAP EXISTS
!
! The domain capability audit found one requirement the re-rooted
! domain law could not satisfy. type(graph) by design stores no name -
! the kernel type is branch(2) and a private token - yet production
! code names a derived subset after the domain it is derived from:
!
!     sg = subset_set(dom % name(), global_carrier, members(1:n))
!
! That is metadata, and metadata belongs beside the mathematics rather
! than inside it. Putting a label in graph_fractal would add to the
! core a component it does not need; putting one in set_representation
! would make two descriptions of one set disagree about its name. So
! the label is stored here, orthogonal to both, and a consumer that
! names nothing never allocates this map.
!
! This is a LABEL map and not an attribute system. One noun, one role.
! When a second kind of metadata is required it will be measured and
! named separately, not added through a generalization.
!
!                        THE UNNAMED RESULT
!
! An unbound graph is not an error - it is a set with no name.
! label_of returns the empty string, matching the metadata convention
! the carriers already follow (name() returns '' when no name was
! set). Binding is the strict half: a set is named ONCE, and an
! undeclared token is rejected, because both would leave two results
! for one query.
!
!                          THE STORAGE LAW
!
! Rows key on type(token), copied at bind. This map references no
! graph object to recognise it, so it may outlive every variable that
! built it, and bind requires no TARGET.
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

  !===================================================================!
  ! The labels run parallel to the table's keys: labels(at) is the
  ! label of the set whose token the table stores at row at.
  !===================================================================!

  type :: label_map

     type(identity_rows)         , private :: rows
     type(string), allocatable, private :: labels(:)

   contains

     procedure :: bind    => bind_label
     procedure :: labelled
     procedure :: label_of

  end type label_map

contains

  !===================================================================!
  ! Name a set. A set is named once: a second binding would leave two
  ! results for one query, and the later one would replace the first
  ! without an error.
  !===================================================================!

  subroutine bind_label(this, element, label)

    class(label_map), intent(inout) :: this
    type(graph)     , intent(in)    :: element
    character(len=*), intent(in)    :: label

    type(string), allocatable :: extended_labels(:)
    integer :: n, at

    at = this % rows % append(element % id(), &
         & 'map_label: a label map is keyed on assigned identity', &
         & 'map_label: a set is named once')

    ! type(string) is finalizable, so the array-constructor grow is
    ! not admitted under -std=f2023; the payload grows by move_alloc,
    ! doubling its capacity so a bind costs amortised constant time.
    if (.not. allocated(this % labels)) allocate(this % labels(max(at, 8)))
    if (at > size(this % labels)) then
       n = size(this % labels)
       allocate(extended_labels(2 * n))
       extended_labels(1:n) = this % labels(1:n)
       call move_alloc(extended_labels, this % labels)
    end if
    this % labels(at) = string(label)

  end subroutine bind_label

  pure logical function labelled(this, element)

    class(label_map), intent(in) :: this
    type(graph)     , intent(in) :: element

    labelled = this % rows % position(element % id()) /= 0

  end function labelled

  !===================================================================!
  ! The set's name, as a value; '' when no name was bound. The
  ! unnamed result is not a rejection, because an absent name is not
  ! a contradiction - metadata may be absent.
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
