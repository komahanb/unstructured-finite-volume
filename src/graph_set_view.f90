!=====================================================================!
! SET VIEW
!
! A graph read as a finite set, through a representation the map holds:
!
!     set_size(g, m)              how many members
!     set_member(g, m, k)         which member stands at position k
!     set_members(g, m, v)        all of them
!     set_has(g, m, v)            is v one of them
!     set_local_index(g, m, v)    where does v stand, 0 for an outsider
!     set_defined(g, m)           is g described at all
!     set_equivalent(a, b, m)     do two sets have equal members
!
!                       WHAT THE BRANCHES HOLD
!
! Nothing. A set graph's branches may both stay NULL, and under this
! view they do: the extension lives in a representation, and encoding a
! member list, a count or a local numbering into branch(1) would put
! O(N_extent) semantic objects where O(1) is correct.
!
! The kernel PERMITS graph -> branch -> graph. It does not require
! every member to become a graph. That boundary is what makes the
! ontology usable at exascale, and this view is where it is drawn.
!
! UNKNOWN stays unspent. A candidate meaning exists - a set whose
! extent is not known HERE, a halo before exchange - and it is not
! earned, so it is not encoded. Spending UNKNOWN on "the extension is
! stored elsewhere" would say nothing, since that is true of every set
! under this view.
!
!                    IDENTITY IS NOT EQUALITY
!
! set_equivalent compares MEMBERS. same_as compares SETS. Two
! independently declared sets over 1..8 are equivalent and are not the
! same set, and no amount of shared extension makes them one.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module graph_set_view

  use graph_fractal , only : graph
  use map_set , only : set_map

  implicit none

  private
  public :: set_defined, set_size, set_member, set_members
  public :: set_has, set_local_index, set_equivalent

contains

  logical function set_defined(g, m)

    type(graph)  , intent(in) :: g
    type(set_map), intent(in) :: m

    set_defined = m % describes(g)

  end function set_defined

  integer function set_size(g, m)

    type(graph)  , intent(in) :: g
    type(set_map), intent(in) :: m

    set_size = m % size_of(g)

  end function set_size

  integer function set_member(g, m, position)

    type(graph)  , intent(in) :: g
    type(set_map), intent(in) :: m
    integer      , intent(in) :: position

    set_member = m % member_of(g, position)

  end function set_member

  subroutine set_members(g, m, values)

    type(graph)         , intent(in)  :: g
    type(set_map)       , intent(in)  :: m
    integer, allocatable, intent(out) :: values(:)

    call m % members_of(g, values)

  end subroutine set_members

  logical function set_has(g, m, value)

    type(graph)  , intent(in) :: g
    type(set_map), intent(in) :: m
    integer      , intent(in) :: value

    set_has = m % has_in(g, value)

  end function set_has

  integer function set_local_index(g, m, value)

    type(graph)  , intent(in) :: g
    type(set_map), intent(in) :: m
    integer      , intent(in) :: value

    set_local_index = m % index_in(g, value)

  end function set_local_index

  !===================================================================!
  ! Extensional equality, defined here and never mistaken for identity.
  ! Equal members, position for position - the enumeration order is the
  ! representation's, and two representations of one set enumerate it
  ! the same way or they are not representations of one set.
  !===================================================================!

  logical function set_equivalent(a, b, m) result(equal)

    type(graph)  , intent(in) :: a, b
    type(set_map), intent(in) :: m

    integer, allocatable :: va(:), vb(:)

    call m % members_of(a, va)
    call m % members_of(b, vb)

    equal = size(va) .eq. size(vb)
    if (equal) equal = all(va .eq. vb)

  end function set_equivalent

end module graph_set_view
