!=====================================================================!
! INCLUSION MAP
!
! The declared embedding of one set in another:
!
!     S c--> A
!
! keyed on graph identity, stored outside both. It records a SEMANTIC
! ASSOCIATION and nothing else.
!
!                  WHY NO MEMBER TRANSLATION IS STORED
!
! A subset reuses its ambient's member VALUES - the inclusion relation
! this repository already builds holds (s, s), the same value on both
! sides. So the value map is the identity and needs no storage:
!
!     inclusion_value(s) = s
!
! What differs between S and A is POSITION, and position belongs to
! each set's own representation. Two coordinate systems, and only the
! second is representation-local.
!
!               WHAT THE EXTENSION CANNOT TELL YOU
!
! Two ambients may have identical extensions and be two domains. Then
! two subsets over the same members - one carved from A, one from B -
! have IDENTICAL inclusion tuples and different ambients. The tuples
! are derivable; the association is not. That is the whole reason this
! map exists rather than being computed.
!
!        DECLARED SUBOBJECT IS NOT EXTENSIONAL CONTAINMENT
!
!     extensional subset   every member of S belongs to A
!     declared subobject   an inclusion path S -> ... -> A exists
!
! The first is a question about two extents; the second is a question
! about what was declared. S = {2,5,6} carved from A = 1..8 is
! extensionally inside every set that holds 2, 5 and 6, and is a
! declared subobject of A alone. This module answers only the second,
! and never infers an edge from the first.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module graph_inclusion_map

  use fractal_graph, only : graph

  implicit none

  private
  public :: inclusion_map, declared_subobject

  !===================================================================!
  ! One declared embedding.
  !===================================================================!

  type :: inclusion_row
     type(graph), pointer :: part    => null()
     type(graph), pointer :: ambient => null()
  end type inclusion_row

  type :: inclusion_map

     type(inclusion_row), allocatable, private :: rows(:)

   contains

     procedure :: include_in
     procedure :: included
     procedure :: ambient_of

  end type inclusion_map

contains

  !===================================================================!
  ! Declare S c--> A. A set is carved from ONE domain, so a second
  ! declaration for the same S is refused: two ambients would be two
  ! answers to one question, and the chain would fork.
  !===================================================================!

  subroutine include_in(this, part, ambient)

    class(inclusion_map), intent(inout)      :: this
    type(graph)         , intent(in), target :: part
    type(graph)         , intent(in), target :: ambient

    type(inclusion_row), allocatable :: grown(:)
    integer                          :: n

    ! An undeclared token does not match itself.
    if (.not. part % same_as(part) .or. .not. ambient % same_as(ambient)) then
       error stop 'graph_inclusion_map: an inclusion is keyed on assigned identity'
    end if

    if (part % same_as(ambient)) then
       error stop 'graph_inclusion_map: a set is not declared into itself'
    end if

    if (this % included(part)) then
       error stop 'graph_inclusion_map: a set is carved from one domain'
    end if

    if (.not. allocated(this % rows)) allocate(this % rows(0))
    n = size(this % rows)
    allocate(grown(n + 1))
    grown(1:n) = this % rows
    grown(n + 1) % part    => part
    grown(n + 1) % ambient => ambient
    call move_alloc(grown, this % rows)

  end subroutine include_in

  pure integer function row_of(this, part) result(at)

    class(inclusion_map), intent(in) :: this
    type(graph)         , intent(in) :: part

    integer :: k

    at = 0
    if (.not. allocated(this % rows)) return

    do k = 1, size(this % rows)
       if (this % rows(k) % part % same_as(part)) then
          at = k
          return
       end if
    end do

  end function row_of

  pure logical function included(this, part)

    class(inclusion_map), intent(in) :: this
    type(graph)         , intent(in) :: part

    included = row_of(this, part) /= 0

  end function included

  !===================================================================!
  ! The declared ambient, as a copy - which is to say, as the same
  ! declared set: a copy carries the token. Nothing is lent.
  !===================================================================!

  function ambient_of(this, part) result(host)

    class(inclusion_map), intent(in) :: this
    type(graph)         , intent(in) :: part
    type(graph)                      :: host

    integer :: at

    at = row_of(this, part)
    if (at == 0) then
       error stop 'graph_inclusion_map: that set was not declared into anything'
    end if

    host = this % rows(at) % ambient

  end function ambient_of

  !===================================================================!
  ! THE SUBOBJECT ORDER: reflexive, and transitive along declared
  ! inclusions.
  !
  !     S <= S
  !     S c--> A  and  A <= B   =>   S <= B
  !
  ! The walk is bounded by the number of declared inclusions, because a
  ! chain that revisits a set is a cycle and no set is carved from
  ! itself twice removed.
  !===================================================================!

  logical function declared_subobject(part, ancestor, m) result(below)

    type(graph)        , intent(in) :: part
    type(graph)        , intent(in) :: ancestor
    type(inclusion_map), intent(in) :: m

    type(graph) :: here
    integer     :: steps, bound

    below = part % same_as(ancestor)
    if (below) return

    bound = 0
    if (allocated(m % rows)) bound = size(m % rows)

    here = part
    do steps = 1, bound
       if (.not. m % included(here)) return
       here  = m % ambient_of(here)
       below = here % same_as(ancestor)
       if (below) return
    end do

    if (m % included(here)) then
       error stop 'graph_inclusion_map: an inclusion chain is finite'
    end if

  end function declared_subobject

end module graph_inclusion_map
