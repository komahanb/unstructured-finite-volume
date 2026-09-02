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
! this repository already builds contains (s, s), the same value on both
! sides. So the value map is the identity and needs no storage:
!
!     inclusion_value(s) = s
!
! What differs between S and A is POSITION, and position belongs to
! each set's own representation. Two coordinate systems, and only the
! second is representation-local.
!
!               WHAT THE EXTENSION DOES NOT DETERMINE
!
! Two ambients may have identical extensions and be two domains. Then
! two subsets over the same members - one declared into A, one into B -
! have IDENTICAL inclusion tuples and different ambients. The tuples
! are derivable; the association is not. That is the whole reason this
! map exists rather than being computed.
!
!        DECLARED SUBOBJECT IS NOT EXTENSIONAL CONTAINMENT
!
!     extensional subset   every member of S belongs to A
!     declared subobject   an inclusion path S -> ... -> A exists
!
! The first is a predicate on two extents; the second is a predicate
! on what was declared. S = {2,5,6} declared into A = 1..8 is
! extensionally inside every set that contains 2, 5 and 6, and is a
! declared subobject of A alone. This module evaluates only the
! second, and never infers an edge from the first.
!
!               THE MAP OWNS ITS KEYS, AND ONLY ITS KEYS
!
! A declared embedding is nothing but a directed pair of identities,
!
!     id(S) -> id(A)
!
! so a row stores two identities BY VALUE. The map once stored two
! graph pointers, and that made the map a borrower of the objects it
! was an association between - stored outside both, and referencing
! both. With the declaring graphs deallocated, the traversal read freed
! storage and returned correct values regardless; valgrind counted the
! reads.
!
!     identity map owns its keys by value;
!     it references no graph object in order to identify it.
!
! Every query here - included, declared_into, and the transitive
! order - is a comparison of identities, so all of them are computable
! from the stored pairs alone. The traversal reads only this map.
!
!                  WHY THERE IS NO ambient_of
!
! An ambient_of(S) -> type(graph) would have to rebuild a graph OBJECT
! from a stored identity, and an identity is not a graph: the identity
! states which graph, not the graph's current contents. Reconstructing
! the object needs a registry of every declared graph, and a global
! table costs more than the query requires.
!
! So the requirement was checked rather than assumed. No production
! caller required the object; the one caller in the tower evaluated
!
!     host = m % ambient_of(s);  host % same_as(a)
!
! which is not a request for a graph - it is the identity predicate,
! written in two steps. declared_into evaluates it in one, and the
! operation that would have needed a registry is removed rather than
! implemented.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module map_inclusion

  use graph_fractal , only : graph
  use token_identity, only : token
  use map_token_rows, only : identity_rows

  implicit none

  private
  public :: inclusion_map, declared_subobject

  !===================================================================!
  ! The ambients run parallel to the table's keys: the part is the key
  ! at row at, ambients(at) the ambient it was declared into. Both are
  ! copied tokens.
  !===================================================================!

  type :: inclusion_map

     type(identity_rows)      , private :: rows
     type(token), allocatable , private :: ambients(:)

   contains

     procedure :: include_in
     procedure :: included
     procedure :: declared_into

  end type inclusion_map

contains

  !===================================================================!
  ! Declare S c--> A. A set is declared into ONE ambient, so a second
  ! declaration for the same S is rejected: two ambients would be two
  ! values for one key, and the chain would branch.
  !===================================================================!

  subroutine include_in(this, part, ambient)

    class(inclusion_map), intent(inout) :: this
    type(graph)         , intent(in)    :: part
    type(graph)         , intent(in)    :: ambient

    type(token), allocatable :: grown(:)
    type(token) :: below, above
    integer     :: n, at

    below = part % id()
    above = ambient % id()

    if (.not. below % declared() .or. .not. above % declared()) then
       error stop 'map_inclusion: an inclusion is keyed on assigned identity'
    end if

    if (below % matches(above)) then
       error stop 'map_inclusion: a set is not declared into itself'
    end if

    if (this % rows % position(below) /= 0) then
       error stop 'map_inclusion: a set is declared into one ambient'
    end if

    at = this % rows % append(below)

    if (.not. allocated(this % ambients)) allocate(this % ambients(0))
    n = size(this % ambients)
    allocate(grown(n + 1))
    grown(1:n)   = this % ambients
    grown(n + 1) = above
    call move_alloc(grown, this % ambients)

  end subroutine include_in

  pure logical function included(this, part)

    class(inclusion_map), intent(in) :: this
    type(graph)         , intent(in) :: part

    included = this % rows % position(part % id()) /= 0

  end function included

  !===================================================================!
  ! The declared edge itself: was S declared into exactly this A. One
  ! step, not the transitive order - S c--> S' c--> A returns false
  ! here and true from declared_subobject, and the difference is the
  ! reason both exist.
  !===================================================================!

  pure logical function declared_into(this, part, ambient) result(declared)

    class(inclusion_map), intent(in) :: this
    type(graph)         , intent(in) :: part
    type(graph)         , intent(in) :: ambient

    type(token) :: below, above
    integer     :: at

    declared = .false.

    below = part % id()
    at    = this % rows % position(below)
    if (at == 0) return

    above  = ambient % id()
    declared = this % ambients(at) % matches(above)

  end function declared_into

  !===================================================================!
  ! THE SUBOBJECT ORDER: reflexive, and transitive along declared
  ! inclusions.
  !
  !     S <= S
  !     S c--> A  and  A <= B   =>   S <= B
  !
  ! The traversal is bounded by the number of declared inclusions,
  ! because a chain that revisits a set is a cycle and no set is
  ! declared into itself twice removed.
  !
  ! The traversal steps in IDENTITIES. Each step reads one stored row
  ! and compares two tokens, so the closure of the order needs no graph
  ! object other than the two the caller named - and needs those only
  ! to read their tokens.
  !===================================================================!

  logical function declared_subobject(part, ancestor, m) result(below)

    type(graph)        , intent(in) :: part
    type(graph)        , intent(in) :: ancestor
    type(inclusion_map), intent(in) :: m

    type(token) :: here, target_id
    integer     :: steps, bound, at

    here      = part % id()
    target_id = ancestor % id()

    below = here % matches(target_id)
    if (below) return

    bound = m % rows % num_rows()

    do steps = 1, bound
       at = m % rows % position(here)
       if (at == 0) return
       here  = m % ambients(at)
       below = here % matches(target_id)
       if (below) return
    end do

    if (m % rows % position(here) /= 0) then
       error stop 'map_inclusion: an inclusion chain is finite'
    end if

  end function declared_subobject

end module map_inclusion
