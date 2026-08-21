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
! The first is a question about two extents; the second is a question
! about what was declared. S = {2,5,6} declared into A = 1..8 is
! extensionally inside every set that holds 2, 5 and 6, and is a
! declared subobject of A alone. This module answers only the second,
! and never infers an edge from the first.
!
!               THE MAP OWNS ITS KEYS, AND ONLY ITS KEYS
!
! A declared embedding is nothing but a directed pair of identities,
!
!     id(S) -> id(A)
!
! so a row stores two identities BY VALUE. It once stored two graph
! pointers, and that made the map a borrower of the very objects it
! was supposed to be an association between - stored outside both, and
! reaching back into both. With the declaring graphs destroyed, the
! walk read freed storage and answered correctly anyway; valgrind
! counted the reads.
!
!     identity map owns its keys by value;
!     it borrows no graph object merely to recognize it.
!
! Every question here - included, declared_into, and the transitive
! order - is a comparison of identities, so all of them are answerable
! from the stored pairs alone. The walk never leaves this map.
!
!                  WHY THERE IS NO ambient_of
!
! An ambient_of(S) -> type(graph) would have to rebuild a graph OBJECT
! from a stored identity, and an identity is not a graph: it says which
! one, not what it is now. Reconstructing the object needs a registry
! of every declared graph, and a global table is a heavier thing than
! the question deserves.
!
! So the question was checked rather than assumed. No production caller
! wanted the object; the one caller in the tower asked
!
!     host = m % ambient_of(s);  host % same_as(a)
!
! which is not a request for a graph at all - it is the identity
! predicate, written in two steps. declared_into answers it in one, and
! the operation that would have needed a registry is gone rather than
! served.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module map_inclusion

  use graph_fractal , only : graph
  use token_identity, only : token, index_of

  implicit none

  private
  public :: inclusion_map, declared_subobject

  !===================================================================!
  ! One declared embedding: which part, which ambient, by value.
  !===================================================================!

  type :: inclusion_pair
     type(token) :: part
     type(token) :: ambient
  end type inclusion_pair

  type :: inclusion_map

     type(inclusion_pair), allocatable, private :: rows(:)

   contains

     procedure :: include_in
     procedure :: included
     procedure :: declared_into

  end type inclusion_map

contains

  !===================================================================!
  ! Declare S c--> A. A set is declared into ONE ambient, so a second
  ! declaration for the same S is refused: two ambients would be two
  ! answers to one question, and the chain would fork.
  !===================================================================!

  subroutine include_in(this, part, ambient)

    class(inclusion_map), intent(inout) :: this
    type(graph)         , intent(in)    :: part
    type(graph)         , intent(in)    :: ambient

    type(inclusion_pair), allocatable :: grown(:)
    type(token)                      :: below, above
    integer                          :: n

    below = part % id()
    above = ambient % id()

    ! An undeclared token does not match itself.
    if (.not. below % matches(below) .or. .not. above % matches(above)) then
       error stop 'map_inclusion: an inclusion is keyed on assigned identity'
    end if

    if (below % matches(above)) then
       error stop 'map_inclusion: a set is not declared into itself'
    end if

    if (row_at(this, below) /= 0) then
       error stop 'map_inclusion: a set is declared into one ambient'
    end if

    if (.not. allocated(this % rows)) allocate(this % rows(0))
    n = size(this % rows)
    allocate(grown(n + 1))
    grown(1:n) = this % rows
    grown(n + 1) % part    = below
    grown(n + 1) % ambient = above
    call move_alloc(grown, this % rows)

  end subroutine include_in

  !===================================================================!
  ! Where the row for an identity is, or zero. The whole of lookup:
  ! the walk below steps in identities, so it needs no other form.
  !===================================================================!

  pure integer function row_at(this, below) result(at)

    class(inclusion_map), intent(in) :: this
    type(token)         , intent(in) :: below

    at = 0
    if (.not. allocated(this % rows)) return

    at = index_of(this % rows % part, below)

  end function row_at

  pure logical function included(this, part)

    class(inclusion_map), intent(in) :: this
    type(graph)         , intent(in) :: part

    type(token) :: below

    below = part % id()
    included = row_at(this, below) /= 0

  end function included

  !===================================================================!
  ! The declared edge itself: was S declared into exactly this A. One
  ! step, not the transitive order - S c--> S' c--> A answers false
  ! here and true to declared_subobject, and the difference is the
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
    at    = row_at(this, below)
    if (at == 0) return

    above  = ambient % id()
    declared = this % rows(at) % ambient % matches(above)

  end function declared_into

  !===================================================================!
  ! THE SUBOBJECT ORDER: reflexive, and transitive along declared
  ! inclusions.
  !
  !     S <= S
  !     S c--> A  and  A <= B   =>   S <= B
  !
  ! The walk is bounded by the number of declared inclusions, because a
  ! chain that revisits a set is a cycle and no set is declared into
  ! itself twice removed.
  !
  ! It steps in IDENTITIES. Each step reads one stored row and compares
  ! two tokens, so the closure of the order needs no graph object other
  ! than the two the caller named - and needs those only to ask for
  ! their tokens.
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

    bound = 0
    if (allocated(m % rows)) bound = size(m % rows)

    do steps = 1, bound
       at = row_at(m, here)
       if (at == 0) return
       here  = m % rows(at) % ambient
       below = here % matches(target_id)
       if (below) return
    end do

    if (row_at(m, here) /= 0) then
       error stop 'map_inclusion: an inclusion chain is finite'
    end if

  end function declared_subobject

end module map_inclusion
