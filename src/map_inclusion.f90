!=====================================================================!
! Module that states the inclusion map: the declared embedding of one
! set in another,
!
!     S c--> A
!
! keyed on graph identity and stored outside both. It records a
! semantic association and nothing else.
!
! Why no member translation is stored. A subset reuses its ambient's
! member values: the inclusion relation this repository already builds
! contains (s, s), the same value on both sides. So the value map is
! the identity and needs no storage,
!
!     inclusion_value(s) = s
!
! What differs between S and A is position, and position belongs to
! each set's own representation. There are two coordinate systems, and
! only the second is local to a representation.
!
! What the extension does not determine. Two ambients may have
! identical extensions and be two domains. Then two subsets over the
! same members, one declared into A and one into B, have identical
! inclusion tuples and different ambients. The tuples are derivable;
! the association is not. That is the whole reason this map exists
! rather than being computed.
!
! A declared subobject is not extensional containment:
!
!     extensional subset   every member of S belongs to A
!     declared subobject   an inclusion path S -> ... -> A exists
!
! The first is a predicate on two extents; the second is a predicate
! on what was declared. S = {2,5,6} declared into A = 1..8 is
! extensionally inside every set that contains 2, 5 and 6, and is a
! declared subobject of A alone. This module evaluates only the
! second, and never infers an edge from the first:
!
!   S ---c--> A ---c--> B        S <= S,  S <= A,  S <= B
!   T ---c--> A                  T <= A,  T <= B;  S and T unordered
!
! The map owns its keys, and only its keys. A declared embedding is
! nothing but a directed pair of identities,
!
!     id(S) -> id(A)
!
! so a row stores two identities by value. The map once stored two
! graph pointers, which made it a borrower of the objects it was an
! association between: stored outside both, and referencing both. With
! the declaring graphs deallocated, the traversal read freed storage
! and returned correct values regardless; valgrind counted the reads.
!
!     an identity map owns its keys by value;
!     it references no graph object in order to identify it.
!
! The one query here, the transitive order, is a comparison of
! identities, so it is computable from the stored pairs alone. The
! traversal reads only this map.
!
! Why there is no ambient_of. An ambient_of(S) -> type(graph) would
! have to rebuild a graph object from a stored identity, and an
! identity is not a graph: it states which graph, not the graph's
! current contents. Rebuilding the object needs a registry of every
! declared graph, and a global table costs more than the query
! requires.
!
! So the requirement was checked rather than assumed. No production
! caller required the object; the one caller in the tower evaluated
!
!     host = m % ambient_of(s);  host % same_as(a)
!
! which is not a request for a graph but the identity predicate,
! written in two steps. declared_subobject evaluates it in one, and the
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

  !-------------------------------------------------------------------!
  ! Type that models the inclusion map. The ambients run parallel to
  ! the table's keys: the part is the key at row at, and ambients(at)
  ! the ambient it was declared into. Both are copied tokens.
  !
  !   rows                 ambients
  !   at 1   id(S)   ->    id(A)
  !   at 2   id(A)   ->    id(B)
  !   ...
  !-------------------------------------------------------------------!

  type :: inclusion_map

     type(identity_rows)     , private :: rows
     type(token), allocatable, private :: ambients(:)

   contains

     procedure :: include_in

  end type inclusion_map

contains

  !===================================================================!
  ! Declares S c--> A. A set is declared into one ambient, so a second
  ! declaration for the same S is rejected: two ambients would be two
  ! values for one key, and the chain would branch. A set cannot be
  ! declared into itself.
  !===================================================================!

  subroutine include_in(this, part, ambient)

    class(inclusion_map), intent(inout) :: this
    type(graph)         , intent(in)    :: part
    type(graph)         , intent(in)    :: ambient

    type(token)              :: ambient_id
    type(token), allocatable :: extended_ambients(:)
    integer                  :: at

    ambient_id = ambient % id()

    if (.not. ambient_id % declared()) then
       error stop 'map_inclusion: include_in requires ambient to have an assigned identity, &
            &but ambient_id % declared() is false'
    end if

    if (part % same_as(ambient)) then
       error stop 'map_inclusion: part and ambient are the same set (part % same_as(ambient) &
            &is true); a set cannot be declared into itself'
    end if

    at = this % rows % append(part % id(), &
         & 'map_inclusion: append requires the added set to have an assigned identity', &
         & 'map_inclusion: this set is already declared into an ambient; it cannot be declared twice')

    ! The ambients double their capacity, so an inclusion costs
    ! amortised constant time
    if (.not. allocated(this % ambients)) allocate(this % ambients(max(at, 8)))
    if (at > size(this % ambients)) then
       allocate(extended_ambients(2 * size(this % ambients)))
       extended_ambients(1:at - 1) = this % ambients(1:at - 1)
       call move_alloc(extended_ambients, this % ambients)
    end if
    this % ambients(at) = ambient_id

  end subroutine include_in

  !===================================================================!
  ! Returns whether part lies below ancestor in the subobject order,
  ! which is reflexive and transitive along declared inclusions:
  !
  !     S <= S
  !     S c--> A  and  A <= B   =>   S <= B
  !
  ! The traversal is bounded by the number of declared inclusions,
  ! because a chain that revisits a set is a cycle and no set is
  ! declared into itself twice removed; a chain still going at the
  ! bound stops the program.
  !
  ! The traversal steps in identities. Each step reads one stored row
  ! and compares two tokens, so the closure of the order needs no graph
  ! object other than the two the caller named, and needs those only to
  ! read their tokens:
  !
  !   current = id(part)
  !      |  current == id(ancestor)?  yes: true
  !      v  no: current = the ambient stored at current's row
  !   repeat; a current with no row: false
  !===================================================================!

  logical function declared_subobject(part, ancestor, m) result(is_subobject)

    type(graph)        , intent(in) :: part
    type(graph)        , intent(in) :: ancestor
    type(inclusion_map), intent(in) :: m

    type(token)        :: current_id, target_id
    integer            :: steps, bound, at
    character(len=250) :: message

    current_id = part % id()
    target_id  = ancestor % id()

    is_subobject = current_id % matches(target_id)
    if (is_subobject) return

    bound = m % rows % num_rows()

    do steps = 1, bound
       at = m % rows % position(current_id)
       if (at == 0) return
       current_id   = m % ambients(at)
       is_subobject = current_id % matches(target_id)
       if (is_subobject) return
    end do

    if (m % rows % position(current_id) /= 0) then
       write(message,'(a,i0,a)') 'map_inclusion: the inclusion chain still names a declared set &
            &after ', bound, ' step(s) (one per declared row), so it revisits a set - a cycle'
       error stop trim(message)
    end if

  end function declared_subobject

end module map_inclusion
