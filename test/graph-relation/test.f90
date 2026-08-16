!=====================================================================!
! The relation suite: the laws of the finite-arity relation
! (AGENTS.md, level 1, phase 2).
!
! A relation is a named subset of a cartesian product, first-class:
! no graph stands anywhere in this file. The checks are the
! acceptance laws of sections 61 and 62 that phase 2 owes - domain
! validity, arity, membership, verbatim tuples, signature identity -
! plus the design commitments the reviews demanded: carriers shared
! across relations without an owning graph, adjacency and incidence
! as readings of one primitive, the boundary written with no
! imaginary far-side member, SET semantics at the interface (a
! duplicate tuple collapses), and a signature genuinely generic
! over the REPRESENTATION - proved by mixing a counted set with a
! listed one in one relation.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_graph_relation

  use graph_identity    , only : token
  use fractal_graph        , only : set_graph => graph
  use graph_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use graph_set_map        , only : set_map
  use graph_label_map      , only : label_map
  use graph_relation    , only : stored_relation

  implicit none

  integer :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "graph relation suite (AGENTS phase 2)"
  write(*,'(1x,a)') "============================================="

  call check_incidence_contract(nfail)
  call check_set_semantics(nfail)
  call check_carrier_set_semantics(nfail)
  call check_mixed_carriers(nfail)
  call check_adjacency_is_a_reading(nfail)
  call check_ternary_endpoint(nfail)
  call check_shared_carriers(nfail)
  call check_relation_identity(nfail)
  call check_empty_relation(nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all relation checks passed"
  else
     write(*,'(1x,a,i0,a)') "FAILED: ", nfail, " relation check(s)"
     error stop
  end if

contains

  subroutine report(ok, label, nfail)

    logical         , intent(in)    :: ok
    character(len=*), intent(in)    :: label
    integer         , intent(inout) :: nfail

    if (ok) then
       write(*,'(1x,a,a)') "PASS : ", label
    else
       write(*,'(1x,a,a)') "FAIL : ", label
       nfail = nfail + 1
    end if

  end subroutine report

  !===================================================================!
  ! The contract on a cells x faces incidence: identity, arity,
  ! membership, tuples verbatim, and a signature that answers by
  ! structural identity.
  !===================================================================!

  subroutine check_incidence_contract(nfail)

    integer, intent(inout) :: nfail

    type(set_graph)              :: cells, faces
    type(stored_relation)          :: r
    type(set_graph) :: d
    type(token)                    :: tk
    integer, allocatable           :: t(:,:)
    type(set_map)     :: sets

    call cells % declare()
    call sets % bind(cells, counted_set_representation(4))
    call faces % declare()
    call sets % bind(faces, counted_set_representation(5))

    r = stored_relation('touches', [cells, faces], &
         & reshape([1,1,  1,2,  2,2,  3,4], [2, 4]), sets)

    tk = r % id()
    call report(tk % declared(), &
         & "a declared relation carries a token that has signed", nfail)
    call report(r % name() == 'touches', &
         & "and the name it was declared with", nfail)
    call report(r % arity() .eq. 2, &
         & "the arity is the signature's length", nfail)
    call report(r % num_tuples() .eq. 4, &
         & "the count is the tuples handed in", nfail)

    call report(r % has([1, 2]), &
         & "has finds a tuple that is there", nfail)
    call report(.not. r % has([2, 1]), &
         & "order matters: the mirrored pair is absent", nfail)
    call report(.not. r % has([4, 5]), &
         & "two valid members do not make a present tuple", nfail)
    call report(.not. r % has([1]), &
         & "a tuple of the wrong length belongs to nothing", nfail)

    call r % tuples(t)
    call report(size(t, 1) .eq. 2 .and. size(t, 2) .eq. 4 .and. &
         &      all(t(:, 3) .eq. [2, 2]) .and. all(t(:, 4) .eq. [3, 4]), &
         & "tuples come back as handed, order kept", nfail)

    d = r % domain(1)
    call report(d % same_as(cells), &
         & "slot one answers the cells domain, by identity", nfail)
    call report(.not. d % same_as(faces), &
         & "and is not the faces domain, whatever the sizes say", nfail)
    d = r % domain(2)
    call report(d % same_as(faces), &
         & "slot two answers the faces domain", nfail)

  end subroutine check_incidence_contract

  !===================================================================!
  ! A relation is a set, not a bag: a tuple handed in twice is in
  ! the relation once, and everything the interface answers - the
  ! count, the table, membership - answers set semantics. Parallel
  ! edges do not live here; they are distinct members of an edge
  ! domain, as the ternary check shows.
  !===================================================================!

  subroutine check_set_semantics(nfail)

    integer, intent(inout) :: nfail

    type(set_graph)     :: cells, faces
    type(stored_relation) :: r
    integer, allocatable  :: t(:,:)
    type(set_map)     :: sets

    call cells % declare()
    call sets % bind(cells, counted_set_representation(4))
    call faces % declare()
    call sets % bind(faces, counted_set_representation(5))

    r = stored_relation('touches', [cells, faces], &
         & reshape([1,1,  2,2,  1,1,  3,3,  2,2], [2, 5]), sets)

    call report(r % num_tuples() .eq. 3, &
         & "a tuple handed in twice is in the relation once", nfail)

    call r % tuples(t)
    call report(size(t, 2) .eq. 3 .and. &
         &      all(t(:, 1) .eq. [1, 1]) .and. &
         &      all(t(:, 2) .eq. [2, 2]) .and. &
         &      all(t(:, 3) .eq. [3, 3]), &
         & "first appearances stand, in their first order", nfail)

    call report(r % has([1, 1]) .and. r % has([2, 2]), &
         & "and membership reads the set, not the handing", nfail)

  end subroutine check_set_semantics

  !===================================================================!
  ! Carriers obey set semantics too: a member handed in twice is in
  ! the domain once, enumeration is injective, and member and
  ! local_index invert each other - on the listed concretion, where
  ! members and positions genuinely differ.
  !===================================================================!

  subroutine check_carrier_set_semantics(nfail)

    integer, intent(inout) :: nfail

    type(set_graph)     :: dup
    integer, allocatable :: idx(:)
    integer              :: i, j
    logical              :: ok
    type(set_map)     :: sets

    call dup % declare()
    call sets % bind(dup, listed_set_representation([10, 20, 10, 30, 20]))

    call report(sets % size_of(dup) .eq. 3, &
         & "a member handed in twice is in the domain once", nfail)

    call sets % members_of(dup, idx)
    call report(all(idx .eq. [10, 20, 30]), &
         & "first appearances stand, in their first order", nfail)

    ok = .true.
    do i = 1, sets % size_of(dup)
       do j = i + 1, sets % size_of(dup)
          ok = ok .and. (sets % member_of(dup, i) /= sets % member_of(dup, j))
       end do
    end do
    call report(ok, &
         & "enumeration is injective: each member once", nfail)

    call report(sets % index_in(dup, 20) .eq. 2 .and. &
         &      sets % index_in(dup, 15) .eq. 0, &
         & "local_index finds the standing, zero for outsiders", nfail)

    ok = .true.
    do i = 1, sets % size_of(dup)
       ok = ok .and. (sets % member_of(dup, sets % index_in(dup, sets % member_of(dup, i))) &
            &         .eq. sets % member_of(dup, i))
       ok = ok .and. (sets % index_in(dup, sets % member_of(dup, i)) .eq. i)
    end do
    call report(ok, &
         & "member and local_index invert each other, both ways", nfail)

  end subroutine check_carrier_set_semantics

  !===================================================================!
  ! The signature is generic over the representation, proved: a
  ! counted set and a listed one stand in one signature, each slot
  ! judged by its own representation and answering its own identity.
  !===================================================================!

  subroutine check_mixed_carriers(nfail)

    integer, intent(inout) :: nfail

    type(set_graph)              :: cells
    type(set_graph)               :: sensors
    type(stored_relation)          :: r
    type(set_graph) :: d
    type(set_map)     :: sets

    call cells % declare()
    call sets % bind(cells, counted_set_representation(3))
    call sensors % declare()
    call sets % bind(sensors, listed_set_representation([10, 20, 30]))

    r = stored_relation('reads', [cells, sensors], &
         & reshape([1,10,  2,30], [2, 2]), sets)

    call report(r % arity() .eq. 2 .and. r % num_tuples() .eq. 2, &
         & "two carrier concretions stand in one signature", nfail)

    call report(r % has([2, 30]) .and. .not. r % has([2, 20]), &
         & "membership crosses the concretions untroubled", nfail)

    d = r % domain(2)
    call report(d % same_as(sensors), &
         & "the listed slot answers the listed domain, by identity", nfail)
    call report(sets % has_in(d, 20) .and. .not. sets % has_in(d, 15), &
         & "and the carrier's own membership law travels with it", nfail)

  end subroutine check_mixed_carriers

  !===================================================================!
  ! Adjacency is a reading, not a primitive: the same carrier stands
  ! in both slots, and the signature says so by identity.
  !===================================================================!

  subroutine check_adjacency_is_a_reading(nfail)

    integer, intent(inout) :: nfail

    type(set_graph)              :: cells
    type(stored_relation)          :: adj
    type(set_graph) :: d1, d2
    type(set_map)     :: sets

    call cells % declare()
    call sets % bind(cells, counted_set_representation(3))

    adj = stored_relation('beside', [cells, cells], &
         & reshape([1,2,  2,3], [2, 2]), sets)

    d1 = adj % domain(1)
    d2 = adj % domain(2)
    call report(d1 % same_as(d2), &
         & "one carrier may stand in both slots", nfail)
    call report(d1 % same_as(cells) .and. d2 % same_as(cells), &
         & "and both are the one declared domain", nfail)

  end subroutine check_adjacency_is_a_reading

  !===================================================================!
  ! The ternary endpoint relation of AGENTS.md 8.1: edges x vertices
  ! x roles. An interior edge holds two tuples, tail and head; the
  ! boundary edge holds one, and NO tuple invents its far side.
  !===================================================================!

  subroutine check_ternary_endpoint(nfail)

    integer, intent(inout) :: nfail

    type(set_graph)     :: edges, verts, roles
    type(stored_relation) :: ends
    integer               :: v
    logical               :: ok
    type(set_map)     :: sets

    call edges % declare()
    call sets % bind(edges, counted_set_representation(3))
    call verts % declare()
    call sets % bind(verts, counted_set_representation(4))
    call roles % declare()                 ! 1 = tail, 2 = head
    call sets % bind(roles, counted_set_representation(2))

    ! 1 --> 2 --> 3, then edge 3 leaves vertex 3 for the wall.
    ends = stored_relation('endpoint', [edges, verts, roles], &
         & reshape([1,1,1,  1,2,2,  2,2,1,  2,3,2,  3,3,1], [3, 5]), sets)

    call report(ends % arity() .eq. 3, &
         & "roles are a slot, not an attribute", nfail)
    call report(ends % has([1, 1, 1]) .and. ends % has([1, 2, 2]), &
         & "an interior edge holds its tail and its head", nfail)
    call report(ends % has([3, 3, 1]), &
         & "the boundary edge holds its one tail", nfail)

    ok = .true.
    do v = 1, 4
       ok = ok .and. .not. ends % has([3, v, 2])
    end do
    call report(ok, &
         & "and no tuple invents a head across the wall", nfail)

  end subroutine check_ternary_endpoint

  !===================================================================!
  ! Carriers are not owned by a graph, or by anything: two relations
  ! declared over the one cells domain agree, slot against slot,
  ! across relations.
  !===================================================================!

  subroutine check_shared_carriers(nfail)

    integer, intent(inout) :: nfail

    type(set_graph)              :: cells, faces
    type(stored_relation)          :: adj, inc
    type(set_graph) :: da, di
    type(set_map)     :: sets

    call cells % declare()
    call sets % bind(cells, counted_set_representation(4))
    call faces % declare()
    call sets % bind(faces, counted_set_representation(5))

    adj = stored_relation('beside' , [cells, cells], &
         & reshape([1,2,  3,4], [2, 2]), sets)
    inc = stored_relation('touches', [cells, faces], &
         & reshape([1,1,  4,5], [2, 2]), sets)

    da = adj % domain(1)
    di = inc % domain(1)
    call report(da % same_as(di), &
         & "two relations share one carrier, by identity", nfail)

  end subroutine check_shared_carriers

  !===================================================================!
  ! Relation identity is independent of signature: two relations
  ! over the same slots coexist, distinct, both answering.
  !===================================================================!

  subroutine check_relation_identity(nfail)

    integer, intent(inout) :: nfail

    type(set_graph)     :: cells, faces
    type(stored_relation) :: physical, coarse, copy
    type(set_map)     :: sets

    call cells % declare()
    call sets % bind(cells, counted_set_representation(4))
    call faces % declare()
    call sets % bind(faces, counted_set_representation(5))

    physical = stored_relation('physical', [cells, faces], &
         & reshape([1,1,  2,2], [2, 2]), sets)
    coarse   = stored_relation('coarse'  , [cells, faces], &
         & reshape([1,1,  3,3], [2, 2]), sets)

    call report(.not. physical % same_as(coarse), &
         & "same signature, two relations - no collision", nfail)
    call report(physical % has([2, 2]) .and. .not. coarse % has([2, 2]), &
         & "and each keeps its own tuples", nfail)

    copy = physical
    call report(copy % same_as(physical), &
         & "a copy is the same declared relation", nfail)

  end subroutine check_relation_identity

  !===================================================================!
  ! The empty relation is a relation: declared, counted at zero,
  ! holding nothing.
  !===================================================================!

  subroutine check_empty_relation(nfail)

    integer, intent(inout) :: nfail

    type(set_graph)     :: cells
    type(stored_relation) :: none
    type(token)           :: tk
    type(set_map)     :: sets

    call cells % declare()
    call sets % bind(cells, counted_set_representation(3))

    none = stored_relation('nothing', [cells, cells], &
         & reshape([integer ::], [2, 0]), sets)

    tk = none % id()
    call report(tk % declared() .and. none % num_tuples() .eq. 0, &
         & "the empty relation is declared and holds nothing", nfail)
    call report(.not. none % has([1, 1]), &
         & "and no tuple is in it", nfail)

  end subroutine check_empty_relation

end program test_graph_relation
