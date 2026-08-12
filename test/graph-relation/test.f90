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
! over member_set - proved by mixing a counted carrier with the
! listed fixture in one relation.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_graph_relation

  use graph_identity    , only : token
  use graph_carrier     , only : counted_set, member_set
  use graph_relation    , only : stored_relation, slot
  use listed_set_fixture, only : listed_set

  implicit none

  integer :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "graph relation suite (AGENTS phase 2)"
  write(*,'(1x,a)') "============================================="

  call check_incidence_contract(nfail)
  call check_set_semantics(nfail)
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

    type(counted_set)              :: cells, faces
    type(stored_relation)          :: r
    class(member_set), allocatable :: d
    type(token)                    :: tk
    integer, allocatable           :: t(:,:)

    cells = counted_set('cells', 4)
    faces = counted_set('faces', 5)

    r = stored_relation('touches', [cells, faces], &
         & reshape([1,1,  1,2,  2,2,  3,4], [2, 4]))

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

    type(counted_set)     :: cells, faces
    type(stored_relation) :: r
    integer, allocatable  :: t(:,:)

    cells = counted_set('cells', 4)
    faces = counted_set('faces', 5)

    r = stored_relation('touches', [cells, faces], &
         & reshape([1,1,  2,2,  1,1,  3,3,  2,2], [2, 5]))

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
  ! The signature is generic over member_set, proved: a counted set
  ! and a listed fixture stand in one signature, each slot judging
  ! membership by its own law and answering its own identity.
  !===================================================================!

  subroutine check_mixed_carriers(nfail)

    integer, intent(inout) :: nfail

    type(counted_set)              :: cells
    type(listed_set)               :: sensors
    type(stored_relation)          :: r
    class(member_set), allocatable :: d

    cells   = counted_set('cells', 3)
    sensors = listed_set('sensors', [10, 20, 30])

    r = stored_relation('reads', [slot(cells), slot(sensors)], &
         & reshape([1,10,  2,30], [2, 2]))

    call report(r % arity() .eq. 2 .and. r % num_tuples() .eq. 2, &
         & "two carrier concretions stand in one signature", nfail)

    call report(r % has([2, 30]) .and. .not. r % has([2, 20]), &
         & "membership crosses the concretions untroubled", nfail)

    d = r % domain(2)
    call report(d % same_as(sensors), &
         & "the listed slot answers the listed domain, by identity", nfail)
    call report(d % has(20) .and. .not. d % has(15), &
         & "and the carrier's own membership law travels with it", nfail)

  end subroutine check_mixed_carriers

  !===================================================================!
  ! Adjacency is a reading, not a primitive: the same carrier stands
  ! in both slots, and the signature says so by identity.
  !===================================================================!

  subroutine check_adjacency_is_a_reading(nfail)

    integer, intent(inout) :: nfail

    type(counted_set)              :: cells
    type(stored_relation)          :: adj
    class(member_set), allocatable :: d1, d2

    cells = counted_set('cells', 3)

    adj = stored_relation('beside', [cells, cells], &
         & reshape([1,2,  2,3], [2, 2]))

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

    type(counted_set)     :: edges, verts, roles
    type(stored_relation) :: ends
    integer               :: v
    logical               :: ok

    edges = counted_set('edges'   , 3)
    verts = counted_set('vertices', 4)
    roles = counted_set('roles'   , 2)     ! 1 = tail, 2 = head

    ! 1 --> 2 --> 3, then edge 3 leaves vertex 3 for the wall.
    ends = stored_relation('endpoint', [edges, verts, roles], &
         & reshape([1,1,1,  1,2,2,  2,2,1,  2,3,2,  3,3,1], [3, 5]))

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

    type(counted_set)              :: cells, faces
    type(stored_relation)          :: adj, inc
    class(member_set), allocatable :: da, di

    cells = counted_set('cells', 4)
    faces = counted_set('faces', 5)

    adj = stored_relation('beside' , [cells, cells], &
         & reshape([1,2,  3,4], [2, 2]))
    inc = stored_relation('touches', [cells, faces], &
         & reshape([1,1,  4,5], [2, 2]))

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

    type(counted_set)     :: cells, faces
    type(stored_relation) :: physical, coarse, copy

    cells = counted_set('cells', 4)
    faces = counted_set('faces', 5)

    physical = stored_relation('physical', [cells, faces], &
         & reshape([1,1,  2,2], [2, 2]))
    coarse   = stored_relation('coarse'  , [cells, faces], &
         & reshape([1,1,  3,3], [2, 2]))

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

    type(counted_set)     :: cells
    type(stored_relation) :: none
    type(token)           :: tk

    cells = counted_set('cells', 3)

    none = stored_relation('nothing', [cells, cells], &
         & reshape([integer ::], [2, 0]))

    tk = none % id()
    call report(tk % declared() .and. none % num_tuples() .eq. 0, &
         & "the empty relation is declared and holds nothing", nfail)
    call report(.not. none % has([1, 1]), &
         & "and no tuple is in it", nfail)

  end subroutine check_empty_relation

end program test_graph_relation
