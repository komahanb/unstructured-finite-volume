!=====================================================================!
! THE RELATIONAL VIEW SUITE
!
! A graph read as (S, P):
!
!     branch(1) = the sequence of member sets
!     branch(2) = the sequence of relations
!
! Each element graph denotes a legacy member set or relation, resolved
! through a binding. The kernel is unchanged.
!
! THE VALIDITY LAW, stated mathematically. G is relationally valid iff
!
!     (i)   S_i /= S_j for i /= j,
!     (ii)  P_i /= P_j for i /= j, and
!     (iii) for every relation P and every j in 1..arity(P),
!           domain_j(P) is a member set of the member-set sequence.
!
! (i) and (ii) say that S and P are SETS: the branches represent them
! as sequences, and a sequence may repeat what a set cannot.
!
! Zero sets and zero relations satisfy all three vacuously. Relations
! without sets do not, unless they have arity zero - and a relation of
! this repository has arity at least one, so a relation always needs
! its domains present.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test

  use fractal_graph        , only : graph, null_branch, known_branch
  use fractal_graph        , only : set_graph => graph
  use graph_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use graph_set_map        , only : set_map
  use graph_relation       , only : relation, stored_relation
  use graph_binary_relation, only : csr_relation
  use graph_relational_view, only : relational_binding, &
       & num_member_sets, member_set_at, num_relations, relation_at, &
       & holds_set, relational_valid

  implicit none

  integer :: failures = 0

  write(*,'(1x,a)') "relational view suite"

  !===================================================================!
  ! Zero sets, zero relations. Both sequences empty: NULL on both
  ! branches. Vacuously valid.
  !===================================================================!

  vacuous_block: block

    type(graph), target      :: g
    type(relational_binding) :: b

    call g % declare()
    g % branch(1) = null_branch()
    g % branch(2) = null_branch()

    call check('zero sets, zero relations: |S| = 0, |P| = 0', &
         & num_member_sets(g) .eq. 0 .and. num_relations(g) .eq. 0)
    call check('and it is vacuously valid', relational_valid(g, b))

  end block vacuous_block

  !===================================================================!
  ! Sets only. Two member sets, no relations.
  !===================================================================!

  sets_only_block: block

    type(graph), target        :: g, cell(2), elem(2)
    type(relational_binding)   :: b
    type(set_graph)          :: e, v
    type(set_graph), pointer :: s
    integer                    :: i
    type(set_map)     :: sets

    call g % declare()
    do i = 1, 2
       call cell(i) % declare(); call elem(i) % declare()
    end do

    call e % declare()
    call sets % bind(e, counted_set_representation(3))
    call v % declare()
    call sets % bind(v, counted_set_representation(4))

    call b % bind_set(elem(1), e)
    call b % bind_set(elem(2), v)

    cell(1) % branch(1) = known_branch(elem(1))
    cell(1) % branch(2) = known_branch(cell(2))
    cell(2) % branch(1) = known_branch(elem(2))
    cell(2) % branch(2) = null_branch()

    g % branch(1) = known_branch(cell(1))
    g % branch(2) = null_branch()

    call check('sets only: |S| = 2, |P| = 0', &
         & num_member_sets(g) .eq. 2 .and. num_relations(g) .eq. 0)

    s => member_set_at(g, b, 1)
    call check('member_set_at(1) is E', s % same_as(e))
    s => member_set_at(g, b, 2)
    call check('member_set_at(2) is V', s % same_as(v))

    call check('holds_set answers for both, and refuses none', &
         & holds_set(g, b, e) .and. holds_set(g, b, v))
    call check('sets only is vacuously valid', relational_valid(g, b))

  end block sets_only_block

  !===================================================================!
  ! Valid relation domains. Two sets, one binary relation over them.
  !===================================================================!

  valid_block: block

    type(graph), target        :: g
    type(graph), target        :: scell(2), selem(2), rcell(1), relem(1)
    type(relational_binding)   :: b
    type(set_graph)          :: e, v
    type(stored_relation)      :: t
    class(relation), pointer   :: r
    integer                    :: i
    type(set_map)     :: sets

    call g % declare()
    do i = 1, 2
       call scell(i) % declare(); call selem(i) % declare()
    end do
    call rcell(1) % declare(); call relem(1) % declare()

    call e % declare()
    call sets % bind(e, counted_set_representation(2))
    call v % declare()
    call sets % bind(v, counted_set_representation(3))
    t = stored_relation('T', [e, v], &
         & reshape([1, 1, 2, 2], [2, 2]), sets)

    call b % bind_set(selem(1), e)
    call b % bind_set(selem(2), v)
    call b % bind_relation(relem(1), t)

    scell(1) % branch(1) = known_branch(selem(1))
    scell(1) % branch(2) = known_branch(scell(2))
    scell(2) % branch(1) = known_branch(selem(2))
    scell(2) % branch(2) = null_branch()

    rcell(1) % branch(1) = known_branch(relem(1))
    rcell(1) % branch(2) = null_branch()

    g % branch(1) = known_branch(scell(1))
    g % branch(2) = known_branch(rcell(1))

    call check('valid: |S| = 2, |P| = 1', &
         & num_member_sets(g) .eq. 2 .and. num_relations(g) .eq. 1)

    r => relation_at(g, b, 1)
    call check('relation_at(1) is T, by identity', r % same_as(t))
    call check('and it is a reference into owned storage, stable across calls', &
         & associated(r, relation_at(g, b, 1)))

    call check('every domain of T is a member set of G: valid', &
         & relational_valid(g, b))

  end block valid_block

  !===================================================================!
  ! A relation naming a set the graph does not hold. The view reports
  ! invalidity; it does not refuse. Refusal belongs to a malformed
  ! representation, which this is not.
  !===================================================================!

  foreign_block: block

    type(graph), target      :: g, scell, selem, rcell, relem
    type(relational_binding) :: b
    type(set_graph)        :: e, v, foreign
    type(stored_relation)    :: t
    type(set_map)     :: sets

    call g % declare(); call scell % declare(); call selem % declare()
    call rcell % declare(); call relem % declare()

    call e % declare()
    call sets % bind(e, counted_set_representation(2))
    call v % declare()
    call sets % bind(v, counted_set_representation(3))
    call foreign % declare()
    call sets % bind(foreign, counted_set_representation(3))

    t = stored_relation('T', [e, foreign], &
         & reshape([1, 1, 2, 2], [2, 2]), sets)

    call b % bind_set(selem, e)                  ! only E is held
    call b % bind_relation(relem, t)

    scell % branch(1) = known_branch(selem)
    scell % branch(2) = null_branch()
    rcell % branch(1) = known_branch(relem)
    rcell % branch(2) = null_branch()

    g % branch(1) = known_branch(scell)
    g % branch(2) = known_branch(rcell)

    call check('foreign domain: the graph is well formed', &
         & num_member_sets(g) .eq. 1 .and. num_relations(g) .eq. 1)
    call check('the graph does not hold W', .not. holds_set(g, b, foreign))
    call check('so it is relationally INVALID, and reported, not refused', &
         & .not. relational_valid(g, b))
    call check('E is still held, and V was never bound', &
         & holds_set(g, b, e) .and. .not. holds_set(g, b, v))

  end block foreign_block

  !===================================================================!
  ! Relations with no sets. Admissible as a representation, invalid as
  ! a relational graph: a relation of arity >= 1 needs its domains.
  !===================================================================!

  relations_only_block: block

    type(graph), target      :: g, rcell, relem
    type(relational_binding) :: b
    type(set_graph)        :: e
    type(stored_relation)    :: p
    type(set_map)     :: sets

    call g % declare(); call rcell % declare(); call relem % declare()

    call e % declare()
    call sets % bind(e, counted_set_representation(2))
    p = stored_relation('P', [e], reshape([1, 2], [1, 2]), sets)

    call b % bind_relation(relem, p)

    rcell % branch(1) = known_branch(relem)
    rcell % branch(2) = null_branch()

    g % branch(1) = null_branch()                ! no member sets
    g % branch(2) = known_branch(rcell)

    call check('relations only: |S| = 0, |P| = 1 is representable', &
         & num_member_sets(g) .eq. 0 .and. num_relations(g) .eq. 1)
    call check('but a relation of arity 1 needs its domain: invalid', &
         & .not. relational_valid(g, b))

  end block relations_only_block

  !===================================================================!
  ! S and P are sets. A sequence may repeat what a set cannot, so
  ! repetition is representable and invalid - answered, not refused.
  ! One member set, denoted by two distinct element graphs.
  !===================================================================!

  repeated_set_block: block

    type(graph), target      :: g, scell(2), selem(2)
    type(relational_binding) :: b
    type(set_graph)        :: e
    integer                  :: i
    type(set_map)     :: sets

    call g % declare()
    do i = 1, 2
       call scell(i) % declare(); call selem(i) % declare()
    end do

    call e % declare()
    call sets % bind(e, counted_set_representation(2))

    call b % bind_set(selem(1), e)
    call b % bind_set(selem(2), e)               ! the same member set

    scell(1) % branch(1) = known_branch(selem(1))
    scell(1) % branch(2) = known_branch(scell(2))
    scell(2) % branch(1) = known_branch(selem(2))

    g % branch(1) = known_branch(scell(1))
    g % branch(2) = null_branch()

    call check('a repeated member set is representable: |S| = 2', &
         & num_member_sets(g) .eq. 2)
    call check('but S is a set, so the graph is INVALID', &
         & .not. relational_valid(g, b))

  end block repeated_set_block

  repeated_relation_block: block

    type(graph), target      :: g, scell, selem, rcell(2), relem(2)
    type(relational_binding) :: b
    type(set_graph)        :: e
    type(stored_relation)    :: p
    integer                  :: i
    type(set_map)     :: sets

    call g % declare(); call scell % declare(); call selem % declare()
    do i = 1, 2
       call rcell(i) % declare(); call relem(i) % declare()
    end do

    call e % declare()
    call sets % bind(e, counted_set_representation(2))
    p = stored_relation('P', [e], reshape([1, 2], [1, 2]), sets)

    call b % bind_set(selem, e)
    call b % bind_relation(relem(1), p)
    call b % bind_relation(relem(2), p)          ! the same relation

    scell % branch(1) = known_branch(selem)
    rcell(1) % branch(1) = known_branch(relem(1))
    rcell(1) % branch(2) = known_branch(rcell(2))
    rcell(2) % branch(1) = known_branch(relem(2))

    g % branch(1) = known_branch(scell)
    g % branch(2) = known_branch(rcell(1))

    call check('a repeated relation is representable: |P| = 2', &
         & num_relations(g) .eq. 2)
    call check('but P is a set, so the graph is INVALID', &
         & .not. relational_valid(g, b))

  end block repeated_relation_block

  !===================================================================!
  ! Identity is the address, never the signature: two relations over
  ! the same slots coexist, each keeping its own tuples. A ternary
  ! relation sits beside a binary one, untroubled.
  !===================================================================!

  signature_block: block

    type(graph), target      :: g, scell(2), selem(2), rcell(3), relem(3)
    type(relational_binding) :: b
    type(set_graph)        :: ops, vals
    type(csr_relation)       :: physical, scheduled
    type(stored_relation)    :: flow
    class(relation), pointer :: r1, r2, r3
    integer                  :: i
    type(set_map)     :: sets

    call g % declare()
    do i = 1, 2
       call scell(i) % declare(); call selem(i) % declare()
    end do
    do i = 1, 3
       call rcell(i) % declare(); call relem(i) % declare()
    end do

    call ops % declare()
    call sets % bind(ops, counted_set_representation(3))
    call vals % declare()
    call sets % bind(vals, counted_set_representation(5))

    physical  = csr_relation('feeds' , ops, ops, reshape([1, 2], [2, 1]), sets)
    scheduled = csr_relation('awaits', ops, ops, reshape([2, 3], [2, 1]), sets)
    flow      = stored_relation('flow', [ops, vals, ops], &
         & reshape([1, 2, 1,  2, 4, 3], [3, 2]), sets)

    call b % bind_set(selem(1), ops)
    call b % bind_set(selem(2), vals)
    call b % bind_relation(relem(1), physical)
    call b % bind_relation(relem(2), scheduled)
    call b % bind_relation(relem(3), flow)

    scell(1) % branch(1) = known_branch(selem(1))
    scell(1) % branch(2) = known_branch(scell(2))
    scell(2) % branch(1) = known_branch(selem(2))

    do i = 1, 3
       rcell(i) % branch(1) = known_branch(relem(i))
       if (i .lt. 3) rcell(i) % branch(2) = known_branch(rcell(i + 1))
    end do

    g % branch(1) = known_branch(scell(1))
    g % branch(2) = known_branch(rcell(1))

    r1 => relation_at(g, b, 1)
    r2 => relation_at(g, b, 2)
    r3 => relation_at(g, b, 3)

    call check('same signature, two relations - no collision', &
         & .not. r1 % same_as(r2))
    call check('and each keeps its own tuples', &
         & r1 % has([1, 2]) .and. .not. r2 % has([1, 2]))
    call check('a ternary relation sits beside two binary ones', &
         & r3 % arity() .eq. 3 .and. r1 % arity() .eq. 2)
    call check('and the whole is valid: every domain is held', &
         & relational_valid(g, b))

  end block signature_block

  !===================================================================!
  ! Identity. A view has none of its own.
  !===================================================================!

  identity_block: block

    type(graph), target :: g, h
    type(set_map)     :: sets

    call g % declare(); call h % declare()
    g % branch(1) = null_branch(); g % branch(2) = null_branch()
    h % branch(1) = null_branch(); h % branch(2) = null_branch()

    call check('the view mints no identity: same_as is the graph''s', &
         & g % same_as(g) .and. .not. g % same_as(h))

  end block identity_block

  !===================================================================!

  if (failures .eq. 0) then
     print *, ''
     print *, ' ALL PROPOSITIONS HOLD'
  else
     print *, ''
     print *, ' FAILURES :', failures
     error stop 'test: a proposition failed'
  end if

contains

  subroutine check(label, ok)

    character(len=*), intent(in) :: label
    logical         , intent(in) :: ok

    if (ok) then
       print *, ' PASS : ', label
    else
       print *, ' FAIL : ', label
       failures = failures + 1
    end if

  end subroutine check

end program test
