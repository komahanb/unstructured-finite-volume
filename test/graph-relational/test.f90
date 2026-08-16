!=====================================================================!
! THE RELATIONAL VIEW SUITE
!
! A graph read as (S, P):
!
!     branch(1) = the sequence of member sets
!     branch(2) = the sequence of relations
!
! Each element graph denotes a legacy member set or relation, resolved
! through a binding. The kernel is unchanged; graph_structure is not
! used here at all.
!
! THE VALIDITY LAW, stated mathematically:
!
!     G is relationally valid  iff  for every relation P in the
!     relation sequence and every j in 1..arity(P), domain_j(P) is a
!     member set of the member-set sequence of G.
!
! Zero sets and zero relations satisfy it vacuously. Relations without
! sets do not, unless they have arity zero - and a relation of this
! repository has arity at least one, so a relation always needs its
! domains present.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test

  use fractal_graph        , only : graph, null_branch, known_branch
  use graph_carrier        , only : member_set, counted_set
  use graph_relation       , only : relation, stored_relation, slot
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
    type(counted_set)          :: e, v
    class(member_set), pointer :: s
    integer                    :: i

    call g % declare()
    do i = 1, 2
       call cell(i) % declare(); call elem(i) % declare()
    end do

    e = counted_set('E', 3)
    v = counted_set('V', 4)

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
    type(counted_set)          :: e, v
    type(stored_relation)      :: t
    class(relation), pointer   :: r
    integer                    :: i

    call g % declare()
    do i = 1, 2
       call scell(i) % declare(); call selem(i) % declare()
    end do
    call rcell(1) % declare(); call relem(1) % declare()

    e = counted_set('E', 2)
    v = counted_set('V', 3)
    t = stored_relation('T', [slot(e), slot(v)], &
         & reshape([1, 1, 2, 2], [2, 2]))

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
    type(counted_set)        :: e, v, foreign
    type(stored_relation)    :: t

    call g % declare(); call scell % declare(); call selem % declare()
    call rcell % declare(); call relem % declare()

    e       = counted_set('E', 2)
    v       = counted_set('V', 3)
    foreign = counted_set('W', 3)

    t = stored_relation('T', [slot(e), slot(foreign)], &
         & reshape([1, 1, 2, 2], [2, 2]))

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
    type(counted_set)        :: e
    type(stored_relation)    :: p

    call g % declare(); call rcell % declare(); call relem % declare()

    e = counted_set('E', 2)
    p = stored_relation('P', [slot(e)], reshape([1, 2], [1, 2]))

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
  ! Identity. A view has none of its own.
  !===================================================================!

  identity_block: block

    type(graph), target :: g, h

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
