!=====================================================================!
! TIME INTEGRATION TOWER . LEVEL 3 . RELATIONAL GRAPH
!
! The level answers one question: CAN TIME STRUCTURE AND THE STATE
! DOMAIN COEXIST AS ONE RELATIONAL STRUCTURE WITHOUT BEING
! CONFLATED.
!
!      G_time = ( { Q, T, E }, { Tail, Head, A1, A2 } )
!
! Three member sets, four relations, one container. The signature
! validity law is the price of admission: every slot of every owned
! relation must answer one of the graph's own carriers, or the
! graph refuses to exist.
!
!                    THE LEVEL'S OWN TRUTH
!
! Q IS OWNED, AND NO RELATION NAMES IT.
!
! That is lawful, and the law is worth stating precisely because it
! is asymmetric: the container validates RELATIONS against its SETS,
! and never the reverse. It demands that every relation's slots
! resolve to an owned carrier. It does not demand that every owned
! carrier be named by a relation - and it should not, because a
! relational graph is a collection of member sets and typed
! relations over them, not a connected object.
!
! So the state axis sits inside the same structure as the time axis,
! answering nothing, waiting for Level 5 to put a field on it. The
! alternative - inventing an incidence merely to attach Q to the
! chain - would manufacture exactly the conflation this tower
! exists to refuse, and would do it in the name of tidiness.
!
!      Q is unrelated here.  Unrelated is not absent.
!
! A1 and A2 are materialized once, by the Level-2 derivation, and
! admitted. The derivation itself lives one rung down; this level
! does not re-derive it as a second source of truth, and no
! provenance type is required to say so.
!
! No ordinary graph profile yet. No values. No marching.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program time_level_3

  use time_assert           , only : report, verdict
  use time_assert           , only : T0, T1, T2, E1
  use graph_carrier         , only : counted_set, member_set
  use graph_relation        , only : relation
  use graph_binary_relation , only : csr_relation
  use fractal_graph        , only : graph, known_branch, null_branch
  use graph_relational_view, only : relational_binding, &
       & num_member_sets, member_set_at, num_relations, relation_at, &
       & holds_set
  use time_carriers_fixture , only : time_carriers
  use time_relations_fixture, only : tail_relation, head_relation
  use time_algebra_fixture  , only : derive_one_step_reach, &
       &                             derive_two_step_reach

  implicit none

  type(counted_set)              :: q, t, e
  type(csr_relation), target     :: tail, head, a1
  type(csr_relation)             :: a2
  type(graph)             , target :: g
  type(graph)             , target :: scell(3), selem(3)
  type(graph)             , target :: rcell(4), relem(4)
  type(relational_binding)         :: bnd
  integer                          :: k
  integer                        :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "time integration tower . level 3 . graph"
  write(*,'(1x,a)') "============================================="

  call time_carriers(q, t, e)
  tail = tail_relation(e, t)
  head = head_relation(e, t)
  a1   = derive_one_step_reach(tail, head)
  a2   = derive_two_step_reach(a1)

  ! 'time': (S, P) as one sequence on each branch.
  call g % declare()
  do k = 1, 3
     call scell(k) % declare()
     call selem(k) % declare()
  end do
  do k = 1, 4
     call rcell(k) % declare()
     call relem(k) % declare()
  end do

  call bnd % bind_set(selem(1), q)
  call bnd % bind_set(selem(2), t)
  call bnd % bind_set(selem(3), e)
  call bnd % bind_relation(relem(1), tail)
  call bnd % bind_relation(relem(2), head)
  call bnd % bind_relation(relem(3), a1)
  call bnd % bind_relation(relem(4), a2)

  do k = 1, 3
     scell(k) % branch(1) = known_branch(selem(k))
     if (k .lt. 3) scell(k) % branch(2) = &
          & known_branch(scell(k + 1))
  end do
  do k = 1, 4
     rcell(k) % branch(1) = known_branch(relem(k))
     if (k .lt. 4) rcell(k) % branch(2) = &
          & known_branch(rcell(k + 1))
  end do

  g % branch(1) = known_branch(scell(1))
  g % branch(2) = known_branch(rcell(1))

  call check_ownership(nfail)
  call check_signature_closure(nfail)
  call check_signatures_survive_ownership(nfail)
  call check_state_axis_is_owned_and_unrelated(nfail)

  call verdict(nfail, "level 3")

contains

  !===================================================================!
  ! The graph owns exactly the three carriers and the four
  ! relations it was handed - by structural identity, through the
  ! generators alone.
  !===================================================================!

  subroutine check_ownership(nfail)

    integer, intent(inout) :: nfail

    call report(num_member_sets(g) .eq. 3, &
         & "G_time owns three member sets", nfail)
    call report(num_relations(g) .eq. 4, &
         & "and four relations", nfail)

    call report(holds_set(g, bnd, q) .and. holds_set(g, bnd, t) .and. &
         &      holds_set(g, bnd, e), &
         & "it holds Q, T and E - both axes and the steps between " // &
         & "them, in one structure", nfail)

    call report(owns_relation(tail) .and. owns_relation(head) .and. &
         &      owns_relation(a1) .and. owns_relation(a2), &
         & "and it holds Tail, Head, A1 and A2, each by identity", &
         & nfail)

  end subroutine check_ownership

  !===================================================================!
  ! THE admission law: every slot of every owned relation resolves
  ! to a member set the graph owns. Checked here from outside,
  ! against the graph's own accessors, rather than trusted because
  ! construction did not stop.
  !===================================================================!

  subroutine check_signature_closure(nfail)

    integer, intent(inout) :: nfail

    class(relation)  , pointer     :: r
    class(member_set), allocatable :: d
    integer                        :: k, slot, s
    logical                        :: ok, found

    ok = .true.
    do k = 1, num_relations(g)
       r => relation_at(g, bnd, k)
       do slot = 1, r % arity()
          d = r % domain(slot)
          found = .false.
          do s = 1, num_member_sets(g)
             found = found .or. d % same_as(member_set_at(g, bnd, s))
          end do
          ok = ok .and. found
       end do
    end do

    call report(ok, &
         & "every slot of every owned relation answers a carrier " // &
         & "G_time owns: the signature validity law holds", nfail)

  end subroutine check_signature_closure

  !===================================================================!
  ! Ownership changes nothing about what a relation IS. Tail and
  ! Head remain E x T; A1 and A2 remain T x T. A container that
  ! quietly renormalized its contents would be a different object.
  !===================================================================!

  subroutine check_signatures_survive_ownership(nfail)

    integer, intent(inout) :: nfail

    call report(owned_runs(tail, e, t) .and. owned_runs(head, e, t), &
         & "Tail and Head are still E x T after ownership", nfail)

    call report(owned_runs(a1, t, t) .and. owned_runs(a2, t, t), &
         & "and A1 and A2 are still T x T", nfail)

    call report(holds_fact(tail, [E1, T0]) .and. &
         &      holds_fact(a1, [T0, T1]) .and. &
         &      holds_fact(a2, [T0, T2]), &
         & "and their extensions came through whole - read back " // &
         & "from graph-owned storage, not from the originals", nfail)

  end subroutine check_signatures_survive_ownership

  !===================================================================!
  ! THE rung's truth: Q is owned, and nothing names it. Proved
  ! positively - no slot of any owned relation resolves to Q - so
  ! the claim is a measurement rather than an absence of evidence.
  !===================================================================!

  subroutine check_state_axis_is_owned_and_unrelated(nfail)

    integer, intent(inout) :: nfail

    class(relation)  , pointer     :: r
    class(member_set), allocatable :: d
    integer                        :: k, slot
    logical                        :: named

    named = .false.
    do k = 1, num_relations(g)
       r => relation_at(g, bnd, k)
       do slot = 1, r % arity()
          d = r % domain(slot)
          named = named .or. d % same_as(q)
       end do
    end do

    call report(holds_set(g, bnd, q) .and. .not. named, &
         & "Q IS OWNED AND NO RELATION NAMES IT - the container " // &
         & "validates relations against sets, never sets against " // &
         & "relations", nfail)

    call report(.not. named .and. num_relations(g) .eq. 4, &
         & "and no relation was invented to attach it: a relational " // &
         & "structure need not be connected", nfail)

    call report(holds_set(g, bnd, t) .and. .not. q % same_as(t) .and. &
         &      .not. q % same_as(e), &
         & "both axes live in one structure, and neither has become " // &
         & "the other", nfail)

  end subroutine check_state_axis_is_owned_and_unrelated

  !===================================================================!
  ! Helpers that ask the GRAPH, never the local variables.
  !===================================================================!

  logical function owns_relation(r)

    class(relation), intent(in) :: r

    class(relation), pointer :: held
    integer                  :: k

    owns_relation = .false.
    do k = 1, num_relations(g)
       held => relation_at(g, bnd, k)
       owns_relation = owns_relation .or. held % same_as(r)
    end do

  end function owns_relation

  logical function owned_runs(selector, from, into)

    class(relation)  , intent(in) :: selector
    class(member_set), intent(in) :: from, into

    class(relation)  , pointer     :: held
    class(member_set), allocatable :: d
    integer                        :: k

    owned_runs = .false.
    do k = 1, num_relations(g)
       held => relation_at(g, bnd, k)
       if (held % same_as(selector)) then
          d = held % domain(1)
          owned_runs = d % same_as(from)
          d = held % domain(2)
          owned_runs = owned_runs .and. d % same_as(into)
          owned_runs = owned_runs .and. (held % arity() .eq. 2)
       end if
    end do

  end function owned_runs

  logical function holds_fact(selector, tuple)

    class(relation), intent(in) :: selector
    integer        , intent(in) :: tuple(:)

    class(relation), pointer :: held
    integer                  :: k

    holds_fact = .false.
    do k = 1, num_relations(g)
       held => relation_at(g, bnd, k)
       if (held % same_as(selector)) holds_fact = held % has(tuple)
    end do

  end function holds_fact

end program time_level_3
