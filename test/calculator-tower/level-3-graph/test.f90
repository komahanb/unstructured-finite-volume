!=====================================================================!
! CALCULATOR TOWER . LEVEL 3 . THE RELATIONAL GRAPH
!
! The level answers one question: HOW CARRIERS AND RELATIONS COEXIST
! AS ONE STRUCTURE. The persistent calculator becomes
!
!      GAMMA = ( { X, O, P }, { T_flow, D } )
!
! where D is ADMITTED, not redefined: it is derived here once more
! by the approved Level-2 road - restrict, project, compose - and
! the graph receives the materialized result. The pair (+, x) is
! never written on the way in. No provenance machinery rides the
! graph; the derivation stands in this file, as the level contract
! asks.
!
! The strongest negative truth of the level is the import list: the
! calculator compiles as a relational graph WITHOUT the ordinary
! profile - no vertex, no edge, no tail, no head anywhere.
!
! A second, smaller graph witnesses the coexistence law: D and its
! own restriction to the times-only subset - the empty relation
! over O x O - share one signature and remain two citizens.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program calculator_level_3

  use calculator_assert, only : report, verdict
  use calculator_assert, only : SLOT_A, SLOT_B, SLOT_C, SLOT_D, SLOT_E
  use calculator_assert, only : OP_PLUS, OP_TIMES
  use calculator_assert, only : PORT_IN1, PORT_IN2, PORT_OUT
  use graph_fractal        , only : graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set        , only : set_map
  use map_inclusion  , only : inclusion_map, declared_subobject
  use relation_finitary   , only : stored_relation, relation
  use relation_algebra, only : restrict_slot, project_slots, &
       &                             compose_binary
  use graph_fractal        , only : graph, known_branch, null_branch
  use view_relational, only : relational_binding, &
       & num_member_sets, member_set_at, num_relations, relation_at, &
       & has_set, relational_valid

  implicit none

  type(graph)              :: x, o, p
  type(graph)               :: p_out, p_in
  type(stored_relation)          :: flow, t_out3, t_in3
  type(stored_relation)          :: produces, consumes
  class(relation), allocatable   :: d
  type(graph)             , target :: g
  type(graph)             , target :: scell(3), selem(3)
  type(graph)             , target :: rcell(2), relem(2)
  type(relational_binding)         :: bnd
  integer                          :: k
  integer                        :: table(3, 6)
  integer                        :: nfail
  type(set_map)     :: sets
  type(inclusion_map)     :: inclusions

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "calculator tower . level 3 . relational graph"
  write(*,'(1x,a)') "============================================="

  call x % declare()
  call sets % bind(x, counted_set_representation(5))
  call o % declare()
  call sets % bind(o, counted_set_representation(2))
  call p % declare()
  call sets % bind(p, counted_set_representation(3))

  table(:, 1) = [OP_PLUS , SLOT_A, PORT_IN1]
  table(:, 2) = [OP_PLUS , SLOT_B, PORT_IN2]
  table(:, 3) = [OP_PLUS , SLOT_C, PORT_OUT]
  table(:, 4) = [OP_TIMES, SLOT_C, PORT_IN1]
  table(:, 5) = [OP_TIMES, SLOT_D, PORT_IN2]
  table(:, 6) = [OP_TIMES, SLOT_E, PORT_OUT]

  flow = stored_relation('flow', [o, x, p], table, sets)

  ! The approved Level-2 road, walked once more; the graph admits
  ! what the algebra derived.
  call p_out % declare()
  call sets       % bind(p_out, listed_set_representation([PORT_OUT]))
  call inclusions % include_in(p_out, p)
  call p_in % declare()
  call sets       % bind(p_in, listed_set_representation([PORT_IN1, PORT_IN2]))
  call inclusions % include_in(p_in, p)
  t_out3   = restrict_slot(flow, 3, p_out, sets, inclusions)
  t_in3    = restrict_slot(flow, 3, p_in, sets, inclusions)
  produces = project_slots(t_out3, [1, 2], sets)
  consumes = project_slots(t_in3 , [2, 1], sets)
  d        = compose_binary(produces, consumes, sets)

  ! 'calculator': (S, P) as one sequence on each branch.
  call g % declare()
  do k = 1, 3
     call scell(k) % declare()
     call selem(k) % declare()
  end do
  do k = 1, 2
     call rcell(k) % declare()
     call relem(k) % declare()
  end do

  call bnd % bind_set(selem(1), x)
  call bnd % bind_set(selem(2), o)
  call bnd % bind_set(selem(3), p)
  call bnd % bind_relation(relem(1), flow)
  call bnd % bind_relation(relem(2), d)

  do k = 1, 3
     scell(k) % branch(1) = known_branch(selem(k))
     if (k .lt. 3) scell(k) % branch(2) = &
          & known_branch(scell(k + 1))
  end do
  do k = 1, 2
     rcell(k) % branch(1) = known_branch(relem(k))
     if (k .lt. 2) rcell(k) % branch(2) = &
          & known_branch(rcell(k + 1))
  end do

  g % branch(1) = known_branch(scell(1))
  g % branch(2) = known_branch(rcell(1))

  call check_ownership(nfail)
  call check_signature_closure(nfail)
  call check_ternary_stands(nfail)
  call check_coexistence(nfail)
  call check_validity(nfail)

  call verdict(nfail, "level 3")

contains

  !===================================================================!
  ! The graph owns exactly the three carriers and the two relations
  ! it was handed - by structural identity, through the generators
  ! alone: has_set for the carriers, and a scan of relation_at
  ! identities for the relations, composed here rather than added
  ! to the contract.
  !===================================================================!

  subroutine check_ownership(nfail)

    integer, intent(inout) :: nfail

    call report(num_member_sets(g) .eq. 3, &
         & "the graph owns three member sets", nfail)
    call report(num_relations(g) .eq. 2, &
         & "and two relations", nfail)

    call report(has_set(g, bnd, x) .and. has_set(g, bnd, o) .and. &
         &      has_set(g, bnd, p), &
         & "X, O and P are its own, by identity", nfail)

    call report(graph_holds_relation(g, bnd, flow), &
         & "the flow is owned, the same declared relation", nfail)
    call report(graph_holds_relation(g, bnd, d), &
         & "and so is the derived dependency", nfail)

  end subroutine check_ownership

  !===================================================================!
  ! Signature closure: every slot of every owned relation resolves
  ! to a carrier the graph holds.
  !===================================================================!

  subroutine check_signature_closure(nfail)

    integer, intent(inout) :: nfail

    class(relation), pointer       :: rp
    type(graph) :: dom
    integer                        :: k, s
    logical                        :: ok

    ok = .true.
    do k = 1, num_relations(g)
       rp => relation_at(g, bnd, k)
       do s = 1, rp % arity()
          dom = rp % domain(s)
          ok = ok .and. has_set(g, bnd, dom)
       end do
    end do
    call report(ok, &
         & "every relation slot resolves to an owned carrier", nfail)

  end subroutine check_signature_closure

  !===================================================================!
  ! The ternary relation remains first-class inside the graph: the
  ! owned flow answers arity three and its six tuples; the owned
  ! dependency answers its one. Nothing asked it to become binary.
  !===================================================================!

  subroutine check_ternary_stands(nfail)

    integer, intent(inout) :: nfail

    class(relation), pointer :: rp
    integer                  :: k

    do k = 1, num_relations(g)
       rp => relation_at(g, bnd, k)
       if (rp % same_as(flow)) then
          call report(rp % arity() .eq. 3 .and. &
               &      rp % num_tuples() .eq. 6 .and. &
               &      rp % has([OP_TIMES, SLOT_E, PORT_OUT]), &
               & "the owned flow is still ternary, six tuples strong", &
               & nfail)
       else if (rp % same_as(d)) then
          call report(rp % num_tuples() .eq. 1 .and. &
               &      rp % has([OP_PLUS, OP_TIMES]), &
               & "the owned dependency still holds its one derived pair", &
               & nfail)
       end if
    end do

  end subroutine check_ternary_stands

  !===================================================================!
  ! Two relations of one signature are two citizens: D, and D
  ! restricted to the times-only subset - the empty relation over
  ! O x O - coexist in one graph with their own identities and
  ! their own extensions.
  !===================================================================!

  subroutine check_coexistence(nfail)

    integer, intent(inout) :: nfail

    type(graph)               :: times_only
    type(stored_relation)          :: d_empty
    type(graph)             , target :: aux
    type(graph)             , target :: scell2(1), selem2(1)
    type(graph)             , target :: rcell2(2), relem2(2)
    type(relational_binding)         :: bnd2
    integer                          :: k2
    class(relation), pointer       :: r1, r2
    type(graph) :: da, db
    logical                        :: ok

    call times_only % declare()
    call sets       % bind(times_only, listed_set_representation([OP_TIMES]))
    call inclusions % include_in(times_only, o)
    d_empty    = restrict_slot(d, 1, times_only, sets, inclusions)

    call report(d_empty % num_tuples() .eq. 0, &
         & "restricting D to times-only leaves the empty relation", nfail)

    ! 'coexistence': (S, P) as one sequence on each branch.
    call aux % declare()
    do k2 = 1, 1
       call scell2(k2) % declare()
       call selem2(k2) % declare()
    end do
    do k2 = 1, 2
       call rcell2(k2) % declare()
       call relem2(k2) % declare()
    end do

    call bnd2 % bind_set(selem2(1), o)
    call bnd2 % bind_relation(relem2(1), d)
    call bnd2 % bind_relation(relem2(2), d_empty)

    do k2 = 1, 1
       scell2(k2) % branch(1) = known_branch(selem2(k2))
       if (k2 .lt. 1) scell2(k2) % branch(2) = &
            & known_branch(scell2(k2 + 1))
    end do
    do k2 = 1, 2
       rcell2(k2) % branch(1) = known_branch(relem2(k2))
       if (k2 .lt. 2) rcell2(k2) % branch(2) = &
            & known_branch(rcell2(k2 + 1))
    end do

    aux % branch(1) = known_branch(scell2(1))
    aux % branch(2) = known_branch(rcell2(1))

    r1 => relation_at(aux, bnd2, 1)
    r2 => relation_at(aux, bnd2, 2)

    call report(.not. r1 % same_as(r2), &
         & "one signature, two citizens - no collision", nfail)

    da = r1 % domain(1)
    db = r2 % domain(1)
    ok = da % same_as(db)
    da = r1 % domain(2)
    db = r2 % domain(2)
    ok = ok .and. da % same_as(db)
    call report(ok, &
         & "over the very same signature, slot for slot", nfail)

    call report(r1 % has([OP_PLUS, OP_TIMES]) .and. &
         &      r2 % num_tuples() .eq. 0, &
         & "and each keeps its own extension", nfail)

  end subroutine check_coexistence

  !===================================================================!
  ! Is this declared relation among the graph's own - composed from
  ! the generators, as the generation rule asks.
  !===================================================================!

  !===================================================================!
  ! Validity, ANSWERED. The specimen is valid; three malformations of
  ! it are representable and invalid. These were refusals while a
  ! constructor enforced them; a view over a graph it did not build
  ! reports instead.
  !
  !     foreign   a relation over a carrier the graph does not hold
  !     dupset    one domain seated twice
  !     duprel    one relation seated twice
  !===================================================================!

  subroutine check_validity(nfail)

    integer, intent(inout) :: nfail

    type(graph), target      :: h, scell(2), selem(2), rcell(2), relem(2)
    type(relational_binding) :: b
    type(stored_relation)    :: uses
    integer                  :: j

    call report(relational_valid(g, bnd), &
         & "the specimen is a valid relational graph", nfail)

    uses = stored_relation('uses', [o, x], &
         & reshape([OP_PLUS, SLOT_A], [2, 1]), sets)

    ! O alone, and a relation reaching X: a foreign domain.
    call h % declare()
    do j = 1, 2
       call scell(j) % declare(); call selem(j) % declare()
       call rcell(j) % declare(); call relem(j) % declare()
    end do

    call b % bind_set(selem(1), o)
    call b % bind_set(selem(2), o)          ! and the same domain twice
    call b % bind_relation(relem(1), uses)
    call b % bind_relation(relem(2), uses)  ! and the same relation twice

    scell(1) % branch(1) = known_branch(selem(1))
    rcell(1) % branch(1) = known_branch(relem(1))
    h % branch(1) = known_branch(scell(1))
    h % branch(2) = known_branch(rcell(1))

    call report(.not. relational_valid(h, b), &
         & "a relation over a carrier the graph lacks: INVALID", nfail)

    scell(1) % branch(2) = known_branch(scell(2))
    scell(2) % branch(1) = known_branch(selem(2))
    call report(.not. relational_valid(h, b), &
         & "one domain seated twice: INVALID, S is a set", nfail)

    rcell(1) % branch(2) = known_branch(rcell(2))
    rcell(2) % branch(1) = known_branch(relem(2))
    call report(.not. relational_valid(h, b), &
         & "one relation seated twice: INVALID, P is a set", nfail)

  end subroutine check_validity

  !===================================================================!

  logical function graph_holds_relation(g, b, r)

    type(graph)             , intent(in) :: g
    type(relational_binding), intent(in) :: b
    class(relation)         , intent(in) :: r

    class(relation), pointer :: rp
    integer                  :: k

    graph_holds_relation = .false.
    do k = 1, num_relations(g)
       rp => relation_at(g, b, k)
       if (rp % same_as(r)) graph_holds_relation = .true.
    end do

  end function graph_holds_relation

end program calculator_level_3
