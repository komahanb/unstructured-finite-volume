!=====================================================================!
! LEARNING TOWER . LEVEL 3 . THE RELATIONAL GRAPH
!
! The level answers one question: CAN THE LEARNING CARRIERS AND
! THEIR DERIVED RELATIONS COEXIST AS ONE STRUCTURE. The persistent
! object becomes
!
!      GAMMA = ( { V, O, P }, { T_flow, D } )
!
! where D is DERIVED once more by the approved Level-2 road and
! then ADMITTED - the container stores structure, it infers no
! learning semantics and re-derives nothing. The most important
! truth of the rung: the TERNARY flow survives graph ownership
! unchanged - nothing asks the native learning schema to become a
! binary graph merely because it now lives inside something called
! a graph. D remains a plain relation in a graph: no directed edge,
! no source, no walk - interpretation is Level 4's business. And
! still no number anywhere: no data, no parameter, no law.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program learning_level_3

  use learning_assert, only : report, verdict
  use learning_assert, only : SLOT_W, SLOT_X, SLOT_YHAT, SLOT_Y, SLOT_E
  use learning_assert, only : OP_PREDICT, OP_ERROR
  use learning_assert, only : PORT_IN1, PORT_IN2, PORT_OUT
  use graph_fractal        , only : graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set        , only : set_map
  use map_inclusion  , only : inclusion_map, declared_subobject
  use relation_finitary , only : stored_relation, relation
  use relation_algebra, only : restrict_slot, project_slots, &
       &                             compose_binary
  use graph_fractal        , only : graph, known_branch, null_branch
  use view_relational, only : relational_binding, &
       & num_member_sets, member_set_at, num_relations, relation_at, &
       & has_set, relational_valid

  implicit none

  type(graph)              :: v, o, p
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
  write(*,'(1x,a)') "learning tower . level 3 . relational graph"
  write(*,'(1x,a)') "============================================="

  call v % declare()
  call sets % bind(v, counted_set_representation(5))
  call o % declare()
  call sets % bind(o, counted_set_representation(2))
  call p % declare()
  call sets % bind(p, counted_set_representation(3))

  table(:, 1) = [OP_PREDICT, SLOT_W   , PORT_IN1]
  table(:, 2) = [OP_PREDICT, SLOT_X   , PORT_IN2]
  table(:, 3) = [OP_PREDICT, SLOT_YHAT, PORT_OUT]
  table(:, 4) = [OP_ERROR  , SLOT_YHAT, PORT_IN1]
  table(:, 5) = [OP_ERROR  , SLOT_Y   , PORT_IN2]
  table(:, 6) = [OP_ERROR  , SLOT_E   , PORT_OUT]
  flow = stored_relation('flow', [o, v, p], table, sets)

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

  ! 'learning': (S, P) as one sequence on each branch.
  call g % declare()
  do k = 1, 3
     call scell(k) % declare()
     call selem(k) % declare()
  end do
  do k = 1, 2
     call rcell(k) % declare()
     call relem(k) % declare()
  end do

  call bnd % bind_set(selem(1), v)
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
  call check_ternary_survives(nfail)
  call check_dependency_survives(nfail)
  call check_signature_closure(nfail)
  call check_graph_identity(nfail)
  call check_validity(nfail)

  call verdict(nfail, "level 3")

contains

  !===================================================================!
  ! Three carriers and two relations, owned by identity - the
  ! relation question composed locally from the generators, as the
  ! generation rule asks.
  !===================================================================!

  subroutine check_ownership(nfail)

    integer, intent(inout) :: nfail

    call report(num_member_sets(g) .eq. 3 .and. &
         &      num_relations(g) .eq. 2, &
         & "the graph owns three member sets and two relations", nfail)

    call report(has_set(g, bnd, v) .and. has_set(g, bnd, o) .and. &
         &      has_set(g, bnd, p), &
         & "V, O and P are its own, by identity", nfail)

    call report(graph_holds_relation(g, bnd, flow), &
         & "the flow is owned, the same declared relation", nfail)
    call report(graph_holds_relation(g, bnd, d), &
         & "and so is the derived dependency", nfail)

  end subroutine check_ownership

  !===================================================================!
  ! The ternary flow survives ownership: arity three, six tuples,
  ! representative members of BOTH operations - the native learning
  ! schema, unchanged inside the container.
  !===================================================================!

  subroutine check_ternary_survives(nfail)

    integer, intent(inout) :: nfail

    class(relation), pointer :: rp
    integer                  :: k

    do k = 1, num_relations(g)
       rp => relation_at(g, bnd, k)
       if (rp % same_as(flow)) then
          call report(rp % arity() .eq. 3 .and. &
               &      rp % num_tuples() .eq. 6 .and. &
               &      rp % has([OP_PREDICT, SLOT_W   , PORT_IN1]) .and. &
               &      rp % has([OP_PREDICT, SLOT_YHAT, PORT_OUT]) .and. &
               &      rp % has([OP_ERROR  , SLOT_YHAT, PORT_IN1]) .and. &
               &      rp % has([OP_ERROR  , SLOT_E   , PORT_OUT]), &
               & "the owned flow is still ternary, six tuples strong", &
               & nfail)
       end if
    end do

  end subroutine check_ternary_survives

  !===================================================================!
  ! The derived dependency survives too: binary, one tuple, both
  ! slots the operations - admitted, never recreated.
  !===================================================================!

  subroutine check_dependency_survives(nfail)

    integer, intent(inout) :: nfail

    class(relation), pointer       :: rp
    type(graph) :: dom
    integer                        :: k

    do k = 1, num_relations(g)
       rp => relation_at(g, bnd, k)
       if (rp % same_as(d)) then
          call report(rp % arity() .eq. 2 .and. &
               &      rp % num_tuples() .eq. 1 .and. &
               &      rp % has([OP_PREDICT, OP_ERROR]), &
               & "the owned dependency still holds its one derived pair", &
               & nfail)
          dom = rp % domain(1)
          call report(dom % same_as(o), &
               & "its first slot is the operations", nfail)
          dom = rp % domain(2)
          call report(dom % same_as(o), &
               & "and so is its second", nfail)
       end if
    end do

  end subroutine check_dependency_survives

  !===================================================================!
  ! Signature closure: every slot of every owned relation resolves
  ! to an owned carrier - no equal-sized foreign domain sneaks in.
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
  ! GAMMA is a declared graph: itself, and not an identically built
  ! twin - extension equality never collapses identity.
  !===================================================================!

  subroutine check_graph_identity(nfail)

    integer, intent(inout) :: nfail

    type(graph)             , target :: g2
    type(graph)             , target :: scell2(3), selem2(3)
    type(relational_binding)         :: bnd2
    integer                          :: k2

    call report(g % same_as(g), &
         & "the learning graph is itself", nfail)

    ! 'learning again': (S, P) as one sequence on each branch.
    call g2 % declare()
    do k2 = 1, 3
       call scell2(k2) % declare()
       call selem2(k2) % declare()
    end do

    call bnd2 % bind_set(selem2(1), v)
    call bnd2 % bind_set(selem2(2), o)
    call bnd2 % bind_set(selem2(3), p)

    do k2 = 1, 3
       scell2(k2) % branch(1) = known_branch(selem2(k2))
       if (k2 .lt. 3) scell2(k2) % branch(2) = &
            & known_branch(scell2(k2 + 1))
    end do

    g2 % branch(1) = known_branch(scell2(1))
    g2 % branch(2) = null_branch()
    call report(.not. g % same_as(g2), &
         & "and no identically stocked twin is it", nfail)

  end subroutine check_graph_identity

  !===================================================================!
  ! Is this declared relation among the graph's own - composed from
  ! num_relations and relation_at, no convenience API.
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

    uses = stored_relation('uses', [o, v], &
         & reshape([OP_PREDICT, SLOT_W], [2, 1]), sets)

    ! O alone, and a relation reaching V: a foreign domain.
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

end program learning_level_3
