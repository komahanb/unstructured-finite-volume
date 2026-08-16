!=====================================================================!
! DERIVATIVE ACTION TOWER . LEVEL 3 . THE RELATIONAL GRAPH
!
! The level answers one question: CAN THE SYMBOLIC DERIVATIVE
! SPECIMEN LIVE AS ONE OWNED RELATIONAL STRUCTURE. The persistent
! object becomes
!
!      GAMMA = ( { V, O, P }, { T_flow, D } )
!
! where D is DERIVED once more by the certified Level-2 road and
! then ADMITTED - the container stores structure, it infers nothing
! and re-derives nothing. The ternary flow survives ownership
! unchanged: nothing asks the native computation schema to become
! binary merely because it now lives inside something called a
! graph. What ownership contributes here is structural closure and
! identity - every relation slot resolves to an owned carrier, and
! the owned citizens are the SAME declared relations, not copies of
! their extensions. Nothing derivative-specific asked the container
! for anything: no tangent seat, no adjoint seat, no derivative
! metadata. Interpretation is Level 4's business; values are
! Level 5's; and nothing differentiates anything.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program derivative_level_3

  use derivative_assert, only : report, verdict
  use derivative_assert, only : SLOT_X, SLOT_Y, SLOT_U, SLOT_Z
  use derivative_assert, only : OP_PRODUCT, OP_SUM
  use derivative_assert, only : PORT_IN1, PORT_IN2, PORT_OUT
  use graph_carrier    , only : counted_set, subset_set, member_set
  use graph_relation   , only : stored_relation, relation
  use graph_relation_algebra, only : restrict_slot, project_slots, &
       &                             compose_binary
  use fractal_graph        , only : graph, known_branch, null_branch
  use graph_relational_view, only : relational_binding, &
       & num_member_sets, member_set_at, num_relations, relation_at, &
       & holds_set

  implicit none

  type(counted_set)              :: v, o, p
  type(subset_set)               :: p_out, p_in
  type(stored_relation)          :: flow
  class(relation), allocatable   :: d
  type(graph)             , target :: g
  type(graph)             , target :: scell(3), selem(3)
  type(graph)             , target :: rcell(2), relem(2)
  type(relational_binding)         :: bnd
  integer                          :: k
  integer                        :: table(3, 6)
  integer                        :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "derivative action tower . level 3 . graph"
  write(*,'(1x,a)') "============================================="

  v = counted_set('value-slots', 4)
  o = counted_set('operations' , 2)
  p = counted_set('ports'      , 3)

  table(:, 1) = [OP_PRODUCT, SLOT_X, PORT_IN1]
  table(:, 2) = [OP_PRODUCT, SLOT_Y, PORT_IN2]
  table(:, 3) = [OP_PRODUCT, SLOT_U, PORT_OUT]
  table(:, 4) = [OP_SUM    , SLOT_U, PORT_IN1]
  table(:, 5) = [OP_SUM    , SLOT_Y, PORT_IN2]
  table(:, 6) = [OP_SUM    , SLOT_Z, PORT_OUT]
  flow = stored_relation('flow', [o, v, p], table)

  ! The certified Level-2 road, walked once more; the graph admits
  ! what the algebra derived.
  p_out = subset_set('output-port', p, [PORT_OUT])
  p_in  = subset_set('input-ports', p, [PORT_IN1, PORT_IN2])
  d = compose_binary( &
       & project_slots(restrict_slot(flow, 3, p_out), [1, 2]), &
       & project_slots(restrict_slot(flow, 3, p_in ), [2, 1]))

  ! 'derivative specimen': (S, P) as one sequence on each branch.
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

  call verdict(nfail, "level 3")

contains

  !===================================================================!
  ! Three carriers and two relations, owned by identity - the
  ! relation question composed locally from the generators.
  !===================================================================!

  subroutine check_ownership(nfail)

    integer, intent(inout) :: nfail

    call report(num_member_sets(g) .eq. 3 .and. &
         &      num_relations(g) .eq. 2, &
         & "the graph owns three member sets and two relations", nfail)

    call report(holds_set(g, bnd, v) .and. holds_set(g, bnd, o) .and. &
         &      holds_set(g, bnd, p), &
         & "V, O and P are its own, by identity", nfail)

    call report(graph_holds_relation(g, bnd, flow), &
         & "the flow is owned, the same declared relation", nfail)
    call report(graph_holds_relation(g, bnd, d), &
         & "and so is the derived dependency", nfail)

  end subroutine check_ownership

  !===================================================================!
  ! The ternary flow survives ownership: arity three, six tuples,
  ! representative members of BOTH operations.
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
               &      rp % has([OP_PRODUCT, SLOT_X, PORT_IN1]) .and. &
               &      rp % has([OP_PRODUCT, SLOT_U, PORT_OUT]) .and. &
               &      rp % has([OP_SUM    , SLOT_Y, PORT_IN2]) .and. &
               &      rp % has([OP_SUM    , SLOT_Z, PORT_OUT]), &
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
    class(member_set), allocatable :: dom
    integer                        :: k

    do k = 1, num_relations(g)
       rp => relation_at(g, bnd, k)
       if (rp % same_as(d)) then
          call report(rp % arity() .eq. 2 .and. &
               &      rp % num_tuples() .eq. 1 .and. &
               &      rp % has([OP_PRODUCT, OP_SUM]), &
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
    class(member_set), allocatable :: dom
    integer                        :: k, s
    logical                        :: ok

    ok = .true.
    do k = 1, num_relations(g)
       rp => relation_at(g, bnd, k)
       do s = 1, rp % arity()
          dom = rp % domain(s)
          ok = ok .and. holds_set(g, bnd, dom)
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
         & "the specimen graph is itself", nfail)

    ! 'derivative specimen again': (S, P) as one sequence on each branch.
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

end program derivative_level_3
