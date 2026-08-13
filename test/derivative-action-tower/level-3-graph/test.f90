!=====================================================================!
! DERIVATIVE ACTION TOWER . LEVEL 3 . THE RELATIONAL GRAPH
!
! The level answers one question: CAN THE SYMBOLIC DERIVATIVE
! SPECIMEN LIVE AS ONE OWNED RELATIONAL STRUCTURE. The persistent
! object becomes
!
!      G = ( { V, O, P }, { R_flow, D } )
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
  use graph_structure  , only : relational_graph, held_set, held_relation

  implicit none

  type(counted_set)              :: v, o, p
  type(subset_set)               :: p_out, p_in
  type(stored_relation)          :: flow
  class(relation), allocatable   :: d
  type(relational_graph), target :: g
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

  g = relational_graph('derivative specimen', &
       & [held_set(v), held_set(o), held_set(p)], &
       & [held_relation(flow), held_relation(d)])

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

    call report(g % num_member_sets() .eq. 3 .and. &
         &      g % num_relations() .eq. 2, &
         & "the graph owns three member sets and two relations", nfail)

    call report(g % holds_set(v) .and. g % holds_set(o) .and. &
         &      g % holds_set(p), &
         & "V, O and P are its own, by identity", nfail)

    call report(graph_holds_relation(g, flow), &
         & "the flow is owned, the same declared relation", nfail)
    call report(graph_holds_relation(g, d), &
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

    do k = 1, g % num_relations()
       rp => g % relation_at(k)
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

    do k = 1, g % num_relations()
       rp => g % relation_at(k)
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
    do k = 1, g % num_relations()
       rp => g % relation_at(k)
       do s = 1, rp % arity()
          dom = rp % domain(s)
          ok = ok .and. g % holds_set(dom)
       end do
    end do
    call report(ok, &
         & "every relation slot resolves to an owned carrier", nfail)

  end subroutine check_signature_closure

  !===================================================================!
  ! G is a declared graph: itself, and not an identically built
  ! twin - extension equality never collapses identity.
  !===================================================================!

  subroutine check_graph_identity(nfail)

    integer, intent(inout) :: nfail

    type(relational_graph) :: g2
    type(held_relation)    :: none(0)

    call report(g % same_as(g), &
         & "the specimen graph is itself", nfail)

    g2 = relational_graph('derivative specimen again', &
         & [held_set(v), held_set(o), held_set(p)], none)
    call report(.not. g % same_as(g2), &
         & "and no identically stocked twin is it", nfail)

  end subroutine check_graph_identity

  !===================================================================!
  ! Is this declared relation among the graph's own - composed from
  ! num_relations and relation_at, no convenience API.
  !===================================================================!

  logical function graph_holds_relation(g, r)

    type(relational_graph), target, intent(in) :: g
    class(relation)               , intent(in) :: r

    class(relation), pointer :: rp
    integer                  :: k

    graph_holds_relation = .false.
    do k = 1, g % num_relations()
       rp => g % relation_at(k)
       if (rp % same_as(r)) graph_holds_relation = .true.
    end do

  end function graph_holds_relation

end program derivative_level_3
