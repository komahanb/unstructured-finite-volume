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
  use graph_carrier  , only : counted_set, subset_set, member_set
  use graph_relation , only : stored_relation, relation
  use graph_relation_algebra, only : restrict_slot, project_slots, &
       &                             compose_binary
  use graph_structure, only : relational_graph, held_set, held_relation

  implicit none

  type(counted_set)              :: v, o, p
  type(subset_set)               :: p_out, p_in
  type(stored_relation)          :: flow, t_out3, t_in3
  type(stored_relation)          :: produces, consumes
  class(relation), allocatable   :: d
  type(relational_graph), target :: g
  integer                        :: table(3, 6)
  integer                        :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "learning tower . level 3 . relational graph"
  write(*,'(1x,a)') "============================================="

  v = counted_set('value-slots', 5)
  o = counted_set('operations' , 2)
  p = counted_set('ports'      , 3)

  table(:, 1) = [OP_PREDICT, SLOT_W   , PORT_IN1]
  table(:, 2) = [OP_PREDICT, SLOT_X   , PORT_IN2]
  table(:, 3) = [OP_PREDICT, SLOT_YHAT, PORT_OUT]
  table(:, 4) = [OP_ERROR  , SLOT_YHAT, PORT_IN1]
  table(:, 5) = [OP_ERROR  , SLOT_Y   , PORT_IN2]
  table(:, 6) = [OP_ERROR  , SLOT_E   , PORT_OUT]
  flow = stored_relation('flow', [o, v, p], table)

  ! The approved Level-2 road, walked once more; the graph admits
  ! what the algebra derived.
  p_out    = subset_set('output-port', p, [PORT_OUT])
  p_in     = subset_set('input-ports', p, [PORT_IN1, PORT_IN2])
  t_out3   = restrict_slot(flow, 3, p_out)
  t_in3    = restrict_slot(flow, 3, p_in)
  produces = project_slots(t_out3, [1, 2])
  consumes = project_slots(t_in3 , [2, 1])
  d        = compose_binary(produces, consumes)

  g = relational_graph('learning', &
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
  ! relation question composed locally from the generators, as the
  ! generation rule asks.
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
  ! representative members of BOTH operations - the native learning
  ! schema, unchanged inside the container.
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
    class(member_set), allocatable :: dom
    integer                        :: k

    do k = 1, g % num_relations()
       rp => g % relation_at(k)
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
  ! GAMMA is a declared graph: itself, and not an identically built
  ! twin - extension equality never collapses identity.
  !===================================================================!

  subroutine check_graph_identity(nfail)

    integer, intent(inout) :: nfail

    type(relational_graph) :: g2
    type(held_relation)    :: none(0)

    call report(g % same_as(g), &
         & "the learning graph is itself", nfail)

    g2 = relational_graph('learning again', &
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

end program learning_level_3
