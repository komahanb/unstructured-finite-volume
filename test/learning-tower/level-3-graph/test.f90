!=====================================================================!
! LEARNING TOWER . LEVEL 3 . THE RELATED GRAPH
!
! The level answers one question: CAN THE LEARNING SETS AND
! THEIR DERIVED RELATIONS COEXIST AS ONE STRUCTURE. The persistent
! object becomes
!
!      G = ( { V, O, P }, { R_flow, D } )
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
  use graph_set  , only : index_set, subset, set, unrelated_graph
  use graph_relation , only : stored_relation, relation
  use graph_relation_algebra, only : restrict_slot, project_slots, &
       &                             compose_binary
  use graph_structure, only : related_graph, declared_set, declared_relation

  implicit none

  type(index_set)              :: v, o, p
  type(subset)               :: p_out, p_in
  type(stored_relation)          :: flow, r_out3, r_in3
  type(stored_relation)          :: produces, consumes
  class(relation), allocatable   :: d
  type(related_graph), target :: g
  integer                        :: table(3, 6)
  integer                        :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "learning tower . level 3 . related graph"
  write(*,'(1x,a)') "============================================="

  v = index_set('value-slots', 5)
  o = index_set('operations' , 2)
  p = index_set('ports'      , 3)

  table(:, 1) = [OP_PREDICT, SLOT_W   , PORT_IN1]
  table(:, 2) = [OP_PREDICT, SLOT_X   , PORT_IN2]
  table(:, 3) = [OP_PREDICT, SLOT_YHAT, PORT_OUT]
  table(:, 4) = [OP_ERROR  , SLOT_YHAT, PORT_IN1]
  table(:, 5) = [OP_ERROR  , SLOT_Y   , PORT_IN2]
  table(:, 6) = [OP_ERROR  , SLOT_E   , PORT_OUT]
  flow = stored_relation('flow', [o, v, p], table)

  ! The approved Level-2 road, walked once more; the graph admits
  ! what the algebra derived.
  p_out    = subset('output-port', p, [PORT_OUT])
  p_in     = subset('input-ports', p, [PORT_IN1, PORT_IN2])
  r_out3   = restrict_slot(flow, 3, p_out)
  r_in3    = restrict_slot(flow, 3, p_in)
  produces = project_slots(r_out3, [1, 2])
  consumes = project_slots(r_in3 , [2, 1])
  d        = compose_binary(produces, consumes)

  g = related_graph('learning', &
       & [declared_set(v), declared_set(o), declared_set(p)], &
       & [declared_relation(flow), declared_relation(d)])

  call check_ownership(nfail)
  call check_ternary_survives(nfail)
  call check_dependency_survives(nfail)
  call check_signature_closure(nfail)
  call check_graph_identity(nfail)

  call verdict(nfail, "level 3")

contains

  !===================================================================!
  ! Three sets and two relations, owned by identity - the
  ! relation question composed locally from the generators, as the
  ! generation rule asks.
  !===================================================================!

  subroutine check_ownership(nfail)

    integer, intent(inout) :: nfail

    call report(g % num_sets() .eq. 3 .and. &
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
       if (rp % equals(flow)) then
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
    class(set), allocatable :: dom
    integer                        :: k

    do k = 1, g % num_relations()
       rp => g % relation_at(k)
       if (rp % equals(d)) then
          call report(rp % arity() .eq. 2 .and. &
               &      rp % num_tuples() .eq. 1 .and. &
               &      rp % has([OP_PREDICT, OP_ERROR]), &
               & "the owned dependency still holds its one derived pair", &
               & nfail)
          dom = rp % domain(1)
          call report(dom % equals(o), &
               & "its first slot is the operations", nfail)
          dom = rp % domain(2)
          call report(dom % equals(o), &
               & "and so is its second", nfail)
       end if
    end do

  end subroutine check_dependency_survives

  !===================================================================!
  ! Signature closure: every slot of every owned relation resolves
  ! to an owned set - no equal-sized unheld domain sneaks in.
  !===================================================================!

  subroutine check_signature_closure(nfail)

    integer, intent(inout) :: nfail

    class(relation), pointer       :: rp
    class(set), allocatable :: dom
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
         & "every relation slot resolves to an owned set", nfail)

  end subroutine check_signature_closure

  !===================================================================!
  ! G is a declared graph: itself, and not an identically built
  ! twin - extension equality never collapses identity.
  !===================================================================!

  subroutine check_graph_identity(nfail)

    integer, intent(inout) :: nfail

    type(unrelated_graph) :: g2

    call report(g % equals(g), &
         & "the learning graph is itself", nfail)

    g2 = unrelated_graph('learning again', &
         & [declared_set(v), declared_set(o), declared_set(p)])
    call report(.not. g % equals(g2), &
         & "and no identically stocked twin is it", nfail)

  end subroutine check_graph_identity

  !===================================================================!
  ! Is this declared relation among the graph's own - composed from
  ! num_relations and relation_at, no convenience API.
  !===================================================================!

  logical function graph_holds_relation(g, r)

    type(related_graph), target, intent(in) :: g
    class(relation)               , intent(in) :: r

    class(relation), pointer :: rp
    integer                  :: k

    graph_holds_relation = .false.
    do k = 1, g % num_relations()
       rp => g % relation_at(k)
       if (rp % equals(r)) graph_holds_relation = .true.
    end do

  end function graph_holds_relation

end program learning_level_3
