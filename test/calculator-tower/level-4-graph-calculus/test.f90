!=====================================================================!
! CALCULATOR TOWER . LEVEL 4 . GRAPH CALCULUS
!
! The level answers one question: WHAT GRAPH-THEORETIC QUESTIONS CAN
! BE ASKED. The dependency road remains one road,
!
!      T_flow --algebra--> D --admitted--> GAMMA --interpreted-->
!      directed adjacency --walked--> sources, sinks, reach, order
!
! and no second dependency is ever built: the view borrows the
! GRAPH-OWNED relation - the external selector is deallocated the
! moment the view exists, and every answer still stands, which pins
! the borrow where the ownership law says it lives.
!
! The architect-owned truths: sources = {+}, sinks = {x},
! + reaches x and never the reverse, and the one topological walk
! is [+, x]. Still no values, still no arithmetic: D says only
! that + structurally precedes x.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program calculator_level_4

  use calculator_assert, only : report, verdict
  use calculator_assert, only : SLOT_A, SLOT_B, SLOT_C, SLOT_D, SLOT_E
  use calculator_assert, only : OP_PLUS, OP_TIMES
  use calculator_assert, only : PORT_IN1, PORT_IN2, PORT_OUT
  use graph_carrier    , only : counted_set, subset_set, member_set
  use graph_relation   , only : stored_relation, relation
  use graph_relation_algebra, only : restrict_slot, project_slots, &
       &                             compose_binary
  use graph_structure  , only : relational_graph, held_set, held_relation
  use graph_profile    , only : directed_adjacency_view
  use graph_algorithms , only : sources, sinks, reachable, &
       &                        topological_order
  use fractal_graph        , only : graph
  use graph_relational_view, only : relational_binding
  use relational_fixture   , only : fractal_fixture

  implicit none

  type(fractal_fixture)             :: fx_
  type(graph)             , pointer :: fg_
  type(relational_binding), pointer :: fb_

  type(counted_set)              :: x, o, p
  type(subset_set)               :: p_out, p_in
  type(stored_relation)          :: flow, t_out3, t_in3
  type(stored_relation)          :: produces, consumes
  class(relation), allocatable   :: d
  type(relational_graph), target :: g
  type(directed_adjacency_view)  :: view
  integer                        :: table(3, 6)
  integer                        :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "calculator tower . level 4 . graph calculus"
  write(*,'(1x,a)') "============================================="

  x = counted_set('value-slots', 5)
  o = counted_set('operations' , 2)
  p = counted_set('ports'      , 3)

  table(:, 1) = [OP_PLUS , SLOT_A, PORT_IN1]
  table(:, 2) = [OP_PLUS , SLOT_B, PORT_IN2]
  table(:, 3) = [OP_PLUS , SLOT_C, PORT_OUT]
  table(:, 4) = [OP_TIMES, SLOT_C, PORT_IN1]
  table(:, 5) = [OP_TIMES, SLOT_D, PORT_IN2]
  table(:, 6) = [OP_TIMES, SLOT_E, PORT_OUT]

  flow = stored_relation('flow', [o, x, p], table)

  p_out    = subset_set('output-port', p, [PORT_OUT])
  p_in     = subset_set('input-ports', p, [PORT_IN1, PORT_IN2])
  t_out3   = restrict_slot(flow, 3, p_out)
  t_in3    = restrict_slot(flow, 3, p_in)
  produces = project_slots(t_out3, [1, 2])
  consumes = project_slots(t_in3 , [2, 1])
  d        = compose_binary(produces, consumes)

  g = relational_graph('calculator', &
       & [held_set(x), held_set(o), held_set(p)], &
       & [held_relation(flow), held_relation(d)])

  ! The interpretation reads the GRAPH-OWNED dependency; the
  ! selector has served and may die.
  call fx_ % to_fractal(g, fg_, fb_)
  view = directed_adjacency_view(fg_, fb_, d)
  deallocate(d)

  call check_sources_and_sinks(nfail)
  call check_reachability(nfail)
  call check_walk(nfail)

  call verdict(nfail, "level 4")

contains

  !===================================================================!
  ! sources(D) = {+} and sinks(D) = {x} - as subobjects of the
  ! operations, proved after the selector's death, so the borrow
  ! demonstrably lives in the graph.
  !===================================================================!

  subroutine check_sources_and_sinks(nfail)

    integer, intent(inout) :: nfail

    type(subset_set) :: src, snk

    src = sources(view)
    snk = sinks(view)

    call report(src % size() .eq. 1 .and. src % has(OP_PLUS), &
         & "sources(D) = { + }, the selector long dead", nfail)
    call report(snk % size() .eq. 1 .and. snk % has(OP_TIMES), &
         & "sinks(D) = { x }", nfail)

    call report(src % is_subobject_of(o), &
         & "the sources stand embedded in the operations", nfail)
    call report(snk % is_subobject_of(o), &
         & "and so do the sinks", nfail)

  end subroutine check_sources_and_sinks

  !===================================================================!
  ! + reaches x, never the reverse - and each reaches itself by
  ! the zero-length path.
  !===================================================================!

  subroutine check_reachability(nfail)

    integer, intent(inout) :: nfail

    call report(reachable(view, OP_PLUS, OP_TIMES), &
         & "reachable(+, x) = true", nfail)
    call report(.not. reachable(view, OP_TIMES, OP_PLUS), &
         & "reachable(x, +) = false", nfail)
    call report(reachable(view, OP_PLUS, OP_PLUS), &
         & "and + reaches itself by the zero-length path", nfail)

  end subroutine check_reachability

  !===================================================================!
  ! One walk, one order: [+, x], and nothing else.
  !===================================================================!

  subroutine check_walk(nfail)

    integer, intent(inout) :: nfail

    integer, allocatable :: order(:)

    call topological_order(view, order)

    call report(size(order) .eq. 2 .and. &
         &      order(1) .eq. OP_PLUS .and. order(2) .eq. OP_TIMES, &
         & "the topological walk is [+, x], exactly", nfail)

  end subroutine check_walk

end program calculator_level_4
