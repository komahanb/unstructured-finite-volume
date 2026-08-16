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
  use fractal_graph        , only : set_graph => graph
  use graph_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use graph_set_map        , only : set_map
  use graph_label_map      , only : label_map
  use graph_inclusion_map  , only : inclusion_map, declared_subobject
  use graph_relation   , only : stored_relation, relation
  use graph_relation_algebra, only : restrict_slot, project_slots, &
       &                             compose_binary
  use graph_profile    , only : directed_adjacency_view
  use graph_algorithms , only : sources, sinks, reachable, &
       &                        topological_order
  use fractal_graph        , only : graph, known_branch, null_branch
  use graph_relational_view, only : relational_binding, &
       & num_member_sets, member_set_at, num_relations, relation_at, &
       & holds_set

  implicit none


  type(set_graph)              :: x, o, p
  type(set_graph)               :: p_out, p_in
  type(stored_relation)          :: flow, t_out3, t_in3
  type(stored_relation)          :: produces, consumes
  class(relation), allocatable   :: d
  type(graph)             , target :: g
  type(graph)             , target :: scell(3), selem(3)
  type(graph)             , target :: rcell(2), relem(2)
  type(relational_binding)         :: bnd
  integer                          :: kcell
  type(directed_adjacency_view)  :: view
  integer                        :: table(3, 6)
  integer                        :: nfail
  type(set_map)     :: sets
  type(inclusion_map)     :: inclusions

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "calculator tower . level 4 . graph calculus"
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
  do kcell = 1, 3
     call scell(kcell) % declare()
     call selem(kcell) % declare()
  end do
  do kcell = 1, 2
     call rcell(kcell) % declare()
     call relem(kcell) % declare()
  end do

  call bnd % bind_set(selem(1), x)
  call bnd % bind_set(selem(2), o)
  call bnd % bind_set(selem(3), p)
  call bnd % bind_relation(relem(1), flow)
  call bnd % bind_relation(relem(2), d)

  do kcell = 1, 3
     scell(kcell) % branch(1) = known_branch(selem(kcell))
     if (kcell .lt. 3) scell(kcell) % branch(2) = &
          & known_branch(scell(kcell + 1))
  end do
  do kcell = 1, 2
     rcell(kcell) % branch(1) = known_branch(relem(kcell))
     if (kcell .lt. 2) rcell(kcell) % branch(2) = &
          & known_branch(rcell(kcell + 1))
  end do

  g % branch(1) = known_branch(scell(1))
  g % branch(2) = known_branch(rcell(1))

  ! The interpretation reads the GRAPH-OWNED dependency; the
  ! selector has served and may die.
  view = directed_adjacency_view(g, bnd, sets, d)
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

    type(set_graph) :: src, snk
    type(label_map)     :: labels

    call sources(view, sets, labels, inclusions, src)
    call sinks(view, sets, labels, inclusions, snk)

    call report(sets % size_of(src) .eq. 1 .and. sets % has_in(src, OP_PLUS), &
         & "sources(D) = { + }, the selector long dead", nfail)
    call report(sets % size_of(snk) .eq. 1 .and. sets % has_in(snk, OP_TIMES), &
         & "sinks(D) = { x }", nfail)

    call report(declared_subobject(src, o, inclusions), &
         & "the sources stand embedded in the operations", nfail)
    call report(declared_subobject(snk, o, inclusions), &
         & "and so do the sinks", nfail)

  end subroutine check_sources_and_sinks

  !===================================================================!
  ! + reaches x, never the reverse - and each reaches itself by
  ! the zero-length path.
  !===================================================================!

  subroutine check_reachability(nfail)

    integer, intent(inout) :: nfail

    call report(reachable(view, sets, OP_PLUS, OP_TIMES), &
         & "reachable(+, x) = true", nfail)
    call report(.not. reachable(view, sets, OP_TIMES, OP_PLUS), &
         & "reachable(x, sets, +) = false", nfail)
    call report(reachable(view, sets, OP_PLUS, OP_PLUS), &
         & "and + reaches itself by the zero-length path", nfail)

  end subroutine check_reachability

  !===================================================================!
  ! One walk, one order: [+, x], and nothing else.
  !===================================================================!

  subroutine check_walk(nfail)

    integer, intent(inout) :: nfail

    integer, allocatable :: order(:)

    call topological_order(view, sets, order)

    call report(size(order) .eq. 2 .and. &
         &      order(1) .eq. OP_PLUS .and. order(2) .eq. OP_TIMES, &
         & "the topological walk is [+, x], exactly", nfail)

  end subroutine check_walk

end program calculator_level_4
