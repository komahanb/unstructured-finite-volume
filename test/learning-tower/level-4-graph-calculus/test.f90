!=====================================================================!
! LEARNING TOWER . LEVEL 4 . GRAPH CALCULUS
!
! The level answers one question: WHAT GRAPH-THEORETIC MEANING MAY
! BE READ from the derived dependency. The relation predict->error
! has existed since Level 2; only THIS rung chooses to read it as
! directed execution - interpretation is explicit, never automatic:
! a binary relation carries the directed reading itself.
!
!      graph-owned D --interpreted--> directed adjacency over O
!      sources = {predict}   sinks = {error}   walk = [predict, error]
!
! owned stable relation, never the selector. And still: execution
! order has meaning; operation laws do not - predict does not yet
! multiply, nothing is evaluated, nothing is trained, and the word
! backprop appears nowhere. No neuron, no layer, no edge carrier:
! the operations themselves are the domain walked.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program learning_level_4

  use learning_assert, only : report, verdict
  use learning_assert, only : SLOT_W, SLOT_X, SLOT_YHAT, SLOT_Y, SLOT_E
  use learning_assert, only : OP_PREDICT, OP_ERROR
  use learning_assert, only : PORT_IN1, PORT_IN2, PORT_OUT
  use graph_fractal        , only : set_graph => graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set        , only : set_map
  use map_label      , only : label_map
  use map_inclusion  , only : inclusion_map, declared_subobject
  use relation_finitary , only : stored_relation, relation
  use relation_binary, only : binary_relation
  use relation_algebra, only : restrict_slot, project_slots, &
       &                             compose_binary
  use graph_algorithms, only : sources, sinks, reachable, &
       &                       topological_order
  use graph_fractal        , only : graph, known_branch, null_branch
  use graph_relational_view, only : relational_binding, &
       & num_member_sets, member_set_at, num_relations, relation_at, &
       & holds_set

  implicit none


  type(set_graph)              :: v, o, p
  type(set_graph)               :: p_out, p_in
  type(stored_relation)          :: flow, t_out3, t_in3
  type(stored_relation)          :: produces, consumes
  class(relation), allocatable   :: d
  type(graph)             , target :: g
  type(graph)             , target :: scell(3), selem(3)
  type(graph)             , target :: rcell(2), relem(2)
  type(relational_binding)         :: bnd
  integer                          :: kcell
  integer                        :: table(3, 6)
  integer                        :: nfail
  type(set_map)     :: sets
  type(inclusion_map)     :: inclusions

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "learning tower . level 4 . graph calculus"
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
  do kcell = 1, 3
     call scell(kcell) % declare()
     call selem(kcell) % declare()
  end do
  do kcell = 1, 2
     call rcell(kcell) % declare()
     call relem(kcell) % declare()
  end do

  call bnd % bind_set(selem(1), v)
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

  ! The interpretive jump, made explicitly - and the selector dies
  ! the moment the reading exists.

  call check_view_domain(nfail)
  call check_sources_and_sinks(nfail)
  call check_reachability(nfail)
  call check_execution_order(nfail)

  call verdict(nfail, "level 4")

contains

  !===================================================================!
  ! The interpretation runs over the operations themselves - no new
  ! vertex carrier, no manufactured edge members.
  !===================================================================!

  subroutine check_view_domain(nfail)

    integer, intent(inout) :: nfail

    type(set_graph) :: dom

    select type (d)
    class is (binary_relation)
       dom = d % source()
    end select
    call report(dom % same_as(o) .and. sets % size_of(dom) .eq. 2, &
         & "the view walks the operations, and nothing invented", nfail)

  end subroutine check_view_domain

  !===================================================================!
  ! Sources and sinks - subobjects of O, proved after the selector's
  ! death, so the borrow demonstrably lives in the graph.
  !===================================================================!

  subroutine check_sources_and_sinks(nfail)

    integer, intent(inout) :: nfail

    type(set_graph) :: src, snk
    type(label_map)     :: labels

    call sources(d, sets, labels, inclusions, src)
    call sinks(d, sets, labels, inclusions, snk)

    call report(sets % size_of(src) .eq. 1 .and. sets % has_in(src, OP_PREDICT) .and. &
         &      .not. sets % has_in(src, OP_ERROR), &
         & "sources = { predict }, the selector long dead", nfail)
    call report(sets % size_of(snk) .eq. 1 .and. sets % has_in(snk, OP_ERROR) .and. &
         &      .not. sets % has_in(snk, OP_PREDICT), &
         & "sinks = { error }", nfail)

    call report(declared_subobject(src, o, inclusions) .and. declared_subobject(snk, o, inclusions), &
         & "both answers stand embedded in the operations", nfail)

  end subroutine check_sources_and_sinks

  !===================================================================!
  ! Directional dependency: predict reaches error, never the
  ! reverse - and each reaches itself by the zero-length path.
  !===================================================================!

  subroutine check_reachability(nfail)

    integer, intent(inout) :: nfail

    call report(reachable(d, sets, OP_PREDICT, OP_ERROR), &
         & "reachable(predict, sets, error) = true", nfail)
    call report(.not. reachable(d, sets, OP_ERROR, OP_PREDICT), &
         & "reachable(error, sets, predict) = false", nfail)
    call report(reachable(d, sets, OP_PREDICT, OP_PREDICT) .and. &
         &      reachable(d, sets, OP_ERROR, OP_ERROR), &
         & "each operation reaches itself by the zero-length path", nfail)
    call report(.not. reachable(d, sets, 7, OP_ERROR), &
         & "an outsider reaches nothing", nfail)

  end subroutine check_reachability

  !===================================================================!
  ! One walk, one order: [predict, error], exactly - execution
  ! order has meaning; operation laws still do not.
  !===================================================================!

  subroutine check_execution_order(nfail)

    integer, intent(inout) :: nfail

    integer, allocatable :: order(:)

    call topological_order(d, sets, order)

    call report(size(order) .eq. 2 .and. &
         &      order(1) .eq. OP_PREDICT .and. order(2) .eq. OP_ERROR, &
         & "the execution order is [predict, error], exactly", nfail)

  end subroutine check_execution_order

end program learning_level_4
