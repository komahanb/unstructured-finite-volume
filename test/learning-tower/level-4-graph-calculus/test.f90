!=====================================================================!
! LEARNING TOWER . LEVEL 4 . GRAPH CALCULUS
!
! The level answers one question: WHAT GRAPH-THEORETIC MEANING MAY
! BE READ from the derived dependency. The relation predict->error
! has existed since Level 2; only THIS rung chooses to read it as
! directed execution - interpretation is explicit, never automatic:
! a binary relation is not a directed graph until a view says so.
!
!      graph-owned D --interpreted--> directed adjacency over O
!      sources = {predict}   sinks = {error}   walk = [predict, error]
!
! The external selector D is DEALLOCATED the moment the view exists
! and every algorithm still answers: the view borrows the graph's
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
  use graph_carrier  , only : counted_set, subset_set, member_set
  use graph_relation , only : stored_relation, relation
  use graph_relation_algebra, only : restrict_slot, project_slots, &
       &                             compose_binary
  use graph_profile  , only : directed_adjacency_view
  use graph_algorithms, only : sources, sinks, reachable, &
       &                       topological_order
  use fractal_graph        , only : graph, known_branch, null_branch
  use graph_relational_view, only : relational_binding, &
       & num_member_sets, member_set_at, num_relations, relation_at, &
       & holds_set

  implicit none


  type(counted_set)              :: v, o, p
  type(subset_set)               :: p_out, p_in
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

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "learning tower . level 4 . graph calculus"
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

  p_out    = subset_set('output-port', p, [PORT_OUT])
  p_in     = subset_set('input-ports', p, [PORT_IN1, PORT_IN2])
  t_out3   = restrict_slot(flow, 3, p_out)
  t_in3    = restrict_slot(flow, 3, p_in)
  produces = project_slots(t_out3, [1, 2])
  consumes = project_slots(t_in3 , [2, 1])
  d        = compose_binary(produces, consumes)

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
  view = directed_adjacency_view(g, bnd, d)
  deallocate(d)

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

    class(member_set), allocatable :: dom

    dom = view % domain()
    call report(dom % same_as(o) .and. dom % size() .eq. 2, &
         & "the view walks the operations, and nothing invented", nfail)

  end subroutine check_view_domain

  !===================================================================!
  ! Sources and sinks - subobjects of O, proved after the selector's
  ! death, so the borrow demonstrably lives in the graph.
  !===================================================================!

  subroutine check_sources_and_sinks(nfail)

    integer, intent(inout) :: nfail

    type(subset_set) :: src, snk

    src = sources(view)
    snk = sinks(view)

    call report(src % size() .eq. 1 .and. src % has(OP_PREDICT) .and. &
         &      .not. src % has(OP_ERROR), &
         & "sources = { predict }, the selector long dead", nfail)
    call report(snk % size() .eq. 1 .and. snk % has(OP_ERROR) .and. &
         &      .not. snk % has(OP_PREDICT), &
         & "sinks = { error }", nfail)

    call report(src % is_subobject_of(o) .and. snk % is_subobject_of(o), &
         & "both answers stand embedded in the operations", nfail)

  end subroutine check_sources_and_sinks

  !===================================================================!
  ! Directional dependency: predict reaches error, never the
  ! reverse - and each reaches itself by the zero-length path.
  !===================================================================!

  subroutine check_reachability(nfail)

    integer, intent(inout) :: nfail

    call report(reachable(view, OP_PREDICT, OP_ERROR), &
         & "reachable(predict, error) = true", nfail)
    call report(.not. reachable(view, OP_ERROR, OP_PREDICT), &
         & "reachable(error, predict) = false", nfail)
    call report(reachable(view, OP_PREDICT, OP_PREDICT) .and. &
         &      reachable(view, OP_ERROR, OP_ERROR), &
         & "each operation reaches itself by the zero-length path", nfail)
    call report(.not. reachable(view, 7, OP_ERROR), &
         & "an outsider reaches nothing", nfail)

  end subroutine check_reachability

  !===================================================================!
  ! One walk, one order: [predict, error], exactly - execution
  ! order has meaning; operation laws still do not.
  !===================================================================!

  subroutine check_execution_order(nfail)

    integer, intent(inout) :: nfail

    integer, allocatable :: order(:)

    call topological_order(view, order)

    call report(size(order) .eq. 2 .and. &
         &      order(1) .eq. OP_PREDICT .and. order(2) .eq. OP_ERROR, &
         & "the execution order is [predict, error], exactly", nfail)

  end subroutine check_execution_order

end program learning_level_4
