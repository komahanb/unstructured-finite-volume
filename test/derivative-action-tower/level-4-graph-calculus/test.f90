!=====================================================================!
! DERIVATIVE ACTION TOWER . LEVEL 4 . GRAPH CALCULUS
!
! The level answers one question: WHEN THE OPERATION RELATION IS
! INTERPRETED DIRECTIONALLY, WHAT EXECUTION STRUCTURE APPEARS. The
! relation product->sum has existed since Level 2; only THIS rung
! chooses to read it as directed dependency - interpretation is
! explicit, never automatic:
!
!      graph-owned D --interpreted--> directed adjacency over O
!      sources = {product}   sinks = {sum}   walk = [product, sum]
!
! And the three concepts stay three:
!
!      dependency relation  /=  directed interpretation
!                           /=  derivative traversal
!
! This walk is computation order and nothing more. No tangent
! travels it, no cotangent travels it backward, and the words
! forward mode and reverse mode attach to NOTHING here - a future
! gate must earn them over finer structure than this order.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program derivative_level_4

  use derivative_assert, only : report, verdict
  use derivative_assert, only : SLOT_X, SLOT_Y, SLOT_U, SLOT_Z
  use derivative_assert, only : OP_PRODUCT, OP_SUM
  use derivative_assert, only : PORT_IN1, PORT_IN2, PORT_OUT
  use graph_set    , only : index_set, subset, set
  use graph_relation   , only : stored_relation, relation
  use graph_relation_algebra, only : restrict_slot, project_slots, &
       &                             compose_binary
  use graph_structure  , only : related_graph, declared_set, declared_relation
  use graph_interpretation    , only : directed_adjacency_view
  use graph_algorithms , only : sources, sinks, reachable, &
       &                        topological_order

  implicit none

  type(index_set)              :: v, o, p
  type(subset)               :: p_out, p_in
  type(stored_relation)          :: flow
  class(relation), allocatable   :: d
  type(related_graph), target :: g
  type(directed_adjacency_view)  :: view
  integer                        :: table(3, 6)
  integer                        :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "derivative action tower . level 4 . calculus"
  write(*,'(1x,a)') "============================================="

  v = index_set('value-slots', 4)
  o = index_set('operations' , 2)
  p = index_set('ports'      , 3)

  table(:, 1) = [OP_PRODUCT, SLOT_X, PORT_IN1]
  table(:, 2) = [OP_PRODUCT, SLOT_Y, PORT_IN2]
  table(:, 3) = [OP_PRODUCT, SLOT_U, PORT_OUT]
  table(:, 4) = [OP_SUM    , SLOT_U, PORT_IN1]
  table(:, 5) = [OP_SUM    , SLOT_Y, PORT_IN2]
  table(:, 6) = [OP_SUM    , SLOT_Z, PORT_OUT]
  flow = stored_relation('flow', [o, v, p], table)

  p_out = subset('output-port', p, [PORT_OUT])
  p_in  = subset('input-ports', p, [PORT_IN1, PORT_IN2])
  d = compose_binary( &
       & project_slots(restrict_slot(flow, 3, p_out), [1, 2]), &
       & project_slots(restrict_slot(flow, 3, p_in ), [2, 1]))

  g = related_graph('derivative specimen', &
       & [declared_set(v), declared_set(o), declared_set(p)], &
       & [declared_relation(flow), declared_relation(d)])

  ! The interpretive jump, made explicitly.
  view = directed_adjacency_view(g, d)

  call check_view_domain(nfail)
  call check_sources_and_sinks(nfail)
  call check_reachability(nfail)
  call check_execution_order(nfail)

  call verdict(nfail, "level 4")

contains

  !===================================================================!
  ! The interpretation runs over the operations themselves - no new
  ! vertex set, no manufactured edge members.
  !===================================================================!

  subroutine check_view_domain(nfail)

    integer, intent(inout) :: nfail

    class(set), allocatable :: dom

    dom = view % domain()
    call report(dom % equals(o) .and. dom % size() .eq. 2, &
         & "the view walks the operations, and nothing invented", nfail)

  end subroutine check_view_domain

  !===================================================================!
  ! Sources and sinks - subobjects of O.
  !===================================================================!

  subroutine check_sources_and_sinks(nfail)

    integer, intent(inout) :: nfail

    type(subset) :: src, snk

    src = sources(view)
    snk = sinks(view)

    call report(src % size() .eq. 1 .and. src % has(OP_PRODUCT) .and. &
         &      .not. src % has(OP_SUM), &
         & "sources = { product }", nfail)
    call report(snk % size() .eq. 1 .and. snk % has(OP_SUM) .and. &
         &      .not. snk % has(OP_PRODUCT), &
         & "sinks = { sum }", nfail)

    call report(src % is_subobject_of(o) .and. snk % is_subobject_of(o), &
         & "both answers stand embedded in the operations", nfail)

  end subroutine check_sources_and_sinks

  !===================================================================!
  ! Directional dependency: product reaches sum, never the reverse -
  ! and each reaches itself by the zero-length path.
  !===================================================================!

  subroutine check_reachability(nfail)

    integer, intent(inout) :: nfail

    call report(reachable(view, OP_PRODUCT, OP_SUM), &
         & "reachable(product, sum) = true", nfail)
    call report(.not. reachable(view, OP_SUM, OP_PRODUCT), &
         & "reachable(sum, product) = false", nfail)
    call report(reachable(view, OP_PRODUCT, OP_PRODUCT) .and. &
         &      reachable(view, OP_SUM, OP_SUM), &
         & "each operation reaches itself by the zero-length path", nfail)

  end subroutine check_reachability

  !===================================================================!
  ! One walk, one order: [product, sum], exactly - computation
  ! order, and NOT derivative traversal: no tangent rides it, no
  ! cotangent rides it backward.
  !===================================================================!

  subroutine check_execution_order(nfail)

    integer, intent(inout) :: nfail

    integer, allocatable :: order(:)

    call topological_order(view, order)

    call report(size(order) .eq. 2 .and. &
         &      order(1) .eq. OP_PRODUCT .and. order(2) .eq. OP_SUM, &
         & "the execution order is [product, sum], exactly", nfail)

  end subroutine check_execution_order

end program derivative_level_4
