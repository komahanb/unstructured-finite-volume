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
  use graph_fractal        , only : set_graph => graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set        , only : set_map
  use map_label      , only : label_map
  use map_inclusion  , only : inclusion_map, declared_subobject
  use relation_finitary   , only : stored_relation, relation
  use relation_binary, only : binary_relation
  use relation_algebra, only : restrict_slot, project_slots, &
       &                             compose_binary
  use graph_algorithms , only : sources, sinks, reachable, &
       &                        topological_order
  use graph_fractal        , only : graph, known_branch, null_branch
  use graph_relational_view, only : relational_binding, &
       & num_member_sets, member_set_at, num_relations, relation_at, &
       & holds_set

  implicit none


  type(set_graph)              :: v, o, p
  type(set_graph)               :: p_out, p_in
  type(stored_relation)          :: flow
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
  write(*,'(1x,a)') "derivative action tower . level 4 . calculus"
  write(*,'(1x,a)') "============================================="

  call v % declare()
  call sets % bind(v, counted_set_representation(4))
  call o % declare()
  call sets % bind(o, counted_set_representation(2))
  call p % declare()
  call sets % bind(p, counted_set_representation(3))

  table(:, 1) = [OP_PRODUCT, SLOT_X, PORT_IN1]
  table(:, 2) = [OP_PRODUCT, SLOT_Y, PORT_IN2]
  table(:, 3) = [OP_PRODUCT, SLOT_U, PORT_OUT]
  table(:, 4) = [OP_SUM    , SLOT_U, PORT_IN1]
  table(:, 5) = [OP_SUM    , SLOT_Y, PORT_IN2]
  table(:, 6) = [OP_SUM    , SLOT_Z, PORT_OUT]
  flow = stored_relation('flow', [o, v, p], table, sets)

  call p_out % declare()
  call sets       % bind(p_out, listed_set_representation([PORT_OUT]))
  call inclusions % include_in(p_out, p)
  call p_in % declare()
  call sets       % bind(p_in, listed_set_representation([PORT_IN1, PORT_IN2]))
  call inclusions % include_in(p_in, p)
  d = compose_binary( &
       & project_slots(restrict_slot(flow, 3, p_out, sets, inclusions), [1, 2], sets), &
       & project_slots(restrict_slot(flow, 3, p_in , sets, inclusions), [2, 1], sets), sets)

  ! 'derivative specimen': (S, P) as one sequence on each branch.
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

  ! The interpretive jump, made explicitly.

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
         & "the relation walks the operations, and nothing invented", nfail)

  end subroutine check_view_domain

  !===================================================================!
  ! Sources and sinks - subobjects of O.
  !===================================================================!

  subroutine check_sources_and_sinks(nfail)

    integer, intent(inout) :: nfail

    type(set_graph) :: src, snk
    type(label_map)     :: labels

    call sources(d, sets, labels, inclusions, src)
    call sinks(d, sets, labels, inclusions, snk)

    call report(sets % size_of(src) .eq. 1 .and. sets % has_in(src, OP_PRODUCT) .and. &
         &      .not. sets % has_in(src, OP_SUM), &
         & "sources = { product }", nfail)
    call report(sets % size_of(snk) .eq. 1 .and. sets % has_in(snk, OP_SUM) .and. &
         &      .not. sets % has_in(snk, OP_PRODUCT), &
         & "sinks = { sum }", nfail)

    call report(declared_subobject(src, o, inclusions) .and. declared_subobject(snk, o, inclusions), &
         & "both answers stand embedded in the operations", nfail)

  end subroutine check_sources_and_sinks

  !===================================================================!
  ! Directional dependency: product reaches sum, never the reverse -
  ! and each reaches itself by the zero-length path.
  !===================================================================!

  subroutine check_reachability(nfail)

    integer, intent(inout) :: nfail

    call report(reachable(d, sets, OP_PRODUCT, OP_SUM), &
         & "reachable(product, sets, sum) = true", nfail)
    call report(.not. reachable(d, sets, OP_SUM, OP_PRODUCT), &
         & "reachable(sum, sets, product) = false", nfail)
    call report(reachable(d, sets, OP_PRODUCT, OP_PRODUCT) .and. &
         &      reachable(d, sets, OP_SUM, OP_SUM), &
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

    call topological_order(d, sets, order)

    call report(size(order) .eq. 2 .and. &
         &      order(1) .eq. OP_PRODUCT .and. order(2) .eq. OP_SUM, &
         & "the execution order is [product, sum], exactly", nfail)

  end subroutine check_execution_order

end program derivative_level_4
