!=====================================================================!
! DERIVATIVE ACTION TOWER . LEVEL 2 . RELATION ALGEBRA
!
! The level answers one question: WHAT OPERATION DEPENDS ON WHAT
! OPERATION. The answer is DERIVED, never written: restrict the flow
! to the output port and to the input ports, project to O x V and
! V x O so the middle domain aligns, compose - and one pair remains,
!
!      D = { (product, sum) }
!
! because product PRODUCES u and sum CONSUMES u. Operation
! dependency is ordinary relation algebra: nothing derivative-
! specific has appeared, and the derivation is exactly the road the
! calculator and learning towers walked before this one. Order of
! writing means nothing: the six facts handed backwards derive the
! same D, proved extensionally in both directions.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program derivative_level_2

  use derivative_assert, only : report, verdict
  use derivative_assert, only : SLOT_X, SLOT_Y, SLOT_U, SLOT_Z
  use derivative_assert, only : OP_PRODUCT, OP_SUM
  use derivative_assert, only : PORT_IN1, PORT_IN2, PORT_OUT
  use fractal_graph        , only : set_graph => graph
  use graph_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use graph_set_map        , only : set_map
  use graph_inclusion_map  , only : inclusion_map, declared_subobject
  use graph_relation   , only : stored_relation, relation
  use graph_relation_algebra, only : restrict_slot, project_slots, &
       &                             compose_binary

  implicit none

  type(set_graph)            :: v, o, p
  type(set_graph)             :: p_out, p_in
  type(stored_relation)        :: flow, backwards
  type(stored_relation)        :: t_out3, t_in3, produces, consumes
  class(relation), allocatable :: d, d2
  integer                      :: table(3, 6)
  integer                      :: nfail
  type(set_map)     :: sets
  type(inclusion_map)     :: inclusions

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "derivative action tower . level 2 . algebra"
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

  ! The structural selectors: which ports mean leaving, arriving.
  call p_out % declare()
  call sets       % bind(p_out, listed_set_representation([PORT_OUT]))
  call inclusions % include_in(p_out, p)
  call p_in % declare()
  call sets       % bind(p_in, listed_set_representation([PORT_IN1, PORT_IN2]))
  call inclusions % include_in(p_in, p)

  call check_restrictions(nfail)
  call check_projections(nfail)
  call check_dependency(nfail)
  call check_order_invariance(nfail)

  call verdict(nfail, "level 2")

contains

  !===================================================================!
  ! The two restrictions, pinned by extension, not count alone.
  !===================================================================!

  subroutine check_restrictions(nfail)

    integer, intent(inout) :: nfail

    t_out3 = restrict_slot(flow, 3, p_out, sets, inclusions)
    t_in3  = restrict_slot(flow, 3, p_in, sets, inclusions)

    call report(t_out3 % num_tuples() .eq. 2 .and. &
         &      t_out3 % has([OP_PRODUCT, SLOT_U, PORT_OUT]) .and. &
         &      t_out3 % has([OP_SUM    , SLOT_Z, PORT_OUT]), &
         & "two tuples pass the output port - exactly these two", nfail)

    call report(t_in3 % num_tuples() .eq. 4 .and. &
         &      t_in3 % has([OP_PRODUCT, SLOT_X, PORT_IN1]) .and. &
         &      t_in3 % has([OP_PRODUCT, SLOT_Y, PORT_IN2]) .and. &
         &      t_in3 % has([OP_SUM    , SLOT_U, PORT_IN1]) .and. &
         &      t_in3 % has([OP_SUM    , SLOT_Y, PORT_IN2]), &
         & "four tuples pass the input ports - exactly these four", nfail)

  end subroutine check_restrictions

  !===================================================================!
  ! The two projections: produces <= O x V, consumes <= V x O - the
  ! [2,1] reversal is structural, so the middle domain V aligns for
  ! the composition to come. consumes carries y TWICE: once toward
  ! product, once toward sum.
  !===================================================================!

  subroutine check_projections(nfail)

    integer, intent(inout) :: nfail

    type(set_graph) :: dom

    produces = project_slots(t_out3, [1, 2], sets)
    consumes = project_slots(t_in3 , [2, 1], sets)

    dom = produces % domain(1)
    call report(dom % same_as(o), &
         & "produces runs from the operations", nfail)
    dom = produces % domain(2)
    call report(dom % same_as(v), &
         & "into the value slots", nfail)
    call report(produces % num_tuples() .eq. 2 .and. &
         &      produces % has([OP_PRODUCT, SLOT_U]) .and. &
         &      produces % has([OP_SUM    , SLOT_Z]), &
         & "produces = { (product, u), (sum, z) } - exactly", nfail)

    dom = consumes % domain(1)
    call report(dom % same_as(v), &
         & "consumes runs from the value slots", nfail)
    dom = consumes % domain(2)
    call report(dom % same_as(o), &
         & "into the operations", nfail)
    call report(consumes % num_tuples() .eq. 4 .and. &
         &      consumes % has([SLOT_X, OP_PRODUCT]) .and. &
         &      consumes % has([SLOT_Y, OP_PRODUCT]) .and. &
         &      consumes % has([SLOT_U, OP_SUM    ]) .and. &
         &      consumes % has([SLOT_Y, OP_SUM    ]), &
         & "consumes = { (x,product), (y,product), (u,sum), " // &
         & "(y,sum) } - y toward both", nfail)

  end subroutine check_projections

  !===================================================================!
  ! The dependency, derived and exact - with its structural witness
  ! stated first: product produces u, sum consumes u, and that is
  ! the whole reason (product, sum) exists.
  !===================================================================!

  subroutine check_dependency(nfail)

    integer, intent(inout) :: nfail

    type(set_graph) :: dom

    call report(produces % has([OP_PRODUCT, SLOT_U]) .and. &
         &      consumes % has([SLOT_U, OP_SUM]), &
         & "the witness: product produces u and sum consumes it", &
         & nfail)

    d = compose_binary(produces, consumes, sets)

    dom = d % domain(1)
    call report(dom % same_as(o), &
         & "D runs from the operations", nfail)
    dom = d % domain(2)
    call report(dom % same_as(o), &
         & "back into the operations", nfail)

    call report(d % num_tuples() .eq. 1, &
         & "|D| = 1", nfail)
    call report(d % has([OP_PRODUCT, OP_SUM]), &
         & "D.has(product, sum) = true", nfail)
    call report(.not. d % has([OP_SUM, OP_PRODUCT]), &
         & "D.has(sum, product) = false", nfail)
    call report(.not. d % has([OP_PRODUCT, OP_PRODUCT]), &
         & "D.has(product, product) = false", nfail)
    call report(.not. d % has([OP_SUM, OP_SUM]), &
         & "D.has(sum, sum) = false", nfail)

  end subroutine check_dependency

  !===================================================================!
  ! The six facts handed backwards derive the same dependency -
  ! equal as sets, both directions, domains agreeing slot for slot;
  ! never a table comparison.
  !===================================================================!

  subroutine check_order_invariance(nfail)

    integer, intent(inout) :: nfail

    type(stored_relation)          :: out2, in2, prod2, cons2
    type(set_graph) :: da, db
    integer                        :: rev(3, 6), j
    integer, allocatable           :: dt(:,:), dt2(:,:)
    logical                        :: ok

    do j = 1, 6
       rev(:, j) = table(:, 7 - j)
    end do
    backwards = stored_relation('flow backwards', [o, v, p], rev, sets)

    out2  = restrict_slot(backwards, 3, p_out, sets, inclusions)
    in2   = restrict_slot(backwards, 3, p_in, sets, inclusions)
    prod2 = project_slots(out2, [1, 2], sets)
    cons2 = project_slots(in2 , [2, 1], sets)
    d2    = compose_binary(prod2, cons2, sets)

    call report(d2 % num_tuples() .eq. d % num_tuples(), &
         & "|D1| = |D2|", nfail)

    call d % tuples(dt)
    call d2 % tuples(dt2)
    ok = .true.
    do j = 1, size(dt, 2)
       ok = ok .and. d2 % has(dt(:, j))
    end do
    do j = 1, size(dt2, 2)
       ok = ok .and. d % has(dt2(:, j))
    end do
    call report(ok, &
         & "each holds every tuple of the other: equal as sets", nfail)

    da = d % domain(1)
    db = d2 % domain(1)
    call report(da % same_as(db), &
         & "the first slots are one domain", nfail)
    da = d % domain(2)
    db = d2 % domain(2)
    call report(da % same_as(db), &
         & "and so are the second", nfail)

  end subroutine check_order_invariance

end program derivative_level_2
