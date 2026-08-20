!=====================================================================!
! LEARNING TOWER . LEVEL 2 . RELATION ALGEBRA
!
! The level answers one question: WHAT CAN BE DERIVED FROM T_flow.
! The answer is the operation dependency - and it is DERIVED, never
! written: restrict the flow to the output port and to the input
! ports, project to O x V and V x O so the middle domain aligns,
! compose - and one pair remains,
!
!      D = { (predict, error) }
!
! because predict PRODUCES yhat and error CONSUMES yhat: the
! structural witness, and nothing numerical anywhere - Level 2
! still does not know that yhat will one day be w*x. D is held as
! class(relation): the learning client leans on the relation
! abstraction, never the storage behind it. Order of writing means
! nothing: the six facts handed backwards derive the same D,
! proved extensionally in both directions.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program learning_level_2

  use learning_assert, only : report, verdict
  use learning_assert, only : SLOT_W, SLOT_X, SLOT_YHAT, SLOT_Y, SLOT_E
  use learning_assert, only : OP_PREDICT, OP_ERROR
  use learning_assert, only : PORT_IN1, PORT_IN2, PORT_OUT
  use graph_fractal        , only : set_graph => graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set        , only : set_map
  use map_inclusion  , only : inclusion_map, declared_subobject
  use relation_finitary , only : stored_relation, relation
  use relation_algebra, only : restrict_slot, project_slots, &
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
  write(*,'(1x,a)') "learning tower . level 2 . relation algebra"
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
         &      t_out3 % has([OP_PREDICT, SLOT_YHAT, PORT_OUT]) .and. &
         &      t_out3 % has([OP_ERROR  , SLOT_E   , PORT_OUT]), &
         & "two tuples pass the output port - exactly these two", nfail)

    call report(t_in3 % num_tuples() .eq. 4 .and. &
         &      t_in3 % has([OP_PREDICT, SLOT_W   , PORT_IN1]) .and. &
         &      t_in3 % has([OP_PREDICT, SLOT_X   , PORT_IN2]) .and. &
         &      t_in3 % has([OP_ERROR  , SLOT_YHAT, PORT_IN1]) .and. &
         &      t_in3 % has([OP_ERROR  , SLOT_Y   , PORT_IN2]), &
         & "four tuples pass the input ports - exactly these four", nfail)

  end subroutine check_restrictions

  !===================================================================!
  ! The two projections: produces <= O x V, consumes <= V x O - the
  ! [2,1] reversal is structural, so the middle domain V aligns for
  ! the composition to come.
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
         &      produces % has([OP_PREDICT, SLOT_YHAT]) .and. &
         &      produces % has([OP_ERROR  , SLOT_E   ]), &
         & "produces = { (predict, yhat), (error, e) } - exactly", nfail)

    dom = consumes % domain(1)
    call report(dom % same_as(v), &
         & "consumes runs from the value slots", nfail)
    dom = consumes % domain(2)
    call report(dom % same_as(o), &
         & "into the operations", nfail)
    call report(consumes % num_tuples() .eq. 4 .and. &
         &      consumes % has([SLOT_W   , OP_PREDICT]) .and. &
         &      consumes % has([SLOT_X   , OP_PREDICT]) .and. &
         &      consumes % has([SLOT_YHAT, OP_ERROR  ]) .and. &
         &      consumes % has([SLOT_Y   , OP_ERROR  ]), &
         & "consumes = { (w,predict), (x,predict), (yhat,error), " // &
         & "(y,error) } - exactly", nfail)

  end subroutine check_projections

  !===================================================================!
  ! The dependency, derived and exact - with its structural witness
  ! stated first: predict produces yhat, error consumes yhat, and
  ! that is the whole reason (predict, error) exists.
  !===================================================================!

  subroutine check_dependency(nfail)

    integer, intent(inout) :: nfail

    type(set_graph) :: dom

    call report(produces % has([OP_PREDICT, SLOT_YHAT]) .and. &
         &      consumes % has([SLOT_YHAT, OP_ERROR]), &
         & "the witness: predict produces yhat and error consumes it", &
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
    call report(d % has([OP_PREDICT, OP_ERROR]), &
         & "D.has(predict, error) = true", nfail)
    call report(.not. d % has([OP_ERROR, OP_PREDICT]), &
         & "D.has(error, predict) = false", nfail)
    call report(.not. d % has([OP_PREDICT, OP_PREDICT]), &
         & "D.has(predict, predict) = false", nfail)
    call report(.not. d % has([OP_ERROR, OP_ERROR]), &
         & "D.has(error, error) = false", nfail)

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

end program learning_level_2
