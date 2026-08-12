!=====================================================================!
! CALCULATOR TOWER . LEVEL 2 . RELATION ALGEBRA
!
! The level answers one question: WHAT NEW RELATIONS CAN BE DERIVED.
! The dependency of the calculator's two operations is not stored
! anywhere - it is DERIVED from the six flow tuples, by the three
! primitives this level earned:
!
!      R_out3   = restrict R_flow to the output port
!      R_in3    = restrict R_flow to the input ports
!      produces = project R_out3 onto O x X
!      consumes = project R_in3  onto X x O
!      D        = consumes o produces  <=  O x O
!
! and the whole answer is one tuple:  D = { (+, x) }. The pair
! (+, x) is never written down on the derivation path - that is the
! point of the level. Order of writing means nothing: the same
! extension handed backwards derives the same D, proved
! extensionally, never by comparing enumeration tables.
!
! Still true, still forbidden: no values on the slots, no
! arithmetic on the operations, no graph anywhere. D says only that
! + structurally precedes x.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program calculator_level_2

  use calculator_assert, only : report, verdict
  use calculator_assert, only : SLOT_A, SLOT_B, SLOT_C, SLOT_D, SLOT_E
  use calculator_assert, only : OP_PLUS, OP_TIMES
  use calculator_assert, only : PORT_IN1, PORT_IN2, PORT_OUT
  use graph_carrier    , only : counted_set, subset_set, member_set
  use graph_relation   , only : stored_relation, relation
  use graph_relation_algebra, only : restrict_slot, project_slots, &
       &                             compose_binary

  implicit none

  type(counted_set)            :: x, o, p
  type(subset_set)             :: p_out, p_in
  type(stored_relation)        :: flow, backwards
  type(stored_relation)        :: r_out3, r_in3, produces, consumes
  class(relation), allocatable :: d, d2
  integer                      :: table(3, 6)
  integer                      :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "calculator tower . level 2 . relation algebra"
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

  ! The structural selections: which ports mean leaving, which mean
  ! arriving. Subobjects of P - not fields, not values.
  p_out = subset_set('output-port', p, [PORT_OUT])
  p_in  = subset_set('input-ports', p, [PORT_IN1, PORT_IN2])

  call check_restrictions(nfail)
  call check_projections(nfail)
  call check_collapse(nfail)
  call check_dependency(nfail)
  call check_order_invariance(nfail)

  call verdict(nfail, "level 2")

contains

  !===================================================================!
  ! The two restrictions: two tuples leave, four arrive.
  !===================================================================!

  subroutine check_restrictions(nfail)

    integer, intent(inout) :: nfail

    r_out3 = restrict_slot(flow, 3, p_out)
    r_in3  = restrict_slot(flow, 3, p_in)

    call report(r_out3 % num_tuples() .eq. 2, &
         & "two tuples pass the output port", nfail)
    call report(r_in3 % num_tuples() .eq. 4, &
         & "four tuples pass the input ports", nfail)

  end subroutine check_restrictions

  !===================================================================!
  ! The two projections: produces <= O x X reads who writes which
  ! slot; consumes <= X x O reads which slot feeds whom - the slot
  ! order [2,1] reversed on purpose, so the middle domain lines up.
  !===================================================================!

  subroutine check_projections(nfail)

    integer, intent(inout) :: nfail

    class(member_set), allocatable :: dom

    produces = project_slots(r_out3, [1, 2])
    consumes = project_slots(r_in3 , [2, 1])

    dom = produces % domain(1)
    call report(dom % same_as(o), &
         & "produces runs from the operations", nfail)
    dom = produces % domain(2)
    call report(dom % same_as(x), &
         & "into the value slots", nfail)
    call report(produces % num_tuples() .eq. 2 .and. &
         &      produces % has([OP_PLUS , SLOT_C]) .and. &
         &      produces % has([OP_TIMES, SLOT_E]), &
         & "produces = { (+, c), (x, e) } - exactly", nfail)

    dom = consumes % domain(1)
    call report(dom % same_as(x), &
         & "consumes runs from the value slots", nfail)
    dom = consumes % domain(2)
    call report(dom % same_as(o), &
         & "into the operations", nfail)
    call report(consumes % num_tuples() .eq. 4 .and. &
         &      consumes % has([SLOT_A, OP_PLUS ]) .and. &
         &      consumes % has([SLOT_B, OP_PLUS ]) .and. &
         &      consumes % has([SLOT_C, OP_TIMES]) .and. &
         &      consumes % has([SLOT_D, OP_TIMES]), &
         & "consumes = { (a,+), (b,+), (c,x), (d,x) } - exactly", nfail)

  end subroutine check_projections

  !===================================================================!
  ! Projection genuinely collapses: the four input tuples name only
  ! two operations, and the image onto O alone holds two members,
  ! not four entries.
  !===================================================================!

  subroutine check_collapse(nfail)

    integer, intent(inout) :: nfail

    type(stored_relation) :: ops_used

    ops_used = project_slots(r_in3, [1])

    call report(ops_used % num_tuples() .eq. 2 .and. &
         &      ops_used % has([OP_PLUS]) .and. &
         &      ops_used % has([OP_TIMES]), &
         & "four input tuples project to two operations, not four", nfail)

  end subroutine check_collapse

  !===================================================================!
  ! The dependency, derived and exact: D = consumes o produces,
  ! one tuple, (+, x), both slots the operations - and the pair was
  ! never written down on the way here.
  !===================================================================!

  subroutine check_dependency(nfail)

    integer, intent(inout) :: nfail

    class(member_set), allocatable :: dom

    d = compose_binary(produces, consumes)

    dom = d % domain(1)
    call report(dom % same_as(o), &
         & "D runs from the operations", nfail)
    dom = d % domain(2)
    call report(dom % same_as(o), &
         & "back into the operations", nfail)

    call report(d % num_tuples() .eq. 1, &
         & "|D| = 1", nfail)
    call report(d % has([OP_PLUS, OP_TIMES]), &
         & "D.has(+, x) = true", nfail)
    call report(.not. d % has([OP_TIMES, OP_PLUS]), &
         & "D.has(x, +) = false", nfail)
    call report(.not. d % has([OP_PLUS, OP_PLUS]), &
         & "D.has(+, +) = false", nfail)
    call report(.not. d % has([OP_TIMES, OP_TIMES]), &
         & "D.has(x, x) = false", nfail)

  end subroutine check_dependency

  !===================================================================!
  ! The set does not remember how it was written down: the same six
  ! facts handed backwards derive the same dependency, proved
  ! extensionally - equal count, every tuple of one in the other,
  ! domains agreeing slot for slot - never by comparing tables.
  !===================================================================!

  subroutine check_order_invariance(nfail)

    integer, intent(inout) :: nfail

    type(stored_relation)          :: out2, in2, prod2, cons2
    class(member_set), allocatable :: da, db
    integer                        :: rev(3, 6)
    integer, allocatable           :: dt(:,:)
    integer                        :: j
    logical                        :: ok

    do j = 1, 6
       rev(:, j) = table(:, 7 - j)
    end do
    backwards = stored_relation('flow backwards', [o, x, p], rev)

    out2  = restrict_slot(backwards, 3, p_out)
    in2   = restrict_slot(backwards, 3, p_in)
    prod2 = project_slots(out2, [1, 2])
    cons2 = project_slots(in2 , [2, 1])
    d2    = compose_binary(prod2, cons2)

    call report(d2 % num_tuples() .eq. d % num_tuples(), &
         & "|D1| = |D2|", nfail)

    call d % tuples(dt)
    ok = .true.
    do j = 1, size(dt, 2)
       ok = ok .and. d2 % has(dt(:, j))
    end do
    call report(ok, &
         & "every tuple of D1 stands in D2: equal as sets", nfail)

    da = d % domain(1)
    db = d2 % domain(1)
    call report(da % same_as(db), &
         & "the first slots are one domain", nfail)
    da = d % domain(2)
    db = d2 % domain(2)
    call report(da % same_as(db), &
         & "and so are the second", nfail)

  end subroutine check_order_invariance

end program calculator_level_2
