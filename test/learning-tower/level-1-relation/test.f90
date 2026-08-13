!=====================================================================!
! LEARNING TOWER . LEVEL 1 . THE RELATION
!
! The level answers one question: HOW ARE THE SYMBOLIC SLOTS
! STRUCTURALLY RELATED. The three carriers stand as at Level 0, and
! one ternary relation joins them,
!
!      R_flow  <=  O x V x P
!
!      (predict, w,    in1)     (predict, x, in2)   (predict, yhat, out)
!      (error,   yhat, in1)     (error,   y, in2)   (error,   e,    out)
!
! Six tuples, exactly - handed to the constructor SEVEN times, one
! duplicated on purpose, because a relation is a set and this tower
! proves it independently of the calculator. The relation states
! only that predict consumes w and x and produces yhat, and that
! error consumes yhat and y and produces e: predict does NOT yet
! mean w*x, error does NOT yet mean yhat-y - structure before
! meaning, the learning tower's first real proof of it. No neuron,
! no layer, no network object anywhere: the import list IS the
! negative truth - learning_assert, graph_carrier, graph_relation,
! and nothing above.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program learning_level_1

  use learning_assert, only : report, verdict
  use learning_assert, only : SLOT_W, SLOT_X, SLOT_YHAT, SLOT_Y, SLOT_E
  use learning_assert, only : OP_PREDICT, OP_ERROR
  use learning_assert, only : PORT_IN1, PORT_IN2, PORT_OUT
  use graph_carrier  , only : counted_set, member_set
  use graph_relation , only : stored_relation

  implicit none

  type(counted_set)     :: v, o, p
  type(stored_relation) :: flow
  integer               :: table(3, 7)
  integer               :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "learning tower . level 1 . relation"
  write(*,'(1x,a)') "============================================="

  v = counted_set('value-slots', 5)
  o = counted_set('operations' , 2)
  p = counted_set('ports'      , 3)

  ! The six facts of the model's flow - and the first handed twice.
  table(:, 1) = [OP_PREDICT, SLOT_W   , PORT_IN1]
  table(:, 2) = [OP_PREDICT, SLOT_X   , PORT_IN2]
  table(:, 3) = [OP_PREDICT, SLOT_YHAT, PORT_OUT]
  table(:, 4) = [OP_ERROR  , SLOT_YHAT, PORT_IN1]
  table(:, 5) = [OP_ERROR  , SLOT_Y   , PORT_IN2]
  table(:, 6) = [OP_ERROR  , SLOT_E   , PORT_OUT]
  table(:, 7) = [OP_PREDICT, SLOT_W   , PORT_IN1]

  flow = stored_relation('flow', [o, v, p], table)

  call check_signature(nfail)
  call check_complete_extension(nfail)
  call check_meaningful_absences(nfail)

  call verdict(nfail, "level 1")

contains

  !===================================================================!
  ! Arity three, and the ordered signature answers the three
  ! DECLARED domains by identity - equal sizes buy nothing.
  !===================================================================!

  subroutine check_signature(nfail)

    integer, intent(inout) :: nfail

    class(member_set), allocatable :: d

    call report(flow % arity() .eq. 3, &
         & "the flow is genuinely ternary", nfail)

    d = flow % domain(1)
    call report(d % same_as(o), &
         & "slot one is the operations, by identity", nfail)
    d = flow % domain(2)
    call report(d % same_as(v), &
         & "slot two is the value slots", nfail)
    d = flow % domain(3)
    call report(d % same_as(p), &
         & "slot three is the ports", nfail)

  end subroutine check_signature

  !===================================================================!
  ! The complete extension: |R| = 6 although seven tuples were
  ! handed in - the duplicate collapsed - and every one of the six
  ! expected tuples is a member. Six present in a six-element set:
  ! no seventh exists.
  !===================================================================!

  subroutine check_complete_extension(nfail)

    integer, intent(inout) :: nfail

    call report(flow % num_tuples() .eq. 6, &
         & "seven handed, six held: a relation is a set", nfail)

    call report(flow % has([OP_PREDICT, SLOT_W   , PORT_IN1]) .and. &
         &      flow % has([OP_PREDICT, SLOT_X   , PORT_IN2]) .and. &
         &      flow % has([OP_PREDICT, SLOT_YHAT, PORT_OUT]) .and. &
         &      flow % has([OP_ERROR  , SLOT_YHAT, PORT_IN1]) .and. &
         &      flow % has([OP_ERROR  , SLOT_Y   , PORT_IN2]) .and. &
         &      flow % has([OP_ERROR  , SLOT_E   , PORT_OUT]), &
         & "all six expected tuples are members - the extension is exact", &
         & nfail)

  end subroutine check_complete_extension

  !===================================================================!
  ! Structurally meaningful absences: the target never feeds
  ! predict, the parameter never feeds error, and no operation
  ! consumes its own output slot on an input port.
  !===================================================================!

  subroutine check_meaningful_absences(nfail)

    integer, intent(inout) :: nfail

    call report(.not. flow % has([OP_PREDICT, SLOT_Y, PORT_IN1]), &
         & "the target y never feeds predict", nfail)
    call report(.not. flow % has([OP_ERROR, SLOT_W, PORT_IN1]), &
         & "the parameter w never feeds error", nfail)
    call report(.not. flow % has([OP_ERROR, SLOT_E, PORT_IN1]), &
         & "error does not consume its own output", nfail)

    call report(.not. flow % has([OP_PREDICT, SLOT_W]), &
         & "a tuple of the wrong length belongs to nothing", nfail)

  end subroutine check_meaningful_absences

end program learning_level_1
