!=====================================================================!
! CALCULATOR TOWER . LEVEL 1 . THE RELATION
!
! The level answers one question: HOW MAY MEMBERS BE RELATED. The
! three carriers stand as at Level 0, and one ternary relation
! joins them,
!
!      T_flow  <=  O x X x P
!
!      (+, a, in1)     (+, b, in2)     (+, c, out)
!      (x, c, in1)     (x, d, in2)     (x, e, out)
!
! Six tuples, exactly - handed to the constructor SEVEN times, one
! duplicated on purpose, because a relation is a set and this level
! proves it independently. The + here is only a member of O; the
! slots carry no values; no graph reads these tuples yet. The
! import list IS the negative truth: calculator_assert,
! graph_carrier, graph_relation, and nothing above (CALCULATOR.md
! 8; the ternary shape is the point - no binary reduction, no
! vertex, no edge).
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program calculator_level_1

  use calculator_assert, only : report, verdict
  use calculator_assert, only : SLOT_A, SLOT_B, SLOT_C, SLOT_D, SLOT_E
  use calculator_assert, only : OP_PLUS, OP_TIMES
  use calculator_assert, only : PORT_IN1, PORT_IN2, PORT_OUT
  use graph_carrier    , only : counted_set, member_set
  use graph_relation   , only : stored_relation

  implicit none

  type(counted_set)     :: x, o, p
  type(stored_relation) :: flow
  integer               :: table(3, 7)
  integer               :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "calculator tower . level 1 . relation"
  write(*,'(1x,a)') "============================================="

  x = counted_set('value-slots', 5)
  o = counted_set('operations' , 2)
  p = counted_set('ports'      , 3)

  ! The six facts of the flow - and the first of them handed twice.
  table(:, 1) = [OP_PLUS , SLOT_A, PORT_IN1]
  table(:, 2) = [OP_PLUS , SLOT_B, PORT_IN2]
  table(:, 3) = [OP_PLUS , SLOT_C, PORT_OUT]
  table(:, 4) = [OP_TIMES, SLOT_C, PORT_IN1]
  table(:, 5) = [OP_TIMES, SLOT_D, PORT_IN2]
  table(:, 6) = [OP_TIMES, SLOT_E, PORT_OUT]
  table(:, 7) = [OP_PLUS , SLOT_A, PORT_IN1]

  flow = stored_relation('flow', [o, x, p], table)

  call check_signature(nfail)
  call check_complete_extension(nfail)
  call check_representative_membership(nfail)

  call verdict(nfail, "level 1")

contains

  !===================================================================!
  ! Arity three, and the ordered signature answers the three
  ! declared domains by identity: O, then X, then P.
  !===================================================================!

  subroutine check_signature(nfail)

    integer, intent(inout) :: nfail

    class(member_set), allocatable :: d

    call report(flow % arity() .eq. 3, &
         & "the flow is genuinely ternary", nfail)

    d = flow % domain(1)
    call report(d % same_as(o), &
         & "slot one is the operations", nfail)
    d = flow % domain(2)
    call report(d % same_as(x), &
         & "slot two is the value slots", nfail)
    d = flow % domain(3)
    call report(d % same_as(p), &
         & "slot three is the ports", nfail)

  end subroutine check_signature

  !===================================================================!
  ! The complete extension: |R| = 6 although seven tuples were
  ! handed in - the duplicate collapsed, for a relation is a set -
  ! and every one of the six expected tuples is a member. Six
  ! members present in a six-element set: no seventh exists.
  !===================================================================!

  subroutine check_complete_extension(nfail)

    integer, intent(inout) :: nfail

    call report(flow % num_tuples() .eq. 6, &
         & "seven handed, six held: a relation is a set", nfail)

    call report(flow % has([OP_PLUS , SLOT_A, PORT_IN1]) .and. &
         &      flow % has([OP_PLUS , SLOT_B, PORT_IN2]) .and. &
         &      flow % has([OP_PLUS , SLOT_C, PORT_OUT]) .and. &
         &      flow % has([OP_TIMES, SLOT_C, PORT_IN1]) .and. &
         &      flow % has([OP_TIMES, SLOT_D, PORT_IN2]) .and. &
         &      flow % has([OP_TIMES, SLOT_E, PORT_OUT]), &
         & "all six expected tuples are members - the extension is exact", &
         & nfail)

  end subroutine check_complete_extension

  !===================================================================!
  ! The representative truths of CALCULATOR.md 8, stated in their
  ! own words, plus the absences that keep membership honest.
  !===================================================================!

  subroutine check_representative_membership(nfail)

    integer, intent(inout) :: nfail

    call report(flow % has([OP_PLUS, SLOT_A, PORT_IN1]), &
         & "R.has(+, a, in1) = true", nfail)
    call report(flow % has([OP_PLUS, SLOT_C, PORT_OUT]), &
         & "R.has(+, c, out) = true", nfail)
    call report(flow % has([OP_TIMES, SLOT_C, PORT_IN1]), &
         & "R.has(x, c, in1) = true", nfail)
    call report(.not. flow % has([OP_TIMES, SLOT_A, PORT_IN1]), &
         & "R.has(x, a, in1) = false", nfail)

    call report(.not. flow % has([OP_PLUS, SLOT_A]), &
         & "a tuple of the wrong length belongs to nothing", nfail)

  end subroutine check_representative_membership

end program calculator_level_1
