!=====================================================================!
! DERIVATIVE ACTION TOWER . LEVEL 1 . THE RELATION
!
! The level answers one question: HOW IS THE SYMBOLIC COMPUTATION
! STRUCTURALLY WIRED. The three sets stand as at Level 0, and
! one ternary relation joins them,
!
!      R_flow  <=  O x V x P
!
!      (product, x, in1)   (product, y, in2)   (product, u, out)
!      (sum,     u, in1)   (sum,     y, in2)   (sum,     z, out)
!
! Six tuples, exactly - handed to the constructor SEVEN times, one
! duplicated on purpose, because a relation is a set and this tower
! proves it independently of its siblings. The load-bearing wiring
! fact: y is consumed TWICE - by product and by sum - the fan-out
! that will one day make reverse accumulation nontrivial. But at
! this level derivative potential is completely LATENT: R_flow
! contains computation structure, not derivative metadata, and it
! looks exactly like ordinary computation structure because it IS
! ordinary computation structure. product does not yet multiply;
! sum does not yet add; nothing differentiates anything. The import
! list is the negative truth - derivative_assert, graph_set,
! graph_relation, and nothing above.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program derivative_level_1

  use derivative_assert, only : report, verdict
  use derivative_assert, only : SLOT_X, SLOT_Y, SLOT_U, SLOT_Z
  use derivative_assert, only : OP_PRODUCT, OP_SUM
  use derivative_assert, only : PORT_IN1, PORT_IN2, PORT_OUT
  use graph_set    , only : index_set, set
  use graph_relation   , only : stored_relation

  implicit none

  type(index_set)     :: v, o, p
  type(stored_relation) :: flow
  integer               :: table(3, 7)
  integer               :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "derivative action tower . level 1 . relation"
  write(*,'(1x,a)') "============================================="

  v = index_set('value-slots', 4)
  o = index_set('operations' , 2)
  p = index_set('ports'      , 3)

  ! The six facts of the specimen's flow - and the first handed twice.
  table(:, 1) = [OP_PRODUCT, SLOT_X, PORT_IN1]
  table(:, 2) = [OP_PRODUCT, SLOT_Y, PORT_IN2]
  table(:, 3) = [OP_PRODUCT, SLOT_U, PORT_OUT]
  table(:, 4) = [OP_SUM    , SLOT_U, PORT_IN1]
  table(:, 5) = [OP_SUM    , SLOT_Y, PORT_IN2]
  table(:, 6) = [OP_SUM    , SLOT_Z, PORT_OUT]
  table(:, 7) = [OP_PRODUCT, SLOT_X, PORT_IN1]

  flow = stored_relation('flow', [o, v, p], table)

  call check_signature(nfail)
  call check_complete_extension(nfail)
  call check_fan_out(nfail)
  call check_meaningful_absences(nfail)

  call verdict(nfail, "level 1")

contains

  !===================================================================!
  ! Arity three, and the ordered signature answers the three
  ! DECLARED domains by identity - equal sizes buy nothing.
  !===================================================================!

  subroutine check_signature(nfail)

    integer, intent(inout) :: nfail

    class(set), allocatable :: d

    call report(flow % arity() .eq. 3, &
         & "the flow is genuinely ternary", nfail)

    d = flow % domain(1)
    call report(d % equals(o), &
         & "slot one is the operations, by identity", nfail)
    d = flow % domain(2)
    call report(d % equals(v), &
         & "slot two is the value slots", nfail)
    d = flow % domain(3)
    call report(d % equals(p), &
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

    call report(flow % has([OP_PRODUCT, SLOT_X, PORT_IN1]) .and. &
         &      flow % has([OP_PRODUCT, SLOT_Y, PORT_IN2]) .and. &
         &      flow % has([OP_PRODUCT, SLOT_U, PORT_OUT]) .and. &
         &      flow % has([OP_SUM    , SLOT_U, PORT_IN1]) .and. &
         &      flow % has([OP_SUM    , SLOT_Y, PORT_IN2]) .and. &
         &      flow % has([OP_SUM    , SLOT_Z, PORT_OUT]), &
         & "all six expected tuples are members - the extension is exact", &
         & nfail)

  end subroutine check_complete_extension

  !===================================================================!
  ! The wiring fact this specimen was chosen for: y feeds BOTH
  ! operations. Two consumers of one slot - stated here as plain
  ! structure, and worth nothing derivative yet.
  !===================================================================!

  subroutine check_fan_out(nfail)

    integer, intent(inout) :: nfail

    call report(flow % has([OP_PRODUCT, SLOT_Y, PORT_IN2]) .and. &
         &      flow % has([OP_SUM    , SLOT_Y, PORT_IN2]), &
         & "y is consumed twice - by product and by sum", nfail)
    call report(flow % has([OP_PRODUCT, SLOT_X, PORT_IN1]) .and. &
         &      .not. flow % has([OP_SUM, SLOT_X, PORT_IN1]) .and. &
         &      .not. flow % has([OP_SUM, SLOT_X, PORT_IN2]), &
         & "x is consumed once - by product alone", nfail)

  end subroutine check_fan_out

  !===================================================================!
  ! Structurally meaningful absences: no operation consumes its own
  ! output, z feeds nothing, and a wrong-length tuple belongs to
  ! nothing.
  !===================================================================!

  subroutine check_meaningful_absences(nfail)

    integer, intent(inout) :: nfail

    call report(.not. flow % has([OP_PRODUCT, SLOT_U, PORT_IN1]) .and. &
         &      .not. flow % has([OP_PRODUCT, SLOT_U, PORT_IN2]), &
         & "product does not consume its own output", nfail)
    call report(.not. flow % has([OP_PRODUCT, SLOT_Z, PORT_OUT]), &
         & "z is not product's output", nfail)
    call report(.not. flow % has([OP_SUM, SLOT_Z, PORT_IN1]) .and. &
         &      .not. flow % has([OP_PRODUCT, SLOT_Z, PORT_IN1]), &
         & "the response feeds nothing", nfail)

    call report(.not. flow % has([OP_PRODUCT, SLOT_X]), &
         & "a tuple of the wrong length belongs to nothing", nfail)

  end subroutine check_meaningful_absences

end program derivative_level_1
