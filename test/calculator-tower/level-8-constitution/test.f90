!=====================================================================!
! CALCULATOR TOWER . LEVEL 8 . CONSTITUTION
!
! The level answers one question: WHAT LAWS DO THE SYMBOLS OBEY.
! Only here, in one small table, + becomes addition and x becomes
! multiplication - the table is TEST-LOCAL, because one calculator
! has not earned a universal production constitution. Everything
! else was already true below:
!
!      R_flow chooses each operation's in1, in2 and out slots
!      L locates each residual row at its computed slot
!      K = { d, a, b } and U = { e, c } hold the values
!      the law supplies only the arithmetic
!
! and the generated residual on Y reproduces the very system Level
! 7 solved - proved at several off-solution states, never only at
! the root - while the constitution-generated dependency equals the
! Level-6 structural supports exactly: meaning added, topology
! untouched. The evaluator names no slot and no row: structure
! chooses the operands; the constitution chooses the law.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module arithmetic_constitution_fixture

  ! The one place where the symbols mean something.

  use iso_fortran_env  , only : dp => REAL64
  use calculator_assert, only : OP_PLUS, OP_TIMES

  implicit none

  private
  public :: apply_law

contains

  real(dp) function apply_law(op, x, y) result(z)

    integer , intent(in) :: op
    real(dp), intent(in) :: x, y

    select case (op)
    case (OP_PLUS)
       z = x + y
    case (OP_TIMES)
       z = x * y
    case default
       error stop 'constitution: no law binds this operation symbol'
    end select

  end function apply_law

end module arithmetic_constitution_fixture

program calculator_level_8

  use iso_fortran_env  , only : dp => REAL64
  use calculator_assert, only : report, verdict
  use calculator_assert, only : SLOT_A, SLOT_B, SLOT_C, SLOT_D, SLOT_E
  use calculator_assert, only : OP_PLUS, OP_TIMES
  use calculator_assert, only : PORT_IN1, PORT_IN2, PORT_OUT
  use graph_carrier    , only : counted_set, subset_set, member_set
  use graph_relation   , only : stored_relation, relation
  use graph_relation_algebra, only : restrict_slot, project_slots, &
       &                             compose_binary
  use class_graph_field, only : field
  use arithmetic_constitution_fixture, only : apply_law

  implicit none

  integer, parameter :: ROW_C = 1
  integer, parameter :: ROW_E = 2

  type(counted_set)     :: x, o, p, y
  type(subset_set)      :: k, u, p_out
  type(stored_relation) :: flow, located
  type(field)           :: qk, qu
  integer               :: table(3, 6)
  integer               :: nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "calculator tower . level 8 . constitution"
  write(*,'(1x,a)') "============================================="

  x = counted_set('value-slots'  , 5)
  o = counted_set('operations'   , 2)
  p = counted_set('ports'        , 3)
  y = counted_set('residual-rows', 2)

  table(:, 1) = [OP_PLUS , SLOT_A, PORT_IN1]
  table(:, 2) = [OP_PLUS , SLOT_B, PORT_IN2]
  table(:, 3) = [OP_PLUS , SLOT_C, PORT_OUT]
  table(:, 4) = [OP_TIMES, SLOT_C, PORT_IN1]
  table(:, 5) = [OP_TIMES, SLOT_D, PORT_IN2]
  table(:, 6) = [OP_TIMES, SLOT_E, PORT_OUT]
  flow = stored_relation('flow', [o, x, p], table)

  located = stored_relation('located', [y, x], &
       & reshape([ROW_C, SLOT_C,  ROW_E, SLOT_E], [2, 2]))

  p_out = subset_set('output-port', p, [PORT_OUT])

  k = subset_set('known'  , x, [SLOT_D, SLOT_A, SLOT_B])
  u = subset_set('unknown', x, [SLOT_E, SLOT_C])
  qk = field('q known', k)
  call qk % set_real_vector([4.0_dp, 2.0_dp, 3.0_dp])
  qu = field('q unknown', u)

  call check_coverage(nfail)
  call check_laws(nfail)
  call check_generated_map(nfail)
  call check_solution(nfail)
  call check_topology_preserved(flow, nfail)
  call check_order_invariance(nfail)

  call verdict(nfail, "level 8")

contains

  !===================================================================!
  ! K and U are disjoint and together cover X.
  !===================================================================!

  subroutine check_coverage(nfail)

    integer, intent(inout) :: nfail

    integer :: i, m
    logical :: ok

    ok = .true.
    do i = 1, x % size()
       m = x % member(i)
       ok = ok .and. (count([k % has(m), u % has(m)]) .eq. 1)
    end do
    call report(ok, &
         & "K and U are disjoint and cover X together", nfail)

  end subroutine check_coverage

  !===================================================================!
  ! The new meaning itself - forbidden below this rung.
  !===================================================================!

  subroutine check_laws(nfail)

    integer, intent(inout) :: nfail

    call report(abs(apply_law(OP_PLUS, 2.0_dp, 3.0_dp) - 5.0_dp) &
         & < 1.0d-14, "law(+, 2, 3) = 5", nfail)
    call report(abs(apply_law(OP_TIMES, 5.0_dp, 4.0_dp) - 20.0_dp) &
         & < 1.0d-14, "law(x, 5, 4) = 20", nfail)

  end subroutine check_laws

  !===================================================================!
  ! One value, one home: known slots answer the known field,
  ! unknown slots the unknown state, outsiders are refused.
  !===================================================================!

  real(dp) function q_of(slot, ustate) result(v)

    integer , intent(in) :: slot
    real(dp), intent(in) :: ustate(:)

    real(dp), allocatable :: kv(:)

    if (k % has(slot)) then
       call qk % get_real_vector(kv)
       v = kv(k % local_index(slot))
    else if (u % has(slot)) then
       v = ustate(u % local_index(slot))
    else
       error stop 'constitution: a slot with no home holds no value'
    end if

  end function q_of

  !===================================================================!
  ! The one tuple the flow holds for (operation, port) - refused
  ! from the fixture if zero or many. The calculator relation
  ! already keeps this law; nothing in production needs to.
  !===================================================================!

  integer function the_slot(r_flow, op, port) result(slot)

    type(stored_relation), intent(in) :: r_flow
    integer, intent(in) :: op, port

    integer :: i, n

    n = 0
    do i = 1, x % size()
       if (r_flow % has([op, x % member(i), port])) then
          n = n + 1
          slot = x % member(i)
       end if
    end do
    if (n /= 1) error stop 'constitution: one port, one slot - or refusal'

  end function the_slot

  !===================================================================!
  ! The generated residual: structure chooses every slot, the law
  ! supplies every number. No row and no operand is named here.
  !===================================================================!

  subroutine generated_residual(r_flow, ustate, r)

    type(stored_relation), intent(in)  :: r_flow
    real(dp), intent(in)               :: ustate(:)
    real(dp), allocatable, intent(out) :: r(:)

    integer  :: i, j, row, x_out, op, in1, in2, nop
    real(dp) :: produced

    allocate(r(y % size()))
    do i = 1, y % size()
       row = y % member(i)

       ! L locates the row's computed slot.
       x_out = 0
       do j = 1, x % size()
          if (located % has([row, x % member(j)])) x_out = x % member(j)
       end do
       if (x_out == 0) error stop 'constitution: an unlocated residual row'

       ! The unique operation whose out port names that slot.
       op  = 0
       nop = 0
       do j = 1, o % size()
          if (r_flow % has([o % member(j), x_out, PORT_OUT])) then
             nop = nop + 1
             op  = o % member(j)
          end if
       end do
       if (nop /= 1) error stop 'constitution: one producer per slot - or refusal'

       in1 = the_slot(r_flow, op, PORT_IN1)
       in2 = the_slot(r_flow, op, PORT_IN2)

       produced = apply_law(op, q_of(in1, ustate), q_of(in2, ustate))
       r(y % local_index(row)) = q_of(x_out, ustate) - produced
    end do

  end subroutine generated_residual

  !===================================================================!
  ! The generated map IS the Level-7 system - probed at several
  ! off-solution states, the expected formulas living only here in
  ! the assertions.
  !===================================================================!

  subroutine check_generated_map(nfail)

    integer, intent(inout) :: nfail

    real(dp), allocatable :: r(:)

    ! U order is {e, c}.
    call generated_residual(flow, [0.0_dp, 0.0_dp], r)
    call report(abs(r(y % local_index(ROW_C)) + 5.0_dp) < 1.0d-14 .and. &
         &      abs(r(y % local_index(ROW_E)) - 0.0_dp) < 1.0d-14, &
         & "at q(e)=0, q(c)=0: r = (-5, 0)", nfail)

    call generated_residual(flow, [0.0_dp, 1.0_dp], r)
    call report(abs(r(y % local_index(ROW_C)) + 4.0_dp) < 1.0d-14 .and. &
         &      abs(r(y % local_index(ROW_E)) + 4.0_dp) < 1.0d-14, &
         & "at q(e)=0, q(c)=1: r = (-4, -4)", nfail)

    call generated_residual(flow, [1.0_dp, 0.0_dp], r)
    call report(abs(r(y % local_index(ROW_C)) + 5.0_dp) < 1.0d-14 .and. &
         &      abs(r(y % local_index(ROW_E)) - 1.0_dp) < 1.0d-14, &
         & "at q(e)=1, q(c)=0: r = (-5, 1)", nfail)

  end subroutine check_generated_map

  !===================================================================!
  ! At the Level-7 solution the generated residual vanishes: with K
  ! this is exactly q = (2, 3, 5, 4, 20), and no full field is ever
  ! manufactured to say so.
  !===================================================================!

  subroutine check_solution(nfail)

    integer, intent(inout) :: nfail

    real(dp), allocatable :: r(:)

    call generated_residual(flow, [20.0_dp, 5.0_dp], r)
    call report(abs(r(1)) < 1.0d-14 .and. abs(r(2)) < 1.0d-14, &
         & "at q(c)=5, q(e)=20 the constitution's residual vanishes", &
         & nfail)

  end subroutine check_solution

  !===================================================================!
  ! Meaning added, topology unchanged: the slots the generated
  ! equation touches equal the Level-6 J supports, row for row -
  ! dependency read structurally, never by perturbation.
  !===================================================================!

  subroutine check_topology_preserved(r_flow, nfail)

    type(stored_relation), intent(in) :: r_flow
    integer, intent(inout) :: nfail

    type(stored_relation)        :: q_, a_
    class(relation), allocatable :: b_, jac
    integer :: i, row, x_out, op, j
    logical :: ok

    ! The Level-6 road, walked again: J = A o (Q o L).
    q_  = project_slots(restrict_slot(r_flow, 3, p_out), [2, 1])
    b_  = compose_binary(located, q_)
    a_  = project_slots(r_flow, [1, 2])
    jac = compose_binary(b_, a_)

    ok = .true.
    do i = 1, y % size()
       row = y % member(i)
       ! the constitution's slots for this row
       x_out = 0
       do j = 1, x % size()
          if (located % has([row, x % member(j)])) x_out = x % member(j)
       end do
       op = merge(OP_PLUS, OP_TIMES, &
            &     r_flow % has([OP_PLUS, x_out, PORT_OUT]))
       do j = 1, x % size()
          ! every slot the generated equation uses is in J, and
          ! every slot J names is used: equality both ways.
          ok = ok .and. ( jac % has([row, x % member(j)]) .eqv. &
               & ( x % member(j) == x_out .or. &
               &   x % member(j) == the_slot(r_flow, op, PORT_IN1) .or. &
               &   x % member(j) == the_slot(r_flow, op, PORT_IN2) ) )
       end do
    end do
    call report(ok, &
         & "the constitution's dependency equals the level-6 supports", &
         & nfail)

  end subroutine check_topology_preserved

  !===================================================================!
  ! The constitution consumes a relation SET: the same six facts
  ! handed backwards generate the same residuals and supports.
  !===================================================================!

  subroutine check_order_invariance(nfail)

    integer, intent(inout) :: nfail

    type(stored_relation) :: backwards
    integer               :: rev(3, 6), j
    real(dp), allocatable :: r1(:), r2(:)

    do j = 1, 6
       rev(:, j) = table(:, 7 - j)
    end do
    backwards = stored_relation('flow backwards', [o, x, p], rev)

    call generated_residual(flow     , [7.0_dp, -2.0_dp], r1)
    call generated_residual(backwards, [7.0_dp, -2.0_dp], r2)
    call report(all(abs(r1 - r2) < 1.0d-14), &
         & "the tuples' order means nothing to the constitution", nfail)

  end subroutine check_order_invariance

end program calculator_level_8
