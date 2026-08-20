!=====================================================================!
! CALCULATOR TOWER . LEVEL 8 . CONSTITUTION
!
! The level answers one question: WHAT LAWS DO THE SYMBOLS OBEY.
! Only here, in one small table, + becomes addition and x becomes
! multiplication - the table is TEST-LOCAL, because one calculator
! has not earned a universal production constitution. Everything
! else was already true below:
!
!      T_flow chooses each operation's in1, in2 and out slots
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

program calculator_level_8

  use iso_fortran_env  , only : dp => REAL64
  use calculator_assert, only : report, verdict
  use calculator_assert, only : SLOT_A, SLOT_B, SLOT_C, SLOT_D, SLOT_E
  use calculator_assert, only : OP_PLUS, OP_TIMES
  use calculator_assert, only : PORT_IN1, PORT_IN2, PORT_OUT
  use graph_fractal        , only : set_graph => graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set        , only : set_map
  use map_inclusion  , only : inclusion_map, declared_subobject
  use relation_finitary   , only : stored_relation, relation
  use relation_algebra, only : restrict_slot, project_slots, &
       &                             compose_binary
  use field_stored, only : field
  use arithmetic_constitution_fixture, only : apply_law, &
       &  generated_residual, constitution_support

  implicit none

  integer, parameter :: ROW_C = 1
  integer, parameter :: ROW_E = 2

  type(set_graph)     :: x, o, p, y
  type(set_graph)      :: k, u, p_out
  type(stored_relation) :: flow, located
  type(field)           :: qk, qu
  integer               :: table(3, 6)
  integer               :: nfail
  type(set_map)     :: sets
  type(inclusion_map)     :: inclusions

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "calculator tower . level 8 . constitution"
  write(*,'(1x,a)') "============================================="

  call x % declare()
  call sets % bind(x, counted_set_representation(5))
  call o % declare()
  call sets % bind(o, counted_set_representation(2))
  call p % declare()
  call sets % bind(p, counted_set_representation(3))
  call y % declare()
  call sets % bind(y, counted_set_representation(2))

  table(:, 1) = [OP_PLUS , SLOT_A, PORT_IN1]
  table(:, 2) = [OP_PLUS , SLOT_B, PORT_IN2]
  table(:, 3) = [OP_PLUS , SLOT_C, PORT_OUT]
  table(:, 4) = [OP_TIMES, SLOT_C, PORT_IN1]
  table(:, 5) = [OP_TIMES, SLOT_D, PORT_IN2]
  table(:, 6) = [OP_TIMES, SLOT_E, PORT_OUT]
  flow = stored_relation('flow', [o, x, p], table, sets)

  located = stored_relation('located', [y, x], &
       & reshape([ROW_C, SLOT_C,  ROW_E, SLOT_E], [2, 2]), sets)

  call p_out % declare()
  call sets       % bind(p_out, listed_set_representation([PORT_OUT]))
  call inclusions % include_in(p_out, p)

  call k % declare()
  call sets       % bind(k, listed_set_representation([SLOT_D, SLOT_A, SLOT_B]))
  call inclusions % include_in(k, x)
  call u % declare()
  call sets       % bind(u, listed_set_representation([SLOT_E, SLOT_C]))
  call inclusions % include_in(u, x)
  qk = field('q known', k, sets % size_of(k))
  call qk % set_real_vector([4.0_dp, 2.0_dp, 3.0_dp])
  qu = field('q unknown', u, sets % size_of(u))

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
    do i = 1, sets % size_of(x)
       m = sets % member_of(x, i)
       ok = ok .and. (count([sets % has_in(k, m), sets % has_in(u, m)]) .eq. 1)
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
  ! The generated map IS the Level-7 system - probed at several
  ! off-solution states, the expected formulas living only here in
  ! the assertions.
  !===================================================================!

  subroutine check_generated_map(nfail)

    integer, intent(inout) :: nfail

    real(dp), allocatable :: r(:)

    ! U order is {e, c}.
    call gen(flow, [0.0_dp, 0.0_dp], r)
    call report(abs(r(sets % index_in(y, ROW_C)) + 5.0_dp) < 1.0d-14 .and. &
         &      abs(r(sets % index_in(y, ROW_E)) - 0.0_dp) < 1.0d-14, &
         & "at q(e)=0, q(c)=0: r = (-5, 0)", nfail)

    call gen(flow, [0.0_dp, 1.0_dp], r)
    call report(abs(r(sets % index_in(y, ROW_C)) + 4.0_dp) < 1.0d-14 .and. &
         &      abs(r(sets % index_in(y, ROW_E)) + 4.0_dp) < 1.0d-14, &
         & "at q(e)=0, q(c)=1: r = (-4, -4)", nfail)

    call gen(flow, [1.0_dp, 0.0_dp], r)
    call report(abs(r(sets % index_in(y, ROW_C)) + 5.0_dp) < 1.0d-14 .and. &
         &      abs(r(sets % index_in(y, ROW_E)) - 1.0_dp) < 1.0d-14, &
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

    call gen(flow, [20.0_dp, 5.0_dp], r)
    call report(abs(r(1)) < 1.0d-14 .and. abs(r(2)) < 1.0d-14, &
         & "at q(c)=5, q(e)=20 the constitution's residual vanishes", &
         & nfail)

  end subroutine check_solution

  !===================================================================!
  ! Meaning added, topology unchanged: the slots the generated
  ! equation touches equal the Level-6 J supports, row for row -
  ! dependency read structurally, never by perturbation.
  !===================================================================!

  subroutine check_topology_preserved(t_flow, nfail)

    type(stored_relation), intent(in) :: t_flow
    integer, intent(inout) :: nfail

    type(stored_relation)        :: q_, a_
    class(relation), allocatable :: b_, jac
    integer, allocatable :: mine(:)
    integer :: i, row, j
    logical :: ok

    ! The Level-6 road, walked again: J = A o (Q o L).
    q_  = project_slots(restrict_slot(t_flow, 3, p_out, sets, inclusions), [2, 1], sets)
    b_  = compose_binary(located, q_, sets)
    a_  = project_slots(t_flow, [1, 2], sets)
    jac = compose_binary(b_, a_, sets)

    ok = .true.
    do i = 1, sets % size_of(y)
       row = sets % member_of(y, i)
       call constitution_support(t_flow, located, x, o, sets, row, mine)
       do j = 1, sets % size_of(x)
          ok = ok .and. ( jac % has([row, sets % member_of(x, j)]) .eqv. &
               &          any(mine == sets % member_of(x, j)) )
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
    backwards = stored_relation('flow backwards', [o, x, p], rev, sets)

    call gen(flow     , [7.0_dp, -2.0_dp], r1)
    call gen(backwards, [7.0_dp, -2.0_dp], r2)
    call report(all(abs(r1 - r2) < 1.0d-14), &
         & "the tuples' order means nothing to the constitution", nfail)

    block
      integer, allocatable :: s1(:), s2(:)
      integer :: i, jj
      logical :: ok
      ok = .true.
      do i = 1, sets % size_of(y)
         call constitution_support(flow     , located, x, o, sets, sets % member_of(y, i), s1)
         call constitution_support(backwards, located, x, o, sets, sets % member_of(y, i), s2)
         ok = ok .and. (size(s1) .eq. size(s2))
         do jj = 1, size(s1)
            ok = ok .and. any(s2 == s1(jj))
         end do
      end do
      call report(ok, &
           & "and neither does it to the generated supports", nfail)
    end block

  end subroutine check_order_invariance

  !===================================================================!
  ! The fixture, applied with this calculator's own objects.
  !===================================================================!

  subroutine gen(t_flow, ustate, r)

    type(stored_relation), intent(in)  :: t_flow
    real(dp), intent(in)               :: ustate(:)
    real(dp), allocatable, intent(out) :: r(:)

    real(dp), allocatable :: kv(:)

    call qk % get_real_vector(kv)
    call generated_residual(t_flow, located, x, o, y, sets, k, kv, u, ustate, r)

  end subroutine gen

end program calculator_level_8
