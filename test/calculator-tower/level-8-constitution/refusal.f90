!=====================================================================!
! The level-8 refusal, EXPECTED TO DIE: a residual row located at
! TWO value slots violates the constitution's one-location schema
! law, and the fixture refuses it - a test-local law, no production
! functional-relation abstraction earned.
!=====================================================================!
program constitution_refusal

  use iso_fortran_env  , only : dp => REAL64
  use calculator_assert, only : SLOT_A, SLOT_B, SLOT_C, SLOT_D, SLOT_E
  use calculator_assert, only : OP_PLUS, OP_TIMES
  use calculator_assert, only : PORT_IN1, PORT_IN2, PORT_OUT
  use graph_carrier    , only : counted_set, subset_set
  use graph_relation   , only : stored_relation
  use arithmetic_constitution_fixture, only : generated_residual

  implicit none

  type(counted_set)     :: x, o, y
  type(subset_set)      :: k, u
  type(stored_relation) :: flow, doubled
  real(dp), allocatable :: r(:)
  integer               :: table(3, 6)

  x = counted_set('value-slots'  , 5)
  o = counted_set('operations'   , 2)
  y = counted_set('residual-rows', 2)

  table(:, 1) = [OP_PLUS , SLOT_A, PORT_IN1]
  table(:, 2) = [OP_PLUS , SLOT_B, PORT_IN2]
  table(:, 3) = [OP_PLUS , SLOT_C, PORT_OUT]
  table(:, 4) = [OP_TIMES, SLOT_C, PORT_IN1]
  table(:, 5) = [OP_TIMES, SLOT_D, PORT_IN2]
  table(:, 6) = [OP_TIMES, SLOT_E, PORT_OUT]
  flow = stored_relation('flow', [o, x, counted_set('ports', 3)], table)

  ! Row 1 located at TWO slots: the schema law must refuse.
  doubled = stored_relation('doubled', [y, x], &
       & reshape([1, SLOT_C,  1, SLOT_D,  2, SLOT_E], [2, 3]))

  k = subset_set('known'  , x, [SLOT_D, SLOT_A, SLOT_B])
  u = subset_set('unknown', x, [SLOT_E, SLOT_C])

  call generated_residual(flow, doubled, x, o, y, &
       & k, [4.0_dp, 2.0_dp, 3.0_dp], u, [0.0_dp, 0.0_dp], r)

  write(*,'(1x,a)') "REACHED PAST THE REFUSAL"

end program constitution_refusal
