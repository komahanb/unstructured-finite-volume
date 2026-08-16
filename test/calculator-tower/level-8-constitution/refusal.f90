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
  use fractal_graph        , only : set_graph => graph
  use graph_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use graph_set_map        , only : set_map
  use graph_inclusion_map  , only : inclusion_map, declared_subobject
  use graph_relation   , only : stored_relation
  use arithmetic_constitution_fixture, only : generated_residual

  implicit none

  type(set_graph)     :: x, o, y
  type(set_graph)     :: ports
  type(set_graph)      :: k, u
  type(stored_relation) :: flow, doubled
  real(dp), allocatable :: r(:)
  integer               :: table(3, 6)
  type(set_map)     :: sets
  type(inclusion_map)     :: inclusions

  call x % declare()
  call sets % bind(x, counted_set_representation(5))
  call o % declare()
  call sets % bind(o, counted_set_representation(2))
  call y % declare()
  call sets % bind(y, counted_set_representation(2))

  table(:, 1) = [OP_PLUS , SLOT_A, PORT_IN1]
  table(:, 2) = [OP_PLUS , SLOT_B, PORT_IN2]
  table(:, 3) = [OP_PLUS , SLOT_C, PORT_OUT]
  table(:, 4) = [OP_TIMES, SLOT_C, PORT_IN1]
  table(:, 5) = [OP_TIMES, SLOT_D, PORT_IN2]
  table(:, 6) = [OP_TIMES, SLOT_E, PORT_OUT]
  call ports % declare()
  call sets % bind(ports, counted_set_representation(3))
  flow = stored_relation('flow', [o, x, ports], table, sets)

  ! Row 1 located at TWO slots: the schema law must refuse.
  doubled = stored_relation('doubled', [y, x], &
       & reshape([1, SLOT_C,  1, SLOT_D,  2, SLOT_E], [2, 3]), sets)

  call k % declare()
  call sets       % bind(k, listed_set_representation([SLOT_D, SLOT_A, SLOT_B]))
  call inclusions % include_in(k, x)
  call u % declare()
  call sets       % bind(u, listed_set_representation([SLOT_E, SLOT_C]))
  call inclusions % include_in(u, x)

  call generated_residual(flow, doubled, x, o, y, sets, &
       & k, [4.0_dp, 2.0_dp, 3.0_dp], u, [0.0_dp, 0.0_dp], r)

  write(*,'(1x,a)') "REACHED PAST THE REFUSAL"

end program constitution_refusal
