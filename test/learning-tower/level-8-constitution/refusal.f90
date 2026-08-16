!=====================================================================!
! LEARNING TOWER . LEVEL 8 . REFUSALS
!
! Each case must die, and die for its stated reason:
!
!   unbound-law ..... a symbol no law binds cannot be evaluated
!   starved-input ... an operation scheduled before its inputs
!                     exist cannot fabricate them - computed values
!                     are produced by laws or not at all
!
! The starved case hands the evaluator a deliberately WRONG order -
! [error, predict] - which no derivation would ever answer: yhat is
! demanded before any law has produced it.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program learning_level_8_refusal

  use iso_fortran_env, only : dp => REAL64
  use learning_assert, only : SLOT_W, SLOT_X, SLOT_YHAT, SLOT_Y, SLOT_E
  use learning_assert, only : OP_PREDICT, OP_ERROR
  use learning_assert, only : PORT_IN1, PORT_IN2, PORT_OUT
  use fractal_graph        , only : set_graph => graph
  use graph_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use graph_set_map        , only : set_map
  use graph_inclusion_map  , only : inclusion_map, declared_subobject
  use graph_relation , only : stored_relation
  use learning_constitution_fixture, only : apply_law, generated_residual

  implicit none

  integer, parameter :: ROW_R = 1

  type(set_graph)     :: v, o, p, y
  type(set_graph)      :: k, theta, u
  type(stored_relation) :: flow, located
  integer               :: table(3, 6)
  real(dp)              :: r(1), z
  character(len=32)     :: which
  type(set_map)     :: sets
  type(inclusion_map)     :: inclusions

  if (command_argument_count() .lt. 1) then
     error stop 'usage: refusal <case>'
  end if
  call get_command_argument(1, which)

  select case (trim(which))

  case ('unbound-law')
     z = apply_law(9, 1.0_dp, 1.0_dp)
     write(*,*) 'a lawless symbol computed', z

  case ('starved-input')
     call v % declare()
     call sets % bind(v, counted_set_representation(5))
     call o % declare()
     call sets % bind(o, counted_set_representation(2))
     call p % declare()
     call sets % bind(p, counted_set_representation(3))
     call y % declare()
     call sets % bind(y, counted_set_representation(1))
     call k % declare()
     call sets       % bind(k, listed_set_representation([SLOT_Y, SLOT_X]))
     call inclusions % include_in(k, v)
     call theta % declare()
     call sets       % bind(theta, listed_set_representation([SLOT_W]))
     call inclusions % include_in(theta, v)
     call u % declare()
     call sets       % bind(u, listed_set_representation([SLOT_E, SLOT_YHAT]))
     call inclusions % include_in(u, v)

     table(:, 1) = [OP_PREDICT, SLOT_W   , PORT_IN1]
     table(:, 2) = [OP_PREDICT, SLOT_X   , PORT_IN2]
     table(:, 3) = [OP_PREDICT, SLOT_YHAT, PORT_OUT]
     table(:, 4) = [OP_ERROR  , SLOT_YHAT, PORT_IN1]
     table(:, 5) = [OP_ERROR  , SLOT_Y   , PORT_IN2]
     table(:, 6) = [OP_ERROR  , SLOT_E   , PORT_OUT]
     flow = stored_relation('flow', [o, v, p], table, sets)

     located = stored_relation('located', [y, v], &
          & reshape([ROW_R, SLOT_E], [2, 1]), sets)

     ! The wrong order, on purpose: error before predict.
     call generated_residual(flow, located, v, sets, y, &
          & k, [6.0_dp, 2.0_dp], theta, [0.0_dp], u, &
          & [OP_ERROR, OP_PREDICT], r)
     write(*,*) 'a starved operation computed', r

  case default
     error stop 'usage: refusal <case>'

  end select

end program learning_level_8_refusal
