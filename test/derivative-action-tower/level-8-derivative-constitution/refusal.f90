!=====================================================================!
! DERIVATIVE ACTION TOWER . LEVEL 8 . REFUSALS
!
! Each case must die, and die for its stated reason:
!
!   unbound-law ......... a symbol no law binds cannot be evaluated
!                         or linearized
!   missing-port ........ an operation without exactly one slot on a
!                         required port cannot be differentiated
!   primal-starvation ... a law cannot run before its primal inputs
!                         exist - computed values are produced by
!                         laws or not at all
!   tangent-starvation .. a tangent cannot be read before it is
!                         established - derivative state is earned,
!                         never fabricated
!
! The starvation cases hand the evaluators a deliberately WRONG
! order - [sum, product] - which no derivation would ever answer.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program derivative_level_8_refusal

  use iso_fortran_env  , only : dp => REAL64
  use derivative_assert, only : SLOT_X, SLOT_Y, SLOT_U, SLOT_Z
  use derivative_assert, only : OP_PRODUCT, OP_SUM
  use derivative_assert, only : PORT_IN1, PORT_IN2, PORT_OUT
  use fractal_graph        , only : set_graph => graph
  use graph_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use graph_set_map        , only : set_map
  use graph_inclusion_map  , only : inclusion_map, declared_subobject
  use graph_relation   , only : stored_relation
  use derivative_constitution_fixture, only : apply_law, &
       & slot_for_port, primal_execution, tangent_action

  implicit none

  type(set_graph)     :: v, o, p
  type(set_graph)      :: x_dom, c
  type(stored_relation) :: flow, lame
  integer               :: table(3, 6), short(3, 5)
  real(dp)              :: vals(4), dots(4), zz
  logical               :: got(4), dgot(4)
  integer               :: found
  character(len=32)     :: which
  type(set_map)     :: sets
  type(inclusion_map)     :: inclusions

  if (command_argument_count() .lt. 1) then
     error stop 'usage: refusal <case>'
  end if
  call get_command_argument(1, which)

  call v % declare()
  call sets % bind(v, counted_set_representation(4))
  call o % declare()
  call sets % bind(o, counted_set_representation(2))
  call p % declare()
  call sets % bind(p, counted_set_representation(3))
  call x_dom % declare()
  call sets       % bind(x_dom, listed_set_representation([SLOT_Y, SLOT_X]))
  call inclusions % include_in(x_dom, v)
  call c % declare()
  call sets       % bind(c, listed_set_representation([SLOT_U, SLOT_Z]))
  call inclusions % include_in(c, v)

  table(:, 1) = [OP_PRODUCT, SLOT_X, PORT_IN1]
  table(:, 2) = [OP_PRODUCT, SLOT_Y, PORT_IN2]
  table(:, 3) = [OP_PRODUCT, SLOT_U, PORT_OUT]
  table(:, 4) = [OP_SUM    , SLOT_U, PORT_IN1]
  table(:, 5) = [OP_SUM    , SLOT_Y, PORT_IN2]
  table(:, 6) = [OP_SUM    , SLOT_Z, PORT_OUT]

  select case (trim(which))

  case ('unbound-law')
     zz = apply_law(9, 1.0_dp, 1.0_dp)
     write(*,*) 'a lawless symbol computed', zz

  case ('missing-port')
     ! Five facts only: sum's in2 is never declared.
     short(:, 1:3) = table(:, 1:3)
     short(:, 4)   = table(:, 4)
     short(:, 5)   = table(:, 6)
     lame = stored_relation('flow missing a port', [o, v, p], short, sets)
     found = slot_for_port(lame, v, sets, OP_SUM, PORT_IN2)
     write(*,*) 'a missing port named slot', found

  case ('primal-starvation')
     flow = stored_relation('flow', [o, v, p], table, sets)
     ! The wrong order, on purpose: sum before product.
     call primal_execution(flow, v, sets, x_dom, [3.0_dp, 2.0_dp], c, &
          & [OP_SUM, OP_PRODUCT], vals, got)
     write(*,*) 'a starved operation computed', vals

  case ('tangent-starvation')
     flow = stored_relation('flow', [o, v, p], table, sets)
     call primal_execution(flow, v, sets, x_dom, [3.0_dp, 2.0_dp], c, &
          & [OP_PRODUCT, OP_SUM], vals, got)
     ! Primal is honest; the tangent order is not.
     call tangent_action(flow, v, sets, x_dom, [1.0_dp, 0.0_dp], c, &
          & [OP_SUM, OP_PRODUCT], vals, dots, dgot)
     write(*,*) 'a starved tangent computed', dots

  case default
     error stop 'usage: refusal <case>'

  end select

end program derivative_level_8_refusal
