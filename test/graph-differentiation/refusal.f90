!=====================================================================!
! The refusals that must die on the differentiation tower, one per
! invocation:
!
!      bdforder      an order the table shelf does not carry
!      bdfcount      a step count disagreeing with the order
!      bdfstep       a nonpositive time step
!      dupslot       two channels naming one input slot
!      badslot       a channel naming a slot the statement lacks
!      negdegree     a negative derivative degree
!      pastcalculus  a composition needing a partial past the
!                    statement's declared calculus
!      hugedegree    a multinomial count past the integer range
!      statepath     a caller-supplied path on the state slot
!      orderzero     a directional walk of order zero
!      blindreverse  a governed reverse walk with no trajectory
!      unfrozen      an exact tangent taken at no state
!      eagerclock    an adaptive walk over no duration
!
! Every case must error stop; a case that returns is a failure of
! the suite.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program refusal

  use iso_fortran_env     , only : dp => REAL64
  use graph_field_calculus, only : graph_field
  use fractal_graph       , only : set_graph => graph
  use class_graph         , only : directed_stored_graph
  use class_graph_field   , only : field
  use class_graph_step    , only : step_operator
  use class_graph_chain_rule, only : chain_rule, chain_channel
  use class_graph_exact_linearization, only : exact_linearization
  use class_graph_marcher , only : marcher, MARCH_BACKWARD
  use class_graph_step_policy, only : halving_policy
  use toy_differentiable_forms, only : quartic_form, equilibrium_law, &
       & scalar_pair, fill_scalar_channel

  implicit none

  type(quartic_form)    :: quartic
  type(equilibrium_law) :: equil

  type(directed_stored_graph) :: lone
  type(set_graph)             :: cells
  type(step_operator)         :: statement
  type(chain_rule)            :: composer
  type(chain_channel)         :: channels(2)
  type(exact_linearization)   :: tangent
  type(marcher)               :: clock
  type(halving_policy)        :: policy

  type(field) :: inputs(2), direction
  class(graph_field), allocatable :: answer
  real(dp), allocatable :: sensitivities(:,:,:), taken(:)
  real(dp) :: trajectory(1,3), lambda(1), q(1)
  logical  :: completed

  character(len=32) :: which

  call get_command_argument(1, which)

  lone  = directed_stored_graph(1, tails=[integer ::], heads=[integer ::])
  cells = lone % vertex_set()

  clock % rule = MARCH_BACKWARD
  trajectory   = 1.0_dp

  select case (trim(which))

  case ('bdforder')

     call statement % set_bdf(3, [1.0_dp, 1.0_dp, 1.0_dp])

  case ('bdfcount')

     call statement % set_bdf(2, [1.0_dp])

  case ('bdfstep')

     call statement % set_bdf(1, [-1.0_dp])

  case ('dupslot')

     call scalar_pair(1.0_dp, 2.0_dp, cells, inputs)
     call fill_scalar_channel(channels(1), 1, [1.0_dp], cells)
     call fill_scalar_channel(channels(2), 1, [1.0_dp], cells)
     call composer % assemble(quartic, lone, inputs, 1, channels, answer)

  case ('badslot')

     call scalar_pair(1.0_dp, 2.0_dp, cells, inputs)
     call fill_scalar_channel(channels(1), 1, [1.0_dp], cells)
     call fill_scalar_channel(channels(2), 3, [1.0_dp], cells)
     call composer % assemble(quartic, lone, inputs, 1, channels, answer)

  case ('negdegree')

     call scalar_pair(1.0_dp, 2.0_dp, cells, inputs)
     call fill_scalar_channel(channels(1), 1, [1.0_dp], cells)
     call composer % assemble(quartic, lone, inputs, -1, channels(1:1), &
          & answer)

  case ('pastcalculus')

     call scalar_pair(1.0_dp, 2.0_dp, cells, inputs)
     call fill_scalar_channel(channels(1), 1, &
          & [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], cells)
     call composer % assemble(quartic, lone, inputs, 5, channels(1:1), &
          & answer)

  case ('hugedegree')

     call scalar_pair(1.0_dp, 2.0_dp, cells, inputs)
     call fill_scalar_channel(channels(1), 1, [1.0_dp], cells)
     call composer % assemble(quartic, lone, inputs, 21, channels(1:1), &
          & answer)

  case ('statepath')

     call fill_scalar_channel(channels(1), 1, [1.0_dp], cells)
     call clock % march_directional(equil, lone, 2, trajectory, 1, &
          & sensitivities, channels=channels(1:1))

  case ('orderzero')

     call clock % march_directional(equil, lone, 2, trajectory, 0, &
          & sensitivities)

  case ('blindreverse')

     lambda = [1.0_dp]
     call clock % march_adjoint(equil, lone, lambda, 2)

  case ('unfrozen')

     tangent = exact_linearization(quartic)
     direction = field('v', cells, 1, ncomp=1)
     call direction % set_real_vector([1.0_dp])
     call tangent % apply(lone, [direction], answer)

  case ('eagerclock')

     q = [2.0_dp]
     call clock % march_adaptive(equil, lone, q, 0.0_dp, policy, 3, &
          & taken, completed)

  case default

     error stop 'refusal: unknown case'

  end select

  write(*,*) 'refusal case survived: ', trim(which)

end program refusal
