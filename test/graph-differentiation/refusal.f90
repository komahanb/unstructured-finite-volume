!=====================================================================!
! Invalid-input cases for the differentiation stack, one per
! invocation. The case is selected by the first command-line
! argument; each must terminate in error stop with the message
! run.sh expects, and a case that returns normally is reported as
! a failure by run.sh.
!
!      bdforder      set_bdf with order 3 (only 1 and 2 have
!                    tables)
!      bdfcount      a step count that disagrees with the order
!      bdfstep       a nonpositive time step
!      dupslot       two paths naming one input slot
!      badslot       a path naming an input slot that does not
!                    exist
!      negdegree     a negative derivative degree
!      pastcalculus  a degree needing more derivative slots than
!                    the operation's max_degree
!      hugedegree    a multinomial coefficient past the int64
!                    range (degree 21)
!      statepath     a caller-supplied path on the state slot,
!                    which march_directional computes itself
!      orderzero     march_directional with order 0
!      blindreverse  an implicit-rule march_adjoint without the
!                    action and trajectory it needs
!      unfrozen      an exact tangent applied before freeze
!      eagerclock    march_adaptive with duration 0
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program refusal

  use iso_fortran_env     , only : dp => REAL64
  use graph_field_calculus, only : graph_field
  use graph_fractal       , only : set_graph => graph
  use class_graph         , only : directed_stored_graph
  use class_graph_field   , only : field
  use class_graph_step    , only : step_operator
  use class_graph_chain_rule, only : chain_rule, argument_path
  use class_graph_linearization, only : exact_linearization
  use class_graph_marcher , only : marcher, MARCH_BACKWARD
  use class_graph_step_policy, only : halving_policy
  use toy_differentiable_forms, only : quartic_form, equilibrium_law, &
       & scalar_pair, fill_path

  implicit none

  type(quartic_form)    :: quartic
  type(equilibrium_law) :: equil

  type(directed_stored_graph) :: lone
  type(set_graph)             :: cells
  type(step_operator)         :: statement
  type(chain_rule)            :: composer
  type(argument_path)         :: paths(2)
  type(exact_linearization)   :: tangent
  type(marcher)               :: clock
  type(halving_policy)        :: policy

  type(field) :: inputs(2), direction
  class(graph_field), allocatable :: output
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
     call fill_path(paths(1), 1, [1.0_dp], cells)
     call fill_path(paths(2), 1, [1.0_dp], cells)
     call composer % assemble(quartic, lone, inputs, 1, paths, output)

  case ('badslot')

     call scalar_pair(1.0_dp, 2.0_dp, cells, inputs)
     call fill_path(paths(1), 1, [1.0_dp], cells)
     call fill_path(paths(2), 3, [1.0_dp], cells)
     call composer % assemble(quartic, lone, inputs, 1, paths, output)

  case ('negdegree')

     call scalar_pair(1.0_dp, 2.0_dp, cells, inputs)
     call fill_path(paths(1), 1, [1.0_dp], cells)
     call composer % assemble(quartic, lone, inputs, -1, paths(1:1), &
          & output)

  case ('pastcalculus')

     call scalar_pair(1.0_dp, 2.0_dp, cells, inputs)
     call fill_path(paths(1), 1, &
          & [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], cells)
     call composer % assemble(quartic, lone, inputs, 5, paths(1:1), &
          & output)

  case ('hugedegree')

     call scalar_pair(1.0_dp, 2.0_dp, cells, inputs)
     call fill_path(paths(1), 1, [1.0_dp], cells)
     call composer % assemble(quartic, lone, inputs, 21, paths(1:1), &
          & output)

  case ('statepath')

     call fill_path(paths(1), 1, [1.0_dp], cells)
     call clock % march_directional(equil, lone, 2, trajectory, 1, &
          & sensitivities, paths=paths(1:1))

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
     call tangent % apply(lone, [direction], output)

  case ('eagerclock')

     q = [2.0_dp]
     call clock % march_adaptive(equil, lone, q, 0.0_dp, policy, 3, &
          & taken, completed)

  case default

     error stop 'refusal: unknown case'

  end select

  write(*,*) 'refusal case survived: ', trim(which)

end program refusal
