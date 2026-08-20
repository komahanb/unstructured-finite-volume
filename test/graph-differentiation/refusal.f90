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
!      dupslot       two paths naming one argument
!      badslot       an argument the operation does not declare
!      foreignpath   a path naming another operation's argument
!      foreignvariation  a variation on another operation's
!                    argument, of the same position
!      undeclared    an argument asked of an operation built
!                    without its constructor
!      historyreach  history(2) of a reach-1 statement
!      historyshape  a supplied history state of another storage
!                    shape than the state
!      negdegree     a negative derivative degree
!      pastcalculus  a degree needing more derivative slots than
!                    the operation's max_degree
!      hugedegree    a multinomial coefficient past the int64
!                    range (degree 21)
!      statepath     a caller-supplied path on the state slot,
!                    which march_directional computes itself
!      orderzero     march_directional with order 0
!      unfrozen      an exact tangent applied before freeze
!      eagerclock    march_adaptive with duration 0
!      flatcalculus  a partial action requested from an operation
!                    that declares none
!      shallowcalculus  march_directional past the action's
!                    max_degree
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program refusal

  use iso_fortran_env     , only : dp => REAL64
  use field_calculus, only : field
  use graph_fractal       , only : graph
  use view_directed_stored         , only : stored_directed_graph
  use field_stored   , only : stored_field
  use operation_action  , only : variation
  use operation_step    , only : scheme, backward_euler
  use operation_chain_rule, only : chain_rule, argument_path
  use operation_linearization, only : linearization, tangent_of
  use operation_marching , only : marcher, MARCH_BACKWARD
  use operation_step_policy, only : halving_policy
  use toy_differentiable_forms, only : quartic_form, equilibrium_law, &
       & linear_law, scalar_pair, fill_path

  implicit none

  type(quartic_form)    :: quartic
  type(equilibrium_law) :: equil
  type(linear_law)      :: lin
  type(linear_law)      :: bare

  type(stored_directed_graph) :: lone
  type(graph)             :: cells
  type(scheme)         :: statement
  type(chain_rule)            :: composer
  type(argument_path)         :: paths(2)
  type(linearization)         :: tangent
  type(marcher)               :: clock
  type(halving_policy)        :: policy

  type(stored_field) :: inputs(2), direction, wide
  class(field), allocatable :: output
  real(dp), allocatable :: sensitivities(:,:,:), taken(:)
  real(dp) :: trajectory(1,3), q(1)
  logical  :: completed

  character(len=32) :: which

  call get_command_argument(1, which)

  lone  = stored_directed_graph(1, tails=[integer ::], heads=[integer ::])
  cells = lone % vertex_set()

  quartic = quartic_form()
  equil   = equilibrium_law()
  lin     = linear_law()

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
     call fill_path(paths(1), quartic % argument(1), [1.0_dp], cells)
     call fill_path(paths(2), quartic % argument(1), [1.0_dp], cells)
     call composer % assemble(quartic, lone, inputs, 1, paths, output)

  case ('badslot')

     call fill_path(paths(1), quartic % argument(3), [1.0_dp], cells)

  case ('foreignpath')

     ! equil's first argument is not quartic's, though both are position 1
     call scalar_pair(1.0_dp, 2.0_dp, cells, inputs)
     call fill_path(paths(1), equil % argument(1), [1.0_dp], cells)
     call composer % assemble(quartic, lone, inputs, 1, paths(1:1), output)

  case ('foreignvariation')

     call scalar_pair(1.0_dp, 2.0_dp, cells, inputs)
     direction = stored_field('v', cells, 1, num_components=1)
     call direction % set_real_vector([1.0_dp])
     call quartic % partial_action(lone, inputs, &
          & [variation(equil % argument(1), direction)], output)

  case ('undeclared')

     ! bare was never built by linear_law(), so it owns no arguments
     call fill_path(paths(1), bare % argument(1), [1.0_dp], cells)

  case ('historyshape')

     statement = backward_euler(quartic, 0.5_dp)
     call scalar_pair(1.0_dp, 2.0_dp, cells, inputs)
     wide = stored_field('qold', cells, 1, num_components=2)
     call wide % set_real_vector([1.0_dp, 1.0_dp])
     call statement % apply(lone, [inputs(1), inputs(2), wide], output)

  case ('historyreach')

     statement = backward_euler(quartic, 1.0_dp)
     call fill_path(paths(1), statement % history(2), [1.0_dp], cells)

  case ('negdegree')

     call scalar_pair(1.0_dp, 2.0_dp, cells, inputs)
     call fill_path(paths(1), quartic % argument(1), [1.0_dp], cells)
     call composer % assemble(quartic, lone, inputs, -1, paths(1:1), &
          & output)

  case ('pastcalculus')

     call scalar_pair(1.0_dp, 2.0_dp, cells, inputs)
     call fill_path(paths(1), quartic % argument(1), &
          & [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], cells)
     call composer % assemble(quartic, lone, inputs, 5, paths(1:1), &
          & output)

  case ('hugedegree')

     call scalar_pair(1.0_dp, 2.0_dp, cells, inputs)
     call fill_path(paths(1), quartic % argument(1), [1.0_dp], cells)
     call composer % assemble(quartic, lone, inputs, 21, paths(1:1), &
          & output)

  case ('statepath')

     call fill_path(paths(1), equil % argument(1), [1.0_dp], cells)
     call clock % march_directional(equil, lone, 2, trajectory, 1, &
          & sensitivities, paths=paths(1:1))

  case ('orderzero')

     call clock % march_directional(equil, lone, 2, trajectory, 0, &
          & sensitivities)

  case ('unfrozen')

     tangent = tangent_of(quartic)
     direction = stored_field('v', cells, 1, num_components=1)
     call direction % set_real_vector([1.0_dp])
     call tangent % apply(lone, [direction], output)

  case ('eagerclock')

     q = [2.0_dp]
     call clock % march_adaptive(equil, lone, q, 0.0_dp, policy, 3, &
          & taken, completed)

  case ('flatcalculus')

     direction = stored_field('v', cells, 1, num_components=1)
     call direction % set_real_vector([1.0_dp])
     call lin % partial_action(lone, [direction], &
          & [variation(lin % argument(1), direction)], output)

  case ('shallowcalculus')

     call clock % march_directional(lin, lone, 2, trajectory, 1, &
          & sensitivities)

  case default

     error stop 'refusal: unknown case'

  end select

  write(*,*) 'refusal case survived: ', trim(which)

end program refusal
