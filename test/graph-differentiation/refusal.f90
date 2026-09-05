!=====================================================================!
! Invalid-input cases for the differentiation stack, one per
! invocation. The case is selected by the first command-line
! argument; each must terminate in error stop with the message
! run.sh expects, and a case that returns normally is reported as
! a failure by run.sh.
!
!      dupslot       two paths naming one argument
!      badslot       an argument the operation does not declare
!      foreignpath   a path naming another operation's argument
!      foreignvariation  a variation on another operation's
!                    argument, of the same position
!      undeclared    an argument requested of an operation built
!                    without its constructor
!      negdegree     a negative derivative degree
!      pastcalculus  a degree needing more derivative factors than
!                    the operation's max_degree
!      hugedegree    a multinomial coefficient past the int64
!                    range (degree 21)
!      unfrozen      an exact tangent applied before freeze
!      flatcalculus  a partial action requested from an operation
!                    that declares none
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
  use operation_chain_rule, only : total_derivative, argument_path
  use operation_linearization, only : linearization, tangent_of
  use toy_differentiable_forms, only : quartic_form, equilibrium_law, &
       & linear_law, scalar_pair, fill_path

  implicit none

  type(quartic_form)    :: quartic
  type(equilibrium_law) :: equil
  type(linear_law)      :: lin
  type(linear_law)      :: bare

  type(stored_directed_graph) :: lone
  type(graph)             :: cells
  type(total_derivative)      :: composer
  type(argument_path)         :: paths(2)
  type(linearization)         :: tangent

  type(stored_field) :: inputs(2), direction
  class(field), allocatable :: output

  character(len=32) :: which

  call get_command_argument(1, which)

  lone  = stored_directed_graph(1, tails=[integer ::], heads=[integer ::])
  cells = lone % vertex_set()

  quartic = quartic_form()
  equil   = equilibrium_law()
  lin     = linear_law()

  select case (trim(which))

  case ('dupslot')

     composer = total_derivative(1)
     call scalar_pair(1.0_dp, 2.0_dp, cells, inputs)
     call fill_path(paths(1), quartic % argument(1), [1.0_dp], cells)
     call fill_path(paths(2), quartic % argument(1), [1.0_dp], cells)
     call composer % assemble(quartic, lone, quartic % bind(inputs), 1, paths, output)

  case ('badslot')

     call fill_path(paths(1), quartic % argument(3), [1.0_dp], cells)

  case ('foreignpath')

     ! equil's first argument is not quartic's, though both are position 1
     composer = total_derivative(1)
     call scalar_pair(1.0_dp, 2.0_dp, cells, inputs)
     call fill_path(paths(1), equil % argument(1), [1.0_dp], cells)
     call composer % assemble(quartic, lone, quartic % bind(inputs), 1, paths(1:1), output)

  case ('foreignvariation')

     call scalar_pair(1.0_dp, 2.0_dp, cells, inputs)
     direction = stored_field('v', cells, 1, num_components=1)
     call direction % set_real_vector([1.0_dp])
     call quartic % partial_action(lone, quartic % bind(inputs), &
          & [variation(equil % argument(1), direction)], output)

  case ('undeclared')

     ! bare was never built by linear_law(), so it owns no arguments
     call fill_path(paths(1), bare % argument(1), [1.0_dp], cells)

  case ('negdegree')

     composer = total_derivative(1)
     call scalar_pair(1.0_dp, 2.0_dp, cells, inputs)
     call fill_path(paths(1), quartic % argument(1), [1.0_dp], cells)
     call composer % assemble(quartic, lone, quartic % bind(inputs), -1, paths(1:1), &
          & output)

  case ('pastcalculus')

     composer = total_derivative(5)
     call scalar_pair(1.0_dp, 2.0_dp, cells, inputs)
     call fill_path(paths(1), quartic % argument(1), &
          & [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], cells)
     call composer % assemble(quartic, lone, quartic % bind(inputs), 5, paths(1:1), &
          & output)

  case ('hugedegree')

     composer = total_derivative(21)

  case ('unfrozen')

     tangent = tangent_of(quartic)
     direction = stored_field('v', cells, 1, num_components=1)
     call direction % set_real_vector([1.0_dp])
     call tangent % apply(lone, tangent % bind([direction]), output)

  case ('flatcalculus')

     direction = stored_field('v', cells, 1, num_components=1)
     call direction % set_real_vector([1.0_dp])
     call lin % partial_action(lone, lin % bind([direction]), &
          & [variation(lin % argument(1), direction)], output)

  case default

     error stop 'refusal: unknown case'

  end select

  write(*,*) 'refusal case survived: ', trim(which)

end program refusal
