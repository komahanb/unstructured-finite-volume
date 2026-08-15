!=====================================================================!
! ADJOINT TOWER . LEVEL 7 . THE SAME-SIZE REFUSALS
!
! Gate A's theorem, made numerical and made to bite. Every field
! offered below has exactly the RIGHT SIZE and the WRONG IDENTITY,
! so nothing but equals can reject it:
!
!   primal-rhs-on-Q ..... a right-hand side on Q handed to a solver
!                         whose residual domain is Y
!   adjoint-rhs-on-Y .... a right-hand side on Y handed to a solver
!                         whose residual domain is Q
!   primal-state-on-Y ... a state on Y handed to the primal equation
!   adjoint-covector-on-Q  a covector on Q handed to the adjoint one
!
! The first two are refused by PRODUCTION - graph_minimization
! already knows that a right-hand side lives on the residual
! domain - and the last two by the test-local equations themselves.
! Both layers must hold: a solver that trusted sizes, or an
! equation that read whatever it was handed, would return a
! plausible wrong answer rather than stopping.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program adjoint_level_7_refusal

  use iso_fortran_env  , only : dp => REAL64
  use adjoint_assert   , only : VAR_U, VAR_V, TGT_R1, TGT_R2
  use graph_set    , only : index_set, subset
  use graph_grammar    , only : graph_field
  use class_graph      , only : stored_graph
  use class_graph_field, only : field
  use class_graph_gmres, only : gmres
  use opaque_equation_fixture, only : opaque_primal, opaque_adjoint

  implicit none

  type(index_set)               :: v, t
  type(subset)                :: q_dom, y_dom
  type(stored_graph)              :: host
  type(opaque_primal)             :: primal_eq
  type(opaque_adjoint)            :: adjoint_eq
  type(gmres)                     :: solver
  type(field)                     :: wrong, state
  class(graph_field), allocatable :: answer
  character(len=32)               :: which

  if (command_argument_count() .lt. 1) then
     error stop 'usage: refusal <case>'
  end if
  call get_command_argument(1, which)

  v = index_set('variables', 3)
  t = index_set('targets'  , 3)
  q_dom = subset('state'   , v, [VAR_U, VAR_V])
  y_dom = subset('residual', t, [TGT_R1, TGT_R2])
  host  = stored_graph(5, tails=[1,2,3,4], heads=[2,3,4,5])

  primal_eq  = opaque_primal(q_dom, y_dom)
  adjoint_eq = opaque_adjoint(y_dom, q_dom)

  select case (trim(which))

  case ('primal-rhs-on-Q')
     ! Residual domain is Y; the right-hand side is offered on Q.
     call solver % attach(primal_eq, host, q_dom)
     wrong = field('rhs on the wrong domain', q_dom)
     call wrong % set_real_vector([8.0_dp, 22.0_dp])
     call solver % apply(host, [wrong], answer)
     write(*,*) 'a right-hand side on Q was accepted by a Y-residual solver'

  case ('adjoint-rhs-on-Y')
     ! Residual domain is Q; the right-hand side is offered on Y.
     call solver % attach(adjoint_eq, host, y_dom)
     wrong = field('rhs on the wrong domain', y_dom)
     call wrong % set_real_vector([1.0_dp, 2.0_dp])
     call solver % apply(host, [wrong], answer)
     write(*,*) 'a right-hand side on Y was accepted by a Q-residual solver'

  case ('primal-state-on-Y')
     state = field('state on the wrong domain', y_dom)
     call state % set_real_vector([2.0_dp, 4.0_dp])
     call primal_eq % apply(host, [state], answer)
     write(*,*) 'the primal equation read a state on Y'

  case ('adjoint-covector-on-Q')
     state = field('covector on the wrong domain', q_dom)
     call state % set_real_vector([-0.4_dp, 0.6_dp])
     call adjoint_eq % apply(host, [state], answer)
     write(*,*) 'the adjoint equation read a covector on Q'

  case default
     error stop 'usage: refusal <case>'

  end select

end program adjoint_level_7_refusal
