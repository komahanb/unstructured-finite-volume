!=====================================================================!
! The two-stage, third-order SDIRK coefficient system can be written
! as this abstract syntax tree (AST).
!
! For two-stage DIRK, the unknown vertices are
! \[
!   W_{11},W_{21},W_{22},B_1,B_2,\delta_1,\delta_2,
! \]
! and the prescribed input is \(h\).
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_graph_algebra

  use util_precision
  use operation_expression
  use view_expression,        only : expression_view
  use operation_residual,     only : residual_operator
  use operation_stencil,      only : stencil
  use operation_newton,       only : newton
  use operation_dense_direct, only : dense_direct
  use operation_minimization, only : solve_result
  use view_directed_stored,   only : stored_directed_graph
  use field_stored,           only : stored_field

  implicit none

  type(expression) :: W11, W21, W22, B1, B2, delta1, delta2
  type(expression) :: h, C(7)

  type(expression)            :: L
  type(residual_operator)     :: equations
  type(stencil)               :: zero_stencil
  type(newton)                :: solver
  type(solve_result)          :: outcome
  type(stored_directed_graph) :: domain, parameter_domain
  type(stored_field)          :: step
  type(expression_view)       :: view

  real(dp) :: theta(7), rhs(7), dt, achieved
  character(len=6), parameter :: names(7) = &
       & [character(len=6) :: 'W11', 'W21', 'W22', 'B1', 'B2', 'delta1', 'delta2']
  integer :: i

  W11    = unknown(1)
  W21    = unknown(2)
  W22    = unknown(3)
  B1     = unknown(4)
  B2     = unknown(5)
  delta1 = unknown(6)
  delta2 = unknown(7)

  h = design() ! supplied physical step size

  C(1) = W11 - delta1
  C(2) = W21 + W22 - delta2
  C(3) = B1 + B2 - h
  C(4) = B1*delta1 + B2*delta2 - h**2/2.0_dp
  C(5) = B1*delta1**2 + B2*delta2**2 - h**3/3.0_dp
  C(6) = B1*W11*delta1 + B2*(W21*delta1 + W22*delta2) - h**3/6.0_dp
  C(7) = W11 - W22


  ! Lagrangian: each constraint C(i) paired with multiplier(i)
  L = multiplier(1)*C(1)
  do i = 2, 7
     L = L + multiplier(i)*C(i)
  end do

  view = expression_view(L)
  print '(a)', view % tree()
  print '(a)', view % diagram()
  print '(a)', view % formula()
  print '(a)', view % brackets()
  print '(a)', view % digraph()

  ! Seven algebraic unknowns; no differential stencil contribution.
  zero_stencil = stencil( &
       [integer ::], [integer ::], [real(dp) ::], &
       spread(0.0_dp, 1, 7), 'zero linear contribution')

  equations = residual_operator( &
       zero_stencil, L, at=[0], unknowns=7, degrees=7, &
       primary=[0,1,2,3,4,5,6], &
       fixed_rows=[integer ::], fixed=[real(dp) ::])

  domain = equations % unknown_graph()
  parameter_domain = equations % point_domain()

  ! Supplied physical time step.
  dt = 0.1_dp
  step = stored_field('physical step', parameter_domain % vertex_set(), 1)
  call step % set_real_vector([dt])

  ! Existing nonlinear and linear solvers.
  allocate(solver % inner, source=dense_direct())
  solver % tolerance = 1.0e-12_dp
  solver % max_iterations = 30

  call solver % state(equations, domain, domain % vertex_set(), 7, &
       stored_inputs=[step])

  ! Initial estimate selects the Crouzeix root.
  theta = dt * [0.8_dp, -0.6_dp, 0.8_dp, &
       0.5_dp, 0.5_dp, 0.8_dp, 0.2_dp]
  rhs = 0.0_dp

  call solver % solve(rhs, theta, achieved)

  outcome = solver % result()
  if (.not. outcome % converged()) then
     print *, outcome % description()
     error stop 'coefficient solve did not converge'
  end if

  do i = 1, 7
     write(*, '(a6, " = ", es24.16)') names(i), theta(i)
  end do

end program test_graph_algebra
