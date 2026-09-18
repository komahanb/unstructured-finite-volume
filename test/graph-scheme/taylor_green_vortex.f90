!=====================================================================!
! THE TAYLOR-GREEN VORTEX.
!
! Let B = (R / 2 pi Z)^2 be the periodic box and Omega = [0, T] x B.
! Let (u, v, p) : Omega -> R^3 satisfy, with nu > 0,
!
!      u_t + u u_x + v u_y + p_x - nu (u_xx + u_yy) = 0
!      v_t + u v_x + v v_y + p_y - nu (v_xx + v_yy) = 0
!      p_xx + p_yy + u_x^2 + 2 u_y v_x + v_y^2     = 0
!
! the third equation being the divergence of the first two under
! u_x + v_y = 0. The equations determine the flow up to its initial
! state and p up to a function of t; two conditions select one
! solution:
!
!      u - u* = 0,   v - v* = 0        on dOmega = {0} x B
!      integral over B of p(t, .) = 0  on tau = [0, T]
!
! where u* = sin x cos y e^(-2 nu t), v* = -cos x sin y e^(-2 nu t),
! p* = (cos 2x + cos 2y) e^(-4 nu t) / 4 is the solution.
!
! With r the residual of the equations on Omega, g the residual of
! the initial data on dOmega, gauge the residual of the mean pressure
! on tau, and the multipliers lambda : dOmega -> R^2, mu : tau -> R,
! the Lagrangian
!
!      L = r on Omega  +  lambda . g on dOmega  +  mu gauge on tau
!
! has as its stationarities the equations on Omega, the two
! conditions, and the reactions lambda and mu added to the equations
! on each multiplier's manifold. The zero of L is the solution together
! with lambda and mu.
!
! Let Omega_h = {t_0, ..., t_n} x C, n = 10, C the cells of a mesh of
! B, and dOmega_h = {t_0} x C. The discrete Lagrangian L_h is L
! evaluated at every point of Omega_h and dOmega_h with d/dt
! approximated by a chain of schemes over the instants, dirk(2) from
! the first, bdf(2) from the fifth, adams(2) from the eighth, and
! d/dx, d/dy, d^2/dx^2, ... by finite differences on each cell's
! neighbourhood, exact on polynomials of degree 2.
! The discrete solution is the zero of L_h; the printed number is
! || (u, v, p)_h - (u*, v*, p*)|Omega_h ||.
!=====================================================================!

program taylor_green_vortex

  use util_precision             , only : dp
  use operation_manifold         , only : continuous_manifold, discrete_manifold, interval, region, instants, mesh
  use operation_field            , only : continuous_field, discrete_field, integral
  use operation_field            , only : sin, cos, exp, operator(+), operator(-), operator(*), operator(/), operator(**)
  use operation_residual         , only : continuous_residual, discrete_residual, operator(+), operator(*)
  use operation_family           , only : dirk, bdf, adams, chain
  use operation_finite_difference, only : finite_difference
  use view_paraview_writer       , only : paraview
  use util_verbosity             , only : set_verbosity
  use view_expression            , only : expression_view

  implicit none

  real(dp), parameter :: T_final = 1.0_dp, nu = 0.01_dp

  type(continuous_manifold) :: omega, domega, tau
  type(discrete_manifold)   :: omega_h

  type(continuous_field)    :: t, x, y, u_exact, v_exact, p_exact, exact
  type(continuous_field)    :: u, v, p, lambda, mu, momentum_x, momentum_y, pressure
  type(discrete_field)      :: estimate, solution, error

  type(continuous_residual) :: r, g, gauge, L
  type(discrete_residual)   :: L_h
  type(expression_view)     :: view
  integer                   :: k, i

  ! Omega = [0, T] x B, B the region of box.geo; dOmega = {0} x B;
  ! tau = [0, T]; Omega_h = {t_0..t_n} x C, C the cells of box.msh
  omega   = continuous_manifold(time=interval(0.0_dp, T_final), space=region('box.geo'))
  domega  = omega % boundary(time=0.0_dp)
  tau     = omega % time()
  omega_h = omega % discretize(time=instants(10), space=mesh('box.msh'))

  ! the coordinate functions t, x, y : Omega -> R
  t = omega % coordinate('t')
  x = omega % coordinate('x')
  y = omega % coordinate('y')

  ! u*, v*, p* : Omega -> R
  u_exact =  sin(x)*cos(y)*exp(-2.0_dp*nu*t)
  v_exact = -cos(x)*sin(y)*exp(-2.0_dp*nu*t)
  p_exact = (cos(2.0_dp*x) + cos(2.0_dp*y))*exp(-4.0_dp*nu*t)/4.0_dp
  exact   = continuous_field(omega, [u_exact, v_exact, p_exact])

  ! the unknown functions u, v, p : Omega -> R; the multipliers
  ! lambda : dOmega -> R^2 of the initial data and mu : [0, T] -> R
  ! of the gauge
  u      = omega % unknown('u')
  v      = omega % unknown('v')
  p      = omega % unknown('p')
  lambda = domega % unknown('lambda', components=2)
  mu     = tau % unknown('mu')

  ! the three equations as functions Omega -> R of the jet of (u, v, p);
  ! derivative([x, x]) is d^2/dx^2, derivative([x, y]) is d^2/dx dy
  momentum_x = u % derivative([t]) + u*u % derivative([x]) + v*u % derivative([y]) + p % derivative([x]) &
       &     - nu*(u % derivative([x, x]) + u % derivative([y, y]))
  momentum_y = v % derivative([t]) + u*v % derivative([x]) + v*v % derivative([y]) + p % derivative([y]) &
       &     - nu*(v % derivative([x, x]) + v % derivative([y, y]))
  pressure   = p % derivative([x, x]) + p % derivative([y, y]) &
       &     + u % derivative([x])**2 + 2.0_dp*u % derivative([y])*v % derivative([x]) + v % derivative([y])**2

  ! r on Omega, g on dOmega, the gauge on tau, and the Lagrangian
  ! L = r + lambda . g + mu gauge
  r     = continuous_residual(omega,  [momentum_x, momentum_y, pressure])
  g     = continuous_residual(domega, [u - u_exact, v - v_exact])
  gauge = continuous_residual(tau,    [integral(p, over=omega % space())])
  L     = r + lambda*g + mu*gauge

  ! L_h on Omega_h and dOmega_h: d/dt by the chain of schemes over the
  ! instants of Omega_h, dirk(2) from the first, bdf(2) from the fifth,
  ! adams(2) from the eighth, each reading the instants before it;
  ! d/dx, d/dy and their second derivatives by finite differences
  ! exact on polynomials of degree 2
  L_h = L % discretize(omega_h, time=chain([dirk(2), bdf(2), adams(2)], from=[1, 5, 8]), &
       &                        space=finite_difference(degree=2))

  ! ONE NEWTON SOLVE PER INSTANT instead of per block: every block
  ! holds one new instant and the history its family reads. The
  ! smallest system is BDF's, one moment of unknowns (the history
  ! instants enter as fixed rows); a DIRK block of one step holds its
  ! stages beside the instant. Both were run on this case: the first
  ! gives 2.3264E-02 in nine solves, the second 2.3266E-02 in ten.
  !
  !   L_h = L % discretize(omega_h, time=chain([dirk(2), (bdf(2), k = 3, 10)], from=[1, (k, k = 3, 10)]), &
  !        &                        space=finite_difference(degree=2))
  !
  !   L_h = L % discretize(omega_h, time=chain([(dirk(2), k = 1, 10)], from=[(k, k = 1, 10)]), &
  !        &                        space=finite_difference(degree=2))
  !
  ! One solve per DIRK stage is not a chain: a block is bounded by
  ! instants, and the stages of a step are solved with it.

  ! the abstract syntax trees: each equation of each term of L as a
  ! formula, q_j the j-th unknown of the term's manifold with its
  ! derivatives as suffixes (q1x = du/dx, q1xx = d^2u/dx^2); then the
  ! rule of L_h, the sum over the equations of Omega weighted by the
  ! row multipliers lambda_j, as a formula and as a tree
  do k = 1, L % num_terms()
     do i = 1, size(L % term(k) % equation)
        view = expression_view(L % term(k) % equation(i) % graph(1))
        print '(a, i0, a, i0, a, a)', 'L term ', k, ' equation ', i, ': ', view % formula()
     end do
  end do
  view = expression_view(L_h % rule)
  print '(a, a)', 'L_h rule: ', view % formula()
  print '(a)', view % tree()

  ! (u, v, p, lambda, mu)_h = the zero of L_h, from the initial
  ! estimate (u*, v*, p*) restricted to Omega_h and lambda = mu = 0;
  ! the convergence of every Newton step and linear solve is printed
  estimate = exact % discretize(omega_h)

  if (this_image() .eq. 1) then
     call set_verbosity(1)
  end if

  call L_h % minimize(estimate, solution)

  ! || (u, v, p)_h - (u*, v*, p*)|Omega_h ||
  error = solution % fields(['u', 'v', 'p']) - estimate
  print '(a, es12.4)', 'error against the exact solution ', error % norm()
  call paraview(solution, 'taylor_green_vortex')

end program taylor_green_vortex
