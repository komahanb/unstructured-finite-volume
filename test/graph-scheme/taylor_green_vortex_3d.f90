!=====================================================================!
! THE TAYLOR-GREEN VORTEX IN THREE DIMENSIONS: the statement of
! taylor_green_vortex.f90 on the periodic box [0, 2 pi]^3, with the
! velocity (u, v, w) and the pressure p as the unknowns. The exact
! solution is the two-dimensional vortex extended along z,
!
!      u* =  sin x cos y e^(-2 nu t)      v* = -cos x sin y e^(-2 nu t)
!      w* =  0                             p* = (cos 2x + cos 2y) e^(-4 nu t) / 4,
!
! which satisfies the three-dimensional equations exactly, so that
! the error of the discrete solution is measured against it as in
! two dimensions. The residual is the three momentum equations and
! the pressure equation
!
!      laplacian p + u_x^2 + v_y^2 + w_z^2
!                  + 2 (u_y v_x + u_z w_x + v_z w_y) = 0,
!
! the initial data of (u, v, w) on the face t = 0 paired with a
! multiplier of three components, the gauge on the time factor with
! one, and the Lagrangian L = r + lambda . g + mu gauge minimised
! block by block, one instant per block, by dirk(2).
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program taylor_green_vortex_3d

  use util_precision             , only : dp
  use util_verbosity             , only : set_verbosity
  use operation_manifold         , only : continuous_manifold, discrete_manifold, interval, region, instants, mesh
  use operation_field            , only : continuous_field, discrete_field, integral, sin, cos, exp, &
       &                                  operator(+), operator(-), operator(*), operator(/), operator(**)
  use operation_residual         , only : continuous_residual, discrete_residual, operator(+), operator(*)
  use operation_family           , only : dirk, chain
  use operation_finite_volume    , only : finite_volume
  use view_paraview_writer       , only : paraview

  implicit none

  real(dp), parameter :: T_final = 0.5_dp, nu = 0.01_dp
  integer , parameter :: num_instants = 5

  type(continuous_manifold) :: omega, domega, tau
  type(discrete_manifold)   :: omega_h
  type(continuous_field)    :: t, x, y, z, u_exact, v_exact, w_exact, p_exact, exact
  type(continuous_field)    :: u, v, w, p, lambda, mu, momentum_x, momentum_y, momentum_z, pressure
  type(discrete_field)      :: estimate, solution, error
  type(continuous_residual) :: r, g, gauge, L
  type(discrete_residual)   :: L_h
  integer                   :: k

  ! Omega = [0, T] x B, B the box of box_3d.py; dOmega = {0} x B;
  ! tau = [0, T]; Omega_h = {t_0..t_n} x C, C the cells of box_3d.msh
  omega   = continuous_manifold(time=interval(0.0_dp, T_final), space=region('box_3d.py'))
  domega  = omega % boundary(time=0.0_dp)
  tau     = omega % time()
  omega_h = omega % discretize(time=instants(num_instants), space=mesh('box_3d.msh'))

  ! the coordinate functions t, x, y, z : Omega -> R
  t = omega % coordinate('t')
  x = omega % coordinate('x')
  y = omega % coordinate('y')
  z = omega % coordinate('z')

  ! u*, v*, w*, p* : Omega -> R
  u_exact =  sin(x)*cos(y)*exp(-2.0_dp*nu*t)
  v_exact = -cos(x)*sin(y)*exp(-2.0_dp*nu*t)
  w_exact =  0.0_dp*z
  p_exact = (cos(2.0_dp*x) + cos(2.0_dp*y))*exp(-4.0_dp*nu*t)/4.0_dp
  exact   = continuous_field(omega, [u_exact, v_exact, w_exact, p_exact])

  ! the unknown functions u, v, w, p : Omega -> R; the multipliers
  ! lambda : dOmega -> R^3 of the initial data and mu : [0, T] -> R
  ! of the gauge
  u      = omega % unknown('u')
  v      = omega % unknown('v')
  w      = omega % unknown('w')
  p      = omega % unknown('p')
  lambda = domega % unknown('lambda', components=3)
  mu     = tau % unknown('mu')

  ! the four equations as functions Omega -> R of the jet of (u, v, w, p)
  momentum_x = u % derivative([t]) + u*u % derivative([x]) + v*u % derivative([y]) + w*u % derivative([z]) &
       &     + p % derivative([x]) - nu*(u % derivative([x, x]) + u % derivative([y, y]) + u % derivative([z, z]))
  momentum_y = v % derivative([t]) + u*v % derivative([x]) + v*v % derivative([y]) + w*v % derivative([z]) &
       &     + p % derivative([y]) - nu*(v % derivative([x, x]) + v % derivative([y, y]) + v % derivative([z, z]))
  momentum_z = w % derivative([t]) + u*w % derivative([x]) + v*w % derivative([y]) + w*w % derivative([z]) &
       &     + p % derivative([z]) - nu*(w % derivative([x, x]) + w % derivative([y, y]) + w % derivative([z, z]))
  pressure   = p % derivative([x, x]) + p % derivative([y, y]) + p % derivative([z, z]) &
       &     + u % derivative([x])**2 + v % derivative([y])**2 + w % derivative([z])**2 &
       &     + 2.0_dp*(u % derivative([y])*v % derivative([x]) + u % derivative([z])*w % derivative([x]) &
       &             + v % derivative([z])*w % derivative([y]))

  ! r on Omega, g on dOmega, the gauge on tau, and the Lagrangian
  ! L = r + lambda . g + mu gauge
  r     = continuous_residual(omega,  [momentum_x, momentum_y, momentum_z, pressure])
  g     = continuous_residual(domega, [u - u_exact, v - v_exact, w - w_exact])
  gauge = continuous_residual(tau,    [integral(p, over=omega % space())])
  L     = r + lambda*g + mu*gauge

  ! L_h on Omega_h and dOmega_h: d/dt by dirk(2), one block per
  ! instant, each block reading the instant before it; d/dx, d/dy,
  ! d/dz and their second derivatives by finite volumes of order 2
  L_h = L % discretize(omega_h, time=chain([(dirk(2), k = 1, num_instants)], from=[(k, k = 1, num_instants)]), &
       &                        space=finite_volume(order=2))

  ! (u, v, w, p, lambda, mu)_h = the zero of L_h, from the initial
  ! estimate (u*, v*, w*, p*) restricted to Omega_h and lambda = mu = 0;
  ! the convergence of every Newton step and linear solve is printed
  ! by image 1
  estimate = exact % discretize(omega_h)
  if (this_image() == 1) call set_verbosity(1)
  call L_h % minimize(estimate, solution)

  ! || (u, v, w, p)_h - (u*, v*, w*, p*)|Omega_h ||
  error = solution % fields(['u', 'v', 'w', 'p']) - estimate
  print '(a, es12.4)', 'error against the exact solution ', error % norm()
  call paraview(solution, 'taylor_green_vortex_3d')

end program taylor_green_vortex_3d
