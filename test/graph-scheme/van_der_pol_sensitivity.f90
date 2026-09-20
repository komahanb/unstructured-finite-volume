!=====================================================================!
! THE VAN DER POL OSCILLATOR AND ITS DESIGN SENSITIVITIES, as the
! derivatives of the solution along a coordinate. The manifold is
! Omega = [0, T] x {nu}: the time t and the parameter nu are its two
! coordinates, the second a design coordinate, a point of the design
! space rather than an interval. The unknown u is a function of both,
! u = u(t, nu), stated as such; the oscillator
!
!      r(u, nu) = u_tt - nu (1 - u^2) u_t + u = 0,   u(0) = 2,  u_t(0) = 0,
!
! holds for every nu, so its derivatives along nu - the tangent
! equations of every order - are the derivatives of the residual along
! that coordinate, exact, and the sensitivities of the energy
!
!      J(nu) = integral over [0, T] of (u^2 + u_t^2) / 2
!
! are the derivatives of J along nu.
!
! THE LAGRANGIAN pairs every equality with a multiplier:
!
!      L = J  +  lambda_r . r on Omega  +  lambda . g on dOmega
!
! Its stationarity in lambda_r is the oscillator, in lambda the
! initial data, and in u the adjoint equation, whose solution
! lambda_r(t) is the sensitivity of J to a forcing of the oscillator
! at the instant t; the components of lambda are the sensitivities of
! J to the initial position and the initial velocity. Without J the
! multipliers are zero and the zero of L is the forward solution.
!
! Along t the discretization is a chain of families over the
! instants, as in van_der_pol.f90; along nu it is the Taylor expansion
! of an order at the design value, so that the discrete solution
! carries u and d^k u / d nu^k, k up to the order, at every instant,
! and the discrete energy carries J and its derivatives to the same
! order.
!
! The statement of every dependence is explicit: an unknown lists its
! arguments among the coordinates, a derivative along a coordinate an
! unknown is not a function of is the zero field, and a multiplier is
! declared by the term it pairs with, one component per equation, its
! arguments the coordinates of the term's manifold.
!
! The condition nu - nu_design = 0 on the design factor is paired
! too, with the multiplier kappa: the stationarity of L in nu gives
! kappa = -dJ/dnu by the adjoint applied to the
! partial derivative of the equations in nu, which the expansion's
! dJ/dnu must equal: kappa = 0.65683913482447476 against dJ/dnu =
! -0.65683913482447487 at 21 instants, order 3. The derivatives of
! kappa along nu, the adjoint's own expansion, give the higher
! derivatives of J by the adjoint: -1.5485964513054384 and
! -1.0209571320096991 against d^2J/dnu^2 = 1.5485964513054380 and
! d^3J/dnu^3 = 1.0209571320096964 by the expansion of the state.
!
! At 21 instants the jets agree with central differences of the
! program run at nu = 1 +- 0.01, +- 0.02, to the differences' own
! truncation error: du(T)/dnu = 1.2502856 against 1.2502527,
! d^2u/dnu^2 = -0.598952 against -0.598900, d^3u/dnu^3 = -1.9712
! against -1.9710; dJ/dnu = -0.6568391 against -0.6568221,
! d^2J/dnu^2 = 1.548596 against 1.548462, d^3J/dnu^3 = 1.0210
! against 1.0213.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program van_der_pol_sensitivity

  use util_precision    , only : dp
  use util_verbosity    , only : set_verbosity
  use operation_manifold, only : continuous_manifold, discrete_manifold, interval, instants, parameter, expansion
  use operation_field   , only : continuous_field, discrete_field, integral, &
       &                         operator(+), operator(-), operator(*), operator(/), operator(**)
  use operation_residual, only : continuous_residual, discrete_residual, operator(+), operator(*)
  use operation_family  , only : dirk, bdf, adams, chain

  implicit none

  real(dp), parameter :: T_final = 2.0_dp, nu_design = 1.0_dp, u_initial = 2.0_dp, v_initial = 0.0_dp
  integer , parameter :: num_instants = 21, order = 3

  type(continuous_manifold) :: omega, domega
  type(discrete_manifold)   :: omega_h
  type(continuous_field)    :: t, nu, u, lambda_r, lambda, kappa, oscillator, u0, v0, energy, rest
  type(discrete_field)      :: estimate, solution, u_h, lambda_h, kappa_h, sensitivity
  type(continuous_residual) :: r, g, d, L
  type(discrete_residual)   :: L_h
  real(dp)                  :: values(0:order)
  integer                   :: k

  ! Omega = [0, T] x {nu}, the design coordinate nu at the value 1;
  ! dOmega = {0} x {nu}; Omega_h = {t_0..t_n} x the expansion of
  ! order 3 at nu = 1
  omega   = continuous_manifold(time=interval(0.0_dp, T_final), design=parameter('nu', nu_design))
  domega  = omega % boundary(time=0.0_dp)
  omega_h = omega % discretize(time=instants(num_instants), design=expansion(order=order))

  ! the coordinate functions t, nu : Omega -> R
  t  = omega % coordinate('t')
  nu = omega % coordinate('nu')

  ! the unknown function u(t, nu) : Omega -> R
  u = omega % unknown('u', [t, nu])

  ! the oscillator as a function Omega -> R of the jet of u along t,
  ! reading nu as a coordinate
  oscillator = u % derivative([t, t]) - nu*(1.0_dp - u*u)*u % derivative([t]) + u

  ! the initial data, functions dOmega -> R: the value and the first
  ! derivative of u at t = 0
  u0 = u - u_initial
  v0 = u % derivative([t]) - v_initial

  ! the energy J, a function of nu alone: the integral over the time
  ! factor
  energy = integral((u*u + u % derivative([t])**2)/2.0_dp, over=omega % time())

  ! r on Omega with its multiplier lambda_r(t, nu) : Omega -> R, the
  ! adjoint; g on dOmega with its multiplier lambda(nu) : dOmega -> R^2,
  ! one component per condition, a function of the design alone; the
  ! condition d on the design factor {nu} fixing the design at its
  ! value, with its multiplier kappa, a number; and the Lagrangian
  ! L = J + lambda_r . r + lambda . g + kappa d, every equality paired
  r        = continuous_residual(omega,  [oscillator])
  g        = continuous_residual(domega, [u0, v0])
  d        = continuous_residual(omega % design(), [nu - nu_design])
  lambda_r = r % multiplier('adjoint', [t, nu])
  lambda   = g % multiplier('lambda', [nu])
  kappa    = d % multiplier('kappa')
  L        = energy + lambda_r*r + lambda*g + kappa*d

  ! L_h on Omega_h: d/dt by the chain dirk(2), bdf(2), adams(2), one
  ! block per instant, each reading the instants before it; along nu
  ! the derivatives of the equations are exact, and need no
  ! approximation to be named
  L_h = L % discretize(omega_h, time=chain([(dirk(2), k = 1, 3), (bdf(2), k = 4, 7), &
       &                                    (adams(2), k = 8, num_instants)], from=[(k, k = 1, num_instants)]))

  ! (u, lambda_r, lambda)_h with their derivatives along nu to order 3
  ! = the zero of L_h and of its derivatives along nu, from the
  ! estimate u = u(0) at every instant; the convergence of every
  ! solve is printed by image 1
  rest     = continuous_field(omega, [u_initial])
  estimate = rest % discretize(omega_h)
  if (this_image() == 1) call set_verbosity(1)
  call L_h % minimize(estimate, solution)

  ! u_h and its derivatives along nu at the last instant
  u_h = solution % fields(['u'])
  print '(a, es24.16)', 'u(T)                 = ', u_h % value(num_instants, 1)
  do k = 1, order
     u_h = u_h % derivative([nu])
     print '(a, i0, a, i0, a, es24.16)', 'd^', k, ' u(T) / d nu^', k, ' = ', u_h % value(num_instants, 1)
  end do

  ! the multiplier of the initial data: the sensitivities of J to the
  ! initial position and the initial velocity; and the multiplier of
  ! the design condition with its derivatives along nu: -dJ/dnu,
  ! -d^2J/dnu^2, ... by the adjoint and its expansion, against the
  ! derivatives of J by the expansion below
  lambda_h = solution % fields(['lambda'])
  kappa_h  = solution % fields(['kappa'])
  print '(a, es24.16)', '-dJ / du(0)          = ', lambda_h % value(1, 1)
  print '(a, es24.16)', '-dJ / du_t(0)        = ', lambda_h % value(1, 2)
  print '(a, es24.16)', 'kappa = -dJ / d nu   = ', kappa_h % value(1, 1)
  do k = 1, order - 1
     kappa_h = kappa_h % derivative([nu])
     print '(a, i0, a, i0, a, i0, a, es24.16)', 'd^', k, ' kappa / d nu^', k, ' = -d^', k + 1, ' J / d nu = ', &
          & kappa_h % value(1, 1)
  end do

  ! the energy and its derivatives along nu: the discrete energy of
  ! the discrete solution, a field on the design coordinate alone
  sensitivity = energy % at(solution)
  call sensitivity % values(values)
  print '(a, es24.16)', 'J                    = ', values(0)
  do k = 1, order
     print '(a, i0, a, i0, a, es24.16)', 'd^', k, ' J / d nu^', k, '      = ', values(k)
  end do

end program van_der_pol_sensitivity
