!=====================================================================!
! THE VAN DER POL OSCILLATOR AND ITS DESIGN SENSITIVITIES, as the
! derivatives of the solution along design coordinates. The manifold
! is Omega = [0, T] x {(nu, mu)}: the time t and the parameters nu
! and mu are its coordinates, the last two design coordinates, a
! point of the design space rather than an interval. The unknown u
! is a function of all three, u = u(t, nu, mu), stated as such; the
! oscillator
!
!      r(u, nu, mu) = u_tt - nu (1 - u^2) u_t + mu u = 0,
!      u(0) = 2,  u_t(0) = 0,
!
! holds for every (nu, mu), so its derivatives along each design
! coordinate - the tangent equations of every order - are the
! derivatives of the residual along that coordinate, exact, and the
! sensitivities of the energy
!
!      J(nu, mu) = integral over [0, T] of (u^2 + u_t^2) / 2
!
! are the derivatives of J along nu and along mu.
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
! instants, as in van_der_pol.f90; along each design coordinate it is
! the Taylor expansion of an order at the design value, so that the
! discrete solution stores u and d^k u / d nu^k, d^k u / d mu^k, k up
! to the order, at every instant, and the discrete energy stores J
! and its derivatives along each coordinate to the same order. The
! condition on the design factor fixes every design coordinate, one
! equation per coordinate, and its multiplier kappa is the vector
! -(dJ/dnu, dJ/dmu) by the adjoint; the derivative of kappa along a
! coordinate is a row of the Hessian of J, so that the mixed
! derivative is read twice, as -d(kappa_nu)/dmu and -d(kappa_mu)/dnu.
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
! dJ/dnu must equal: kappa = 0.65656201984799478 against dJ/dnu =
! -0.65656201984799489 at 21 instants, order 3. The derivatives of
! kappa along nu, the adjoint's own expansion, give the higher
! derivatives of J by the adjoint: -1.5484062330664192 and
! -1.0211138729502853 against d^2J/dnu^2 = 1.5484062330664188 and
! d^3J/dnu^3 = 1.0211138729502878 by the expansion of the state.
!
! With mu the second coordinate, kappa_mu = -1.4383847623474970
! against dJ/dmu = 1.4383847623474963 by the expansion; the mixed
! derivative -d^2J/dnu dmu read as d(kappa_nu)/dmu = 3.5336347345361871
! and as d(kappa_mu)/dnu = 3.5336347345361867.
!
! The energy is integrated by the schemes' own quadrature: the
! tableau's weights at the stages of a DIRK step, the interpolatory
! rule over the instants of a multistep step. At 21 instants the jets
! agree with central differences of the program run at nu = 1 +-
! 0.01, +- 0.02, to the differences' own truncation error: du(T)/dnu
! = 1.2502856 against 1.2502527, d^2u/dnu^2 = -0.598952 against
! -0.598900, d^3u/dnu^3 = -1.9712 against -1.9710; dJ/dnu =
! -0.6565620 against -0.6565450, d^2J/dnu^2 = 1.548406 against
! 1.548272, d^3J/dnu^3 = 1.0211 against 1.0215.
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

  real(dp), parameter :: T_final = 2.0_dp, nu_design = 1.0_dp, mu_design = 1.0_dp
  real(dp), parameter :: u_initial = 2.0_dp, v_initial = 0.0_dp
  integer , parameter :: num_instants = 21, order = 3

  type(continuous_manifold) :: omega, domega
  type(discrete_manifold)   :: omega_h
  type(continuous_field)    :: t, nu, mu, u, lambda_r, lambda, kappa, oscillator, u0, v0, energy, rest
  type(discrete_field)      :: estimate, solution, u_h, lambda_h, kappa_h, kappa_nu, kappa_mu, sensitivity
  type(continuous_residual) :: r, g, d, L
  type(discrete_residual)   :: L_h
  real(dp)                  :: values(1 + 2 * order)
  integer                   :: k

  ! Omega = [0, T] x {(nu, mu)}, the design coordinates nu and mu at
  ! the values 1; dOmega = {0} x {(nu, mu)}; Omega_h = {t_0..t_n} x the
  ! expansion of order 3 at (1, 1) along each coordinate
  omega   = continuous_manifold(time=interval(0.0_dp, T_final), &
       &                        design=[parameter('nu', nu_design), parameter('mu', mu_design)])
  domega  = omega % boundary(time=0.0_dp)
  omega_h = omega % discretize(time=instants(num_instants), design=expansion(order=order))

  ! the coordinate functions t, nu, mu : Omega -> R
  t  = omega % coordinate('t')
  nu = omega % coordinate('nu')
  mu = omega % coordinate('mu')

  ! the unknown function u(t, nu, mu) : Omega -> R
  u = omega % unknown('u', [t, nu, mu])

  ! the oscillator as a function Omega -> R of the jet of u along t,
  ! reading nu and mu as coordinates
  oscillator = u % derivative([t, t]) - nu*(1.0_dp - u*u)*u % derivative([t]) + mu*u

  ! the initial data, functions dOmega -> R: the value and the first
  ! derivative of u at t = 0
  u0 = u - u_initial
  v0 = u % derivative([t]) - v_initial

  ! the energy J, a function of the design alone: the integral over
  ! the time factor
  energy = integral((u*u + u % derivative([t])**2)/2.0_dp, over=omega % time())

  ! r on Omega with its multiplier lambda_r(t, nu, mu) : Omega -> R,
  ! the adjoint; g on dOmega with its multiplier lambda(nu, mu) :
  ! dOmega -> R^2, one component per condition, a function of the
  ! design alone; the condition d on the design factor {(nu, mu)}
  ! fixing each design coordinate at its value, with its multiplier
  ! kappa, a vector of two; and the Lagrangian L = J + lambda_r . r
  ! + lambda . g + kappa . d, every equality paired
  r        = continuous_residual(omega,  [oscillator])
  g        = continuous_residual(domega, [u0, v0])
  d        = continuous_residual(omega % design(), [nu - nu_design, mu - mu_design])
  lambda_r = r % multiplier('adjoint', [t, nu, mu])
  lambda   = g % multiplier('lambda', [nu, mu])
  kappa    = d % multiplier('kappa')
  L        = energy + lambda_r*r + lambda*g + kappa*d

  ! L_h on Omega_h: d/dt by the chain dirk(2), bdf(2), adams(2), one
  ! block per instant, each reading the instants before it; along nu
  ! the derivatives of the equations are exact, and need no
  ! approximation to be named
  L_h = L % discretize(omega_h, time=chain([(dirk(2), k = 1, 3), (bdf(2), k = 4, 7), &
       &                                    (adams(2), k = 8, num_instants)], from=[(k, k = 1, num_instants)]))

  ! (u, lambda_r, lambda, kappa)_h with their derivatives along nu and
  ! along mu to order 3 = the zero of L_h and of its derivatives along
  ! each coordinate, from the estimate u = u(0) at every instant; the
  ! convergence of every solve is printed by image 1
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
  ! initial position and the initial velocity; the multiplier of the
  ! design conditions, the vector -(dJ/dnu, dJ/dmu) by the adjoint,
  ! with its derivatives along nu: -d^2J/dnu^2, ... and the mixed
  ! derivative -d^2J/dnu dmu read from both rows of the Hessian,
  ! against the derivatives of J by the expansion below
  lambda_h = solution % fields(['lambda'])
  kappa_h  = solution % fields(['kappa'])
  print '(a, es24.16)', '-dJ / du(0)          = ', lambda_h % value(1, 1)
  print '(a, es24.16)', '-dJ / du_t(0)        = ', lambda_h % value(1, 2)
  print '(a, es24.16)', 'kappa_nu = -dJ / d nu = ', kappa_h % value(1, 1)
  print '(a, es24.16)', 'kappa_mu = -dJ / d mu = ', kappa_h % value(1, 2)
  kappa_mu = kappa_h % derivative([mu])
  kappa_nu = kappa_h % derivative([nu])
  print '(a, 2es24.16)', '-d^2 J / d nu d mu = d kappa_nu / d mu, d kappa_mu / d nu = ', &
       & kappa_mu % value(1, 1), kappa_nu % value(1, 2)
  do k = 1, order - 1
     if (k > 1) kappa_nu = kappa_nu % derivative([nu])
     print '(a, i0, a, i0, a, i0, a, es24.16)', 'd^', k, ' kappa_nu / d nu^', k, ' = -d^', k + 1, ' J / d nu = ', &
          & kappa_nu % value(1, 1)
  end do

  ! the energy and its derivatives along nu and along mu: the discrete
  ! energy of the discrete solution, a field on the design factor:
  ! J, then the derivatives along nu of order 1 to 3, then along mu
  sensitivity = energy % at(solution)
  call sensitivity % values(values)
  print '(a, es24.16)', 'J                    = ', values(1)
  do k = 1, order
     print '(a, i0, a, i0, a, es24.16)', 'd^', k, ' J / d nu^', k, '      = ', values(1 + k)
  end do
  do k = 1, order
     print '(a, i0, a, i0, a, es24.16)', 'd^', k, ' J / d mu^', k, '      = ', values(1 + order + k)
  end do

end program van_der_pol_sensitivity
