!=====================================================================!
! THE AUDIT OF THE EXPANSION ALONG THE DESIGN, on the van der Pol
! oscillator with two design coordinates,
!
!      u'' - nu (1 - u^2) u' + mu u = 0,   u(0) = 2,   u'(0) = nu / 4,
!      J = integral over [0, T] of (u^2 + mu u_t^2) / 2,
!
! the initial velocity and the objective reading the design. The
! expansion is of order two along nu, along mu and along the
! direction w = nu + mu. Three laws are audited, each by a difference
! that is zero to rounding, and the program stops when one exceeds
! 1e-10.
!
! THE TWO HALVES OF L AGREE. Along each expansion, of seed s, the
! derivatives of J by the expansion of the state against those by the
! adjoint, dJ = -kappa . s and d^2J = -(dkappa) . s. The multiplier
! of the design condition takes, beside the equations' products, the
! reaction of the initial data times the derivative of the datum:
! without that term dJ/dnu is 0.10 in error, lambda_2 / 4.
!
! THE HESSIAN IS SYMMETRIC. The mixed derivative -d^2J/dnu dmu read
! as dkappa_nu/dmu and as dkappa_mu/dnu, and against the polarization
! of the expansions, -(J_ww - J_nunu - J_mumu) / 2.
!
! A CONDITION NEEDS NO MULTIPLIER. The same Lagrangian with the
! initial data and the design condition summed bare has the same
! state, adjoint and energy, and its solution has no multiplier of
! either among its fields.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program van_der_pol_design_audit

  use util_precision     , only : dp
  use operation_manifold , only : continuous_manifold, discrete_manifold, interval, points, parameter
  use operation_expansion, only : expansion
  use operation_field    , only : continuous_field, discrete_field, integral, &
       &                          operator(+), operator(-), operator(*), operator(/), operator(**)
  use operation_residual , only : continuous_residual, discrete_residual, operator(+), operator(*)
  use operation_family   , only : dirk, bdf, adams, chain

  implicit none

  real(dp), parameter :: T_final = 2.0_dp, tolerance = 1.0e-10_dp
  integer , parameter :: num_instants = 21, order = 2
  ! the seeds of the three expansions: along nu, along mu, along nu + mu
  real(dp), parameter :: seed(2, 3) = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], [2, 3])

  type(continuous_manifold) :: omega
  type(discrete_manifold)   :: omega_h
  type(continuous_field)    :: t, nu, mu, w, u, energy, lambda_r, lambda, kappa, rest
  type(continuous_residual) :: r, g, d, L
  type(discrete_residual)   :: L_h
  type(discrete_field)      :: estimate, solution, bare, kappa_h, slope, energy_h, bare_energy
  type(continuous_field)    :: along(3)
  real(dp) :: J(1 + 3 * order), by_adjoint, mixed, apart, largest
  integer  :: k, e

  ! THE FIELD on Omega = [0, T] x {(nu, mu)}
  omega = continuous_manifold(time=interval(0.0_dp, T_final), design=[parameter('nu', 1.0_dp), parameter('mu', 1.0_dp)])
  t  = omega % coordinate('t')
  nu = omega % coordinate('nu')
  mu = omega % coordinate('mu')
  w  = omega % direction(seed(:, 3))
  u  = omega % unknown('u', [t, nu, mu])

  ! THE RESIDUALS: the oscillator, the initial data whose velocity
  ! reads nu, the design fixed at (1, 1); THE OBJECTIVE, reading mu
  r = continuous_residual(omega, [u % derivative([t, t]) - nu*(1.0_dp - u*u)*u % derivative([t]) + mu*u])
  g = continuous_residual(omega % boundary(t, 0.0_dp), [u - 2.0_dp, u % derivative([t]) - 0.25_dp*nu])
  d = continuous_residual(omega % design(), [nu - 1.0_dp, mu - 1.0_dp])
  energy = integral((u*u + mu*u % derivative([t])**2)/2.0_dp, over=omega % time())

  ! THE LAGRANGIAN, every condition paired
  lambda_r = r % multiplier('adjoint', [t, nu, mu])
  lambda   = g % multiplier('lambda', [nu, mu])
  kappa    = d % multiplier('kappa')
  L        = energy + lambda_r*r + lambda*g + kappa*d

  ! THE DISCRETE IMAGE AND THE ZERO
  omega_h  = omega % discretize(time=points(num_instants), design=expansion(order, along=seed))
  L_h      = L % discretize(omega_h, time=chain([(dirk(2), k = 1, 3), (bdf(2), k = 4, 7), &
       &                                         (adams(2), k = 8, num_instants)], from=[(k, k = 1, num_instants)]))
  rest     = continuous_field(omega, [2.0_dp])
  estimate = rest % discretize(omega_h)
  call L_h % minimize(estimate, solution)

  energy_h = energy % at(solution)
  call energy_h % values(J)
  kappa_h  = solution % fields(['kappa'])
  along(1) = nu
  along(2) = mu
  along(3) = w
  largest  = 0.0_dp

  ! THE TWO HALVES OF L AGREE, along each expansion and at each order
  print '(a)', 'expansion, order, the derivative of J by the state, by the adjoint, their difference'
  do e = 1, 3
     by_adjoint = -dot_product(kappa_h % value(1, :), seed(:, e))
     apart      = J(1 + (e - 1) * order + 1) - by_adjoint
     largest    = max(largest, abs(apart))
     print '(2i3, 2es24.16, es12.2)', e, 1, J(1 + (e - 1) * order + 1), by_adjoint, apart
     slope      = kappa_h % derivative([along(e)])
     by_adjoint = -dot_product(slope % value(1, :), seed(:, e))
     apart      = J(1 + (e - 1) * order + 2) - by_adjoint
     largest    = max(largest, abs(apart))
     print '(2i3, 2es24.16, es12.2)', e, 2, J(1 + (e - 1) * order + 2), by_adjoint, apart
  end do

  ! THE HESSIAN IS SYMMETRIC
  mixed = (J(1 + 2 * order + 2) - J(1 + 2) - J(1 + order + 2)) / 2.0_dp
  slope = kappa_h % derivative([mu])
  print '(a, es24.16)', '-d kappa_nu / d mu                = ', -slope % value(1, 1)
  apart = mixed + slope % value(1, 1)
  slope = kappa_h % derivative([nu])
  print '(a, es24.16)', '-d kappa_mu / d nu                = ', -slope % value(1, 2)
  print '(a, es24.16)', '(J_ww - J_nunu - J_mumu) / 2      = ', mixed
  largest = max(largest, abs(apart), abs(mixed + slope % value(1, 2)))

  ! A CONDITION NEEDS NO MULTIPLIER: the initial data and the design
  ! condition bare
  L   = energy + lambda_r*r + g + d
  L_h = L % discretize(omega_h, time=chain([(dirk(2), k = 1, 3), (bdf(2), k = 4, 7), &
       &                                    (adams(2), k = 8, num_instants)], from=[(k, k = 1, num_instants)]))
  call L_h % minimize(estimate, bare)
  bare_energy = energy % at(bare)
  apart = maxval(abs(bare % value(:, 1:2) - solution % value(:, 1:2)))
  print '(a, *(a, 1x))', 'the fields of the solution with the conditions bare: ', &
       & (trim(bare % name(k)), k = 1, size(bare % name))
  print '(a, es12.2, a, es12.2)', 'its state and adjoint apart by ', apart, ', its energy by ', &
       & abs(bare_energy % value(1, 1) - J(1))
  largest = max(largest, apart, abs(bare_energy % value(1, 1) - J(1)))
  if (size(bare % name) /= 2) error stop 'van_der_pol_design_audit: a bare condition returned a multiplier'

  print '(a, es12.2)', 'the largest difference of the audit = ', largest
  if (largest > tolerance) error stop 'van_der_pol_design_audit: a law of the expansion along the design fails'

end program van_der_pol_design_audit
