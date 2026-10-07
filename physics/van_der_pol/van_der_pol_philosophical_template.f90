!=====================================================================!
! THE VAN DER POL OSCILLATOR AS THE TEMPLATE OF A PROBLEM: the
! second-order equation
!
!      u'' - nu (1 - u^2) u' + u = 0,   u(0) = u_0,   u'(0) = v_0,
!
! on the interval [0, T] at nu = nu_0, with the energy
!
!      J = integral over [0, T] of (u^2 + u_t^2) / 2
!
! as the objective. A problem is five statements, each made once, in
! the order in which each reads the ones before it:
!
!   THE FIELD       the unknown u, a function on the manifold
!                   Omega = [0, T] x {nu} of the coordinates t and nu
!   THE RESIDUALS   every law of the problem as a function that
!                   vanishes on a domain of its own: the oscillator r
!                   on Omega, the initial data g on the face t = 0,
!                   the design d = nu - nu_0 on the factor {nu}
!   THE OBJECTIVE   a functional of the field, J
!   THE LAGRANGIAN  the objective and every residual paired with a
!                   multiplier of its own,
!
!                      L = J + lambda_r . r + lambda . g + kappa . d
!
!   THE ZERO        of the discrete image L_h of L
!
! The stationarity of L in the multipliers is the oscillator with its
! data and its design; in u it is the adjoint equation; in nu it is
! the gradient. So lambda_r(t) is the sensitivity of J to a forcing at
! the instant t, lambda = -(dJ/du(0), dJ/du_t(0)) the reaction of the
! initial data, and kappa = -dJ/dnu the reaction of the design.
!
! THE DISCRETE IMAGE is chosen coordinate by coordinate, after the
! statement and without a change to it: along t the instants and a
! chain of families, dirk(2) over the first three, bdf(2) over the
! next four, adams(2) to the end; along nu the expansion of an order
! at nu_0. An extension of the problem is such a choice, and restates
! nothing. The expansion of order zero is the oscillator at nu_0
! alone: its state, adjoint, energy and reactions are those of
! van_der_pol.f90 in every digit, as they are at every order above.
! The expansion of order m adds the derivatives of u, of the
! multipliers and of J along nu to that order: at 21 instants and
! order 3 the derivatives of J along nu are those of
! van_der_pol_sensitivity.f90 in every digit.
!
! THE AUDIT is the agreement of the two halves of L: d^k J / dnu^k by
! the expansion of the state, and by the expansion of the adjoint,
! the derivative of order k - 1 of -kappa. At eleven instants and
! order three they are apart by 5.6e-16, 4.4e-16 and 1.1e-15.
!
! The instant count, the initial data and the order of the expansion
! are read from the command line (default 11, 2, 0 and 1).
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program van_der_pol_philosophical_template

  use util_precision     , only : dp
  use operation_manifold , only : continuous_manifold, discrete_manifold, interval, points, parameter
  use operation_expansion, only : expansion
  use operation_field    , only : continuous_field, discrete_field, integral
  use operation_residual , only : continuous_residual, discrete_residual
  use operation_family   , only : dirk, bdf, adams, chain

  implicit none

  real(dp), parameter :: T_final = 2.0_dp, nu_0 = 1.0_dp
  real(dp) :: u_0, v_0
  integer  :: num_instants, order

  type(continuous_manifold) :: omega
  type(discrete_manifold)   :: omega_h
  type(continuous_field)    :: t, nu, u, energy, lambda_r, lambda, kappa, rest
  type(continuous_residual) :: r, g, d, L
  type(discrete_residual)   :: L_h
  type(discrete_field)      :: estimate, solution, u_h, adjoint_h, lambda_h, kappa_h, energy_h
  integer                   :: k

  num_instants = nint(argument(1, 11.0_dp))
  u_0          = argument(2, 2.0_dp)
  v_0          = argument(3, 0.0_dp)
  order        = nint(argument(4, 1.0_dp))

  ! THE FIELD. Omega = [0, T] x {nu}, the design coordinate nu at the
  ! value nu_0; the coordinate functions t, nu : Omega -> R; the
  ! unknown function u(t, nu) : Omega -> R
  omega = continuous_manifold(time=interval(0.0_dp, T_final), design=parameter('nu', nu_0))
  t     = omega % coordinate('t')
  nu    = omega % coordinate('nu')
  u     = omega % unknown('u', [t, nu])

  ! THE RESIDUALS, each on its own domain: on Omega the oscillator, a
  ! function of the jet of u to order two along t, reading nu as a
  ! coordinate; on the face t = 0 the initial data, the value and the
  ! first derivative of u; on the factor {nu} the design, fixed at nu_0
  r = continuous_residual(omega, [u % derivative([t, t]) - nu*(1.0_dp - u*u)*u % derivative([t]) + u])
  g = continuous_residual(omega % boundary(t, 0.0_dp), [u - u_0, u % derivative([t]) - v_0])
  d = continuous_residual(omega % design(), [nu - nu_0])

  ! THE OBJECTIVE. The energy J, the integral over the time factor
  energy = integral((u*u + u % derivative([t])**2)/2.0_dp, over=omega % time())

  ! THE LAGRANGIAN. Every residual declares its multiplier, one
  ! component per equation, a function of the coordinates of its
  ! domain: lambda_r(t, nu) the adjoint, lambda(nu) the reaction of
  ! the initial data, kappa that of the design; and L = J + lambda_r . r
  ! + lambda . g + kappa . d
  lambda_r = r % multiplier('adjoint', [t, nu])
  lambda   = g % multiplier('lambda', [nu])
  kappa    = d % multiplier('kappa')
  L        = energy + lambda_r*r + lambda*g + kappa*d

  ! THE DISCRETE IMAGE, coordinate by coordinate: along t the instants
  ! {t_1..t_n} and the chain dirk(2), bdf(2), adams(2), one block per
  ! instant, each block reading the instants before it; along nu the
  ! expansion of the order at nu_0, the derivatives of the equations
  ! along nu exact
  omega_h = omega % discretize(time=points(num_instants), design=expansion(order=order))
  L_h     = L % discretize(omega_h, time=chain([(dirk(2), k = 1, 3), (bdf(2), k = 4, 7), &
       &                                        (adams(2), k = 8, num_instants)], from=[(k, k = 1, num_instants)]))

  ! THE ZERO of L_h, from the estimate u = u_0 at every instant: u,
  ! lambda_r, lambda and kappa at the points, each with its
  ! derivatives along nu to the order
  rest     = continuous_field(omega, [u_0])
  estimate = rest % discretize(omega_h)
  call L_h % minimize(estimate, solution)

  ! the state and the adjoint at every tenth of the interval
  u_h       = solution % fields(['u'])
  adjoint_h = solution % fields(['adjoint'])
  do k = 1, num_instants
     if (mod(k - 1, max(1, (num_instants - 1) / 10)) == 0) then
        print '(a, f8.4, a, es16.8, a, es16.8)', 't = ', omega_h % instant(k), '   u = ', u_h % value(k, 1), &
             & '   adjoint = ', adjoint_h % value(k, 1)
     end if
  end do

  ! the energy of the discrete solution, a field on the design factor:
  ! J, then its derivatives along nu of order 1 to the order; and the
  ! reactions, -dJ/du(0), -dJ/du_t(0) and -dJ/dnu
  energy_h = energy % at(solution)
  lambda_h = solution % fields(['lambda'])
  kappa_h  = solution % fields(['kappa'])
  print '(a, es24.16)', 'J             = ', energy_h % value(1, 1)
  print '(a, es24.16)', '-dJ / du(0)   = ', lambda_h % value(1, 1)
  print '(a, es24.16)', '-dJ / du_t(0) = ', lambda_h % value(1, 2)
  print '(a, es24.16)', '-dJ / dnu     = ', kappa_h % value(1, 1)

  ! THE AUDIT. d^k J / dnu^k by the expansion of the state, against
  ! the derivative of order k - 1 of -kappa, the expansion of the
  ! adjoint
  if (order >= 1) print '(a)', 'k, d^k J / dnu^k by the expansion of the state, by that of the adjoint, their difference'
  do k = 1, order
     if (k > 1) kappa_h = kappa_h % derivative([nu])
     print '(i3, 2es24.16, es12.2)', k, energy_h % value(1, 1 + k), -kappa_h % value(1, 1), &
          & energy_h % value(1, 1 + k) + kappa_h % value(1, 1)
  end do

contains

  ! the k-th number of the command line, or the default without one
  real(dp) function argument(k, default)
    integer , intent(in) :: k
    real(dp), intent(in) :: default
    character(len=32) :: item
    argument = default
    if (command_argument_count() >= k) then
       call get_command_argument(k, item)
       read(item, *) argument
    end if
  end function argument

end program van_der_pol_philosophical_template
