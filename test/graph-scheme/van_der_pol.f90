!=====================================================================!
! THE VAN DER POL OSCILLATOR, the second-order equation
!
!      u'' - nu (1 - u^2) u' + u = 0
!
! on the interval [0, T] with u(0) = 2, u'(0) = 0 and nu = 1, stated
! on a manifold of time alone, with the energy
!
!      J = integral over [0, T] of (u^2 + u_t^2) / 2
!
! as the objective. The Lagrangian pairs every equality with a
! multiplier,
!
!      L = J  +  lambda_r . r on [0, T]  +  lambda . g on {0},
!
! r the residual of the equation, g the initial data, lambda_r the
! adjoint and lambda the reaction of the initial data with one
! component per condition. Its stationarity in the multipliers is
! the equation with its data; in u it is the adjoint equation, whose
! solution gives the sensitivities of J: lambda_r(t) to a forcing at
! the instant t, and lambda = -(dJ/du(0), dJ/du_t(0)) to the initial
! data. L is discretized by a heterogeneous chain of families over
! the instants - dirk(2) over the first three, bdf(2) over the next
! four, adams(2) to the end, one block per instant - and minimised
! block by block forward for the state, then in reverse for the
! multipliers. The instant count and the initial data are read from
! the command line (default 11, 2, 0); the solution converges with
! the count: u(2) = 0.32331 at 401 and at 1601 instants.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program van_der_pol

  use util_precision    , only : dp
  use util_verbosity    , only : set_verbosity
  use operation_manifold, only : continuous_manifold, discrete_manifold, interval, instants
  use operation_field   , only : continuous_field, discrete_field, integral, &
       &                         operator(+), operator(-), operator(*), operator(/), operator(**)
  use operation_residual, only : continuous_residual, discrete_residual, operator(+), operator(*)
  use operation_family  , only : dirk, bdf, adams, chain

  implicit none

  real(dp), parameter :: T_final = 2.0_dp, nu = 1.0_dp
  real(dp) :: u_initial = 2.0_dp, v_initial = 0.0_dp
  integer  :: num_instants = 11

  type(continuous_manifold) :: omega, domega
  type(discrete_manifold)   :: omega_h
  type(continuous_field)    :: t, u, lambda_r, lambda, oscillator, u0, v0, energy, rest
  type(discrete_field)      :: estimate, solution, u_h, lambda_h, adjoint_h, functional
  type(continuous_residual) :: r, g, L
  type(discrete_residual)   :: L_h
  real(dp), allocatable     :: values(:)
  real(dp)                  :: J(1)
  character(len=32)         :: item
  integer                   :: k

  if (command_argument_count() >= 1) then
     call get_command_argument(1, item)
     read(item, *) num_instants
  end if
  if (command_argument_count() >= 2) then
     call get_command_argument(2, item)
     read(item, *) u_initial
  end if
  if (command_argument_count() >= 3) then
     call get_command_argument(3, item)
     read(item, *) v_initial
  end if

  ! Omega = [0, T], a manifold of time alone; dOmega = {0};
  ! Omega_h = {t_0..t_n}
  omega   = continuous_manifold(time=interval(0.0_dp, T_final))
  domega  = omega % boundary(time=0.0_dp)
  omega_h = omega % discretize(time=instants(num_instants))

  ! the coordinate function t and the unknown function u(t) : Omega -> R
  t = omega % coordinate('t')
  u = omega % unknown('u', [t])

  ! the equation as a function Omega -> R of the jet of u to order
  ! two along t
  oscillator = u % derivative([t, t]) - nu*(1.0_dp - u*u)*u % derivative([t]) + u

  ! the initial data, functions dOmega -> R: the value and the first
  ! derivative of u at t = 0
  u0 = u - u_initial
  v0 = u % derivative([t]) - v_initial

  ! the energy J, the integral over the time factor
  energy = integral((u*u + u % derivative([t])**2)/2.0_dp, over=omega % time())

  ! r on Omega with its multiplier lambda_r(t) : Omega -> R, the
  ! adjoint; g on dOmega with its multiplier lambda : dOmega -> R^2,
  ! one component per condition, a function on a point; and the
  ! Lagrangian L = J + lambda_r . r + lambda . g
  r        = continuous_residual(omega,  [oscillator])
  g        = continuous_residual(domega, [u0, v0])
  lambda_r = r % multiplier('adjoint', [t])
  lambda   = g % multiplier('lambda')
  L        = energy + lambda_r*r + lambda*g

  ! L_h on Omega_h: d/dt by the chain, one block per instant, each
  ! block reading the instants before it
  L_h = L % discretize(omega_h, time=chain([(dirk(2), k = 1, 3), (bdf(2), k = 4, 7), &
       &                                    (adams(2), k = 8, num_instants)], from=[(k, k = 1, num_instants)]))

  ! (u, lambda_r, lambda)_h = the zero of L_h, from the estimate u(0)
  ! at every instant: the state block by block forward, the
  ! multipliers block by block in reverse; the convergence of every
  ! solve is printed by image 1 where the instants are few
  rest     = continuous_field(omega, [u_initial])
  estimate = rest % discretize(omega_h)
  if (this_image() == 1 .and. num_instants <= 21) call set_verbosity(1)
  call L_h % minimize(estimate, solution)

  ! the state and the adjoint at every tenth of the interval
  allocate(values(num_instants))
  u_h = solution % fields(['u'])
  call u_h % values(values)
  adjoint_h = solution % fields(['adjoint'])
  do k = 1, num_instants
     if (mod(k - 1, max(1, (num_instants - 1) / 10)) == 0) then
        print '(a, f8.4, a, es16.8, a, es16.8)', 't = ', omega_h % instant(k), '   u = ', values(k), &
             & '   adjoint = ', adjoint_h % value(k, 1)
     end if
  end do

  ! the energy of the discrete solution, and the reactions of the
  ! initial data: -dJ/du(0) and -dJ/du_t(0)
  functional = energy % at(solution)
  call functional % values(J)
  lambda_h = solution % fields(['lambda'])
  print '(a, es24.16)', 'J             = ', J(1)
  print '(a, es24.16)', '-dJ / du(0)   = ', lambda_h % value(1, 1)
  print '(a, es24.16)', '-dJ / du_t(0) = ', lambda_h % value(1, 2)

end program van_der_pol
