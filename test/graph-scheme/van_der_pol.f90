!=====================================================================!
! THE VAN DER POL OSCILLATOR, u'' - mu (1 - u^2) u' + u = 0, as the
! first-order system
!
!      u' = v            v' = mu (1 - u^2) v - u,
!
! on the interval [0, T] with u(0) = 2, v(0) = 0 and mu = 1, stated
! on a manifold of time alone: the residual r on [0, T], the initial
! data g on its boundary {0} paired with a multiplier of two
! components, the Lagrangian L = r + lambda . g discretized by a
! heterogeneous chain of families over the instants - dirk(2) over
! the first three, bdf(2) over the next four, adams(2) to the end,
! one block per instant - and minimised block by block. The instant
! count is read from the command line (default 11); the solution
! converges with it: u(2) = 0.32331 at 401 and at 1601 instants.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program van_der_pol

  use util_precision    , only : dp
  use util_verbosity    , only : set_verbosity
  use operation_manifold, only : continuous_manifold, discrete_manifold, interval, instants
  use operation_field   , only : continuous_field, discrete_field, &
       &                         operator(+), operator(-), operator(*), operator(/), operator(**)
  use operation_residual, only : continuous_residual, discrete_residual, operator(+), operator(*)
  use operation_family  , only : dirk, bdf, adams, chain

  implicit none

  real(dp), parameter :: T_final = 2.0_dp, mu = 1.0_dp, u_initial = 2.0_dp, v_initial = 0.0_dp
  integer  :: num_instants = 11

  type(continuous_manifold) :: omega, domega
  type(discrete_manifold)   :: omega_h
  type(continuous_field)    :: t, u, v, lambda, first, second, u0, v0, rest
  type(discrete_field)      :: estimate, solution
  type(continuous_residual) :: r, g, L
  type(discrete_residual)   :: L_h
  real(dp), allocatable     :: values(:)
  character(len=32)         :: item
  integer                   :: k

  if (command_argument_count() >= 1) then
     call get_command_argument(1, item)
     read(item, *) num_instants
  end if

  ! Omega = [0, T], a manifold of time alone; dOmega = {0};
  ! Omega_h = {t_0..t_n}
  omega   = continuous_manifold(time=interval(0.0_dp, T_final))
  domega  = omega % boundary(time=0.0_dp)
  omega_h = omega % discretize(time=instants(num_instants))

  ! the coordinate function t and the unknown functions u(t), v(t) :
  ! Omega -> R; the multiplier lambda : dOmega -> R^2 of the initial
  ! data, a function on a point
  t      = omega % coordinate('t')
  u      = omega % unknown('u', [t])
  v      = omega % unknown('v', [t])
  lambda = domega % unknown('lambda', components=2)

  ! the two equations as functions Omega -> R of the jet of (u, v)
  first  = u % derivative([t]) - v
  second = v % derivative([t]) - (mu*(1.0_dp - u*u)*v - u)

  ! the initial data, functions dOmega -> R
  u0 = u - u_initial
  v0 = v - v_initial

  ! r on Omega, g on dOmega, and the Lagrangian L = r + lambda . g
  r = continuous_residual(omega,  [first, second])
  g = continuous_residual(domega, [u0, v0])
  L = r + lambda*g

  ! L_h on Omega_h: d/dt by the chain, one block per instant, each
  ! block reading the instants before it
  L_h = L % discretize(omega_h, time=chain([(dirk(2), k = 1, 3), (bdf(2), k = 4, 7), &
       &                                    (adams(2), k = 8, num_instants)], from=[(k, k = 1, num_instants)]))

  ! (u, v, lambda)_h = the zero of L_h, from the estimate (u(0), v(0))
  ! at every instant; the convergence of every solve is printed by
  ! image 1 where the instants are few
  rest     = continuous_field(omega, [u_initial, v_initial])
  estimate = rest % discretize(omega_h)
  if (this_image() == 1 .and. num_instants <= 21) call set_verbosity(1)
  call L_h % minimize(estimate, solution)

  ! the solution at every tenth of the interval
  allocate(values(2 * num_instants))
  call solution % values(values)
  do k = 1, num_instants
     if (mod(k - 1, max(1, (num_instants - 1) / 10)) == 0) then
        print '(a, f8.4, a, 2es16.8)', 't = ', omega_h % instant(k), '   u, v = ', values(2*k-1), values(2*k)
     end if
  end do

end program van_der_pol
