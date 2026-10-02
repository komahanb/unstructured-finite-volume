!=====================================================================!
! THE STIFFNESS OF A BEAM IDENTIFIED BY THE STATIONARITY OF ONE
! LAGRANGIAN.
!
! The clamped Euler-Bernoulli beam of euler_bernoulli_beam.f90,
!
!      rhoA w_tt + d_xx( EI w_xx ) - f = 0     on Omega = [0, L] x [0, T],
!
! with a stiffness of two design coordinates,
!
!      EI(x) = EI_0 (1 + s1 x / L + s2 (x / L)^2),
!
! and a load f built so that at the true stiffness (s1, s2) = (0.5, -0.3)
! the deflection is w_true = sin(pi x / L)^2 cos(omega t). The design
! is sought as the minimizer of the misfit
!
!      J(s) = integral over Omega of (w - w_true)^2 / 2.
!
! The Lagrangian pairs every law with its multiplier. No condition
! fixes the design:
!
!      L = J + <adjoint, r> + <lambda, ic> + <lambda_0, bc_0>
!            + <lambda_L, bc_L>
!
! Its stationarity in the multipliers is the beam, in the state the
! adjoint, and in the design dJ/ds = 0. With the condition
! s - s_given = 0 added, paired with kappa, the same Lagrangian is the
! sensitivity at s_given, kappa = -dJ/ds, and the design does not
! move. Without it, minimize finds the design at which that reaction
! vanishes, by Newton's iteration on s from the value the manifold
! states, s = 0: at each step the beam with the design held at the
! iterate, one forward and one reverse sweep, the gradient from kappa
! and the hessian from the expansion of order one along each design
! coordinate. The objective is a choice of this example: an
! identification, which has an optimum at the true design up to the
! error of the discretization.
!
! Measured, the iterations printed by minimize at verbosity one: from
! s = 0 the gradient falls as 1.19e-1, 2.19e-2, 1.55e-3, 2.66e-5,
! 4.56e-7, 3.10e-11 (quadratic, the derivatives exact), and the design
! found approaches the true one as the grid is refined: s1 =
! 0.50961827, 0.50205024, 0.50054911 on 21 x 17, 41 x 33 and 81 x 33
! instants by points. At the design found at 21 x 17, dJ/ds by the
! expansion of the state is (-2.47e-11, -1.82e-11).
!=====================================================================!

program euler_bernoulli_beam_optimization

  use util_precision             , only : dp
  use operation_manifold         , only : continuous_manifold, discrete_manifold, interval, region, points, parameter
  use operation_expansion        , only : expansion
  use operation_field            , only : continuous_field, discrete_field, integral
  use operation_field            , only : sin, cos, operator(+), operator(-), operator(*), operator(/), operator(**)
  use operation_residual         , only : continuous_residual, discrete_residual, operator(+), operator(*)
  use operation_family           , only : dirk, bdf, adams, chain
  use operation_finite_difference, only : finite_difference

  implicit none

  real(dp), parameter :: pi = acos(-1.0_dp)
  real(dp), parameter :: L = 2.0_dp, T_final = 1.0_dp, EI_0 = 2.0_dp, rhoA_0 = 1.5_dp
  real(dp), parameter :: s_true(2) = [0.5_dp, -0.3_dp]

  integer  :: num_instants = 21, num_points = 17
  character(len=32) :: item

  type(continuous_manifold) :: omega, initial, left, right
  type(discrete_manifold)   :: omega_h
  type(continuous_field)    :: t, x, s1, s2, w, w_true, exact
  type(continuous_field)    :: EI, EI_true, rhoA, moment, moment_true, load, misfit
  type(continuous_field)    :: adjoint, lambda, lambda_0, lambda_L
  type(discrete_field)      :: estimate, solution, functional, s_h
  type(continuous_residual) :: r, ic, bc_0, bc_L, Lagrangian
  type(discrete_residual)   :: L_h
  real(dp) :: frequency, values(3)
  integer  :: k

  if (command_argument_count() >= 1) then
     call get_command_argument(1, item)
     read(item, *) num_instants
  end if
  if (command_argument_count() >= 2) then
     call get_command_argument(2, item)
     read(item, *) num_points
  end if

  ! Omega = [0, L] x [0, T] x {(s1, s2)}, the design from (0, 0); its
  ! faces at t = 0, x = 0 and x = L
  omega   = continuous_manifold(time=interval(0.0_dp, T_final), space=region(interval(0.0_dp, L)), &
       &                        design=[parameter('s1', 0.0_dp), parameter('s2', 0.0_dp)])
  t       = omega % coordinate('t')
  x       = omega % coordinate('x')
  s1      = omega % coordinate('s1')
  s2      = omega % coordinate('s2')
  initial = omega % boundary(t, 0.0_dp)
  left    = omega % boundary(x, 0.0_dp)
  right   = omega % boundary(x, L)
  omega_h = omega % discretize(time=points(num_instants), space=points(num_points), design=expansion(order=1))

  ! the stiffness at the design, the true stiffness, the mass per length
  rhoA    = rhoA_0 + 0.0_dp * x
  EI      = EI_0 * (1.0_dp + s1 * x / L + s2 * (x / L)**2)
  EI_true = EI_0 * (1.0_dp + s_true(1) * x / L + s_true(2) * (x / L)**2)

  ! the deflection at the true stiffness, and the load it needs
  frequency   = (pi / L)**2 * sqrt(EI_0 / rhoA_0)
  w_true      = sin(pi / L * x)**2 * cos(frequency * t)
  exact       = continuous_field(omega, [w_true])
  moment_true = EI_true * w_true % derivative([x, x])
  load        = rhoA * w_true % derivative([t, t]) + moment_true % derivative([x, x])

  ! the physics in flux form, the initial data and the clamped ends;
  ! no condition on the design
  w      = omega % unknown('w')
  moment = EI * w % derivative([x, x])
  r      = continuous_residual(omega, [rhoA * w % derivative([t, t]) + moment % derivative([x, x]) - load])
  ic     = continuous_residual(initial, [w - w_true, w % derivative([t])])
  bc_0   = continuous_residual(left , [w, w % derivative([x])])
  bc_L   = continuous_residual(right, [w, w % derivative([x])])

  ! the misfit, and the Lagrangian
  misfit   = integral((w - w_true)**2 / 2.0_dp, over=omega)
  adjoint  = r    % multiplier('adjoint')
  lambda   = ic   % multiplier('lambda')
  lambda_0 = bc_0 % multiplier('lambda_0')
  lambda_L = bc_L % multiplier('lambda_L')
  Lagrangian = misfit + adjoint*r + lambda*ic + lambda_0*bc_0 + lambda_L*bc_L

  L_h = Lagrangian % discretize(omega_h, &
       & time=chain([(dirk(2), k = 1, 4), (bdf(2), k = 5, 7), (adams(2), k = 8, num_instants)], &
       &            from=[(k, k = 1, num_instants)]), space=finite_difference(order=6))

  estimate = exact % discretize(omega_h)
  call L_h % minimize(estimate, solution)

  ! J at the design found, dJ/ds there by the expansion of the state,
  ! and the design
  functional = misfit % at(solution)
  call functional % values(values)
  if (this_image() == 1) print '(a, 3es12.4)', 'J, dJ/ds    = ', values
  s_h = solution % fields(['s1', 's2'])
  if (this_image() == 1) print '(a, 2f12.8)', 's found     = ', s_h % value(1, 1), s_h % value(1, 2)
  if (this_image() == 1) print '(a, 2f12.8)', 'true s      = ', s_true

end program euler_bernoulli_beam_optimization
