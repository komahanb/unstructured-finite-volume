!=====================================================================!
! THE EULER-BERNOULLI BEAM.
!
! Let Omega = [0, L] x [0, T], the beamline times the timeline, and
! let w : Omega -> R be the deflection of a beam of bending stiffness
! EI(x) and mass per length rhoA(x) under the load f(t, x):
!
!      rhoA w_tt + d_xx( EI w_xx ) - f = 0          on Omega
!
! with the initial deflection and velocity on {t = 0}, and two
! conditions at each end of the beamline on {x = 0} and {x = L}:
!
!      simply supported   w = 0,  w_xx = 0
!      clamped            w = 0,  w_x  = 0
!
! Every law is a residual on its own manifold, paired with its own
! multiplier. The derivatives of w are components of its jet, and the
! representation laws state what each component is,
!
!      r_t = ( d_t w - w_t,  d_t w_t - w_tt )                  on Omega
!      r_x = ( d_x w - w_x,  ..,  d_x w_xxx - w_xxxx )          on Omega
!
! so that with J the energy integrated over Omega the Lagrangian is
!
!      L = J + <adjoint, r> + <mu_t, r_t> + <mu_x, r_x>
!            + <lambda, ic> + <lambda_0, bc_0> + <lambda_L, bc_L>
!
! Omega_h is n instants of [0, T] times np points of [0, L], ends
! included; no mesh file is read. d_t is discretized by the chain of
! families of the other examples, d_x by finite differences of order
! 6, which give w_xxxx at order 4, or by finite elements. A condition
! at an end of the beamline takes, under finite differences, the row
! of the physics at a point next to that end, so the physics holds on
! the points away from the ends; under finite elements, the row of the
! coefficient it fixes at the end.
!
! Eight cases:
!
!   1. simply supported, uniform EI and rhoA, no load: the mode
!      w = sin(pi x / L) cos(omega t), omega = (pi / L)^2 sqrt(EI / rhoA),
!      whose energy is EI (pi / L)^4 L / 4 at every instant, so that
!      J = T EI (pi / L)^4 L / 4;
!   2. clamped, EI(x) = EI_0 (1 + x / (2 L)) and
!      rhoA(x) = rhoA_0 (1 + cos(pi x / L) / 4), loaded so that
!      w = sin(pi x / L)^2 cos(omega t) is the solution; the load is
!      the law applied to that w, the derivative of EI w_xx formed by
!      the chain rule when the law is discretized;
!   3. case 2 on the open line: d_x by the trapezoid rule on the
!      segments of the line of points, whose relations leave four rows
!      open at the lower end for the four conditions of both ends; of
!      order two, and stable where the fitted differences of order 6,
!      one-sided at the ends, are not (801 instants, 129 points: case 1
!      diverges at order 6, while the open line gives 8.64e-5 and
!      case 3 6.53e-4, falling by 4.0 per halving of the step);
!   4. case 1 on the open line;
!   5. case 2 with the Kelvin-Voigt damping of the bending moment,
!
!      rhoA w_tt + d_xx( EI w_xx + c w_xxt ) - f = 0,
!
!      c(x) = c_0 (1 + (x / L)^2): the law reads the mixed components
!      w_xxt, w_xxxt, w_xxxxt, each the difference along x of w_t, the
!      path the concretion chooses. Measured, the error falls as
!      3.06e-2, 4.07e-4, 9.03e-6 over 9 to 33 points at 801 instants,
!      and as 8.91e-4, 1.51e-4, 3.13e-5, 8.33e-6 over 21 to 161 instants
!      at 65 points; one Newton step per block; with c_0 the design, the
!      misfit's dJ/dc_0 at c_0 = 0.5 on 21 x 17 is 0.1096767998 by the
!      tangent, 0.1096767982 by the adjoint, 0.1096764971 by the central
!      difference of step 1e-3 (its error h^2/6 J''' = -3e-7);
!   6. case 2 by finite elements: the law in flux form tested against
!      the Hermite cubics on the segments of the line, the derivative
!      d_xx of the bending moment moved onto the multiplier, so that
!      the elements read w_xx alone; each point stores the coefficients
!      w and w_x, the clamped ends fix them, and the initial data are
!      projected onto the cubics that satisfy the clamped ends.
!      Measured at 801 instants, the error falls as 1.37e-4, 8.47e-6,
!      4.70e-7 over 9 to 33 points, by 16 per halving: order 4; at 65
!      points as 1.69e-3, 2.50e-4, 4.14e-5, 8.15e-6 over 21 to 161
!      instants;
!   7. case 5 by finite elements, the damping's w_xxt formed by the
!      basis from the jets in time of w and w_x: 4.52e-5, 3.49e-6,
!      3.57e-7 over 9 to 33 points at 801 instants;
!   8. case 6 on a hub that turns: the hub's angle theta, of inertia I
!      under a torque tau and a spring k, an unknown of time alone, and
!      the beam read in the hub's turning axes, a point at the radius
!      r0 + x. The action
!
!      int ( I theta_t^2 / 2 - k theta^2 / 2 + tau theta ) dt
!        + int ( rhoA ((w_t + theta_t (r0 + x))^2 + theta_t^2 w^2) / 2
!                - EI w_xx^2 / 2 + f w ) dx dt
!
!      is stationary in w where
!
!      rhoA (w_tt + theta_tt (r0 + x) - theta_t^2 w) + d_xx( EI w_xx )
!        - f = 0,
!
!      and in theta where the moment of momentum of the hub and of the
!      beam about the axis changes at the torque,
!
!      I theta_tt + k theta + int_0^L rhoA ((r0 + x) (w_tt + theta_tt
!        (r0 + x)) + theta_tt w^2 + 2 theta_t w w_t) dx - tau = 0,
!
!      the equation of an unknown of time alone: stated once per
!      instant, its integral over the beam at the quadrature points of
!      the elements. The chosen solution is case 6's deflection and
!      theta = 3 t / 10 + sin(2 t) / 4, the load and the torque the
!      laws applied to them; each block starts from the instant before
!      (continued). Measured at 801 instants over 9, 17 and 33 points,
!      the error of w falls as 1.51e-4, 9.47e-6, 5.40e-7, order 4, and
!      of theta as 4.18e-5, 2.66e-6, 5.10e-7; with theta's law reading
!      no integral the error of w is case 6's to every printed digit.
!      The same laws formed as the stationarity of the action, along t
!      and x in w and along t in theta, give the same errors, theta
!      the same to 2e-11. With the inertia I (1 + s), s = 0.2, the
!      misfit's dJ/ds is 3.79228278020e-4 by the adjoint, the same by
!      the expansion to 9e-16, and the central differences of steps
!      2e-4 and 1e-4 are within 2.0e-11 and 3.5e-12.
!
! In cases 1 and 4 the multipliers mu_x also give an estimate of the
! error of J: the representation laws evaluated on the jet that
! differences of order 8 form from the values of w, their mismatch
! weighted by mu_x,
!
!      eta = - sum over the rows of mu_x (D_8 w - w_h's jet),
!
! a relation of the open line evaluated on that jet as a whole. On the
! open line eta / (J - J_h) is 1.009 at 41 x 17, and with differences
! of order 6 1.0435, 1.0101, 1.0018, 1.0007 at 801 instants over 9 to
! 65 points. Under the fitted differences of order 6 it is 0.78 at
! 41 x 17 and 0.41, 0.35 at 801 x 17 and 33: next to each end a
! condition takes the row of a representation law, whose mismatch then
! has no multiplier (at 801 x 17 the points beside the ends contribute
! 1.7e-6 and the next ones -3.8e-4).
!
! The multipliers of the conditions are returned under their names;
! the sum over the instants of the multiplier of w = 0 at an end is
! -dJ/de for the datum e of w there at every instant (measured against
! the central difference of J on 21 x 17: 7.9479677 and 7.9479451).
!
! The instant count and the point count are read from the command
! line (default 41 and 17). Measured at 801 instants from 9 to 33
! points, the error of case 1 falls as 3.95e-3, 2.03e-4, 5.67e-6 and
! of case 2 as 6.68e-2, 8.72e-4, 2.14e-5; at 65 points from 21 to 161
! instants, case 1 as 2.36e-3, 4.82e-4, 1.06e-4, 2.47e-5 and case 2
! as 1.70e-3, 2.57e-4, 4.98e-5, 1.73e-5: order 4 or more in x and 2
! in t, case 2 reaching 9.5e-6 at 65 points.
!
! THE MINIMIZER STACKS, paired with the rows of every block as a
! multiplier is with its condition (operation_minimization), the
! minimization of each case printed before its solve: newton over
! GMRES(60) to a relative residual of 1e-3 preconditioned by four
! Gauss-Seidel sweeps, each sweep's blocks the rows' points - the whole
! line at one moment on the open line; on the hub newton to 1e-13 over
! the elimination of the rows its unknowns of time alone tie and an
! equilibrated direct solve of the rest.
!=====================================================================!

program euler_bernoulli_beam

  use util_precision             , only : dp
  use util_verbosity             , only : set_verbosity
  use operation_manifold         , only : continuous_manifold, discrete_manifold, interval, region, points
  use operation_field            , only : continuous_field, discrete_field, integral, representation
  use operation_field            , only : sin, cos, operator(+), operator(-), operator(*), operator(/), operator(**)
  use operation_lagrangian         , only : continuous_lagrangian, discrete_lagrangian, operator(+), operator(*)
  use operation_family           , only : dirk, bdf, adams, chain
  use operation_finite_difference, only : finite_difference
  use operation_finite_element   , only : finite_element
  use operation_segment_trapezoid, only : segment_trapezoid
  use operation_newton           , only : newton
  use operation_gmres            , only : gmres
  use operation_gauss_seidel     , only : gauss_seidel
  use operation_elimination      , only : elimination
  use operation_dense_direct     , only : dense_direct
  use operation_minimization     , only : minimization

  implicit none

  real(dp), parameter :: pi = acos(-1.0_dp)
  real(dp), parameter :: L = 2.0_dp, T_final = 1.0_dp, EI_0 = 2.0_dp, rhoA_0 = 1.5_dp
  integer  :: num_instants = 41, num_points = 17
  character(len=32) :: item
  type(newton) :: stack, hub_stack

  if (command_argument_count() >= 1) then
     call get_command_argument(1, item)
     read(item, *) num_instants
  end if
  if (command_argument_count() >= 2) then
     call get_command_argument(2, item)
     read(item, *) num_points
  end if

  ! each case under its name, the convergence of every solve printed
  ! by image 1
  if (this_image() == 1) call set_verbosity(1)

  ! THE MINIMIZER STACKS, paired with the rows of every block: newton
  ! over GMRES(60) to a relative residual of 1e-3, preconditioned by
  ! four Gauss-Seidel sweeps over the blocks of the rows' points; on the
  ! hub, whose unknowns of time alone tie their rows, newton to 1e-13
  ! over the elimination of the tying rows and an equilibrated direct
  ! solve of the rest
  stack     = newton(gmres(gauss_seidel(4), restart=60, tolerance=1.0e-3_dp))
  hub_stack = newton(elimination(dense_direct(equilibrated=.true.)), tolerance=1.0e-13_dp)

  call simply_supported(open_line=.false.)
  call clamped(open_line=.false.)
  call clamped(open_line=.true.)
  call simply_supported(open_line=.true.)
  call clamped(open_line=.false., damping_scale=0.3_dp)
  call clamped(open_line=.false., elements=.true.)
  call clamped(open_line=.false., damping_scale=0.3_dp, elements=.true.)
  call on_hub()

contains

  !===================================================================!
  ! Cases 1 and 4: the mode of the simply supported beam, and the
  ! estimate of the error of its energy.
  !===================================================================!

  subroutine simply_supported(open_line)

    logical, intent(in) :: open_line

    type(continuous_manifold) :: omega, initial, left, right
    type(discrete_manifold)   :: omega_h
    type(continuous_field)    :: t, x, w, w_exact, exact, energy
    type(continuous_field)    :: adjoint, mu_t, mu_x, lambda, lambda_0, lambda_L
    type(discrete_field)      :: estimate, solution, error, functional, at_0, at_L
    type(continuous_lagrangian) :: r, r_t, r_x, ic, bc_0, bc_L, Lagrangian
    type(discrete_lagrangian)   :: L_h
    type(minimization)        :: graph
    real(dp) :: frequency, J(1), J_exact
    integer  :: k

    ! Omega = [0, L] x [0, T]; its faces at t = 0, x = 0 and x = L
    omega   = continuous_manifold(time=interval(0.0_dp, T_final), space=region(interval(0.0_dp, L)))
    t       = omega % coordinate('t')
    x       = omega % coordinate('x')
    initial = omega % boundary(t, 0.0_dp)
    left    = omega % boundary(x, 0.0_dp)
    right   = omega % boundary(x, L)
    omega_h = omega % discretize(time=points(num_instants), space=points(num_points))

    frequency = (pi / L)**2 * sqrt(EI_0 / rhoA_0)
    w_exact   = sin(pi / L * x) * cos(frequency * t)
    exact     = continuous_field(omega, [w_exact])

    ! the physics, the representation laws, the initial data, and the
    ! conditions at the two ends
    w    = omega % unknown('w', [t, x])
    r    = continuous_lagrangian(omega, [rhoA_0 * w % derivative([t, t]) + EI_0 * w % derivative([x, x, x, x])])
    r_t  = continuous_lagrangian(omega, representation(w, t, 2))
    r_x  = continuous_lagrangian(omega, representation(w, x, 4))
    ic   = continuous_lagrangian(initial, [w - w_exact, w % derivative([t])])
    bc_0 = continuous_lagrangian(left , [w, w % derivative([x, x])])
    bc_L = continuous_lagrangian(right, [w, w % derivative([x, x])])

    ! the energy, and the Lagrangian
    energy   = integral((rhoA_0 * w % derivative([t])**2 + EI_0 * w % derivative([x, x])**2) / 2.0_dp, over=omega)
    adjoint  = r    % multiplier('adjoint', [t, x])
    mu_t     = r_t  % multiplier('mu_t', [t, x])
    mu_x     = r_x  % multiplier('mu_x', [t, x])
    lambda   = ic   % multiplier('lambda', [x])
    lambda_0 = bc_0 % multiplier('lambda_0', [t])
    lambda_L = bc_L % multiplier('lambda_L', [t])
    Lagrangian = energy + adjoint*r + mu_t*r_t + mu_x*r_x + lambda*ic + lambda_0*bc_0 + lambda_L*bc_L

    if (open_line) then
       L_h = Lagrangian % discretize(omega_h, &
            & time=chain([(dirk(2), k = 1, 4), (bdf(2), k = 5, 7), (adams(2), k = 8, num_instants)], &
            &            from=[(k, k = 1, num_instants)]), space=segment_trapezoid())
    else
       L_h = Lagrangian % discretize(omega_h, &
            & time=chain([(dirk(2), k = 1, 4), (bdf(2), k = 5, 7), (adams(2), k = 8, num_instants)], &
            &            from=[(k, k = 1, num_instants)]), space=finite_difference(order=6))
    end if

    if (open_line) then
       if (this_image() == 1) print '(a)', 'simply supported, uniform EI and rhoA, the open line'
    else
       if (this_image() == 1) print '(a)', 'simply supported, uniform EI and rhoA'
    end if
    estimate = exact % discretize(omega_h)
    graph = minimization(L_h, stack)
    if (this_image() == 1) call graph % view(L_h)
    call L_h % minimize(estimate, solution, by=stack)

    error      = solution % fields(['w']) - estimate
    functional = energy % at(solution)
    call functional % values(J)
    J_exact    = T_final * EI_0 * (pi / L)**4 * L / 4.0_dp
    if (this_image() == 1) print '(a, es12.4)', '   error against the exact solution  ', error % norm()
    if (this_image() == 1) print '(a, es24.16, a, es24.16)', '   J = ', J(1), '   exact ', J_exact
    if (this_image() == 1) print '(a, es12.4, a, es12.4)', '   J_exact - J = ', J_exact - J(1), &
         & '   estimate by mu_x and differences of order 8 ', &
         & sum(L_h % disagreement(solution, finite_difference(order=8)))

    ! the multipliers of w = 0 at the two ends, summed over the
    ! instants: -dJ/de for the datum e of w at that end at every instant
    at_0 = solution % fields(['lambda_0'])
    at_L = solution % fields(['lambda_L'])
    if (this_image() == 1) print '(a, 2es24.16)', '   -dJ/de at x = 0 and at x = L ', &
         & sum([(at_0 % value((k - 1) * num_points + 1, 1), k = 1, num_instants)]), &
         & sum([(at_L % value((k - 1) * num_points + 1, 1), k = 1, num_instants)])

  end subroutine simply_supported

  !===================================================================!
  ! Cases 2, 3, 5, 6 and 7: the clamped beam of varying stiffness and
  ! mass, loaded so that a chosen deflection is the solution; the
  ! derivatives along x by the fitted differences, on the open line by
  ! the trapezoid rule on its segments, or with elements by the Hermite
  ! cubics; with a damping scale, the Kelvin-Voigt damping of the
  ! bending moment.
  !===================================================================!

  subroutine clamped(open_line, damping_scale, elements)

    logical , intent(in) :: open_line
    real(dp), intent(in), optional :: damping_scale
    logical , intent(in), optional :: elements

    type(continuous_manifold) :: omega, initial, left, right
    type(discrete_manifold)   :: omega_h
    type(continuous_field)    :: t, x, w, w_exact, exact, EI, rhoA, damping, moment, moment_exact, load
    type(continuous_field)    :: lambda, lambda_0, lambda_L
    type(discrete_field)      :: estimate, solution, error
    type(continuous_lagrangian) :: r, ic, bc_0, bc_L, Lagrangian
    type(discrete_lagrangian)   :: L_h
    type(minimization)        :: graph
    real(dp) :: frequency
    integer  :: k
    logical  :: by_elements

    by_elements = .false.
    if (present(elements)) by_elements = elements
    omega   = continuous_manifold(time=interval(0.0_dp, T_final), space=region(interval(0.0_dp, L)))
    t       = omega % coordinate('t')
    x       = omega % coordinate('x')
    initial = omega % boundary(t, 0.0_dp)
    left    = omega % boundary(x, 0.0_dp)
    right   = omega % boundary(x, L)
    omega_h = omega % discretize(time=points(num_instants), space=points(num_points))

    ! the stiffness and the mass per length, fields on Omega as w is
    EI   = EI_0   * (1.0_dp + x / (2.0_dp * L))
    rhoA = rhoA_0 * (1.0_dp + cos(pi / L * x) / 4.0_dp)

    ! the chosen solution, clamped at both ends, and the bending moment
    ! of it and of the unknown, damped when a damping scale is given
    frequency    = (pi / L)**2 * sqrt(EI_0 / rhoA_0)
    w_exact      = sin(pi / L * x)**2 * cos(frequency * t)
    exact        = continuous_field(omega, [w_exact])
    w            = omega % unknown('w', [t, x])
    moment_exact = EI * w_exact % derivative([x, x])
    moment       = EI * w % derivative([x, x])
    if (present(damping_scale)) then
       damping      = damping_scale * (1.0_dp + (x / L)**2)
       moment_exact = moment_exact + damping * w_exact % derivative([x, x, t])
       moment       = moment + damping * w % derivative([x, x, t])
    end if

    ! the load the chosen solution needs, and the physics in flux form,
    ! d_xx of the bending moment
    load   = rhoA * w_exact % derivative([t, t]) + moment_exact % derivative([x, x])
    r      = continuous_lagrangian(omega, [rhoA * w % derivative([t, t]) + moment % derivative([x, x]) - load])
    ic     = continuous_lagrangian(initial, [w - w_exact, w % derivative([t])])
    bc_0   = continuous_lagrangian(left , [w, w % derivative([x])])
    bc_L   = continuous_lagrangian(right, [w, w % derivative([x])])

    lambda   = ic   % multiplier('lambda', [x])
    lambda_0 = bc_0 % multiplier('lambda_0', [t])
    lambda_L = bc_L % multiplier('lambda_L', [t])
    Lagrangian = r + lambda*ic + lambda_0*bc_0 + lambda_L*bc_L

    if (open_line) then
       L_h = Lagrangian % discretize(omega_h, &
            & time=chain([(dirk(2), k = 1, 4), (bdf(2), k = 5, 7), (adams(2), k = 8, num_instants)], &
            &            from=[(k, k = 1, num_instants)]), space=segment_trapezoid())
    else if (by_elements) then
       L_h = Lagrangian % discretize(omega_h, &
            & time=chain([(dirk(2), k = 1, 4), (bdf(2), k = 5, 7), (adams(2), k = 8, num_instants)], &
            &            from=[(k, k = 1, num_instants)]), space=finite_element(degree=3))
    else
       L_h = Lagrangian % discretize(omega_h, &
            & time=chain([(dirk(2), k = 1, 4), (bdf(2), k = 5, 7), (adams(2), k = 8, num_instants)], &
            &            from=[(k, k = 1, num_instants)]), space=finite_difference(order=6))
    end if

    if (open_line) then
       if (this_image() == 1) print '(a)', 'clamped, EI(x) and rhoA(x) varying, manufactured load, the open line'
    else if (present(damping_scale) .and. by_elements) then
       if (this_image() == 1) print '(a)', 'clamped, EI(x) and rhoA(x) varying, manufactured load, Kelvin-Voigt damping &
            &c(x) w_xxt, finite elements'
    else if (present(damping_scale)) then
       if (this_image() == 1) print '(a)', 'clamped, EI(x) and rhoA(x) varying, manufactured load, Kelvin-Voigt damping c(x) w_xxt'
    else if (by_elements) then
       if (this_image() == 1) print '(a)', 'clamped, EI(x) and rhoA(x) varying, manufactured load, finite elements'
    else
       if (this_image() == 1) print '(a)', 'clamped, EI(x) and rhoA(x) varying, manufactured load'
    end if
    estimate = exact % discretize(omega_h)
    graph = minimization(L_h, stack)
    if (this_image() == 1) call graph % view(L_h)
    call L_h % minimize(estimate, solution, by=stack)

    error = solution % fields(['w']) - estimate
    if (this_image() == 1) print '(a, es12.4)', '   error against the exact solution  ', error % norm()

  end subroutine clamped

  !===================================================================!
  ! Case 8: the clamped beam of case 6 on a hub that turns, the hub's
  ! angle an unknown of time alone whose law reads the integral of the
  ! beam's moment of momentum over the beam.
  !===================================================================!

  subroutine on_hub()

    real(dp), parameter :: r0 = 0.5_dp, inertia = 3.0_dp, spring = 4.0_dp

    type(continuous_manifold) :: omega, initial, left, right
    type(discrete_manifold)   :: omega_h
    type(continuous_field)    :: t, x, w, theta, w_exact, theta_exact, exact, EI, rhoA, moment, moment_exact
    type(continuous_field)    :: across, across_exact, turning, turning_exact, load, torque
    type(continuous_field)    :: lambda, lambda_0, lambda_L, w_only, theta_only
    type(discrete_field)      :: estimate, solution, error_w, error_theta
    type(continuous_lagrangian) :: r, ic, bc_0, bc_L, Lagrangian
    type(discrete_lagrangian)   :: L_h
    type(minimization)        :: graph
    real(dp) :: frequency
    integer  :: k

    omega   = continuous_manifold(time=interval(0.0_dp, T_final), space=region(interval(0.0_dp, L)))
    t       = omega % coordinate('t')
    x       = omega % coordinate('x')
    initial = omega % boundary(t, 0.0_dp)
    left    = omega % boundary(x, 0.0_dp)
    right   = omega % boundary(x, L)
    omega_h = omega % discretize(time=points(num_instants), space=points(num_points))

    EI   = EI_0   * (1.0_dp + x / (2.0_dp * L))
    rhoA = rhoA_0 * (1.0_dp + cos(pi / L * x) / 4.0_dp)

    ! the chosen motion: the hub's angle, and the deflection in its axes
    frequency   = (pi / L)**2 * sqrt(EI_0 / rhoA_0)
    w_exact     = sin(pi / L * x)**2 * cos(frequency * t)
    theta_exact = 0.3_dp * t + sin(2.0_dp * t) / 4.0_dp
    exact       = continuous_field(omega, [w_exact, theta_exact])

    ! w a function of (t, x); theta of t alone
    w     = omega % unknown('w')
    theta = omega % unknown('theta', [t])

    ! the acceleration across the beam, the bending moment, and the
    ! rate of the moment of momentum of a length of the beam about the
    ! axis, of the unknowns and of the chosen motion
    across        = w % derivative([t, t]) + theta % derivative([t, t]) * (r0 + x) &
         &        - theta % derivative([t])**2 * w
    across_exact  = w_exact % derivative([t, t]) + theta_exact % derivative([t, t]) * (r0 + x) &
         &        - theta_exact % derivative([t])**2 * w_exact
    moment        = EI * w % derivative([x, x])
    moment_exact  = EI * w_exact % derivative([x, x])
    turning       = rhoA * ((r0 + x) * (w % derivative([t, t]) + theta % derivative([t, t]) * (r0 + x)) &
         &        + theta % derivative([t, t]) * w**2 + 2.0_dp * theta % derivative([t]) * w * w % derivative([t]))
    turning_exact = rhoA * ((r0 + x) * (w_exact % derivative([t, t]) + theta_exact % derivative([t, t]) * (r0 + x)) &
         &        + theta_exact % derivative([t, t]) * w_exact**2 &
         &        + 2.0_dp * theta_exact % derivative([t]) * w_exact * w_exact % derivative([t]))

    ! the load and the torque the chosen motion needs, and the laws
    load   = rhoA * across_exact + moment_exact % derivative([x, x])
    torque = inertia * theta_exact % derivative([t, t]) + spring * theta_exact &
         & + integral(turning_exact, over=omega % space())
    r      = continuous_lagrangian(omega, [rhoA * across + moment % derivative([x, x]) - load, &
         &                               inertia * theta % derivative([t, t]) + spring * theta &
         &                               + integral(turning, over=omega % space()) - torque])
    ic     = continuous_lagrangian(initial, [w - w_exact, w % derivative([t]), theta - theta_exact, &
         &                                 theta % derivative([t]) - theta_exact % derivative([t])])
    bc_0   = continuous_lagrangian(left , [w, w % derivative([x])])
    bc_L   = continuous_lagrangian(right, [w, w % derivative([x])])

    lambda   = ic   % multiplier('lambda')
    lambda_0 = bc_0 % multiplier('lambda_0')
    lambda_L = bc_L % multiplier('lambda_L')
    Lagrangian = r + lambda*ic + lambda_0*bc_0 + lambda_L*bc_L

    L_h = Lagrangian % discretize(omega_h, &
         & time=chain([(dirk(2), k = 1, 4), (bdf(2), k = 5, 7), (adams(2), k = 8, num_instants)], &
         &            from=[(k, k = 1, num_instants)]), space=finite_element(degree=3))

    if (this_image() == 1) print '(a)', 'clamped, case 6 on a turning hub, the hub''s angle an unknown of time alone'
    estimate = exact % discretize(omega_h)
    graph = minimization(L_h, hub_stack)
    if (this_image() == 1) call graph % view(L_h)
    call L_h % minimize(estimate, solution, continued=.true., by=hub_stack)

    w_only      = continuous_field(omega, [w_exact])
    theta_only  = continuous_field(omega, [theta_exact])
    error_w     = solution % fields(['w']) - w_only % discretize(omega_h)
    error_theta = solution % fields(['theta']) - theta_only % discretize(omega_h)
    if (this_image() == 1) print '(a, es12.4)', '   error against the exact solution  ', error_w % norm()
    if (this_image() == 1) print '(a, es12.4)', '   error of the hub''s angle          ', error_theta % norm()

  end subroutine on_hub

end program euler_bernoulli_beam
