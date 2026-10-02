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
! 6, which give w_xxxx at order 4. A condition at an end of the
! beamline takes the row of the physics at a point next to that end,
! so the physics holds on the points away from the ends.
!
! Five cases:
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
!      difference of step 1e-3 (its error h^2/6 J''' = -3e-7).
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
!=====================================================================!

program euler_bernoulli_beam

  use util_precision             , only : dp
  use operation_manifold         , only : continuous_manifold, discrete_manifold, interval, region, points
  use operation_field            , only : continuous_field, discrete_field, integral, representation
  use operation_field            , only : sin, cos, operator(+), operator(-), operator(*), operator(/), operator(**)
  use operation_residual         , only : continuous_residual, discrete_residual, operator(+), operator(*)
  use operation_family           , only : dirk, bdf, adams, chain
  use operation_finite_difference, only : finite_difference
  use operation_segment_trapezoid, only : segment_trapezoid

  implicit none

  real(dp), parameter :: pi = acos(-1.0_dp)
  real(dp), parameter :: L = 2.0_dp, T_final = 1.0_dp, EI_0 = 2.0_dp, rhoA_0 = 1.5_dp
  integer  :: num_instants = 41, num_points = 17
  character(len=32) :: item

  if (command_argument_count() >= 1) then
     call get_command_argument(1, item)
     read(item, *) num_instants
  end if
  if (command_argument_count() >= 2) then
     call get_command_argument(2, item)
     read(item, *) num_points
  end if

  call simply_supported(open_line=.false.)
  call clamped(open_line=.false.)
  call clamped(open_line=.true.)
  call simply_supported(open_line=.true.)
  call clamped(open_line=.false., damping_scale=0.3_dp)

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
    type(continuous_residual) :: r, r_t, r_x, ic, bc_0, bc_L, Lagrangian
    type(discrete_residual)   :: L_h
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
    r    = continuous_residual(omega, [rhoA_0 * w % derivative([t, t]) + EI_0 * w % derivative([x, x, x, x])])
    r_t  = continuous_residual(omega, representation(w, t, 2))
    r_x  = continuous_residual(omega, representation(w, x, 4))
    ic   = continuous_residual(initial, [w - w_exact, w % derivative([t])])
    bc_0 = continuous_residual(left , [w, w % derivative([x, x])])
    bc_L = continuous_residual(right, [w, w % derivative([x, x])])

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

    estimate = exact % discretize(omega_h)
    call L_h % minimize(estimate, solution)

    error      = solution % fields(['w']) - estimate
    functional = energy % at(solution)
    call functional % values(J)
    J_exact    = T_final * EI_0 * (pi / L)**4 * L / 4.0_dp
    if (open_line) then
       if (this_image() == 1) print '(a)', 'simply supported, uniform EI and rhoA, the open line'
    else
       if (this_image() == 1) print '(a)', 'simply supported, uniform EI and rhoA'
    end if
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
  ! Cases 2, 3 and 5: the clamped beam of varying stiffness and mass,
  ! loaded so that a chosen deflection is the solution; the derivatives
  ! along x by the fitted differences, or on the open line by the
  ! trapezoid rule on its segments; with a damping scale, the
  ! Kelvin-Voigt damping of the bending moment.
  !===================================================================!

  subroutine clamped(open_line, damping_scale)

    logical , intent(in) :: open_line
    real(dp), intent(in), optional :: damping_scale

    type(continuous_manifold) :: omega, initial, left, right
    type(discrete_manifold)   :: omega_h
    type(continuous_field)    :: t, x, w, w_exact, exact, EI, rhoA, damping, moment, moment_exact, load
    type(continuous_field)    :: lambda, lambda_0, lambda_L
    type(discrete_field)      :: estimate, solution, error
    type(continuous_residual) :: r, ic, bc_0, bc_L, Lagrangian
    type(discrete_residual)   :: L_h
    real(dp) :: frequency
    integer  :: k

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
    r      = continuous_residual(omega, [rhoA * w % derivative([t, t]) + moment % derivative([x, x]) - load])
    ic     = continuous_residual(initial, [w - w_exact, w % derivative([t])])
    bc_0   = continuous_residual(left , [w, w % derivative([x])])
    bc_L   = continuous_residual(right, [w, w % derivative([x])])

    lambda   = ic   % multiplier('lambda', [x])
    lambda_0 = bc_0 % multiplier('lambda_0', [t])
    lambda_L = bc_L % multiplier('lambda_L', [t])
    Lagrangian = r + lambda*ic + lambda_0*bc_0 + lambda_L*bc_L

    if (open_line) then
       L_h = Lagrangian % discretize(omega_h, &
            & time=chain([(dirk(2), k = 1, 4), (bdf(2), k = 5, 7), (adams(2), k = 8, num_instants)], &
            &            from=[(k, k = 1, num_instants)]), space=segment_trapezoid())
    else
       L_h = Lagrangian % discretize(omega_h, &
            & time=chain([(dirk(2), k = 1, 4), (bdf(2), k = 5, 7), (adams(2), k = 8, num_instants)], &
            &            from=[(k, k = 1, num_instants)]), space=finite_difference(order=6))
    end if

    estimate = exact % discretize(omega_h)
    call L_h % minimize(estimate, solution)

    error = solution % fields(['w']) - estimate
    if (open_line) then
       if (this_image() == 1) print '(a)', 'clamped, EI(x) and rhoA(x) varying, manufactured load, the open line'
    else if (present(damping_scale)) then
       if (this_image() == 1) print '(a)', 'clamped, EI(x) and rhoA(x) varying, manufactured load, Kelvin-Voigt damping c(x) w_xxt'
    else
       if (this_image() == 1) print '(a)', 'clamped, EI(x) and rhoA(x) varying, manufactured load'
    end if
    if (this_image() == 1) print '(a, es12.4)', '   error against the exact solution  ', error % norm()

  end subroutine clamped

end program euler_bernoulli_beam
