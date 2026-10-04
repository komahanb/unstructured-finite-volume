!=====================================================================!
! THE STIFFNESS OF A BEAM IDENTIFIED BY THE STATIONARITY OF ONE
! LAGRANGIAN.
!
! The clamped Euler-Bernoulli beam of euler_bernoulli_beam.f90,
!
!      rhoA w_tt + d_xx( EI w_xx ) - f = 0     on Omega = [0, L] x [0, T],
!
! discretized along x by finite elements, the Hermite cubics of its
! case 6, and along t by its chain of families,
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
! s = 0 the gradient falls as 1.20e-1, 2.20e-2, 1.57e-3, 3.74e-5,
! 7.13e-7, 2.00e-11 (quadratic, the derivatives exact), and the design
! found approaches the true one as the grid is refined: s1 =
! 0.51178039, 0.50183579, 0.50033454 on 21 x 17, 41 x 33 and 81 x 33
! instants by points. At the design found at 21 x 17, dJ/ds by the
! expansion of the state is (-1.58e-11, -1.23e-11). At the true design
! on 21 x 17, dJ/ds by the adjoint and by the expansion agree to
! 4e-12 relative, and the central differences of step 1e-4 to 1e-6,
! their own error.
!
! A third argument bounds the design, s1 <= s1_max: the inequality is
! paired with mu >= 0,
!
!      L + mu (s1 - s1_max),
!
! and each Newton step on s is the minimizer of the quadratic model of
! J under the inequality linearized at the iterate. With s1 <= 0.45 on
! 21 x 17 the design found is (0.45, -0.22396067) in six iterations,
! mu = 2.4205141386e-5: -dJ/ds1 there by the adjoint of the design
! fixed at it, to every digit, and by the central differences of steps
! 1e-3 and 1e-4 to 1.05e-7 and 1.05e-9, their own error; dJ/ds2 =
! -4.2e-13. On 2 and 4 images mu is the same to 2.5e-8: dJ/ds is
! computed to 1e-13, as its central differences of steps 1e-9 to 1e-6
! show, and the design found moves within that, s2 by 2.2e-12.
!
! A third argument misfit=m asks instead for the least mean stiffness
! whose misfit is at most m, an inequality on a functional of the
! state paired with mu >= 0:
!
!      C(s) = integral over Omega of EI / (EI_0 L T) = 1 + s1/2 + s2/3,
!      L = C + <adjoint, r> + ... + mu (J - m),
!
! J - m one integral over Omega. C is linear in s, so the model's
! hessian is the Lagrangian's, mu times J's by the reverse sweep with
! J weighted by mu, and the first steps, before mu is known, are those
! of the shifted model from s = 0, where J exceeds m. With m = 1e-4 on
! 21 x 17 the design found is (-0.01538570, 0.40832545) in 13
! iterations, J - m = 3.0e-15, and with m = 1e-6 (0.49081178,
! -0.28174700) in 32, J - m = -1.1e-16: in both mu = -dC/ds_i / dJ/ds_i,
! 104.84257061 and 4135.9907, along s1 and along s2 alike, dJ/ds by
! the expansion of the state.
!
! A third argument scenarios asks for the stiffness from two
! deflections, at the frequency and at twice it, two scenarios each
! with its own state on one design: the sum of their misfits least,
! minimize_together over the two Lagrangians. On 21 x 17 the design
! found is (0.5612926640, -0.3552221853) in six iterations, the
! derivatives of the two misfits along the design by the expansion
! summing to -5.0e-11 and -3.9e-11 against their own 4.1e-3 and
! 2.9e-3.
!
! A fourth argument bfgs discretizes the design by the expansion of
! order zero: the gradients are those of the reverse sweeps alone, J's
! and each functional's of its own, and the hessian of the model their
! damped BFGS update, its cost apart from the number of design
! coordinates. The bound s1 <= 0.45 is reached in 19 iterations where
! the expansion took 6, the misfit at most 1e-4 in 15 and at most 1e-6
! in 45 where it took 13 and 32: the same designs, the multipliers to
! 7e-9 and the mean stiffness to 6e-13 of the expansion's.
!=====================================================================!

program euler_bernoulli_beam_optimization

  use util_precision             , only : dp
  use operation_manifold         , only : continuous_manifold, discrete_manifold, interval, region, points, parameter
  use operation_expansion        , only : expansion
  use operation_field            , only : continuous_field, discrete_field, integral
  use operation_field            , only : sin, cos, operator(+), operator(-), operator(*), operator(/), operator(**)
  use operation_residual         , only : continuous_residual, discrete_residual, operator(+), operator(*), &
       &                                  minimize_together
  use operation_family           , only : dirk, bdf, adams, chain
  use operation_finite_element   , only : finite_element

  implicit none

  real(dp), parameter :: pi = acos(-1.0_dp)
  real(dp), parameter :: L = 2.0_dp, T_final = 1.0_dp, EI_0 = 2.0_dp, rhoA_0 = 1.5_dp
  real(dp), parameter :: s_true(2) = [0.5_dp, -0.3_dp]

  integer  :: num_instants = 21, num_points = 17
  real(dp) :: s1_max, misfit_max = 0.0_dp
  logical  :: bounded = .false., scenarios = .false., quasi = .false.
  character(len=32) :: item

  type(continuous_field)    :: misfit, stiffness, misfits(2), stiffnesses(2)
  type(discrete_field)      :: estimate, solution, functional, s_h, mu_h, estimates(2), solutions(2)
  type(discrete_residual)   :: L_h, L_hs(2)
  real(dp) :: values(3), total(3)
  integer  :: k

  if (command_argument_count() >= 1) then
     call get_command_argument(1, item)
     read(item, *) num_instants
  end if
  if (command_argument_count() >= 2) then
     call get_command_argument(2, item)
     read(item, *) num_points
  end if
  if (command_argument_count() >= 3) then
     call get_command_argument(3, item)
     if (item(1:7) == 'misfit=') then
        read(item(8:), *) misfit_max
     else if (trim(item) == 'scenarios') then
        scenarios = .true.
     else
        read(item, *) s1_max
        bounded = .true.
     end if
  end if
  ! with a fourth argument bfgs, the hessian of the design's iteration
  ! by the BFGS update of the reverse sweeps' gradients
  if (command_argument_count() >= 4) then
     call get_command_argument(4, item)
     quasi = trim(item) == 'bfgs'
  end if

  ! two scenarios on one design: the deflections at the frequency and
  ! at twice it, the sum of their misfits least
  if (scenarios) then
     do k = 1, 2
        call scenario(real(k, dp), L_hs(k), estimates(k), misfits(k), stiffnesses(k))
     end do
     call minimize_together(L_hs, estimates, solutions)
     total = 0.0_dp
     do k = 1, 2
        functional = misfits(k) % at(solutions(k))
        call functional % values(values)
        total = total + values
        if (this_image() == 1) print '(a, i0, a, 3es24.16)', 'J, dJ/ds of scenario ', k, ' = ', values
     end do
     if (this_image() == 1) print '(a, 2es24.16)', 'the sum of dJ/ds          = ', total(2:3)
     s_h = solutions(1) % fields(['s1', 's2'])
     if (this_image() == 1) print '(a, 2f14.10)', 's found     = ', s_h % value(1, 1), s_h % value(1, 2)
     stop
  end if

  call scenario(1.0_dp, L_h, estimate, misfit, stiffness)
  call L_h % minimize(estimate, solution)

  ! J at the design found, dJ/ds there by the expansion of the state,
  ! and the design
  functional = misfit % at(solution)
  if (quasi) then
     call functional % values(values(1:1))
     values(2:3) = 0.0_dp
  else
     call functional % values(values)
  end if
  if (this_image() == 1) print '(a, 3es12.4)', 'J, dJ/ds    = ', values
  s_h = solution % fields(['s1', 's2'])
  if (this_image() == 1) print '(a, 2f12.8)', 's found     = ', s_h % value(1, 1), s_h % value(1, 2)
  if (this_image() == 1) print '(a, 2f12.8)', 'true s      = ', s_true
  if (command_argument_count() >= 3) then
     mu_h = solution % fields(['mu'])
     if (this_image() == 1) print '(a, es16.8)', 'mu          = ', mu_h % value(1, 1)
  end if
  ! with the misfit bounded: J and the mean stiffness, each with its
  ! derivatives by the expansion, and the stationarity of C + mu J
  if (misfit_max > 0.0_dp) then
     if (this_image() == 1) print '(a, 3es24.16)', 'J, dJ/ds    = ', values
     functional = stiffness % at(solution)
     if (quasi) then
        call functional % values(values(1:1))
     else
        call functional % values(values(1:3))
     end if
     if (this_image() == 1) print '(a, 3es24.16)', 'C, dC/ds    = ', values
  end if

contains

  !===================================================================!
  ! THE BEAM AND ITS LAGRANGIAN at the frequency times the harmonic
  ! given: its manifold, its residual discretized, the estimate, the
  ! misfit and the mean stiffness.
  !===================================================================!

  subroutine scenario(harmonic, L_h, estimate, misfit, stiffness)

    real(dp)               , intent(in)  :: harmonic
    type(discrete_residual), intent(out) :: L_h
    type(discrete_field)   , intent(out) :: estimate
    type(continuous_field) , intent(out) :: misfit, stiffness

    type(continuous_manifold) :: omega, initial, left, right
    type(discrete_manifold)   :: omega_h
    type(continuous_field)    :: t, x, s1, s2, w, w_true, exact
    type(continuous_field)    :: EI, EI_true, rhoA, moment, moment_true, load
    type(continuous_field)    :: adjoint, lambda, lambda_0, lambda_L, mu
    type(continuous_residual) :: r, ic, bc_0, bc_L, bound, Lagrangian
    real(dp) :: frequency
    integer  :: k

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
    omega_h = omega % discretize(time=points(num_instants), space=points(num_points), &
         &                       design=expansion(order=merge(0, 1, quasi)))

    ! the stiffness at the design, the true stiffness, the mass per length
    rhoA    = rhoA_0 + 0.0_dp * x
    EI      = EI_0 * (1.0_dp + s1 * x / L + s2 * (x / L)**2)
    EI_true = EI_0 * (1.0_dp + s_true(1) * x / L + s_true(2) * (x / L)**2)

    ! the deflection at the true stiffness, and the load it needs
    frequency   = harmonic * (pi / L)**2 * sqrt(EI_0 / rhoA_0)
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
    misfit    = integral((w - w_true)**2 / 2.0_dp, over=omega)
    stiffness = integral(EI / (EI_0 * L * T_final), over=omega)
    adjoint   = r    % multiplier('adjoint')
    lambda    = ic   % multiplier('lambda')
    lambda_0  = bc_0 % multiplier('lambda_0')
    lambda_L  = bc_L % multiplier('lambda_L')
    if (misfit_max > 0.0_dp) then
       ! the mean stiffness, its misfit at most the bound given: the
       ! misfit less the bound one integral over Omega
       Lagrangian = stiffness + adjoint*r + lambda*ic + lambda_0*bc_0 + lambda_L*bc_L
       bound      = continuous_residual(omega % design(), &
            & [integral((w - w_true)**2 / 2.0_dp - misfit_max / (L * T_final), over=omega)], inequality=.true.)
       mu         = bound % multiplier('mu')
       Lagrangian = Lagrangian + mu*bound
    else
       Lagrangian = misfit + adjoint*r + lambda*ic + lambda_0*bc_0 + lambda_L*bc_L
    end if

    ! the bound on s1 when given, paired with mu >= 0
    if (bounded) then
       bound      = continuous_residual(omega % design(), [s1 - s1_max], inequality=.true.)
       mu         = bound % multiplier('mu')
       Lagrangian = Lagrangian + mu*bound
    end if

    L_h = Lagrangian % discretize(omega_h, &
         & time=chain([(dirk(2), k = 1, 4), (bdf(2), k = 5, 7), (adams(2), k = 8, num_instants)], &
         &            from=[(k, k = 1, num_instants)]), space=finite_element(degree=3))
    estimate = exact % discretize(omega_h)

  end subroutine scenario

end program euler_bernoulli_beam_optimization
