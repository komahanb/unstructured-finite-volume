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
! The Lagrangian pairs every law with its multiplier, the condition
! fixing the design with kappa:
!
!      L = J + <adjoint, r> + <lambda, ic> + <lambda_0, bc_0>
!            + <lambda_L, bc_L> + <kappa, s - s_given>
!
! Its stationarity in the multipliers is the beam at the given design,
! in the state the adjoint, and in the design kappa = -dJ/ds; the
! expansion of order one along each design coordinate gives the
! derivatives of kappa, -d^2J/ds^2. The stationarity dL/ds = 0 of the
! design is solved by Newton's iteration on s with that gradient and
! that hessian, from s = 0, one forward and one reverse sweep per
! iteration. The objective is a choice of this example: an
! identification, which has an optimum at the true design up to the
! error of the discretization.
!
! Measured: from s = 0 the gradient falls as 1.19e-1, 2.19e-2, 1.55e-3,
! 2.66e-5, 4.58e-7, 8.32e-12 (quadratic, the derivatives exact), and
! the design found approaches the true one as the grid is refined:
! s1 = 0.50961829, 0.50205024, 0.50054911 on 21 x 17, 41 x 33 and
! 81 x 33 instants by points.
!=====================================================================!

program euler_bernoulli_beam_optimization

  use util_precision             , only : dp
  use operation_manifold         , only : continuous_manifold, discrete_manifold, interval, region, points
  use operation_manifold         , only : parameter, expansion
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
  real(dp) :: s(2), J, gradient(2), hessian(2,2)
  integer  :: iteration
  character(len=32) :: item

  if (command_argument_count() >= 1) then
     call get_command_argument(1, item)
     read(item, *) num_instants
  end if
  if (command_argument_count() >= 2) then
     call get_command_argument(2, item)
     read(item, *) num_points
  end if

  ! Newton's iteration on the stationarity of L in the design, from
  ! s = 0, until the gradient vanishes to the arithmetic's precision
  s = 0.0_dp
  do iteration = 1, 8
     call evaluate_misfit(s, J, gradient, hessian)
     if (this_image() == 1) print '(a, i0, a, 2f12.8, a, es12.4, a, es12.4)', 'iteration ', iteration, '  s = ', s, &
          & '  J = ', J, '  |dJ/ds| = ', norm2(gradient)
     if (norm2(gradient) < 1.0e-10_dp) exit
     s = s - linear_solution(hessian, gradient)
  end do
  if (this_image() == 1) print '(a, 2f12.8)', 'true s      = ', s_true

contains

  !===================================================================!
  ! The solution x of A x = b for a matrix of order two.
  !===================================================================!

  function linear_solution(A, b) result(x)

    real(dp), intent(in) :: A(2,2), b(2)
    real(dp) :: x(2)

    real(dp) :: determinant

    determinant = A(1,1) * A(2,2) - A(1,2) * A(2,1)
    x = [A(2,2) * b(1) - A(1,2) * b(2), A(1,1) * b(2) - A(2,1) * b(1)] / determinant

  end function linear_solution

  !===================================================================!
  ! The misfit J at the design s, its gradient -kappa, and its hessian
  ! -dkappa/ds: one forward and one reverse sweep of the Lagrangian,
  ! the design expanded to order one along each coordinate.
  !===================================================================!

  subroutine evaluate_misfit(s, J, gradient, hessian)

    real(dp), intent(in)  :: s(2)
    real(dp), intent(out) :: J, gradient(2), hessian(2,2)

    type(continuous_manifold) :: omega, initial, left, right
    type(discrete_manifold)   :: omega_h
    type(continuous_field)    :: t, x, s1, s2, w, w_true, exact
    type(continuous_field)    :: EI, EI_true, rhoA, moment, moment_true, load, misfit
    type(continuous_field)    :: adjoint, lambda, lambda_0, lambda_L, kappa
    type(discrete_field)      :: estimate, solution, functional, kappa_h, kappa_1, kappa_2
    type(continuous_residual) :: r, ic, bc_0, bc_L, d, Lagrangian
    type(discrete_residual)   :: L_h
    real(dp) :: frequency, values(3)
    integer  :: k

    ! Omega = [0, L] x [0, T] at the design point (s1, s2); its faces at
    ! t = 0, x = 0 and x = L
    omega   = continuous_manifold(time=interval(0.0_dp, T_final), space=region(interval(0.0_dp, L)), &
         &                        design=[parameter('s1', s(1)), parameter('s2', s(2))])
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

    ! the physics in flux form, the initial data, the clamped ends, and
    ! the condition fixing the design
    w      = omega % unknown('w', [t, x, s1, s2])
    moment = EI * w % derivative([x, x])
    r      = continuous_residual(omega, [rhoA * w % derivative([t, t]) + moment % derivative([x, x]) - load])
    ic     = continuous_residual(initial, [w - w_true, w % derivative([t])])
    bc_0   = continuous_residual(left , [w, w % derivative([x])])
    bc_L   = continuous_residual(right, [w, w % derivative([x])])
    d      = continuous_residual(omega % design(), [s1 - s(1), s2 - s(2)])

    ! the misfit, and the Lagrangian
    misfit   = integral((w - w_true)**2 / 2.0_dp, over=omega)
    adjoint  = r    % multiplier('adjoint', [t, x, s1, s2])
    lambda   = ic   % multiplier('lambda', [x, s1, s2])
    lambda_0 = bc_0 % multiplier('lambda_0', [t, s1, s2])
    lambda_L = bc_L % multiplier('lambda_L', [t, s1, s2])
    kappa    = d    % multiplier('kappa')
    Lagrangian = misfit + adjoint*r + lambda*ic + lambda_0*bc_0 + lambda_L*bc_L + kappa*d

    L_h = Lagrangian % discretize(omega_h, &
         & time=chain([(dirk(2), k = 1, 4), (bdf(2), k = 5, 7), (adams(2), k = 8, num_instants)], &
         &            from=[(k, k = 1, num_instants)]), space=finite_difference(order=6))

    estimate = exact % discretize(omega_h)
    call L_h % minimize(estimate, solution)

    ! J, its gradient -kappa, and its hessian from the derivatives of
    ! kappa along each design coordinate
    functional = misfit % at(solution)
    call functional % values(values)
    J = values(1)

    kappa_h  = solution % fields(['kappa'])
    kappa_1  = kappa_h % derivative([s1])
    kappa_2  = kappa_h % derivative([s2])
    gradient = -[kappa_h % value(1, 1), kappa_h % value(1, 2)]
    hessian(:, 1) = -[kappa_1 % value(1, 1), kappa_1 % value(1, 2)]
    hessian(:, 2) = -[kappa_2 % value(1, 1), kappa_2 % value(1, 2)]

  end subroutine evaluate_misfit

end program euler_bernoulli_beam_optimization
