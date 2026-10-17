!=====================================================================!
! THE POINTS OF THE BEAMLINE AS A DESIGN.
!
! The clamped Euler-Bernoulli beam of
! physics/beam/euler_bernoulli_beam.f90 on the open line, uniform EI
! and rhoA, loaded so that
!
!      w = sin(pi (x / L)^2)^2 cos(omega t)
!
! is the solution, its curvature growing along the beamline. The goal
! is J = integral over Omega of w^2 / 2. The np points of the beamline
! are graded by one design coordinate s,
!
!      x_i = L (xi_i + s xi_i (1 - xi_i)),   xi_i = i / (np - 1),
!
! the ends fixed and the points ascending for |s| < 1; s > 0 moves the
! points towards x = L.
!
! The multipliers mu_x of the representation laws estimate the error
! of J, one contribution eta_ik per point and instant (section 15 of
! doc/unified-residuals.md). The objective of the grid is
!
!      Phi(s) = sum over i, k of eta_ik^2,
!
! in which contributions of opposite sign do not cancel, as they do in
! eta = sum of eta_ik; the sum of |eta_ik| has its minimum near the
! same grid but is not differentiable where a contribution changes
! sign. Newton's iteration on dPhi/ds = 0 moves s from 0, the
! derivatives of Phi by central differences of step 0.02: the gradient
! of eta in the positions would need a second adjoint, since eta reads
! the multipliers, and the derivatives of the weights of the
! concretion in the positions.
!
! Measured at 101 instants and 17 points, four iterations from s = 0 to
! s = 0.372: Phi 1.46e-5 to 1.03e-5, the error of J 2.34e-2 to -7.4e-4
! and of w 5.49e-2 to 3.14e-2. At 201 instants the same grid is found
! (0.379). Scanned at 201 instants, the estimate eta has its zero near
! s = 0.41 and the error of J near s = 0.37, and the error of w is
! least near s = 0.3. On grids graded towards x = 0 (s < 0) the
! estimate loses its accuracy: eta / (J - J_h) = 0.37 at s = -0.6.
!=====================================================================!

program euler_bernoulli_beam_grid

  use util_precision             , only : dp
  use operation_manifold         , only : continuous_manifold, discrete_manifold, interval, region, points
  use operation_field            , only : continuous_field, discrete_field, integral, representation
  use operation_field            , only : sin, cos, operator(+), operator(-), operator(*), operator(/), operator(**)
  use operation_lagrangian         , only : continuous_lagrangian, discrete_lagrangian, operator(+), operator(*)
  use operation_family           , only : dirk, bdf, adams, chain
  use operation_finite_difference, only : finite_difference
  use operation_segment_trapezoid, only : segment_trapezoid

  implicit none

  real(dp), parameter :: pi = acos(-1.0_dp)
  real(dp), parameter :: L = 2.0_dp, T_final = 1.0_dp, EI = 2.0_dp, rhoA = 1.5_dp
  real(dp), parameter :: step = 0.02_dp

  integer  :: num_instants = 101, num_points = 17
  real(dp) :: s, phi(-1:1), eta, J_error, w_error, gradient, curvature, change
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

  ! Newton's iteration on dPhi/ds = 0 from the uniform grid, the step
  ! of s bounded by 0.2 where the quadratic model is not trusted, and s
  ! within 0.9 so that the points ascend
  s = 0.0_dp
  do iteration = 1, 8
     call evaluate_grid(s - step, phi(-1))
     call evaluate_grid(s + step, phi(1))
     call evaluate_grid(s, phi(0), eta, J_error, w_error)
     if (this_image() == 1) print '(a, i0, a, f9.5, a, es11.4, a, es11.4, a, es11.4, a, es11.4)', 'iteration ', iteration, &
          & '  s = ', s, '  Phi = ', phi(0), '  eta = ', eta, '  J - J_h = ', J_error, &
          & '  |w - w_exact| = ', w_error

     gradient  = (phi(1) - phi(-1)) / (2.0_dp * step)
     curvature = (phi(1) - 2.0_dp * phi(0) + phi(-1)) / step**2
     if (curvature <= 0.0_dp) then
        error stop 'euler_bernoulli_beam_grid: Phi is not convex at s; Newton''s step does not descend'
     end if
     change = -gradient / curvature
     if (abs(change) < 1.0e-3_dp) exit
     change = max(-0.2_dp, min(0.2_dp, change))
     s = max(-0.9_dp, min(0.9_dp, s + change))
  end do

contains

  !===================================================================!
  ! The objective Phi of the grid graded by s, and optionally the
  ! estimate eta, the error of J and the error of w against the
  ! chosen solution: one forward and one reverse sweep.
  !===================================================================!

  subroutine evaluate_grid(s, phi, eta, J_error, w_error)

    real(dp), intent(in)            :: s
    real(dp), intent(out)           :: phi
    real(dp), intent(out), optional :: eta, J_error, w_error

    type(continuous_manifold) :: omega, initial, left, right
    type(discrete_manifold)   :: omega_h
    type(continuous_field)    :: t, x, w, w_exact, exact, load, goal
    type(continuous_field)    :: adjoint, mu_x, lambda, lambda_0, lambda_L
    type(discrete_field)      :: estimate, solution, functional, error
    type(continuous_lagrangian) :: r, r_x, ic, bc_0, bc_L, Lagrangian
    type(discrete_lagrangian)   :: L_h
    real(dp), allocatable :: xi(:), contribution(:)
    real(dp) :: frequency, J(1)
    integer  :: k, i

    ! Omega = [0, L] x [0, T], the beamline graded by s
    xi      = [(real(i, dp) / (num_points - 1), i = 0, num_points - 1)]
    omega   = continuous_manifold(time=interval(0.0_dp, T_final), space=region(interval(0.0_dp, L)))
    t       = omega % coordinate('t')
    x       = omega % coordinate('x')
    initial = omega % boundary(t, 0.0_dp)
    left    = omega % boundary(x, 0.0_dp)
    right   = omega % boundary(x, L)
    omega_h = omega % discretize(time=points(num_instants), space=points(L * (xi + s * xi * (1.0_dp - xi))))

    ! the chosen solution, clamped at both ends, and the load it needs
    frequency = (pi / L)**2 * sqrt(EI / rhoA)
    w_exact   = sin(pi * (x / L)**2)**2 * cos(frequency * t)
    exact     = continuous_field(omega, [w_exact])
    load      = rhoA * w_exact % derivative([t, t]) + EI * w_exact % derivative([x, x, x, x])

    ! the physics, the representation laws along x, the initial data,
    ! and the clamped ends
    w    = omega % unknown('w', [t, x])
    r    = continuous_lagrangian(omega, [rhoA * w % derivative([t, t]) + EI * w % derivative([x, x, x, x]) - load])
    r_x  = continuous_lagrangian(omega, representation(w, x, 4))
    ic   = continuous_lagrangian(initial, [w - w_exact, w % derivative([t])])
    bc_0 = continuous_lagrangian(left , [w, w % derivative([x])])
    bc_L = continuous_lagrangian(right, [w, w % derivative([x])])

    ! the goal, and the Lagrangian
    goal     = integral(w**2 / 2.0_dp, over=omega)
    adjoint  = r    % multiplier('adjoint', [t, x])
    mu_x     = r_x  % multiplier('mu_x', [t, x])
    lambda   = ic   % multiplier('lambda', [x])
    lambda_0 = bc_0 % multiplier('lambda_0', [t])
    lambda_L = bc_L % multiplier('lambda_L', [t])
    Lagrangian = goal + adjoint*r + mu_x*r_x + lambda*ic + lambda_0*bc_0 + lambda_L*bc_L

    L_h = Lagrangian % discretize(omega_h, &
         & time=chain([(dirk(2), k = 1, 4), (bdf(2), k = 5, 7), (adams(2), k = 8, num_instants)], &
         &            from=[(k, k = 1, num_instants)]), space=segment_trapezoid())

    estimate = exact % discretize(omega_h)
    call L_h % minimize(estimate, solution)

    ! the contributions to the estimate, the differences of order 6
    ! forming the second jet
    contribution = L_h % disagreement(solution, finite_difference(order=6))
    phi = sum(contribution**2)
    if (present(eta)) eta = sum(contribution)

    if (present(J_error)) then
       functional = goal % at(solution)
       call functional % values(J)
       J_error = goal_exact(frequency) - J(1)
    end if
    if (present(w_error)) then
       error   = solution % fields(['w']) - estimate
       w_error = error % norm()
    end if

  end subroutine evaluate_grid

  !===================================================================!
  ! J of the chosen solution: the integral of cos^2 over [0, T] times
  ! half the integral of sin(pi (x / L)^2)^4 over [0, L], the latter
  ! by Simpson's rule on 20000 intervals.
  !===================================================================!

  function goal_exact(frequency) result(J)

    real(dp), intent(in) :: frequency
    real(dp) :: J

    integer, parameter :: intervals = 20000
    real(dp) :: h, along_x
    integer  :: i

    h = L / intervals
    along_x = 0.0_dp
    do i = 0, intervals
       along_x = along_x + merge(1.0_dp, merge(4.0_dp, 2.0_dp, mod(i, 2) == 1), i == 0 .or. i == intervals) &
            & * sin(pi * (i * h / L)**2)**4
    end do
    along_x = along_x * h / 3.0_dp
    J = (T_final / 2.0_dp + sin(2.0_dp * frequency * T_final) / (4.0_dp * frequency)) * along_x / 2.0_dp

  end function goal_exact

end program euler_bernoulli_beam_grid
