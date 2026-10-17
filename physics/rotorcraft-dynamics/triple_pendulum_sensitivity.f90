!=====================================================================!
! THE SENSITIVITIES OF THE TRIPLE PENDULUM, by the adjoint and by the
! expansion along the design: the pendulum of triple_pendulum.f90 over
! one second, 101 instants, with two design coordinates. s1 scales the
! mass of A's translation, (s1 - 1) times A's kinetic energy of
! translation added to the action; s2 moves the fixed point of joint 1
! along z by (s2 - 1) 0.01 (1 - cos 2t), a condition stated at the
! first instant by its second derivative along time, where s2 enters.
! The objective is J, the integral over time of the square of C's fall,
! and the design conditions s - s_0 pair with kappa = -dJ/ds.
!
! At s = (1, 1), J = 3.8528427173363906: dJ/ds by the adjoint
! -2.6054981030793269e-2 and -8.7722777922432382e-3, by the expansion
! the same to 2.9e-13 and 2.3e-14, 7.6e-15 and 2e-16, where J itself
! varies by 6e-15 between designs 1e-9 apart; by central differences
! of step 1e-4 the same to 1.5e-9 and 8.3e-10, their own error. The
! same digits on 1, 2 and 4 images.
! The design values and the order of the expansion are read from the
! command line (default 1, 1, 1).
!
! THE MINIMIZER STACK, paired with the rows of every block as a
! multiplier is with its condition (operation_minimization), the
! minimization printed before the march: newton to 1e-13 by the chord
! method, every step of a block against the factors of its first,
! over the elimination of the rows the family's relations, the data and
! the joints' unknowns tie and a dense direct solve of the rest; the
! expansion and the adjoint of each block against the factors of its
! jacobian at its solution.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program triple_pendulum_sensitivity

  use util_precision     , only : dp
  use util_verbosity     , only : set_verbosity
  use operation_manifold , only : continuous_manifold, discrete_manifold, interval, points, parameter
  use operation_expansion, only : expansion
  use operation_field    , only : continuous_field, discrete_field, integral, cos, &
       &                          operator(+), operator(-), operator(*), operator(/), operator(**)
  use operation_residual , only : continuous_residual, discrete_residual, operator(+), operator(*)
  use operation_family   , only : dirk, bdf, chain
  use multibody_dynamics , only : rigid_body, mechanism, joint_conditions, spherical, revolute
  use operation_newton      , only : newton
  use operation_elimination , only : elimination
  use operation_dense_direct, only : dense_direct
  use operation_minimization, only : minimization

  implicit none

  real(dp), parameter :: T_final = 1.0_dp, width = 0.1_dp, turning(3) = [0.0_dp, 0.0_dp, 1.0_dp]
  integer  :: num_instants = 101, order = 1
  real(dp) :: s1_value = 1.0_dp, s2_value = 1.0_dp

  type(continuous_manifold) :: omega, domega
  type(discrete_manifold)   :: omega_h
  type(continuous_field)    :: t, s1, s2, motion, objective, lambda_r, lambda, kappa
  type(joint_conditions)    :: g1
  type(discrete_field)      :: estimate, solution, kappa_h, functional
  type(continuous_residual) :: r, g, d, L
  type(discrete_residual)   :: L_h
  type(newton)              :: stack
  type(minimization)        :: graph
  type(mechanism)           :: pendulum
  type(rigid_body)          :: a, b, c
  real(dp), allocatable     :: values(:)
  character(len=32)         :: item
  integer :: k

  if (command_argument_count() >= 1) then
     call get_command_argument(1, item); read(item, *) s1_value
  end if
  if (command_argument_count() >= 2) then
     call get_command_argument(2, item); read(item, *) s2_value
  end if
  if (command_argument_count() >= 3) then
     call get_command_argument(3, item); read(item, *) order
  end if

  omega   = continuous_manifold(time=interval(0.0_dp, T_final), &
       &                        design=[parameter('s1', s1_value), parameter('s2', s2_value)])
  domega  = omega % boundary(omega % coordinate('t'), 0.0_dp)
  omega_h = omega % discretize(time=points(num_instants), design=expansion(order=order))
  t  = omega % coordinate('t')
  s1 = omega % coordinate('s1')
  s2 = omega % coordinate('s2')

  pendulum = mechanism(t)
  call pendulum % body(omega, a, 'a', [0.5_dp, 0.0_dp, 0.0_dp], 1.0_dp, bar(1.0_dp, 1.0_dp), angular_velocity=turning)
  call pendulum % body(omega, b, 'b', [2.0_dp, 0.0_dp, 0.0_dp], 2.0_dp, bar(2.0_dp, 2.0_dp), angular_velocity=turning)
  call pendulum % body(omega, c, 'c', [4.5_dp, 0.0_dp, 0.0_dp], 3.0_dp, bar(3.0_dp, 3.0_dp), angular_velocity=turning)

  ! the design s1 in the mass of A: (s1 - 1) times A's kinetic energy of translation added
  pendulum % density = pendulum % density + (s1 - 1.0_dp) * (1.0_dp / 2.0_dp) * (a % u(1) % derivative([t])**2 &
       & + a % u(2) % derivative([t])**2 + a % u(3) % derivative([t])**2)

  ! the design s2 in joint 1: its point moved along z by (s2 - 1) 0.01 (1 - cos 2t),
  ! whose rate at the first instant is zero
  g1 = spherical(a, [0.0_dp, 0.0_dp, 0.0_dp])
  g1 % field(3) = g1 % field(3) - (s2 - 1.0_dp) * 0.01_dp * (1.0_dp - cos(2.0_dp * t))
  call pendulum % joint(omega, 'joint1', g1)
  call pendulum % joint(omega, 'joint2', revolute(a, [1.0_dp, 0.0_dp, 0.0_dp], [0.0_dp, 1.0_dp, 0.0_dp], b))
  call pendulum % joint(omega, 'joint3', revolute(b, [3.0_dp, 0.0_dp, 0.0_dp], [0.0_dp, 1.0_dp, 0.0_dp], c))

  ! J, the integral of the square of C's fall
  objective = integral(c % u(3) * c % u(3), over=omega % time())

  r = continuous_residual(omega,  pendulum % laws())
  g = continuous_residual(domega, pendulum % initial)
  d = continuous_residual(omega % design(), [s1 - s1_value, s2 - s2_value])
  lambda_r = r % multiplier('adjoint')
  lambda   = g % multiplier('lambda')
  kappa    = d % multiplier('kappa')
  L = objective + lambda_r * r + lambda * g + kappa * d
  L_h = L % discretize(omega_h, time=chain([dirk(2), bdf(1), (bdf(2), k = 3, num_instants)], &
       &                                   from=[(k, k = 1, num_instants)]))

  ! the motion from the estimate, the convergence of every solve
  ! printed by image 1
  motion   = pendulum % estimate(omega)
  estimate = motion % discretize(omega_h)
  if (this_image() == 1) call set_verbosity(1)
  stack = newton(elimination(dense_direct()), tolerance=1.0e-13_dp, chord=.true.)
  graph = minimization(L_h, stack)
  if (this_image() == 1) call graph % view(L_h)
  call L_h % minimize(estimate, solution, continued=.true., by=stack)

  kappa_h = solution % fields(['kappa'])
  functional = objective % at(solution)
  allocate(values(1 + 2 * order))
  call functional % values(values)
  print '(a, es24.16)', 'J                        = ', values(1)
  print '(a, 2es24.16)', 'dJ/ds by the adjoint     = ', -kappa_h % value(1, 1), -kappa_h % value(1, 2)
  if (order >= 1) print '(a, 2es24.16)', 'dJ/ds by the expansion   = ', values(2), values(2 + order)

contains

  function bar(mass, length) result(inertia)
    real(dp), intent(in) :: mass, length
    real(dp) :: inertia(3,3)
    inertia = 0.0_dp
    inertia(1, 1) = mass * (width**2 + width**2) / 12.0_dp
    inertia(2, 2) = mass * (length**2 + width**2) / 12.0_dp
    inertia(3, 3) = mass * (length**2 + width**2) / 12.0_dp
  end function bar

end program triple_pendulum_sensitivity
