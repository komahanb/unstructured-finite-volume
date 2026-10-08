!=====================================================================!
! THE TRIPLE PENDULUM of the dissertation (section 5.1): three rigid
! bars, A, B and C, of masses 1, 2 and 3 and lengths 1, 2 and 3, each
! 0.1 wide and thick, A on a spherical joint at the origin, B hinged
! to A and C to B about the y axis. The bars lie along x at the start
! and turn together about the vertical at 1 rad/s under gravity.
!
! Each bar is a rigid body, its coordinates the displacement of its
! centre of mass and the Euler parameters of its rotation; each joint
! is a list of conditions paired with their reactions. The Lagrangian
! density is the bodies' kinetic energy less their potential energy,
! with the conditions times their reactions, and its stationarity in
! each unknown is that unknown's law: the equations of motion and the
! conditions of the joints, an index-three system marched by bdf(2).
! At the first instant the initial data fix the positions and the
! velocities, and each condition is stated by its second derivative
! along time, which determines the reactions there.
!
! The energy is conserved by the motion and dissipated by the
! scheme at its order; the conditions hold at every instant. The
! instant count is read from the command line (default 201).
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program triple_pendulum

  use util_precision    , only : dp
  use util_verbosity    , only : set_verbosity
  use operation_manifold, only : continuous_manifold, discrete_manifold, interval, points
  use operation_field   , only : continuous_field, discrete_field, operator(+), operator(-)
  use operation_residual, only : continuous_residual, discrete_residual, operator(+)
  use operation_family  , only : dirk, bdf, chain
  use multibody_dynamics, only : rigid_body, mechanism, spherical, revolute

  implicit none

  real(dp), parameter :: T_final = 2.0_dp, width = 0.1_dp, turning(3) = [0.0_dp, 0.0_dp, 1.0_dp]
  integer  :: num_instants = 201

  type(continuous_manifold) :: omega, domega
  type(discrete_manifold)   :: omega_h
  type(continuous_field)    :: t, energy, joints, motion
  type(discrete_field)      :: estimate, solution, energy_h, joints_h
  type(continuous_residual) :: r, g, L
  type(discrete_residual)   :: L_h
  type(mechanism)           :: pendulum
  type(rigid_body)          :: a, b, c
  real(dp), allocatable     :: e(:), defect(:,:)
  character(len=32)         :: item
  integer :: k

  if (command_argument_count() >= 1) then
     call get_command_argument(1, item)
     read(item, *) num_instants
  end if

  omega   = continuous_manifold(time=interval(0.0_dp, T_final))
  domega  = omega % boundary(omega % coordinate('t'), 0.0_dp)
  omega_h = omega % discretize(time=points(num_instants))
  t       = omega % coordinate('t')

  ! the bars, turning about the vertical through the origin
  pendulum = mechanism(t)
  call pendulum % body(omega, a, 'a', [0.5_dp, 0.0_dp, 0.0_dp], 1.0_dp, bar(1.0_dp, 1.0_dp), &
       & angular_velocity=turning)
  call pendulum % body(omega, b, 'b', [2.0_dp, 0.0_dp, 0.0_dp], 2.0_dp, bar(2.0_dp, 2.0_dp), &
       & angular_velocity=turning)
  call pendulum % body(omega, c, 'c', [4.5_dp, 0.0_dp, 0.0_dp], 3.0_dp, bar(3.0_dp, 3.0_dp), &
       & angular_velocity=turning)

  ! the joints: A on a fixed point, B hinged to A, C to B
  call pendulum % joint(omega, 'joint1', spherical(a, [0.0_dp, 0.0_dp, 0.0_dp]))
  call pendulum % joint(omega, 'joint2', revolute(a, [1.0_dp, 0.0_dp, 0.0_dp], [0.0_dp, 1.0_dp, 0.0_dp], b))
  call pendulum % joint(omega, 'joint3', revolute(b, [3.0_dp, 0.0_dp, 0.0_dp], [0.0_dp, 1.0_dp, 0.0_dp], c))

  ! the laws on [0, T], the initial data on {0}
  r = continuous_residual(omega,  pendulum % laws())
  g = continuous_residual(domega, pendulum % initial)
  L = r + g
  L_h = L % discretize(omega_h, time=chain([dirk(2), bdf(1), (bdf(2), k = 3, num_instants)], &
       &                                   from=[(k, k = 1, num_instants)]))

  ! the motion from the estimate, the convergence of every solve
  ! printed by image 1
  motion   = pendulum % estimate(omega)
  estimate = motion % discretize(omega_h)
  if (this_image() == 1) call set_verbosity(1)
  call L_h % minimize(estimate, solution, continued=.true.)

  ! the energy at every instant, and the conditions of the joints
  energy = a % kinetic_energy(t) + b % kinetic_energy(t) + c % kinetic_energy(t) &
       & + a % potential_energy(pendulum % gravity) + b % potential_energy(pendulum % gravity) &
       & + c % potential_energy(pendulum % gravity)
  energy_h = energy % at(solution)
  allocate(e(num_instants))
  call energy_h % values(e)
  joints = continuous_field(omega, [spherical(a, [0.0_dp, 0.0_dp, 0.0_dp]), &
       & revolute(a, [1.0_dp, 0.0_dp, 0.0_dp], [0.0_dp, 1.0_dp, 0.0_dp], b), &
       & revolute(b, [3.0_dp, 0.0_dp, 0.0_dp], [0.0_dp, 1.0_dp, 0.0_dp], c)])
  joints_h = joints % at(solution)
  defect = joints_h % value

  print '(a, i0, a)', 'the triple pendulum at ', num_instants, ' instants:'
  print '(a)', '      t       energy         change    largest condition'
  do k = 1, num_instants, max(1, (num_instants - 1) / 10)
     print '(f8.3, es16.8, 2es12.3)', omega_h % instant(k), e(k), e(k) - e(1), maxval(abs(defect(k, :)))
  end do

contains

  ! the inertia of a bar along x of the given mass and length about
  ! its centre of mass
  function bar(mass, length) result(inertia)
    real(dp), intent(in) :: mass, length
    real(dp) :: inertia(3,3)
    inertia = 0.0_dp
    inertia(1, 1) = mass * (width**2 + width**2) / 12.0_dp
    inertia(2, 2) = mass * (length**2 + width**2) / 12.0_dp
    inertia(3, 3) = mass * (length**2 + width**2) / 12.0_dp
  end function bar

end program triple_pendulum
