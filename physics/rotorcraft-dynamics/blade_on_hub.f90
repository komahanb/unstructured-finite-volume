!=====================================================================!
! A BLADE ON A TURNING HUB: the blade of the rotor, of span 1.56 from
! the radius 0.44, chord 0.121, thickness 0.02 and density 2500, fixed
! to a hub that its driver turns at 109.12 rad/s about z, over a tenth
! of a revolution by bdf(1), under gravity.
!
! As a beam (multibody_dynamics, beam) the blade is a field of time and
! of its span in the hub's axes: the displacement along its span,
! chord and thickness and its twist, by the hermite cubics on the span,
! the hub's coordinates and the driver's reactions unknowns of time
! alone, whose laws read the integral over the span of the blade's
! density. As a rigid plate it is a body of the plate's mass and
! inertia linked to the hub. Both begin at rest in the hub's axes.
!
!      blade_on_hub beam  [steps] [span points] [E]
!      blade_on_hub rest  [steps] [span points] [E]
!      blade_on_hub rigid [steps]
!
! (default 36 steps, 5 span points, E = 70e9). Printed: the reactions
! of the driver's three rotational conditions at the last instant, and
! for the beam the displacements and twist at the tip. At rest the hub
! is still and the initial data hold the rates and the accelerations,
! zero (mechanism % initial([1, 2])): the configuration at the first
! instant is solved, the blade's droop under its weight, and the tip's
! displacements there are printed beside the droop of a clamped beam,
! q L^4 / (8 EI).
!
! Measured: as E grows the beam's reactions approach the plate's. The
! two about x and y, the moment of the blade's weight, are -225.9008
! and -225.91163 at E = 7e15 and 7e17 against the plate's -225.91174,
! the remainder falling as 1/E: the axial force of the turning on the
! sag. The one about z falls with the step for both, 3531.5, 1764.7,
! 882.1 at 36, 72 and 144 steps for the beam at E = 7e15, 6446.2,
! 3223.3, 1611.7 for the plate: the dissipation of bdf(1), which the
! driver makes up, of the first order, each formulation its own. At
! E = 70e9 the tip has moved by 7.55e-4 along the span, -4.15e-3 along
! the chord and -1.66e-4 along the thickness. At rest the droop at the
! first instant is the clamped beam's to every printed digit where the
! beam is linear, -7.78109356e-8 at E = 7e15 by the hermite cubics,
! exact for a uniform load; at E = 70e9 the axial strain of the sag
! stiffens it, -7.77976430e-3, -7.78105888e-3 and -7.78109283e-3 at 3,
! 5 and 9 points against -7.78109356e-3, the tip drawn in along the
! span by 2.2e-5; the march then stays at the droop.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program blade_on_hub

  use util_precision          , only : dp
  use operation_manifold      , only : continuous_manifold, discrete_manifold, interval, region, points
  use operation_field         , only : continuous_field, discrete_field
  use operation_residual      , only : continuous_residual, discrete_residual
  use operation_family        , only : dirk, bdf, chain
  use operation_finite_element, only : finite_element
  use multibody_dynamics      , only : rigid_body, beam_blade, solid_rectangle, mechanism, revolute_driver, rigid_link

  implicit none

  real(dp), parameter :: pi = acos(-1.0_dp), speed = 109.12_dp, ez(3) = [0.0_dp, 0.0_dp, 1.0_dp]
  real(dp), parameter :: density = 2500.0_dp, root = 0.44_dp, tip = 2.0_dp, chord = 0.121_dp, height = 0.9_dp
  real(dp), parameter :: thickness = 0.02_dp, poisson = 0.3_dp
  real(dp), parameter :: axes(3,3) = reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
       &                                      1.0_dp], [3, 3])

  character(len=16) :: model = 'beam'
  character(len=32) :: item
  real(dp) :: modulus = 70.0e9_dp, inertia(3,3), plate(3,3), mass, rate
  integer  :: steps = 36, span_points = 5, k, num_instants

  type(continuous_manifold) :: omega, initial, left
  type(discrete_manifold)   :: omega_h
  type(continuous_field)    :: t, x, motion
  type(discrete_field)      :: estimate, solution, reactions, tip_h
  type(continuous_residual) :: L
  type(discrete_residual)   :: L_h
  type(mechanism)           :: spinner
  type(rigid_body)          :: hub, slab
  type(beam_blade)          :: blade
  character(len=32), allocatable :: names(:)

  if (command_argument_count() >= 1) call get_command_argument(1, model)
  if (command_argument_count() >= 2) then
     call get_command_argument(2, item)
     read(item, *) steps
  end if
  if (command_argument_count() >= 3) then
     call get_command_argument(3, item)
     read(item, *) span_points
  end if
  if (command_argument_count() >= 4) then
     call get_command_argument(4, item)
     read(item, *) modulus
  end if
  num_instants = steps + 1

  ! a tenth of a revolution; for the beam, of time and of the span
  select case (trim(model))
  case ('beam', 'rest')
     omega = continuous_manifold(time=interval(0.0_dp, 0.2_dp * pi / speed), space=region(interval(0.0_dp, tip - root)))
     x     = omega % coordinate('x')
     left  = omega % boundary(x, 0.0_dp)
  case ('rigid')
     omega = continuous_manifold(time=interval(0.0_dp, 0.2_dp * pi / speed))
  case default
     error stop 'blade_on_hub: the blade is a beam, a beam at rest or a rigid plate'
  end select
  t       = omega % coordinate('t')
  initial = omega % boundary(t, 0.0_dp)

  ! the hub, turning at its rate about z, and the blade
  inertia = 0.0_dp
  inertia(1, 1) = 0.5_dp
  inertia(2, 2) = 0.5_dp
  inertia(3, 3) = 1.0_dp
  ! at rest the hub still, its rate zero
  rate = speed
  if (trim(model) == 'rest') rate = 0.0_dp
  spinner = mechanism(t)
  call spinner % body(omega, hub, 'hub', [0.0_dp, 0.0_dp, height], 20.0_dp, inertia, angular_velocity=rate * ez)
  if (trim(model) /= 'rigid') then
     call spinner % beam(omega, blade, hub, 'blade', x, [root, 0.0_dp, height], axes, &
          & solid_rectangle(density, modulus, poisson, chord, continuous_field(omega, [thickness])))
  else
     mass  = density * thickness * (tip - root) * chord
     plate = 0.0_dp
     plate(1, 1) = mass * (chord**2 + thickness**2) / 12.0_dp
     plate(2, 2) = mass * ((tip - root)**2 + thickness**2) / 12.0_dp
     plate(3, 3) = mass * ((tip - root)**2 + chord**2) / 12.0_dp
     call spinner % body(omega, slab, 'blade', [(root + tip) / 2.0_dp, 0.0_dp, height], mass, plate, &
          & angular_velocity=speed * ez)
     call spinner % joint(omega, 'link', rigid_link(hub, slab))
  end if
  call spinner % joint(omega, 'driver', revolute_driver(hub, ez, rate, t))

  ! the laws, the initial data - at rest the rates and the
  ! accelerations, the configuration then solved - and for the beam
  ! its root fixed to the hub
  if (trim(model) == 'rest') then
     L = continuous_residual(omega, spinner % laws()) + continuous_residual(initial, spinner % initial([1, 2]))
  else
     L = continuous_residual(omega, spinner % laws()) + continuous_residual(initial, spinner % initial([0, 1]))
  end if
  if (trim(model) /= 'rigid') then
     L = L + continuous_residual(left, spinner % root)
     omega_h = omega % discretize(time=points(num_instants), space=points(span_points))
     L_h = L % discretize(omega_h, time=chain([dirk(2), (bdf(1), k = 2, num_instants)], &
          &                                   from=[(k, k = 1, num_instants)]), space=finite_element(degree=3))
  else
     omega_h = omega % discretize(time=points(num_instants))
     L_h = L % discretize(omega_h, time=chain([dirk(2), (bdf(1), k = 2, num_instants)], &
          &                                   from=[(k, k = 1, num_instants)]))
  end if

  motion   = spinner % estimate(omega)
  estimate = motion % discretize(omega_h)
  call L_h % minimize(estimate, solution, continued=.true.)

  names = [character(len=32) :: 'driver_4', 'driver_5', 'driver_6']
  reactions = solution % fields(names)
  print '(a, a, a, i0, a)', 'the blade as a ', trim(model), ', ', steps, ' steps:'
  print '(a, 3es24.16)', '   the driver''s rotational reactions at the last instant ', &
       & reactions % value(size(reactions % value, 1), :)
  if (trim(model) /= 'rigid') then
     names = [character(len=32) :: 'blade_u', 'blade_v', 'blade_w', 'blade_phi']
     tip_h = solution % fields(names)
     print '(a, 4es16.8)', '   the tip''s displacements and twist ', tip_h % value(size(tip_h % value, 1), :)
  end if
  ! at rest, the tip's displacements at the first instant, and the
  ! droop of a clamped beam under its weight, q L^4 / (8 EI)
  if (trim(model) == 'rest') then
     print '(a, 4es16.8)', '   at the first instant                ', tip_h % value(span_points, :)
     print '(a, es16.8)', '   q L^4 / (8 EI)                      ', -density * chord * thickness * 9.81_dp &
          & * (tip - root)**4 / (8.0_dp * modulus * chord * thickness**3 / 12.0_dp)
  end if

end program blade_on_hub
