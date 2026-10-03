!=====================================================================!
! THE DYNAMICS OF A FOUR-BLADED ROTOR (dissertation, section 5.4): the
! case of komahanb/rotorcraft-optimization, its 26 bodies and 34
! joints as collective.py, longitudinal_cyclic.py and
! lateral_cyclic.py state them.
!
! The hub turns at 109.12 rad/s about z. Three push rods move the
! lower swashplate along z by a sin(w t), w a quarter of the hub's
! rate; the upper swashplate turns with the hub on a hinge about z,
! and through the lower and upper pitch links and the blade caps it
! sets the pitch of the blades. The base plate is fixed, the sphere
! slides along the shaft's axis, and two horns hinge the lower
! swashplate to the base plate.
!
! The bodies are rigid, their table in rotor_bodies.txt. The blades
! are plates of span 1.56 from the radius 0.44, chord 0.121 and the
! thickness of the setup, density 2500, E = 70e9 and nu = 0.3; each
! lies along the azimuth of its cap, 0, 90, 180 and 270 degrees. The
! blade model is chosen at setup: rigid, each blade a rigid body of
! the plate's mass and inertia linked to its cap; beam, each blade a
! beam in the axes of its cap (multibody_dynamics, beam), fixed to it
! at the root, a field of time and of the span by the hermite cubics,
! the bodies' coordinates and the reactions unknowns of time alone;
! plate and solid are not yet stated.
!
! The Lagrangian density is the bodies' kinetic energy less their
! potential energy, with every condition of the joints times its
! reaction; its stationarity in each unknown is that unknown's law.
! The motion is marched by bdf(1) over the revolution, each instant
! from the solution at the one before. The energy at the first instant
! is 4.99144494e5 J, and over the revolution of 360 steps it falls to
! 4.94255e5, the dissipation of bdf(1); the hub's angle holds to
! 1e-16 and every condition of the joints to 2e-14 at every instant.
! With the blades as beams the march is solved with the rows that tie
! the copies of the unknowns of time alone eliminated, the rest by a
! direct solve; at the tenth step of 360 the tip of the blade at zero
! degrees has moved by 7.22238e-4 along the span, and along the chord
! by -7.10262e-4, -7.07935e-4, -7.07820e-4, -7.07813e-4 at 3, 5, 9
! and 17 points of the span, order four.
!
! THE SETUP, rotorcraft.nml:
!
!      blade        rigid or beam (plate and solid to come)
!      scenario     collective, longitudinal or lateral: the push
!                   rods' amplitudes, 0.05 each, or 0.025, 0.025 and
!                   0.05, or 0.025, 0.05 and 0.025 at 90, 180 and 270
!                   degrees
!      thickness    of the blades
!      revolutions  of the hub
!      steps        per revolution
!      meshes       the directory of the bodies' surface meshes, the
!                   case's data/mesh; none, no drawing
!      every        the stride of the instants drawn
!      span_points  of a blade as a beam
!
! THE DRAWING: every rigid body is drawn with the surface of its mesh,
! each blade with the box of its plate, moved with the body to each
! instant drawn; a blade as a beam with each point of the box moved
! with its section of the beam, and the beam's deflection a vector at
! the point: rotor_<instant>.vtu and rotor.pvd for paraview.
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program rotorcraft_dynamics

  use util_precision    , only : dp
  use util_verbosity    , only : set_verbosity
  use operation_manifold, only : continuous_manifold, discrete_manifold, interval, region, points
  use operation_field   , only : continuous_field, discrete_field
  use operation_residual, only : continuous_residual, discrete_residual, operator(+)
  use operation_family  , only : dirk, bdf, chain
  use operation_finite_element, only : finite_element
  use multibody_dynamics, only : rigid_body, beam_blade, mechanism, body_table, surface, spherical, revolute, &
       &                         cylindrical, fixed, rigid_link, revolute_driver, motion_driver
  use view_bdf_loader   , only : bdf_surface

  implicit none

  real(dp), parameter :: pi = acos(-1.0_dp), speed = 109.12_dp
  real(dp), parameter :: ex(3) = [1.0_dp, 0.0_dp, 0.0_dp], ey(3) = [0.0_dp, 1.0_dp, 0.0_dp], &
       &                 ez(3) = [0.0_dp, 0.0_dp, 1.0_dp], origin(3) = 0.0_dp
  real(dp), parameter :: density = 2500.0_dp, root = 0.44_dp, tip = 2.0_dp, chord = 0.121_dp, height = 0.9_dp
  real(dp), parameter :: modulus = 70.0e9_dp, poisson = 0.3_dp

  ! the setup
  character(len=16)  :: blade = 'rigid', scenario = 'collective'
  character(len=256) :: meshes = ''
  real(dp)           :: thickness = 0.02_dp
  integer            :: revolutions = 1, steps = 360, every = 10, span_points = 9
  namelist /setup/ blade, scenario, thickness, revolutions, steps, meshes, every, span_points

  type(continuous_manifold) :: omega, domega, left
  type(discrete_manifold)   :: omega_h
  type(continuous_field)    :: t, x, motion, energy
  type(discrete_field)      :: estimate, solution, energy_h, hub_h, tip_h
  type(continuous_residual) :: r, g, L
  type(discrete_residual)   :: L_h
  type(mechanism)           :: rotor
  type(body_table)          :: table
  type(rigid_body), allocatable :: part(:)
  type(rigid_body)          :: blades(4)
  type(beam_blade)          :: beams(4)
  real(dp)                  :: amplitude(3)
  real(dp), allocatable     :: e(:), hub_eps(:)
  character(len=32), allocatable :: names(:)
  integer :: unit, k, n, num_instants, lower, upper, cap, rod, np

  open(newunit=unit, file='rotorcraft.nml', status='old', action='read')
  read(unit, nml=setup)
  close(unit)
  select case (trim(scenario))
  case ('collective');   amplitude = [0.05_dp, 0.05_dp, 0.05_dp]
  case ('longitudinal'); amplitude = [0.025_dp, 0.025_dp, 0.05_dp]
  case ('lateral');      amplitude = [0.025_dp, 0.05_dp, 0.025_dp]
  case default
     error stop 'rotorcraft_dynamics: the scenario is collective, longitudinal or lateral'
  end select

  ! the revolutions of the hub, a manifold of time alone; with the
  ! blades as beams, of time and of the span from the root to the tip
  num_instants = revolutions * steps + 1
  np = 1
  select case (trim(blade))
  case ('rigid')
     omega   = continuous_manifold(time=interval(0.0_dp, revolutions * 2.0_dp * pi / speed))
     omega_h = omega % discretize(time=points(num_instants))
  case ('beam')
     omega   = continuous_manifold(time=interval(0.0_dp, revolutions * 2.0_dp * pi / speed), &
          &                        space=region(interval(0.0_dp, tip - root)))
     omega_h = omega % discretize(time=points(num_instants), space=points(span_points))
     x       = omega % coordinate('x')
     left    = omega % boundary(x, 0.0_dp)
     np      = span_points
  case default
     error stop 'rotorcraft_dynamics: the blade model is rigid or beam; plate and solid are not yet stated'
  end select
  domega  = omega % boundary(omega % coordinate('t'), 0.0_dp)
  t       = omega % coordinate('t')

  ! THE BODIES: the table's, each turning at its rate about z
  rotor = mechanism(t)
  table = body_table('rotor_bodies.txt')
  allocate(part(size(table % name)))
  do k = 1, size(table % name)
     call rotor % body(omega, part(k), trim(table % name(k)), table % centre(:, k), table % mass(k), &
          & table % inertia(:, :, k), angular_velocity=table % rate(k) * ez)
  end do

  ! the blades, each along the azimuth of its cap and turning with it:
  ! rigid plates linked to the caps, or beams in the caps' axes from
  ! the root, the span along the azimuth, the chord across it and the
  ! thickness along z
  select case (trim(blade))
  case ('rigid')
     do k = 1, 4
        call rotor % body(omega, blades(k), 'blade_' // trim(degrees(k - 1)), &
             & turned([(root + tip) / 2.0_dp, 0.0_dp, height], k - 1), &
             & density * thickness * (tip - root) * chord, plate(k - 1), angular_velocity=speed * ez)
     end do
  case ('beam')
     do k = 1, 4
        call rotor % beam(omega, beams(k), part(table % row('cap_' // trim(degrees(k - 1)))), &
             & 'blade_' // trim(degrees(k - 1)), x, turned([root, 0.0_dp, height], k - 1), &
             & reshape([turned(ex, k - 1), turned(ey, k - 1), ez], [3, 3]), density, modulus, poisson, chord, thickness)
     end do
  end select

  ! THE JOINTS
  associate (hub => part(table % row('hub')), lsp => part(table % row('lower_swashplate')), &
       &     usp => part(table % row('upper_swashplate')), sphere => part(table % row('sphere')), &
       &     bsp => part(table % row('base_plate')), lph => part(table % row('lower_horn')), &
       &     uph => part(table % row('upper_horn')))

    ! the shaft turns at its rate; the swashplates hinge about z, the
    ! lower one on the sphere, which slides along the shaft's axis
    call rotor % joint(omega, 'hub_driver', revolute_driver(hub, ez, speed, t))
    call rotor % joint(omega, 'swashplates', revolute(lsp, [0.0_dp, 0.0_dp, 0.4_dp], ez, usp))
    call rotor % joint(omega, 'swashplate_sphere', spherical(lsp, [0.0_dp, 0.0_dp, 0.4_dp], sphere))
    call rotor % joint(omega, 'sphere_shaft', cylindrical(sphere, [0.0_dp, 0.0_dp, 0.4_dp], ez))

    ! at each of the four azimuths: the upper swashplate to the lower
    ! pitch link, the lower to the upper pitch link, the upper pitch
    ! link hinged to the shaft, the cap and the blade linked to it
    do k = 1, 4
       lower = table % row('lower_link_' // trim(degrees(k - 1, 30)))
       upper = table % row('upper_link_' // trim(degrees(k - 1, 30)))
       cap   = table % row('cap_' // trim(degrees(k - 1)))
       call rotor % joint(omega, 'swashplate_link_' // trim(degrees(k - 1)), &
            & spherical(usp, turned([0.278167359696_dp, 0.1606_dp, 0.43_dp], k - 1), part(lower)))
       call rotor % joint(omega, 'links_' // trim(degrees(k - 1)), &
            & spherical(part(lower), turned([0.278167359696_dp, 0.1606_dp, height], k - 1), part(upper)))
       call rotor % joint(omega, 'shaft_link_' // trim(degrees(k - 1)), &
            & revolute(hub, turned([0.231_dp, 0.0_dp, height], k - 1), abs(turned(ex, k - 1)), part(upper)))
       call rotor % joint(omega, 'link_cap_' // trim(degrees(k - 1)), rigid_link(part(upper), part(cap)))
       if (trim(blade) == 'rigid') call rotor % joint(omega, 'cap_blade_' // trim(degrees(k - 1)), &
            & rigid_link(part(cap), blades(k)))
    end do

    ! the push rods on the lower swashplate, driven along z
    do k = 1, 3
       rod = table % row('push_rod_' // trim(degrees(k)))
       call rotor % joint(omega, 'rod_swashplate_' // trim(degrees(k)), &
            & spherical(part(rod), turned([0.4004_dp, 0.0_dp, 0.38_dp], k), lsp))
       call rotor % joint(omega, 'rod_driver_' // trim(degrees(k)), &
            & motion_driver(part(rod), amplitude(k) * ez, speed / 4.0_dp, t))
    end do

    ! the base plate fixed, the horns hinged to it and to each other
    ! about y, the upper horn on the lower swashplate
    call rotor % joint(omega, 'base', fixed(bsp))
    call rotor % joint(omega, 'base_horn', revolute(lph, [0.4004_dp, 0.0_dp, 0.05_dp], ey, bsp))
    call rotor % joint(omega, 'horns', revolute(lph, [0.5654_dp, 0.0_dp, 0.215_dp], ey, uph))
    call rotor % joint(omega, 'horn_swashplate', spherical(uph, [0.4004_dp, 0.0_dp, 0.38_dp], lsp))

  end associate

  ! the surfaces the bodies are drawn with
  if (len_trim(meshes) > 0) then
     do k = 1, size(table % name)
        call rotor % draw_with(part(k), mesh_surface(trim(meshes) // '/' // trim(table % mesh(k))))
     end do
     do k = 1, 4
        if (trim(blade) == 'rigid') then
           call rotor % draw_with(blades(k), plate_surface(k - 1))
        else
           call rotor % draw_beam_with(beams(k), part(table % row('cap_' // trim(degrees(k - 1)))), &
                & plate_surface(k - 1))
        end if
     end do
  end if

  ! the laws on the revolution, the initial data on {0}, and with the
  ! blades as beams their conditions at the root, the derivatives
  ! along the span by the hermite cubics
  r = continuous_residual(omega,  rotor % laws())
  g = continuous_residual(domega, rotor % initial)
  L = r + g
  if (trim(blade) == 'rigid') then
     L_h = L % discretize(omega_h, time=chain([dirk(2), (bdf(1), k = 2, num_instants)], &
          &                                   from=[(k, k = 1, num_instants)]))
  else
     L = L + continuous_residual(left, rotor % root)
     L_h = L % discretize(omega_h, time=chain([dirk(2), (bdf(1), k = 2, num_instants)], &
          &                                   from=[(k, k = 1, num_instants)]), space=finite_element(degree=3))
  end if

  motion   = rotor % estimate(omega)
  estimate = motion % discretize(omega_h)
  if (this_image() == 1 .and. num_instants <= 11) call set_verbosity(1)
  call L_h % minimize(estimate, solution, continued=.true.)

  ! the energy at every instant - of the rigid bodies alone when the
  ! blades are beams - and the hub's angle; the first point of each
  ! instant, and with beams the tip of the blade at zero degrees
  energy   = rotor % energy
  energy_h = energy % at(solution)
  allocate(e(num_instants * np), hub_eps(num_instants * np))
  call energy_h % values(e)
  hub_h = solution % fields(['hub_eps3'])
  call hub_h % values(hub_eps)
  print '(a, a, a, a, a, i0, a)', 'the rotor, ', trim(blade), ' blades, ', trim(scenario), ', ', num_instants, &
       & ' instants:'
  n = max(1, (num_instants - 1) / (12 * revolutions))
  if (trim(blade) == 'rigid') then
     print '(a)', '   angle         energy   hub eps3 - sin(w t/2)'
     do k = 1, num_instants, n
        print '(f8.1, es16.8, es14.3)', speed * omega_h % instant(k) * 180.0_dp / pi, e(k), &
             & hub_eps(k) - sin(speed * omega_h % instant(k) / 2.0_dp)
     end do
  else
     names = [character(len=32) :: 'blade_0_u', 'blade_0_v', 'blade_0_w', 'blade_0_phi']
     tip_h = solution % fields(names)
     print '(a)', '   angle   bodies'' energy   hub eps3 - sin(w t/2)    the tip of the blade at 0: u, v, w, phi'
     do k = 1, num_instants, n
        print '(f8.1, es16.8, es14.3, 4es14.5)', speed * omega_h % instant(k) * 180.0_dp / pi, e((k - 1) * np + 1), &
             & hub_eps((k - 1) * np + 1) - sin(speed * omega_h % instant(k) / 2.0_dp), tip_h % value(k * np, :)
     end do
  end if

  ! the time history
  if (len_trim(meshes) > 0) call rotor % draw(solution, 'rotor', every)

contains

  ! the surface of a body's mesh
  function mesh_surface(file) result(s)
    character(len=*), intent(in) :: file
    type(surface) :: s
    call bdf_surface(file, s % points, s % corners, s % count, s % kinds)
  end function mesh_surface

  ! the surface of the plate of the blade of k quarter turns: the box
  ! of its span, chord and thickness, 48 cells along the span and 4
  ! along the chord on each face
  function plate_surface(k) result(s)
    integer, intent(in) :: k
    type(surface) :: s
    integer, parameter :: ns = 48, nc = 4
    real(dp) :: x(3, 8), corner(3)
    integer :: i, j, f, face(4, 6), n1, n2, p, c
    ! the corners of the box along the azimuth zero, and its faces
    do i = 1, 8
       corner = [merge(tip, root, btest(i - 1, 0)), merge(chord, -chord, btest(i - 1, 1)) / 2.0_dp, &
            & height + merge(thickness, -thickness, btest(i - 1, 2)) / 2.0_dp]
       x(:, i) = corner
    end do
    face = reshape([1, 2, 4, 3,  5, 6, 8, 7,  1, 2, 6, 5,  3, 4, 8, 7,  1, 3, 7, 5,  2, 4, 8, 6], [4, 6])
    allocate(s % points(3, 0), s % corners(4, 0), s % count(0), s % kinds(0))
    do f = 1, 6
       ! the faces of the span and the chord, of the span and the
       ! thickness, of the chord and the thickness
       n1 = merge(ns, nc, f <= 4)
       n2 = merge(nc, 1, f <= 2)
       p = size(s % points, 2)
       do j = 0, n2
          do i = 0, n1
             s % points = reshape([s % points, turned(x(:, face(1, f)) + (x(:, face(2, f)) - x(:, face(1, f))) &
                  & * real(i, dp) / n1 + (x(:, face(4, f)) - x(:, face(1, f))) * real(j, dp) / n2, k)], &
                  & [3, size(s % points, 2) + 1])
          end do
       end do
       do j = 1, n2
          do i = 1, n1
             c = p + (j - 1) * (n1 + 1) + i
             s % corners = reshape([s % corners, c, c + 1, c + n1 + 2, c + n1 + 1], [4, size(s % count) + 1])
             s % count   = [s % count, 4]
             s % kinds   = [s % kinds, 3]
          end do
       end do
    end do
  end function plate_surface

  ! the inertia of the blade of the azimuth of k quarter turns about its
  ! centre of mass: a plate of the span, the chord and the thickness
  function plate(k) result(inertia)
    integer, intent(in) :: k
    real(dp) :: inertia(3,3), mass, span
    span    = tip - root
    mass    = density * thickness * span * chord
    inertia = 0.0_dp
    inertia(1, 1) = mass * (chord**2 + thickness**2) / 12.0_dp
    inertia(2, 2) = mass * (span**2 + thickness**2) / 12.0_dp
    inertia(3, 3) = mass * (span**2 + chord**2) / 12.0_dp
    if (mod(k, 2) == 1) inertia(1:2, 1:2) = reshape([inertia(2, 2), 0.0_dp, 0.0_dp, inertia(1, 1)], [2, 2])
  end function plate

  ! the point x turned about z by k quarter turns
  function turned(x, k) result(y)
    real(dp), intent(in) :: x(3)
    integer , intent(in) :: k
    real(dp) :: y(3)
    integer :: i
    y = x
    do i = 1, modulo(k, 4)
       y = [-y(2), y(1), y(3)]
    end do
  end function turned

  ! the azimuth of k quarter turns in degrees, offset by the degrees
  ! given, as the bodies are named
  function degrees(k, offset) result(name)
    integer, intent(in) :: k
    integer, intent(in), optional :: offset
    character(len=8) :: name
    integer :: angle
    angle = 90 * k
    if (present(offset)) angle = angle + offset
    write(name, '(i0)') angle
  end function degrees

end program rotorcraft_dynamics
