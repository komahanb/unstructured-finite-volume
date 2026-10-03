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
! thickness of the setup, density 2500; each lies along the azimuth of
! its cap, 0, 90, 180 and 270 degrees, and is linked to that cap. The
! blade model is chosen at setup: rigid, each blade a rigid body of
! the plate's mass and inertia; beam, plate and solid are not yet
! stated.
!
! The Lagrangian density is the bodies' kinetic energy less their
! potential energy, with every condition of the joints times its
! reaction; its stationarity in each unknown is that unknown's law.
! The motion is marched by bdf(1) over the revolution, each instant
! from the solution at the one before. The energy at the first instant
! is 4.99144494e5 J, and over the revolution of 360 steps it falls to
! 4.94255e5, the dissipation of bdf(1); the hub's angle holds to
! 1e-16 and every condition of the joints to 2e-14 at every instant.
!
! THE SETUP, rotorcraft.nml:
!
!      blade        rigid (beam, plate and solid to come)
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
!
! THE DRAWING: every rigid body is drawn with the surface of its mesh,
! each blade with the box of its plate, moved with the body to each
! instant drawn: rotor_<instant>.vtu and rotor.pvd for paraview.
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program rotorcraft_dynamics

  use util_precision    , only : dp
  use util_verbosity    , only : set_verbosity
  use operation_manifold, only : continuous_manifold, discrete_manifold, interval, points
  use operation_field   , only : continuous_field, discrete_field
  use operation_residual, only : continuous_residual, discrete_residual, operator(+)
  use operation_family  , only : dirk, bdf, chain
  use multibody_dynamics, only : rigid_body, mechanism, body_table, surface, spherical, revolute, cylindrical, &
       &                         fixed, rigid_link, revolute_driver, motion_driver
  use view_bdf_loader   , only : bdf_surface

  implicit none

  real(dp), parameter :: pi = acos(-1.0_dp), speed = 109.12_dp
  real(dp), parameter :: ex(3) = [1.0_dp, 0.0_dp, 0.0_dp], ey(3) = [0.0_dp, 1.0_dp, 0.0_dp], &
       &                 ez(3) = [0.0_dp, 0.0_dp, 1.0_dp], origin(3) = 0.0_dp
  real(dp), parameter :: density = 2500.0_dp, root = 0.44_dp, tip = 2.0_dp, chord = 0.121_dp, height = 0.9_dp

  ! the setup
  character(len=16)  :: blade = 'rigid', scenario = 'collective'
  character(len=256) :: meshes = ''
  real(dp)           :: thickness = 0.02_dp
  integer            :: revolutions = 1, steps = 360, every = 10
  namelist /setup/ blade, scenario, thickness, revolutions, steps, meshes, every

  type(continuous_manifold) :: omega, domega
  type(discrete_manifold)   :: omega_h
  type(continuous_field)    :: t, motion, energy
  type(discrete_field)      :: estimate, solution, energy_h, hub_h
  type(continuous_residual) :: r, g, L
  type(discrete_residual)   :: L_h
  type(mechanism)           :: rotor
  type(body_table)          :: table
  type(rigid_body), allocatable :: part(:)
  type(rigid_body)          :: blades(4)
  real(dp)                  :: amplitude(3)
  real(dp), allocatable     :: e(:), hub_eps(:)
  integer :: unit, k, n, num_instants, lower, upper, cap, rod

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

  ! one revolution of the hub, a manifold of time alone
  num_instants = revolutions * steps + 1
  omega   = continuous_manifold(time=interval(0.0_dp, revolutions * 2.0_dp * pi / speed))
  domega  = omega % boundary(omega % coordinate('t'), 0.0_dp)
  omega_h = omega % discretize(time=points(num_instants))
  t       = omega % coordinate('t')

  ! THE BODIES: the table's, each turning at its rate about z
  rotor = mechanism(t)
  table = body_table('rotor_bodies.txt')
  allocate(part(size(table % name)))
  do k = 1, size(table % name)
     call rotor % body(omega, part(k), trim(table % name(k)), table % centre(:, k), table % mass(k), &
          & table % inertia(:, :, k), angular_velocity=table % rate(k) * ez)
  end do

  ! the blades, each along the azimuth of its cap and turning with it
  select case (trim(blade))
  case ('rigid')
     do k = 1, 4
        call rotor % body(omega, blades(k), 'blade_' // trim(degrees(k - 1)), &
             & turned([(root + tip) / 2.0_dp, 0.0_dp, height], k - 1), &
             & density * thickness * (tip - root) * chord, plate(k - 1), angular_velocity=speed * ez)
     end do
  case default
     error stop 'rotorcraft_dynamics: the blade model is rigid; beam, plate and solid are not yet stated'
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
       call rotor % joint(omega, 'cap_blade_' // trim(degrees(k - 1)), rigid_link(part(cap), blades(k)))
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
        call rotor % draw_with(blades(k), plate_surface(k - 1))
     end do
  end if

  ! the laws on the revolution, the initial data on {0}
  r = continuous_residual(omega,  rotor % laws())
  g = continuous_residual(domega, rotor % initial)
  L = r + g
  L_h = L % discretize(omega_h, time=chain([dirk(2), (bdf(1), k = 2, num_instants)], &
       &                                   from=[(k, k = 1, num_instants)]))

  motion   = rotor % estimate(omega)
  estimate = motion % discretize(omega_h)
  if (this_image() == 1 .and. num_instants <= 11) call set_verbosity(1)
  call L_h % minimize(estimate, solution, continued=.true.)

  ! the energy at every instant, and the hub's angle
  energy   = rotor % energy
  energy_h = energy % at(solution)
  allocate(e(num_instants), hub_eps(num_instants))
  call energy_h % values(e)
  hub_h = solution % fields(['hub_eps3'])
  call hub_h % values(hub_eps)
  print '(a, a, a, a, a, i0, a)', 'the rotor, ', trim(blade), ' blades, ', trim(scenario), ', ', num_instants, &
       & ' instants:'
  print '(a)', '   angle         energy   hub eps3 - sin(w t/2)'
  n = max(1, (num_instants - 1) / (12 * revolutions))
  do k = 1, num_instants, n
     print '(f8.1, es16.8, es14.3)', speed * omega_h % instant(k) * 180.0_dp / pi, e(k), &
          & hub_eps(k) - sin(speed * omega_h % instant(k) / 2.0_dp)
  end do

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
