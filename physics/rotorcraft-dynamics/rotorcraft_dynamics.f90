!=====================================================================!
! THE DYNAMICS OF A FOUR-BLADED ROTOR (dissertation, section 5.4): the
! case of komahanb/rotorcraft-optimization, its 26 bodies and 34
! joints as collective.py, longitudinal_cyclic.py and
! lateral_cyclic.py state them, and a 35th, the rotating scissors.
!
! The hub turns at 109.12 rad/s about z. Three push rods move the
! lower swashplate along z by a sin(w t), w a quarter of the hub's
! rate; the upper swashplate hinges about z on the lower one, the
! rotating scissors turns it with the hub, and through the lower and
! upper pitch links and the blade caps it sets the pitch of the
! blades. The base plate is fixed, the sphere slides along the
! shaft's axis, and two horns hinge the lower swashplate to the base
! plate. Without the scissors, as in TACS, the upper swashplate's turn
! is free, held to the hub by the pitch links alone as they tilt: its
! lag behind the hub after a revolution is the march's error in its
! rate, 64.4 degrees at 36 steps of bdf(1), 10.1 at 360, 8.60 at 36
! steps of bdf(2) and 2.14 at 72, the pitch links tilted by 43.3,
! 6.1, 5.09 and 0.77 degrees. The bodies begin at the velocities the
! drivers' rates require (multibody_dynamics, joint), the push rods'
! a w along z where the table states them at rest.
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
! The motion is marched by bdf(2) over the revolution, from one step
! of bdf(1). The energy at the first instant is 4.99432169e5 J,
! 287.675 J above that of the table's rates, the motion the push rods'
! rates require; at the end of a revolution of 360 steps it is
! 4.99195025e5, the drivers' work and the march's error between, and
! the hub's angle holds to 2.4e-18 at every instant printed. With the
! blades as beams the march is solved with the rows that tie the
! copies of the unknowns of time alone eliminated, the rest by a
! direct solve; at 30 degrees of 144 steps the tip of the blade at
! zero degrees has moved by 7.82649e-4 along the span, about the
! static stretch of the turning blade, the integral over the span of
! the centrifugal stress rho w^2 (R^2 - r^2)/2 over E, 7.6583e-4, and
! along the chord by -2.13686e-3, -2.13816e-3, -2.13845e-3 and
! -2.13847e-3 at 3, 5, 9 and 17 points of the span, the differences
! falling by 4.5 and 15.
!
! THE SETUP, rotorcraft.nml:
!
!      blade        rigid or beam (plate and solid to come)
!      scenario     collective, longitudinal or lateral: the push
!                   rods' amplitudes, 0.05 each, or 0.025, 0.025 and
!                   0.05, or 0.025, 0.05 and 0.025 at 90, 180 and 270
!                   degrees; several listed, each is stated on a
!                   manifold of its own and marched in turn
!      thickness    of the blades; with strips, of each strip
!      revolutions  of the hub
!      steps        per revolution
!      meshes       the directory of the bodies' surface meshes, the
!                   case's data/mesh; none, no drawing
!      every        the stride of the instants drawn
!      span_points  of a blade as a beam
!      strips       of the span of a blade as a beam, the thickness of
!                   each a design coordinate, the same for the four
!                   blades; none, the thickness of the setup throughout
!      yield        the yield stress, 350e6
!      ks_weight    the weight of the smooth maximum of the stress
!      sensitivities with strips, the functions' derivatives along each
!                   thickness (.true.), or the functions alone
!      optimize     with strips, the thickness of the least mass whose
!                   smooth maximum of the von Mises function is at
!                   most ks_bound^2 in every scenario listed, each
!                   strip's within bounds and adjacent strips' within
!                   smoothness of each other: one design of the
!                   scenarios' Lagrangians (minimize_together), the
!                   first stating the mass and the inequalities on
!                   the design alone
!      ks_bound     the bound of the von Mises ratio, 0.25
!      bounds       the least and the greatest thickness, 0.01 and 0.02
!      smoothness   the largest difference of adjacent strips, 1e-3
!      order        of the backward differences of the march, 1 or 2
!                   (2), from one step of bdf(1)
!      lagrangian   the Lagrangian of each scenario written as a
!                   formula, one equation a line, to
!                   rotor_lagrangian.txt (.true.), or not (.false.)
!
! THE FUNCTIONS, of the blades as beams: the blades' mass over the
! revolutions, and the smooth maximum (operation_field, maximum) of
! the von Mises function of their stress over the four blades, the
! corners and the middles of the long sides of their sections and the
! revolutions; with strips, their derivatives along the thickness of
! each strip by the expansion of the motion, and the smooth maximum's
! by the reverse sweep, -kappa, its objective. At 144 steps, 5 points
! of the span and eight strips at 2 cm the mass is 2.1737794328870783
! kg s, TACS's 2.173779432887122, and the smooth maximum at ks_weight =
! 100 2.5151758681472021e-2, the same on 1 and 2 images: the stress is
! the centrifugal stress of the turning blade, rho w^2 (R^2 - r^2)/2,
! 0.162 of the yield at the root at every thickness, its square there
! 0.0262; the smooth maximum is 2.51151e-2, 2.51518e-2, 2.51569e-2 and
! 2.51574e-2 at 3, 5, 9 and 17 points of the span, and 3.961e-2,
! 2.515e-2 and 2.504e-2 at 72, 144 and 288 steps. Along the root
! strip's thickness its derivative by the expansion is
! -1.2065199020683000, a thinner root strip under the outer strips'
! load, by the reverse sweep the same to 2.2e-15 (along the eight
! strips to 5.0e-15), their sum 1.13e-3, a uniform thickness leaving
! the stress as it is; central differences of steps 5e-6, 1e-5 and
! 2e-5 miss it by their own error, 7.681e-7, 3.073e-6 and 1.229e-5,
! extrapolated in the step's second and fourth powers to within
! 2.0e-10.
!
! THE LEAST MASS (optimize), at 144 steps, 5 points of the span and
! eight strips from 2 cm, with ks_bound = 0.25: the smooth maximum, the
! square of the centrifugal stress, is below the bound's 0.0625 at
! every uniform thickness from 1 to 2 cm, and the design found in two
! iterations is every strip at its bound of 1 cm, the mass
! 1.0868897164435392 kg s, the multiplier of the stress zero and that
! of each strip's bound 13.586121455544628, the mass's derivative along
! it, 4 rho c (L/8) 2 pi/w = 13.5861. The three scenarios together,
! their smooth maxima at 1 cm 2.51338e-2, 2.51312e-2 and 2.51287e-2:
! the same design in two iterations. At 36 steps of bdf(1) without the
! scissors the smooth maximum was 1.676 at 2 cm, the stress of the
! upper swashplate's lag, and the least mass under the bound 1.25 in
! its place left the root strip at 1.02 cm.
!
! THE DRAWING: every rigid body is drawn with the surface of its mesh,
! each blade with the box of its plate, its thickness that of each
! strip, moved with the body to each instant drawn; a blade as a beam
! with each point of the box moved with its section of the beam, and
! the beam's deflection a vector at the point: rotor_<instant>.vtu and
! rotor.pvd for paraview, rotor_<scenario>_... where several scenarios
! are listed. To optimize, the initial design is marched, reported and
! drawn first, rotor_initial, and the design found drawn with the
! thicknesses found, rotor_optimized. README.md states how to run and
! draw.
!
! THE LAGRANGIAN (lagrangian): L written as a formula before its
! march, one equation a line: the rigid rotor's 356 laws, the
! stationarity of the density in each unknown, and the 364 initial
! data of its 26 bodies, without multipliers, 720 lines of 876 kB;
! the blades as beams with eight strips, the smooth maximum and every
! law times its multiplier, 689 lines of 2.0 MB. Its tree is not
! written: that of the density alone is 13 921 lines of 10 MB.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program rotorcraft_dynamics

  use util_precision    , only : dp
  use util_verbosity    , only : set_verbosity
  use operation_manifold, only : continuous_manifold, discrete_manifold, interval, region, points, parameter
  use operation_expansion, only : expansion
  use operation_field   , only : continuous_field, discrete_field, integral, maximum, heaviside, &
       &                         operator(+), operator(-), operator(*)
  use operation_lagrangian, only : continuous_lagrangian, discrete_lagrangian, operator(+), operator(*), &
       &                           minimize_together
  use operation_family  , only : family, dirk, bdf, chain
  use operation_finite_element, only : finite_element
  use multibody_dynamics, only : rigid_body, beam_blade, mechanism, body_table, surface, spherical, revolute, &
       &                         cylindrical, fixed, rigid_link, revolute_driver, motion_driver, scissors
  use view_bdf_loader   , only : bdf_surface
  use view_residual     , only : formula, tree
  use operation_newton      , only : newton
  use operation_elimination , only : elimination
  use operation_dense_direct, only : dense_direct
  use operation_minimization, only : discrete_minimizer

  implicit none

  real(dp), parameter :: pi = acos(-1.0_dp), speed = 109.12_dp
  real(dp), parameter :: ex(3) = [1.0_dp, 0.0_dp, 0.0_dp], ey(3) = [0.0_dp, 1.0_dp, 0.0_dp], &
       &                 ez(3) = [0.0_dp, 0.0_dp, 1.0_dp], origin(3) = 0.0_dp
  real(dp), parameter :: density = 2500.0_dp, root = 0.44_dp, tip = 2.0_dp, chord = 0.121_dp, height = 0.9_dp
  real(dp), parameter :: modulus = 70.0e9_dp, poisson = 0.3_dp

  ! the setup
  character(len=16)  :: blade = 'rigid', scenario(3) = [character(len=16) :: 'collective', '', '']
  character(len=256) :: meshes = ''
  real(dp)           :: thickness(48) = 0.02_dp, yield = 350.0e6_dp, ks_weight = 100.0_dp
  real(dp)           :: ks_bound = 0.25_dp, bounds(2) = [0.01_dp, 0.02_dp], smoothness = 1.0e-3_dp
  integer            :: revolutions = 1, steps = 360, every = 10, span_points = 9, strips = 0, order = 2
  logical            :: sensitivities = .true., optimize = .false., lagrangian = .false.
  namelist /setup/ blade, scenario, thickness, revolutions, steps, meshes, every, span_points, strips, yield, &
       & ks_weight, sensitivities, optimize, ks_bound, bounds, smoothness, order, lagrangian

  type(discrete_manifold) :: omega_h
  type(discrete_lagrangian) :: L_h(size(scenario))
  type(discrete_field)    :: estimates(size(scenario)), solutions(size(scenario))
  type(mechanism)         :: rotors(size(scenario))
  type(continuous_field)  :: masses(size(scenario)), kss(size(scenario))
  integer :: unit, j, ns, num_instants, np
  logical :: optimizing
  type(discrete_minimizer) :: M_h

  open(newunit=unit, file='rotorcraft.nml', status='old', action='read')
  read(unit, nml=setup)
  close(unit)

  ! the scenarios listed, each on a manifold of its own
  ns = count(len_trim(scenario) > 0)
  if (ns == 0) error stop 'rotorcraft_dynamics: the setup lists a scenario at least'
  num_instants = revolutions * steps + 1
  np = 1
  if (trim(blade) == 'beam') np = span_points
  optimizing = strips > 0 .and. optimize .and. trim(blade) == 'beam'

  ! the convergence of every solve and of the design's iteration
  ! printed by image 1
  if (this_image() == 1) call set_verbosity(1)

  ! THE MINIMIZER M_h, paired with each L_h once L_h is stated
  ! (operation_minimization).
  ! Each block of instants is solved by Newton's method to 1e-13.
  ! The chord method reuses the LU factors of the first jacobian of a
  ! block for every Newton step of that block.
  ! A Newton step first eliminates the tying rows.
  ! A tying row fixes one unknown from unknowns before it: by the
  ! history, by the relations of the time scheme, or as a copy.
  ! A dense LU solve gives the remaining unknowns.
  ! For blades on cells, the rows mix positions, rotations and
  ! reactions, of different units.
  ! The LU solve then scales the rows, then the columns, to a largest
  ! entry of one (equilibrated).

  ! to optimize with the bodies drawn, the initial design marched,
  ! reported and drawn first, beside the design found
  if (optimizing .and. len_trim(meshes) > 0) then
     do j = 1, ns
        call stated(scenario(j), j == 1, L_h(j), estimates(j), rotors(j), masses(j), kss(j), omega_h, initial=.true.)
        if (this_image() == 1) print '(a, a, a, a, a)', 'the march of the rotor, ', trim(blade), ' blades, ', &
             & trim(scenario(j)), ', the initial design'
        M_h = discrete_minimizer(L_h(j), newton(elimination(dense_direct(equilibrated=trim(blade) /= 'rigid')), &
             &                                  tolerance=1.0e-13_dp, chord=.true.))
        call L_h(j) % minimize(estimates(j), solutions(j), continued=.true., by=M_h)
        call reported(scenario(j), ns == 1, rotors(j), masses(j), kss(j), solutions(j), 'initial')
     end do
  end if
  do j = 1, ns
     call stated(scenario(j), j == 1, L_h(j), estimates(j), rotors(j), masses(j), kss(j), omega_h)
  end do

  ! the minimizer of every block: the march over the blocks, and the
  ! minimizer of a block
  M_h = discrete_minimizer(L_h(1), newton(elimination(dense_direct(equilibrated=trim(blade) /= 'rigid')), &
       &                               tolerance=1.0e-13_dp, chord=.true.))
  if (this_image() == 1) print '(a)', tree(M_h)

  ! each marched in turn; to optimize, the one design of them all
  if (optimizing) then
     if (this_image() == 1) print '(*(a))', 'the design of the rotor, ', trim(blade), ' blades, of the scenarios', &
          & (' ' // trim(scenario(j)), j = 1, ns)
     call minimize_together(L_h(1:ns), estimates(1:ns), solutions(1:ns), continued=.true., by=M_h)
     do j = 1, ns
        call reported(scenario(j), ns == 1, rotors(j), masses(j), kss(j), solutions(j), 'optimized')
     end do
  else
     do j = 1, ns
        if (this_image() == 1) print '(a, a, a, a)', 'the march of the rotor, ', trim(blade), ' blades, ', &
             & trim(scenario(j))
        call L_h(j) % minimize(estimates(j), solutions(j), continued=.true., by=M_h)
     end do
     do j = 1, ns
        call reported(scenario(j), ns == 1, rotors(j), masses(j), kss(j), solutions(j))
     end do
  end if

contains

  !===================================================================!
  ! ONE SCENARIO STATED on a manifold of its own: the revolutions, the
  ! bodies, the blades and the joints with the push rods' amplitudes of
  ! the scenario, the functions, the Lagrangian discretized and the
  ! estimate. To optimize, every scenario bounds its stress, and the
  ! first states the objective and the inequalities on the design
  ! alone; initial, the march at the setup's thickness alone.
  !===================================================================!

  subroutine stated(name, first, L_h, estimate, rotor, mass, ks, omega_h, initial)

    character(len=*)       , intent(in)  :: name
    logical                , intent(in)  :: first
    type(discrete_lagrangian), intent(out) :: L_h
    type(discrete_field)   , intent(out) :: estimate
    type(mechanism)        , intent(out) :: rotor
    type(continuous_field) , intent(out) :: mass, ks
    type(discrete_manifold), intent(out) :: omega_h
    logical                , intent(in), optional :: initial

    type(continuous_manifold) :: omega, domega, left
    type(continuous_field)    :: t, x, motion, h, adjoint, lambda, nu, kappa, mu_ks, mu_bounds, mu_smooth, bounded
    type(continuous_field), allocatable :: stresses(:), thicknesses(:), limits(:)
    type(continuous_lagrangian) :: r, g, b, d, L, stressed, within, smooth
    type(chain)               :: march
    type(body_table)          :: table
    type(rigid_body), allocatable :: part(:)
    type(rigid_body)          :: blades(4)
    type(beam_blade)          :: beams(4)
    real(dp)                  :: amplitude(3)
    type(points)              :: instants
    character(len=:), allocatable :: file
    integer :: k, lower, upper, cap, rod, unit
    logical :: marched

    marched = .false.
    if (present(initial)) marched = initial

    select case (trim(name))
    case ('collective');   amplitude = [0.05_dp, 0.05_dp, 0.05_dp]
    case ('longitudinal'); amplitude = [0.025_dp, 0.025_dp, 0.05_dp]
    case ('lateral');      amplitude = [0.025_dp, 0.05_dp, 0.025_dp]
    case default
       error stop 'rotorcraft_dynamics: the scenario is collective, longitudinal or lateral'
    end select

    ! the revolutions of the hub, a manifold of time alone; with the
    ! blades as beams, of time and of the span from the root to the tip
    instants = points(num_instants)

    select case (trim(blade))
    case ('rigid')
       omega   = continuous_manifold(time=interval(0.0_dp, revolutions * 2.0_dp * pi / speed))
       omega_h = omega % discretize(time=instants)
    case ('beam')
       if (strips > 0) then
          omega   = continuous_manifold(time=interval(0.0_dp, revolutions * 2.0_dp * pi / speed), &
               &                        space=region(interval(0.0_dp, tip - root)), &
               &                        design=[(parameter(numbered('h', k), thickness(k)), k = 1, strips)])
          ! the derivatives by the expansion of order one; to optimize, by
          ! the reverse sweeps alone, the hessian by its BFGS update
          if (sensitivities .and. .not. optimize) then
             omega_h = omega % discretize(time=instants, space=points(span_points), &
                  &                       design=expansion(order=1))
          else
             omega_h = omega % discretize(time=instants, space=points(span_points), &
                  &                       design=expansion(order=0))
          end if
       else
          omega   = continuous_manifold(time=interval(0.0_dp, revolutions * 2.0_dp * pi / speed), &
               &                        space=region(interval(0.0_dp, tip - root)))
          omega_h = omega % discretize(time=instants, space=points(span_points))
       end if
       x       = omega % coordinate('x')
       left    = omega % boundary(x, 0.0_dp)
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
               & density * thickness(1) * (tip - root) * chord, plate(k - 1), angular_velocity=speed * ez)
       end do
    case ('beam')
       ! the thickness: the setup's, or with strips of the span the
       ! thickness of each strip a design coordinate, the same for the
       ! four blades, h = h_1 + sum over j > 1 of (h_j - h_(j-1)) H(x -
       ! x_j), x_j the inner edge of strip j
       if (strips > 0) then
          h = omega % coordinate(numbered('h', 1))
          do k = 2, strips
             h = h + (omega % coordinate(numbered('h', k)) - omega % coordinate(numbered('h', k - 1))) &
                  & * heaviside(x - real(k - 1, dp) * (tip - root) / real(strips, dp))
          end do
       else
          h = continuous_field(omega, [thickness(1)])
       end if
       do k = 1, 4
          call rotor % beam(omega, beams(k), part(table % row('cap_' // trim(degrees(k - 1)))), &
               & 'blade_' // trim(degrees(k - 1)), x, turned([root, 0.0_dp, height], k - 1), &
               & reshape([turned(ex, k - 1), turned(ey, k - 1), ez], [3, 3]), density, modulus, poisson, chord, h)
       end do
    end select

    ! THE JOINTS
    associate (hub => part(table % row('hub')), lsp => part(table % row('lower_swashplate')), &
         &     usp => part(table % row('upper_swashplate')), sphere => part(table % row('sphere')), &
         &     bsp => part(table % row('base_plate')), lph => part(table % row('lower_horn')), &
         &     uph => part(table % row('upper_horn')))

      ! the shaft turns at its rate; the swashplates hinge about z, the
      ! upper one turned with the shaft by the rotating scissors, the
      ! lower one on the sphere, which slides along the shaft's axis
      call rotor % joint(omega, 'hub_driver', revolute_driver(hub, ez, speed, t))
      call rotor % joint(omega, 'swashplates', revolute(lsp, [0.0_dp, 0.0_dp, 0.4_dp], ez, usp))
      call rotor % joint(omega, 'scissors', scissors(hub, ez, usp))
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
             call rotor % draw_with(blades(k), plate_surface(k - 1, thickness(1:1)))
          else
             call rotor % draw_beam_with(beams(k), part(table % row('cap_' // trim(degrees(k - 1)))), &
                  & plate_surface(k - 1, thickness(1:max(strips, 1))))
          end if
       end do
    end if

    ! with the blades as beams, the blades' mass over the revolutions and
    ! the smooth maximum of the von Mises function of their stress over
    ! the blades, the corners and the middles of the long sides of their
    ! sections and the revolutions
    if (trim(blade) == 'beam') then
       mass     = integral(4.0_dp * density * chord * h, over=omega)
       stresses = [(beams(k) % stress(yield), k = 1, 4)]
       ks       = maximum(continuous_field(omega, stresses), over=omega, weight=ks_weight)
    end if

    ! the march: bdf of the setup's order, from one step of bdf(1)
    if (order == 1) then
       march = chain([dirk(2), (bdf(1), k = 2, num_instants)], from=[(k, k = 1, num_instants)])
    else
       march = chain([dirk(2), bdf(1), (bdf(order), k = 3, num_instants)], from=[(k, k = 1, num_instants)])
    end if

    ! the laws on the revolution, the initial data on {0}, and with the
    ! blades as beams their conditions at the root, the derivatives
    ! along the span by the hermite cubics. With strips, every condition
    ! paired with its multiplier: each thickness at the setup's and the
    ! smooth maximum the objective, so that -kappa is its derivative
    ! along each thickness by the reverse sweep; or, to optimize, the
    ! least mass over the revolutions whose smooth maximum of the von
    ! Mises function is at most ks_bound^2, the thickness of each strip
    ! within the bounds and adjacent strips within the smoothness, no
    ! condition fixing the thicknesses
    r = continuous_lagrangian(omega,  rotor % laws())
    g = continuous_lagrangian(domega, rotor % initial)
    if (trim(blade) == 'rigid') then
       L   = r + g
       L_h = L % discretize(omega_h, time=march)
    else
       b = continuous_lagrangian(left, rotor % root)
       if (strips > 0) then
          thicknesses = [(omega % coordinate(numbered('h', k)), k = 1, strips)]
          d = continuous_lagrangian(omega % design(), [(thicknesses(k) - thickness(k), k = 1, strips)])
       end if
       if (strips > 0 .and. optimize .and. .not. marched) then
          allocate(limits(size(stresses)))
          do k = 1, size(stresses)
             limits(k) = stresses(k) - ks_bound**2
          end do
          bounded  = maximum(continuous_field(omega, limits), over=omega, weight=ks_weight)
          stressed = continuous_lagrangian(omega % design(), [bounded], inequality=.true.)
          adjoint  = r % multiplier('adjoint')
          lambda   = g % multiplier('lambda')
          nu       = b % multiplier('nu')
          mu_ks    = stressed % multiplier('mu_ks')
          L = adjoint*r + lambda*g + nu*b + mu_ks*stressed
          if (first) then
             within    = continuous_lagrangian(omega % design(), [(bounds(1) - thicknesses(k), k = 1, strips), &
                  &                                             (thicknesses(k) - bounds(2), k = 1, strips)], &
                  &                            inequality=.true.)
             smooth    = continuous_lagrangian(omega % design(), [(thicknesses(k) - thicknesses(k + 1) - smoothness, &
                  &                                              k = 1, strips - 1), (thicknesses(k + 1) &
                  &                                              - thicknesses(k) - smoothness, k = 1, strips - 1)], &
                  &                            inequality=.true.)
             mu_bounds = within % multiplier('mu_bounds')
             mu_smooth = smooth % multiplier('mu_smooth')
             L = mass + L + mu_bounds*within + mu_smooth*smooth
          end if
       else if (strips > 0 .and. sensitivities .and. .not. marched) then
          adjoint = r % multiplier('adjoint')
          lambda  = g % multiplier('lambda')
          nu      = b % multiplier('nu')
          kappa   = d % multiplier('kappa')
          L = ks + adjoint*r + lambda*g + nu*b + kappa*d
       else if (strips > 0) then
          L = r + g + b + d
       else
          L = r + g + b
       end if
       L_h = L % discretize(omega_h, time=march, space=finite_element(degree=3))
    end if

    ! THE LAGRANGIAN, with lagrangian: the objective and every law with
    ! its multiplier where it is paired, one equation a line, written
    ! by image 1 to rotor_lagrangian.txt, the scenario's name between
    ! where several are listed
    if (lagrangian .and. .not. marched .and. this_image() == 1) then
       file = 'rotor_lagrangian'
       if (ns > 1) file = file // '_' // trim(name)
       open(newunit=unit, file=file // '.txt', status='replace', action='write')
       write(unit, '(a)') formula(L)
       close(unit)
    end if

    motion   = rotor % estimate(omega)
    estimate = motion % discretize(omega_h)

  end subroutine stated

  !===================================================================!
  ! ONE SCENARIO REPORTED: the energy at every instant - of the rigid
  ! bodies alone when the blades are beams - and the hub's angle, the
  ! first point of each instant; with beams the tip of the blade at
  ! zero degrees and the functions; and the time history drawn, named
  ! after the scenario where several are listed, and to optimize after
  ! the design, initial or optimized.
  !===================================================================!

  subroutine reported(name, alone, rotor, mass, ks, solution, design)

    character(len=*)      , intent(in) :: name
    logical               , intent(in) :: alone
    type(mechanism)       , intent(in) :: rotor
    type(continuous_field), intent(in) :: mass, ks
    type(discrete_field)  , intent(in) :: solution
    character(len=*)      , intent(in), optional :: design

    type(continuous_field) :: energy
    type(discrete_field)   :: energy_h, hub_h, tip_h, mass_h, ks_h, kappa_h, found
    type(mechanism)        :: drawn
    real(dp), allocatable  :: e(:), hub_eps(:), found_thickness(:)
    character(len=32), allocatable :: names(:)
    character(len=:), allocatable :: drawing
    integer :: k, n

    ! the energy at every instant - of the rigid bodies alone when the
    ! blades are beams - and the hub's angle; the first point of each
    ! instant, and with beams the tip of the blade at zero degrees
    energy   = rotor % energy
    energy_h = energy % at(solution)
    allocate(e(num_instants * np), hub_eps(num_instants * np))
    call energy_h % values(e)
    hub_h = solution % fields(['hub_eps3'])
    call hub_h % values(hub_eps)
    if (present(design)) then
       print '(a, a, a, a, a, i0, a, a, a)', 'the rotor, ', trim(blade), ' blades, ', trim(name), ', ', &
            & num_instants, ' instants, the ', design, ' design:'
    else
       print '(a, a, a, a, a, i0, a)', 'the rotor, ', trim(blade), ' blades, ', trim(name), ', ', num_instants, &
            & ' instants:'
    end if
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

       ! the functions, and along each strip's thickness the mass's and
       ! the smooth maximum's derivatives by the expansion, and the smooth
       ! maximum's by the reverse sweep
       mass_h = mass % at(solution)
       ks_h   = ks % at(solution)
       print '(a, es24.16)', '   the blades'' mass over the revolutions, kg s     ', mass_h % value(1, 1)
       print '(a, es24.16)', '   the smooth maximum of their von Mises function   ', ks_h % value(1, 1)
       if (strips > 0 .and. optimize) then
          if (present(design)) then
             if (design == 'optimized') then
                names = [(numbered('h', k), k = 1, strips)]
                found = solution % fields(names)
                found_thickness = found % value(1, :)
                print '(a)', '   the thickness found, strip by strip, and the multiplier of the stress'
                print '(8es16.8)', found_thickness
                found = solution % fields(['mu_ks'])
                print '(a, es24.16)', '   mu_ks ', found % value(1, 1)
             end if
          end if
       else if (strips > 0 .and. sensitivities) then
          kappa_h = solution % fields(['kappa'])
          print '(a)', '   along the thickness of each strip: the mass''s, the smooth maximum''s by the expansion and &
               &by the reverse sweep'
          do k = 1, strips
             print '(i8, 3es24.16)', k, mass_h % value(1, 1 + k), ks_h % value(1, 1 + k), -kappa_h % value(1, k)
          end do
       end if
    end if

    ! the time history; of the design found, each blade drawn with the
    ! thickness found of each strip
    if (len_trim(meshes) > 0) then
       drawing = 'rotor'
       if (.not. alone) drawing = drawing // '_' // trim(name)
       if (present(design)) drawing = drawing // '_' // design
       drawn = rotor
       if (allocated(found_thickness)) then
          do k = 1, size(drawn % beam_surfaces)
             drawn % beam_surfaces(k) = plate_surface(k - 1, found_thickness)
          end do
       end if
       call drawn % draw(solution, drawing, every)
    end if

  end subroutine reported

  ! the surface of a body's mesh
  function mesh_surface(file) result(s)
    character(len=*), intent(in) :: file
    type(surface) :: s
    call bdf_surface(file, s % points, s % corners, s % count, s % kinds)
  end function mesh_surface

  ! the surface of the plate of the blade of k quarter turns: the box
  ! of its span, chord and thickness, 48 cells along the span and 4
  ! along the chord on each face, its thickness at each point that of
  ! the strip of the span it lies in, h(1) to h(size(h)) from the root
  function plate_surface(k, h) result(s)
    integer , intent(in) :: k
    real(dp), intent(in) :: h(:)
    type(surface) :: s
    integer, parameter :: ns = 48, nc = 4
    real(dp) :: x(3, 8), corner(3), y(3)
    integer :: i, j, f, face(4, 6), n1, n2, p, c, m
    ! the corners of the box along the azimuth zero, of unit thickness,
    ! and its faces
    do i = 1, 8
       corner = [merge(tip, root, btest(i - 1, 0)), merge(chord, -chord, btest(i - 1, 1)) / 2.0_dp, &
            & height + merge(1.0_dp, -1.0_dp, btest(i - 1, 2)) / 2.0_dp]
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
             y = x(:, face(1, f)) + (x(:, face(2, f)) - x(:, face(1, f))) * real(i, dp) / n1 &
                  & + (x(:, face(4, f)) - x(:, face(1, f))) * real(j, dp) / n2
             m    = max(1, min(size(h), 1 + int((y(1) - root) / (tip - root) * size(h))))
             y(3) = height + (y(3) - height) * h(m)
             s % points = reshape([s % points, turned(y, k)], [3, size(s % points, 2) + 1])
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
    mass    = density * thickness(1) * span * chord
    inertia = 0.0_dp
    inertia(1, 1) = mass * (chord**2 + thickness(1)**2) / 12.0_dp
    inertia(2, 2) = mass * (span**2 + thickness(1)**2) / 12.0_dp
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

  ! a name and a number: the design coordinate of strip k
  function numbered(prefix, k) result(name)
    character(len=*), intent(in) :: prefix
    integer         , intent(in) :: k
    character(len=8) :: name
    write(name, '(a, i0)') prefix, k
  end function numbered

end program rotorcraft_dynamics
