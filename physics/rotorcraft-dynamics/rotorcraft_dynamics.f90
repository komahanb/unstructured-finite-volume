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
! plate and solid are not yet stated. The section of a beam is chosen
! at setup too: the plate's, a uniform solid rectangle; or given by
! its properties per length, as TACS's TimoshenkoStiffness takes them,
! by default those of TACS's uniform HART-II case (its examples,
! rotating_beam): 0.95 kg/m, the moments of inertia 7.1828e-4,
! 1.7e-5 and 7.0128e-4 kg m, EA = 1.17e7 N, EI_f = 250 and EI_c =
! 5200 N m^2, GJ = 160 N m^2. A section given so states no stress.
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
!      section_kind of a blade as a beam: solid, the plate of the
!                   setup's thickness; or given, by section_mass,
!                   section_moments (about the span, the chord and the
!                   thickness), section_axial, section_flap (EI_f),
!                   section_lag (EI_c) and section_torsion (GJ), per
!                   length, constant along the span; the smooth
!                   maximum, the strips and the optimization ask for
!                   the solid's stress
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
!      blade_meshes the directory of the blades' surface meshes,
!                   blade0.bdf to blade270.bdf - the HART II blades of
!                   the geometry's hart2 preset - drawn in place of the
!                   boxes of their plates
!      every        the stride of the instants drawn
!      span_points  of a blade as a beam
!      precone      of the blades, degrees: each span inclined upward
!                   about the chord at the root; 0 by default, 2.5 HART
!                   II's
!      twist        of the blades as beams, degrees over the blade's
!                   radius, linear, nose up positive; 0 by default, -8
!                   HART II's
!      twist_zero   the radius of zero twist, 1.5 m, HART II's
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
!      objective    the functional of the Lagrangian - mass, stress,
!                   strain_energy, flap, lag, twist, pitch_link,
!                   push_rod, hub_force or hub_moment: to optimize,
!                   the one minimized, the mass unless named; with
!                   strips and sensitivities, the one the reverse
!                   sweep differentiates too, the stress unless named
!      constraints  to optimize, the functionals bounded in every
!                   scenario, the stress by default
!      limits       their limits: each functional at most its limit,
!                   the stress's ks_bound^2 unless given
!      bounds       the least and the greatest thickness, 0.01 and 0.02
!      smoothness   the largest difference of adjacent strips, 1e-3
!      order        of the backward differences of the march, 1 or 2,
!                   from one step of bdf(1)
!
! THE FUNCTIONS, of the blades as beams: the blades' mass over the
! revolutions, and with the solid section the smooth maximum
! (operation_field, maximum) of
! the von Mises function of their stress over the four blades, the
! corners and the middles of the long sides of their sections and the
! revolutions; and the further objectives (pieces_of), the means
! (operation_field, mean) over the revolutions, the span and the
! components of their fields: the blades' strain energy per length;
! the squares of their flap and lag displacements and of their elastic
! twist; the squares of the forces of the four pitch links on the
! upper swashplate, of the three push rods' drivers and of the hub's
! driver, and of the reactions of the hub driver's conditions on its
! tilt about x and y. A constraint is the same reduction of the field
! less its limit, at most zero. With strips, the derivatives
! of every one along the thickness of each strip by the expansion of
! the motion, and the setup's objective's by the reverse sweep too,
! -kappa, the smooth maximum unless another is named. At 12 steps, 5
! points of the span and two strips of 2 and 1.5 cm, the reverse
! sweep of each objective in turn gives its derivatives as the
! expansion does to 2.5e-14 to 9.4e-13, the mass's exactly, and the
! hub driver's tilt reactions' to 1.5e-10, their mean square some
! 1e-6 of the hub force's; central differences of step 1e-5 miss
! them by their own error, 4e-9 to 1.5e-6. At 36 steps, 5 points of
! the span and eight strips at 2 cm the mass is 2.1737794328871449 kg
! s, TACS's 2.173779432887122, and the smooth maximum at ks_weight =
! 100 1.6760974217836406, the same on 1 and 2 images; along the root
! strip's thickness its derivative by the expansion is
! -166.45582953687486, by the reverse sweep the same to 5.4e-14 (along
! the eight strips to 4e-14 to 2.5e-13), and central differences of
! steps 5e-6, 1e-5 and 2e-5 miss it by their own error, 3.29e-5,
! 1.317e-4 and 5.270e-4, extrapolated in the step's second and fourth
! powers to within 1.4e-9.
!
! THE LEAST MASS (optimize), at 36 steps, 5 points of the span and
! eight strips from 2 cm, with ks_bound = 1.25 - the stress of this
! march is 1.29 times the yield at 2 cm, and no uniform blade meets
! 0.25 at any march of the note - is found in six iterations: the root
! strip 1.0203272556788471 cm and the seven others at their bound of
! 1 cm, the mass 1.0896514020886789 kg s, the smooth maximum on its
! bound to 5e-15 and its multiplier 0.037891294249. The three
! scenarios together with ks_bound = 1.5 - their smooth maxima at 2 cm
! 1.676, 2.106 and 2.991, at 1 cm 1.638, 1.811 and 2.323 - in six
! iterations: the root strip 1.0308636977683841 cm and the seven others
! at 1 cm, the mass 1.0910828959080567 kg s, the lateral's smooth
! maximum on its bound to 7e-15 and its multiplier 0.0578856747211,
! the others' 1.525 and 1.739 and zero; the lateral's derivative along
! the root strip by central differences of step 1e-6, -234.706109775,
! gives the multiplier as the mass's over it to 5e-9.
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
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program rotorcraft_dynamics

  use util_precision    , only : dp
  use util_verbosity    , only : set_verbosity
  use operation_manifold, only : continuous_manifold, discrete_manifold, interval, region, points, parameter
  use operation_expansion, only : expansion
  use operation_field   , only : continuous_field, discrete_field, integral, mean, maximum, heaviside, &
       &                         operator(+), operator(-), operator(*)
  use operation_residual, only : continuous_residual, discrete_residual, operator(+), operator(*), &
       &                         minimize_together
  use operation_family  , only : dirk, bdf, chain
  use operation_finite_element, only : finite_element
  use multibody_dynamics, only : rigid_body, beam_blade, section, solid_rectangle, mechanism, body_table, surface, &
       &                         spherical, revolute, cylindrical, fixed, rigid_link, revolute_driver, motion_driver
  use view_bdf_loader   , only : bdf_surface

  implicit none

  real(dp), parameter :: pi = acos(-1.0_dp), speed = 109.12_dp
  real(dp), parameter :: ex(3) = [1.0_dp, 0.0_dp, 0.0_dp], ey(3) = [0.0_dp, 1.0_dp, 0.0_dp], &
       &                 ez(3) = [0.0_dp, 0.0_dp, 1.0_dp], origin(3) = 0.0_dp
  real(dp), parameter :: density = 2500.0_dp, root = 0.44_dp, tip = 2.0_dp, chord = 0.121_dp, height = 0.9_dp
  real(dp), parameter :: modulus = 70.0e9_dp, poisson = 0.3_dp

  ! the setup
  character(len=16)  :: blade = 'rigid', scenario(3) = [character(len=16) :: 'collective', '', '']
  character(len=256) :: meshes = '', blade_meshes = ''
  real(dp)           :: thickness(48) = 0.02_dp, yield = 350.0e6_dp, ks_weight = 100.0_dp
  ! the precone of the blades, degrees: each span inclined upward about
  ! the chord at the root
  real(dp)           :: precone = 0.0_dp
  ! the twist of the blades as beams: degrees over the blade's radius,
  ! linear, nose up positive, zero at the radius twist_zero
  real(dp)           :: twist = 0.0_dp, twist_zero = 1.5_dp
  real(dp)           :: ks_bound = 0.25_dp, bounds(2) = [0.01_dp, 0.02_dp], smoothness = 1.0e-3_dp
  integer            :: revolutions = 1, steps = 360, every = 10, span_points = 9, strips = 0, order = 1
  logical            :: sensitivities = .true., optimize = .false.
  ! the section of a blade as a beam: the solid plate's, or given by its
  ! properties per length, TACS's uniform HART-II case by default
  character(len=16)  :: section_kind = 'solid'
  ! the objective of the Lagrangian: to optimize, the functional
  ! minimized, the mass unless named; with sensitivities, the functional
  ! the reverse sweep differentiates too, the stress unless named. To
  ! optimize, the functionals bounded in every scenario, and their
  ! limits, the stress's ks_bound^2 unless given
  character(len=16)  :: objective = ''
  character(len=16)  :: constraints(10) = [character(len=16) :: 'stress', '', '', '', '', '', '', '', '', '']
  real(dp)           :: limits(10) = -huge(1.0_dp)
  real(dp)           :: section_mass = 0.95_dp, section_moments(3) = [7.1828e-4_dp, 1.7e-5_dp, 7.0128e-4_dp]
  real(dp)           :: section_axial = 1.17e7_dp, section_flap = 2.5e2_dp, section_lag = 5.2e3_dp
  real(dp)           :: section_torsion = 1.6e2_dp
  namelist /setup/ blade, scenario, thickness, revolutions, steps, meshes, every, span_points, strips, yield, &
       & ks_weight, sensitivities, optimize, ks_bound, bounds, smoothness, order, section_kind, section_mass, &
       & section_moments, section_axial, section_flap, section_lag, section_torsion, objective, constraints, limits, &
       & precone, twist, twist_zero, blade_meshes
  ! the objectives besides the mass and the smooth maximum, of the blades
  ! as beams (objectives_of)
  character(len=16), parameter :: further(8) = [character(len=16) :: 'strain_energy', 'flap', 'lag', 'twist', &
       & 'pitch_link', 'push_rod', 'hub_force', 'hub_moment']

  type(discrete_manifold) :: omega_h
  type(discrete_residual) :: L_h(size(scenario))
  type(discrete_field)    :: estimates(size(scenario)), solutions(size(scenario))
  type(mechanism)         :: rotors(size(scenario))
  type(continuous_field)  :: masses(size(scenario)), kss(size(scenario)), others(size(further), size(scenario))
  integer :: unit, j, ns, num_instants, np
  logical :: optimizing

  open(newunit=unit, file='rotorcraft.nml', status='old', action='read')
  read(unit, nml=setup)
  close(unit)

  ! the scenarios listed, each on a manifold of its own
  ns = count(len_trim(scenario) > 0)
  if (ns == 0) error stop 'rotorcraft_dynamics: the setup lists a scenario at least'
  if (trim(section_kind) /= 'solid' .and. trim(section_kind) /= 'given') then
     error stop 'rotorcraft_dynamics: the section of a blade as a beam is solid or given'
  end if
  if (len_trim(objective) == 0) objective = merge('mass  ', 'stress', optimize)
  if (.not. functional_named(objective)) then
     error stop 'rotorcraft_dynamics: the objective is mass, stress, strain_energy, flap, lag, twist, pitch_link, &
          &push_rod, hub_force or hub_moment'
  end if
  do j = 1, count(len_trim(constraints) > 0)
     if (.not. functional_named(constraints(j))) then
        error stop 'rotorcraft_dynamics: a constraint is one of the functionals an objective may be'
     end if
     if (limits(j) == -huge(1.0_dp)) then
        if (trim(constraints(j)) /= 'stress') error stop 'rotorcraft_dynamics: each constraint but the stress takes its limit'
        limits(j) = ks_bound**2
     end if
  end do
  if (trim(blade) == 'beam' .and. trim(section_kind) == 'given' .and. strips > 0) then
     error stop 'rotorcraft_dynamics: strips design the thickness of the solid section; a section given by its &
          &properties has none, and states no stress'
  end if
  num_instants = revolutions * steps + 1
  np = 1
  if (trim(blade) == 'beam') np = span_points
  optimizing = strips > 0 .and. optimize .and. trim(blade) == 'beam'
  if (this_image() == 1 .and. num_instants <= 11) call set_verbosity(1)

  ! to optimize with the bodies drawn, the initial design marched,
  ! reported and drawn first, beside the design found
  if (optimizing .and. len_trim(meshes) > 0) then
     do j = 1, ns
        call stated(scenario(j), j == 1, L_h(j), estimates(j), rotors(j), masses(j), kss(j), others(:, j), omega_h, &
             & initial=.true.)
        call L_h(j) % minimize(estimates(j), solutions(j), continued=.true.)
        call reported(scenario(j), ns == 1, rotors(j), masses(j), kss(j), others(:, j), solutions(j), 'initial')
     end do
  end if
  do j = 1, ns
     call stated(scenario(j), j == 1, L_h(j), estimates(j), rotors(j), masses(j), kss(j), others(:, j), omega_h)
  end do

  ! each marched in turn; to optimize, the one design of them all
  if (optimizing) then
     call minimize_together(L_h(1:ns), estimates(1:ns), solutions(1:ns), continued=.true.)
     do j = 1, ns
        call reported(scenario(j), ns == 1, rotors(j), masses(j), kss(j), others(:, j), solutions(j), 'optimized')
     end do
  else
     do j = 1, ns
        call L_h(j) % minimize(estimates(j), solutions(j), continued=.true.)
     end do
     do j = 1, ns
        call reported(scenario(j), ns == 1, rotors(j), masses(j), kss(j), others(:, j), solutions(j))
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

  subroutine stated(name, first, L_h, estimate, rotor, mass, ks, others, omega_h, initial)

    character(len=*)       , intent(in)  :: name
    logical                , intent(in)  :: first
    type(discrete_residual), intent(out) :: L_h
    type(discrete_field)   , intent(out) :: estimate
    type(mechanism)        , intent(out) :: rotor
    type(continuous_field) , intent(out) :: mass, ks, others(:)
    type(discrete_manifold), intent(out) :: omega_h
    logical                , intent(in), optional :: initial

    type(continuous_manifold) :: omega, domega, left
    type(continuous_field)    :: t, x, motion, h, adjoint, lambda, nu, kappa, mu_ks, mu_bounds, mu_smooth, bounded
    type(continuous_field), allocatable :: stresses(:), thicknesses(:), shifted(:)
    class(section)        , allocatable :: s
    type(continuous_field) :: moments(3), pieces(4, size(further)), mu
    type(surface)          :: blade_surface
    integer :: counts(size(further)), c, j
    type(continuous_residual) :: r, g, b, d, L, stressed, within, smooth
    type(chain)               :: march
    type(body_table)          :: table
    type(rigid_body), allocatable :: part(:)
    type(rigid_body)          :: blades(4)
    type(beam_blade)          :: beams(4)
    real(dp)                  :: amplitude(3)
    integer :: k, lower, upper, cap, rod
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
    select case (trim(blade))
    case ('rigid')
       omega   = continuous_manifold(time=interval(0.0_dp, revolutions * 2.0_dp * pi / speed))
       omega_h = omega % discretize(time=points(num_instants))
    case ('beam')
       if (strips > 0) then
          omega   = continuous_manifold(time=interval(0.0_dp, revolutions * 2.0_dp * pi / speed), &
               &                        space=region(interval(0.0_dp, tip - root)), &
               &                        design=[(parameter(numbered('h', k), thickness(k)), k = 1, strips)])
          ! the derivatives by the expansion of order one; to optimize, by
          ! the reverse sweeps alone, the hessian by its BFGS update
          if (sensitivities .and. .not. optimize) then
             omega_h = omega % discretize(time=points(num_instants), space=points(span_points), &
                  &                       design=expansion(order=1))
          else
             omega_h = omega % discretize(time=points(num_instants), space=points(span_points), &
                  &                       design=expansion(order=0))
          end if
       else
          omega   = continuous_manifold(time=interval(0.0_dp, revolutions * 2.0_dp * pi / speed), &
               &                        space=region(interval(0.0_dp, tip - root)))
          omega_h = omega % discretize(time=points(num_instants), space=points(span_points))
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
               & turned(coned([(root + tip) / 2.0_dp, 0.0_dp, height]), k - 1), &
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
       ! the section: the solid plate of that thickness, or given by its
       ! properties; the moments an array variable, not an array
       ! constructor, whose temporary's allocatable components the
       ! compiler in use does not release
       if (trim(section_kind) == 'solid') then
          allocate(s, source=solid_rectangle(density, modulus, poisson, chord, h))
       else
          do k = 1, 3
             moments(k) = continuous_field(omega, [section_moments(k)])
          end do
          allocate(s, source=section(continuous_field(omega, [section_mass]), moments, &
               & continuous_field(omega, [section_axial]), continuous_field(omega, [section_flap]), &
               & continuous_field(omega, [section_lag]), continuous_field(omega, [section_torsion])))
       end if
       do k = 1, 4
          if (twist /= 0.0_dp) then
             call rotor % beam(omega, beams(k), part(table % row('cap_' // trim(degrees(k - 1)))), &
                  & 'blade_' // trim(degrees(k - 1)), x, turned([root, 0.0_dp, height], k - 1), blade_axes(k - 1), s, &
                  & pitch=(twist * pi / (180.0_dp * tip)) * (x + (root - twist_zero)))
          else
             call rotor % beam(omega, beams(k), part(table % row('cap_' // trim(degrees(k - 1)))), &
                  & 'blade_' // trim(degrees(k - 1)), x, turned([root, 0.0_dp, height], k - 1), blade_axes(k - 1), s)
          end if
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
          if (len_trim(blade_meshes) > 0) then
             blade_surface = mesh_surface(trim(blade_meshes) // '/blade' // trim(degrees(k - 1)) // '.bdf')
          else if (trim(blade) == 'rigid') then
             blade_surface = plate_surface(k - 1, thickness(1:1))
          else
             blade_surface = plate_surface(k - 1, thickness(1:max(strips, 1)))
          end if
          if (trim(blade) == 'rigid') then
             call rotor % draw_with(blades(k), blade_surface)
          else
             call rotor % draw_beam_with(beams(k), part(table % row('cap_' // trim(degrees(k - 1)))), blade_surface)
          end if
       end do
    end if

    ! with the blades as beams, the blades' mass over the revolutions and,
    ! where their section states its stress, the smooth maximum of the
    ! von Mises function of their stress over the blades, the corners and
    ! the middles of the long sides of their sections and the revolutions
    if (trim(blade) == 'beam') then
       mass     = integral(4.0_dp * s % mass, over=omega)
       stresses = [(beams(k) % stress(yield), k = 1, 4)]
       if (size(stresses) > 0) ks = maximum(continuous_field(omega, stresses), over=omega, weight=ks_weight)
       call pieces_of(rotor, beams, pieces, counts)
       do j = 1, size(further)
          others(j) = mean_of(pieces(1:counts(j), j), omega, 0.0_dp)
       end do
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
    r = continuous_residual(omega,  rotor % laws())
    g = continuous_residual(domega, rotor % initial)
    if (trim(blade) == 'rigid') then
       L   = r + g
       L_h = L % discretize(omega_h, time=march)
    else
       b = continuous_residual(left, rotor % root)
       if (strips > 0) then
          thicknesses = [(omega % coordinate(numbered('h', k)), k = 1, strips)]
          d = continuous_residual(omega % design(), [(thicknesses(k) - thickness(k), k = 1, strips)])
       end if
       if (strips > 0 .and. optimize .and. .not. marched) then
          adjoint  = r % multiplier('adjoint')
          lambda   = g % multiplier('lambda')
          nu       = b % multiplier('nu')
          L = adjoint*r + lambda*g + nu*b
          ! each constraint the functional's reduction of its field less
          ! the limit, at most zero, with its multiplier; the mass's, of
          ! the design alone, stated in the first scenario
          do c = 1, count(len_trim(constraints) > 0)
             select case (trim(constraints(c)))
             case ('stress')
                allocate(shifted(size(stresses)))
                do k = 1, size(stresses)
                   shifted(k) = stresses(k) - limits(c)
                end do
                bounded  = maximum(continuous_field(omega, shifted), over=omega, weight=ks_weight)
                deallocate(shifted)
                stressed = continuous_residual(omega % design(), [bounded], inequality=.true.)
                mu_ks    = stressed % multiplier('mu_ks')
                L = L + mu_ks*stressed
             case ('mass')
                if (.not. first) cycle
                bounded  = integral(4.0_dp * s % mass - limits(c) / omega % measure(), over=omega)
                stressed = continuous_residual(omega % design(), [bounded], inequality=.true.)
                mu       = stressed % multiplier('mu_mass')
                L = L + mu*stressed
             case default
                j        = findloc(further, constraints(c), 1)
                bounded  = mean_of(pieces(1:counts(j), j), omega, limits(c))
                stressed = continuous_residual(omega % design(), [bounded], inequality=.true.)
                mu       = stressed % multiplier('mu_' // trim(constraints(c)))
                L = L + mu*stressed
             end select
          end do
          ! the objective: of the state, in every scenario; the mass, of
          ! the design alone, in the first, with the inequalities on the
          ! design alone
          select case (trim(objective))
          case ('mass')
          case ('stress')
             L = ks + L
          case default
             L = others(findloc(further, objective, 1)) + L
          end select
          if (first) then
             within    = continuous_residual(omega % design(), [(bounds(1) - thicknesses(k), k = 1, strips), &
                  &                                           (thicknesses(k) - bounds(2), k = 1, strips)], &
                  &                          inequality=.true.)
             smooth    = continuous_residual(omega % design(), [(thicknesses(k) - thicknesses(k + 1) - smoothness, &
                  &                                            k = 1, strips - 1), (thicknesses(k + 1) &
                  &                                            - thicknesses(k) - smoothness, k = 1, strips - 1)], &
                  &                          inequality=.true.)
             mu_bounds = within % multiplier('mu_bounds')
             mu_smooth = smooth % multiplier('mu_smooth')
             if (trim(objective) == 'mass') then
                L = mass + L + mu_bounds*within + mu_smooth*smooth
             else
                L = L + mu_bounds*within + mu_smooth*smooth
             end if
          end if
       else if (strips > 0 .and. sensitivities .and. .not. marched) then
          adjoint = r % multiplier('adjoint')
          lambda  = g % multiplier('lambda')
          nu      = b % multiplier('nu')
          kappa   = d % multiplier('kappa')
          select case (trim(objective))
          case ('stress')
             L = ks + adjoint*r + lambda*g + nu*b + kappa*d
          case ('mass')
             L = mass + adjoint*r + lambda*g + nu*b + kappa*d
          case default
             L = others(findloc(further, objective, 1)) + adjoint*r + lambda*g + nu*b + kappa*d
          end select
       else if (strips > 0) then
          L = r + g + b + d
       else
          L = r + g + b
       end if
       L_h = L % discretize(omega_h, time=march, space=finite_element(degree=3))
    end if

    motion   = rotor % estimate(omega)
    estimate = motion % discretize(omega_h)

  end subroutine stated

  !===================================================================!
  ! THE FURTHER OBJECTIVES of the blades as beams, each the mean over
  ! the manifold - the revolutions and the span - and over its
  ! components (operation_field, mean) of a field: the blades' strain
  ! energy per length; the squares of their flap and lag displacements
  ! and of their elastic twist; the squares of the forces of the four
  ! pitch links on the upper swashplate, of the three push rods'
  ! drivers and of the hub's driver, and of the reactions of the hub
  ! driver's conditions on its tilt about x and y. The components of
  ! objective j are pieces(1:counts(j), j), the blades', the links',
  ! the rods' or the hub's; a reaction, a function of time alone, has
  ! over the manifold its mean over the revolutions.
  !===================================================================!

  subroutine pieces_of(rotor, beams, pieces, counts)

    type(mechanism)       , intent(in)  :: rotor
    type(beam_blade)      , intent(in)  :: beams(4)
    type(continuous_field), intent(out) :: pieces(:,:)
    integer               , intent(out) :: counts(:)

    type(continuous_field), allocatable :: r(:)
    integer :: k

    counts = [4, 4, 4, 4, 4, 3, 1, 1]
    do k = 1, 4
       pieces(k, 1) = beams(k) % strain_energy()
       pieces(k, 2) = beams(k) % w * beams(k) % w
       pieces(k, 3) = beams(k) % v * beams(k) % v
       pieces(k, 4) = beams(k) % phi * beams(k) % phi
       r = rotor % reactions('swashplate_link_' // trim(degrees(k - 1)))
       pieces(k, 5) = r(1) * r(1) + r(2) * r(2) + r(3) * r(3)
    end do
    do k = 1, 3
       r = rotor % reactions('rod_driver_' // trim(degrees(k)))
       pieces(k, 6) = r(1) * r(1) + r(2) * r(2) + r(3) * r(3)
    end do
    r = rotor % reactions('hub_driver')
    pieces(1, 7) = r(1) * r(1) + r(2) * r(2) + r(3) * r(3)
    pieces(1, 8) = r(4) * r(4) + r(5) * r(5)

  end subroutine pieces_of

  ! the mean over the manifold and the components of the field of these
  ! components, each less the shift: the functional less the shift
  function mean_of(pieces, omega, shift) result(f)

    type(continuous_field)   , intent(in) :: pieces(:)
    type(continuous_manifold), intent(in) :: omega
    real(dp)                 , intent(in) :: shift
    type(continuous_field) :: f

    type(continuous_field), allocatable :: shifted(:)
    integer :: k

    if (shift == 0.0_dp) then
       f = mean(continuous_field(omega, pieces), over=omega)
       return
    end if
    allocate(shifted(size(pieces)))
    do k = 1, size(pieces)
       shifted(k) = pieces(k) - shift
    end do
    f = mean(continuous_field(omega, shifted), over=omega)

  end function mean_of

  ! whether the name is one of the functionals an objective may be
  logical function functional_named(name)
    character(len=*), intent(in) :: name
    functional_named = trim(name) == 'mass' .or. trim(name) == 'stress' .or. any(further == name)
  end function functional_named

  !===================================================================!
  ! ONE SCENARIO REPORTED: the energy at every instant - of the rigid
  ! bodies alone when the blades are beams - and the hub's angle, the
  ! first point of each instant; with beams the tip of the blade at
  ! zero degrees and the functions; and the time history drawn, named
  ! after the scenario where several are listed, and to optimize after
  ! the design, initial or optimized.
  !===================================================================!

  subroutine reported(name, alone, rotor, mass, ks, others, solution, design)

    character(len=*)      , intent(in) :: name
    logical               , intent(in) :: alone
    type(mechanism)       , intent(in) :: rotor
    type(continuous_field), intent(in) :: mass, ks, others(:)
    type(discrete_field)  , intent(in) :: solution
    character(len=*)      , intent(in), optional :: design

    type(continuous_field) :: energy
    type(discrete_field)   :: energy_h, hub_h, tip_h, mass_h, ks_h, kappa_h, found, other_h
    type(mechanism)        :: drawn
    real(dp), allocatable  :: e(:), hub_eps(:), found_thickness(:)
    character(len=32), allocatable :: names(:)
    character(len=:), allocatable :: drawing
    integer :: k, n, j, c

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
       print '(a, es24.16)', '   the blades'' mass over the revolutions, kg s     ', mass_h % value(1, 1)
       if (trim(section_kind) == 'solid') then
          ks_h = ks % at(solution)
          print '(a, es24.16)', '   the smooth maximum of their von Mises function   ', ks_h % value(1, 1)
       end if
       if (strips > 0 .and. optimize) then
          if (present(design)) then
             if (design == 'optimized') then
                names = [(numbered('h', k), k = 1, strips)]
                found = solution % fields(names)
                found_thickness = found % value(1, :)
                print '(a)', '   the thickness found, strip by strip, and the multipliers of the constraints'
                print '(8es16.8)', found_thickness
                if (any(constraints == 'stress')) then
                   found = solution % fields(['mu_ks'])
                   print '(a, es24.16)', '   mu_ks ', found % value(1, 1)
                end if
                do c = 1, count(len_trim(constraints) > 0)
                   if (trim(constraints(c)) == 'stress') cycle
                   if (trim(constraints(c)) == 'mass' .and. .not. alone .and. name /= scenario(1)) cycle
                   names = [character(len=32) :: 'mu_' // trim(constraints(c))]
                   found = solution % fields(names)
                   print '(a, a, es24.16)', '   mu_', trim(constraints(c)), found % value(1, 1)
                end do
             end if
          end if
       else if (strips > 0 .and. sensitivities) then
          kappa_h = solution % fields(['kappa'])
          if (trim(objective) == 'stress') then
             print '(a)', '   along the thickness of each strip: the mass''s, the smooth maximum''s by the expansion and &
                  &by the reverse sweep'
             do k = 1, strips
                print '(i8, 3es24.16)', k, mass_h % value(1, 1 + k), ks_h % value(1, 1 + k), -kappa_h % value(1, k)
             end do
          else
             print '(a)', '   along the thickness of each strip: the mass''s, the smooth maximum''s by the expansion'
             do k = 1, strips
                print '(i8, 2es24.16)', k, mass_h % value(1, 1 + k), ks_h % value(1, 1 + k)
             end do
          end if
       end if

       ! the further objectives, and with the strips' thicknesses design
       ! coordinates their derivatives along each by the expansion; the
       ! objective of the reverse sweep's by it too, beside them
       print '(a)', '   the further objectives, means over the revolutions:'
       do j = 1, size(further)
          other_h = others(j) % at(solution)
          print '(3x, a16, es24.16)', further(j), other_h % value(1, 1)
          if (strips > 0 .and. sensitivities .and. .not. optimize) then
             print '(5x, a, *(es24.16))', 'along the strips', other_h % value(1, 2:1 + strips)
             if (trim(objective) == trim(further(j))) then
                print '(5x, a, *(es24.16))', 'by the reverse  ', -kappa_h % value(1, 1:strips)
             end if
          end if
       end do
       if (strips > 0 .and. sensitivities .and. .not. optimize .and. trim(objective) == 'mass') then
          print '(a, *(es24.16))', '   the mass''s by the reverse sweep along the strips', -kappa_h % value(1, 1:strips)
       end if
    end if

    ! the time history; of the design found, each blade drawn with the
    ! thickness found of each strip
    if (len_trim(meshes) > 0) then
       drawing = 'rotor'
       if (.not. alone) drawing = drawing // '_' // trim(name)
       if (present(design)) drawing = drawing // '_' // design
       drawn = rotor
       if (allocated(found_thickness) .and. len_trim(blade_meshes) == 0) then
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
             s % points = reshape([s % points, turned(coned(y), k)], [3, size(s % points, 2) + 1])
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
  ! centre of mass: a plate of the span, the chord and the thickness, its
  ! moments about the blade's axes
  function plate(k) result(inertia)
    integer, intent(in) :: k
    real(dp) :: inertia(3,3), mass, span, moments(3), axes(3,3)
    integer :: i, j
    span    = tip - root
    mass    = density * thickness(1) * span * chord
    inertia = 0.0_dp
    inertia(1, 1) = mass * (chord**2 + thickness(1)**2) / 12.0_dp
    inertia(2, 2) = mass * (span**2 + thickness(1)**2) / 12.0_dp
    inertia(3, 3) = mass * (span**2 + chord**2) / 12.0_dp
    if (precone /= 0.0_dp) then
       moments = [inertia(1, 1), inertia(2, 2), inertia(3, 3)]
       axes    = blade_axes(k)
       do i = 1, 3
          do j = 1, 3
             inertia(i, j) = sum(moments * axes(i, :) * axes(j, :))
          end do
       end do
       return
    end if
    if (mod(k, 2) == 1) inertia(1:2, 1:2) = reshape([inertia(2, 2), 0.0_dp, 0.0_dp, inertia(1, 1)], [2, 2])
  end function plate

  ! the axes of the blade of k quarter turns - span, chord, thickness -
  ! the span along the azimuth inclined upward by the precone
  function blade_axes(k) result(axes)
    integer, intent(in) :: k
    real(dp) :: axes(3,3), angle
    if (precone == 0.0_dp) then
       axes = reshape([turned(ex, k), turned(ey, k), ez], [3, 3])
       return
    end if
    angle = precone * pi / 180.0_dp
    axes  = reshape([turned([cos(angle), 0.0_dp, sin(angle)], k), turned(ey, k), &
         &           turned([-sin(angle), 0.0_dp, cos(angle)], k)], [3, 3])
  end function blade_axes

  ! a point of the blade along the azimuth zero, coned about the chord at
  ! the root by the precone
  function coned(p) result(q)
    real(dp), intent(in) :: p(3)
    real(dp) :: q(3), angle
    q = p
    if (precone == 0.0_dp) return
    angle = precone * pi / 180.0_dp
    q(1) = root + (p(1) - root) * cos(angle) - (p(3) - height) * sin(angle)
    q(3) = height + (p(1) - root) * sin(angle) + (p(3) - height) * cos(angle)
  end function coned

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
