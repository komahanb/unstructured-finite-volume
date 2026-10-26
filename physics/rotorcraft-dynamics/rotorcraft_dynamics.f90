!=====================================================================!
! THE DYNAMICS OF A FOUR-BLADED ROTOR
!
! Consider changing the material thickness along a rotor blade. Its
! mass, resistance to deformation and capacity to carry load change
! together. To choose a thickness distribution, we need the blade's
! motion, the stresses it experiences whilst in motion, and their
! dependence on that distribution. We shall use the same rotor to
! examine these questions with a model chosen for each blade, then
! allow the beam blades to share their thicknesses or vary them
! separately.
!
! The rotor is the case of komahanb/rotorcraft-optimization
! (dissertation, section 5.4). With rigid blades it has 26 bodies and
! the 34 joints of collective.py, longitudinal_cyclic.py and
! lateral_cyclic.py, with a 35th joint, the rotating scissors. The hub
! turns at 109.12 rad/s about z. Three push rods move the lower
! swashplate along z by a sin(w t), w a quarter of the hub's rate.
! The upper swashplate hinges about z on the lower one; the scissors
! turns it with the hub. Through the lower and upper pitch links and
! the blade caps, this motion sets the pitch of the blades. The base
! plate is fixed, the sphere slides along the shaft's axis, and two
! horns hinge the lower swashplate to the base plate. The bodies
! begin at the velocities the drivers' rates require
! (multibody_dynamics, joint), the push rods' a w along z where
! rotor_bodies.txt states them at rest.
!
! Each blade extends along its cap's azimuth, 0, 90, 180 or 270
! degrees, with span 1.56 from radius 0.44 and chord 0.121. The
! material has density = 2500, modulus = 70e9 and poisson = 0.3.
! The setup key blades gives four models in that azimuth order:
! 'rigid' or 'beam'.
! If blades is not given, blade chooses the model of all four, as
! before. 'plate' and 'solid' are refused by error stop.
!
! To locate a beam's deformation, we need both an instant and a
! position along its span. With any beam present, omega is therefore
! time x span, with the design coordinates added when strips > 0.
! Each beam is fixed at its root in its cap's axes, its displacement
! and twist fields approximated along the span by hermite cubics.
! A rigid blade is a rigid body on the same omega, with the plate's
! mass and inertia at thickness(1). Its unknowns depend on time
! alone, and rigid_link joins it to its cap through cap_blade_<deg>,
! just as when all four blades are rigid. Only beams contribute
! span fields and root conditions. With all blades rigid, omega is
! time alone. The bodies' kinetic energy less their potential energy,
! with each joint condition paired with its reaction, supplies their
! Lagrangian density; the beams add kinetic energy less strain and
! potential energies over the span.
! Stationarity gives the motion laws, with initial data and the
! beams' root conditions completing the problem.
!
! Now let the beam thickness vary from one strip of the span to the
! next. With each_blade = .false., the design coordinates h1 through
! h<strips> give one distribution shared by every beam. With
! each_blade = .true., each beam has its own coordinates h<deg>_<k>,
! named by strip: for example, h90_1 through h270_8 for eight strips
! on the beams at 90, 180 and 270 degrees. Strip k of every beam
! begins at thickness(k). Without strips, every beam has the uniform
! thickness(1). The rigid blades keep their constant thickness and
! mass throughout the design calculation.
!
! To measure the cost of these thicknesses, mass integrates the
! beams' mass over the revolutions, in kg s; the rigid blades' mass
! is not counted. To measure their stress, ks is the smooth maximum
! (operation_field, maximum) of the von Mises function over all the
! beams, the corners and the middles of the long sides of their
! sections, and the revolutions. With each_blade, blade_ks also gives
! each beam's own smooth maximum, in the order of flexible.
!
! We can first hold the thicknesses at their setup values and ask
! how these quantities change with each design coordinate. The
! expansion of the motion gives their derivatives. For a reverse
! sweep, every law is paired with its multiplier and ks is the
! objective; the multiplier kappa of the conditions fixing the
! thicknesses gives -kappa as the derivative of ks. This objective
! remains the smooth maximum over all beams with each_blade: a sum
! of smooth maxima is not a functional the library forms. The
! derivatives of each blade_ks are read from the expansion.
!
! To find thicknesses of least mass, optimize releases those fixed
! values and makes mass the objective. With shared strips, ks must
! be at most ks_bound^2. With each_blade, each beam's blade_ks must
! be at most ks_bound^2, and mu_ks has one component per beam. Every
! strip lies within bounds; smoothness limits the difference
! between adjacent strips of the same blade. Several scenarios
! retain one set of design coordinates: each scenario has its own
! motion and stress bounds, while the first states mass and the
! inequalities on thickness alone (minimize_together).
!
! THE SETUP, rotorcraft.nml:
!
!      blade        'rigid' or 'beam' for all blades if blades is not
!                   given; 'plate' and 'solid' are refused
!      blades       four models, 'rigid' or 'beam', at 0, 90, 180 and
!                   270 degrees; omitted, every blade takes blade
!      scenario     collective, longitudinal or lateral: the push
!                   rods' amplitudes, 0.05 each, or 0.025, 0.025 and
!                   0.05, or 0.025, 0.05 and 0.025 at 90, 180 and 270
!                   degrees; several listed, each is stated on a
!                   manifold of its own and marched in turn
!      thickness    thickness(1) for a rigid blade or a beam without
!                   strips; with strips, thickness(k) initializes
!                   strip k of every beam
!      revolutions  of the hub
!      steps        per revolution
!      meshes       the directory of the bodies' surface meshes, the
!                   case's data/mesh; none, no drawing
!      every        the stride of the instants drawn
!      span_points  of a blade as a beam
!      strips       of each beam's span, each thickness a design
!                   coordinate; 0, thickness(1) throughout
!      each_blade   .false. by default: shared strips; .true., each
!                   beam's own strips and reported blade_ks, with
!                   optimize bounding each blade_ks separately
!      yield        the yield stress, 350e6
!      ks_weight    the weight of the smooth maximum of the stress
!      sensitivities with strips, the functions' derivatives along each
!                   thickness (.true.), or the functions alone
!      optimize     with strips, find the least mass under the stress,
!                   thickness and smoothness bounds (.false. by default)
!      ks_bound     the bound of the von Mises ratio, 0.25; its square
!                   bounds the smooth maximum in each scenario
!      bounds       the least and the greatest thickness, 0.01 and 0.02
!      smoothness   the largest difference of adjacent strips of the
!                   same blade, 1e-3
!      order        of the backward differences of the march, 1 or 2
!                   (2), from one step of bdf(1)
!      lagrangian   the Lagrangian of each scenario written as a
!                   formula, one equation a line, to
!                   rotor_lagrangian.txt (.true.), or not (.false.)
!
! We can check the motion first with all blades rigid. At 360 steps
! of bdf(2), from one step of bdf(1), the energy is 4.99432169e5 J
! initially, 287.675 J above that of the table's rates because of the
! motion the push rods' rates require. At the end of a revolution it
! is 4.99195025e5, the drivers' work and the march's error between;
! the hub's angle holds to 2.4e-18 at every instant printed. Without
! the scissors, as in TACS, the upper swashplate's turn is free, held
! to the hub by the pitch links alone as they tilt. Its lag after a
! revolution is the march's error in its rate: 64.4 degrees at 36
! steps of bdf(1), 10.1 at 360, 8.60 at 36 steps of bdf(2) and 2.14
! at 72, with pitch-link tilts of 43.3, 6.1, 5.09 and 0.77 degrees.
!
! With all four blades as beams, the solve eliminates the rows that
! tie copies of unknowns of time alone, then solves the remaining
! rows directly. At 30 degrees of 144 steps, the tip of the blade at
! zero degrees has moved by 7.82649e-4 along the span. The static
! stretch, the span integral of centrifugal stress
! rho w^2 (R^2 - r^2)/2 over modulus, is 7.6583e-4. Along the chord, the
! displacements are -2.13686e-3, -2.13816e-3, -2.13845e-3 and
! -2.13847e-3 at 3, 5, 9 and 17 span points; the differences fall by
! 4.5 and 15.
!
! For these four beams with shared strips, at 144 steps, 5 span
! points and eight strips at 2 cm, mass is 2.1737794328870783 kg s,
! TACS's 2.173779432887122. At ks_weight = 100, ks is
! 2.5151758681472021e-2, the same on 1 and 2 images. The stress is
! the centrifugal stress of the turning blade, 0.162 of the yield
! at the root at every uniform thickness, its square there 0.0262.
! The smooth maximum is 2.51151e-2, 2.51518e-2, 2.51569e-2 and
! 2.51574e-2 at 3, 5, 9 and 17 span points, and 3.961e-2, 2.515e-2
! and 2.504e-2 at 72, 144 and 288 steps. Its derivative along the
! root strip is -1.2065199020683000 by the expansion: thinning that
! strip leaves it carrying the outer strips' load. The reverse
! sweep agrees to 2.2e-15 along the root strip and 5.0e-15 along the
! eight strips. The derivatives sum to 1.13e-3, a uniform thickness
! leaving the stress as it is. Central differences with steps 5e-6,
! 1e-5 and 2e-5 miss it by their own error, 7.681e-7, 3.073e-6 and
! 1.229e-5, extrapolated in the step's second and fourth powers to
! within 2.0e-10.
!
! At 144 steps and 5 span points, the design question has a simple
! answer. From eight shared strips at 2 cm, with ks_bound = 0.25,
! ks is below the bound's 0.0625 at every uniform thickness from 1 to 2 cm.
! optimize finds every strip at its bound of 1 cm in two iterations:
! mass is 1.0868897164435392 kg s, mu_ks is zero and each strip's
! bound has multiplier 13.586121455544628, its mass derivative,
! 4 rho c (L/8) 2 pi/w = 13.5861. The three scenarios together give
! the same design in two iterations, with smooth maxima at 1 cm of
! 2.51338e-2, 2.51312e-2 and 2.51287e-2. At 36 steps of bdf(1)
! without the scissors, the smooth maximum was 1.676 at 2 cm, from
! the upper swashplate's lag; the least mass under the bound 1.25
! in its place left the root strip at 1.02 cm.
!
! Suppose instead that the blade at zero degrees is rigid and the
! other three are beams, each with its own strips. The following
! measurements use 12 steps of bdf(2), from one step of bdf(1),
! 5 span points and 8 strips initially at 2 cm. At 12 steps the
! march is coarse: ks is near 5 where 144 steps give 2.5e-2. These
! runs verify the statement and its derivatives, not the physics.
! At the same 12-step resolution, blade = 'beam' with blades and
! each_blade omitted prints every line as before, both with
! sensitivities and with optimize at ks_bound = 2.24; only the
! times differ.
!
! With blades = 'rigid', 'beam', 'beam', 'beam' and
! each_blade = .true., there are 24 design coordinates. The beams'
! mass is 1.6303345746653470 kg s, three quarters of the four beams'
! 2.1737794328871334 at 12 steps. Its derivative along each strip
! is 3.3965303638861 kg s/m, a quarter of the four beams' shared-strip
! derivative 13.586121455544571. Over the three beams, ks is
! 4.9951328594446212; blade_ks at 90, 180 and 270 degrees is
! 4.9841664260969605, 4.9837199257344640 and 4.9845371374955736.
! The derivatives of ks by the expansion and the reverse sweep
! agree to 1.4e-13 relative along all 24 strips. Each blade_ks has
! derivative at most 1.1e-14 in magnitude along another blade's
! strip: the hub and swashplates are driven, so one blade's thickness
! does not change another blade's motion. Along their own root
! strips the derivatives are -513.08, -513.03 and -513.11, compared
! with -513.09 for the four beams' shared root strip. The march
! with its 24 expansions and its adjoint took 4077 s and 633 MB.
!
! With optimize and ks_bound = 2.24, the same mixed rotor reaches
! every strip's bound of 1 cm in two iterations. The beams' mass is
! 8.1516728733267363e-1 kg s, each stress bound is inactive, and
! every component of mu_ks is zero. Tightening ks_bound to 2.2
! makes each stress bound active at 4.84. At the five design iterates,
! the three blade_ks values are above the bound by up to roughly
! 1.45e-1, below it by up to 8.8e-2, then above it by up to
! 1.1e-3, 1.8e-7 and 1.3e-14. Each blade_ks reaches 4.84 to
! 1.4e-14. In azimuth order, the root strips are 1.00420219e-2,
! 1.00415906e-2 and 1.00423813e-2 m; every other strip is at 1 cm.
! The components of mu_ks are 3.4085047050319067e-3,
! 3.4083674569569708e-3 and 3.4086315756461419e-3. The beams' mass
! is 8.1559522908988735e-1 kg s; this design took 4 min 26 s.
!
! To read these results, the title names the models in azimuth order
! when they differ, for example 'rigid, beam, beam, beam blades'.
! The motion report gives the rigid bodies' energy and the hub's angle.
! With beams, it adds the first beam's tip and its azimuth, and reports
! mass and ks, and with each_blade each blade_ks. The derivative
! table includes every design coordinate; with each_blade its name
! is the first column, followed by the mass derivative, the ks
! derivatives by expansion and reverse sweep, and each blade_ks
! derivative by expansion in azimuth order. The optimized report
! gives every thickness found and every component of mu_ks.
!
! The drawing lets us compare these thickness distributions with
! their motion. The table's rigid bodies use their surface meshes;
! a rigid blade uses a box at thickness(1), moved with the body.
! Each beam uses a box whose points move with their beam sections,
! with deflection stored as a vector at each point. The files are
! rotor_<instant>.vtu and rotor.pvd for ParaView, or
! rotor_<scenario>_... for several scenarios. With meshes and
! optimize, the initial design is marched, reported and drawn first
! as rotor_initial_*, then the design found as rotor_optimized_*.
! Each beam is drawn at the thicknesses found for its own strips
! with each_blade. This was also checked for the mixed rotor: both
! drawings contain the rigid blade's box and each beam's box.
! README.md gives the run and drawing instructions.
!
! With lagrangian, L is written as a formula before its march, one
! equation a line. The all-rigid rotor has 356 laws from stationarity
! of the density and 364 initial data for its 26 bodies, without
! multipliers: 720 lines of 876 kB. The four beams with eight shared
! strips, ks and every law times its multiplier give 689 lines of
! 2.0 MB. The tree of L is not written; the density's tree alone is
! 13 921 lines of 10 MB.
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
  character(len=16)  :: blade = 'rigid', blades(4) = '', scenario(3) = [character(len=16) :: 'collective', '', '']
  character(len=256) :: meshes = ''
  real(dp)           :: thickness(48) = 0.02_dp, yield = 350.0e6_dp, ks_weight = 100.0_dp
  real(dp)           :: ks_bound = 0.25_dp, bounds(2) = [0.01_dp, 0.02_dp], smoothness = 1.0e-3_dp
  integer            :: revolutions = 1, steps = 360, every = 10, span_points = 9, strips = 0, order = 2
  logical            :: sensitivities = .true., optimize = .false., lagrangian = .false., each_blade = .false.
  namelist /setup/ blade, blades, scenario, thickness, revolutions, steps, meshes, every, span_points, strips, &
       & yield, ks_weight, sensitivities, optimize, ks_bound, bounds, smoothness, order, lagrangian, each_blade

  type(discrete_manifold) :: omega_h
  type(discrete_lagrangian) :: L_h(size(scenario))
  type(discrete_field)    :: estimates(size(scenario)), solutions(size(scenario))
  type(mechanism)         :: rotors(size(scenario))
  type(continuous_field)  :: masses(size(scenario)), kss(size(scenario)), blade_kss(4, size(scenario))
  integer, allocatable :: flexible(:)
  character(len=8), allocatable :: designs(:)
  character(len=:), allocatable :: models
  integer :: unit, j, ns, num_instants, np, sets, k
  logical :: optimizing
  type(discrete_minimizer) :: M_h

  open(newunit=unit, file='rotorcraft.nml', status='old', action='read')
  read(unit, nml=setup)
  close(unit)

  ! Choose the models at 0, 90, 180 and 270 degrees. With blades
  ! omitted, blade retains its choice for all four. flexible lists the
  ! beams in azimuth order; models names the choice in the run's titles.
  if (all(len_trim(blades) == 0)) blades = blade
  do k = 1, 4
     if (trim(blades(k)) /= 'rigid' .and. trim(blades(k)) /= 'beam') &
          & error stop 'rotorcraft_dynamics: the blade model is rigid or beam; plate and solid are not yet stated'
  end do
  flexible = pack([(k, k = 1, 4)], blades == 'beam')
  models   = trim(blades(1))
  if (any(blades /= blades(1))) models = trim(blades(1)) // ', ' // trim(blades(2)) // ', ' // trim(blades(3)) &
       & // ', ' // trim(blades(4))

  ! To vary thickness, give every beam the shared h1..h<strips>, or
  ! with each_blade give each beam its own h<deg>_<k> through strip.
  ! designs orders the strips from root to tip, then the beams by
  ! azimuth. Rigid blades add no design coordinates.
  sets = 1
  if (each_blade) sets = size(flexible)
  allocate(designs(0))
  if (strips > 0 .and. size(flexible) > 0) then
     if (each_blade) then
        designs = [((strip(flexible(j), k), k = 1, strips), j = 1, sets)]
     else
        designs = [(numbered('h', k), k = 1, strips)]
     end if
  end if

  ! the scenarios listed, each on a manifold of its own
  ns = count(len_trim(scenario) > 0)
  if (ns == 0) error stop 'rotorcraft_dynamics: the setup lists a scenario at least'
  num_instants = revolutions * steps + 1
  np = 1
  if (size(flexible) > 0) np = span_points
  optimizing = strips > 0 .and. optimize .and. size(flexible) > 0

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
        call stated(scenario(j), j == 1, L_h(j), estimates(j), rotors(j), masses(j), kss(j), blade_kss(:, j), omega_h, &
             & initial=.true.)
        if (this_image() == 1) print '(a, a, a, a, a)', 'the march of the rotor, ', models, ' blades, ', &
             & trim(scenario(j)), ', the initial design'
        M_h = discrete_minimizer(L_h(j), newton(elimination(dense_direct(equilibrated=size(flexible) > 0)), &
             &                                  tolerance=1.0e-13_dp, chord=.true.))
        call L_h(j) % minimize(estimates(j), solutions(j), continued=.true., by=M_h)
        call reported(scenario(j), ns == 1, rotors(j), masses(j), kss(j), blade_kss(:, j), solutions(j), 'initial')
     end do
  end if
  do j = 1, ns
     call stated(scenario(j), j == 1, L_h(j), estimates(j), rotors(j), masses(j), kss(j), blade_kss(:, j), omega_h)
  end do

  ! the minimizer of every block: the march over the blocks, and the
  ! minimizer of a block
  M_h = discrete_minimizer(L_h(1), newton(elimination(dense_direct(equilibrated=size(flexible) > 0)), &
       &                               tolerance=1.0e-13_dp, chord=.true.))
  if (this_image() == 1) print '(a)', tree(M_h)

  ! each marched in turn; to optimize, the one design of them all
  if (optimizing) then
     if (this_image() == 1) print '(*(a))', 'the design of the rotor, ', models, ' blades, of the scenarios', &
          & (' ' // trim(scenario(j)), j = 1, ns)
     call minimize_together(L_h(1:ns), estimates(1:ns), solutions(1:ns), continued=.true., by=M_h)
     do j = 1, ns
        call reported(scenario(j), ns == 1, rotors(j), masses(j), kss(j), blade_kss(:, j), solutions(j), 'optimized')
     end do
  else
     do j = 1, ns
        if (this_image() == 1) print '(a, a, a, a)', 'the march of the rotor, ', models, ' blades, ', &
             & trim(scenario(j))
        call L_h(j) % minimize(estimates(j), solutions(j), continued=.true., by=M_h)
     end do
     do j = 1, ns
        call reported(scenario(j), ns == 1, rotors(j), masses(j), kss(j), blade_kss(:, j), solutions(j))
     end do
  end if

contains

  !===================================================================!
  ! ONE SCENARIO STATED on a manifold of its own: the revolutions, the
  ! bodies, the blades and the joints with the push rods' amplitudes of
  ! the scenario, the functions, the Lagrangian discretized and the
  ! estimate. To optimize, every scenario bounds ks, or with each_blade
  ! each blade_ks; the first states mass and the inequalities on the
  ! thicknesses alone. initial requests a march at the setup's
  ! thicknesses so that its motion can be compared with the design found.
  !===================================================================!

  subroutine stated(name, first, L_h, estimate, rotor, mass, ks, blade_ks, omega_h, initial)

    character(len=*)       , intent(in)  :: name
    logical                , intent(in)  :: first
    type(discrete_lagrangian), intent(out) :: L_h
    type(discrete_field)   , intent(out) :: estimate
    type(mechanism)        , intent(out) :: rotor
    type(continuous_field) , intent(out) :: mass, ks, blade_ks(4)
    type(discrete_manifold), intent(out) :: omega_h
    logical                , intent(in), optional :: initial

    type(continuous_manifold) :: omega, domega, left
    type(continuous_field)    :: t, x, motion, h(4), adjoint, lambda, nu, kappa, mu_ks, mu_bounds, mu_smooth, line_mass
    type(continuous_field), allocatable :: stresses(:), thicknesses(:), limits(:), bounded(:)
    type(continuous_lagrangian) :: r, g, b, d, L, stressed, within, smooth
    type(chain)               :: march
    type(body_table)          :: table
    type(rigid_body), allocatable :: part(:)
    type(rigid_body)          :: rigid(4)
    type(beam_blade)          :: beams(4)
    real(dp)                  :: amplitude(3)
    type(points)              :: instants
    character(len=:), allocatable :: file
    integer :: k, s, o, lower, upper, cap, rod, unit
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

    ! A beam's motion needs an instant and a position along its span.
    ! With any beam, omega is time x span, with design coordinates for
    ! the strips; otherwise it is time alone. The rigid bodies retain
    ! unknowns of time alone on this same manifold.
    instants = points(num_instants)

    if (size(flexible) == 0) then
       omega   = continuous_manifold(time=interval(0.0_dp, revolutions * 2.0_dp * pi / speed))
       omega_h = omega % discretize(time=instants)
    else
       if (strips > 0) then
          ! Repeat thickness(k) for strip k of each beam when its
          ! thicknesses are separate; the shared set uses it once.
          omega   = continuous_manifold(time=interval(0.0_dp, revolutions * 2.0_dp * pi / speed), &
               &                        space=region(interval(0.0_dp, tip - root)), &
               &                        design=[(parameter(designs(k), thickness(1 + mod(k - 1, strips))), &
               &                                 k = 1, size(designs))])
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
    end if
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

    ! Build the beam thickness along the span before stating its motion.
    ! h(1) is shared by every beam unless each_blade gives each its own
    ! h(s), in flexible order. Without strips it is thickness(1).
    ! Within each set, h = h_1 + sum over j > 1 of
    ! (h_j - h_(j-1)) H(x - x_j), with x_j the inner edge of strip j.
    if (size(flexible) > 0) then
       if (strips > 0) then
          do s = 1, sets
             o = (s - 1) * strips
             h(s) = omega % coordinate(designs(o + 1))
             do k = 2, strips
                h(s) = h(s) + (omega % coordinate(designs(o + k)) - omega % coordinate(designs(o + k - 1))) &
                     & * heaviside(x - real(k - 1, dp) * (tip - root) / real(strips, dp))
             end do
          end do
       else
          h(1) = continuous_field(omega, [thickness(1)])
       end if
    end if
    ! Each blade lies along its cap's azimuth, with chord across it and
    ! thickness along z. A rigid blade keeps thickness(1), its plate
    ! mass and inertia, and unknowns of time alone. Only a beam adds
    ! fields of the span and conditions fixing its root to its cap.
    s = 0
    do k = 1, 4
       select case (trim(blades(k)))
       case ('rigid')
          call rotor % body(omega, rigid(k), 'blade_' // trim(degrees(k - 1)), &
               & turned([(root + tip) / 2.0_dp, 0.0_dp, height], k - 1), &
               & density * thickness(1) * (tip - root) * chord, plate(k - 1), angular_velocity=speed * ez)
       case ('beam')
          if (each_blade .and. strips > 0) s = s + 1
          call rotor % beam(omega, beams(k), part(table % row('cap_' // trim(degrees(k - 1)))), &
               & 'blade_' // trim(degrees(k - 1)), x, turned([root, 0.0_dp, height], k - 1), &
               & reshape([turned(ex, k - 1), turned(ey, k - 1), ez], [3, 3]), density, modulus, poisson, chord, &
               & h(max(s, 1)))
       end select
    end do

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

      ! At each azimuth, the pitch links connect the swashplate, shaft
      ! and cap. Every rigid blade is attached by cap_blade_<deg> using
      ! rigid_link, also when beams share omega; a beam's attachment is
      ! already stated by its root conditions.
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
         if (trim(blades(k)) == 'rigid') call rotor % joint(omega, 'cap_blade_' // trim(degrees(k - 1)), &
              & rigid_link(part(cap), rigid(k)))
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

    ! Draw the table's bodies with their meshes and each blade with a
    ! box: thickness(1) for a rigid blade, the initial strip thicknesses
    ! for a beam. reported replaces the beam surfaces at the design found.
    if (len_trim(meshes) > 0) then
       do k = 1, size(table % name)
          call rotor % draw_with(part(k), mesh_surface(trim(meshes) // '/' // trim(table % mesh(k))))
       end do
       do k = 1, 4
          if (trim(blades(k)) == 'rigid') then
             call rotor % draw_with(rigid(k), plate_surface(k - 1, thickness(1:1)))
          else
             call rotor % draw_beam_with(beams(k), part(table % row('cap_' // trim(degrees(k - 1)))), &
                  & plate_surface(k - 1, thickness(1:max(strips, 1))))
          end if
       end do
    end if

    ! To compare thickness distributions, mass integrates the beams'
    ! mass over the revolutions. A rigid blade's constant mass is excluded.
    ! ks measures the smooth maximum of the beams' von Mises function
    ! over their sections, span and revolutions. Each beam contributes
    ! six stress fields, at the corners and middles of the long sides;
    ! with each_blade, blade_ks(s) reduces those of flexible(s) alone.
    if (size(flexible) > 0) then
       if (each_blade .and. strips > 0) then
          line_mass = density * chord * h(1)
          do s = 2, sets
             line_mass = line_mass + density * chord * h(s)
          end do
          mass = integral(line_mass, over=omega)
       else
          mass = integral(real(size(flexible), dp) * density * chord * h(1), over=omega)
       end if
       stresses = [(beams(flexible(s)) % stress(yield), s = 1, size(flexible))]
       ks       = maximum(continuous_field(omega, stresses), over=omega, weight=ks_weight)
       if (each_blade) then
          do s = 1, size(flexible)
             blade_ks(s) = maximum(continuous_field(omega, stresses(6 * s - 5:6 * s)), over=omega, weight=ks_weight)
          end do
       end if
    end if

    ! the march: bdf of the setup's order, from one step of bdf(1)
    if (order == 1) then
       march = chain([dirk(2), (bdf(1), k = 2, num_instants)], from=[(k, k = 1, num_instants)])
    else
       march = chain([dirk(2), bdf(1), (bdf(order), k = 3, num_instants)], from=[(k, k = 1, num_instants)])
    end if

    ! The motion laws r and initial data g apply to the whole mechanism;
    ! b contains only the beams' root conditions. For sensitivities, d
    ! holds every strip at its setup thickness. Pairing these laws with
    ! their multipliers and taking ks as objective gives -kappa as its
    ! derivative along each design coordinate. To optimize, release d
    ! and minimize mass under stress, thickness and smoothness bounds.
    r = continuous_lagrangian(omega,  rotor % laws())
    g = continuous_lagrangian(domega, rotor % initial)
    if (size(flexible) == 0) then
       L   = r + g
       L_h = L % discretize(omega_h, time=march)
    else
       b = continuous_lagrangian(left, rotor % root)
       if (strips > 0) then
          thicknesses = [(omega % coordinate(designs(k)), k = 1, size(designs))]
          d = continuous_lagrangian(omega % design(), [(thicknesses(k) - thickness(1 + mod(k - 1, strips)), &
               &                                        k = 1, size(designs))])
       end if
       if (strips > 0 .and. optimize .and. .not. marched) then
          allocate(limits(size(stresses)))
          do k = 1, size(stresses)
             limits(k) = stresses(k) - ks_bound**2
          end do
          if (each_blade) then
             ! One stress inequality per beam: blade_ks(s) <= ks_bound^2.
             ! mu_ks consequently has one component for each flexible(s).
             bounded = [(maximum(continuous_field(omega, limits(6 * s - 5:6 * s)), over=omega, weight=ks_weight), &
                  &      s = 1, size(flexible))]
          else
             bounded = [maximum(continuous_field(omega, limits), over=omega, weight=ks_weight)]
          end if
          stressed = continuous_lagrangian(omega % design(), bounded, inequality=.true.)
          adjoint  = r % multiplier('adjoint')
          lambda   = g % multiplier('lambda')
          nu       = b % multiplier('nu')
          mu_ks    = stressed % multiplier('mu_ks')
          L = adjoint*r + lambda*g + nu*b + mu_ks*stressed
          if (first) then
             ! Bound every strip, but compare neighbours only within
             ! its own set. No smoothness condition joins two blades.
             within    = continuous_lagrangian(omega % design(), [(bounds(1) - thicknesses(k), k = 1, size(designs)), &
                  &                                             (thicknesses(k) - bounds(2), k = 1, size(designs))], &
                  &                            inequality=.true.)
             smooth    = continuous_lagrangian(omega % design(), &
                  & [((thicknesses(o + k) - thicknesses(o + k + 1) - smoothness, k = 1, strips - 1), &
                  &   o = 0, size(designs) - strips, strips), &
                  &  ((thicknesses(o + k + 1) - thicknesses(o + k) - smoothness, k = 1, strips - 1), &
                  &   o = 0, size(designs) - strips, strips)], inequality=.true.)
             mu_bounds = within % multiplier('mu_bounds')
             mu_smooth = smooth % multiplier('mu_smooth')
             L = mass + L + mu_bounds*within + mu_smooth*smooth
          end if
       else if (strips > 0 .and. sensitivities .and. .not. marched) then
          ! Keep ks over all beams as the one sensitivity objective.
          ! The library does not form a sum of smooth maxima as a
          ! functional; each blade_ks derivative comes from the expansion.
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
  ! ONE SCENARIO REPORTED: the models, the rigid bodies' energy and
  ! hub's angle, then the first beam's tip identified by azimuth and
  ! the functions. With each_blade, report each blade_ks and name
  ! every design coordinate in the derivative table. Draw the
  ! time history with the scenario's name where several are listed,
  ! and with initial or optimized when comparing designs.
  !===================================================================!

  subroutine reported(name, alone, rotor, mass, ks, blade_ks, solution, design)

    character(len=*)      , intent(in) :: name
    logical               , intent(in) :: alone
    type(mechanism)       , intent(in) :: rotor
    type(continuous_field), intent(in) :: mass, ks, blade_ks(4)
    type(discrete_field)  , intent(in) :: solution
    character(len=*)      , intent(in), optional :: design

    type(continuous_field) :: energy
    type(discrete_field)   :: energy_h, hub_h, tip_h, mass_h, ks_h, kappa_h, found
    type(discrete_field), allocatable :: blade_ks_h(:)
    type(mechanism)        :: drawn
    real(dp), allocatable  :: e(:), hub_eps(:), found_thickness(:)
    character(len=32), allocatable :: names(:)
    character(len=:), allocatable :: drawing, tip
    integer :: k, n, s, o

    ! The rigid bodies' energy includes rigid blades; it excludes beams.
    ! Read it and the hub's angle at the first span point of each
    ! instant. flexible(1) chooses the first beam's tip, whose azimuth
    ! is printed whether or not each_blade is set.
    energy   = rotor % energy
    energy_h = energy % at(solution)
    allocate(e(num_instants * np), hub_eps(num_instants * np))
    call energy_h % values(e)
    hub_h = solution % fields(['hub_eps3'])
    call hub_h % values(hub_eps)
    if (present(design)) then
       print '(a, a, a, a, a, i0, a, a, a)', 'the rotor, ', models, ' blades, ', trim(name), ', ', &
            & num_instants, ' instants, the ', design, ' design:'
    else
       print '(a, a, a, a, a, i0, a)', 'the rotor, ', models, ' blades, ', trim(name), ', ', num_instants, &
            & ' instants:'
    end if
    n = max(1, (num_instants - 1) / (12 * revolutions))
    if (size(flexible) == 0) then
       print '(a)', '   angle         energy   hub eps3 - sin(w t/2)'
       do k = 1, num_instants, n
          print '(f8.1, es16.8, es14.3)', speed * omega_h % instant(k) * 180.0_dp / pi, e(k), &
               & hub_eps(k) - sin(speed * omega_h % instant(k) / 2.0_dp)
       end do
    else
       tip   = 'blade_' // trim(degrees(flexible(1) - 1))
       names = [character(len=32) :: tip // '_u', tip // '_v', tip // '_w', tip // '_phi']
       tip_h = solution % fields(names)
       print '(a, a, a)', '   angle   bodies'' energy   hub eps3 - sin(w t/2)    the tip of the blade at ', &
            & trim(degrees(flexible(1) - 1)), ': u, v, w, phi'
       do k = 1, num_instants, n
          print '(f8.1, es16.8, es14.3, 4es14.5)', speed * omega_h % instant(k) * 180.0_dp / pi, e((k - 1) * np + 1), &
               & hub_eps((k - 1) * np + 1) - sin(speed * omega_h % instant(k) / 2.0_dp), tip_h % value(k * np, :)
       end do

       ! Report mass and ks over the beams, then with each_blade the
       ! individual blade_ks values in flexible order. Their expansion
       ! derivatives follow the mass and the two ks derivatives in each
       ! row of the table; designs supplies every coordinate's name.
       mass_h = mass % at(solution)
       ks_h   = ks % at(solution)
       print '(a, es24.16)', '   the blades'' mass over the revolutions, kg s     ', mass_h % value(1, 1)
       print '(a, es24.16)', '   the smooth maximum of their von Mises function   ', ks_h % value(1, 1)
       if (each_blade) then
          allocate(blade_ks_h(size(flexible)))
          do s = 1, size(flexible)
             blade_ks_h(s) = blade_ks(s) % at(solution)
             print '(a, a, a, es24.16)', '   the smooth maximum of the blade at ', degrees(flexible(s) - 1), '  ', &
                  & blade_ks_h(s) % value(1, 1)
          end do
       end if
       if (strips > 0 .and. optimize) then
          if (present(design)) then
             if (design == 'optimized') then
                names = designs
                found = solution % fields(names)
                found_thickness = found % value(1, :)
                print '(a)', '   the thickness found, strip by strip, and the multiplier of the stress'
                print '(8es16.8)', found_thickness
                ! One component for the shared stress bound, or one
                ! per beam in flexible order with each_blade.
                found = solution % fields(['mu_ks'])
                print '(a, *(es24.16))', '   mu_ks ', found % value(1, :)
             end if
          end if
       else if (strips > 0 .and. sensitivities) then
          kappa_h = solution % fields(['kappa'])
          print '(a)', '   along the thickness of each strip: the mass''s, the smooth maximum''s by the expansion and &
               &by the reverse sweep'
          do k = 1, size(designs)
             if (each_blade) then
                print '(a8, *(es24.16))', designs(k), mass_h % value(1, 1 + k), ks_h % value(1, 1 + k), &
                     & -kappa_h % value(1, k), (blade_ks_h(s) % value(1, 1 + k), s = 1, size(flexible))
             else
                print '(i8, 3es24.16)', k, mass_h % value(1, 1 + k), ks_h % value(1, 1 + k), -kappa_h % value(1, k)
             end if
          end do
       end if
    end if

    ! Draw each beam with its own strips at the design found. flexible
    ! locates its azimuth, and o locates its set in found_thickness;
    ! without each_blade all beams use the shared set. The rigid blade
    ! surfaces retain thickness(1).
    if (len_trim(meshes) > 0) then
       drawing = 'rotor'
       if (.not. alone) drawing = drawing // '_' // trim(name)
       if (present(design)) drawing = drawing // '_' // design
       drawn = rotor
       if (allocated(found_thickness)) then
          do k = 1, size(drawn % beam_surfaces)
             o = 0
             if (each_blade) o = (k - 1) * strips
             drawn % beam_surfaces(k) = plate_surface(flexible(k) - 1, found_thickness(o + 1:o + strips))
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

  ! Name strip k by its blade's azimuth: h<deg>_<k>, with b the blade
  ! index at b - 1 quarter turns (for example strip(2, 1) is h90_1).
  function strip(b, k) result(name)
    integer, intent(in) :: b, k
    character(len=8) :: name
    write(name, '(a, a, a, i0)') 'h', trim(degrees(b - 1)), '_', k
  end function strip

end program rotorcraft_dynamics
