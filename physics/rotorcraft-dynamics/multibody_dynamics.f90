!=====================================================================!
! RIGID BODIES AND THEIR JOINTS, stated as fields of time, from which
! a program forms the Lagrangian of a mechanism (dissertation, section
! 4.1.4): the action of the bodies, and the conditions of the joints
! paired with their reactions.
!
! A BODY is its centre of mass r0 in the reference configuration, its
! mass m and its inertia J about the centre of mass in the reference
! axes. Its coordinates, unknowns of time, are the displacement u of
! the centre of mass and the Euler parameters (eta, eps) of its
! rotation; C = I - 2 eta eps^x + 2 eps^x eps^x takes the reference
! axes to the body's, so that the point x of the reference
! configuration is at r0 + u + C^T (x - r0) and a direction e of the
! body is C^T e. The angular velocity in the body's axes is
!
!      omega = 2 (eta eps' - eta' eps - eps x eps'),
!
! and the action is the kinetic energy m u'.u'/2 + omega.J omega/2
! less the potential energy - m g.u in the gravity g, with the
! condition eta^2 + eps.eps - 1 = 0 paired with a multiplier of its
! own, an unknown of time.
!
! A JOINT is a list of conditions on the coordinates of the bodies it
! joins, each paired with its reaction, an unknown of time:
!
!      spherical(a, x, b)          the point x of a and of b coincide;
!                                  of a and the fixed point x when b
!                                  is absent
!      revolute(a, x, e, b)        spherical, and the direction e of a
!                                  normal to the two directions of b
!                                  normal to e: a hinge about e; of a
!                                  to the ground when b is absent
!      line(a, e)                  the centre of a on the fixed line
!                                  through it along e
!      cylindrical(a, x, e)        the point x of a on the fixed line
!                                  through x along e, and the direction
!                                  e of a along that line
!      fixed(a)                    a neither moves nor turns
!      rigid(a, b)                 b moves and turns with a
!      revolute_driver(a, e, w, t) a turns about e at the rate w, its
!                                  centre fixed
!      motion_driver(a, d, w, t)   the centre of a moves by d sin(w t)
!
! The two directions normal to e are those of TACS: with i the axis of
! the least component of e, e2 = e x i and e1 = e2 x e. Two rotations
! agree where the vector part of the rotation between them vanishes,
! eta_a eps_b - eta_b eps_a - eps_a x eps_b for the Euler parameters
! of a and b: three conditions independent at every orientation,
! which equal parameters eps_a = eps_b are not where eta = 0.
!
! THE MECHANISM collects the bodies and the joints: the Lagrangian
! density, the sum of the bodies' actions and of the products of the
! conditions with their reactions; the unknowns in the order
! declared; and the initial data, each body moving at the velocity v
! and turning at the angular velocity w about the origin,
!
!      u = 0,  eta = 1,  eps = 0,  u' = v + w x r0,  eta' = 0,
!      eps' = w/2.
!
! Its laws are the stationarities of the density in each unknown: the
! equation of motion of each coordinate, and the condition paired
! with each reaction. Its estimate is each body continuing its initial
! motion, every reaction zero; its energy the bodies' kinetic and
! potential energies summed, and the blades' as beams their density of
! energy over the span (spread_energy).
!
! A DRIVER fixes an actuator - the turn of a body about an axis,
! rotation(a, e), or the displacement of its centre along an axis,
! displacement(a, e) - at the image of the controls u(t) by its row of
! a map: its condition c(q, u) = 0 is paired with its reaction mu as a
! joint's is. The efforts conjugate to the controls are
! e = -(dc/du)^T mu (efforts), and the power they supply, e^T u', is
! read from the motion as mu (dc/dq) q' (power): by the energy theorem
! the rate of the energy, the bodies' and the blades'.
!
! A BLADE AS A BEAM (beam) is a field of time and of its span x, in
! the axes of the body it is fixed to: on a manifold of time and of
! the span, the bodies' coordinates and the reactions are unknowns of
! time alone. The blade's density is a density over the span, and the
! law of an unknown of time alone it reads gains the integral over the
! span of its stationarity there; the laws of the blade's fields are
! the stationarity of that density along t and x, the outer
! derivatives along x read by the elements as fluxes. Its section
! (section) is its properties per length of the span, fields of the
! span and of the design: given, as TACS's TimoshenkoStiffness takes
! them, or implied by a uniform solid rectangle (solid_rectangle) of a
! chord, a thickness and a material, the thickness a field of the span
! and of the design. Its stress (stress) is read where its section
! states it: the rectangle's, the von Mises function at its corners
! and at the middles of its long sides; a section given by its
! properties states none.
!
! A TABLE OF BODIES is read from a file of one body a line, a line
! that begins with # aside: the name, the name of a surface mesh, the
! mass, the centre of mass, the inertia about it in the order xx xy xz
! yy yz zz, and the rate at which the body turns about z at the first
! instant.
!
! THE DRAWING. A body is drawn with a surface, its points in the
! reference configuration and its cells; on a solution, every surface
! moved with its body to an instant, x = r0 + u + C^T (X - r0), is
! one paraview file of that instant, the cells marked by their body,
! and the files of the instants a time series. A blade as a beam is
! drawn with a surface of the blade, each point moved with its
! section of the beam and then with the body the blade is fixed to,
! its deflection a vector at the point.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module multibody_dynamics

  use util_precision      , only : dp
  use util_string         , only : string
  use relation_binary     , only : ragged
  use operation_manifold  , only : continuous_manifold, discrete_manifold, window, region, parameter, points
  use operation_expansion , only : expansion
  use operation_finite_element, only : finite_element
  use operation_field     , only : continuous_field, discrete_field, integral, mean, maximum, sin, cos
  use operation_residual  , only : continuous_residual
  use operation_expression, only : expression
  use view_paraview_writer, only : paraview_writer, time_series

  implicit none

  private
  public :: rigid_body, beam_blade, section, solid_rectangle, mechanism, body_table, table, surface
  public :: spherical, revolute, cylindrical, fixed, rigid, revolute_driver, motion_driver
  public :: line, rotation, displacement
  public :: assembly, attachment, beam

  type :: rigid_body
     character(len=:), allocatable :: name
     real(dp) :: centre(3)    = 0.0_dp
     real(dp) :: mass         = 0.0_dp
     real(dp) :: inertia(3,3) = 0.0_dp
     type(continuous_field) :: u(3), eta, eps(3), norm
   contains
     procedure :: position
     procedure :: moved
     procedure :: direction
     procedure :: angular_velocity
     procedure :: kinetic_energy
     procedure :: potential_energy
  end type rigid_body

  !===================================================================!
  ! AN ACTUATOR: a coordinate of a body that a driver fixes at the
  ! image of the controls - its turn about an axis, fixed by the
  ! agreement of its orientation with the turn by the target angle, or
  ! the displacement of its centre along an axis.
  !===================================================================!

  type :: actuator
     integer          :: kind = 0
     type(rigid_body) :: body
     real(dp)         :: axis(3) = 0.0_dp
     ! of an assembly, the name of the body, which its mechanism forms
     character(len=32) :: name = ''
  end type actuator

  integer, parameter :: TURN = 1, SHIFT = 2

  !===================================================================!
  ! THE SECTION OF A BLADE AS A BEAM, per length of its span: its mass
  ! m; its moments of inertia about the span, the chord and the
  ! thickness; its stiffnesses EA along the span, EI_f in bending along
  ! the thickness, EI_c in bending along the chord and GJ in twist -
  ! fields of the span and of the design. Given so, as TACS's
  ! TimoshenkoStiffness takes them (mA; m11, m22, m33; EA; EI22, EI33;
  ! GJ), a section states no point where its stress is read.
  !===================================================================!

  type :: section
     type(continuous_field) :: mass, axial, flap, lag, torsion
     ! an allocatable array: a fixed array of fields in a section copied
     ! as a class is released twice by the compiler in use
     type(continuous_field), allocatable :: moments(:)
     ! the strength its stress is measured against, where it states one
     real(dp) :: strength = 0.0_dp
   contains
     procedure :: stress   => section_stress
     procedure :: pullback => section_pullback
  end type section

  interface section
     module procedure given_section
  end interface section

  !===================================================================!
  ! A UNIFORM SOLID RECTANGLE of chord c and thickness h, of the density
  ! rho and the moduli E and nu, a section whose properties are implied:
  !
  !      m = rho c h,  rho (c h^3 + h c^3)/12,  rho c h^3/12,  rho h c^3/12,
  !      EA = E c h,  EI_f = E c h^3/12,  EI_c = E h c^3/12,
  !      GJ = G (1 - 0.63 h/c) c h^3/3,
  !
  ! G = E/(2 (1 + nu)), GJ of the thin rectangle; the chord a number,
  ! the thickness a field of the span and of the design. Its stress is
  ! read at its four corners and at the middles of its long sides.
  !===================================================================!

  type, extends(section) :: solid_rectangle
     real(dp) :: chord = 0.0_dp, modulus = 0.0_dp, shear = 0.0_dp
     type(continuous_field) :: thickness
   contains
     procedure :: stress   => rectangle_stress
     procedure :: pullback => rectangle_pullback
  end type solid_rectangle

  interface solid_rectangle
     module procedure uniform_solid_rectangle
     module procedure rectangle_of_strength
  end interface solid_rectangle

  ! a blade as a beam in the axes of the body it is fixed to: its root
  ! and its axes - span, chord, thickness - in that body's reference
  ! axes, its section, and its fields of time and span
  type :: beam_blade
     character(len=:), allocatable :: name
     real(dp) :: root(3) = 0.0_dp, axes(3,3) = 0.0_dp
     class(section), allocatable :: section
     type(continuous_field) :: x, u, v, w, phi
     ! the section's pitch along the span, nose up positive, when it is
     ! twisted: its axes turned about the span from the chord and the
     ! thickness
     logical :: pitched = .false.
     type(continuous_field) :: pitch
   contains
     procedure :: stress
     procedure :: strain_energy
     procedure :: curvatures
  end type beam_blade

  ! the points of a surface in the reference configuration, and its
  ! cells: their points, padded to four, their counts and their gmsh
  ! numbers
  type :: surface
     real(dp), allocatable :: points(:,:)
     integer , allocatable :: corners(:,:), count(:), kinds(:)
  end type surface

  type :: mechanism
     type(continuous_field) :: t
     real(dp) :: gravity(3) = [0.0_dp, 0.0_dp, -9.81_dp]
     type(continuous_field) :: density, energy
     logical :: stated = .false.
     type(continuous_field), allocatable :: unknowns(:)
     ! the data on t = 0 of each coordinate the bodies and the blades
     ! declare, in their order: its value, its rate, the coordinate, and
     ! the member it belongs to
     type(continuous_field), allocatable :: datum_value(:), datum_rate(:), datum_field(:)
     integer               , allocatable :: datum_member(:)
     integer :: members = 0
     type(continuous_field), allocatable :: guess(:)
     ! the density over the span of the blades as beams, the span
     ! coordinate and the space factor, the unknowns of time and span,
     ! and the conditions at the root
     type(continuous_field) :: distributed, x
     type(continuous_manifold) :: space
     logical :: spread = .false.
     ! the energy density of the blades as beams over the span, their
     ! kinetic, strain and gravitational energy per length
     type(continuous_field) :: spread_energy
     ! THE DRIVERS: of each, its actuator, its row of the map from the
     ! controls to the actuators, its target - the row applied to the
     ! controls - and its reaction
     type(actuator)        , allocatable :: drivers(:)
     real(dp)              , allocatable :: driver_row(:,:)
     type(continuous_field), allocatable :: driver_target(:), driver_reaction(:)
     logical, allocatable :: of_space(:)
     type(continuous_field), allocatable :: root(:)
     ! the bodies drawn, and the surface of each; the blades as beams
     ! drawn, the body each is fixed to and its surface
     type(rigid_body), allocatable :: drawn(:)
     type(surface)   , allocatable :: surfaces(:)
     type(beam_blade), allocatable :: drawn_beams(:)
     type(rigid_body), allocatable :: beam_bodies(:)
     type(surface)   , allocatable :: beam_surfaces(:)
     ! the conditions of its joints and drivers, in their order
     type(continuous_field), allocatable :: conditions(:)
     ! OF AN ASSEMBLY: the manifold it is stated on and its blades as
     ! beams
     type(continuous_manifold) :: omega
     type(beam_blade), allocatable :: blades(:)
   contains
     procedure :: body
     procedure :: beam => beam_member
     procedure :: joint
     procedure :: driver
     procedure :: efforts
     procedure :: power
     procedure :: laws
     procedure :: initial
     generic   :: estimate => estimate_on, estimate_image
     procedure, private :: estimate_on, estimate_image
     procedure :: stress => mechanism_stress
     procedure :: attachments
     procedure :: audit
     procedure :: reactions
     procedure, private :: time_unknown
     procedure :: draw_with
     procedure :: draw_beam_with
     procedure :: draw
  end type mechanism

  interface mechanism
     module procedure create_mechanism
     module procedure mechanism_of_assembly
  end interface mechanism

  type :: body_table
     character(len=32), allocatable :: name(:), mesh(:)
     real(dp), allocatable :: mass(:), centre(:,:), inertia(:,:,:), rate(:)
   contains
     procedure :: row
  end type body_table

  !===================================================================!
  ! A COUPLING of bodies by their names: the kind of its joint, the
  ! bodies - the second absent for the ground - and its point and axis.
  ! A mechanism of an assembly forms its conditions on its bodies.
  !===================================================================!

  type :: coupling
     integer           :: kind = 0
     character(len=32) :: first = '', second = ''
     real(dp)          :: point(3) = 0.0_dp, axis(3) = 0.0_dp
  end type coupling

  integer, parameter :: SPHERICAL_COUPLING = 1, REVOLUTE_COUPLING = 2, CYLINDRICAL_COUPLING = 3, &
       & FIXED_COUPLING = 4, RIGID_COUPLING = 5, LINE_COUPLING = 6

  ! A FRAME of a body by its name: a point of it and three axes, the
  ! first along the span of the member attached to it there
  type :: frame
     character(len=32) :: body = ''
     real(dp)          :: point(3) = 0.0_dp, axes(3,3) = 0.0_dp
  end type frame

  interface attachment
     module procedure frame_of
  end interface attachment

  ! A MEMBER'S MODEL: a beam of a section, its pitch along the span
  ! where it is twisted
  type :: beam_model
     class(section), allocatable :: section
     logical :: pitched = .false.
     type(continuous_field) :: pitch
  end type beam_model

  ! A FLEXIBLE MEMBER of an assembly: its name, the region and the
  ! points of its domain, its attachment and its model
  type :: flexible_member
     character(len=32) :: name = ''
     type(region)      :: domain
     integer           :: points = 0
     type(frame)       :: attachment
     type(beam_model), allocatable :: model
  end type flexible_member

  !===================================================================!
  ! THE ASSEMBLY: the members and their domains, the joints and the
  ! drivers, by name, and the design. Its manifold Sigma is the
  ! members' domain times the design - the flexible members of one
  ! domain, the rigid ones of none - on which the design's fields, the
  ! sections and the pitches, are stated; a mechanism forms them on a
  ! product of time with Sigma by their pullback.
  !===================================================================!

  type :: assembly
     type(body_table)  :: bodies
     real(dp)          :: gravity(3) = [0.0_dp, 0.0_dp, -9.81_dp]
     type(continuous_manifold) :: sigma
     type(flexible_member), allocatable :: members(:)
     type(coupling)   , allocatable :: couplings(:)
     character(len=64), allocatable :: coupling_name(:)
     type(actuator)   , allocatable :: actuators(:)
     character(len=64), allocatable :: driver_name(:)
     real(dp)         , allocatable :: driver_row(:,:)
   contains
     procedure :: manifold   => assembly_manifold
     procedure :: mass       => assembly_mass
     procedure :: discretize => assembly_discretize
     procedure :: member     => assembly_member
     procedure :: model      => assembly_model
     procedure :: joint      => assembly_joint
     procedure :: driver     => assembly_driver
     generic   :: coordinate => design_coordinate, member_coordinate
     procedure, private :: design_coordinate, member_coordinate
  end type assembly

  interface assembly
     module procedure create_assembly
  end interface assembly

  interface table
     module procedure read_table
  end interface table

  ! THE JOINTS AND THE ACTUATORS of bodies, by the bodies themselves -
  ! their conditions, or the actuator - or by their names - a coupling,
  ! or the actuator of the named body, which a mechanism of an assembly
  ! forms on its bodies by the same procedures
  interface spherical
     module procedure spherical_of_bodies, spherical_of_names
  end interface spherical

  interface revolute
     module procedure revolute_of_bodies, revolute_of_names
  end interface revolute

  interface cylindrical
     module procedure cylindrical_of_body, cylindrical_of_name
  end interface cylindrical

  interface fixed
     module procedure fixed_body, fixed_name
  end interface fixed

  interface rigid
     module procedure rigid_of_bodies, rigid_of_names
  end interface rigid

  interface line
     module procedure line_of_body, line_of_name
  end interface line

  interface rotation
     module procedure rotation_of_body, rotation_of_name
  end interface rotation

  interface displacement
     module procedure displacement_of_body, displacement_of_name
  end interface displacement

contains

  function create_mechanism(t, gravity) result(this)

    type(continuous_field), intent(in) :: t
    real(dp)              , intent(in), optional :: gravity(3)
    type(mechanism) :: this

    this % t = t
    if (present(gravity)) this % gravity = gravity
    allocate(this % unknowns(0), this % datum_value(0), this % datum_rate(0), this % datum_field(0), &
         & this % datum_member(0), this % guess(0), this % drawn(0), this % surfaces(0))
    allocate(this % of_space(0), this % root(0), this % drawn_beams(0), this % beam_bodies(0), this % beam_surfaces(0))
    allocate(this % conditions(0), this % blades(0))

  end function create_mechanism

  !===================================================================!
  ! A BODY of the mechanism: its coordinates and the multiplier of
  ! their norm declared on the manifold under the body's name, its
  ! action added to the density, its initial data and its estimate.
  !===================================================================!

  subroutine body(this, omega, b, name, centre, mass, inertia, velocity, angular_velocity)

    class(mechanism)         , intent(inout) :: this
    type(continuous_manifold), intent(inout) :: omega
    type(rigid_body)         , intent(out)   :: b
    character(len=*)         , intent(in)    :: name
    real(dp)                 , intent(in)    :: centre(3), mass, inertia(3,3)
    real(dp)                 , intent(in), optional :: velocity(3), angular_velocity(3)

    type(continuous_field) :: action, zero, unit, theta
    real(dp) :: v(3), w(3), rate, axis(3), start(3), turned(3)
    character(len=1), parameter :: digit(3) = ['1', '2', '3']
    integer :: k

    v = 0.0_dp
    w = 0.0_dp
    if (present(velocity))         v = velocity
    if (present(angular_velocity)) w = angular_velocity
    b % name    = name
    b % centre  = centre
    b % mass    = mass
    b % inertia = inertia
    do k = 1, 3
       b % u(k) = this % time_unknown(omega, name // '_u' // digit(k))
    end do
    b % eta = this % time_unknown(omega, name // '_eta')
    do k = 1, 3
       b % eps(k) = this % time_unknown(omega, name // '_eps' // digit(k))
    end do
    b % norm = this % time_unknown(omega, name // '_norm')
    this % unknowns = [this % unknowns, b % u, b % eta, b % eps, b % norm]
    this % of_space = [this % of_space, (.false., k = 1, 8)]

    action = b % kinetic_energy(this % t) - b % potential_energy(this % gravity) &
         & + b % norm * (b % eta * b % eta + b % eps(1) * b % eps(1) + b % eps(2) * b % eps(2) &
         & + b % eps(3) * b % eps(3) - 1.0_dp)
    call add(this, action)
    if (this % stated .and. allocated(this % energy % component)) then
       this % energy = this % energy + b % kinetic_energy(this % t) + b % potential_energy(this % gravity)
    else
       this % energy = b % kinetic_energy(this % t) + b % potential_energy(this % gravity)
    end if

    ! the initial data: at rest in the reference configuration but for
    ! the velocity of the centre and the rate of the rotation
    start = v + cross_product(w, centre)
    this % members      = this % members + 1
    this % datum_value  = [this % datum_value, (b % u(k), k = 1, 3), b % eta - 1.0_dp, (b % eps(k), k = 1, 3)]
    this % datum_rate   = [this % datum_rate, (b % u(k) % derivative([this % t]) - start(k), k = 1, 3), &
         & b % eta % derivative([this % t]), (b % eps(k) % derivative([this % t]) - w(k) / 2.0_dp, k = 1, 3)]
    this % datum_field  = [this % datum_field, (b % u(k), k = 1, 3), b % eta, (b % eps(k), k = 1, 3)]
    this % datum_member = [this % datum_member, (this % members, k = 1, 7)]

    ! the estimate: the initial motion continued, the rotation by the
    ! angle |w| t about w and the translation v t
    zero = continuous_field(omega, [0.0_dp])
    unit = continuous_field(omega, [1.0_dp])
    rate = norm2(w)
    if (rate > 0.0_dp) then
       axis   = w / rate
       turned = cross_product(axis, centre)
       theta  = rate * this % t
       this % guess = [this % guess, (v(k) * this % t + (cos(theta) - 1.0_dp) * (centre(k) - axis(k) &
            & * dot_product(axis, centre)) + turned(k) * sin(theta), k = 1, 3), &
            & cos(theta / 2.0_dp), (axis(k) * sin(theta / 2.0_dp), k = 1, 3), zero]
    else
       this % guess = [this % guess, (v(k) * this % t, k = 1, 3), unit, zero, zero, zero, zero]
    end if

  end subroutine body

  !===================================================================!
  ! A JOINT of the mechanism: a reaction declared for each condition
  ! under the joint's name, the products of the conditions with their
  ! reactions added to the density.
  !===================================================================!

  subroutine joint(this, omega, name, conditions)

    class(mechanism)         , intent(inout) :: this
    type(continuous_manifold), intent(inout) :: omega
    character(len=*)         , intent(in)    :: name
    type(continuous_field)   , intent(in)    :: conditions(:)

    type(continuous_field) :: reaction
    character(len=8) :: label
    integer :: k

    do k = 1, size(conditions)
       write(label, '(i0)') k
       reaction = this % time_unknown(omega, name // '_' // trim(label))
       this % unknowns = [this % unknowns, reaction]
       this % of_space = [this % of_space, .false.]
       this % guess    = [this % guess, continuous_field(omega, [0.0_dp])]
       call add(this, reaction * conditions(k))
    end do
    this % conditions = [this % conditions, conditions]

  end subroutine joint

  !===================================================================!
  ! A DRIVER fixes an actuator at its target, the row of the map from
  ! the controls applied to them: a turn by the condition
  ! cos(psi/2) (eps.e) - sin(psi/2) eta = 0, the axial component of the
  ! vector part of the rotation between the body and the turn by psi
  ! about e, which agrees at every turn; a displacement by u.e - psi =
  ! 0. Its condition and reaction are a joint's, the reaction named
  ! name_1.
  !===================================================================!

  subroutine driver(this, omega, name, act, row, controls)

    class(mechanism)         , intent(inout) :: this
    type(continuous_manifold), intent(inout) :: omega
    character(len=*)         , intent(in)    :: name
    type(actuator)           , intent(in)    :: act
    real(dp)                 , intent(in)    :: row(:)
    type(continuous_field)   , intent(in)    :: controls(:)

    type(continuous_field) :: target, condition
    logical :: started
    integer :: i

    if (size(row) /= size(controls)) error stop 'multibody_dynamics: a driver''s row has one entry per control'
    started = .false.
    do i = 1, size(row)
       if (row(i) == 0.0_dp) cycle
       if (started) then
          target = target + row(i) * controls(i)
       else
          target  = row(i) * controls(i)
          started = .true.
       end if
    end do
    if (.not. started) target = 0.0_dp * controls(1)
    select case (act % kind)
    case (TURN)
       condition = cos(target / 2.0_dp) * dotted(act % body % eps, act % axis) - sin(target / 2.0_dp) * act % body % eta
    case (SHIFT)
       condition = dotted(act % body % u, act % axis) - target
    case default
       error stop 'multibody_dynamics: a driver fixes a turn or a displacement'
    end select
    call this % joint(omega, name, [condition])
    if (.not. allocated(this % drivers)) then
       allocate(this % drivers(0), this % driver_target(0), this % driver_reaction(0), this % driver_row(0, size(row)))
    end if
    if (size(this % driver_row, 2) /= size(row)) error stop 'multibody_dynamics: every driver reads the same controls'
    this % drivers         = [this % drivers, act]
    this % driver_target   = [this % driver_target, target]
    this % driver_reaction = [this % driver_reaction, this % unknowns(size(this % unknowns))]
    this % driver_row      = reshape([transpose(this % driver_row), row], [size(row), size(this % drivers)])
    this % driver_row      = transpose(this % driver_row)

  end subroutine driver

  !===================================================================!
  ! THE EFFORTS conjugate to the controls, e = -(dc/du)^T mu summed over
  ! the drivers: of a displacement its reaction times its row; of a
  ! turn its reaction times (sin(psi/2) (eps.e) + cos(psi/2) eta) / 2
  ! times its row, the reaction over two at agreement.
  !===================================================================!

  function efforts(this) result(e)

    class(mechanism), intent(in) :: this
    type(continuous_field), allocatable :: e(:)

    type(continuous_field) :: term
    logical, allocatable :: started(:)
    integer :: d, i

    allocate(e(size(this % driver_row, 2)), started(size(this % driver_row, 2)))
    started = .false.
    do d = 1, size(this % drivers)
       associate (act => this % drivers(d), psi => this % driver_target(d))
         if (act % kind == TURN) then
            term = this % driver_reaction(d) * (sin(psi / 2.0_dp) * dotted(act % body % eps, act % axis) &
                 & + cos(psi / 2.0_dp) * act % body % eta) / 2.0_dp
         else
            term = this % driver_reaction(d)
         end if
       end associate
       do i = 1, size(this % driver_row, 2)
          if (this % driver_row(d, i) == 0.0_dp) cycle
          if (started(i)) then
             e(i) = e(i) + this % driver_row(d, i) * term
          else
             e(i) = this % driver_row(d, i) * term
             started(i) = .true.
          end if
       end do
    end do
    do i = 1, size(e)
       if (.not. started(i)) e(i) = 0.0_dp * this % t
    end do

  end function efforts

  !===================================================================!
  ! THE POWER the controls supply, e^T u', read from the motion: each
  ! driver's reaction times the rate of its condition through the
  ! coordinates, mu (dc/dq) q', which equals -mu (dc/du) u' where the
  ! condition is zero along the motion - for a displacement u'.e, for
  ! a turn cos(psi/2) (eps'.e) - sin(psi/2) eta'; no rate of a control
  ! is read.
  !===================================================================!

  function power(this) result(p)

    class(mechanism), intent(in) :: this
    type(continuous_field) :: p

    type(continuous_field) :: rate
    integer :: d, i3

    p = 0.0_dp * this % t
    do d = 1, size(this % drivers)
       associate (act => this % drivers(d), psi => this % driver_target(d))
         if (act % kind == TURN) then
            rate = cos(psi / 2.0_dp) * dotted([(act % body % eps(i3) % derivative([this % t]), i3 = 1, 3)], act % axis) &
                 & - sin(psi / 2.0_dp) * act % body % eta % derivative([this % t])
         else
            rate = dotted([(act % body % u(i3) % derivative([this % t]), i3 = 1, 3)], act % axis)
         end if
       end associate
       p = p + this % driver_reaction(d) * rate
    end do

  end function power

  !===================================================================!
  ! A BLADE AS A BEAM in the axes of the body it is fixed to: the
  ! displacement u along its span, v along its chord and w along its
  ! thickness, and its twist phi about its span, fields of time and of
  ! the span x, declared on the manifold under the blade's name. The
  ! centre of the section at x is root + (x + u) e_s + v e_c + w e_f
  ! in the body's reference axes, which the body moves; its section s
  ! is its properties per length (section). The blade's density over
  ! the span is the kinetic energy of the section - of its centre, and
  ! of its turning at the body's angular velocity and the twist's rate
  ! about e_s with the section's moments of inertia - less its strain
  ! energy and its potential energy:
  !
  !      m |p_t|^2/2 + w.J w/2 - EA eps^2/2 - EI_f w_xx^2/2
  !        - EI_c v_xx^2/2 - GJ phi_x^2/2 + m g.p,
  !
  ! the axial strain to second order in the slopes, eps = u_x + (v_x^2
  ! + w_x^2)/2, which stiffens the turning blade. At the root the blade
  ! is fixed to the body, u, v, w, phi, v_x and w_x zero; it begins at
  ! rest in the body's axes, and its estimate is zero.
  !
  ! A TWISTED BLADE: of the pitch theta(x), nose up positive, its
  ! section's axes are e_c cos theta + e_f sin theta and -e_c sin theta
  ! + e_f cos theta; its moments of inertia turn with them, and its
  ! bending stiffnesses take the curvatures along them,
  ! k_f = w_xx cos theta - v_xx sin theta and k_c = v_xx cos theta +
  ! w_xx sin theta (curvatures).
  !===================================================================!

  subroutine beam_member(this, omega, blade, b, name, x, root, axes, s, pitch)

    class(mechanism)         , intent(inout) :: this
    type(continuous_manifold), intent(inout) :: omega
    type(beam_blade)         , intent(out)   :: blade
    type(rigid_body)         , intent(in)    :: b
    character(len=*)         , intent(in)    :: name
    type(continuous_field)   , intent(in)    :: x
    real(dp)                 , intent(in)    :: root(3), axes(3,3)
    class(section)           , intent(in)    :: s
    type(continuous_field)   , intent(in), optional :: pitch

    type(continuous_field) :: centre(3), p(3), rate, turning(3), body_rate(3), slope_v, slope_w, strain, energy
    type(continuous_field) :: zero, inertia, k_f, k_c, c, sn, turned_moments(3), kinetic
    logical  :: stated
    integer  :: k, i, j

    blade % name = name
    blade % root = root
    blade % axes = axes
    blade % x    = x
    allocate(blade % section, source=s)
    if (present(pitch)) then
       blade % pitched = .true.
       blade % pitch   = pitch
    end if
    blade % u   = omega % unknown(name // '_u')
    blade % v   = omega % unknown(name // '_v')
    blade % w   = omega % unknown(name // '_w')
    blade % phi = omega % unknown(name // '_phi')
    this % unknowns = [this % unknowns, blade % u, blade % v, blade % w, blade % phi]
    this % of_space = [this % of_space, .true., .true., .true., .true.]
    this % x     = x
    this % space = omega % space()

    ! the centre of the section and where the body moves it; its
    ! kinetic energy and that of the section's turning
    do k = 1, 3
       centre(k) = root(k) + axes(k, 1) * (x + blade % u) + axes(k, 2) * blade % v + axes(k, 3) * blade % w
    end do
    p = b % moved(centre)
    rate   = p(1) % derivative([this % t])
    energy = rate * rate
    do k = 2, 3
       rate   = p(k) % derivative([this % t])
       energy = energy + rate * rate
    end do
    energy    = (s % mass / 2.0_dp) * energy
    body_rate = b % angular_velocity(this % t)
    rate      = blade % phi % derivative([this % t])
    do k = 1, 3
       turning(k) = body_rate(k)
       if (axes(k, 1) /= 0.0_dp) turning(k) = turning(k) + axes(k, 1) * rate
    end do
    ! the moments of inertia in the body's reference axes, J_ij the sum
    ! of the moments about the axes k with a_ik a_jk nonzero; of a
    ! twisted blade, about the span and the chord and thickness turned
    ! by the pitch, the product of the two the moments' difference times
    ! cos theta sin theta
    if (blade % pitched) then
       c  = cos(pitch)
       sn = sin(pitch)
       turned_moments(1) = s % moments(1)
       turned_moments(2) = s % moments(2) * c * c + s % moments(3) * sn * sn
       turned_moments(3) = s % moments(2) * sn * sn + s % moments(3) * c * c
       do i = 1, 3
          do j = i, 3
             inertia = turned_moments(1) * (axes(i, 1) * axes(j, 1)) + turned_moments(2) * (axes(i, 2) * axes(j, 2)) &
                  & + turned_moments(3) * (axes(i, 3) * axes(j, 3)) &
                  & + ((s % moments(2) - s % moments(3)) * c * sn) * (axes(i, 2) * axes(j, 3) + axes(i, 3) * axes(j, 2))
             if (i == j) then
                energy = energy + (inertia / 2.0_dp) * turning(i) * turning(i)
             else
                energy = energy + inertia * turning(i) * turning(j)
             end if
          end do
       end do
    end if
    do i = 1, 3
       if (blade % pitched) exit
       do j = i, 3
          stated = .false.
          do k = 1, 3
             if (axes(i, k) * axes(j, k) == 0.0_dp) cycle
             if (stated) then
                inertia = inertia + s % moments(k) * (axes(i, k) * axes(j, k))
             else
                inertia = s % moments(k) * (axes(i, k) * axes(j, k))
                stated  = .true.
             end if
          end do
          if (.not. stated) cycle
          if (i == j) then
             energy = energy + (inertia / 2.0_dp) * turning(i) * turning(i)
          else
             energy = energy + inertia * turning(i) * turning(j)
          end if
       end do
    end do

    ! less the strain energy, and the potential energy in gravity
    kinetic = energy
    slope_v = blade % v % derivative([x])
    slope_w = blade % w % derivative([x])
    strain  = blade % u % derivative([x]) + (slope_v * slope_v + slope_w * slope_w) / 2.0_dp
    if (blade % pitched) then
       call blade % curvatures(k_f, k_c)
       energy  = energy - (s % axial / 2.0_dp) * strain * strain &
            & - (s % flap / 2.0_dp) * k_f**2 - (s % lag / 2.0_dp) * k_c**2 &
            & - (s % torsion / 2.0_dp) * blade % phi % derivative([x])**2 + s % mass * combination(this % gravity, p)
    else
       energy  = energy - (s % axial / 2.0_dp) * strain * strain &
            & - (s % flap / 2.0_dp) * blade % w % derivative([x, x])**2 &
            & - (s % lag / 2.0_dp) * blade % v % derivative([x, x])**2 &
            & - (s % torsion / 2.0_dp) * blade % phi % derivative([x])**2 + s % mass * combination(this % gravity, p)
    end if
    ! the energy density, kinetic and potential, 2 T - (T - V)
    if (this % spread) then
       this % distributed   = this % distributed + energy
       this % spread_energy = this % spread_energy + (2.0_dp * kinetic - energy)
    else
       this % distributed   = energy
       this % spread_energy = 2.0_dp * kinetic - energy
       this % spread        = .true.
    end if

    ! at rest in the body's axes at the first instant; fixed to the
    ! body at the root; the estimate zero
    this % members      = this % members + 1
    this % datum_value  = [this % datum_value, blade % u, blade % v, blade % w, blade % phi]
    this % datum_rate   = [this % datum_rate, blade % u % derivative([this % t]), blade % v % derivative([this % t]), &
         & blade % w % derivative([this % t]), blade % phi % derivative([this % t])]
    this % datum_field  = [this % datum_field, blade % u, blade % v, blade % w, blade % phi]
    this % datum_member = [this % datum_member, (this % members, k = 1, 4)]
    this % root = [this % root, blade % u, blade % v, blade % w, blade % phi, slope_v, slope_w]
    zero = continuous_field(omega, [0.0_dp])
    this % guess = [this % guess, zero, zero, zero, zero]

  end subroutine beam_member

  !===================================================================!
  ! THE STRESS OF A BLADE AS A BEAM, fields of time and of the span:
  ! the von Mises function (sigma^2 + 3 tau^2) / sigma_y^2 of the yield
  ! stress sigma_y at each point of its section where the section
  ! states it, from the axial strain eps, the curvatures v_xx and w_xx
  ! and the twist's rate phi_x; none when the section is given by its
  ! properties alone.
  !===================================================================!

  function stress(blade, yield) result(f)

    class(beam_blade), intent(in) :: blade
    ! the yield stress, or the strength its section states
    real(dp)         , intent(in), optional :: yield
    type(continuous_field), allocatable :: f(:)

    type(continuous_field) :: slope_v, slope_w, strain

    type(continuous_field) :: k_f, k_c
    real(dp) :: strength

    strength = blade % section % strength
    if (present(yield)) strength = yield
    if (.not. strength > 0.0_dp) error stop 'multibody_dynamics: the stress of a beam is measured against a strength'

    slope_v = blade % v % derivative([blade % x])
    slope_w = blade % w % derivative([blade % x])
    strain  = blade % u % derivative([blade % x]) + (slope_v * slope_v + slope_w * slope_w) / 2.0_dp
    if (blade % pitched) then
       call blade % curvatures(k_f, k_c)
       f = blade % section % stress(strain, k_c, k_f, blade % phi % derivative([blade % x]), strength)
       return
    end if
    f = blade % section % stress(strain, blade % v % derivative([blade % x, blade % x]), &
         & blade % w % derivative([blade % x, blade % x]), blade % phi % derivative([blade % x]), strength)

  end function stress

  ! the strain energy of a blade as a beam per length of its span, a
  ! field of time and of the span: (EA eps^2 + EI_f w_xx^2 + EI_c v_xx^2
  ! + GJ phi_x^2)/2, the strain energy the beam's density states
  function strain_energy(blade) result(e)
    class(beam_blade), intent(in) :: blade
    type(continuous_field) :: e
    type(continuous_field) :: slope_v, slope_w, strain, k_f, k_c
    slope_v = blade % v % derivative([blade % x])
    slope_w = blade % w % derivative([blade % x])
    strain  = blade % u % derivative([blade % x]) + (slope_v * slope_v + slope_w * slope_w) / 2.0_dp
    if (blade % pitched) then
       call blade % curvatures(k_f, k_c)
       e = (blade % section % axial / 2.0_dp) * strain * strain + (blade % section % flap / 2.0_dp) * k_f**2 &
            & + (blade % section % lag / 2.0_dp) * k_c**2 &
            & + (blade % section % torsion / 2.0_dp) * blade % phi % derivative([blade % x])**2
       return
    end if
    e = (blade % section % axial / 2.0_dp) * strain * strain &
         & + (blade % section % flap / 2.0_dp) * blade % w % derivative([blade % x, blade % x])**2 &
         & + (blade % section % lag / 2.0_dp) * blade % v % derivative([blade % x, blade % x])**2 &
         & + (blade % section % torsion / 2.0_dp) * blade % phi % derivative([blade % x])**2
  end function strain_energy

  ! the curvatures of a twisted blade along its section's axes: k_f in
  ! bending along the turned thickness, k_c along the turned chord
  subroutine curvatures(blade, k_f, k_c)
    class(beam_blade)     , intent(in)  :: blade
    type(continuous_field), intent(out) :: k_f, k_c
    type(continuous_field) :: v_xx, w_xx, c, sn
    v_xx = blade % v % derivative([blade % x, blade % x])
    w_xx = blade % w % derivative([blade % x, blade % x])
    c    = cos(blade % pitch)
    sn   = sin(blade % pitch)
    k_f  = w_xx * c - v_xx * sn
    k_c  = v_xx * c + w_xx * sn
  end subroutine curvatures

  ! a section given by its properties, as TACS's TimoshenkoStiffness
  ! takes them
  function given_section(mass, moments, axial, flap, lag, torsion) result(this)
    type(continuous_field), intent(in) :: mass, moments(3), axial, flap, lag, torsion
    type(section) :: this
    this % mass    = mass
    this % moments = moments
    this % axial   = axial
    this % flap    = flap
    this % lag     = lag
    this % torsion = torsion
  end function given_section

  ! a section given by its properties states no point where its stress
  ! is read
  function section_stress(this, strain, v_xx, w_xx, phi_x, yield) result(f)
    class(section)        , intent(in) :: this
    type(continuous_field), intent(in) :: strain, v_xx, w_xx, phi_x
    real(dp)              , intent(in) :: yield
    type(continuous_field), allocatable :: f(:)
    allocate(f(0))
  end function section_stress

  ! the uniform solid rectangle's section, its properties implied
  function uniform_solid_rectangle(density, modulus, poisson, chord, thickness) result(this)
    real(dp)              , intent(in) :: density, modulus, poisson, chord
    type(continuous_field), intent(in) :: thickness
    type(solid_rectangle) :: this
    type(continuous_field) :: area
    this % chord     = chord
    this % modulus   = modulus
    this % shear     = modulus / (2.0_dp * (1.0_dp + poisson))
    this % thickness = thickness
    area = chord * thickness
    allocate(this % moments(3))
    this % mass       = density * area
    this % moments(1) = density * ((chord * thickness**3 + thickness * chord**3) / 12.0_dp)
    this % moments(2) = density * (chord * thickness**3 / 12.0_dp)
    this % moments(3) = density * (thickness * chord**3 / 12.0_dp)
    this % axial      = modulus * area
    this % flap       = modulus * chord * thickness**3 / 12.0_dp
    this % lag        = modulus * thickness * chord**3 / 12.0_dp
    this % torsion    = this % shear * (1.0_dp - 0.63_dp * thickness / chord) * chord * thickness**3 / 3.0_dp
  end function uniform_solid_rectangle

  ! the same of a strength, which its stress is measured against
  function rectangle_of_strength(density, modulus, poisson, strength, chord, thickness) result(this)
    real(dp)              , intent(in) :: density, modulus, poisson, strength, chord
    type(continuous_field), intent(in) :: thickness
    type(solid_rectangle) :: this
    this = uniform_solid_rectangle(density, modulus, poisson, chord, thickness)
    this % strength = strength
  end function rectangle_of_strength

  !===================================================================!
  ! THE PULLBACK of a section stated on the base of a product of time
  ! with it - an assembly's - to the product: its fields read there,
  ! its numbers as they are.
  !===================================================================!

  function section_pullback(this, omega) result(pulled)
    class(section)           , intent(in) :: this
    type(continuous_manifold), intent(in) :: omega
    class(section), allocatable :: pulled
    type(section) :: s
    call section_fields(this, omega, s)
    allocate(pulled, source=s)
  end function section_pullback

  function rectangle_pullback(this, omega) result(pulled)
    class(solid_rectangle)   , intent(in) :: this
    type(continuous_manifold), intent(in) :: omega
    class(section), allocatable :: pulled
    type(solid_rectangle) :: r
    call section_fields(this, omega, r)
    r % chord     = this % chord
    r % modulus   = this % modulus
    r % shear     = this % shear
    r % thickness = omega % pullback(this % thickness)
    allocate(pulled, source=r)
  end function rectangle_pullback

  ! the fields every section states, pulled back, and its strength
  subroutine section_fields(this, omega, s)
    class(section)           , intent(in)    :: this
    type(continuous_manifold), intent(in)    :: omega
    class(section)           , intent(inout) :: s
    integer :: k
    s % mass    = omega % pullback(this % mass)
    s % axial   = omega % pullback(this % axial)
    s % flap    = omega % pullback(this % flap)
    s % lag     = omega % pullback(this % lag)
    s % torsion = omega % pullback(this % torsion)
    allocate(s % moments(size(this % moments)))
    do k = 1, size(this % moments)
       s % moments(k) = omega % pullback(this % moments(k))
    end do
    s % strength = this % strength
  end subroutine section_fields

  !===================================================================!
  ! THE STRESS OF A UNIFORM SOLID RECTANGLE at its four corners and at
  ! the middles of its long sides: sigma = E (eps - y v_xx - z w_xx) the
  ! normal stress at the offsets y along the chord and z along the
  ! thickness, tau = G h phi_x the shear of the twist, at the middles of
  ! the long sides, zero at the corners.
  !===================================================================!

  function rectangle_stress(this, strain, v_xx, w_xx, phi_x, yield) result(f)

    class(solid_rectangle), intent(in) :: this
    type(continuous_field), intent(in) :: strain, v_xx, w_xx, phi_x
    real(dp)              , intent(in) :: yield
    type(continuous_field), allocatable :: f(:)

    type(continuous_field) :: lag, flap, shear

    lag   = (this % chord / 2.0_dp) * v_xx
    flap  = (this % thickness / 2.0_dp) * w_xx
    shear = this % shear * this % thickness * phi_x
    allocate(f(6))
    f(1) = (this % modulus * (strain - lag - flap))**2 / yield**2
    f(2) = (this % modulus * (strain + lag - flap))**2 / yield**2
    f(3) = (this % modulus * (strain - lag + flap))**2 / yield**2
    f(4) = (this % modulus * (strain + lag + flap))**2 / yield**2
    f(5) = ((this % modulus * (strain - flap))**2 + 3.0_dp * shear**2) / yield**2
    f(6) = ((this % modulus * (strain + flap))**2 + 3.0_dp * shear**2) / yield**2

  end function rectangle_stress

  subroutine add(this, term)
    type(mechanism)       , intent(inout) :: this
    type(continuous_field), intent(in)    :: term
    if (this % stated) then
       this % density = this % density + term
    else
       this % density = term
       this % stated  = .true.
    end if
  end subroutine add

  ! the laws: the stationarity of the density in each unknown, in the
  ! order declared, the density's derivatives along time read once as
  ! the components of the jet. Of a field of time and span, the
  ! stationarity of the density over the span along t and x; of an
  ! unknown of time alone the density over the span reads, the
  ! integral over the span of its stationarity is added
  function laws(this) result(equations)
    class(mechanism), intent(in) :: this
    type(continuous_field), allocatable :: equations(:)
    type(continuous_field) :: density
    type(expression) :: spread_form
    integer :: k
    density = this % density % jet_form()
    if (this % spread) spread_form = this % distributed % graph(1)
    allocate(equations(size(this % unknowns)))
    do k = 1, size(this % unknowns)
       if (this % of_space(k)) then
          equations(k) = this % distributed % stationarity(this % unknowns(k), [this % t, this % x])
          cycle
       end if
       equations(k) = density % stationarity(this % unknowns(k), this % t)
       if (.not. this % spread) cycle
       if (spread_form % highest_order_read(this % unknowns(k) % index(1)) < 0) cycle
       equations(k) = equations(k) + integral(this % distributed % stationarity(this % unknowns(k), this % t), &
            & over=this % space)
    end do
  end function laws

  ! the data on t = 0 of the orders along time given, member by member
  ! and in each member order by order: of order 0 each coordinate less
  ! its reference value, of order 1 its rate less the rate the member
  ! is given, of order 2 its acceleration; orders 0 and 1 the data of
  ! an initial value problem, 1 and 2 a start at rest, the values then
  ! solved at the first instant
  function initial(this, orders) result(data)
    class(mechanism), intent(in) :: this
    integer         , intent(in) :: orders(:)
    type(continuous_field), allocatable :: data(:)
    integer :: first, last, i, k, n
    n = size(this % datum_member)
    allocate(data(0))
    first = 1
    do while (first <= n)
       last = first
       do while (last < n)
          if (this % datum_member(last + 1) /= this % datum_member(first)) exit
          last = last + 1
       end do
       do i = 1, size(orders)
          select case (orders(i))
          case (0)
             data = [data, this % datum_value(first:last)]
          case (1)
             data = [data, this % datum_rate(first:last)]
          case (2)
             data = [data, (this % datum_field(k) % derivative([this % t, this % t]), k = first, last)]
          case default
             error stop 'multibody_dynamics: the data on t = 0 are of the orders 0, 1 and 2 along time'
          end select
       end do
       first = last + 1
    end do
  end function initial

  ! an unknown of the bodies and the joints: of time alone on a
  ! manifold of time and space, and of the design when the manifold
  ! has one
  function time_unknown(this, omega, name) result(u)
    class(mechanism)         , intent(in)    :: this
    type(continuous_manifold), intent(inout) :: omega
    character(len=*)         , intent(in)    :: name
    type(continuous_field) :: u
    type(continuous_field), allocatable :: arguments(:)
    integer :: k
    if (.not. omega % with_space) then
       u = omega % unknown(name)
       return
    end if
    allocate(arguments(1))
    arguments(1) = this % t
    if (omega % with_design) then
       do k = 1, size(omega % design_parameter)
          arguments = [arguments, omega % coordinate(trim(omega % design_parameter(k) % name))]
       end do
    end if
    u = omega % unknown(name, arguments)
  end function time_unknown

  ! the estimate: the tuple of the unknowns' estimates
  ! the reactions of the joint of a name, the multipliers of its
  ! conditions in their order, unknowns of time alone; a name no joint
  ! has stops the program
  function reactions(this, name) result(r)
    class(mechanism), intent(in) :: this
    character(len=*), intent(in) :: name
    type(continuous_field), allocatable :: r(:)
    character(len=8) :: label
    integer :: k, j
    allocate(r(0))
    do k = 1, size(this % unknowns)
       write(label, '(i0)') size(r) + 1
       do j = 1, size(this % unknowns)
          if (this % unknowns(j) % name == name // '_' // trim(label)) then
             r = [r, this % unknowns(j)]
             exit
          end if
       end do
       if (j > size(this % unknowns)) exit
    end do
    if (size(r) == 0) error stop 'multibody_dynamics: no joint of this name states a reaction: ' // name
  end function reactions

  function estimate_on(this, omega) result(tuple)
    class(mechanism)         , intent(in) :: this
    type(continuous_manifold), intent(in) :: omega
    type(continuous_field) :: tuple
    tuple = continuous_field(omega, this % guess)
  end function estimate_on

  ! of a mechanism of an assembly, the estimate on the image of its
  ! manifold
  function estimate_image(this, omega_h) result(estimates)
    class(mechanism)       , intent(in) :: this
    type(discrete_manifold), intent(in) :: omega_h
    type(discrete_field) :: estimates
    type(continuous_field) :: tuple
    tuple     = this % estimate_on(this % omega)
    estimates = tuple % discretize(omega_h)
  end function estimate_image

  !===================================================================!
  ! THE ASSEMBLY: its bodies, the table's, under gravity, and its
  ! design; its members, joints and drivers are added by name. Sigma is
  ! a manifold of the design until a flexible member gives it the
  ! member's domain; a further member of another domain stops the
  ! program, members of domains of their own not yet stated.
  !===================================================================!

  function create_assembly(bodies, gravity, design) result(this)
    type(body_table), intent(in) :: bodies
    real(dp)        , intent(in) :: gravity(3)
    type(parameter) , intent(in), optional :: design
    type(assembly) :: this
    this % bodies  = bodies
    this % gravity = gravity
    if (present(design)) then
       this % sigma = continuous_manifold(design=design)
    else
       this % sigma = continuous_manifold()
    end if
    allocate(this % members(0), this % couplings(0), this % coupling_name(0), this % actuators(0), this % driver_name(0))
  end function create_assembly

  function assembly_manifold(this) result(sigma)
    class(assembly), intent(in) :: this
    type(continuous_manifold) :: sigma
    sigma = this % sigma
  end function assembly_manifold

  ! a coordinate of the design by its name, one field of its components
  function design_coordinate(this, name) result(c)
    class(assembly) , intent(in) :: this
    character(len=*), intent(in) :: name
    type(continuous_field) :: c
    c = this % sigma % coordinate(name)
  end function design_coordinate

  ! a coordinate of a member's domain by its name
  function member_coordinate(this, member, name) result(c)
    class(assembly) , intent(in) :: this
    character(len=*), intent(in) :: member, name
    type(continuous_field) :: c
    integer :: k
    k = member_index(this, member)
    c = this % sigma % coordinate(name)
  end function member_coordinate

  ! a flexible member on a domain of its own points, attached to a frame
  ! of a body at the first point of its domain
  subroutine assembly_member(this, name, domain, points, attachment)
    class(assembly) , intent(inout) :: this
    character(len=*), intent(in)    :: name
    type(region)    , intent(in)    :: domain
    integer         , intent(in)    :: points
    type(frame)     , intent(in)    :: attachment
    type(flexible_member) :: m
    if (.not. allocated(domain % extent)) then
       error stop 'multibody_dynamics: a flexible member is on an interval of its span'
    end if
    if (size(this % members) == 0) then
       this % sigma % with_space = .true.
       this % sigma % geometry   = domain
    else
       if (domain % extent % a /= this % members(1) % domain % extent % a .or. &
            & domain % extent % b /= this % members(1) % domain % extent % b .or. points /= this % members(1) % points) then
          error stop 'multibody_dynamics: the flexible members of an assembly share one domain and its points; &
               &members of domains of their own are not yet stated'
       end if
    end if
    m % name       = name
    m % domain     = domain
    m % points     = points
    m % attachment = attachment
    this % members = [this % members, m]
  end subroutine assembly_member

  subroutine assembly_model(this, name, model)
    class(assembly) , intent(inout) :: this
    character(len=*), intent(in)    :: name
    type(beam_model), intent(in)    :: model
    integer :: k
    k = member_index(this, name)
    this % members(k) % model = model
  end subroutine assembly_model

  subroutine assembly_joint(this, name, c)
    class(assembly) , intent(inout) :: this
    character(len=*), intent(in)    :: name
    type(coupling)  , intent(in)    :: c
    this % couplings     = [this % couplings, c]
    this % coupling_name = [character(len=64) :: this % coupling_name, name]
  end subroutine assembly_joint

  ! a driver of an actuator by its row of the map from the controls
  subroutine assembly_driver(this, name, act, row)
    class(assembly) , intent(inout) :: this
    character(len=*), intent(in)    :: name
    type(actuator)  , intent(in)    :: act
    real(dp)        , intent(in)    :: row(:)
    real(dp), allocatable :: rows(:,:)
    integer :: n
    n = size(this % actuators)
    if (n == 0) allocate(this % driver_row(0, size(row)))
    if (size(this % driver_row, 2) /= size(row)) error stop 'multibody_dynamics: every driver reads the same controls'
    allocate(rows(n + 1, size(row)))
    rows(1:n, :)  = this % driver_row
    rows(n + 1, :) = row
    call move_alloc(rows, this % driver_row)
    this % actuators   = [this % actuators, act]
    this % driver_name = [character(len=64) :: this % driver_name, name]
  end subroutine assembly_driver

  ! THE MASS of the assembly's members, kg: the rigid members', the
  ! table's, and the integral of each flexible member's mass per length
  ! over its domain, a field of the design
  function assembly_mass(this) result(m)
    class(assembly), intent(in) :: this
    type(continuous_field) :: m
    type(continuous_field) :: density
    integer :: k
    density = continuous_field(this % sigma, [sum(this % bodies % mass) / this % sigma % measure()])
    do k = 1, size(this % members)
       if (.not. allocated(this % members(k) % model)) error stop 'multibody_dynamics: a flexible member is given its model'
       density = density + this % members(k) % model % section % mass
    end do
    m = integral(density, over=this % sigma)
  end function assembly_mass

  ! THE IMAGE OF SIGMA: the members' points, their beams' hermite
  ! cubics, and the design's expansion
  function assembly_discretize(this, design) result(sigma_h)
    class(assembly), intent(in) :: this
    type(expansion), intent(in) :: design
    type(discrete_manifold) :: sigma_h
    if (size(this % members) > 0) then
       sigma_h = this % sigma % discretize(space=points(this % members(1) % points), design=design, &
            & approximation=finite_element(degree=3))
    else
       sigma_h = this % sigma % discretize(design=design)
    end if
  end function assembly_discretize

  integer function member_index(this, name)
    class(assembly) , intent(in) :: this
    character(len=*), intent(in) :: name
    do member_index = 1, size(this % members)
       if (trim(this % members(member_index) % name) == trim(name)) return
    end do
    error stop 'multibody_dynamics: the assembly has no member named ' // trim(name)
  end function member_index

  !===================================================================!
  ! THE MECHANISM OF AN ASSEMBLY on a product of time with its Sigma,
  ! driven by the controls: the bodies at rest in their reference
  ! configuration, the flexible members as beams of their sections
  ! pulled back, the couplings' conditions, and the drivers on the
  ! controls, in that order.
  !===================================================================!

  function mechanism_of_assembly(omega, layout, controls) result(this)

    type(continuous_manifold), intent(inout) :: omega
    type(assembly)           , intent(in)    :: layout
    type(continuous_field)   , intent(in)    :: controls(:)
    type(mechanism) :: this

    type(rigid_body), allocatable :: part(:)
    type(beam_blade) :: blade
    type(continuous_field) :: x
    class(section), allocatable :: s
    type(actuator) :: act
    integer :: k

    this = create_mechanism(omega % coordinate('t'), layout % gravity)
    allocate(part(size(layout % bodies % name)))
    do k = 1, size(part)
       call this % body(omega, part(k), trim(layout % bodies % name(k)), layout % bodies % centre(:, k), &
            & layout % bodies % mass(k), layout % bodies % inertia(:, :, k))
    end do
    if (size(layout % members) > 0) x = omega % coordinate('x')
    do k = 1, size(layout % members)
       associate (m => layout % members(k))
         if (.not. allocated(m % model)) error stop 'multibody_dynamics: a flexible member is given its model'
         s = m % model % section % pullback(omega)
         if (m % model % pitched) then
            call this % beam(omega, blade, part(layout % bodies % row(m % attachment % body)), trim(m % name), x, &
                 & m % attachment % point, m % attachment % axes, s, pitch=omega % pullback(m % model % pitch))
         else
            call this % beam(omega, blade, part(layout % bodies % row(m % attachment % body)), trim(m % name), x, &
                 & m % attachment % point, m % attachment % axes, s)
         end if
         this % blades = [this % blades, blade]
       end associate
    end do
    do k = 1, size(layout % couplings)
       call this % joint(omega, trim(layout % coupling_name(k)), coupling_conditions(layout % couplings(k), part, &
            & layout % bodies))
    end do
    do k = 1, size(layout % actuators)
       act = layout % actuators(k)
       act % body = part(layout % bodies % row(act % name))
       call this % driver(omega, trim(layout % driver_name(k)), act, layout % driver_row(k, :), controls)
    end do
    this % omega = omega

  end function mechanism_of_assembly

  ! THE STRESS of the flexible members: the squared ratio of the stress
  ! to the strength at every point of their sections, one field of as
  ! many components
  function mechanism_stress(this) result(f)
    class(mechanism), intent(in) :: this
    type(continuous_field) :: f
    type(continuous_field), allocatable :: ratios(:)
    integer :: k
    allocate(ratios(0))
    do k = 1, size(this % blades)
       ratios = [ratios, this % blades(k) % stress()]
    end do
    if (size(ratios) == 0) error stop 'multibody_dynamics: a mechanism of no flexible member states no stress'
    f = continuous_field(this % omega, ratios)
  end function mechanism_stress

  ! THE ATTACHMENTS of the flexible members: their conditions at the
  ! first point of their domain, a residual on that face
  function attachments(this) result(r)
    class(mechanism), intent(in) :: this
    type(continuous_residual) :: r
    type(continuous_field) :: x
    x = this % omega % coordinate('x')
    r = continuous_residual(this % omega % boundary(x, 0.0_dp), this % root)
  end function attachments

  !===================================================================!
  ! THE AUDIT of a scenario's solution:
  !   the start - the efforts at the first instant, the weights the
  !     drivers carry at rest, and the energy there;
  !   the conditions - the largest of the joints' and the drivers' over
  !     the points;
  !   the stress - the largest ratio of the stress to the strength on
  !     the window, the root of the largest squared ratio at its nodes,
  !     against the limit;
  !   the energy balance - the work of the efforts over the window, its
  !     duration times the power's mean there, less the energy gained,
  !     the bodies' by their mean over the window's faces and the
  !     blades' by their integral, and its defect.
  !===================================================================!

  subroutine audit(this, solution, stress, observed, limit)

    class(mechanism)      , intent(in) :: this
    type(discrete_field)  , intent(in) :: solution
    type(continuous_field), intent(in) :: stress
    type(window)          , intent(in) :: observed
    real(dp)              , intent(in) :: limit

    type(continuous_field) :: largest, work
    type(continuous_field), allocatable :: e(:)
    type(discrete_field) :: v
    real(dp), allocatable :: values(:), at_start(:)
    real(dp) :: first_energy, last_energy, w, largest_condition
    integer :: k

    e = this % efforts()
    allocate(at_start(size(e)))
    do k = 1, size(e)
       v = e(k) % at(solution)
       call v % values(values)
       at_start(k) = values(1)
    end do
    largest_condition = 0.0_dp
    do k = 1, size(this % conditions)
       v = this % conditions(k) % at(solution)
       call v % values(values)
       largest_condition = max(largest_condition, maxval(abs(values)))
    end do
    print '(a)', '   the audit'
    print '(a, *(es16.7))', '   the efforts at the first instant: ', at_start
    print '(a, es16.7)', '   the energy at the first instant, J: ', &
         & energy_on(this % omega % boundary(this % t, this % omega % span % a))
    print '(a, es10.2)', '   the largest condition of the joints and the drivers: ', largest_condition
    largest = maximum(stress, over=observed)
    v = largest % at(solution)
    print '(a, f10.6, a, f10.6)', '   the largest ratio of the stress to the strength on the window: ', &
         & sqrt(max(v % value(1, 1), 0.0_dp)), ', the limit ', limit
    work = mean(this % power(), over=observed)
    v = work % at(solution)
    w = (observed % window_end - observed % face_value) * v % value(1, 1)
    first_energy = energy_on(observed % boundary(this % t, observed % face_value))
    last_energy  = energy_on(observed % boundary(this % t, observed % window_end))
    print '(a, 3es20.10)', '   the energy balance on the window, J: the work, the energy gained, the defect ', &
         & w, last_energy - first_energy, w - (last_energy - first_energy)

  contains

    ! the energy on a face of an instant: the bodies' by their mean over
    ! the face, the blades' by their integral
    real(dp) function energy_on(face)
      type(continuous_manifold), intent(in) :: face
      type(continuous_field) :: f
      type(discrete_field) :: g
      f = mean(this % energy, over=face)
      g = f % at(solution)
      energy_on = g % value(1, 1)
      if (this % spread) then
         f = integral(this % spread_energy, over=face)
         g = f % at(solution)
         energy_on = energy_on + g % value(1, 1)
      end if
    end function energy_on

  end subroutine audit

  !===================================================================!
  ! THE DRAWING: a body drawn with a surface, and a blade as a beam
  ! with a surface of its blade - a box, or a mesh - in the reference
  ! configuration; and on a solution, at every instant of the stride
  ! given, the surfaces moved, <name>_<instant>.vtu, with <name>.pvd
  ! the time series. A body's surface moves with it; a point X of a
  ! beam's surface, at the span s = (X - root).e_s and the offsets
  ! c and f along the chord and the thickness, moves with its section,
  !
  !      X + u e_s + v e_c + w e_f + phi e_s x (c e_c + f e_f)
  !        - (c v_x + f w_x) e_s,
  !
  ! the fields and their slopes read from the hermite cubics at s, and
  ! then with the body the beam is fixed to. Each cell is marked by its
  ! body, and at each point is the deflection of the beam, in the
  ! fixed axes, zero on a rigid body.
  !===================================================================!

  subroutine draw_with(this, b, s)
    class(mechanism), intent(inout) :: this
    type(rigid_body), intent(in)    :: b
    type(surface)   , intent(in)    :: s
    this % drawn    = [this % drawn, b]
    this % surfaces = [this % surfaces, s]
  end subroutine draw_with

  subroutine draw_beam_with(this, blade, b, s)
    class(mechanism), intent(inout) :: this
    type(beam_blade), intent(in)    :: blade
    type(rigid_body), intent(in)    :: b
    type(surface)   , intent(in)    :: s
    this % drawn_beams   = [this % drawn_beams, blade]
    this % beam_bodies   = [this % beam_bodies, b]
    this % beam_surfaces = [this % beam_surfaces, s]
  end subroutine draw_beam_with

  subroutine draw(this, solution, name, every)

    class(mechanism)    , intent(in) :: this
    type(discrete_field), intent(in) :: solution
    character(len=*)    , intent(in) :: name
    integer             , intent(in) :: every

    type(discrete_field) :: coordinates, bodies, values, slopes
    type(paraview_writer) :: writer
    character(len=256), allocatable :: files(:)
    character(len=32), allocatable :: names(:)
    real(dp), allocatable :: points(:,:), body(:,:), times(:), moved(:,:), deflection(:,:), node(:)
    real(dp), allocatable :: field(:,:,:), slope(:,:,:)
    integer , allocatable :: corners(:,:), count(:), kinds(:)
    real(dp) :: c(3,3), u(3), eta, eps(3), local(4), rate(4), along, chord, across, offset(3), p(3), h, xi
    character(len=4), parameter :: suffix(4) = ['_u  ', '_v  ', '_w  ', '_phi']
    character(len=256) :: file
    integer :: b, k, j, i, np, nc, p0, c0, ns, m, first, nb, nr, f

    if (this_image() /= 1) return
    nr = size(this % drawn)
    nb = size(this % drawn_beams)
    ! the coordinates of the bodies drawn, and of the bodies the beams
    ! are fixed to
    allocate(names(0))
    do b = 1, nr
       names = [character(len=32) :: names, coordinate_names(this % drawn(b) % name)]
    end do
    if (nr > 0) coordinates = solution % fields(names)
    deallocate(names)
    allocate(names(0))
    do b = 1, nb
       names = [character(len=32) :: names, coordinate_names(this % beam_bodies(b) % name)]
    end do
    if (nb > 0) bodies = solution % fields(names)
    np = sum([(size(this % surfaces(b) % points, 2), b = 1, nr)]) &
         & + sum([(size(this % beam_surfaces(b) % points, 2), b = 1, nb)])
    nc = sum([(size(this % surfaces(b) % count), b = 1, nr)]) + sum([(size(this % beam_surfaces(b) % count), b = 1, nb)])
    allocate(points(3, np), deflection(3, np), corners(4, nc), count(nc), kinds(nc), body(nc, 1))
    allocate(files(0), times(0))
    select type (on => solution % on)
    type is (discrete_manifold)
       ! the points of an instant, along the span when the blades are
       ! beams, and the span coordinate of each
       ns = on % num_cells()
       allocate(node(ns), source=0.0_dp)
       if (nb > 0) then
          do j = 1, ns
             associate (x => on % position(on % point_of(1, j)))
               node(j) = x(size(x))
             end associate
          end do
          allocate(field(4, ns * on % num_instants(), nb), slope(4, ns * on % num_instants(), nb))
          do b = 1, nb
             do f = 1, 4
                ! the name in an array of its own first: a constructor
                ! with a type spec passed straight as an argument
                ! corrupts the memory under the compiler in use
                deallocate(names)
                allocate(names(1))
                names(1) = this % drawn_beams(b) % name // trim(suffix(f))
                values = solution % fields(names)
                slopes = values % derivative([this % x])
                field(f, :, b) = values % value(:, 1)
                slope(f, :, b) = slopes % value(:, 1)
             end do
          end do
       end if
       do k = 1, on % num_instants(), every
          first = (k - 1) * ns + 1
          deflection = 0.0_dp
          p0 = 0
          c0 = 0
          do b = 1, nr
             call body_at(coordinates, b, first, u, eta, eps)
             c = orientation(eta, eps)
             associate (s => this % surfaces(b), r0 => this % drawn(b) % centre)
               moved = s % points
               do i = 1, size(moved, 2)
                  moved(:, i) = r0 + u + matmul(transpose(c), s % points(:, i) - r0)
               end do
               call placed(s, moved, real(b, dp))
             end associate
          end do
          do b = 1, nb
             call body_at(bodies, b, first, u, eta, eps)
             c = orientation(eta, eps)
             associate (s => this % beam_surfaces(b), blade => this % drawn_beams(b), r0 => this % beam_bodies(b) % centre)
               moved = s % points
               do i = 1, size(moved, 2)
                  ! the section through the point, and the point's offsets
                  along  = min(max(dot_product(s % points(:, i) - blade % root, blade % axes(:, 1)), node(1)), node(ns))
                  chord  = dot_product(s % points(:, i) - blade % root, blade % axes(:, 2))
                  across = dot_product(s % points(:, i) - blade % root, blade % axes(:, 3))
                  m = 1
                  do j = 1, ns - 1
                     if (along >= node(j)) m = j
                  end do
                  h  = node(m + 1) - node(m)
                  xi = (along - node(m)) / h
                  ! the fields and their slopes by the hermite cubics
                  do f = 1, 4
                     associate (v0 => field(f, first + m - 1, b), v1 => field(f, first + m, b), &
                          &     s0 => slope(f, first + m - 1, b), s1 => slope(f, first + m, b))
                       local(f) = (1.0_dp - 3.0_dp * xi**2 + 2.0_dp * xi**3) * v0 + h * (xi - 2.0_dp * xi**2 + xi**3) * s0 &
                            & + (3.0_dp * xi**2 - 2.0_dp * xi**3) * v1 + h * (xi**3 - xi**2) * s1
                       rate(f)  = (6.0_dp * xi**2 - 6.0_dp * xi) / h * v0 + (1.0_dp - 4.0_dp * xi + 3.0_dp * xi**2) * s0 &
                            & + (6.0_dp * xi - 6.0_dp * xi**2) / h * v1 + (3.0_dp * xi**2 - 2.0_dp * xi) * s1
                     end associate
                  end do
                  offset = chord * blade % axes(:, 2) + across * blade % axes(:, 3)
                  p = local(1) * blade % axes(:, 1) + local(2) * blade % axes(:, 2) + local(3) * blade % axes(:, 3) &
                       & + local(4) * cross_product(blade % axes(:, 1), offset) &
                       & - (chord * rate(2) + across * rate(3)) * blade % axes(:, 1)
                  moved(:, i) = r0 + u + matmul(transpose(c), s % points(:, i) + p - r0)
                  deflection(:, p0 + i) = matmul(transpose(c), p)
               end do
               call placed(s, moved, real(nr + b, dp))
             end associate
          end do
          write(file, '(a, a, i4.4, a)') trim(name), '_', k, '.vtu'
          writer = paraview_writer(points, ragged(corners, count), kinds)
          if (nb > 0) then
             call writer % write_binary(trim(file), body, [string('body')], deflection, 'deflection')
          else
             call writer % write_binary(trim(file), body, [string('body')])
          end if
          files = [character(len=256) :: files, file]
          times = [times, on % instant(k)]
       end do
    class default
       error stop 'multibody_dynamics: a mechanism is drawn on a solution of its manifold'
    end select
    call time_series(trim(name) // '.pvd', files, times)

  contains

    ! a surface's moved points and its cells, marked by the body
    subroutine placed(s, moved, mark)
      type(surface), intent(in) :: s
      real(dp)     , intent(in) :: moved(:,:), mark
      points(:, p0 + 1:p0 + size(moved, 2)) = moved
      corners(:, c0 + 1:c0 + size(s % count)) = s % corners + p0
      where (corners(:, c0 + 1:c0 + size(s % count)) == p0) corners(:, c0 + 1:c0 + size(s % count)) = 0
      count(c0 + 1:c0 + size(s % count)) = s % count
      kinds(c0 + 1:c0 + size(s % count)) = s % kinds
      body(c0 + 1:c0 + size(s % count), 1) = mark
      p0 = p0 + size(moved, 2)
      c0 = c0 + size(s % count)
    end subroutine placed

  end subroutine draw

  ! the names of a body's coordinates, as the body declares them
  function coordinate_names(n) result(names)
    character(len=*), intent(in) :: n
    character(len=32) :: names(7)
    names = [character(len=32) :: n // '_u1', n // '_u2', n // '_u3', n // '_eta', n // '_eps1', n // '_eps2', &
         & n // '_eps3']
  end function coordinate_names

  ! the coordinates of the b-th body of a field of bodies' coordinates
  ! at a point
  subroutine body_at(coordinates, b, point, u, eta, eps)
    type(discrete_field), intent(in)  :: coordinates
    integer             , intent(in)  :: b, point
    real(dp)            , intent(out) :: u(3), eta, eps(3)
    integer :: j, i
    j   = 7 * (b - 1)
    u   = [(coordinates % value(point, j + i), i = 1, 3)]
    eta = coordinates % value(point, j + 4)
    eps = [(coordinates % value(point, j + 4 + i), i = 1, 3)]
  end subroutine body_at

  ! the rotation matrix of the Euler parameters, C = I - 2 eta eps^x +
  ! 2 eps^x eps^x
  pure function orientation(eta, eps) result(c)
    real(dp), intent(in) :: eta, eps(3)
    real(dp) :: c(3,3)
    c(1, :) = [1.0_dp - 2.0_dp * (eps(2)**2 + eps(3)**2), 2.0_dp * (eps(1) * eps(2) + eta * eps(3)), &
         & 2.0_dp * (eps(1) * eps(3) - eta * eps(2))]
    c(2, :) = [2.0_dp * (eps(2) * eps(1) - eta * eps(3)), 1.0_dp - 2.0_dp * (eps(1)**2 + eps(3)**2), &
         & 2.0_dp * (eps(2) * eps(3) + eta * eps(1))]
    c(3, :) = [2.0_dp * (eps(3) * eps(1) + eta * eps(2)), 2.0_dp * (eps(3) * eps(2) - eta * eps(1)), &
         & 1.0_dp - 2.0_dp * (eps(1)**2 + eps(2)**2)]
  end function orientation

  !===================================================================!
  ! THE TABLE OF BODIES read from a file, and the row of a name.
  !===================================================================!

  function read_table(file) result(table)

    character(len=*), intent(in) :: file
    type(body_table) :: table

    character(len=512) :: line
    character(len=32)  :: name, mesh
    real(dp) :: mass, centre(3), j(6), rate
    integer :: unit, status, n

    open(newunit=unit, file=file, status='old', action='read')
    allocate(table % name(0), table % mesh(0), table % mass(0), table % centre(3, 0), &
         & table % inertia(3, 3, 0), table % rate(0))
    n = 0
    do
       read(unit, '(a)', iostat=status) line
       if (status /= 0) exit
       if (len_trim(line) == 0 .or. adjustl(line(1:1)) == '#') cycle
       read(line, *) name, mesh, mass, centre, j, rate
       n = n + 1
       table % name   = [table % name, name]
       table % mesh   = [table % mesh, mesh]
       table % mass   = [table % mass, mass]
       table % rate   = [table % rate, rate]
       table % centre = reshape([table % centre, centre], [3, n])
       table % inertia = reshape([table % inertia, j(1), j(2), j(3), j(2), j(4), j(5), j(3), j(5), j(6)], [3, 3, n])
    end do
    close(unit)

  end function read_table

  integer function row(this, name)
    class(body_table), intent(in) :: this
    character(len=*) , intent(in) :: name
    do row = 1, size(this % name)
       if (trim(this % name(row)) == trim(name)) return
    end do
    error stop 'multibody_dynamics: the table names no such body'
  end function row

  !===================================================================!
  ! THE KINEMATICS OF A BODY: the point x of the reference
  ! configuration, r0 + u + C^T (x - r0); the direction e of the body,
  ! C^T e = e + 2 eta eps x e + 2 eps x (eps x e); the angular velocity
  ! in the body's axes.
  !===================================================================!

  function position(b, x) result(p)
    class(rigid_body), intent(in) :: b
    real(dp)         , intent(in) :: x(3)
    type(continuous_field) :: p(3)
    integer :: k
    if (all(x == b % centre)) then
       do k = 1, 3
          p(k) = b % u(k) + b % centre(k)
       end do
       return
    end if
    p = b % direction(x - b % centre)
    do k = 1, 3
       p(k) = p(k) + b % u(k) + b % centre(k)
    end do
  end function position

  ! the point given as fields of the reference configuration,
  ! r0 + u + C^T (x - r0)
  function moved(b, x) result(p)
    class(rigid_body)     , intent(in) :: b
    type(continuous_field), intent(in) :: x(3)
    type(continuous_field) :: p(3)
    type(continuous_field) :: d(3), s(3), r(3)
    integer :: k
    do k = 1, 3
       d(k) = x(k) - b % centre(k)
    end do
    s = cross(b % eps, d)
    r = cross(b % eps, s)
    do k = 1, 3
       p(k) = d(k) + 2.0_dp * (b % eta * s(k) + r(k)) + b % u(k) + b % centre(k)
    end do
  end function moved

  function direction(b, e) result(d)
    class(rigid_body), intent(in) :: b
    real(dp)         , intent(in) :: e(3)
    type(continuous_field) :: d(3)
    type(continuous_field) :: s(3), r(3)
    integer :: k
    s = crossed(b % eps, e)
    r = cross(b % eps, s)
    do k = 1, 3
       d(k) = 2.0_dp * (b % eta * s(k) + r(k))
       if (e(k) /= 0.0_dp) d(k) = d(k) + e(k)
    end do
  end function direction

  function angular_velocity(b, t) result(w)
    class(rigid_body)     , intent(in) :: b
    type(continuous_field), intent(in) :: t
    type(continuous_field) :: w(3)
    type(continuous_field) :: rate(3), c(3)
    integer :: k
    do k = 1, 3
       rate(k) = b % eps(k) % derivative([t])
    end do
    c = cross(b % eps, rate)
    do k = 1, 3
       w(k) = 2.0_dp * (b % eta * rate(k) - b % eta % derivative([t]) * b % eps(k) - c(k))
    end do
  end function angular_velocity

  ! m u'.u'/2 + omega.J omega/2
  function kinetic_energy(b, t) result(energy)
    class(rigid_body)     , intent(in) :: b
    type(continuous_field), intent(in) :: t
    type(continuous_field) :: energy
    type(continuous_field) :: w(3), v
    integer :: i, j
    v = b % u(1) % derivative([t])
    energy = v * v
    do i = 2, 3
       v = b % u(i) % derivative([t])
       energy = energy + v * v
    end do
    energy = (b % mass / 2.0_dp) * energy
    w = b % angular_velocity(t)
    do i = 1, 3
       if (b % inertia(i, i) /= 0.0_dp) energy = energy + (b % inertia(i, i) / 2.0_dp) * w(i) * w(i)
       do j = i + 1, 3
          if (b % inertia(i, j) /= 0.0_dp) energy = energy + b % inertia(i, j) * w(i) * w(j)
       end do
    end do
  end function kinetic_energy

  ! - m g.u
  function potential_energy(b, gravity) result(energy)
    class(rigid_body), intent(in) :: b
    real(dp)         , intent(in) :: gravity(3)
    type(continuous_field) :: energy
    energy = combination(-b % mass * gravity, b % u)
  end function potential_energy

  !===================================================================!
  ! THE JOINTS, each the list of its conditions.
  !===================================================================!

  function spherical_of_bodies(a, x, b) result(g)
    type(rigid_body), intent(in) :: a
    real(dp)        , intent(in) :: x(3)
    type(rigid_body), intent(in), optional :: b
    type(continuous_field) :: g(3)
    type(continuous_field) :: pb(3)
    integer :: k
    g = a % position(x)
    if (present(b)) then
       pb = b % position(x)
       do k = 1, 3
          g(k) = g(k) - pb(k)
       end do
    else
       do k = 1, 3
          g(k) = g(k) - x(k)
       end do
    end if
  end function spherical_of_bodies

  ! a hinge about e at the point x, of a to b or, without b, to the
  ! ground
  function revolute_of_bodies(a, x, e, b) result(g)
    type(rigid_body), intent(in) :: a
    real(dp)        , intent(in) :: x(3), e(3)
    type(rigid_body), intent(in), optional :: b
    type(continuous_field) :: g(5)
    type(continuous_field) :: ea(3)
    real(dp) :: e1(3), e2(3)
    call normals(e, e1, e2)
    ea = a % direction(e)
    if (present(b)) then
       g(1:3) = spherical(a, x, b)
       g(4)   = dot(ea, b % direction(e1))
       g(5)   = dot(ea, b % direction(e2))
    else
       g(1:3) = spherical(a, x)
       g(4)   = dotted(ea, e1)
       g(5)   = dotted(ea, e2)
    end if
  end function revolute_of_bodies

  ! the centre of a on the line through it along e
  function line_of_body(a, e) result(g)
    type(rigid_body), intent(in) :: a
    real(dp)        , intent(in) :: e(3)
    type(continuous_field) :: g(2)
    real(dp) :: e1(3), e2(3)
    call normals(e, e1, e2)
    g(1) = dotted(a % u, e1)
    g(2) = dotted(a % u, e2)
  end function line_of_body

  ! the actuators: the turn of a about e, and the displacement of its
  ! centre along e
  function rotation_of_body(a, e) result(this)
    type(rigid_body), intent(in) :: a
    real(dp)        , intent(in) :: e(3)
    type(actuator) :: this
    this % kind = TURN
    this % body = a
    this % axis = e / norm2(e)
  end function rotation_of_body

  function displacement_of_body(a, e) result(this)
    type(rigid_body), intent(in) :: a
    real(dp)        , intent(in) :: e(3)
    type(actuator) :: this
    this % kind = SHIFT
    this % body = a
    this % axis = e / norm2(e)
  end function displacement_of_body

  function cylindrical_of_body(a, x, e) result(g)
    type(rigid_body), intent(in) :: a
    real(dp)        , intent(in) :: x(3), e(3)
    type(continuous_field) :: g(4)
    type(continuous_field) :: d(3), ea(3)
    real(dp) :: e1(3), e2(3)
    call normals(e, e1, e2)
    d  = spherical(a, x)
    ea = a % direction(e)
    g(1) = dotted(d, e1)
    g(2) = dotted(d, e2)
    g(3) = dotted(ea, e1)
    g(4) = dotted(ea, e2)
  end function cylindrical_of_body

  function fixed_body(a) result(g)
    type(rigid_body), intent(in) :: a
    type(continuous_field) :: g(6)
    g = [a % u, a % eps]
  end function fixed_body

  function rigid_of_bodies(a, b) result(g)
    type(rigid_body), intent(in) :: a, b
    type(continuous_field) :: g(6)
    type(continuous_field) :: c(3)
    integer :: k
    g(1:3) = a % position(b % centre)
    c = cross(a % eps, b % eps)
    do k = 1, 3
       g(k)     = g(k) - b % u(k) - b % centre(k)
       g(3 + k) = a % eta * b % eps(k) - b % eta * a % eps(k) - c(k)
    end do
  end function rigid_of_bodies

  !===================================================================!
  ! THE JOINTS AND THE ACTUATORS BY THE NAMES of their bodies, of an
  ! assembly, whose mechanism forms them by the procedures above.
  !===================================================================!

  function spherical_of_names(a, x, b) result(c)
    character(len=*), intent(in) :: a
    real(dp)        , intent(in) :: x(3)
    character(len=*), intent(in), optional :: b
    type(coupling) :: c
    c % kind  = SPHERICAL_COUPLING
    c % first = a
    c % point = x
    if (present(b)) c % second = b
  end function spherical_of_names

  function revolute_of_names(a, x, e, b) result(c)
    character(len=*), intent(in) :: a
    real(dp)        , intent(in) :: x(3), e(3)
    character(len=*), intent(in), optional :: b
    type(coupling) :: c
    c % kind  = REVOLUTE_COUPLING
    c % first = a
    c % point = x
    c % axis  = e
    if (present(b)) c % second = b
  end function revolute_of_names

  function cylindrical_of_name(a, x, e) result(c)
    character(len=*), intent(in) :: a
    real(dp)        , intent(in) :: x(3), e(3)
    type(coupling) :: c
    c % kind  = CYLINDRICAL_COUPLING
    c % first = a
    c % point = x
    c % axis  = e
  end function cylindrical_of_name

  function fixed_name(a) result(c)
    character(len=*), intent(in) :: a
    type(coupling) :: c
    c % kind  = FIXED_COUPLING
    c % first = a
  end function fixed_name

  function rigid_of_names(a, b) result(c)
    character(len=*), intent(in) :: a, b
    type(coupling) :: c
    c % kind   = RIGID_COUPLING
    c % first  = a
    c % second = b
  end function rigid_of_names

  function line_of_name(a, e) result(c)
    character(len=*), intent(in) :: a
    real(dp)        , intent(in) :: e(3)
    type(coupling) :: c
    c % kind  = LINE_COUPLING
    c % first = a
    c % axis  = e
  end function line_of_name

  function rotation_of_name(a, e) result(this)
    character(len=*), intent(in) :: a
    real(dp)        , intent(in) :: e(3)
    type(actuator) :: this
    this % kind = TURN
    this % name = a
    this % axis = e / norm2(e)
  end function rotation_of_name

  function displacement_of_name(a, e) result(this)
    character(len=*), intent(in) :: a
    real(dp)        , intent(in) :: e(3)
    type(actuator) :: this
    this % kind = SHIFT
    this % name = a
    this % axis = e / norm2(e)
  end function displacement_of_name

  ! the conditions of a coupling on the bodies of its names
  function coupling_conditions(c, part, bodies) result(g)
    type(coupling)  , intent(in) :: c
    type(rigid_body), intent(in) :: part(:)
    type(body_table), intent(in) :: bodies
    type(continuous_field), allocatable :: g(:)
    integer :: a, b
    a = bodies % row(c % first)
    b = 0
    if (len_trim(c % second) > 0) b = bodies % row(c % second)
    select case (c % kind)
    case (SPHERICAL_COUPLING)
       if (b > 0) then
          g = spherical(part(a), c % point, part(b))
       else
          g = spherical(part(a), c % point)
       end if
    case (REVOLUTE_COUPLING)
       if (b > 0) then
          g = revolute(part(a), c % point, c % axis, part(b))
       else
          g = revolute(part(a), c % point, c % axis)
       end if
    case (CYLINDRICAL_COUPLING)
       g = cylindrical(part(a), c % point, c % axis)
    case (FIXED_COUPLING)
       g = fixed(part(a))
    case (RIGID_COUPLING)
       g = rigid(part(a), part(b))
    case (LINE_COUPLING)
       g = line(part(a), c % axis)
    case default
       error stop 'multibody_dynamics: a coupling of no kind'
    end select
  end function coupling_conditions

  ! a frame of a body by its name: a point of it and three axes
  function frame_of(body, point, axes) result(this)
    character(len=*), intent(in) :: body
    real(dp)        , intent(in) :: point(3), axes(3,3)
    type(frame) :: this
    this % body  = body
    this % point = point
    this % axes  = axes
  end function frame_of

  ! a member's model: a beam of a section, twisted by its pitch along
  ! the span
  function beam(s, pitch) result(model)
    class(section)        , intent(in) :: s
    type(continuous_field), intent(in), optional :: pitch
    type(beam_model) :: model
    allocate(model % section, source=s)
    if (present(pitch)) then
       model % pitched = .true.
       model % pitch   = pitch
    end if
  end function beam

  ! the rotation by the angle w t about e: eta = cos(w t/2) and eps =
  ! e sin(w t/2)
  function revolute_driver(a, e, rate, t) result(g)
    type(rigid_body)      , intent(in) :: a
    real(dp)              , intent(in) :: e(3), rate
    type(continuous_field), intent(in) :: t
    type(continuous_field) :: g(6)
    type(continuous_field) :: s, c, turned(3)
    real(dp) :: axis(3)
    integer :: k
    axis   = e / norm2(e)
    s      = sin((rate / 2.0_dp) * t)
    c      = cos((rate / 2.0_dp) * t)
    turned = crossed(a % eps, axis)
    do k = 1, 3
       g(k)     = a % u(k)
       g(3 + k) = c * a % eps(k) + s * turned(k)
       if (axis(k) /= 0.0_dp) g(3 + k) = g(3 + k) - axis(k) * s * a % eta
    end do
  end function revolute_driver

  function motion_driver(a, d, rate, t) result(g)
    type(rigid_body)      , intent(in) :: a
    real(dp)              , intent(in) :: d(3), rate
    type(continuous_field), intent(in) :: t
    type(continuous_field) :: g(3)
    type(continuous_field) :: s
    integer :: k
    s = sin(rate * t)
    do k = 1, 3
       g(k) = a % u(k) - d(k) * s
    end do
  end function motion_driver

  !===================================================================!
  ! THE ALGEBRA OF DIRECTIONS: fields and real vectors, a zero
  ! coefficient omitted from each sum.
  !===================================================================!

  ! the two directions normal to e, as TACS takes them
  subroutine normals(e, e1, e2)
    real(dp), intent(in)  :: e(3)
    real(dp), intent(out) :: e1(3), e2(3)
    real(dp) :: i(3)
    i = 0.0_dp
    i(minloc(abs(e), 1)) = 1.0_dp
    e2 = cross_product(e, i)
    e1 = cross_product(e2, e)
  end subroutine normals

  pure function cross_product(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
    c = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)]
  end function cross_product

  ! the field e x x for the fields e and the real vector x
  function crossed(e, x) result(c)
    type(continuous_field), intent(in) :: e(3)
    real(dp)              , intent(in) :: x(3)
    type(continuous_field) :: c(3)
    c(1) = combination([x(3), -x(2)], [e(2), e(3)])
    c(2) = combination([x(1), -x(3)], [e(3), e(1)])
    c(3) = combination([x(2), -x(1)], [e(1), e(2)])
  end function crossed

  ! the field a x b for the fields a and b
  function cross(a, b) result(c)
    type(continuous_field), intent(in) :: a(3), b(3)
    type(continuous_field) :: c(3)
    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)
  end function cross

  function dot(a, b) result(d)
    type(continuous_field), intent(in) :: a(3), b(3)
    type(continuous_field) :: d
    d = a(1) * b(1) + a(2) * b(2) + a(3) * b(3)
  end function dot

  function dotted(a, x) result(d)
    type(continuous_field), intent(in) :: a(3)
    real(dp)              , intent(in) :: x(3)
    type(continuous_field) :: d
    d = combination(x, a)
  end function dotted

  ! the sum of c(k) f(k) over the nonzero coefficients; a field of
  ! value zero when every coefficient is zero
  function combination(c, f) result(s)
    real(dp)              , intent(in) :: c(:)
    type(continuous_field), intent(in) :: f(:)
    type(continuous_field) :: s
    logical :: begun
    integer :: k
    begun = .false.
    do k = 1, size(c)
       if (c(k) == 0.0_dp) cycle
       if (begun) then
          s = s + c(k) * f(k)
       else
          s = c(k) * f(k)
          begun = .true.
       end if
    end do
    if (.not. begun) s = 0.0_dp * f(1)
  end function combination

end module multibody_dynamics
