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
!                                  normal to e: a hinge about e
!      cylindrical(a, x, e)        the point x of a on the fixed line
!                                  through x along e, and the direction
!                                  e of a along that line
!      fixed(a)                    a neither moves nor turns
!      rigid_link(a, b)            b moves and turns with a
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
! potential energies summed.
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
! and the files of the instants a time series.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module multibody_dynamics

  use util_precision      , only : dp
  use util_string         , only : string
  use relation_binary     , only : ragged
  use operation_manifold  , only : continuous_manifold, discrete_manifold
  use operation_field     , only : continuous_field, discrete_field, sin, cos, &
       &                           operator(+), operator(-), operator(*), operator(/)
  use view_paraview_writer, only : paraview_writer, time_series

  implicit none

  private
  public :: rigid_body, mechanism, body_table, surface
  public :: spherical, revolute, cylindrical, fixed, rigid_link, revolute_driver, motion_driver

  type :: rigid_body
     character(len=:), allocatable :: name
     real(dp) :: centre(3)    = 0.0_dp
     real(dp) :: mass         = 0.0_dp
     real(dp) :: inertia(3,3) = 0.0_dp
     type(continuous_field) :: u(3), eta, eps(3), norm
   contains
     procedure :: position
     procedure :: direction
     procedure :: angular_velocity
     procedure :: kinetic_energy
     procedure :: potential_energy
  end type rigid_body

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
     type(continuous_field), allocatable :: initial(:)
     type(continuous_field), allocatable :: guess(:)
     ! the bodies drawn, and the surface of each
     type(rigid_body), allocatable :: drawn(:)
     type(surface)   , allocatable :: surfaces(:)
   contains
     procedure :: body
     procedure :: joint
     procedure :: laws
     procedure :: estimate
     procedure :: draw_with
     procedure :: draw
  end type mechanism

  interface mechanism
     module procedure create_mechanism
  end interface mechanism

  type :: body_table
     character(len=32), allocatable :: name(:), mesh(:)
     real(dp), allocatable :: mass(:), centre(:,:), inertia(:,:,:), rate(:)
   contains
     procedure :: row
  end type body_table

  interface body_table
     module procedure read_table
  end interface body_table

contains

  function create_mechanism(t, gravity) result(this)

    type(continuous_field), intent(in) :: t
    real(dp)              , intent(in), optional :: gravity(3)
    type(mechanism) :: this

    this % t = t
    if (present(gravity)) this % gravity = gravity
    allocate(this % unknowns(0), this % initial(0), this % guess(0), this % drawn(0), this % surfaces(0))

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
       b % u(k) = omega % unknown(name // '_u' // digit(k), [this % t])
    end do
    b % eta = omega % unknown(name // '_eta', [this % t])
    do k = 1, 3
       b % eps(k) = omega % unknown(name // '_eps' // digit(k), [this % t])
    end do
    b % norm = omega % unknown(name // '_norm', [this % t])
    this % unknowns = [this % unknowns, b % u, b % eta, b % eps, b % norm]

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
    this % initial = [this % initial, (b % u(k), k = 1, 3), b % eta - 1.0_dp, (b % eps(k), k = 1, 3), &
         & (b % u(k) % derivative([this % t]) - start(k), k = 1, 3), b % eta % derivative([this % t]), &
         & (b % eps(k) % derivative([this % t]) - w(k) / 2.0_dp, k = 1, 3)]

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
       reaction = omega % unknown(name // '_' // trim(label), [this % t])
       this % unknowns = [this % unknowns, reaction]
       this % guess    = [this % guess, continuous_field(omega, [0.0_dp])]
       call add(this, reaction * conditions(k))
    end do

  end subroutine joint

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
  ! the components of the jet
  function laws(this) result(equations)
    class(mechanism), intent(in) :: this
    type(continuous_field), allocatable :: equations(:)
    type(continuous_field) :: density
    integer :: k
    density = this % density % jet_form()
    allocate(equations(size(this % unknowns)))
    do k = 1, size(this % unknowns)
       equations(k) = density % stationarity(this % unknowns(k), this % t)
    end do
  end function laws

  ! the estimate: the tuple of the unknowns' estimates
  function estimate(this, omega) result(tuple)
    class(mechanism)         , intent(in) :: this
    type(continuous_manifold), intent(in) :: omega
    type(continuous_field) :: tuple
    tuple = continuous_field(omega, this % guess)
  end function estimate

  !===================================================================!
  ! THE DRAWING: a body drawn with a surface; and on a solution, at
  ! every instant of the stride given, the surfaces moved with their
  ! bodies, <name>_<instant>.vtu, with <name>.pvd the time series.
  !===================================================================!

  subroutine draw_with(this, b, s)
    class(mechanism), intent(inout) :: this
    type(rigid_body), intent(in)    :: b
    type(surface)   , intent(in)    :: s
    this % drawn    = [this % drawn, b]
    this % surfaces = [this % surfaces, s]
  end subroutine draw_with

  subroutine draw(this, solution, name, every)

    class(mechanism)    , intent(in) :: this
    type(discrete_field), intent(in) :: solution
    character(len=*)    , intent(in) :: name
    integer             , intent(in) :: every

    type(discrete_field) :: coordinates
    type(paraview_writer) :: writer
    character(len=256), allocatable :: files(:)
    character(len=32), allocatable :: names(:)
    real(dp), allocatable :: points(:,:), body(:,:), times(:), moved(:,:)
    integer , allocatable :: corners(:,:), count(:), kinds(:)
    real(dp) :: c(3,3), u(3), eta, eps(3)
    character(len=256) :: file
    integer :: b, k, j, i, np, nc, p0, c0

    if (this_image() /= 1) return
    allocate(names(0))
    do b = 1, size(this % drawn)
       associate (n => this % drawn(b) % name)
         names = [character(len=32) :: names, n // '_u1', n // '_u2', n // '_u3', n // '_eta', n // '_eps1', &
              & n // '_eps2', n // '_eps3']
       end associate
    end do
    coordinates = solution % fields(names)
    np = sum([(size(this % surfaces(b) % points, 2), b = 1, size(this % surfaces))])
    nc = sum([(size(this % surfaces(b) % count), b = 1, size(this % surfaces))])
    allocate(points(3, np), corners(4, nc), count(nc), kinds(nc), body(nc, 1))
    allocate(files(0), times(0))
    select type (on => solution % on)
    type is (discrete_manifold)
       do k = 1, on % num_instants(), every
          p0 = 0
          c0 = 0
          do b = 1, size(this % drawn)
             j   = 7 * (b - 1)
             u   = [(coordinates % value(k, j + i), i = 1, 3)]
             eta = coordinates % value(k, j + 4)
             eps = [(coordinates % value(k, j + 4 + i), i = 1, 3)]
             c   = rotation(eta, eps)
             associate (s => this % surfaces(b), r0 => this % drawn(b) % centre)
               moved = s % points
               do i = 1, size(moved, 2)
                  moved(:, i) = r0 + u + matmul(transpose(c), s % points(:, i) - r0)
               end do
               points(:, p0 + 1:p0 + size(moved, 2)) = moved
               corners(:, c0 + 1:c0 + size(s % count)) = s % corners + p0
               where (corners(:, c0 + 1:c0 + size(s % count)) == p0) corners(:, c0 + 1:c0 + size(s % count)) = 0
               count(c0 + 1:c0 + size(s % count)) = s % count
               kinds(c0 + 1:c0 + size(s % count)) = s % kinds
               body(c0 + 1:c0 + size(s % count), 1) = real(b, dp)
               p0 = p0 + size(moved, 2)
               c0 = c0 + size(s % count)
             end associate
          end do
          write(file, '(a, a, i4.4, a)') trim(name), '_', k, '.vtu'
          writer = paraview_writer(points, ragged(corners, count), kinds)
          call writer % write_binary(trim(file), body, [string('body')])
          files = [character(len=256) :: files, file]
          times = [times, on % instant(k)]
       end do
    class default
       error stop 'multibody_dynamics: a mechanism is drawn on a solution of its manifold'
    end select
    call time_series(trim(name) // '.pvd', files, times)

  end subroutine draw

  ! the rotation matrix of the Euler parameters, C = I - 2 eta eps^x +
  ! 2 eps^x eps^x
  pure function rotation(eta, eps) result(c)
    real(dp), intent(in) :: eta, eps(3)
    real(dp) :: c(3,3)
    c(1, :) = [1.0_dp - 2.0_dp * (eps(2)**2 + eps(3)**2), 2.0_dp * (eps(1) * eps(2) + eta * eps(3)), &
         & 2.0_dp * (eps(1) * eps(3) - eta * eps(2))]
    c(2, :) = [2.0_dp * (eps(2) * eps(1) - eta * eps(3)), 1.0_dp - 2.0_dp * (eps(1)**2 + eps(3)**2), &
         & 2.0_dp * (eps(2) * eps(3) + eta * eps(1))]
    c(3, :) = [2.0_dp * (eps(3) * eps(1) + eta * eps(2)), 2.0_dp * (eps(3) * eps(2) - eta * eps(1)), &
         & 1.0_dp - 2.0_dp * (eps(1)**2 + eps(2)**2)]
  end function rotation

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

  function spherical(a, x, b) result(g)
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
  end function spherical

  function revolute(a, x, e, b) result(g)
    type(rigid_body), intent(in) :: a, b
    real(dp)        , intent(in) :: x(3), e(3)
    type(continuous_field) :: g(5)
    type(continuous_field) :: ea(3)
    real(dp) :: e1(3), e2(3)
    call normals(e, e1, e2)
    ea = a % direction(e)
    g(1:3) = spherical(a, x, b)
    g(4)   = dot(ea, b % direction(e1))
    g(5)   = dot(ea, b % direction(e2))
  end function revolute

  function cylindrical(a, x, e) result(g)
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
  end function cylindrical

  function fixed(a) result(g)
    type(rigid_body), intent(in) :: a
    type(continuous_field) :: g(6)
    g = [a % u, a % eps]
  end function fixed

  function rigid_link(a, b) result(g)
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
  end function rigid_link

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
