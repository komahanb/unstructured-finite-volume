!=====================================================================!
! THE CASE of HART II for the template flexible_mechanism
! (flexible_mechanism.f90): the rotor from rest on the ground, its
! assembly, its controls, its breakpoints, the bounds of its design and
! its limit.
!
! THE ROTOR, named as in rotor_bodies.txt: bodies [ ], blades ( ),
! joints on the edges, # the ground.
!
!                azimuth                       collective, longitudinal, lateral
!   # --support-- [hub] <==driver== #          [push_rod_90, _180, _270] <==driver== #
!                   |                                |      |
!               revolute                         spherical  support (line)
!                   |                                |      |
!            [upper_link_30 ...] x4          [lower_swashplate] --spherical-- [sphere] --cylindrical-- #
!              |           |                    |            |
!            rigid     spherical             revolute    spherical
!              |           |                    |            |
!          [cap_0 ...] x4  [lower_link_30 ...] x4 |     [upper_horn] --revolute-- [lower_horn]
!              |           |                    |                                      |
!          attachment  spherical                |                                  revolute
!              |           |                    |                                      |
!         (blade_0 ...) x4 +--------- [upper_swashplate]              # --fixed-- [base_plate]
!
! THE CONTROLS u = (azimuth, collective, longitudinal, lateral):
!
!   the hub's rotation   agrees with the turn by the azimuth: c = sin((theta
!                        - azimuth) / 2), the shaft component of the vector
!                        part of the rotation between them, at every turn
!   push rod k's         displacement along the shaft = collective +
!                        longitudinal x_k + lateral y_k, the height at the
!                        rod's point of the plane its slopes define - linear,
!                        the rows of A
!   efforts              e = -(dc/du)^T mu: the shaft torque, (mu / 2)
!                        cos((theta - azimuth) / 2), mu / 2 at agreement while
!                        the hub's Euler parameters keep the representative they
!                        start with - the antipodal one reverses its sign; the
!                        collective force and the generalized efforts conjugate
!                        to the slopes, the rods' rows of A^T mu
!
! THE BREAKPOINTS: 0; the startup, 1 s, the rate speed s(t / startup)
! from rest, s(z) = z^3 (10 - 15 z + 6 z^2), its first and second
! derivatives zero at 0 and 1; the onset, a revolution, the strokes
! from zero to the scenario's by the same s; the observation, a
! revolution, the strokes held. TACS's three scenarios raise the push
! rods at 90, 180 and 270 degrees by (0.05, 0.05, 0.05), (0.025,
! 0.025, 0.05) and (0.025, 0.05, 0.025) m.
!
! THE DESIGN: the thickness of each of eight strips of the span, 2 cm
! at the start, within 1 to 2 cm, adjacent strips within 6 mm. A blade
! is a beam of its span at nine points, its section the solid rectangle
! of the thickness h(x) on the strips, its pitch HART II's twist,
! attached at x = 0 to its cap in the cap's axes coned by HART II's
! precone. The limit of the ratio of the stress to the strength: 1/4.
!
! At rest nine directions are determined neither by the joints nor by
! the loads, measured as the kernel at the reference: the sphere's spin
! about the shaft, the four pitch links' and the three push rods' spins
! about their axes, and the upper swashplate's turn about the shaft,
! whose pitch links, vertical there, it moves across their lengths, a
! mechanism to the first order; the template's reference phase
! condition holds them at the reference. At the droop the kernel is of
! four, the weights the joints carry stiffening the rest.
!=====================================================================!

module hart2_case

  use util_precision    , only : dp
  use operation_manifold, only : region, interval, parameter
  use operation_field   , only : continuous_field, indicators, heaviside
  use multibody_dynamics, only : assembly, table, beam, solid_rectangle, attachment, revolute, spherical, &
       &                         cylindrical, fixed, rigid, line, rotation, displacement

  implicit none

  private
  public :: machine, controls, bounds, scenarios, breakpoints, observation, limit

  ! HART II and TACS's case: the rate and its revolution; the blade's
  ! root, tip, chord and height; HART II's precone and twist
  real(dp), parameter :: pi = acos(-1.0_dp)
  real(dp), parameter :: speed = 109.12_dp, revolution = 2.0_dp * pi / speed
  real(dp), parameter :: root = 0.44_dp, tip = 2.0_dp, chord = 0.121_dp, height = 0.9_dp
  real(dp), parameter :: precone = 2.5_dp * pi / 180.0_dp, twist = -8.0_dp * pi / 180.0_dp, twist_zero = 1.5_dp
  real(dp), parameter :: ex(3) = [1.0_dp, 0.0_dp, 0.0_dp], ey(3) = [0.0_dp, 1.0_dp, 0.0_dp], &
       &                 ez(3) = [0.0_dp, 0.0_dp, 1.0_dp]

  ! the breakpoints: the start; the end of the startup; the end of the
  ! onset, a revolution later; the end of the observation, a revolution
  ! after that
  real(dp), parameter :: startup = 1.0_dp, onset = revolution
  real(dp), parameter :: breakpoints(4) = [0.0_dp, startup, startup + onset, startup + onset + revolution]
  real(dp), parameter :: observation(2) = breakpoints(3:4)

  ! the design: eight strips of the span, each 2 cm at the start; the
  ! blades' material; gravity; the limit
  real(dp), parameter :: edges(9) = (tip - root) * [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, &
       &                                           7.0_dp, 8.0_dp] / 8.0_dp
  real(dp), parameter :: start(8) = 0.02_dp
  real(dp), parameter :: density = 2500.0_dp, modulus = 70.0e9_dp, poisson = 0.3_dp, strength = 350.0e6_dp
  real(dp), parameter :: gravity(3) = [0.0_dp, 0.0_dp, -9.81_dp]
  real(dp), parameter :: limit = 0.25_dp

  ! the scenarios' strokes - collective, m, longitudinal and lateral,
  ! slopes - the height and the slopes of the plane through the push
  ! rods' points at TACS's three raises
  integer , parameter :: scenarios = 3
  real(dp), parameter :: radius = 0.4004_dp
  real(dp), parameter :: strokes(3, scenarios) = reshape([0.05_dp, 0.0_dp, 0.0_dp, &
       &                                                  0.0375_dp, 0.0125_dp / radius, -0.0125_dp / radius, &
       &                                                  0.025_dp, -0.025_dp / radius, 0.0_dp], [3, scenarios])

contains

  ! THE ASSEMBLY: the bodies of TACS's table under gravity and the
  ! design; the blades; the supports; the passive joints; the drivers
  function machine() result(rotor)

    type(assembly) :: rotor

    type(continuous_field) :: x
    type(continuous_field), allocatable :: h(:)
    character(len=:), allocatable :: blade, cap, lower, upper, rod
    real(dp) :: point(3)
    integer :: k

    rotor = assembly(table('rotor_bodies.txt'), gravity, design=parameter('h', start))
    h     = rotor % coordinate('h')

    ! the blades: each a member on a span of its own at nine points,
    ! attached at x = 0 to its cap at the root point of its azimuth, in
    ! the cap's axes coned by the precone; a beam whose section is the
    ! solid rectangle of the thickness h(x), its pitch the twist
    do k = 0, 3
       blade = 'blade_' // degrees(k)
       call rotor % member(blade, region(interval(0.0_dp, tip - root)), points=9, &
            &              attachment=attachment('cap_' // degrees(k), turned([root, 0.0_dp, height], k), axes(k)))
       x = rotor % coordinate(blade, 'x')
       call rotor % model(blade, beam(solid_rectangle(density, modulus, poisson, strength, chord, &
            &                          thickness=continuous_field(coefficients=h, basis=indicators(x, edges))), &
            &                     pitch=(twist / tip) * (x + root - twist_zero)))
    end do

    ! the supports: the hub's centre held, its axis the shaft; each push
    ! rod's centre on a line along the shaft
    call rotor % joint('hub_support', revolute('hub', [0.0_dp, 0.0_dp, 0.0_dp], ez))
    do k = 1, 3
       call rotor % joint('rod_support_' // degrees(k), line('push_rod_' // degrees(k), ez))
    end do

    ! the swashplates hinged about the shaft, the lower one on the
    ! sphere, which slides along the shaft
    call rotor % joint('swashplates', revolute('lower_swashplate', [0.0_dp, 0.0_dp, 0.4_dp], ez, 'upper_swashplate'))
    call rotor % joint('swashplate_sphere', spherical('lower_swashplate', [0.0_dp, 0.0_dp, 0.4_dp], 'sphere'))
    call rotor % joint('sphere_shaft', cylindrical('sphere', [0.0_dp, 0.0_dp, 0.4_dp], ez))

    ! at each azimuth: the upper swashplate to the lower pitch link, the
    ! lower to the upper pitch link, the upper hinged to the hub, the cap
    ! fixed to it
    do k = 0, 3
       cap   = 'cap_' // degrees(k)
       lower = 'lower_link_' // degrees(k, 30)
       upper = 'upper_link_' // degrees(k, 30)
       call rotor % joint('swashplate_link_' // degrees(k), &
            & spherical('upper_swashplate', turned([0.278167359696_dp, 0.1606_dp, 0.43_dp], k), lower))
       call rotor % joint('links_' // degrees(k), &
            & spherical(lower, turned([0.278167359696_dp, 0.1606_dp, height], k), upper))
       call rotor % joint('shaft_link_' // degrees(k), &
            & revolute('hub', turned([0.231_dp, 0.0_dp, height], k), abs(turned(ex, k)), upper))
       call rotor % joint('link_cap_' // degrees(k), rigid(upper, cap))
    end do

    ! the push rods on the lower swashplate
    do k = 1, 3
       rod = 'push_rod_' // degrees(k)
       call rotor % joint('rod_swashplate_' // degrees(k), &
            & spherical(rod, turned([radius, 0.0_dp, 0.38_dp], k), 'lower_swashplate'))
    end do

    ! the base plate fixed, the horns hinged to it and to each other about
    ! y, the upper horn on the lower swashplate
    call rotor % joint('base', fixed('base_plate'))
    call rotor % joint('base_horn', revolute('lower_horn', [radius, 0.0_dp, 0.05_dp], ey, 'base_plate'))
    call rotor % joint('horns', revolute('lower_horn', [0.5654_dp, 0.0_dp, 0.215_dp], ey, 'upper_horn'))
    call rotor % joint('horn_swashplate', spherical('upper_horn', [radius, 0.0_dp, 0.38_dp], 'lower_swashplate'))

    ! the drivers on (azimuth, collective, longitudinal, lateral), each by
    ! its row: the hub's rotation about the shaft agreeing with the turn by
    ! the azimuth, sin((theta - azimuth) / 2) = 0; each push rod's
    ! displacement along the shaft the height of the plane at its point
    ! (x_k, y_k), linear
    call rotor % driver('hub', rotation('hub', ez), row=[1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp])
    do k = 1, 3
       rod   = 'push_rod_' // degrees(k)
       point = turned([radius, 0.0_dp, 0.0_dp], k)
       call rotor % driver(rod, displacement(rod, ez), row=[0.0_dp, 1.0_dp, point(1), point(2)])
    end do

  end function machine

  ! THE CONTROLS of scenario j, (azimuth, collective, longitudinal,
  ! lateral): the azimuth from rest to the rate over the startup, then
  ! the rate held; the strokes from zero to the scenario's over the
  ! onset, then held; their rates and accelerations zero at t = 0
  function controls(t, j) result(u)
    type(continuous_field), intent(in) :: t
    integer               , intent(in) :: j
    type(continuous_field) :: u(4)
    type(continuous_field) :: rise
    rise = ramp(t, startup, onset)
    u    = [azimuth(t), strokes(1, j) * rise, strokes(2, j) * rise, strokes(3, j) * rise]
  end function controls

  ! the azimuth: speed startup (5 z^4 / 2 - 3 z^5 + z^6) at z = t /
  ! startup, the integral of the rate over the startup; then the rate
  ! held
  function azimuth(t) result(psi)
    type(continuous_field), intent(in) :: t
    type(continuous_field) :: psi, z
    z   = t / startup
    psi = speed * startup * (2.5_dp * z**4 - 3.0_dp * z**5 + z**6) * (1.0_dp - heaviside(t - startup)) &
         & + speed * (t - startup / 2.0_dp) * heaviside(t - startup)
  end function azimuth

  ! zero before a, s((t - a) / tau) from a to a + tau, one after it
  function ramp(t, a, tau) result(s)
    type(continuous_field), intent(in) :: t
    real(dp)              , intent(in) :: a, tau
    type(continuous_field) :: s, z
    z = (t - a) / tau
    s = z**3 * (10.0_dp - 15.0_dp * z + 6.0_dp * z**2) * (heaviside(t - a) - heaviside(t - a - tau)) &
         & + heaviside(t - a - tau)
  end function ramp

  ! the bounds of the design: each thickness within 1 to 2 cm, adjacent
  ! strips within 6 mm, each at most zero
  function bounds(h) result(c)
    type(continuous_field), intent(in) :: h(:)
    type(continuous_field), allocatable :: c(:)
    integer :: n
    n = size(h)
    c = [0.01_dp - h, h - 0.02_dp, h(:n - 1) - h(2:) - 0.006_dp, h(2:) - h(:n - 1) - 0.006_dp]
  end function bounds

  ! a point turned by k quarter turns about the shaft
  function turned(p, k) result(q)
    real(dp), intent(in) :: p(3)
    integer , intent(in) :: k
    real(dp) :: q(3)
    integer  :: i
    q = p
    do i = 1, modulo(k, 4)
       q = [-q(2), q(1), q(3)]
    end do
  end function turned

  ! the axes of the blade of k quarter turns, its span coned up by the
  ! precone
  function axes(k) result(a)
    integer, intent(in) :: k
    real(dp) :: a(3, 3)
    a = reshape([turned([cos(precone), 0.0_dp, sin(precone)], k), turned(ey, k), &
         &       turned([-sin(precone), 0.0_dp, cos(precone)], k)], [3, 3])
  end function axes

  ! the name of the azimuth of k quarter turns, in degrees, offset by
  ! the degrees given
  function degrees(k, offset) result(name)
    integer, intent(in) :: k
    integer, intent(in), optional :: offset
    character(len=:), allocatable :: name
    character(len=8) :: buffer
    integer :: angle_k
    angle_k = 90 * k
    if (present(offset)) angle_k = angle_k + offset
    write(buffer, '(i0)') angle_k
    name = trim(buffer)
  end function degrees

end module hart2_case
