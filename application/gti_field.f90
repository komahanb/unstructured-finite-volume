! THE FIELD: what a field adds to a march over a spatial mesh, and
! what is checked about it.
!
! A block over a field holds every node's components at every moment,
! node by node within a moment, laid out by the constructors in
! gti_march and gti_stage exactly as a single node's block is with
! nodes = 1. What the field adds is the level below: a stencil over
! the nodes carrying minus the framework's diffusion operator, each
! row divided by its cell's area so that the flux balance becomes
! kappa times the laplacian. A block lays it on every moment it
! evaluates its physics at, where it adds to the physics' row. It is
! linear in the state and knows nothing of the design, so it enters
! the apply and the tangent and nothing else. Its order is the form's,
! and the form's degree is given.
!
! The functional over a field is the integral over the domain and the
! duration, so the measure of a node is its cell's area, and the
! chain weights each point by the step times that area.
!
! Three things are checked, each against something known:
!
!   the operator  the balance of the rectangle's mode against kappa
!                 times its laplacian, cell by cell, walls apart
!   the mode      at nu = 0 on a rectangle the field q = cos(pi x / a)
!                 cos(pi y / b) cos(omega t) is exact, with omega^2 =
!                 1 + kappa pi^2 (1/a^2 + 1/b^2), so the last instant
!                 is measured against it, and against the semi-discrete
!                 solution that isolates the time error
!   the ode       at kappa = 0 with a constant field every node is one
!                 node's equation - checked by the program, which
!                 marches the node
!
! and any instant may be written as a vtu file for paraview.
module gti_field

  use util_precision   , only : dp
  use operation_stencil, only : stencil
  use field_calculus   , only : field
  use field_stored     , only : stored_field
  use operation_expression, only : expression
  use gti_configuration, only : worded
  use gti_march        , only : consistent_states
  use gti_space        , only : room, spatial_operator, cartesian, written_paraview

  implicit none

  private
  public :: node_operator, initial_field
  public :: against_the_laplacian, against_the_mode, export_instant

contains

  !-------------------------------------------------------------------!
  ! The level below as a stencil over the nodes: -kappa times the
  ! laplacian, one row per cell. A wall holding a value would enter
  ! as a source, which a block has no place for; the wall here holds
  ! no flux, and the block refuses a constant if one arrives.
  !-------------------------------------------------------------------!

  function node_operator(space, kappa, degree) result(op)

    type(room), intent(in) :: space
    real(dp)  , intent(in) :: kappa
    integer   , intent(in) :: degree
    type(stencil) :: op

    type(stencil) :: balance
    integer , allocatable :: r(:), c(:)
    real(dp), allocatable :: lw(:), w(:), held(:)
    integer :: e, m

    balance = spatial_operator(space, kappa, degree)
    m       = balance % pattern % num_edges()
    call balance % weights % real_vector(lw)
    call balance % constants % real_vector(held)

    allocate(r(m), c(m), w(m))
    do e = 1, m
       r(e) = balance % pattern % edge_head(e)
       c(e) = balance % pattern % edge_tail(e)
       w(e) = -lw(e) / space % volume(r(e))
    end do

    op = stencil(r, c, w, held, 'level below')

  end function node_operator

  !-------------------------------------------------------------------!
  ! The state at the first instant: the components below the highest
  ! at every node - constant from the words given, or the rectangle's
  ! mode, or one plus half of it - and the highest solved from the
  ! physics with the level below laid on, so that the state is
  ! consistent with the equation rather than merely plausible. One
  ! node with no mesh is one node's equation. Invalid input: more
  ! components than lie below the highest; a mode with no rectangle.
  !-------------------------------------------------------------------!

  function initial_field(physics, degrees, kind, initial_state, design, spatial, space, a, b) &
       & result(q)

    type(expression)      , intent(in)           :: physics
    integer               , intent(in)           :: degrees
    character(len=*)      , intent(in)           :: kind, initial_state
    real(dp)              , intent(in)           :: design
    type(stencil)         , intent(in), optional :: spatial
    type(room)            , intent(in), optional :: space
    real(dp)              , intent(in), optional :: a, b
    real(dp), allocatable :: q(:)

    character(len=32), allocatable :: given(:)
    real(dp), allocatable :: lower(:,:)
    integer  :: nodes, i, d

    nodes = 1
    if (present(space)) nodes = space % num_cells
    allocate(lower(degrees - 1, nodes), source=0.0_dp)

    select case (trim(kind))
    case ('constant')
       given = worded(initial_state)
       if (size(given) > degrees - 1) then
          write(*,'(a,i0,a)') ' the initial state holds the ', degrees - 1, &
               & ' components below the highest, which the physics gives.'
          error stop 'gti_field: the initial state is given below the highest derivative'
       end if
       do d = 1, size(given)
          read(given(d), *) lower(d, 1)
       end do
       do i = 2, nodes
          lower(:, i) = lower(:, 1)
       end do
    case ('mode')
       if (.not. present(space)) error stop 'gti_field: the mode is a field over a mesh'
       if (space % geometry /= cartesian) error stop 'gti_field: the mode is the rectangle''s'
       lower(1, :) = mode_shape(space, a, b)
    case ('bump')
       ! one plus half the rectangle's mode, on any geometry: a field
       ! that is not uniform, so the level below has something to do
       if (.not. present(space)) error stop 'gti_field: the bump is a field over a mesh'
       lower(1, :) = 1.0_dp + 0.5_dp * mode_shape(space, a, b)
    case default
       error stop 'gti_field: an initial field is constant, the mode, or the bump'
    end select

    q = consistent_states(physics, degrees, lower, design, spatial)

  end function initial_field

  !-------------------------------------------------------------------!
  ! The rectangle's mode at every cell centre, cos(pi x / a) cos(pi y
  ! / b): the shape every check on the rectangle reads.
  !-------------------------------------------------------------------!

  pure function mode_shape(space, a, b) result(shape)

    type(room), intent(in) :: space
    real(dp)  , intent(in) :: a, b
    real(dp), allocatable :: shape(:)

    real(dp) :: pi
    integer  :: i

    pi = acos(-1.0_dp)
    shape = [(cos(pi * space % centre(1, i) / a) * cos(pi * space % centre(2, i) / b), &
         &    i = 1, space % num_cells)]

  end function mode_shape

  !-------------------------------------------------------------------!
  ! The operator applied to a field over the cells: the integrated
  ! flux balance of that field.
  !-------------------------------------------------------------------!

  subroutine balance_of(space, kappa, degree, values, balanced)

    type(room), intent(in) :: space
    real(dp)  , intent(in) :: kappa, values(:)
    integer   , intent(in) :: degree
    real(dp), allocatable, intent(out) :: balanced(:)

    type(stencil) :: op
    type(stored_field) :: given
    class(field), allocatable :: out

    op    = spatial_operator(space, kappa, degree)
    given = stored_field('values', op % pattern % vertex_set(), size(values))
    call given % set_real_vector(values)
    call op % apply(op % pattern, [given], out)
    call out % real_vector(balanced)

  end subroutine balance_of

  !-------------------------------------------------------------------!
  ! The operator alone against the laplacian of the mode, cell by
  ! cell: the balance over the area against kappa times minus pi^2
  ! (1/a^2 + 1/b^2) times the mode, which has no normal derivative at
  ! any wall. The error is reported over the cells that touch no
  ! wall, those that touch one, and those that touch two, so a wall's
  ! treatment is told apart from the interior's.
  !-------------------------------------------------------------------!

  subroutine against_the_laplacian(space, a, b, kappa, degree)

    type(room), intent(in) :: space
    real(dp)  , intent(in) :: a, b, kappa
    integer   , intent(in) :: degree

    real(dp), allocatable :: shape(:), balanced(:), exact(:)
    real(dp) :: pi, err(0:2), norm(0:2)
    integer  :: i, walls, count(0:2)

    if (space % geometry /= cartesian) then
       error stop 'gti_field: the laplacian check is the rectangle''s'
    end if

    pi    = acos(-1.0_dp)
    shape = mode_shape(space, a, b)
    exact = -kappa * pi ** 2 * (1.0_dp / a ** 2 + 1.0_dp / b ** 2) * shape
    call balance_of(space, kappa, degree, shape, balanced)

    err   = 0.0_dp
    norm  = 0.0_dp
    count = 0
    do i = 1, space % num_cells
       walls = 0
       if (space % cell_ij(1, i) == 1 .or. space % cell_ij(1, i) == space % n2) walls = walls + 1
       if (space % cell_ij(2, i) == 1 .or. space % cell_ij(2, i) == space % n1) walls = walls + 1
       err(walls)   = err(walls)   + (balanced(i) / space % volume(i) - exact(i)) ** 2
       norm(walls)  = norm(walls)  + exact(i) ** 2
       count(walls) = count(walls) + 1
    end do

    write(*,'(a,i0,a,i0,a)') '   the operator against kappa laplacian of the mode, ', &
         & space % num_cells, ' cells, form degree ', degree, ':'
    write(*,'(a,3(a,es10.3))') '   relative rms error', &
         & '   interior ', sqrt(err(0) / max(norm(0), tiny(1.0_dp))), &
         & '   one wall ', sqrt(err(1) / max(norm(1), tiny(1.0_dp))), &
         & '   corner ',   sqrt(err(2) / max(norm(2), tiny(1.0_dp)))

  end subroutine against_the_laplacian

  !-------------------------------------------------------------------!
  ! At nu = 0 on a rectangle, the last instant against the exact mode
  ! and against the semi-discrete mode, whose frequency carries the
  ! operator's eigenvalue as built, by the rayleigh quotient of the
  ! mode. The first error holds space and time, the second time
  ! alone. The instant's components arrive node by node, degrees
  ! within a node. Nothing is said on any other geometry or design.
  !-------------------------------------------------------------------!

  subroutine against_the_mode(space, a, b, kappa, degree, design, t_last, x, degrees)

    type(room), intent(in) :: space
    real(dp)  , intent(in) :: a, b, kappa, design, t_last, x(:)
    integer   , intent(in) :: degree, degrees

    real(dp) :: pi, omega, omega_h, exact, semi, e_exact, e_semi, area, mode
    real(dp), allocatable :: shape(:), balanced(:)
    integer  :: i

    if (space % geometry /= cartesian .or. design /= 0.0_dp) return

    pi    = acos(-1.0_dp)
    omega = sqrt(1.0_dp + kappa * pi ** 2 * (1.0_dp / a ** 2 + 1.0_dp / b ** 2))

    ! minus the mode against its own balance, over the mode against
    ! itself by area
    shape = mode_shape(space, a, b)
    call balance_of(space, kappa, degree, shape, balanced)
    omega_h = sqrt(1.0_dp - dot_product(shape, balanced) / &
         & dot_product(shape, space % volume * shape))

    e_exact = 0.0_dp
    e_semi  = 0.0_dp
    area    = sum(space % volume)

    do i = 1, space % num_cells
       mode  = shape(i)
       exact = mode * cos(omega   * t_last)
       semi  = mode * cos(omega_h * t_last)
       e_exact = e_exact + space % volume(i) * (x((i - 1) * degrees + 1) - exact) ** 2
       e_semi  = e_semi  + space % volume(i) * (x((i - 1) * degrees + 1) - semi) ** 2
    end do

    write(*,'(a,es12.3,a,es12.3,a,f10.6,a,f10.6)') &
         & '      error at the last instant, against the mode ', sqrt(e_exact / area), &
         & '   semi-discrete ', sqrt(e_semi / area), '   omega ', omega, '   omega_h ', omega_h

  end subroutine against_the_mode

  !-------------------------------------------------------------------!
  ! One instant as one vtu file: every degree of every node, the
  ! components arriving node by node, degrees within a node.
  !-------------------------------------------------------------------!

  subroutine export_instant(space, path, degrees, x)

    type(room)      , intent(in) :: space
    character(len=*), intent(in) :: path
    integer         , intent(in) :: degrees
    real(dp)        , intent(in) :: x(:)

    character(len=8), allocatable :: names(:)
    real(dp), allocatable :: values(:,:)
    integer :: i, d

    allocate(names(degrees), values(space % num_cells, degrees))
    do d = 0, degrees - 1
       write(names(d + 1),'(a,i0)') 'q', d
    end do
    do i = 1, space % num_cells
       do d = 0, degrees - 1
          values(i, d + 1) = x((i - 1) * degrees + d + 1)
       end do
    end do

    call written_paraview(space, path, names, values)

  end subroutine export_instant

end module gti_field
