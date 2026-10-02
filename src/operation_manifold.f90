!=====================================================================!
! THE MANIFOLD, CONTINUOUS AND DISCRETE.
!
! A continuous manifold is a product of coordinates with their
! extents: a time interval, a region of space, either, or neither,
! in which case the manifold is a point. Its coordinate functions
! and its unknown functions are fields on it; a part of it - the
! face at one instant, the time factor, the space factor - is a
! manifold of its own that records its parent.
!
! Its discrete image is the product of the points chosen along each
! coordinate: the instants of the interval, the cells of a mesh of
! the region. Nothing else: the approximation of a derivative is a
! property of the residual placed on the points, not of the points.
! The image states which coordinates it has by what it stores - its
! instants, its cells, the expansion of its design - and repeats no
! declaration of the manifold it is the image of.
!
!      coordinate index   1 the time, then the space axes in order
!      point index        instant outer, cell inner
!      position(p)        the instant and the cell centre
!      measure(p)         the instant's share of the interval times the cell's volume
!
! THE DESIGN. The design coordinate nu is a point of the design
! space, fixed at a value nu_0; its discrete image is the expansion
! of an order at that value, which operation_expansion states, so
! that the discrete manifold is the instants (and cells) times the
! orders of the expansion.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_manifold

  use util_precision      , only : dp
  use token_identity      , only : token, next_token
  use operation_expansion , only : expansion
  use operation_field     , only : continuous_support, discrete_support, continuous_field
  use operation_field     , only : unknown_field, coordinate_field
  use operation_field     , only : WHOLE, TIME_FACE, TIME_FACTOR, SPACE_FACTOR, DESIGN_FACTOR, SPACE_FACE
  use view_mesh           , only : mesh, values_of
  use view_mesh_builder   , only : mesh_from_gmsh, mesh_from_points, mesh_of_corners

  implicit none

  private
  public :: continuous_manifold, discrete_manifold
  public :: interval, region, points, mesh, corners, parameter

  integer, parameter :: MAX_NAME = 32

  type :: interval
     real(dp) :: a = 0.0_dp, b = 0.0_dp
  end type interval

  interface interval
     module procedure create_interval
  end interface interval

  ! the region of space: as a gmsh geometry file describes it, or
  ! the interval of one space coordinate
  type :: region
     character(len=:), allocatable :: geometry
     type(interval), allocatable :: extent
  end type region

  interface region
     module procedure create_region
     module procedure create_region_of_interval
  end interface region

  ! the points of an interval of time or of a space coordinate: n
  ! equally spaced, or the list given, ascending from one end of the
  ! interval to the other
  type :: points
     integer :: n = 0
     real(dp), allocatable :: position(:)
  end type points

  interface points
     module procedure create_points
     module procedure create_points_at
  end interface points

  interface mesh
     module procedure mesh_from_file
  end interface mesh

  ! the corners of a mesh as the points of a region, the mesh's cells
  ! the elements of their corners
  interface corners
     module procedure corners_of
  end interface corners

  ! the design: a coordinate fixed at a value, a point of the design
  ! space; its discrete image is the expansion of an order at that
  ! value, so that the discrete solution stores its derivatives along
  ! the design to the order
  type :: parameter
     character(len=MAX_NAME) :: name = ''
     real(dp) :: value = 0.0_dp
  end type parameter

  interface parameter
     module procedure create_parameter
  end interface parameter

  type, extends(continuous_support) :: continuous_manifold

     logical        :: with_time  = .false.
     logical        :: with_space = .false.
     logical        :: with_design = .false.
     type(interval) :: span
     type(region)   :: geometry
     type(parameter), allocatable :: design_parameter(:)
     integer        :: num_unknowns = 0
     character(len=MAX_NAME), allocatable :: unknown_name(:)

   contains

     procedure :: unknown         => manifold_unknown
     procedure :: coordinate      => manifold_coordinate
     procedure :: boundary        => manifold_boundary
     procedure :: time            => manifold_time
     procedure :: space           => manifold_space
     procedure :: design          => manifold_design
     procedure :: direction       => manifold_direction
     procedure :: discretize      => manifold_discretize
     procedure :: num_coordinates => manifold_num_coordinates
     procedure :: is_point
     procedure :: is_time_face
     procedure :: is_space_face
     procedure :: is_time_factor
     procedure :: is_space_factor
     procedure :: is_design_factor

  end type continuous_manifold

  interface continuous_manifold
     module procedure create_manifold
     module procedure create_manifold_of_designs
  end interface continuous_manifold

  type, extends(discrete_support) :: discrete_manifold

     real(dp), allocatable :: instant(:)
     type(mesh) :: cells
     integer :: dimension = 0
     real(dp), allocatable :: centre(:)
     real(dp), allocatable :: volume(:)

   contains

     procedure :: num_points      => discrete_num_points
     procedure :: num_coordinates => discrete_num_coordinates
     procedure :: position        => discrete_position
     procedure :: measure         => discrete_measure
     procedure :: design_factor   => discrete_design_factor
     procedure :: measure_on_face => discrete_measure_on_face
     procedure :: num_instants
     procedure :: num_cells
     procedure :: step
     procedure :: point_of

  end type discrete_manifold

contains

  !===================================================================!
  ! THE CONSTRUCTORS. An interval with b below a, a count of instants
  ! below two, stops the program.
  !===================================================================!

  function create_interval(a, b) result(this)

    real(dp), intent(in) :: a, b
    type(interval) :: this

    character(len=250) :: message

    if (b <= a) then
       write(message,'(a,es12.5,a,es12.5)') 'operation_manifold: an interval requires b > a; a = ', &
            & a, ', b = ', b
       error stop trim(message)
    end if
    this % a = a
    this % b = b

  end function create_interval

  function create_region(geometry) result(this)

    character(len=*), intent(in) :: geometry
    type(region) :: this

    this % geometry = geometry

  end function create_region

  function create_region_of_interval(extent) result(this)

    type(interval), intent(in) :: extent
    type(region) :: this

    this % extent = extent

  end function create_region_of_interval

  function create_points(n) result(this)

    integer, intent(in) :: n
    type(points) :: this

    character(len=250) :: message

    if (n < 2) then
       write(message,'(a,i0)') 'operation_manifold: an interval is discretized by two points at &
            &least; n = ', n
       error stop trim(message)
    end if
    this % n = n

  end function create_points

  function create_points_at(position) result(this)

    real(dp), intent(in) :: position(:)
    type(points) :: this

    integer :: k
    character(len=250) :: message

    if (size(position) < 2) then
       write(message,'(a,i0)') 'operation_manifold: an interval is discretized by two points at &
            &least; size(position) = ', size(position)
       error stop trim(message)
    end if
    do k = 2, size(position)
       if (position(k) <= position(k - 1)) then
          write(message,'(a,i0)') 'operation_manifold: the points ascend; point ', k
          error stop trim(message) // ' does not follow the one before it'
       end if
    end do
    this % n    = size(position)
    this % position = position

  end function create_points_at

  function mesh_from_file(filename) result(cells)

    character(len=*), intent(in) :: filename
    type(mesh) :: cells

    cells = mesh_from_gmsh(filename)

  end function mesh_from_file

  function corners_of(cells) result(points)

    type(mesh), intent(in) :: cells
    type(mesh) :: points

    points = mesh_of_corners(cells)

  end function corners_of

  ! a design parameter named as a coordinate is, other than t, x, y, z
  function create_parameter(name, value) result(this)

    character(len=*), intent(in) :: name
    real(dp)        , intent(in) :: value
    type(parameter) :: this

    select case (trim(name))
    case ('t', 'x', 'y', 'z')
       error stop 'operation_manifold: the design coordinate is named other than t, x, y, z; the name given is ' &
            & // trim(name)
    end select
    if (len_trim(name) < 1 .or. len_trim(name) > MAX_NAME) then
       error stop 'operation_manifold: the design coordinate is named by one to 32 characters'
    end if
    this % name  = name
    this % value = value

  end function create_parameter

  function create_manifold(time, space, design) result(this)

    type(interval) , intent(in), optional :: time
    type(region)   , intent(in), optional :: space
    type(parameter), intent(in), optional :: design
    type(continuous_manifold) :: this

    this % identity = next_token()
    if (present(time)) then
       this % with_time = .true.
       this % span      = time
    end if
    if (present(space)) then
       this % with_space = .true.
       this % geometry   = space
    end if
    if (present(design)) then
       this % with_design = .true.
       this % design_parameter = [design]
    end if
    allocate(this % unknown_name(0))

  end function create_manifold

  ! the same with several design coordinates, each named once
  function create_manifold_of_designs(time, space, design) result(this)

    type(interval) , intent(in), optional :: time
    type(region)   , intent(in), optional :: space
    type(parameter), intent(in) :: design(:)
    type(continuous_manifold) :: this

    integer :: i, j

    this = create_manifold(time, space)
    if (size(design) < 1) then
       error stop 'operation_manifold: a manifold of designs names one design coordinate at least'
    end if
    do i = 2, size(design)
       do j = 1, i - 1
          if (trim(design(i) % name) == trim(design(j) % name)) then
             error stop 'operation_manifold: a design coordinate is named once; the name repeated is ' &
                  & // trim(design(i) % name)
          end if
       end do
    end do
    this % with_design = .true.
    this % design_parameter = design

  end function create_manifold_of_designs

  pure logical function is_point(this)
    class(continuous_manifold), intent(in) :: this
    is_point = .not. (this % with_time .or. this % with_space)
  end function is_point

  pure logical function is_time_face(this)
    class(continuous_manifold), intent(in) :: this
    is_time_face = this % part == TIME_FACE
  end function is_time_face

  pure logical function is_space_face(this)
    class(continuous_manifold), intent(in) :: this
    is_space_face = this % part == SPACE_FACE
  end function is_space_face

  pure logical function is_time_factor(this)
    class(continuous_manifold), intent(in) :: this
    is_time_factor = this % part == TIME_FACTOR
  end function is_time_factor

  pure logical function is_space_factor(this)
    class(continuous_manifold), intent(in) :: this
    is_space_factor = this % part == SPACE_FACTOR
  end function is_space_factor

  pure logical function is_design_factor(this)
    class(continuous_manifold), intent(in) :: this
    is_design_factor = this % part == DESIGN_FACTOR
  end function is_design_factor

  !===================================================================!
  ! The coordinates: the time first when present, then the space
  ! axes. Their number is known once the space is discretized; before
  ! that a name is enough to number them.
  !===================================================================!

  pure integer function manifold_num_coordinates(this)
    class(continuous_manifold), intent(in) :: this
    manifold_num_coordinates = merge(1, 0, this % with_time)
    if (this % with_space) manifold_num_coordinates = manifold_num_coordinates &
         & + merge(1, 3, allocated(this % geometry % extent))
  end function manifold_num_coordinates

  !===================================================================!
  ! An unknown function of the manifold, numbered after those already
  ! declared, one field per component, a function of the coordinates
  ! given as its arguments: a function of every coordinate of the
  ! manifold when none are given, of those alone otherwise, so that
  ! its derivative along any other coordinate is zero. An argument
  ! may be a coordinate of the manifold or of its parent, since the
  ! equations on a part are written in the parent's coordinates.
  ! Invalid input: a name already declared, an argument that is not
  ! a coordinate function, one of another manifold, or one repeated.
  !===================================================================!

  function manifold_unknown(this, name, arguments, components) result(u)

    class(continuous_manifold), intent(inout) :: this
    character(len=*)          , intent(in)    :: name
    type(continuous_field)    , intent(in), optional :: arguments(:)
    integer                   , intent(in), optional :: components
    type(continuous_field) :: u

    integer :: k, n
    character(len=MAX_NAME) :: padded
    character(len=8), allocatable :: names(:)
    character(len=250) :: message

    n = 1
    if (present(components)) n = components
    do k = 1, size(this % unknown_name)
       if (trim(this % unknown_name(k)) == trim(name)) then
          error stop 'operation_manifold: an unknown named ' // trim(name) // ' is already declared'
       end if
    end do

    if (present(arguments)) then
       allocate(names(size(arguments)))
       do k = 1, size(arguments)
          if ((arguments(k) % coordinate < 1 .and. .not. arguments(k) % design_coordinate) &
               & .or. allocated(arguments(k) % direction)) then
             write(message,'(a,i0,a,a,a)') 'operation_manifold: argument ', k, ' of the unknown ', trim(name), &
                  & ' is not a coordinate function'
             error stop trim(message)
          end if
          if (.not. (arguments(k) % on % matches(this % identity) .or. arguments(k) % on % matches(this % parent))) then
             write(message,'(a,i0,a,a,a)') 'operation_manifold: argument ', k, ' of the unknown ', trim(name), &
                  & ' is a coordinate of neither the manifold nor its parent'
             error stop trim(message)
          end if
          if (any(names(1:k - 1) == arguments(k) % name)) then
             write(message,'(a,i0,a,a,a)') 'operation_manifold: argument ', k, ' of the unknown ', trim(name), &
                  & ' repeats a coordinate'
             error stop trim(message)
          end if
          names(k) = arguments(k) % name
       end do
    else
       allocate(names(0))
       if (this % with_time)   names = [character(len=8) :: names, 't']
       if (this % with_space)  names = [character(len=8) :: names, 'x', 'y', 'z']
       if (this % with_design) then
          do k = 1, size(this % design_parameter)
             names = [character(len=8) :: names, this % design_parameter(k) % name(1:8)]
          end do
       end if
    end if

    u = unknown_field(this, name, this % num_unknowns + 1, n, names)
    padded = name
    this % unknown_name = [this % unknown_name, (padded, k = 1, n)]
    this % num_unknowns = this % num_unknowns + n

  end function manifold_unknown

  !===================================================================!
  ! The coordinate function named t, x, y or z, or the design
  ! coordinate by its name. A name the manifold has no coordinate for
  ! stops the program.
  !===================================================================!

  function manifold_coordinate(this, name) result(c)

    class(continuous_manifold), intent(in) :: this
    character(len=*)          , intent(in) :: name
    type(continuous_field) :: c

    integer :: axis, at

    if (this % with_design) then
       do at = 1, size(this % design_parameter)
          if (trim(name) == trim(this % design_parameter(at) % name)) then
             c = coordinate_field(this, name, at, of_design=.true.)
             return
          end if
       end do
    end if
    axis = 0
    select case (trim(name))
    case ('t')
       if (this % with_time) axis = 1
       at = 1
    case ('x')
       at = 1
    case ('y')
       at = 2
    case ('z')
       at = 3
    case default
       error stop 'operation_manifold: a coordinate is named t, x, y or z; the name given is ' // trim(name)
    end select
    if (trim(name) /= 't' .and. this % with_space) then
       if (allocated(this % geometry % extent) .and. at > 1) then
          error stop 'operation_manifold: the region of an interval has the one coordinate x; the name &
               &given is ' // trim(name)
       end if
       axis = merge(1, 0, this % with_time) + at
    end if
    if (axis == 0) then
       error stop 'operation_manifold: the manifold has no coordinate named ' // trim(name)
    end if
    c = coordinate_field(this, name, axis)

  end function manifold_coordinate

  !===================================================================!
  ! THE PARTS: the face where one coordinate takes a value at an end
  ! of its interval, a manifold of the other coordinates that records
  ! its parent, the coordinate and the value; the time factor; the
  ! space factor. Each is a manifold with unknowns of its own. On a
  ! region read from a geometry file a face of space is where the
  ! coordinate takes the value, its points those of the discretization
  ! that lie there. Invalid input: a field that is not a coordinate of
  ! this manifold, a design coordinate, or a value at neither end of an
  ! interval.
  !===================================================================!

  function manifold_boundary(this, along, value) result(face)

    class(continuous_manifold), intent(in) :: this
    type(continuous_field)    , intent(in) :: along
    real(dp)                  , intent(in) :: value
    type(continuous_manifold) :: face

    type(interval) :: span
    logical :: of_time
    character(len=250) :: message

    if (along % coordinate < 1 .or. along % design_coordinate .or. .not. along % on % matches(this % identity)) then
       error stop 'operation_manifold: a face is taken along a time or space coordinate of this manifold'
    end if
    of_time = this % with_time .and. along % coordinate == 1
    if (of_time) then
       span = this % span
    else if (allocated(this % geometry % extent)) then
       span = this % geometry % extent
    else
       span % a = value
       span % b = value
    end if
    if (value /= span % a .and. value /= span % b) then
       write(message,'(a,es12.5,a,es12.5,a,es12.5)') 'operation_manifold: the face must be at an end &
            &of the interval; value = ', value, ', interval = [', span % a, ', ', span % b
       error stop trim(message)
    end if
    face = part_of(this, merge(TIME_FACE, SPACE_FACE, of_time))
    face % with_time  = this % with_time .and. .not. of_time
    if (face % with_time) face % span = this % span
    face % with_space = this % with_space .and. of_time
    face % geometry   = this % geometry
    face % face_axis  = along % coordinate
    face % face_value = value

  end function manifold_boundary

  ! a part of a manifold: a manifold of its own that records its
  ! parent and which part of it it is, with the design coordinates of
  ! the parent and no unknown yet; the part taken states the
  ! coordinates of time and space it retains
  function part_of(this, part) result(sub)

    class(continuous_manifold), intent(in) :: this
    integer                   , intent(in) :: part
    type(continuous_manifold) :: sub

    sub % identity    = next_token()
    sub % with_design = this % with_design
    if (allocated(this % design_parameter)) sub % design_parameter = this % design_parameter
    sub % part        = part
    sub % parent      = this % identity
    allocate(sub % unknown_name(0))

  end function part_of

  function manifold_time(this) result(factor)

    class(continuous_manifold), intent(in) :: this
    type(continuous_manifold) :: factor

    if (.not. this % with_time) then
       error stop 'operation_manifold: the time factor requires a manifold with a time coordinate'
    end if
    factor = part_of(this, TIME_FACTOR)
    factor % with_time = .true.
    factor % span      = this % span

  end function manifold_time

  function manifold_space(this) result(factor)

    class(continuous_manifold), intent(in) :: this
    type(continuous_manifold) :: factor

    if (.not. this % with_space) then
       error stop 'operation_manifold: the space factor requires a manifold with a region'
    end if
    factor = part_of(this, SPACE_FACTOR)
    factor % with_space = .true.
    factor % geometry   = this % geometry

  end function manifold_space

  !===================================================================!
  ! A DIRECTION IN THE DESIGN: a vector over the design coordinates,
  ! as a field of the manifold, so that the derivative of a discrete
  ! field along it - the expansion having been stated along it - is
  ! read as along a coordinate: derivative([w]). A vector of another
  ! length than the design coordinates, or zero, is invalid input.
  !===================================================================!

  function manifold_direction(this, along) result(w)

    class(continuous_manifold), intent(in) :: this
    real(dp)                  , intent(in) :: along(:)
    type(continuous_field) :: w

    if (.not. this % with_design) then
       error stop 'operation_manifold: a direction in the design requires a manifold with a design coordinate'
    end if
    if (size(along) /= size(this % design_parameter)) then
       error stop 'operation_manifold: a direction in the design has one entry per design coordinate'
    end if
    if (all(along == 0.0_dp)) then
       error stop 'operation_manifold: a direction in the design is a nonzero vector'
    end if
    w = coordinate_field(this, 'direction', 0, of_design=.true.)
    w % direction = along

  end function manifold_direction

  !===================================================================!
  ! THE DESIGN FACTOR: the point {nu} of the design space, a manifold
  ! of the design coordinate alone, on which the condition fixing the
  ! design is stated and paired with its multiplier, the sensitivity
  ! of the objective to the design.
  !===================================================================!

  function manifold_design(this) result(factor)

    class(continuous_manifold), intent(in) :: this
    type(continuous_manifold) :: factor

    if (.not. this % with_design) then
       error stop 'operation_manifold: the design factor requires a manifold with a design coordinate'
    end if
    factor = part_of(this, DESIGN_FACTOR)

  end function manifold_design

  !===================================================================!
  ! THE DISCRETE IMAGE: the instants of the interval, equally spaced,
  ! the cells of the mesh, and the expansion of an order along the
  ! design. A coordinate present without its points, or points given
  ! for a coordinate absent, stops the program.
  !===================================================================!

  function manifold_discretize(this, time, space, design) result(image)

    class(continuous_manifold), intent(in) :: this
    type(points)   , intent(in), optional :: time
    class(*)       , intent(in), optional :: space
    type(expansion), intent(in), optional :: design
    type(discrete_manifold) :: image

    integer :: k

    if (this % with_time .neqv. present(time)) then
       error stop 'operation_manifold: the time coordinate is discretized by its points, and only it'
    end if
    if (this % with_space .neqv. present(space)) then
       error stop 'operation_manifold: the region is discretized by a mesh or by points, and only it'
    end if
    if (this % with_design .neqv. present(design)) then
       error stop 'operation_manifold: the design coordinate is discretized by an expansion of an order, &
            &and only it'
    end if
    image % identity = this % identity
    if (present(design)) then
       image % expansion = design
       image % expansion % at   = [(this % design_parameter(k) % value, k = 1, size(this % design_parameter))]
       image % expansion % name = [(this % design_parameter(k) % name, k = 1, size(this % design_parameter))]
       if (allocated(design % along)) then
          if (size(design % along, 1) /= size(this % design_parameter)) then
             error stop 'operation_manifold: a direction of the expansion has one entry per design coordinate'
          end if
       else
          allocate(image % expansion % along(size(this % design_parameter), size(this % design_parameter)), source=0.0_dp)
          do k = 1, size(this % design_parameter)
             image % expansion % along(k, k) = 1.0_dp
          end do
       end if
    end if
    if (present(time)) image % instant = positions_on(time, this % span)
    if (present(space)) then
       select type (space)
       type is (mesh)
          image % cells = space
       type is (points)
          if (.not. allocated(this % geometry % extent)) then
             error stop 'operation_manifold: points discretize the region of an interval; a region &
                  &read from a geometry file is discretized by a mesh'
          end if
          image % cells = mesh_from_points(positions_on(space, this % geometry % extent))
       class default
          error stop 'operation_manifold: a region is discretized by a mesh or by points'
       end select
       image % dimension  = image % cells % dimension
       call values_of(image % cells % cell_centre(), image % centre)
       call values_of(image % cells % cell_volume(), image % volume)
    end if

  end function manifold_discretize

  !===================================================================!
  ! The positions of points on an interval: those listed, which run
  ! from one end to the other, or n equally spaced from end to end.
  !===================================================================!

  function positions_on(given, span) result(position)

    type(points)  , intent(in) :: given
    type(interval), intent(in) :: span
    real(dp), allocatable :: position(:)

    integer :: k

    if (allocated(given % position)) then
       if (abs(given % position(1) - span % a) > 1.0e-12_dp * max(1.0_dp, abs(span % a)) .or. &
            & abs(given % position(given % n) - span % b) > 1.0e-12_dp * max(1.0_dp, abs(span % b))) then
          error stop 'operation_manifold: the points listed run from one end of the interval to the other'
       end if
       position = given % position
    else
       position = [(span % a + (span % b - span % a) * real(k - 1, dp) / real(given % n - 1, dp), k = 1, given % n)]
    end if

  end function positions_on

  pure integer function num_instants(this)
    class(discrete_manifold), intent(in) :: this
    num_instants = 1
    if (allocated(this % instant)) num_instants = size(this % instant)
  end function num_instants

  pure integer function num_cells(this)
    class(discrete_manifold), intent(in) :: this
    num_cells = 1
    if (allocated(this % volume)) num_cells = size(this % volume)
  end function num_cells

  pure integer function discrete_num_points(this)
    class(discrete_manifold), intent(in) :: this
    discrete_num_points = this % num_instants() * this % num_cells()
  end function discrete_num_points

  pure integer function discrete_num_coordinates(this)
    class(discrete_manifold), intent(in) :: this
    discrete_num_coordinates = merge(1, 0, allocated(this % instant)) + this % dimension
  end function discrete_num_coordinates

  !===================================================================!
  ! The discretization of the design factor: the one point of measure
  ! one, on which a functional of the solution takes its value and
  ! its derivatives along the design to the order of the expansion.
  !===================================================================!

  function discrete_design_factor(this) result(factor)
    class(discrete_manifold), intent(in) :: this
    class(discrete_support), allocatable :: factor
    type(discrete_manifold) :: point
    point % identity  = next_token()
    point % expansion = this % expansion
    allocate(factor, source=point)
  end function discrete_design_factor

  !===================================================================!
  ! The measure of a point within the face at an instant: the cell's
  ! volume for a point at that instant, one without a region, zero
  ! for a point at another instant.
  !===================================================================!

  pure real(dp) function discrete_measure_on_face(this, p, time)
    class(discrete_manifold), intent(in) :: this
    integer                 , intent(in) :: p
    real(dp)                , intent(in) :: time
    integer :: k, c
    discrete_measure_on_face = 0.0_dp
    if (.not. allocated(this % instant)) return
    k = (p - 1) / this % num_cells() + 1
    c = p - (k - 1) * this % num_cells()
    if (abs(this % instant(k) - time) > 1.0e-12_dp * max(1.0_dp, abs(time))) return
    discrete_measure_on_face = 1.0_dp
    if (allocated(this % volume)) discrete_measure_on_face = this % volume(c)
  end function discrete_measure_on_face

  ! the point of instant k and cell c: instant outer, cell inner
  pure integer function point_of(this, k, c)
    class(discrete_manifold), intent(in) :: this
    integer                 , intent(in) :: k, c
    point_of = (k - 1) * this % num_cells() + c
  end function point_of

  ! the step ending at instant k; the first instant has none
  pure real(dp) function step(this, k)
    class(discrete_manifold), intent(in) :: this
    integer                 , intent(in) :: k
    step = 0.0_dp
    if (k > 1) step = this % instant(k) - this % instant(k - 1)
  end function step

  pure function discrete_position(this, p) result(x)

    class(discrete_manifold), intent(in) :: this
    integer                 , intent(in) :: p
    real(dp), allocatable :: x(:)

    integer :: k, c, d

    k = (p - 1) / this % num_cells() + 1
    c = p - (k - 1) * this % num_cells()
    d = this % dimension
    allocate(x(this % num_coordinates()))
    if (allocated(this % instant)) x(1) = this % instant(k)
    if (allocated(this % volume)) then
       x(this % num_coordinates() - d + 1:) = this % centre(d * (c - 1) + 1:d * c)
    end if

  end function discrete_position

  !===================================================================!
  ! The measure at a point: the trapezoidal share of the interval at
  ! the instant times the cell's volume; one at a point manifold.
  !===================================================================!

  pure real(dp) function discrete_measure(this, p)

    class(discrete_manifold), intent(in) :: this
    integer                 , intent(in) :: p

    integer :: k, c, n

    k = (p - 1) / this % num_cells() + 1
    c = p - (k - 1) * this % num_cells()
    discrete_measure = 1.0_dp
    if (allocated(this % instant)) then
       n = size(this % instant)
       if (k == 1) then
          discrete_measure = 0.5_dp * this % step(2)
       else if (k == n) then
          discrete_measure = 0.5_dp * this % step(n)
       else
          discrete_measure = 0.5_dp * (this % step(k) + this % step(k + 1))
       end if
    end if
    if (allocated(this % volume)) discrete_measure = discrete_measure * this % volume(c)

  end function discrete_measure

end module operation_manifold
