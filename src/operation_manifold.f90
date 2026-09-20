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
!
!      coordinate index   1 the time, then the space axes in order
!      point index        instant outer, cell inner
!      position(p)        the instant and the cell centre
!      measure(p)         the instant's share of the interval times the cell's volume
!
! THE JET ALONG THE DESIGN. A jet is the derivatives of a quantity
! along a coordinate at one point, to an order: the coefficients of
! its Taylor expansion there, the expansion not summed. The design
! coordinate nu is a point of the design space, fixed at a value
! nu_0, and its discretization is the expansion of an order at that
! value, so that the discrete manifold is the instants (and cells)
! times the orders of the expansion:
!
!             nu_0 ------>  nu  (the design, one point, expanded to order 3)
!
!   t_n  |  [u, u_t, u_tt]  [u, u_t, u_tt]'  [u, u_t, u_tt]''  [u, u_t, u_tt]'''
!    :   |        :
!   t_2  |  [u, u_t, u_tt]  [u, u_t, u_tt]'  [u, u_t, u_tt]''  [u, u_t, u_tt]'''
!   t_1  |  [u, u_t, u_tt]  [u, u_t, u_tt]'  [u, u_t, u_tt]''  [u, u_t, u_tt]'''
!        |
!   time v      order 0          order 1          order 2           order 3
!              (the state)     (d/dnu)          (d^2/dnu^2)       (d^3/dnu^3)
!
! Every point stores the jet of the unknowns along time, which the
! families connect; column 0 is the solution, column 1 its
! derivative along nu at every instant, column 2 the second
! derivative, and so on. The derivative of a discrete field along nu
! is the columns from the next on. The multipliers have the same
! columns: the equations' at every instant, the reactions', the
! gauges', and the design condition's, whose column m is minus the
! (m + 1)-th derivative of the objective.
!
! The columns follow from the equations holding at every nu: with
! F(u(nu), nu) = 0 and A the jacobian of F in u, A u' + dF/dnu = 0
! gives column 1 by one linear solve with A, A u'' + (terms in u',
! u) = 0 gives column 2 by one more solve with the same A, and so
! on; the jet arithmetic of util_derivative_terms evaluates F on the
! jets of u and nu and returns the jet of F, whose coefficient at
! order m is the right side of the m-th solve. Nothing is
! differentiated by hand and nothing is approximated by differences.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_manifold

  use util_precision      , only : dp
  use token_identity      , only : token, next_token
  use operation_field     , only : continuous_support, discrete_support, continuous_field
  use operation_field     , only : unknown_field, coordinate_field
  use operation_field     , only : WHOLE, TIME_FACE, TIME_FACTOR, SPACE_FACTOR, DESIGN_FACTOR
  use view_mesh           , only : mesh, values_of
  use view_mesh_builder   , only : mesh_from_gmsh

  implicit none

  private
  public :: continuous_manifold, discrete_manifold
  public :: interval, region, instants, mesh, parameter, expansion

  integer, parameter :: MAX_NAME = 32

  type :: interval
     real(dp) :: a = 0.0_dp, b = 0.0_dp
  end type interval

  interface interval
     module procedure create_interval
  end interface interval

  ! the region of space, as the gmsh geometry file describes it
  type :: region
     character(len=:), allocatable :: geometry
  end type region

  interface region
     module procedure create_region
  end interface region

  type :: instants
     integer :: n = 0
  end type instants

  interface instants
     module procedure create_instants
  end interface instants

  interface mesh
     module procedure mesh_from_file
  end interface mesh

  ! the design: a coordinate fixed at a value, a point of the design
  ! space; its discretization is the Taylor expansion of an order at
  ! that value, so that the discrete solution stores its derivatives
  ! along the design to the order
  type :: parameter
     character(len=MAX_NAME) :: name = ''
     real(dp) :: value = 0.0_dp
  end type parameter

  interface parameter
     module procedure create_parameter
  end interface parameter

  type :: expansion
     integer :: order = 0
  end type expansion

  interface expansion
     module procedure create_expansion
  end interface expansion

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
     procedure :: discretize      => manifold_discretize
     procedure :: num_coordinates => manifold_num_coordinates
     procedure :: is_point
     procedure :: is_face
     procedure :: is_time_factor
     procedure :: is_space_factor
     procedure :: is_design_factor

  end type continuous_manifold

  interface continuous_manifold
     module procedure create_manifold
     module procedure create_manifold_of_designs
  end interface continuous_manifold

  type, extends(discrete_support) :: discrete_manifold

     logical :: with_time  = .false.
     logical :: with_space = .false.
     logical :: with_design = .false.
     real(dp), allocatable :: instant(:)
     type(mesh) :: cells
     integer :: dimension = 0
     real(dp), allocatable :: centre(:)
     real(dp), allocatable :: volume(:)
     type(parameter), allocatable :: design_parameter(:)
     integer :: order = 0

   contains

     procedure :: num_points      => discrete_num_points
     procedure :: num_coordinates => discrete_num_coordinates
     procedure :: position        => discrete_position
     procedure :: measure         => discrete_measure
     procedure :: design_factor   => discrete_design_factor
     procedure :: design_order    => discrete_design_order
     procedure :: num_designs     => discrete_num_designs
     procedure :: design_values   => discrete_design_values
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

  function create_instants(n) result(this)

    integer, intent(in) :: n
    type(instants) :: this

    character(len=250) :: message

    if (n < 2) then
       write(message,'(a,i0)') 'operation_manifold: an interval is discretized by two instants at &
            &least; n = ', n
       error stop trim(message)
    end if
    this % n = n

  end function create_instants

  function mesh_from_file(filename) result(cells)

    character(len=*), intent(in) :: filename
    type(mesh) :: cells

    cells = mesh_from_gmsh(filename)

  end function mesh_from_file

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

  function create_expansion(order) result(this)

    integer, intent(in) :: order
    type(expansion) :: this

    character(len=250) :: message

    if (order < 0) then
       write(message,'(a,i0)') 'operation_manifold: the expansion along the design is of order zero at &
            &least; order = ', order
       error stop trim(message)
    end if
    this % order = order

  end function create_expansion

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

  pure logical function is_face(this)
    class(continuous_manifold), intent(in) :: this
    is_face = this % part == TIME_FACE
  end function is_face

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
    manifold_num_coordinates = merge(1, 0, this % with_time) + merge(3, 0, this % with_space)
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
          if (arguments(k) % coordinate < 1 .and. .not. arguments(k) % design_coordinate) then
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
    if (trim(name) /= 't') then
       if (this % with_space) axis = merge(1, 0, this % with_time) + at
    end if
    if (axis == 0) then
       error stop 'operation_manifold: the manifold has no coordinate named ' // trim(name)
    end if
    c = coordinate_field(this, name, axis)

  end function manifold_coordinate

  !===================================================================!
  ! THE PARTS: the face at one instant, a manifold of the space alone
  ! that records its parent and the instant; the time factor; the
  ! space factor. Each is a manifold with unknowns of its own.
  !===================================================================!

  function manifold_boundary(this, time) result(face)

    class(continuous_manifold), intent(in) :: this
    real(dp)                  , intent(in) :: time
    type(continuous_manifold) :: face

    character(len=250) :: message

    if (.not. this % with_time) then
       error stop 'operation_manifold: a face at an instant requires a manifold with a time coordinate'
    end if
    if (time /= this % span % a .and. time /= this % span % b) then
       write(message,'(a,es12.5,a,es12.5,a,es12.5)') 'operation_manifold: the face must be at an end &
            &of the interval; time = ', time, ', interval = [', this % span % a, ', ', this % span % b
       error stop trim(message)
    end if
    face % identity   = next_token()
    face % with_space = this % with_space
    face % geometry   = this % geometry
    face % with_design = this % with_design
    if (allocated(this % design_parameter)) face % design_parameter = this % design_parameter
    face % part       = TIME_FACE
    face % face_time  = time
    face % parent     = this % identity
    allocate(face % unknown_name(0))

  end function manifold_boundary

  function manifold_time(this) result(factor)

    class(continuous_manifold), intent(in) :: this
    type(continuous_manifold) :: factor

    if (.not. this % with_time) then
       error stop 'operation_manifold: the time factor requires a manifold with a time coordinate'
    end if
    factor % identity  = next_token()
    factor % with_time = .true.
    factor % span      = this % span
    factor % with_design = this % with_design
    if (allocated(this % design_parameter)) factor % design_parameter = this % design_parameter
    factor % part      = TIME_FACTOR
    factor % parent    = this % identity
    allocate(factor % unknown_name(0))

  end function manifold_time

  function manifold_space(this) result(factor)

    class(continuous_manifold), intent(in) :: this
    type(continuous_manifold) :: factor

    if (.not. this % with_space) then
       error stop 'operation_manifold: the space factor requires a manifold with a region'
    end if
    factor % identity   = next_token()
    factor % with_space = .true.
    factor % geometry   = this % geometry
    factor % with_design = this % with_design
    if (allocated(this % design_parameter)) factor % design_parameter = this % design_parameter
    factor % part       = SPACE_FACTOR
    factor % parent     = this % identity
    allocate(factor % unknown_name(0))

  end function manifold_space

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
    factor % identity    = next_token()
    factor % with_design = .true.
    if (allocated(this % design_parameter)) factor % design_parameter = this % design_parameter
    factor % part        = DESIGN_FACTOR
    factor % parent      = this % identity
    allocate(factor % unknown_name(0))

  end function manifold_design

  !===================================================================!
  ! THE DISCRETE IMAGE: the instants of the interval, equally spaced,
  ! the cells of the mesh, and the expansion of an order along the
  ! design. A coordinate present without its points, or points given
  ! for a coordinate absent, stops the program.
  !===================================================================!

  function manifold_discretize(this, time, space, design) result(image)

    class(continuous_manifold), intent(in) :: this
    type(instants) , intent(in), optional :: time
    type(mesh)     , intent(in), optional :: space
    type(expansion), intent(in), optional :: design
    type(discrete_manifold) :: image

    integer :: k

    if (this % with_time .neqv. present(time)) then
       error stop 'operation_manifold: the time coordinate is discretized by its instants, and only it'
    end if
    if (this % with_space .neqv. present(space)) then
       error stop 'operation_manifold: the region is discretized by a mesh, and only it'
    end if
    if (this % with_design .neqv. present(design)) then
       error stop 'operation_manifold: the design coordinate is discretized by an expansion of an order, &
            &and only it'
    end if
    image % identity = this % identity
    if (present(design)) then
       image % with_design = .true.
       image % design_parameter = this % design_parameter
       image % order       = design % order
    end if
    if (present(time)) then
       image % with_time = .true.
       image % instant = [(this % span % a + (this % span % b - this % span % a) * real(k - 1, dp) &
            & / real(time % n - 1, dp), k = 1, time % n)]
    end if
    if (present(space)) then
       image % with_space = .true.
       image % cells      = space
       image % dimension  = space % dimension
       call values_of(space % cell_centre(), image % centre)
       call values_of(space % cell_volume(), image % volume)
    end if

  end function manifold_discretize

  pure integer function num_instants(this)
    class(discrete_manifold), intent(in) :: this
    num_instants = 1
    if (this % with_time) num_instants = size(this % instant)
  end function num_instants

  pure integer function num_cells(this)
    class(discrete_manifold), intent(in) :: this
    num_cells = 1
    if (this % with_space) num_cells = size(this % volume)
  end function num_cells

  pure integer function discrete_num_points(this)
    class(discrete_manifold), intent(in) :: this
    discrete_num_points = this % num_instants() * this % num_cells()
  end function discrete_num_points

  pure integer function discrete_num_coordinates(this)
    class(discrete_manifold), intent(in) :: this
    discrete_num_coordinates = merge(1, 0, this % with_time) + this % dimension
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
    point % identity    = next_token()
    point % with_design = this % with_design
    if (allocated(this % design_parameter)) point % design_parameter = this % design_parameter
    point % order       = this % order
    allocate(factor, source=point)
  end function discrete_design_factor

  pure integer function discrete_design_order(this)
    class(discrete_manifold), intent(in) :: this
    discrete_design_order = this % order
  end function discrete_design_order

  pure integer function discrete_num_designs(this)
    class(discrete_manifold), intent(in) :: this
    discrete_num_designs = 0
    if (this % with_design) discrete_num_designs = size(this % design_parameter)
  end function discrete_num_designs

  ! the design vector: the value of every design coordinate, or the
  ! one value zero without a design, which no expression reads
  pure function discrete_design_values(this) result(v)
    class(discrete_manifold), intent(in) :: this
    real(dp), allocatable :: v(:)
    integer :: i
    if (.not. this % with_design) then
       allocate(v(1), source=0.0_dp)
       return
    end if
    allocate(v(size(this % design_parameter)))
    do i = 1, size(v)
       v(i) = this % design_parameter(i) % value
    end do
  end function discrete_design_values

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
    if (.not. this % with_time) return
    k = (p - 1) / this % num_cells() + 1
    c = p - (k - 1) * this % num_cells()
    if (abs(this % instant(k) - time) > 1.0e-12_dp * max(1.0_dp, abs(time))) return
    discrete_measure_on_face = 1.0_dp
    if (this % with_space) discrete_measure_on_face = this % volume(c)
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
    if (this % with_time) x(1) = this % instant(k)
    if (this % with_space) then
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
    if (this % with_time) then
       n = size(this % instant)
       if (k == 1) then
          discrete_measure = 0.5_dp * this % step(2)
       else if (k == n) then
          discrete_measure = 0.5_dp * this % step(n)
       else
          discrete_measure = 0.5_dp * (this % step(k) + this % step(k + 1))
       end if
    end if
    if (this % with_space) discrete_measure = discrete_measure * this % volume(c)

  end function discrete_measure

end module operation_manifold
