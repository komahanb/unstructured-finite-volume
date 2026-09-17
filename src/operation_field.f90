!=====================================================================!
! FIELDS ON A MANIFOLD, CONTINUOUS AND DISCRETE.
!
! A continuous field is a function on a manifold, stored as an
! expression graph per component: a coordinate function, an unknown
! function of the manifold, or a formula in these and in constants.
! A discrete field is the same function on the points of a
! discretization of the manifold: one value per point per component.
! The one map between them, discretize, samples the continuous field
! at the points; the residual, defined elsewhere, solves for an unknown
! there instead.
!
! The manifolds themselves are defined above this module. A field
! sees its manifold through the two abstract supports below: the
! continuous one by its identity, the discrete one by its points,
! their positions and their measure. A manifold extends the first,
! its discretization the second.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_field

  use util_precision      , only : dp
  use token_identity      , only : token
  use operation_expression, only : expression, unknown, constant, coordinate, derivative, derivative_along
  use operation_expression, only : FIRST_COORDINATE
  use operation_expression, only : operator(+), operator(-), operator(*), operator(/), operator(**)
  use operation_expression, only : sin, cos, exp, log, sqrt

  implicit none

  private
  public :: continuous_support, discrete_support
  public :: continuous_field, discrete_field
  public :: unknown_field, coordinate_field, integral
  public :: operator(+), operator(-), operator(*), operator(/), operator(**)
  public :: sin, cos, exp, log, sqrt

  !===================================================================!
  ! THE SUPPORTS. A continuous support is known to a field by its
  ! identity alone; a discrete support by its points.
  !===================================================================!

  type, abstract :: continuous_support

     type(token) :: identity

   contains

     procedure :: same_support

  end type continuous_support

  type, abstract :: discrete_support

     type(token) :: identity        ! of the continuous support discretized

   contains

     procedure(count_interface)   , deferred :: num_points
     procedure(count_interface)   , deferred :: num_coordinates
     procedure(position_interface), deferred :: position
     procedure(measure_interface) , deferred :: measure

  end type discrete_support

  abstract interface

     pure integer function count_interface(this)
       import :: discrete_support
       class(discrete_support), intent(in) :: this
     end function count_interface

     ! the coordinates of point p
     pure function position_interface(this, p) result(x)
       import :: discrete_support, dp
       class(discrete_support), intent(in) :: this
       integer                , intent(in) :: p
       real(dp), allocatable :: x(:)
     end function position_interface

     ! the share of the manifold's measure at point p
     pure real(dp) function measure_interface(this, p)
       import :: discrete_support, dp
       class(discrete_support), intent(in) :: this
       integer                , intent(in) :: p
     end function measure_interface

  end interface

  !===================================================================!
  ! THE CONTINUOUS FIELD: one expression per component, on one
  ! support. An unknown has the field number of each component
  ! on its manifold; a coordinate function has the coordinate it
  ! is; an integral has the identity of the factor it runs over.
  !===================================================================!

  type :: continuous_field

     type(token) :: on
     character(len=:), allocatable :: name
     type(expression), allocatable :: component(:)
     integer, allocatable :: index(:)
     integer :: coordinate = 0
     type(token) :: integrated_over

   contains

     procedure :: derivative     => field_derivative
     procedure :: discretize     => field_discretize
     procedure :: num_components => field_num_components
     procedure :: graph          => field_graph
     procedure :: is_unknown
     procedure :: is_coordinate
     procedure :: is_integral

  end type continuous_field

  interface continuous_field
     module procedure create_tuple
     module procedure create_constants
  end interface continuous_field

  !===================================================================!
  ! THE DISCRETE FIELD: its support, and one value per point and
  ! component.
  !===================================================================!

  type :: discrete_field

     class(discrete_support), allocatable :: on
     character(len=32), allocatable :: name(:)
     real(dp), allocatable :: value(:,:)

   contains

     procedure :: values
     procedure :: fields
     procedure :: norm
     procedure :: num_points     => discrete_num_points
     procedure :: num_components => discrete_num_components

  end type discrete_field

  interface discrete_field
     module procedure create_discrete
  end interface discrete_field

  interface operator(+)
     module procedure plus, real_plus, plus_real
  end interface operator(+)

  interface operator(-)
     module procedure minus, real_minus, minus_real, negated, discrete_minus
  end interface operator(-)

  interface operator(*)
     module procedure times, real_times, times_real
  end interface operator(*)

  interface operator(/)
     module procedure over, real_over, over_real
  end interface operator(/)

  interface operator(**)
     module procedure to_integer, to_real
  end interface operator(**)

  interface sin
     module procedure sine_of
  end interface sin

  interface cos
     module procedure cosine_of
  end interface cos

  interface exp
     module procedure exponential_of
  end interface exp

  interface log
     module procedure logarithm_of
  end interface log

  interface sqrt
     module procedure root_of
  end interface sqrt

contains

  pure logical function same_support(this, other)

    class(continuous_support), intent(in) :: this
    type(token)              , intent(in) :: other

    same_support = this % identity % matches(other)

  end function same_support

  !===================================================================!
  ! THE CONSTRUCTORS. An unknown of a manifold: one expression leaf
  ! per component, numbered from first. A coordinate function. A tuple
  ! of fields on one support, each of one component. Constants, one
  ! per component.
  !===================================================================!

  function unknown_field(on, name, first, components) result(this)

    class(continuous_support), intent(in) :: on
    character(len=*)         , intent(in) :: name
    integer                  , intent(in) :: first, components
    type(continuous_field) :: this

    integer :: k

    this % on   = on % identity
    this % name = name
    allocate(this % component(components), this % index(components))
    do k = 1, components
       this % index(k)     = first + k - 1
       this % component(k) = unknown(first + k - 1)
    end do

  end function unknown_field

  function coordinate_field(on, name, c) result(this)

    class(continuous_support), intent(in) :: on
    character(len=*)         , intent(in) :: name
    integer                  , intent(in) :: c
    type(continuous_field) :: this

    this % on         = on % identity
    this % name       = name
    this % coordinate = c
    this % component  = [coordinate(c)]
    this % index      = [0]

  end function coordinate_field

  function create_tuple(on, fields) result(this)

    class(continuous_support), intent(in) :: on
    type(continuous_field)   , intent(in) :: fields(:)
    type(continuous_field) :: this

    integer :: k
    character(len=250) :: message

    this % on   = on % identity
    this % name = ''
    allocate(this % component(size(fields)), this % index(size(fields)))
    do k = 1, size(fields)
       if (.not. on % same_support(fields(k) % on)) then
          write(message,'(a,i0,a)') 'operation_field: a tuple requires every field on the support &
               &given; field ', k, ' is on another'
          error stop trim(message)
       end if
       if (size(fields(k) % component) /= 1) then
          write(message,'(a,i0,a,i0)') 'operation_field: a tuple requires fields of one component; &
               &field ', k, ' has ', size(fields(k) % component)
          error stop trim(message)
       end if
       this % component(k) = fields(k) % component(1)
       this % index(k)     = fields(k) % index(1)
    end do

  end function create_tuple

  function create_constants(on, values) result(this)

    class(continuous_support), intent(in) :: on
    real(dp)                 , intent(in) :: values(:)
    type(continuous_field) :: this

    integer :: k

    this % on   = on % identity
    this % name = ''
    allocate(this % component(size(values)), this % index(size(values)))
    do k = 1, size(values)
       this % component(k) = constant(values(k))
       this % index(k)     = 0
    end do

  end function create_constants

  !===================================================================!
  ! The integral of a field over a factor of its manifold: the field,
  ! marked with the factor's identity. The discrete residual sums it
  ! over the points of that factor.
  !===================================================================!

  function integral(f, over) result(this)

    type(continuous_field)   , intent(in) :: f
    class(continuous_support), intent(in) :: over
    type(continuous_field) :: this

    this = f
    this % integrated_over = over % identity

  end function integral

  pure integer function field_num_components(this)
    class(continuous_field), intent(in) :: this
    field_num_components = size(this % component)
  end function field_num_components

  function field_graph(this, k) result(g)
    class(continuous_field), intent(in) :: this
    integer                , intent(in) :: k
    type(expression) :: g
    g = this % component(k)
  end function field_graph

  pure logical function is_unknown(this)
    class(continuous_field), intent(in) :: this
    is_unknown = all(this % index > 0)
  end function is_unknown

  pure logical function is_coordinate(this)
    class(continuous_field), intent(in) :: this
    is_coordinate = this % coordinate > 0
  end function is_coordinate

  pure logical function is_integral(this)
    class(continuous_field), intent(in) :: this
    is_integral = this % integrated_over % declared()
  end function is_integral

  !===================================================================!
  ! THE DERIVATIVE along a multi-index of coordinate functions: each
  ! entry one order along that coordinate. Invalid input: a field that
  ! is not a scalar unknown, an entry that is not a coordinate of the
  ! same manifold, or two distinct coordinates, since the jet stores
  ! no mixed component.
  !===================================================================!

  function field_derivative(this, along) result(d)

    class(continuous_field), intent(in) :: this
    type(continuous_field) , intent(in) :: along(:)
    type(continuous_field) :: d

    integer :: k, c
    character(len=250) :: message

    if (size(this % component) /= 1 .or. this % index(1) < 1) then
       error stop 'operation_field: derivative requires a scalar unknown of the manifold'
    end if
    if (size(along) < 1) then
       error stop 'operation_field: derivative requires one coordinate at least'
    end if
    c = along(1) % coordinate
    do k = 1, size(along)
       if (along(k) % coordinate < 1) then
          write(message,'(a,i0,a)') 'operation_field: entry ', k, ' of the multi-index is not a coordinate'
          error stop trim(message)
       end if
       if (.not. along(k) % on % matches(this % on)) then
          write(message,'(a,i0,a)') 'operation_field: entry ', k, ' of the multi-index is a coordinate &
               &of another manifold'
          error stop trim(message)
       end if
       if (along(k) % coordinate /= c) then
          error stop 'operation_field: a mixed derivative is not a component of the jet; the multi-index &
               &must repeat one coordinate'
       end if
    end do

    d % on    = this % on
    d % name  = this % name
    d % index = [0]
    if (c == FIRST_COORDINATE) then
       d % component = [derivative(this % component(1), size(along))]
    else
       d % component = [derivative_along(this % component(1), c, size(along))]
    end if

  end function field_derivative

  !===================================================================!
  ! THE DISCRETE IMAGE: the field sampled at the points. Invalid
  ! input: a discretization of another manifold, or a field that
  ! reads an unknown, which has no value to sample.
  !===================================================================!

  function field_discretize(this, on) result(sampled)

    class(continuous_field), intent(in) :: this
    class(discrete_support), intent(in) :: on
    type(discrete_field) :: sampled

    integer :: p, k

    if (.not. on % identity % matches(this % on)) then
       error stop 'operation_field: discretize requires the discretization of the field''s own manifold'
    end if
    allocate(sampled % on, source=on)
    allocate(sampled % name(size(this % component)))
    sampled % name = this % name
    allocate(sampled % value(on % num_points(), size(this % component)))
    do p = 1, on % num_points()
       do k = 1, size(this % component)
          sampled % value(p, k) = this % component(k) % value_at(on % position(p))
       end do
    end do

  end function field_discretize

  !===================================================================!
  ! THE DISCRETE FIELD from values: value(p, k) at point p, component
  ! k, named per component.
  !===================================================================!

  function create_discrete(on, value, name) result(this)

    class(discrete_support), intent(in) :: on
    real(dp)               , intent(in) :: value(:,:)
    character(len=*)       , intent(in), optional :: name(:)
    type(discrete_field) :: this

    character(len=250) :: message

    if (size(value, 1) /= on % num_points()) then
       write(message,'(a,i0,a,i0)') 'operation_field: one row of values is required per point; &
            &size(value, 1) = ', size(value, 1), ', num_points = ', on % num_points()
       error stop trim(message)
    end if
    allocate(this % on, source=on)
    this % value = value
    allocate(this % name(size(value, 2)))
    this % name = ''
    if (present(name)) this % name = name

  end function create_discrete

  pure integer function discrete_num_points(this)
    class(discrete_field), intent(in) :: this
    discrete_num_points = size(this % value, 1)
  end function discrete_num_points

  pure integer function discrete_num_components(this)
    class(discrete_field), intent(in) :: this
    discrete_num_components = size(this % value, 2)
  end function discrete_num_components

  !===================================================================!
  ! The values in the order of the points, each point's components
  ! together. An array of another size than points times components
  ! is invalid input.
  !===================================================================!

  subroutine values(this, v)

    class(discrete_field), intent(in)  :: this
    real(dp)             , intent(out) :: v(:)

    integer :: p, k, nc
    character(len=250) :: message

    nc = size(this % value, 2)
    if (size(v) /= size(this % value, 1) * nc) then
       write(message,'(a,i0,a,i0)') 'operation_field: the values array has one entry per point and &
            &component; size = ', size(v), ', required = ', size(this % value, 1) * nc
       error stop trim(message)
    end if
    do p = 1, size(this % value, 1)
       do k = 1, nc
          v((p - 1) * nc + k) = this % value(p, k)
       end do
    end do

  end subroutine values

  !===================================================================!
  ! The components named, in the order named. A name not present
  ! stops the program.
  !===================================================================!

  function fields(this, names) result(chosen)

    class(discrete_field), intent(in) :: this
    character(len=*)     , intent(in) :: names(:)
    type(discrete_field) :: chosen

    integer :: k, j, at

    allocate(chosen % on, source=this % on)
    allocate(chosen % name(size(names)), chosen % value(size(this % value, 1), size(names)))
    do k = 1, size(names)
       at = 0
       do j = 1, size(this % name)
          if (trim(this % name(j)) == trim(names(k))) at = j
       end do
       if (at == 0) then
          error stop 'operation_field: fields names a component the field does not have: ' // trim(names(k))
       end if
       chosen % name(k)     = this % name(at)
       chosen % value(:, k) = this % value(:, at)
    end do

  end function fields

  !===================================================================!
  ! The root mean square over the manifold: the sum over the points
  ! of the measure at the point times the sum of the squares of the
  ! components, divided by the total measure, under the root.
  !===================================================================!

  pure real(dp) function norm(this)

    class(discrete_field), intent(in) :: this

    real(dp) :: weighted, total
    integer :: p

    weighted = 0.0_dp
    total    = 0.0_dp
    do p = 1, size(this % value, 1)
       weighted = weighted + this % on % measure(p) * sum(this % value(p, :)**2)
       total    = total    + this % on % measure(p)
    end do
    norm = sqrt(weighted / total)

  end function norm

  !===================================================================!
  ! The difference of two discrete fields on one support with equal
  ! component counts; the names are the left field's.
  !===================================================================!

  function discrete_minus(a, b) result(d)

    type(discrete_field), intent(in) :: a, b
    type(discrete_field) :: d

    if (.not. a % on % identity % matches(b % on % identity)) then
       error stop 'operation_field: a difference requires discrete fields of one manifold'
    end if
    if (size(a % value, 2) /= size(b % value, 2) .or. size(a % value, 1) /= size(b % value, 1)) then
       error stop 'operation_field: a difference requires discrete fields of equal extent'
    end if
    allocate(d % on, source=a % on)
    d % name  = a % name
    d % value = a % value - b % value

  end function discrete_minus

  !===================================================================!
  ! THE OPERATORS, on scalar fields of one manifold. A tuple as an
  ! operand, or operands of two manifolds, stops the program.
  !===================================================================!

  subroutine require_scalar(a, symbol)

    type(continuous_field), intent(in) :: a
    character(len=*)      , intent(in) :: symbol

    if (size(a % component) /= 1) then
       error stop 'operation_field: ' // symbol // ' requires fields of one component'
    end if

  end subroutine require_scalar

  subroutine require_one_manifold(a, b, symbol)

    type(continuous_field), intent(in) :: a, b
    character(len=*)      , intent(in) :: symbol

    if (.not. a % on % matches(b % on)) then
       error stop 'operation_field: ' // symbol // ' requires fields of one manifold'
    end if

  end subroutine require_one_manifold

  function composed(a, g) result(this)
    type(continuous_field), intent(in) :: a
    type(expression)      , intent(in) :: g
    type(continuous_field) :: this
    this % on        = a % on
    this % name      = ''
    this % component = [g]
    this % index     = [0]
  end function composed

  function plus(a, b) result(this)
    type(continuous_field), intent(in) :: a, b
    type(continuous_field) :: this
    call require_scalar(a, '+'); call require_scalar(b, '+'); call require_one_manifold(a, b, '+')
    this = composed(a, a % component(1) + b % component(1))
  end function plus

  function real_plus(x, b) result(this)
    real(dp)              , intent(in) :: x
    type(continuous_field), intent(in) :: b
    type(continuous_field) :: this
    call require_scalar(b, '+')
    this = composed(b, x + b % component(1))
  end function real_plus

  function plus_real(a, x) result(this)
    type(continuous_field), intent(in) :: a
    real(dp)              , intent(in) :: x
    type(continuous_field) :: this
    call require_scalar(a, '+')
    this = composed(a, a % component(1) + x)
  end function plus_real

  function minus(a, b) result(this)
    type(continuous_field), intent(in) :: a, b
    type(continuous_field) :: this
    call require_scalar(a, '-'); call require_scalar(b, '-'); call require_one_manifold(a, b, '-')
    this = composed(a, a % component(1) - b % component(1))
  end function minus

  function real_minus(x, b) result(this)
    real(dp)              , intent(in) :: x
    type(continuous_field), intent(in) :: b
    type(continuous_field) :: this
    call require_scalar(b, '-')
    this = composed(b, x - b % component(1))
  end function real_minus

  function minus_real(a, x) result(this)
    type(continuous_field), intent(in) :: a
    real(dp)              , intent(in) :: x
    type(continuous_field) :: this
    call require_scalar(a, '-')
    this = composed(a, a % component(1) - x)
  end function minus_real

  function negated(a) result(this)
    type(continuous_field), intent(in) :: a
    type(continuous_field) :: this
    call require_scalar(a, '-')
    this = composed(a, -a % component(1))
  end function negated

  function times(a, b) result(this)
    type(continuous_field), intent(in) :: a, b
    type(continuous_field) :: this
    call require_scalar(a, '*'); call require_scalar(b, '*'); call require_one_manifold(a, b, '*')
    this = composed(a, a % component(1) * b % component(1))
  end function times

  function real_times(x, b) result(this)
    real(dp)              , intent(in) :: x
    type(continuous_field), intent(in) :: b
    type(continuous_field) :: this
    call require_scalar(b, '*')
    this = composed(b, x * b % component(1))
  end function real_times

  function times_real(a, x) result(this)
    type(continuous_field), intent(in) :: a
    real(dp)              , intent(in) :: x
    type(continuous_field) :: this
    call require_scalar(a, '*')
    this = composed(a, a % component(1) * x)
  end function times_real

  function over(a, b) result(this)
    type(continuous_field), intent(in) :: a, b
    type(continuous_field) :: this
    call require_scalar(a, '/'); call require_scalar(b, '/'); call require_one_manifold(a, b, '/')
    this = composed(a, a % component(1) / b % component(1))
  end function over

  function real_over(x, b) result(this)
    real(dp)              , intent(in) :: x
    type(continuous_field), intent(in) :: b
    type(continuous_field) :: this
    call require_scalar(b, '/')
    this = composed(b, x / b % component(1))
  end function real_over

  function over_real(a, x) result(this)
    type(continuous_field), intent(in) :: a
    real(dp)              , intent(in) :: x
    type(continuous_field) :: this
    call require_scalar(a, '/')
    this = composed(a, a % component(1) / x)
  end function over_real

  function to_integer(a, n) result(this)
    type(continuous_field), intent(in) :: a
    integer               , intent(in) :: n
    type(continuous_field) :: this
    call require_scalar(a, '**')
    this = composed(a, a % component(1) ** n)
  end function to_integer

  function to_real(a, p) result(this)
    type(continuous_field), intent(in) :: a
    real(dp)              , intent(in) :: p
    type(continuous_field) :: this
    call require_scalar(a, '**')
    this = composed(a, a % component(1) ** p)
  end function to_real

  function sine_of(a) result(this)
    type(continuous_field), intent(in) :: a
    type(continuous_field) :: this
    call require_scalar(a, 'sin')
    this = composed(a, sin(a % component(1)))
  end function sine_of

  function cosine_of(a) result(this)
    type(continuous_field), intent(in) :: a
    type(continuous_field) :: this
    call require_scalar(a, 'cos')
    this = composed(a, cos(a % component(1)))
  end function cosine_of

  function exponential_of(a) result(this)
    type(continuous_field), intent(in) :: a
    type(continuous_field) :: this
    call require_scalar(a, 'exp')
    this = composed(a, exp(a % component(1)))
  end function exponential_of

  function logarithm_of(a) result(this)
    type(continuous_field), intent(in) :: a
    type(continuous_field) :: this
    call require_scalar(a, 'log')
    this = composed(a, log(a % component(1)))
  end function logarithm_of

  function root_of(a) result(this)
    type(continuous_field), intent(in) :: a
    type(continuous_field) :: this
    call require_scalar(a, 'sqrt')
    this = composed(a, sqrt(a % component(1)))
  end function root_of

end module operation_field
