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
  use util_derivative_terms, only : derivative_terms, coefficient
  use operation_expression, only : expression, unknown, constant, coordinate, derivative, derivative_along, design
  use operation_expression, only : FIRST_COORDINATE
  use operation_expression, only : operator(+), operator(-), operator(*), operator(/), operator(**)
  use operation_expression, only : sin, cos, exp, log, sqrt

  implicit none

  private
  public :: continuous_support, discrete_support
  public :: continuous_field, discrete_field
  public :: unknown_field, coordinate_field, integral
  public :: WHOLE, TIME_FACE, TIME_FACTOR, SPACE_FACTOR

  ! the parts a support can be of another
  integer, parameter :: WHOLE = 0, TIME_FACE = 1, TIME_FACTOR = 2, SPACE_FACTOR = 3
  public :: operator(+), operator(-), operator(*), operator(/), operator(**)
  public :: sin, cos, exp, log, sqrt

  !===================================================================!
  ! THE SUPPORTS. A continuous support is known to a field by its
  ! identity alone; a discrete support by its points.
  !===================================================================!

  type, abstract :: continuous_support

     type(token) :: identity
     ! the part of another support this one is - the whole, the face
     ! at an instant, the time factor, the space factor - and the
     ! identity of that parent
     integer     :: part = WHOLE
     type(token) :: parent

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
     procedure(factor_interface)  , deferred :: design_factor
     procedure(count_interface)   , deferred :: design_order
     procedure(real_interface)    , deferred :: design_value

  end type discrete_support

  abstract interface

     ! the discretization of the design factor of the support: the
     ! point a functional of the solution takes its values on, with
     ! the order of the expansion along the design and the design's
     ! value
     function factor_interface(this) result(factor)
       import :: discrete_support
       class(discrete_support), intent(in) :: this
       class(discrete_support), allocatable :: factor
     end function factor_interface

     pure real(dp) function real_interface(this)
       import :: discrete_support, dp
       class(discrete_support), intent(in) :: this
     end function real_interface

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
  ! on its manifold and its arguments, the names of the coordinates
  ! it is a function of; a coordinate function has the coordinate it
  ! is, and is its own argument - the design coordinate reads the
  ! design of the point rather than a position; an integral has the
  ! identity of the factor it runs over, which part of its parent
  ! that factor is, and the parent's identity.
  !===================================================================!

  type :: continuous_field

     type(token) :: on
     character(len=:), allocatable :: name
     character(len=8), allocatable :: argument(:)
     type(expression), allocatable :: component(:)
     integer, allocatable :: index(:)
     integer :: coordinate = 0
     logical :: design_coordinate = .false.
     type(token) :: integrated_over
     integer     :: integrated_part = WHOLE
     type(token) :: integrated_parent

   contains

     procedure :: derivative     => field_derivative
     procedure :: discretize     => field_discretize
     procedure :: at             => field_at
     procedure :: differential   => field_differential
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
  ! component. The solution of a residual also stores the jet at
  ! every point - the value and the derivatives of every unknown, as
  ! the residual's rule lays them out, at jet(:, p, 1), and its
  ! derivatives along the design to the order of the expansion at
  ! jet(:, p, k + 1) - so that a functional of the solution reads the
  ! derivatives the residual formed; slot(j) is the row of the jet
  ! component j reads, zero for a component with no jet. The
  ! derivative along the design is the field of the next
  ! coefficients.
  !===================================================================!

  type :: discrete_field

     class(discrete_support), allocatable :: on
     character(len=32), allocatable :: name(:)
     real(dp), allocatable :: value(:,:)
     real(dp), allocatable :: jet(:,:,:)
     integer , allocatable :: slot(:)
     type(expression), allocatable :: rule

   contains

     procedure :: values
     procedure :: fields
     procedure :: derivative      => discrete_derivative
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
  ! per component, numbered from first, a function of the coordinates
  ! named by its arguments. A coordinate function. A tuple of fields
  ! on one support, each of one component. Constants, one per
  ! component.
  !===================================================================!

  function unknown_field(on, name, first, components, arguments) result(this)

    class(continuous_support), intent(in) :: on
    character(len=*)         , intent(in) :: name
    integer                  , intent(in) :: first, components
    character(len=*)         , intent(in) :: arguments(:)
    type(continuous_field) :: this

    integer :: k

    this % on   = on % identity
    this % name = name
    allocate(this % component(components), this % index(components))
    do k = 1, components
       this % index(k)     = first + k - 1
       this % component(k) = unknown(first + k - 1)
    end do
    allocate(this % argument(size(arguments)))
    do k = 1, size(arguments)
       this % argument(k) = arguments(k)
    end do

  end function unknown_field

  function coordinate_field(on, name, c, of_design) result(this)

    class(continuous_support), intent(in) :: on
    character(len=*)         , intent(in) :: name
    integer                  , intent(in) :: c
    logical                  , intent(in), optional :: of_design
    type(continuous_field) :: this

    this % on         = on % identity
    this % name       = name
    this % coordinate = c
    this % index      = [0]
    if (present(of_design)) this % design_coordinate = of_design
    if (this % design_coordinate) then
       this % component = [design()]
    else
       this % component = [coordinate(c)]
    end if
    allocate(this % argument(1))
    this % argument(1) = name

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
  ! marked with the factor's identity, its part and its parent. The
  ! discrete residual sums it over the points of that factor.
  !===================================================================!

  function integral(f, over) result(this)

    type(continuous_field)   , intent(in) :: f
    class(continuous_support), intent(in) :: over
    type(continuous_field) :: this

    this = f
    this % integrated_over   = over % identity
    this % integrated_part   = over % part
    this % integrated_parent = over % parent

  end function integral

  !===================================================================!
  ! THE VALUE OF A FUNCTIONAL ON A SOLUTION: the integral over the
  ! points of the support, with their measure, of the integrand
  ! evaluated on the jet at each point; a field on the design factor
  ! of the support. The differential is its gradient in the jet at
  ! every point, in the layout of the solution's rule: the source of
  ! the adjoint equation. Invalid input: a field that is not a scalar
  ! integral, a solution without a jet, or an integrand reading a
  ! derivative the rule does not store.
  !===================================================================!

  function field_at(this, solution) result(functional)

    class(continuous_field), intent(in) :: this
    type(discrete_field)   , intent(in) :: solution
    type(discrete_field) :: functional

    real(dp), allocatable :: total(:), gradient(:,:)
    character(len=32), allocatable :: names(:)
    integer :: k

    call functional_over(this, solution, total, gradient)
    allocate(names(size(total)))
    names(1) = ''
    if (allocated(this % name)) names(1) = this % name
    do k = 2, size(total)
       write(names(k), '(a,i0)') 'derivative ', k - 1
    end do
    functional = discrete_field(solution % on % design_factor(), reshape(total, [1, size(total)]), names)

  end function field_at

  subroutine field_differential(this, solution, gradient)

    class(continuous_field), intent(in)  :: this
    type(discrete_field)   , intent(in)  :: solution
    real(dp), allocatable  , intent(out) :: gradient(:,:)

    real(dp), allocatable :: total(:)

    call functional_over(this, solution, total, gradient)

  end subroutine field_differential

  !===================================================================!
  ! The value of the functional and its derivatives along the design
  ! to the order of the expansion, total(1:order+1), by the jet
  ! arithmetic: each component of the point's tuple enters as the
  ! quantity whose k-th derivative along the design is the solution's
  ! k-th coefficient, the design as the quantity with derivative one;
  ! and the gradient of the value in the jet at every point.
  !===================================================================!

  subroutine functional_over(this, solution, total, gradient)

    class(continuous_field), intent(in)  :: this
    type(discrete_field)   , intent(in)  :: solution
    real(dp), allocatable  , intent(out) :: total(:)
    real(dp), allocatable  , intent(out) :: gradient(:,:)

    real(dp), allocatable :: q(:), g(:)
    integer , allocatable :: slot(:)
    type(derivative_terms), allocatable :: jets(:)
    type(derivative_terms) :: nu, r
    real(dp) :: value, measure
    integer :: p, k, o, order

    if (.not. this % is_integral()) then
       error stop 'operation_field: a functional is the integral of a field over its manifold'
    end if
    if (size(this % component) /= 1) then
       error stop 'operation_field: a functional is a scalar field'
    end if
    if (.not. (allocated(solution % jet) .and. allocated(solution % rule))) then
       error stop 'operation_field: a functional is evaluated on the solution of a residual, which &
            &stores the jet at every point'
    end if
    order = size(solution % jet, 3) - 1
    call jet_slots(this % component(1), solution % rule, slot)
    allocate(q(0:size(slot) - 1), g(0:size(slot) - 1), jets(0:size(slot) - 1))
    allocate(gradient(size(solution % jet, 1), size(solution % jet, 2)), source=0.0_dp)
    allocate(total(order + 1), source=0.0_dp)
    nu = derivative_terms(solution % on % design_value(), order)
    if (order > 0) call nu % set_symmetric(1, 1.0_dp)
    do p = 1, size(solution % jet, 2)
       do k = 0, size(slot) - 1
          q(k)    = solution % jet(slot(k), p, 1)
          jets(k) = derivative_terms(q(k), order)
          do o = 1, order
             call jets(k) % set_symmetric(o, solution % jet(slot(k), p, o + 1))
          end do
       end do
       measure = solution % on % measure(p)
       call this % component(1) % gradient_at(q, solution % on % design_value(), value, g, &
            & solution % on % position(p))
       r = this % component(1) % at_instant(jets, nu, solution % on % position(p))
       do o = 0, order
          total(o + 1) = total(o + 1) + measure * coefficient(r, 2**o - 1)
       end do
       do k = 0, size(slot) - 1
          gradient(slot(k), p) = gradient(slot(k), p) + measure * g(k)
       end do
    end do

  end subroutine functional_over

  !===================================================================!
  ! THE DERIVATIVE OF A DISCRETE FIELD along the design coordinate,
  ! one order per entry of the multi-index: the field whose values
  ! are the next coefficients of the jet along the design, with the
  ! rest of the jet after them. Invalid input: an entry that is not
  ! the design coordinate, an order past the expansion's, or a
  ! component with no jet, a multiplier's.
  !===================================================================!

  function discrete_derivative(this, along) result(d)

    class(discrete_field)  , intent(in) :: this
    type(continuous_field) , intent(in) :: along(:)
    type(discrete_field) :: d

    integer :: n, k, order, npts
    character(len=250) :: message

    n = size(along)
    do k = 1, n
       if (.not. along(k) % design_coordinate) then
          error stop 'operation_field: a discrete field is differentiated along the design coordinate; &
               &its derivatives along time and space are the residual''s to form'
       end if
    end do
    if (.not. (allocated(this % jet) .and. allocated(this % slot))) then
       error stop 'operation_field: the derivative along the design is read from the solution of a &
            &residual discretized with an expansion along the design'
    end if
    order = size(this % jet, 3) - 1
    if (n > order) then
       write(message,'(a,i0,a,i0)') 'operation_field: the derivative of order ', n, ' along the design &
            &exceeds the order of the expansion, ', order
       error stop trim(message)
    end if
    npts = size(this % value, 1)
    allocate(d % on, source=this % on)
    d % name = this % name
    d % slot = this % slot
    allocate(d % value(npts, size(this % value, 2)))
    do k = 1, size(this % value, 2)
       if (this % slot(k) == 0) then
          error stop 'operation_field: component ' // trim(this % name(k)) // ' has no jet along the design'
       end if
       d % value(:, k) = this % jet(this % slot(k), :, n + 1)
    end do
    allocate(d % jet(size(this % jet, 1), npts, order - n + 1))
    d % jet = this % jet(:, :, n + 1:order + 1)
    if (allocated(this % rule)) d % rule = this % rule

  end function discrete_derivative

  !===================================================================!
  ! The slot in the rule's jet of each component an expression reads:
  ! the value and the derivatives along the first coordinate of each
  ! field, then its derivatives along each later coordinate, the two
  ! expressions laying their tuples out by their own degrees. A field
  ! or a derivative the expression reads and the rule does not store
  ! is invalid input.
  !===================================================================!

  subroutine jet_slots(g, rule, slot)

    type(expression), intent(in)  :: g, rule
    integer, allocatable, intent(out) :: slot(:)

    integer :: f, d, c, o, fields
    character(len=250) :: message

    allocate(slot(0:g % num_components() - 1), source=0)
    fields = rule % num_fields() - rule % num_multipliers()
    do f = 1, g % num_fields() - g % num_multipliers()
       if (f > fields) then
          write(message,'(a,i0,a,i0)') 'operation_field: the functional reads unknown ', f, &
               & ' and the rule stores ', fields
          error stop trim(message)
       end if
       if (g % degree_of_field(f) > rule % degree_of_field(f)) then
          write(message,'(a,i0,a,i0,a,i0)') 'operation_field: the functional reads the derivative of order ', &
               & g % degree_of_field(f), ' of unknown ', f, ' along the first coordinate, and the rule &
               &stores the order ', rule % degree_of_field(f)
          error stop trim(message)
       end if
       do d = 0, g % degree_of_field(f)
          slot(g % offset_of_field(f) + d) = rule % offset_of_field(f) + d + 1
       end do
       do c = FIRST_COORDINATE + 1, g % num_coordinates()
          if (c > rule % num_coordinates()) then
             write(message,'(a,i0)') 'operation_field: the functional reads a derivative along coordinate ', c
             error stop trim(message) // ', which the rule does not read'
          end if
          if (g % degree_along(c) > rule % degree_along(c)) then
             write(message,'(a,i0,a,i0,a,i0)') 'operation_field: the functional reads the derivative of order ', &
                  & g % degree_along(c), ' along coordinate ', c, ', and the rule stores the order ', &
                  & rule % degree_along(c)
             error stop trim(message)
          end if
          do o = 1, g % degree_along(c)
             slot(g % component_at(c, o, f)) = rule % component_at(c, o, f) + 1
          end do
       end do
    end do

  end subroutine jet_slots

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
  ! entry one order along that coordinate. Along a coordinate the
  ! unknown is not a function of, one outside its arguments, the
  ! derivative is the zero field. Invalid input: a field that is not
  ! a scalar unknown, an entry that is not a coordinate of the same
  ! manifold, or two distinct coordinates, since the jet stores no
  ! mixed component.
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
       if (along(k) % design_coordinate) then
          error stop 'operation_field: the derivative of a continuous field along the design coordinate &
               &is not stated; the derivatives along the design are read from the discrete solution, &
               &derivative([nu]) on a discrete field'
       end if
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
    if (allocated(this % argument)) then
       if (.not. any(this % argument == along(1) % name)) then
          d % component = [constant(0.0_dp)]
          return
       end if
    end if
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
  ! The components named, in the order named, every component of a
  ! name in its order - the components of a multiplier share its
  ! name. The jet is retained. A name not present stops the
  ! program.
  !===================================================================!

  function fields(this, names) result(chosen)

    class(discrete_field), intent(in) :: this
    character(len=*)     , intent(in) :: names(:)
    type(discrete_field) :: chosen

    integer :: k, j, count

    allocate(chosen % on, source=this % on)
    count = 0
    do k = 1, size(names)
       count = count + count_named(names(k))
    end do
    allocate(chosen % name(count), chosen % value(size(this % value, 1), count))
    count = 0
    do k = 1, size(names)
       if (count_named(names(k)) == 0) then
          error stop 'operation_field: fields names a component the field does not have: ' // trim(names(k))
       end if
       do j = 1, size(this % name)
          if (trim(this % name(j)) == trim(names(k))) then
             count = count + 1
             chosen % name(count)     = this % name(j)
             chosen % value(:, count) = this % value(:, j)
          end if
       end do
    end do
    if (allocated(this % jet))  chosen % jet  = this % jet
    if (allocated(this % rule)) chosen % rule = this % rule
    if (allocated(this % slot)) then
       allocate(chosen % slot(count))
       count = 0
       do k = 1, size(names)
          do j = 1, size(this % name)
             if (trim(this % name(j)) == trim(names(k))) then
                count = count + 1
                chosen % slot(count) = this % slot(j)
             end if
          end do
       end do
    end if

  contains

    pure integer function count_named(name)
      character(len=*), intent(in) :: name
      integer :: j
      count_named = 0
      do j = 1, size(this % name)
         if (trim(this % name(j)) == trim(name)) count_named = count_named + 1
      end do
    end function count_named

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
