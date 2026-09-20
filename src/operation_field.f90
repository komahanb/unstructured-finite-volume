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
  use util_derivative_terms, only : derivative_terms, coefficient, symmetric_terms
  use operation_expression, only : expression, unknown, constant, coordinate, derivative, derivative_along, design
  use operation_expression, only : partial_derivative, partial_in_design, total_derivative, stationarity
  use operation_expression, only : boundary_terms
  use operation_expression, only : FIRST_COORDINATE
  use operation_expression, only : operator(+), operator(-), operator(*), operator(/), operator(**)
  use operation_expression, only : sin, cos, exp, log, sqrt

  implicit none

  private
  public :: continuous_support, discrete_support
  public :: continuous_field, discrete_field
  public :: unknown_field, coordinate_field, integral
  public :: WHOLE, TIME_FACE, TIME_FACTOR, SPACE_FACTOR, DESIGN_FACTOR

  ! the parts a support can be of another
  integer, parameter :: WHOLE = 0, TIME_FACE = 1, TIME_FACTOR = 2, SPACE_FACTOR = 3, DESIGN_FACTOR = 4
  public :: operator(+), operator(-), operator(*), operator(/), operator(**)
  public :: sin, cos, exp, log, sqrt

  !===================================================================!
  ! THE SUPPORTS. A continuous support is known to a field by its
  ! identity alone; a discrete support by its points.
  !===================================================================!

  type, abstract :: continuous_support

     type(token) :: identity
     ! the part of another support this one is - the whole, the face
     ! at an instant, the time factor, the space factor, the design
     ! factor - the identity of that parent, and the instant of a face
     integer     :: part = WHOLE
     type(token) :: parent
     real(dp)    :: face_time = 0.0_dp

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
     procedure(count_interface)   , deferred :: num_designs
     procedure(values_interface)  , deferred :: design_values
     procedure(count_interface)   , deferred :: num_expansions
     procedure(seed_interface)    , deferred :: expansion_seed
     procedure(face_measure_interface), deferred :: measure_on_face

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

     ! the design vector of the support, one value per design coordinate
     pure function values_interface(this) result(v)
       import :: discrete_support, dp
       class(discrete_support), intent(in) :: this
       real(dp), allocatable :: v(:)
     end function values_interface

     ! the seed of the e-th expansion: the derivative of every design
     ! coordinate along it
     pure function seed_interface(this, e) result(seed)
       import :: discrete_support, dp
       class(discrete_support), intent(in) :: this
       integer                , intent(in) :: e
       real(dp), allocatable :: seed(:)
     end function seed_interface

     ! the measure of point p within the face at the instant given:
     ! zero for a point not on the face
     pure real(dp) function face_measure_interface(this, p, time)
       import :: discrete_support, dp
       class(discrete_support), intent(in) :: this
       integer                , intent(in) :: p
       real(dp)               , intent(in) :: time
     end function face_measure_interface

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
     ! a direction in the design: its vector over the design coordinates
     real(dp), allocatable :: direction(:)
     type(token) :: integrated_over
     integer     :: integrated_part = WHOLE
     type(token) :: integrated_parent
     real(dp)    :: integrated_time = 0.0_dp

   contains

     procedure :: derivative     => field_derivative
     procedure :: partial        => field_partial
     procedure :: stationarity   => field_stationarity
     procedure :: boundary_terms => field_boundary_terms
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
  ! the residual's rule lays them out, at jet(:, p, 1, d), and its
  ! derivatives along the d-th design coordinate to the order of the
  ! expansion at jet(:, p, k + 1, d) - so that a functional of the
  ! solution reads the derivatives the residual formed; slot(j) is
  ! the row of the jet component j reads, zero for a component with
  ! no jet. The derivative along a design coordinate is the field of
  ! the next coefficients along it, which retains the jet along that
  ! coordinate alone, recorded in along: a mixed derivative along two
  ! coordinates is not stored.
  !
  ! THE QUADRATURE NODES. The solution of a residual stores, beside
  ! its points, the nodes of the families' own quadrature: every
  ! instant at every cell with the weight the multistep rules give
  ! it, and every stage of a staged family's steps with the tableau's
  ! weight - node_jet(:, n, k + 1, d) the jet at the node, node_weight
  ! its weight, node_position its position, node_point the point it
  ! is (zero for a stage). A functional of the solution is the sum
  ! over the nodes, so that the objective is discretized by the
  ! scheme that marches the state, not by a rule of the points.
  !===================================================================!

  type :: discrete_field

     class(discrete_support), allocatable :: on
     character(len=32), allocatable :: name(:)
     real(dp), allocatable :: value(:,:)
     real(dp), allocatable :: jet(:,:,:,:)
     integer , allocatable :: slot(:)
     integer :: along = 0
     type(expression), allocatable :: rule
     real(dp), allocatable :: node_jet(:,:,:,:)
     real(dp), allocatable :: node_weight(:)
     real(dp), allocatable :: node_position(:,:)
     integer , allocatable :: node_point(:)
     ! the image that owns the node: a functional sums its owned nodes
     ! and the images' sums once
     integer , allocatable :: node_image(:)

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
       ! a direction, of index zero, is no leaf of an expression
       this % component = [design(max(c, 1))]
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
    this % integrated_time   = over % face_time

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

    real(dp), allocatable :: total(:), gradient(:,:,:,:), partial(:,:,:)
    character(len=32), allocatable :: names(:)
    integer :: k, d, o, order, nd

    call functional_over(this, solution, total, gradient, partial)
    order = size(solution % jet, 3) - 1
    nd    = size(solution % jet, 4)
    allocate(names(size(total)))
    names(1) = ''
    if (allocated(this % name)) names(1) = this % name
    k = 1
    do d = 1, nd
       do o = 1, order
          k = k + 1
          write(names(k), '(a,i0,a,i0)') 'expansion ', d, ' order ', o
       end do
    end do
    functional = discrete_field(solution % on % design_factor(), reshape(total, [1, size(total)]), names)

  end function field_at

  subroutine field_differential(this, solution, gradient, partial)

    class(continuous_field), intent(in)  :: this
    type(discrete_field)   , intent(in)  :: solution
    real(dp), allocatable  , intent(out) :: gradient(:,:,:,:)
    real(dp), allocatable  , intent(out) :: partial(:,:,:)

    real(dp), allocatable :: total(:)

    call functional_over(this, solution, total, gradient, partial)

  end subroutine field_differential

  !===================================================================!
  ! The value of the functional and its derivatives along each design
  ! coordinate to the order of the expansion - total(1) the value,
  ! then for each coordinate d its derivatives of order 1 to the
  ! order, by the jet arithmetic: each component of the point's tuple
  ! enters as the quantity whose k-th derivative along d is the
  ! solution's k-th coefficient along d, the coordinate d as the
  ! quantity with derivative one. The gradient of the value in the
  ! jet at every point, gradient(:, p, 0, d), and its derivatives
  ! along d, gradient(:, p, o, d), by one more direction seeded on the
  ! component varied; the partial derivative of the value in each
  ! design coordinate i at fixed state, partial(i, 0, d), by the
  ! reverse pass, and its derivatives along d, partial(i, o, d), by one
  ! more direction seeded on the coordinate i. Over a face at an
  ! instant the measure is the point's within the face, zero off it:
  ! J = u(T)^2 / 2 is the integral of that field over the face at T.
  !===================================================================!

  subroutine functional_over(this, solution, total, gradient, partial)

    class(continuous_field), intent(in)  :: this
    type(discrete_field)   , intent(in)  :: solution
    real(dp), allocatable  , intent(out) :: total(:)
    real(dp), allocatable  , intent(out) :: gradient(:,:,:,:)
    real(dp), allocatable  , intent(out) :: partial(:,:,:)

    real(dp), allocatable :: q(:), g(:), gd(:), values(:), unit(:), seed(:), position(:)
    integer , allocatable :: slot(:)
    type(derivative_terms), allocatable :: jets(:), varied(:), nu(:), along(:), design(:)
    type(derivative_terms) :: r
    real(dp) :: value, measure
    integer :: p, k, o, i, d, order, ns, nd, nnodes
    logical :: over_nodes

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
    ! over the quadrature nodes of the solution when it stores them,
    ! the families' own rule; over the points with their measure
    ! otherwise
    over_nodes = allocated(solution % node_weight)
    order  = size(solution % jet, 3) - 1
    nd     = size(solution % jet, 4)
    nnodes = size(solution % jet, 2)
    if (over_nodes) nnodes = size(solution % node_weight)
    values = solution % on % design_values()
    if (size(values) < this % component(1) % num_designs()) then
       error stop 'operation_field: the functional reads a design coordinate the manifold does not have'
    end if
    call jet_slots(this % component(1), solution % rule, slot)
    ns = size(slot)
    allocate(q(0:ns - 1), g(0:ns - 1), gd(size(values)), jets(0:ns - 1), varied(0:ns - 1))
    allocate(nu(size(values)), along(size(values)), design(size(values)))
    allocate(gradient(size(solution % jet, 1), nnodes, 0:order, nd), source=0.0_dp)
    allocate(total(1 + nd * order), partial(size(values), 0:order, nd), source=0.0_dp)
    allocate(unit(order), source=0.0_dp)
    if (order > 0) unit(1) = 1.0_dp

    do p = 1, nnodes
       if (over_nodes) then
          if (solution % node_image(p) /= this_image()) cycle
          measure  = solution % node_weight(p)
          position = solution % node_position(:, p)
          if (this % integrated_part == TIME_FACE) then
             measure = 0.0_dp
             if (solution % node_point(p) > 0) then
                measure = solution % on % measure_on_face(solution % node_point(p), this % integrated_time)
             end if
          end if
       else
          measure  = solution % on % measure(p)
          position = solution % on % position(p)
          if (this % integrated_part == TIME_FACE) measure = solution % on % measure_on_face(p, this % integrated_time)
       end if
       if (measure == 0.0_dp) cycle
       ! the value, its gradient in the jet and in the design
       do k = 0, ns - 1
          q(k) = jet_of(slot(k), p, 1, 1)
       end do
       call this % component(1) % gradient_at(q, values, value, g, position, gd)
       total(1) = total(1) + measure * value
       do d = 1, nd
          do k = 0, ns - 1
             gradient(slot(k), p, 0, d) = gradient(slot(k), p, 0, d) + measure * g(k)
          end do
          partial(:, 0, d) = partial(:, 0, d) + measure * gd
       end do
       if (order == 0) cycle
       ! along each expansion d: the jets over the order's directions,
       ! and over one more, the last, seeded on a component or on a
       ! coordinate
       do d = 1, nd
          seed = solution % on % expansion_seed(d)
          do i = 1, size(values)
             nu(i)    = symmetric_terms(values(i), seed(i) * unit, order)
             along(i) = symmetric_terms(values(i), seed(i) * unit, order + 1)
          end do
          do k = 0, ns - 1
             jets(k)   = symmetric_terms(q(k), [(jet_of(slot(k), p, o + 1, d), o = 1, order)], order)
             varied(k) = symmetric_terms(q(k), [(jet_of(slot(k), p, o + 1, d), o = 1, order)], order + 1)
          end do
          r = this % component(1) % at_instant(jets, nu, position)
          do o = 1, order
             total(1 + (d - 1) * order + o) = total(1 + (d - 1) * order + o) + measure * coefficient(r, 2**o - 1)
          end do
          do i = 1, size(values)
             design = along
             call design(i) % set_coefficient(2**order, 1.0_dp)
             r = this % component(1) % at_instant(varied, design, position)
             do o = 1, order
                partial(i, o, d) = partial(i, o, d) + measure * coefficient(r, 2**o - 1 + 2**order)
             end do
          end do
          do k = 0, ns - 1
             call varied(k) % set_coefficient(2**order, 1.0_dp)
             r = this % component(1) % at_instant(varied, along, position)
             call varied(k) % set_coefficient(2**order, 0.0_dp)
             do o = 1, order
                gradient(slot(k), p, o, d) = gradient(slot(k), p, o, d) + measure * coefficient(r, 2**o - 1 + 2**order)
             end do
          end do
       end do
    end do
    ! the sums over the owned nodes, summed over the images once
    if (over_nodes) then
       if (any(solution % node_image /= this_image())) then
          call co_sum(total)
          call co_sum(gradient)
          call co_sum(partial)
       end if
    end if

  contains

    ! the jet coefficient at a node or a point
    pure real(dp) function jet_of(row, n, k, d)
      integer, intent(in) :: row, n, k, d
      if (over_nodes) then
         jet_of = solution % node_jet(row, n, k, d)
      else
         jet_of = solution % jet(row, n, k, d)
      end if
    end function jet_of

  end subroutine functional_over

  ! the quadrature nodes of a solution, copied to a field read from it
  subroutine nodes_copied(from, to)
    type(discrete_field), intent(in)    :: from
    type(discrete_field), intent(inout) :: to
    if (allocated(from % node_weight)) then
       to % node_jet      = from % node_jet
       to % node_weight   = from % node_weight
       to % node_position = from % node_position
       to % node_point    = from % node_point
       to % node_image    = from % node_image
    end if
  end subroutine nodes_copied

  !===================================================================!
  ! THE DERIVATIVE OF A DISCRETE FIELD along the design: the
  ! multi-index names design coordinates or directions, and the
  ! derivative is read from the jets stored - along the expansion
  ! whose seed is the coordinate's unit vector or the direction's
  ! vector, for a multi-index repeating one entry: the field whose
  ! values are the next coefficients, with the rest of that jet after
  ! them and the jets along the other expansions dropped. A mixed
  ! derivative along two coordinates, [nu, mu], is the polarization
  ! (u_ww - u_nunu - u_mumu) / 2 of the second derivatives along nu,
  ! mu and w = nu + mu, and requires the expansions along all three;
  ! it stores no jet of its own. Invalid input: an entry that is not
  ! a design coordinate or direction, an expansion not stored, a
  ! coordinate other than the one the field is already a derivative
  ! along, an order past the expansion's, a mixed derivative of an
  ! order other than two, or a component with no jet.
  !===================================================================!

  function discrete_derivative(this, along) result(d)

    class(discrete_field)  , intent(in) :: this
    type(continuous_field) , intent(in) :: along(:)
    type(discrete_field) :: d

    real(dp), allocatable :: seed(:), other(:), both(:)
    integer :: n, k, order, npts, column, second, third
    logical :: mixed
    character(len=250) :: message

    n = size(along)
    do k = 1, n
       if (.not. along(k) % design_coordinate) then
          error stop 'operation_field: a discrete field is differentiated along a design coordinate or a &
               &direction in the design; its derivatives along time and space are the residual''s to form'
       end if
    end do
    if (.not. (allocated(this % jet) .and. allocated(this % slot))) then
       error stop 'operation_field: the derivative along the design is read from the solution of a &
            &residual discretized with an expansion along the design'
    end if
    order = size(this % jet, 3) - 1
    npts  = size(this % value, 1)
    seed  = seed_of(along(1))

    ! the mixed derivative along two coordinates, by polarization
    mixed = .false.
    if (n == 2) mixed = .not. same_seed(seed, seed_of(along(2)))
    if (mixed) then
       if (this % along /= 0) then
          error stop 'operation_field: a mixed derivative is read from the solution, not from a derivative of it'
       end if
       if (order < 2) then
          error stop 'operation_field: a mixed derivative along two coordinates requires an expansion of order two'
       end if
       other  = seed_of(along(2))
       both   = seed + other
       column = expansion_of(seed)
       second = expansion_of(other)
       third  = expansion_of(both)
       call require_columns(column, second, third)
       allocate(d % on, source=this % on)
       d % name  = this % name
       d % slot  = this % slot
       d % along = 0
       allocate(d % value(npts, size(this % value, 2)))
       do k = 1, size(this % value, 2)
          call require_jet(k)
          d % value(:, k) = (this % jet(this % slot(k), :, 3, third) - this % jet(this % slot(k), :, 3, column) &
               & - this % jet(this % slot(k), :, 3, second)) / 2.0_dp
       end do
       allocate(d % jet(size(this % jet, 1), npts, 1, 1))
       d % jet(:, :, 1, 1) = (this % jet(:, :, 3, third) - this % jet(:, :, 3, column) - this % jet(:, :, 3, second)) / 2.0_dp
       if (allocated(this % rule)) d % rule = this % rule
       call nodes_copied(this, d)
       return
    end if

    do k = 2, n
       if (.not. same_seed(seed, seed_of(along(k)))) then
          error stop 'operation_field: a mixed derivative is along two coordinates, of order two; the &
               &multi-index otherwise repeats one coordinate or direction'
       end if
    end do
    if (this % along /= 0) then
       if (.not. same_seed(seed, this % on % expansion_seed(this % along))) then
          error stop 'operation_field: the field is a derivative along one expansion, and its jet along &
               &another is not stored'
       end if
       column = 1
    else
       column = expansion_of(seed)
       if (column == 0) then
          error stop 'operation_field: no expansion of the solution is along the coordinate or direction named; &
               &state it in expansion(order, along=...)'
       end if
    end if
    if (n > order) then
       write(message,'(a,i0,a,i0)') 'operation_field: the derivative of order ', n, ' along the design &
            &exceeds the order of the expansion, ', order
       error stop trim(message)
    end if
    allocate(d % on, source=this % on)
    d % name  = this % name
    d % slot  = this % slot
    d % along = merge(this % along, column, this % along /= 0)
    allocate(d % value(npts, size(this % value, 2)))
    do k = 1, size(this % value, 2)
       call require_jet(k)
       d % value(:, k) = this % jet(this % slot(k), :, n + 1, column)
    end do
    allocate(d % jet(size(this % jet, 1), npts, order - n + 1, 1))
    d % jet(:, :, :, 1) = this % jet(:, :, n + 1:order + 1, column)
    if (allocated(this % rule)) d % rule = this % rule
    call nodes_copied(this, d)

  contains

    ! the seed a coordinate or a direction names: the unit vector of
    ! the coordinate over the design coordinates, or the vector
    function seed_of(entry) result(s)
      type(continuous_field), intent(in) :: entry
      real(dp), allocatable :: s(:)
      integer :: nd
      if (allocated(entry % direction)) then
         s = entry % direction
         return
      end if
      nd = this % on % num_designs()
      allocate(s(nd), source=0.0_dp)
      if (entry % coordinate < 1 .or. entry % coordinate > nd) then
         error stop 'operation_field: the design coordinate named is not one of the solution''s manifold'
      end if
      s(entry % coordinate) = 1.0_dp
    end function seed_of

    pure logical function same_seed(a, b)
      real(dp), intent(in) :: a(:), b(:)
      same_seed = size(a) == size(b)
      if (same_seed) same_seed = all(abs(a - b) <= 1.0e-12_dp * max(1.0_dp, maxval(abs(a))))
    end function same_seed

    ! the expansion stored along a seed, zero when none is
    integer function expansion_of(s)
      real(dp), intent(in) :: s(:)
      integer :: e
      expansion_of = 0
      do e = 1, this % on % num_expansions()
         if (same_seed(s, this % on % expansion_seed(e))) expansion_of = e
      end do
    end function expansion_of

    subroutine require_columns(a, b, c)
      integer, intent(in) :: a, b, c
      if (a == 0 .or. b == 0 .or. c == 0) then
         error stop 'operation_field: a mixed derivative along nu and mu requires the expansions along nu, &
              &along mu and along nu + mu: expansion(order=2, along=[e_nu, e_mu, e_nu + e_mu])'
      end if
    end subroutine require_columns

    subroutine require_jet(k)
      integer, intent(in) :: k
      if (this % slot(k) == 0) then
         error stop 'operation_field: component ' // trim(this % name(k)) // ' has no jet along the design'
      end if
    end subroutine require_jet

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

    if (size(this % component) /= 1) then
       error stop 'operation_field: derivative requires a scalar field'
    end if
    if (size(along) < 1) then
       error stop 'operation_field: derivative requires one coordinate at least'
    end if
    ! a field that is not an unknown: its total derivative along the
    ! coordinates, by the chain rule over the jets it reads
    if (this % index(1) < 1) then
       d = this
       d % index = [0]
       do k = 1, size(along)
          if (along(k) % design_coordinate .or. along(k) % coordinate < 1) then
             error stop 'operation_field: the total derivative of a field is along a coordinate of time or space'
          end if
          if (.not. along(k) % on % matches(this % on)) then
             error stop 'operation_field: the total derivative is along a coordinate of the field''s manifold'
          end if
          d % component = [total_derivative(d % component(1), along(k) % coordinate)]
       end do
       return
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
  ! THE PARTIAL DERIVATIVE of a scalar field in a component of the
  ! jet of an unknown - u, or u % derivative([t]) - or in a design
  ! coordinate, every other leaf unchanged: the symbolic derivative as a
  ! field on the same manifold. Invalid input: a field of several
  ! components, or an argument that is neither a jet component of an
  ! unknown nor a design coordinate.
  !===================================================================!

  function field_partial(this, in) result(d)

    class(continuous_field), intent(in) :: this
    type(continuous_field) , intent(in) :: in
    type(continuous_field) :: d

    type(expression) :: leaf

    if (size(this % component) /= 1 .or. size(in % component) /= 1) then
       error stop 'operation_field: the partial derivative is of a scalar field in a scalar component'
    end if
    d = this
    d % index = [0]
    if (in % design_coordinate) then
       if (allocated(in % direction)) then
          error stop 'operation_field: the partial derivative is in a design coordinate, not a direction'
       end if
       d % component = [partial_in_design(this % component(1), in % coordinate)]
       return
    end if
    leaf = in % component(1)
    if (leaf % num_vertices() /= 1 .or. .not. leaf % reads_state()) then
       error stop 'operation_field: the partial derivative is in a component of the jet of an unknown, &
            &u or u % derivative([t]), or in a design coordinate'
    end if
    d % component = [partial_derivative(this % component(1), leaf % leaf_field(), leaf % leaf_order(), &
         & leaf % leaf_along())]

  end function field_partial

  !===================================================================!
  ! THE STATIONARITY of a scalar field in an unknown along a
  ! coordinate: the Euler-Lagrange derivative, the sum over the
  ! orders o of (-1)^o times the o-th total derivative along the
  ! coordinate of the partial in the component of order o - the
  ! interior equation of the stationarity of the integral of the
  ! field, the boundary terms of the integration by parts aside.
  ! Invalid input: a field of several components, an argument that
  ! is not a scalar unknown, or a coordinate that is not of the
  ! manifold.
  !===================================================================!

  function field_stationarity(this, in, along) result(d)

    class(continuous_field), intent(in) :: this
    type(continuous_field) , intent(in) :: in, along
    type(continuous_field) :: d

    if (size(this % component) /= 1) then
       error stop 'operation_field: the stationarity is of a scalar field'
    end if
    if (size(in % component) /= 1 .or. in % index(1) < 1) then
       error stop 'operation_field: the stationarity is in a scalar unknown of the manifold'
    end if
    if (along % design_coordinate .or. along % coordinate < 1) then
       error stop 'operation_field: the stationarity integrates by parts along a coordinate of time or space'
    end if
    d = this
    d % index = [0]
    d % component = [stationarity(this % component(1), in % index(1), along % coordinate)]

  end function field_stationarity

  !===================================================================!
  ! THE BOUNDARY TERMS of the stationarity of a scalar field in an
  ! unknown along a coordinate: the fields whose vanishing at the far
  ! end of the coordinate are the natural conditions of the
  ! stationarity of the integral, one per order of the unknown below
  ! the highest, terms(k + 1) the coefficient of the variation of the
  ! k-th derivative. Invalid input as for the stationarity.
  !===================================================================!

  function field_boundary_terms(this, in, along) result(terms)

    class(continuous_field), intent(in) :: this
    type(continuous_field) , intent(in) :: in, along
    type(continuous_field), allocatable :: terms(:)

    type(expression), allocatable :: graphs(:)
    integer :: k

    if (size(this % component) /= 1) then
       error stop 'operation_field: the boundary terms are of a scalar field'
    end if
    if (size(in % component) /= 1 .or. in % index(1) < 1) then
       error stop 'operation_field: the boundary terms are in a scalar unknown of the manifold'
    end if
    if (along % design_coordinate .or. along % coordinate < 1) then
       error stop 'operation_field: the boundary terms are along a coordinate of time or space'
    end if
    graphs = boundary_terms(this % component(1), in % index(1), along % coordinate)
    allocate(terms(size(graphs)))
    do k = 1, size(graphs)
       terms(k) = this
       terms(k) % index = [0]
       terms(k) % component = [graphs(k)]
    end do

  end function field_boundary_terms

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
    chosen % along = this % along
    call nodes_copied(this, chosen)
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
