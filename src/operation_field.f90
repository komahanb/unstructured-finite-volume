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
  use operation_expansion , only : expansion
  use operation_expression, only : expression, unknown, constant, coordinate, component, component_along, design
  use operation_expression, only : derivative, jet_form
  use operation_expression, only : partial_derivative, partial_in_design, total_derivative, stationarity
  use operation_expression, only : boundary_terms
  use operation_expression, only : FIRST_COORDINATE
  use operation_expression, only : operator(+), operator(-), operator(*), operator(/), operator(**)
  use operation_expression, only : sin, cos, exp, log, sqrt, heaviside
  use view_expression     , only : expression_view

  implicit none

  private
  public :: continuous_support, discrete_support
  public :: continuous_field, discrete_field
  public :: unknown_field, coordinate_field, integral, mean, maximum, minimum, representation
  public :: WHOLE, TIME_FACE, TIME_FACTOR, SPACE_FACTOR, DESIGN_FACTOR, SPACE_FACE, TIME_WINDOW

  ! the parts a support can be of another
  integer, parameter :: WHOLE = 0, TIME_FACE = 1, TIME_FACTOR = 2, SPACE_FACTOR = 3, DESIGN_FACTOR = 4, &
       & SPACE_FACE = 5, TIME_WINDOW = 6
  public :: sin, cos, exp, log, sqrt, heaviside

  !===================================================================!
  ! THE SUPPORTS. A continuous support is known to a field by its
  ! identity alone; a discrete support by its points.
  !===================================================================!

  type, abstract :: continuous_support

     type(token) :: identity
     ! the part of another support this one is - the whole, the face
     ! at an instant, the face at a value of a space coordinate, the
     ! time factor, the space factor, the design factor, the window
     ! between two instants - the identity of that parent, and for a
     ! face the coordinate it is taken along and the coordinate's value
     ! there; for a window its first instant, face_value, and its last,
     ! window_end
     integer     :: part = WHOLE
     type(token) :: parent
     integer     :: face_axis  = 0
     real(dp)    :: face_value = 0.0_dp
     real(dp)    :: window_end = 0.0_dp

   contains

     procedure :: same_support
     procedure :: measure => support_measure

  end type continuous_support

  type, abstract :: discrete_support

     type(token) :: identity        ! of the continuous support discretized
     ! the discrete image of the design: the expansion of an order at
     ! the design values, of order zero and no seed without a design
     type(expansion) :: expansion

   contains

     procedure(count_interface)   , deferred :: num_points
     procedure(count_interface)   , deferred :: num_coordinates
     procedure(position_interface), deferred :: position
     procedure(measure_interface) , deferred :: measure
     procedure(factor_interface)  , deferred :: design_factor
     procedure(face_measure_interface), deferred :: measure_on_face

  end type discrete_support

  abstract interface

     ! the discretization of the design factor of the support: the
     ! point a functional of the solution takes its values on, with
     ! the expansion along the design
     function factor_interface(this) result(factor)
       import :: discrete_support
       class(discrete_support), intent(in) :: this
       class(discrete_support), allocatable :: factor
     end function factor_interface

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
     ! the instant of a face, or the first instant of a window and its
     ! last, integrated_end
     real(dp)    :: integrated_time = 0.0_dp
     real(dp)    :: integrated_end  = 0.0_dp
     ! the weight of a smooth maximum over the factor; zero for an
     ! integral
     real(dp)    :: weight = 0.0_dp
     ! A FIELD OF TIME WITH AN INTEGRAL OVER THE REGION, on a manifold
     ! of time and space: its value the component plus the integral over
     ! the region of the integrand, the sum of a field and of integrals
     ! over the space factor; such an integral enters a field linearly,
     ! times fields of time alone
     type(expression), allocatable :: integrand(:)

   contains

     procedure :: derivative     => field_derivative
     procedure, private :: jet_component => field_jet_component
     procedure :: jet_form       => field_jet_form
     procedure :: partial        => field_partial
     procedure, private :: field_stationarity, field_stationarity_along
     generic   :: stationarity   => field_stationarity, field_stationarity_along
     procedure :: boundary_terms => field_boundary_terms
     procedure :: discretize     => field_discretize
     procedure :: at             => field_at
     procedure :: differential   => field_differential
     procedure :: functional     => functional_over
     procedure :: num_components => field_num_components
     procedure :: graph          => field_graph
     procedure :: is_unknown
     procedure :: is_coordinate
     procedure :: is_integral
     ! THE ARITHMETIC of scalar fields, elemental over arrays of them
     generic   :: operator(+)  => plus, real_plus, plus_real
     generic   :: operator(-)  => minus, real_minus, minus_real, negated
     generic   :: operator(*)  => times, real_times, times_real
     generic   :: operator(/)  => over, real_over, over_real
     generic   :: operator(**) => to_integer, to_real
     procedure, private :: plus, minus, times, over, plus_real, minus_real, times_real, over_real, negated
     procedure, private :: to_integer, to_real
     procedure, private, pass(b) :: real_plus, real_minus, real_times, real_over

  end type continuous_field

  interface continuous_field
     module procedure create_tuple
     module procedure create_constants
  end interface continuous_field

  !===================================================================!
  ! THE DISCRETE FIELD: its support, and one value per point and
  ! component. The solution of a residual also stores the jet at
  ! every point - the value and the derivatives of every unknown, as
  ! the residual's rule lays them out - in every column of the
  ! expansion along the design, jet(:, p, c + 1) the column c, so that
  ! a functional of the solution reads the derivatives the residual
  ! formed; slot(j) is the row of the jet component j reads, zero for
  ! a component with no jet. The derivative along a design coordinate
  ! is the field of the next column along it, which retains the
  ! columns of that expansion alone, from the derivative's own on, the
  ! expansion recorded in along: a mixed derivative along two
  ! coordinates is not stored.
  !
  ! THE QUADRATURE NODES. The solution of a residual stores, beside
  ! its points, the nodes of the families' own quadrature: every
  ! instant at every cell with the weight the multistep rules give
  ! it, and every stage of a staged family's steps with the tableau's
  ! weight - node_jet(:, n, c + 1) the jet at the node, node_weight
  ! its weight, node_position its position, node_face_weight its
  ! weight within the face of its instant (zero for a stage, and for
  ! an instant of history). A functional of the solution is the sum
  ! over the nodes, so that the objective is discretized by the
  ! scheme that marches the state, not by a rule of the points; over a
  ! window between two instants, the sum with the weights of the steps
  ! within the window, a multistep step's history included.
  !===================================================================!

  type :: discrete_field

     class(discrete_support), allocatable :: on
     character(len=32), allocatable :: name(:)
     real(dp), allocatable :: value(:,:)
     real(dp), allocatable :: jet(:,:,:)
     integer , allocatable :: slot(:)
     integer :: along = 0
     type(expression), allocatable :: rule
     real(dp), allocatable :: node_jet(:,:,:)
     ! the layout of the nodes' jets where it is not the rule's: under
     ! a projection on the cells, the flux's, with the derivatives along
     ! space its basis forms at the quadrature points
     type(expression), allocatable :: node_rule
     real(dp), allocatable :: node_weight(:)
     real(dp), allocatable :: node_position(:,:)
     real(dp), allocatable :: node_face_weight(:)
     ! the shares of the nodes' weights by the steps of the march: the
     ! step that ends at the instant share_end(s) gives the node
     ! share_node(s) the weight share_weight(s); a node's weight is the
     ! sum of its shares, and within a window the sum of those of the
     ! steps within it
     integer , allocatable :: share_node(:)
     real(dp), allocatable :: share_end(:), share_weight(:)
     ! the image that owns the node: a functional sums its owned nodes
     ! and the images' sums once
     integer , allocatable :: node_image(:)
     ! a solution whose stages were retained over checkpoints alone
     ! integrates the objective of its residual and no other
     ! functional: the value and its derivatives along the design,
     ! summed segment by segment, and the objective's formula
     real(dp), allocatable :: objective_value(:)
     character(len=:), allocatable :: objective_formula

   contains

     procedure :: values
     procedure :: fields
     procedure :: derivative      => discrete_derivative
     procedure :: norm
     procedure :: num_points     => discrete_num_points
     procedure :: num_components => discrete_num_components
     generic   :: operator(-)    => discrete_minus
     procedure, private :: discrete_minus

  end type discrete_field

  interface discrete_field
     module procedure create_discrete
  end interface discrete_field

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

  interface heaviside
     module procedure step_of
  end interface heaviside

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
  ! THE REDUCTIONS of a field over a factor of its manifold - the whole,
  ! a factor or a face - and over its components: the integral, the
  ! mean, the smooth maximum and the smooth minimum. Each is a mean of
  ! one generator,
  !
  !      phi^-1( integral of phi(f) over the support ),
  !
  ! the integral's phi the identity, the mean's the identity over the
  ! measure made one, the smooth maximum's the exponential exp(rho f):
  ! the integral of the semiring in which a sum is log(exp(rho a) +
  ! exp(rho b)) / rho, which is the maximum as rho grows and, of
  ! negative weight, the minimum. So the reductions over a product of
  ! factors are the reductions over each in turn, of one generator,
  ! and a reduction of a reduction of another generator is not yet
  ! stated: a field reduced already stops the program.
  !===================================================================!

  ! the integral over the support and over the components, their sum
  function integral(f, over) result(this)

    type(continuous_field)   , intent(in) :: f
    class(continuous_support), intent(in) :: over
    type(continuous_field) :: this

    this = reduced(summed(f), over)

  end function integral

  ! the mean over the support and over the components: the integral of
  ! the components' mean over the support's measure, so that the mean
  ! of one is one, and a field of time alone on a manifold of time and
  ! space has its mean over the time
  function mean(f, over) result(this)

    type(continuous_field)   , intent(in) :: f
    class(continuous_support), intent(in) :: over
    type(continuous_field) :: this

    this = integral((1.0_dp / (real(size(f % component), dp) * over % measure())) * summed(f), over)

  end function mean

  ! the field marked with the factor's identity, its part and its
  ! parent; the discrete residual sums it over the points of that
  ! factor
  function reduced(f, over) result(this)

    type(continuous_field)   , intent(in) :: f
    class(continuous_support), intent(in) :: over
    type(continuous_field) :: this

    if (f % is_integral() .or. allocated(f % integrand)) then
       error stop 'operation_field: a reduction of a field reduced already is not yet stated; a reduction over &
            &the product of the factors is'
    end if
    this = f
    this % integrated_over   = over % identity
    this % integrated_part   = over % part
    this % integrated_parent = over % parent
    this % integrated_time   = over % face_value
    this % integrated_end    = over % window_end

  end function reduced

  ! the sum of the components of a field, a field of one
  function summed(f) result(this)

    type(continuous_field), intent(in) :: f
    type(continuous_field) :: this

    type(expression) :: total
    integer :: i

    this = f
    if (size(f % component) <= 1) return
    total = f % component(1)
    do i = 2, size(f % component)
       total = total + f % component(i)
    end do
    this % component = [total]
    this % index     = [0]

  end function summed

  ! the measure of a support - the duration of its time times the
  ! measure of its region - which a mean divides by; a support that
  ! does not state it stops the program
  function support_measure(this) result(m)

    class(continuous_support), intent(in) :: this
    real(dp) :: m

    m = 0.0_dp
    error stop 'operation_field: this support does not state its measure'

  end function support_measure

  !===================================================================!
  ! THE SMOOTH MAXIMUM of the components of a field over a factor of
  ! its manifold, of the weight rho > 0, the aggregate of Kreisselmeier
  ! and Steinhauser:
  !
  !      log( integral of sum_i exp(rho f_i) ) / rho,
  !
  ! which tends to the largest f_i over the factor as rho grows. On a
  ! solution it is m + log( integral of sum_i exp(rho (f_i - m)) ) /
  ! rho, m the largest f_i at the nodes, so that no exponential
  ! exceeds one; along the design, the derivative of the logarithm.
  ! THE SMOOTH MINIMUM is the same of the weight -rho, which tends to
  ! the smallest f_i, m then the smallest f_i.
  !===================================================================!

  function maximum(f, over, weight) result(this)

    type(continuous_field)   , intent(in) :: f
    class(continuous_support), intent(in) :: over
    real(dp)                 , intent(in) :: weight
    type(continuous_field) :: this

    if (.not. weight > 0.0_dp) error stop 'operation_field: the weight of a smooth maximum is positive'
    if (allocated(f % integrand)) error stop 'operation_field: a smooth maximum is of a field with no integral'
    this = reduced(f, over)
    this % weight = weight

  end function maximum

  function minimum(f, over, weight) result(this)

    type(continuous_field)   , intent(in) :: f
    class(continuous_support), intent(in) :: over
    real(dp)                 , intent(in) :: weight
    type(continuous_field) :: this

    if (.not. weight > 0.0_dp) error stop 'operation_field: the weight of a smooth minimum is positive'
    if (allocated(f % integrand)) error stop 'operation_field: a smooth minimum is of a field with no integral'
    this = reduced(f, over)
    this % weight = -weight

  end function minimum

  !===================================================================!
  ! THE VALUE OF A FUNCTIONAL ON A SOLUTION: the integral over the
  ! points of the support, with their measure, of the integrand
  ! evaluated on the jet at each point; a field on the design factor
  ! of the support. The differential is its gradient in the jet at
  ! every point, in the layout of the solution's rule: the source of
  ! the adjoint equation. A field that is not an integral is evaluated
  ! on the jet at every point, at the design the jet is expanded
  ! about: a field of the solution's support, one value per point and
  ! component. Invalid input: an integral of several components, a
  ! solution without a jet, or a field reading a derivative the rule
  ! does not store.
  !===================================================================!

  function field_at(this, solution) result(functional)

    class(continuous_field), intent(in) :: this
    type(discrete_field)   , intent(in) :: solution
    type(discrete_field) :: functional

    real(dp), allocatable :: total(:), gradient(:,:,:), partial(:,:)
    character(len=32), allocatable :: names(:)
    integer :: d, o

    if (.not. this % is_integral()) then
       functional = field_on_points(this, solution)
       return
    end if
    if (this % weight /= 0.0_dp) then
       total = maximum_over(this, solution)
    else
       call functional_over(this, solution, total, gradient, partial)
    end if
    allocate(names(size(total)))
    names(1) = ''
    if (allocated(this % name)) names(1) = this % name
    associate (expanded => solution % on % expansion)
      do d = 1, expanded % num_expansions()
         do o = 1, expanded % order
            write(names(1 + expanded % column(d, o)), '(a,i0,a,i0)') 'expansion ', d, ' order ', o
         end do
      end do
    end associate
    functional = discrete_field(solution % on % design_factor(), reshape(total, [1, size(total)]), names)

  end function field_at

  ! a field that is not an integral, each component evaluated on the
  ! jet at every point of the solution
  function field_on_points(this, solution) result(on_points)

    class(continuous_field), intent(in) :: this
    type(discrete_field)   , intent(in) :: solution
    type(discrete_field) :: on_points

    type(expression) :: g
    type(expansion)  :: expanded
    integer , allocatable :: slot(:)
    real(dp), allocatable :: q(:), gradient(:), gd(:), designs(:), values(:,:)
    character(len=32), allocatable :: names(:)
    integer :: j, p, k

    if (.not. (allocated(solution % rule) .and. allocated(solution % jet))) then
       error stop 'operation_field: a field is evaluated on the solution of a residual, which stores the jet &
            &at every point'
    end if
    expanded = jet_expansion(solution)
    designs  = expanded % values()
    allocate(values(size(solution % jet, 2), size(this % component)))
    allocate(names(size(this % component)))
    names = ''
    if (allocated(this % name)) names = this % name
    do j = 1, size(this % component)
       g = jet_form(this % component(j))
       call jet_slots(g, solution % rule, slot)
       if (size(designs) < g % num_designs()) then
          error stop 'operation_field: the field reads a design coordinate the manifold does not have'
       end if
       if (allocated(q)) deallocate(q, gradient, gd)
       allocate(q(0:size(slot) - 1), gradient(0:size(slot) - 1), gd(size(designs)))
       do p = 1, size(solution % jet, 2)
          do k = 0, size(slot) - 1
             q(k) = solution % jet(slot(k), p, 1)
          end do
          call g % gradient_at(q, designs, values(p, j), gradient, solution % on % position(p), gd)
       end do
    end do
    on_points = discrete_field(solution % on, values, names)

  end function field_on_points

  subroutine field_differential(this, solution, gradient, partial)

    class(continuous_field), intent(in)  :: this
    type(discrete_field)   , intent(in)  :: solution
    real(dp), allocatable  , intent(out) :: gradient(:,:,:)
    real(dp), allocatable  , intent(out) :: partial(:,:)

    real(dp), allocatable :: total(:)
    real(dp) :: m, scale
    integer  :: c

    if (this % weight /= 0.0_dp) then
       ! of the smooth maximum m + log(I) / rho: the differential of I
       ! over rho I, and along each expansion the derivative of that
       ! quotient, I's own derivative read in its value's column
       call maximum_parts(this, solution, m, total, gradient, partial)
       scale = 1.0_dp / (this % weight * total(1))
       gradient(:, :, 0) = gradient(:, :, 0) * scale
       partial(:, 0)     = partial(:, 0) * scale
       do c = 1, ubound(gradient, 3)
          gradient(:, :, c) = gradient(:, :, c) * scale - gradient(:, :, 0) * (total(1 + c) / total(1))
          partial(:, c)     = partial(:, c) * scale - partial(:, 0) * (total(1 + c) / total(1))
       end do
       return
    end if
    call functional_over(this, solution, total, gradient, partial)

  end subroutine field_differential

  !===================================================================!
  ! The value of the functional in every column of the expansion -
  ! total(c + 1) the column c: the value, then for each expansion d
  ! its derivatives of order 1 to the order, by the jet arithmetic:
  ! each component of the point's tuple enters as the quantity whose
  ! k-th derivative along d is the solution's column of order k along
  ! d, the coordinates as the quantities with the seed's derivatives.
  ! The gradient of the value in the jet at every point,
  ! gradient(:, p, 0), and its derivatives along d in the columns of
  ! d, by one more direction seeded on the component varied; the
  ! partial derivative of the value in each design coordinate i at
  ! fixed state, partial(i, 0), by the reverse pass, and its
  ! derivatives along d in the columns of d, by one more direction
  ! seeded on the coordinate i. Over a face at an
  ! instant the measure is the point's within the face, zero off it:
  ! J = u(T)^2 / 2 is the integral of that field over the face at T.
  !===================================================================!

  subroutine functional_over(this, solution, total, gradient, partial)

    class(continuous_field), intent(in)  :: this
    type(discrete_field)   , intent(in)  :: solution
    real(dp), allocatable  , intent(out) :: total(:)
    real(dp), allocatable  , intent(out) :: gradient(:,:,:)
    real(dp), allocatable  , intent(out) :: partial(:,:)

    real(dp), allocatable :: q(:), g(:), gd(:), values(:), unit(:), seed(:), position(:), within(:)
    integer , allocatable :: slot(:)
    type(derivative_terms), allocatable :: jets(:), varied(:), nu(:), along(:), design(:)
    type(derivative_terms) :: r
    type(expansion) :: expanded
    real(dp) :: value, measure
    integer :: p, k, o, i, d, order, ns, nnodes, stride, first
    logical :: over_nodes
    type(expression_view) :: view
    type(expression) :: integrand

    if (.not. this % is_integral()) then
       error stop 'operation_field: a functional is the integral of a field over its manifold'
    end if
    if (size(this % component) /= 1) then
       error stop 'operation_field: a functional is a scalar field'
    end if
    if (this % weight /= 0.0_dp) then
       error stop 'operation_field: a smooth maximum is not the integral of its field; at evaluates it'
    end if
    integrand = jet_form(this % component(1))
    if (.not. allocated(solution % rule)) then
       error stop 'operation_field: a functional is evaluated on the solution of a residual, which &
            &stores the jet at every point'
    end if
    ! a solution retained over checkpoints: the objective's value, summed
    ! by its residual segment by segment; no other functional
    if (allocated(solution % objective_value)) then
       view = expression_view(integrand)
       if (view % formula() /= solution % objective_formula) then
          error stop 'operation_field: the solution retained its stages over checkpoints alone, so the &
               &only functional it integrates is the objective of the residual that formed it'
       end if
       total = solution % objective_value
       allocate(gradient(0, 0, 0), partial(0, 0))
       return
    end if
    if (.not. (allocated(solution % jet) .or. allocated(solution % node_jet))) then
       error stop 'operation_field: a functional is evaluated on the solution of a residual, which &
            &stores the jet at every point'
    end if
    ! over the quadrature nodes of the solution when it stores them,
    ! the families' own rule, their jets in the columns of the
    ! support's expansion; over the points with their measure
    ! otherwise
    over_nodes = allocated(solution % node_weight)
    if (over_nodes) then
       stride   = size(solution % node_jet, 1)
       nnodes   = size(solution % node_weight)
       expanded = solution % on % expansion
    else
       stride   = size(solution % jet, 1)
       nnodes   = size(solution % jet, 2)
       expanded = jet_expansion(solution)
    end if
    if (this % integrated_part == TIME_WINDOW) within = window_weights(solution, this % integrated_time, &
         & this % integrated_end)
    order  = expanded % order
    values = expanded % values()
    if (size(values) < integrand % num_designs()) then
       error stop 'operation_field: the functional reads a design coordinate the manifold does not have'
    end if
    if (over_nodes .and. allocated(solution % node_rule)) then
       call jet_slots(integrand, solution % node_rule, slot)
    else
       call jet_slots(integrand, solution % rule, slot)
    end if
    ns = size(slot)
    allocate(q(0:ns - 1), g(0:ns - 1), gd(size(values)), jets(0:ns - 1), varied(0:ns - 1))
    allocate(nu(size(values)), along(size(values)), design(size(values)))
    allocate(gradient(stride, nnodes, 0:expanded % num_columns() - 1), source=0.0_dp)
    allocate(total(expanded % num_columns()), partial(size(values), 0:expanded % num_columns() - 1), source=0.0_dp)
    allocate(unit(order), source=0.0_dp)
    if (order > 0) unit(1) = 1.0_dp

    do p = 1, nnodes
       if (over_nodes) then
          if (solution % node_image(p) /= this_image()) cycle
          measure  = solution % node_weight(p)
          position = solution % node_position(:, p)
          if (this % integrated_part == TIME_FACE) measure = face_weight(solution, p, this % integrated_time)
          if (this % integrated_part == TIME_WINDOW) measure = within(p)
       else
          measure  = solution % on % measure(p)
          position = solution % on % position(p)
          if (this % integrated_part == TIME_FACE) measure = solution % on % measure_on_face(p, this % integrated_time)
       end if
       if (measure == 0.0_dp) cycle
       ! the value, its gradient in the jet and in the design
       do k = 0, ns - 1
          q(k) = jet_of(slot(k), p, 0)
       end do
       call integrand % gradient_at(q, values, value, g, position, gd)
       total(1) = total(1) + measure * value
       do k = 0, ns - 1
          gradient(slot(k), p, 0) = gradient(slot(k), p, 0) + measure * g(k)
       end do
       partial(:, 0) = partial(:, 0) + measure * gd
       if (order == 0) cycle
       ! along each expansion d, whose columns begin at first: the jets
       ! over the order's directions, and over one more, the last,
       ! seeded on a component or on a coordinate
       do d = 1, expanded % num_expansions()
          seed  = expanded % seed(d)
          first = expanded % column(d, 1)
          do i = 1, size(values)
             nu(i)    = symmetric_terms(values(i), seed(i) * unit, order)
             along(i) = symmetric_terms(values(i), seed(i) * unit, order + 1)
          end do
          do k = 0, ns - 1
             jets(k)   = symmetric_terms(q(k), [(jet_of(slot(k), p, first + o - 1), o = 1, order)], order)
             varied(k) = symmetric_terms(q(k), [(jet_of(slot(k), p, first + o - 1), o = 1, order)], order + 1)
          end do
          r = integrand % at_instant(jets, nu, position)
          do o = 1, order
             total(first + o) = total(first + o) + measure * coefficient(r, 2**o - 1)
          end do
          do i = 1, size(values)
             design = along
             call design(i) % set_coefficient(2**order, 1.0_dp)
             r = integrand % at_instant(varied, design, position)
             do o = 1, order
                partial(i, first + o - 1) = partial(i, first + o - 1) + measure * coefficient(r, 2**o - 1 + 2**order)
             end do
          end do
          do k = 0, ns - 1
             call varied(k) % set_coefficient(2**order, 1.0_dp)
             r = integrand % at_instant(varied, along, position)
             call varied(k) % set_coefficient(2**order, 0.0_dp)
             do o = 1, order
                gradient(slot(k), p, first + o - 1) = gradient(slot(k), p, first + o - 1) &
                     & + measure * coefficient(r, 2**o - 1 + 2**order)
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

    ! the jet at a node or a point, in the column c
    pure real(dp) function jet_of(row, n, c)
      integer, intent(in) :: row, n, c
      if (over_nodes) then
         jet_of = solution % node_jet(row, n, c + 1)
      else
         jet_of = solution % jet(row, n, c + 1)
      end if
    end function jet_of

  end subroutine functional_over

  !===================================================================!
  ! THE SMOOTH MAXIMUM ON A SOLUTION in every column of the expansion
  ! of order one: m the largest component at the nodes, or at the
  ! points, of positive measure; I the integral of the sum of exp(rho
  ! (f_i - m)), in every column; the value m + log(I) / rho and its
  ! derivative along each expansion I' / (rho I).
  !===================================================================!

  function maximum_over(this, solution) result(total)

    class(continuous_field), intent(in) :: this
    type(discrete_field)   , intent(in) :: solution
    real(dp), allocatable :: total(:)

    real(dp), allocatable :: gradient(:,:,:), partial(:,:)
    real(dp) :: m

    call maximum_parts(this, solution, m, total, gradient, partial)
    total(2:) = total(2:) / (this % weight * total(1))
    total(1)  = m + log(total(1)) / this % weight

  end function maximum_over

  ! the largest component m - the smallest, of a negative weight - and
  ! the integral I of the shifted exponentials with its differential,
  ! in every column
  subroutine maximum_parts(this, solution, m, total, gradient, partial)

    class(continuous_field), intent(in)  :: this
    type(discrete_field)   , intent(in)  :: solution
    real(dp)               , intent(out) :: m
    real(dp), allocatable  , intent(out) :: total(:), gradient(:,:,:), partial(:,:)

    type(continuous_field) :: shifted
    type(expression) :: summed, g
    real(dp), allocatable :: q(:), w(:), values(:), gd(:), position(:), within(:)
    integer , allocatable :: slot(:)
    real(dp) :: value, measure
    integer  :: i, p, k, nnodes
    logical  :: over_nodes

    if (solution % on % expansion % order > 1) then
       error stop 'operation_field: a smooth maximum is evaluated along the design to the first order'
    end if
    if (.not. (allocated(solution % rule) .and. (allocated(solution % jet) .or. allocated(solution % node_jet)))) then
       error stop 'operation_field: a smooth maximum is evaluated on the solution of a residual, which stores &
            &the jet at every point'
    end if
    ! the largest component, at the nodes of the families' rule or at
    ! the points
    over_nodes = allocated(solution % node_weight)
    if (over_nodes) then
       nnodes = size(solution % node_weight)
    else
       nnodes = size(solution % jet, 2)
    end if
    if (this % integrated_part == TIME_WINDOW) within = window_weights(solution, this % integrated_time, &
         & this % integrated_end)
    values = solution % on % expansion % values()
    m = merge(-huge(1.0_dp), huge(1.0_dp), this % weight > 0.0_dp)
    do i = 1, size(this % component)
       g = jet_form(this % component(i))
       if (over_nodes .and. allocated(solution % node_rule)) then
          call jet_slots(g, solution % node_rule, slot)
       else
          call jet_slots(g, solution % rule, slot)
       end if
       if (allocated(q)) deallocate(q, w, gd)
       allocate(q(0:size(slot) - 1), w(0:size(slot) - 1), gd(size(values)))
       do p = 1, nnodes
          if (over_nodes) then
             if (solution % node_image(p) /= this_image()) cycle
             measure  = solution % node_weight(p)
             position = solution % node_position(:, p)
             if (this % integrated_part == TIME_FACE) measure = face_weight(solution, p, this % integrated_time)
             if (this % integrated_part == TIME_WINDOW) measure = within(p)
             do k = 0, size(slot) - 1
                q(k) = solution % node_jet(slot(k), p, 1)
             end do
          else
             measure  = solution % on % measure(p)
             position = solution % on % position(p)
             if (this % integrated_part == TIME_FACE) measure = solution % on % measure_on_face(p, this % integrated_time)
             do k = 0, size(slot) - 1
                q(k) = solution % jet(slot(k), p, 1)
             end do
          end if
          if (measure == 0.0_dp) cycle
          call g % gradient_at(q, values, value, w, position, gd)
          if (this % weight > 0.0_dp) then
             m = max(m, value)
          else
             m = min(m, value)
          end if
       end do
    end do
    if (over_nodes) then
       if (any(solution % node_image /= this_image())) then
          if (this % weight > 0.0_dp) then
             call co_max(m)
          else
             call co_min(m)
          end if
       end if
    end if

    ! the integral of the shifted exponentials, and the logarithm
    summed = exp(this % weight * (this % component(1) - m))
    do i = 2, size(this % component)
       summed = summed + exp(this % weight * (this % component(i) - m))
    end do
    shifted = this
    shifted % weight    = 0.0_dp
    shifted % component = [summed]
    shifted % index     = [0]
    call functional_over(shifted, solution, total, gradient, partial)

  end subroutine maximum_parts

  !===================================================================!
  ! THE EXPANSION A FIELD'S JET IS ALONG: the support's, to the order
  ! the jet stores; for the derivative along one expansion, that
  ! expansion alone, from the derivative's own column on.
  !===================================================================!

  function jet_expansion(this) result(expanded)

    type(discrete_field), intent(in) :: this
    type(expansion) :: expanded

    real(dp), allocatable :: seed(:)

    expanded = this % on % expansion
    if (this % along /= 0) then
       seed = expanded % seed(this % along)
       expanded % along = reshape(seed, [size(seed), 1])
    end if
    expanded % order = (size(this % jet, 3) - 1) / max(1, expanded % num_expansions())

  end function jet_expansion

  ! the quadrature nodes of a solution, copied to a field read from it
  subroutine nodes_copied(from, to)
    type(discrete_field), intent(in)    :: from
    type(discrete_field), intent(inout) :: to
    if (allocated(from % node_weight)) then
       to % node_jet      = from % node_jet
       if (allocated(from % node_rule)) to % node_rule = from % node_rule
       to % node_weight   = from % node_weight
       to % node_position = from % node_position
       to % node_face_weight = from % node_face_weight
       if (allocated(from % share_node)) then
          to % share_node   = from % share_node
          to % share_end    = from % share_end
          to % share_weight = from % share_weight
       end if
       to % node_image    = from % node_image
    end if
    if (allocated(from % objective_value)) then
       to % objective_value   = from % objective_value
       to % objective_formula = from % objective_formula
    end if
  end subroutine nodes_copied

  ! the weight of node p within the face at the instant given: its
  ! weight within the face of its own instant there, zero elsewhere
  pure real(dp) function face_weight(solution, p, time)
    type(discrete_field), intent(in) :: solution
    integer             , intent(in) :: p
    real(dp)            , intent(in) :: time
    face_weight = 0.0_dp
    if (abs(solution % node_position(1, p) - time) <= 1.0e-12_dp * max(1.0_dp, abs(time))) then
       face_weight = solution % node_face_weight(p)
    end if
  end function face_weight

  ! the weights of the nodes within the window from the instant a to
  ! the instant b: of each node, its shares of the steps that end
  ! after a and at b at the latest, the steps within the window
  function window_weights(solution, a, b) result(within)
    type(discrete_field), intent(in) :: solution
    real(dp)            , intent(in) :: a, b
    real(dp), allocatable :: within(:)
    real(dp) :: tolerance
    integer  :: s
    if (.not. allocated(solution % share_node)) then
       error stop 'operation_field: a window is measured by the steps of the march, on the nodes of a solution &
            &that stores their shares'
    end if
    tolerance = 1.0e-12_dp * max(1.0_dp, abs(a), abs(b))
    allocate(within(size(solution % node_weight)), source=0.0_dp)
    do s = 1, size(solution % share_node)
       if (solution % share_end(s) > a + tolerance .and. solution % share_end(s) <= b + tolerance) then
          within(solution % share_node(s)) = within(solution % share_node(s)) + solution % share_weight(s)
       end if
    end do
  end function window_weights

  !===================================================================!
  ! THE DERIVATIVE ALONG TIME AND SPACE of an unknown read from a
  ! solution: the component its jet stores of the orders the
  ! multi-index names along the first coordinate and along at most
  ! one later coordinate, with its columns along the design. Nothing
  ! is computed. Invalid input: a field of several components, or not
  ! read from a solution, or a multiplier; a solution under a
  ! checkpoint or snapshot schedule, whose jet keeps the values alone;
  ! two later coordinates; an order the jet does not store.
  !===================================================================!

  function stored_component(this, along) result(d)

    class(discrete_field)  , intent(in) :: this
    type(continuous_field) , intent(in) :: along(:)
    type(discrete_field) :: d

    integer :: k, f, field, first, later, coordinate, row
    character(len=32) :: name

    if (size(this % value, 2) /= 1) then
       error stop 'operation_field: a derivative along time or space is read of one unknown'
    end if
    if (.not. (allocated(this % jet) .and. allocated(this % slot) .and. allocated(this % rule))) then
       error stop 'operation_field: a derivative along time or space is read from the jet of a solution'
    end if
    if (this % along /= 0) then
       error stop 'operation_field: a derivative along time or space is read before those along the design'
    end if
    if (allocated(this % objective_value)) then
       error stop 'operation_field: under a checkpoint or snapshot schedule the jet of a solution keeps the &
            &values alone'
    end if
    field = 0
    ! the unknowns of the state, the multipliers the rule reads after
    ! them having no row of the tuple
    do f = 1, this % rule % num_fields() - this % rule % num_multipliers()
       if (this % rule % offset_of_field(f) + 1 == this % slot(1)) field = f
    end do
    if (field == 0) then
       error stop 'operation_field: a derivative along time or space is read of an unknown, whose jet the &
            &solution stores'
    end if
    first      = 0
    later      = 0
    coordinate = 0
    name       = trim(this % name(1)) // '_'
    do k = 1, size(along)
       name = trim(name) // trim(along(k) % name)
       if (along(k) % coordinate == FIRST_COORDINATE) then
          first = first + 1
          cycle
       end if
       if (coordinate /= 0 .and. along(k) % coordinate /= coordinate) then
          error stop 'operation_field: the jet stores the derivatives along one later coordinate at a time'
       end if
       coordinate = along(k) % coordinate
       later      = later + 1
    end do
    if (later == 0) then
       if (first > this % rule % degree_of_field(field)) then
          error stop 'operation_field: the jet does not store this derivative along the first coordinate'
       end if
       row = this % rule % component_at(FIRST_COORDINATE, first, field) + 1
    else
       if (later > this % rule % degree_along(coordinate) .or. first > this % rule % mixed_degree(field)) then
          error stop 'operation_field: the jet does not store this derivative'
       end if
       row = this % rule % component_at(coordinate, later, field, first) + 1
    end if

    allocate(d % on, source=this % on)
    allocate(d % name(1))
    d % name(1) = name
    d % value   = reshape(this % jet(row, :, 1), [size(this % jet, 2), 1])
    d % jet     = this % jet
    d % slot    = [row]
    d % rule    = this % rule
    d % along   = 0
    call nodes_copied(this, d)

  end function stored_component

  !===================================================================!
  ! THE DERIVATIVE OF A DISCRETE FIELD. Along time and space, the
  ! component the jet of a solution stores, read and not computed;
  ! the entries of the multi-index along time and space come before
  ! those along the design, which are then read from that component's
  ! columns, as below.
  !
  ! Along the design: the
  ! multi-index names design coordinates or directions, and the
  ! derivative is read from the columns stored - along the expansion
  ! whose seed is the coordinate's unit vector or the direction's
  ! vector, for a multi-index repeating one entry: the field whose
  ! values are the column of that order, with the columns of that
  ! expansion after it and those of the other expansions dropped. A mixed
  ! derivative along two coordinates, [nu, mu], is the polarization
  ! (u_ww - u_nunu - u_mumu) / 2 of the second derivatives along nu,
  ! mu and w = nu + mu, and requires the expansions along all three;
  ! it stores no jet of its own. Invalid input: an entry that is not
  ! a design coordinate or direction, an expansion not stored, a
  ! coordinate other than the one the field is already a derivative
  ! along, an order past the expansion's, a mixed derivative of an
  ! order other than two, or a component with no jet.
  !===================================================================!

  recursive function discrete_derivative(this, along) result(d)

    class(discrete_field)  , intent(in) :: this
    type(continuous_field) , intent(in) :: along(:)
    type(discrete_field) :: d

    type(expansion) :: expanded
    real(dp), allocatable :: seed(:), other(:), both(:)
    integer :: n, k, order, npts, e, second, third, c, m
    logical :: mixed
    character(len=250) :: message

    ! along time and space: the stored component, then its columns
    ! along the design
    m = count(.not. along % design_coordinate)
    if (m > 0) then
       if (any(along(1:m) % design_coordinate)) then
          error stop 'operation_field: a multi-index names the coordinates of time and space before those of &
               &the design'
       end if
       d = stored_component(this, along(1:m))
       if (m < size(along)) d = d % derivative(along(m + 1:))
       return
    end if

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
    expanded = jet_expansion(this)
    order = expanded % order
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
       e      = expansion_of(seed)
       second = expansion_of(other)
       third  = expansion_of(both)
       call require_columns(e, second, third)
       allocate(d % on, source=this % on)
       d % name  = this % name
       d % slot  = this % slot
       d % along = 0
       allocate(d % value(npts, size(this % value, 2)))
       do k = 1, size(this % value, 2)
          call require_jet(k)
          d % value(:, k) = (this % jet(this % slot(k), :, expanded % column(third, 2) + 1) &
               & - this % jet(this % slot(k), :, expanded % column(e, 2) + 1) &
               & - this % jet(this % slot(k), :, expanded % column(second, 2) + 1)) / 2.0_dp
       end do
       allocate(d % jet(size(this % jet, 1), npts, 1))
       d % jet(:, :, 1) = (this % jet(:, :, expanded % column(third, 2) + 1) - this % jet(:, :, expanded % column(e, 2) + 1) &
            & - this % jet(:, :, expanded % column(second, 2) + 1)) / 2.0_dp
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
       if (.not. same_seed(seed, this % on % expansion % seed(this % along))) then
          error stop 'operation_field: the field is a derivative along one expansion, and its jet along &
               &another is not stored'
       end if
       e = 1
    else
       e = expansion_of(seed)
       if (e == 0) then
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
    d % along = merge(this % along, e, this % along /= 0)
    ! the column of the derivative, and the columns of its expansion
    ! after it
    c = expanded % column(e, n)
    allocate(d % value(npts, size(this % value, 2)))
    do k = 1, size(this % value, 2)
       call require_jet(k)
       d % value(:, k) = this % jet(this % slot(k), :, c + 1)
    end do
    allocate(d % jet(size(this % jet, 1), npts, order - n + 1))
    d % jet = this % jet(:, :, c + 1:c + 1 + order - n)
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
      nd = this % on % expansion % num_designs()
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
      do e = 1, this % on % expansion % num_expansions()
         if (same_seed(s, this % on % expansion % seed(e))) expansion_of = e
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

    integer :: f, d, c, o, a, fields
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
          if (g % mixed_degree(f) > rule % mixed_degree(f)) then
             write(message,'(a,i0,a,i0)') 'operation_field: the functional reads a mixed derivative of order ', &
                  & g % mixed_degree(f), ' along the first coordinate, and the rule stores the order ', &
                  & rule % mixed_degree(f)
             error stop trim(message)
          end if
          do a = 0, g % mixed_degree(f)
             do o = 1, g % degree_along(c)
                slot(g % component_at(c, o, f, a)) = rule % component_at(c, o, f, a) + 1
             end do
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
  ! entry one order along that coordinate, an operator stored in the law
  ! until jet_form replaces it by the jet components it denotes. An
  ! unknown along a coordinate it is not a function of has the zero
  ! derivative. Invalid input: a field of several components, an empty
  ! multi-index, or an entry that is not a time or space coordinate of
  ! the field's manifold.
  !===================================================================!

  function field_derivative(this, along) result(d)

    class(continuous_field), intent(in) :: this
    type(continuous_field) , intent(in) :: along(:)
    type(continuous_field) :: d

    integer :: k

    if (size(this % component) /= 1) then
       error stop 'operation_field: derivative requires a scalar field'
    end if
    if (size(along) < 1) then
       error stop 'operation_field: derivative requires one coordinate at least'
    end if
    do k = 1, size(along)
       if (along(k) % design_coordinate .or. along(k) % coordinate < 1) then
          error stop 'operation_field: the derivative of a field is along a coordinate of time or space'
       end if
       if (.not. along(k) % on % matches(this % on)) then
          error stop 'operation_field: the derivative is along a coordinate of the field''s manifold'
       end if
    end do
    ! an unknown along a coordinate it is not a function of: zero
    if (this % index(1) >= 1 .and. allocated(this % argument)) then
       if (any([(.not. any(this % argument == along(k) % name), k = 1, size(along))])) then
          d = this % jet_component(along)
          return
       end if
    end if
    ! the operator, one order per entry of the multi-index, replaced by
    ! jet components when the law is discretized
    d = this
    d % index = [0]
    do k = 1, size(along)
       d % component = [derivative(d % component(1), along(k) % coordinate, 1)]
    end do
    ! of a field with an integral over the region, along time: the
    ! derivative of its terms and the integral of the integrand's
    if (allocated(this % integrand)) then
       if (any([(along(k) % name /= 't', k = 1, size(along))])) then
          error stop 'operation_field: a field with an integral over the region is differentiated along time alone'
       end if
       do k = 1, size(along)
          d % integrand = [derivative(d % integrand(1), along(k) % coordinate, 1)]
       end do
    end if

  end function field_derivative

  !===================================================================!
  ! The field with every derivative operator of its components replaced
  ! by the jet components it denotes.
  !===================================================================!

  function field_jet_form(this) result(form)

    class(continuous_field), intent(in) :: this
    type(continuous_field) :: form

    integer :: k

    form = this
    do k = 1, size(this % component)
       form % component(k) = jet_form(this % component(k))
    end do
    if (allocated(this % integrand)) then
       do k = 1, size(this % integrand)
          form % integrand(k) = jet_form(this % integrand(k))
       end do
    end if

  end function field_jet_form

  !===================================================================!
  ! THE JET COMPONENT of an unknown along a multi-index of
  ! coordinate functions: each entry one order along that coordinate.
  ! Along a coordinate the unknown is not a function of, one outside
  ! its arguments, the component is the zero field. A multi-index of
  ! the first coordinate and one later coordinate names a mixed
  ! component, the entries commuting. Invalid input: a field that is
  ! not a scalar unknown, an entry that is not a coordinate of the
  ! same manifold, or two distinct later coordinates, whose mixed
  ! components the jet does not store.
  !===================================================================!

  function field_jet_component(this, along) result(d)

    class(continuous_field), intent(in) :: this
    type(continuous_field) , intent(in) :: along(:)
    type(continuous_field) :: d

    integer :: k, c, mixed_order
    character(len=250) :: message

    if (size(this % component) /= 1 .or. this % index(1) < 1) then
       error stop 'operation_field: jet_component requires a scalar unknown'
    end if
    if (size(along) < 1) then
       error stop 'operation_field: jet_component requires one coordinate at least'
    end if
    c = maxval([(along(k) % coordinate, k = 1, size(along))])
    mixed_order = count([(along(k) % coordinate == FIRST_COORDINATE, k = 1, size(along))])
    if (c == FIRST_COORDINATE) mixed_order = 0
    do k = 1, size(along)
       if (along(k) % design_coordinate) then
          error stop 'operation_field: the derivative of a continuous field along the design coordinate &
               &is not stated; the derivatives along the design are read from the discrete solution, &
               &component([nu]) on a discrete field'
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
       if (along(k) % coordinate /= c .and. along(k) % coordinate /= FIRST_COORDINATE) then
          error stop 'operation_field: a mixed derivative of two later coordinates is not a component of &
               &the jet; the multi-index names the first coordinate and one later coordinate'
       end if
    end do

    d % on    = this % on
    d % name  = this % name
    d % index = [0]
    if (allocated(this % argument)) then
       if (any([(.not. any(this % argument == along(k) % name), k = 1, size(along))])) then
          d % component = [constant(0.0_dp)]
          return
       end if
    end if
    if (c == FIRST_COORDINATE) then
       d % component = [component(this % component(1), size(along))]
    else
       d % component = [component_along(this % component(1), c, size(along) - mixed_order, mixed_order)]
    end if

  end function field_jet_component

  !===================================================================!
  ! THE REPRESENTATION LAWS of an unknown along a coordinate, to an
  ! order: for k = 1 to the order, the derivative of order k of u less
  ! the component of its jet the derivative denotes, d^k u - u_(k).
  ! The concretion along the coordinate states each as its rows; stated
  ! as a residual, each is paired with a multiplier. Invalid input: an
  ! order below one.
  !===================================================================!

  function representation(u, along, order) result(laws)

    type(continuous_field), intent(in) :: u, along
    integer               , intent(in) :: order
    type(continuous_field), allocatable :: laws(:)

    integer :: k, j

    if (order < 1) error stop 'operation_field: a representation is of order one or more'
    allocate(laws(order))
    do k = 1, order
       laws(k) = u % derivative([(along, j = 1, k)]) - u % jet_component([(along, j = 1, k)])
    end do

  end function representation

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
       d % component = [partial_in_design(jet_form(this % component(1)), in % coordinate)]
       return
    end if
    ! a derivative of an unknown is the component of its jet it denotes
    leaf = jet_form(in % component(1))
    if (leaf % num_vertices() /= 1 .or. .not. leaf % reads_state()) then
       error stop 'operation_field: the partial derivative is in an unknown or a derivative of one, &
            &u or u % derivative([t]), or in a design coordinate'
    end if
    if (leaf % mixed_degree(leaf % leaf_field()) > 0) then
       error stop 'operation_field: the partial derivative in a mixed component is not taken'
    end if
    d % component = [partial_derivative(jet_form(this % component(1)), leaf % leaf_field(), leaf % leaf_order(), &
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
    d % component = [stationarity(jet_form(this % component(1)), in % index(1), along % coordinate)]

  end function field_stationarity

  !===================================================================!
  ! ALONG SEVERAL COORDINATES: the partial in the value, and for each
  ! coordinate and each order o the density reads along it, (-1)^o
  ! times the o-th derivative along it of the partial in the component
  ! of that order, each derivative an operator of the field, replaced
  ! when the law is discretized - so that a projection on the cells
  ! reads the outer derivatives along space as fluxes. Invalid input
  ! as for one coordinate.
  !===================================================================!

  function field_stationarity_along(this, in, along) result(d)

    class(continuous_field), intent(in) :: this
    type(continuous_field) , intent(in) :: in, along(:)
    type(continuous_field) :: d

    type(continuous_field) :: term
    type(expression) :: form
    integer :: k, o, top

    if (size(this % component) /= 1) then
       error stop 'operation_field: the stationarity is of a scalar field'
    end if
    if (size(in % component) /= 1 .or. in % index(1) < 1) then
       error stop 'operation_field: the stationarity is in a scalar unknown of the manifold'
    end if
    form = jet_form(this % component(1))
    d = this % partial(in)
    do k = 1, size(along)
       if (along(k) % design_coordinate .or. along(k) % coordinate < 1) then
          error stop 'operation_field: the stationarity integrates by parts along a coordinate of time or space'
       end if
       top = form % highest_degree_along(along(k) % coordinate, in % index(1))
       do o = 1, top
          term = this % partial(in % derivative(spread(along(k), 1, o)))
          term = term % derivative(spread(along(k), 1, o))
          if (mod(o, 2) == 0) then
             d = d + term
          else
             d = d - term
          end if
       end do
    end do

  end function field_stationarity_along

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
    graphs = boundary_terms(jet_form(this % component(1)), in % index(1), along % coordinate)
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
  ! name. The jet is retained. A design coordinate of the support
  ! that no component is named after is read as the component whose
  ! value at every point is the design the support is at, and whose
  ! jet, a row after the field's, is that of the coordinate: its
  ! value, the entry of each seed at order one, zero above. A name
  ! present as neither stops the program.
  !===================================================================!

  function fields(this, names) result(chosen)

    class(discrete_field), intent(in) :: this
    character(len=*)     , intent(in) :: names(:)
    type(discrete_field) :: chosen

    integer, allocatable :: slot(:), design(:)
    integer :: k, j, count, nrow, d

    allocate(chosen % on, source=this % on)
    count = 0
    do k = 1, size(names)
       if (count_named(names(k)) == 0 .and. design_index(names(k)) == 0) then
          error stop 'operation_field: fields names a component the field does not have: ' // trim(names(k))
       end if
       count = count + max(count_named(names(k)), 1)
    end do
    allocate(chosen % name(count), chosen % value(size(this % value, 1), count), slot(count))
    allocate(design(0))
    nrow = 0
    if (allocated(this % jet)) nrow = size(this % jet, 1)
    count = 0
    do k = 1, size(names)
       if (count_named(names(k)) == 0) then
          if (this % along /= 0) then
             error stop 'operation_field: the design is read from the solution, not from a derivative of it'
          end if
          count  = count + 1
          design = [design, design_index(names(k))]
          chosen % name(count)     = names(k)
          chosen % value(:, count) = this % on % expansion % at(design(size(design)))
          slot(count) = nrow + size(design)
          cycle
       end if
       do j = 1, size(this % name)
          if (trim(this % name(j)) == trim(names(k))) then
             count = count + 1
             chosen % name(count)     = this % name(j)
             chosen % value(:, count) = this % value(:, j)
             if (allocated(this % slot)) slot(count) = this % slot(j)
          end if
       end do
    end do
    if (allocated(this % jet)) then
       allocate(chosen % jet(nrow + size(design), size(this % jet, 2), size(this % jet, 3)), source=0.0_dp)
       chosen % jet(1:nrow, :, :) = this % jet
       associate (expanded => this % on % expansion)
         do k = 1, size(design)
            chosen % jet(nrow + k, :, 1) = expanded % at(design(k))
            if (expanded % order < 1) cycle
            do d = 1, expanded % num_expansions()
               chosen % jet(nrow + k, :, 1 + expanded % column(d, 1)) = expanded % along(design(k), d)
            end do
         end do
       end associate
    end if
    if (allocated(this % rule)) chosen % rule = this % rule
    chosen % along = this % along
    call nodes_copied(this, chosen)
    if (allocated(this % slot)) chosen % slot = slot

  contains

    pure integer function count_named(name)
      character(len=*), intent(in) :: name
      integer :: j
      count_named = 0
      do j = 1, size(this % name)
         if (trim(this % name(j)) == trim(name)) count_named = count_named + 1
      end do
    end function count_named

    ! the design coordinate of the support named so, zero for none
    pure integer function design_index(name)
      character(len=*), intent(in) :: name
      integer :: i
      design_index = 0
      if (.not. allocated(this % on % expansion % name)) return
      do i = 1, size(this % on % expansion % name)
         if (trim(this % on % expansion % name(i)) == trim(name)) design_index = i
      end do
    end function design_index

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

    class(discrete_field), intent(in) :: a
    type(discrete_field) , intent(in) :: b
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

  !===================================================================!
  ! THE TWO PARTS OF A SCALAR FIELD: the field itself and no integral,
  ! or, for the integral over the space factor of a manifold of time
  ! and space or a field with such an integral, its terms of time and
  ! the integrand over the region. A field with an integral is
  ! composed of both parts.
  !===================================================================!

  logical function with_integral(a)
    type(continuous_field), intent(in) :: a
    with_integral = allocated(a % integrand) .or. (a % is_integral() .and. a % integrated_part == SPACE_FACTOR)
  end function with_integral

  subroutine parts(a, terms, integrand)
    type(continuous_field), intent(in)  :: a
    type(expression)      , intent(out) :: terms, integrand
    if (allocated(a % integrand)) then
       terms     = a % component(1)
       integrand = a % integrand(1)
    else if (with_integral(a)) then
       terms     = constant(0.0_dp)
       integrand = a % component(1)
    else
       terms     = a % component(1)
       integrand = constant(0.0_dp)
    end if
  end subroutine parts

  function composed_parts(a, terms, integrand) result(this)
    type(continuous_field), intent(in) :: a
    type(expression)      , intent(in) :: terms, integrand
    type(continuous_field) :: this
    this = composed(a, terms)
    this % integrand = [integrand]
  end function composed_parts

  subroutine require_linear(a, symbol)
    type(continuous_field), intent(in) :: a
    character(len=*)      , intent(in) :: symbol
    if (with_integral(a)) then
       error stop 'operation_field: an integral over the region enters a field linearly, not through ' // symbol
    end if
  end subroutine require_linear

  impure elemental function plus(a, b) result(this)
    class(continuous_field), intent(in) :: a
    type(continuous_field) , intent(in) :: b
    type(continuous_field) :: this
    type(expression) :: ta, ia, tb, ib
    call require_scalar(a, '+'); call require_scalar(b, '+'); call require_one_manifold(a, b, '+')
    if (with_integral(a) .or. with_integral(b)) then
       call parts(a, ta, ia); call parts(b, tb, ib)
       this = composed_parts(a, ta + tb, ia + ib)
       return
    end if
    this = composed(a, a % component(1) + b % component(1))
  end function plus

  impure elemental function real_plus(x, b) result(this)
    real(dp)               , intent(in) :: x
    class(continuous_field), intent(in) :: b
    type(continuous_field) :: this
    type(expression) :: tb, ib
    call require_scalar(b, '+')
    if (with_integral(b)) then
       call parts(b, tb, ib)
       this = composed_parts(b, x + tb, ib)
       return
    end if
    this = composed(b, x + b % component(1))
  end function real_plus

  impure elemental function plus_real(a, x) result(this)
    class(continuous_field), intent(in) :: a
    real(dp)               , intent(in) :: x
    type(continuous_field) :: this
    call require_scalar(a, '+')
    this = real_plus(x, a)
    if (.not. with_integral(a)) this = composed(a, a % component(1) + x)
  end function plus_real

  impure elemental function minus(a, b) result(this)
    class(continuous_field), intent(in) :: a
    type(continuous_field) , intent(in) :: b
    type(continuous_field) :: this
    type(expression) :: ta, ia, tb, ib
    call require_scalar(a, '-'); call require_scalar(b, '-'); call require_one_manifold(a, b, '-')
    if (with_integral(a) .or. with_integral(b)) then
       call parts(a, ta, ia); call parts(b, tb, ib)
       this = composed_parts(a, ta - tb, ia - ib)
       return
    end if
    this = composed(a, a % component(1) - b % component(1))
  end function minus

  impure elemental function real_minus(x, b) result(this)
    real(dp)               , intent(in) :: x
    class(continuous_field), intent(in) :: b
    type(continuous_field) :: this
    type(expression) :: tb, ib
    call require_scalar(b, '-')
    if (with_integral(b)) then
       call parts(b, tb, ib)
       this = composed_parts(b, x - tb, -ib)
       return
    end if
    this = composed(b, x - b % component(1))
  end function real_minus

  impure elemental function minus_real(a, x) result(this)
    class(continuous_field), intent(in) :: a
    real(dp)               , intent(in) :: x
    type(continuous_field) :: this
    type(expression) :: ta, ia
    call require_scalar(a, '-')
    if (with_integral(a)) then
       call parts(a, ta, ia)
       this = composed_parts(a, ta - x, ia)
       return
    end if
    this = composed(a, a % component(1) - x)
  end function minus_real

  impure elemental function negated(a) result(this)
    class(continuous_field), intent(in) :: a
    type(continuous_field) :: this
    type(expression) :: ta, ia
    call require_scalar(a, '-')
    if (with_integral(a)) then
       call parts(a, ta, ia)
       this = composed_parts(a, -ta, -ia)
       return
    end if
    this = composed(a, -a % component(1))
  end function negated

  ! a field with an integral times a field of time alone: both parts
  ! times it
  impure elemental function times(a, b) result(this)
    class(continuous_field), intent(in) :: a
    type(continuous_field) , intent(in) :: b
    type(continuous_field) :: this
    type(expression) :: ta, ia
    call require_scalar(a, '*'); call require_scalar(b, '*'); call require_one_manifold(a, b, '*')
    if (with_integral(a) .and. with_integral(b)) call require_linear(a, '*')
    if (with_integral(a)) then
       call parts(a, ta, ia)
       this = composed_parts(a, ta * b % component(1), ia * b % component(1))
       return
    else if (with_integral(b)) then
       call parts(b, ta, ia)
       this = composed_parts(b, a % component(1) * ta, a % component(1) * ia)
       return
    end if
    this = composed(a, a % component(1) * b % component(1))
  end function times

  impure elemental function real_times(x, b) result(this)
    real(dp)               , intent(in) :: x
    class(continuous_field), intent(in) :: b
    type(continuous_field) :: this
    type(expression) :: tb, ib
    call require_scalar(b, '*')
    if (with_integral(b)) then
       call parts(b, tb, ib)
       this = composed_parts(b, x * tb, x * ib)
       return
    end if
    this = composed(b, x * b % component(1))
  end function real_times

  impure elemental function times_real(a, x) result(this)
    class(continuous_field), intent(in) :: a
    real(dp)               , intent(in) :: x
    type(continuous_field) :: this
    type(expression) :: ta, ia
    call require_scalar(a, '*')
    if (with_integral(a)) then
       call parts(a, ta, ia)
       this = composed_parts(a, ta * x, ia * x)
       return
    end if
    this = composed(a, a % component(1) * x)
  end function times_real

  impure elemental function over(a, b) result(this)
    class(continuous_field), intent(in) :: a
    type(continuous_field) , intent(in) :: b
    type(continuous_field) :: this
    type(expression) :: ta, ia
    call require_scalar(a, '/'); call require_scalar(b, '/'); call require_one_manifold(a, b, '/')
    call require_linear(b, '/')
    if (with_integral(a)) then
       call parts(a, ta, ia)
       this = composed_parts(a, ta / b % component(1), ia / b % component(1))
       return
    end if
    this = composed(a, a % component(1) / b % component(1))
  end function over

  impure elemental function real_over(x, b) result(this)
    real(dp)               , intent(in) :: x
    class(continuous_field), intent(in) :: b
    type(continuous_field) :: this
    call require_scalar(b, '/'); call require_linear(b, '/')
    this = composed(b, x / b % component(1))
  end function real_over

  impure elemental function over_real(a, x) result(this)
    class(continuous_field), intent(in) :: a
    real(dp)               , intent(in) :: x
    type(continuous_field) :: this
    type(expression) :: ta, ia
    call require_scalar(a, '/')
    if (with_integral(a)) then
       call parts(a, ta, ia)
       this = composed_parts(a, ta / x, ia / x)
       return
    end if
    this = composed(a, a % component(1) / x)
  end function over_real

  impure elemental function to_integer(a, n) result(this)
    class(continuous_field), intent(in) :: a
    integer                , intent(in) :: n
    type(continuous_field) :: this
    call require_scalar(a, '**'); call require_linear(a, '**')
    this = composed(a, a % component(1) ** n)
  end function to_integer

  impure elemental function to_real(a, p) result(this)
    class(continuous_field), intent(in) :: a
    real(dp)               , intent(in) :: p
    type(continuous_field) :: this
    call require_scalar(a, '**'); call require_linear(a, '**')
    this = composed(a, a % component(1) ** p)
  end function to_real

  function sine_of(a) result(this)
    type(continuous_field), intent(in) :: a
    type(continuous_field) :: this
    call require_scalar(a, 'sin'); call require_linear(a, 'sin')
    this = composed(a, sin(a % component(1)))
  end function sine_of

  function cosine_of(a) result(this)
    type(continuous_field), intent(in) :: a
    type(continuous_field) :: this
    call require_scalar(a, 'cos'); call require_linear(a, 'cos')
    this = composed(a, cos(a % component(1)))
  end function cosine_of

  function exponential_of(a) result(this)
    type(continuous_field), intent(in) :: a
    type(continuous_field) :: this
    call require_scalar(a, 'exp'); call require_linear(a, 'exp')
    this = composed(a, exp(a % component(1)))
  end function exponential_of

  function logarithm_of(a) result(this)
    type(continuous_field), intent(in) :: a
    type(continuous_field) :: this
    call require_scalar(a, 'log'); call require_linear(a, 'log')
    this = composed(a, log(a % component(1)))
  end function logarithm_of

  function root_of(a) result(this)
    type(continuous_field), intent(in) :: a
    type(continuous_field) :: this
    call require_scalar(a, 'sqrt'); call require_linear(a, 'sqrt')
    this = composed(a, sqrt(a % component(1)))
  end function root_of

  ! the step H(a), of the coordinates alone
  function step_of(a) result(this)
    type(continuous_field), intent(in) :: a
    type(continuous_field) :: this
    call require_scalar(a, 'heaviside'); call require_linear(a, 'heaviside')
    this = composed(a, heaviside(a % component(1)))
  end function step_of

end module operation_field
