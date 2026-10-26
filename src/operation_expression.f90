!=====================================================================!
! Module that states a continuous law in named coordinates, stored as
! data: a graph whose vertices are typed operators and whose edges are
! the reads between them.
!
! Domain. The law itself has no point graph. Its domain is read from
! its leaves: the fields, the highest derivative order of each field
! along each coordinate, and the multipliers. A discrete domain places
! the law on a directed graph when values are evaluated.
!
!      input graph       the points where the law is evaluated
!      state argument    the state, one field over those points with
!                        num_components() components per point
!      design argument   the design, num_designs() values per point,
!                        so a design that varies over the points needs
!                        no other shape
!      output            one value per point
!
! Vertices. A leaf reads a component of an argument (the state's
! component of degree d, or the design), or a multiplier, or a
! coordinate of the point's position, or stores a constant. Every
! other vertex is an arithmetic operation, a power, or an elementary
! function of the vertices it reads. The vertices are stored in
! evaluation order, every read before its reader, and the root is the
! last one.
!
! A rule is built by the intrinsic operators on expressions, so it
! reads as it is written and the compiler checks it:
!
!      q  = unknown()
!      nu = design()
!      r  = component(q, 2) + nu * component(q, 1) + sin(component(q, 0))
!
! component(q, d) is the component of degree d: along the instants
! the derivatives are unknowns the scheme relates, so the vertex
! selects and does not differentiate. The state r reads stores the
! components of degree 0, 1 and 2: the highest degree any leaf names.
! The rule component(q, 2) + nu * component(q, 1) is stored as
!
!      vertex   kind        reads
!      1        leaf q_2
!      2        leaf nu
!      3        leaf q_1
!      4        product     2, 3
!      5        sum         1, 4     the root, last
!
! Several fields. unknown(i) is the i-th zeroth-order field; unknown()
! is the first. A rule that reads fields 1..m is over m fields, the
! highest index read. Along the first coordinate each field stores its
! own jet, to the highest degree the rule names for it; along every
! later coordinate the fields share the highest degree named there. A
! point's tuple lists the fields in order, each field's components
! together.
!
!      tuple   [ u  u_t  u_tt | u_x  u_xx | u_tx  u_txx ] [ v ... ]
!                along t        along x     mixed           field 2
!
! The Lagrangian. A rule that reads multiplier(j), j = 1..k, is a
! Lagrangian, its j-th multiplier paired with its j-th field, and the
! rules it generates are its partials. euler_lagrange(l, j) is the
! stationarity of l in its j-th multiplier, the coefficient of one
! direction seeded on that multiplier, and at_zero(l) is l with every
! multiplier at zero. A multiplier is not part of the tuple a rule
! reads, so the state contract excludes it; it is supplied inside the
! evaluation, zero in value, with its own direction when its
! stationarity is read. Nothing is rewritten: the same vertices are
! evaluated with one more direction, but for a Lagrangian that is a
! sum of terms, each a multiplier times a law free of multipliers,
! whose stationarity in a multiplier is the vertices of the law it
! multiplies, the same numbers without the vertices of the other
! laws.
!
! The partials. The rule is evaluated over derivative_terms by one
! loop over the vertices, so the value and every mixed partial in the
! directions requested come from one pass and are exact; nothing is
! differentiated symbolically. A variation may name either argument,
! so the partials in the state, in the design and the mixed partials
! all come from one path and to any degree. A Newton block requires
! the partials in the state as a row of numbers, which is one call per
! degree with a direction that is one at that degree and zero
! elsewhere.
!
! What is rejected. A derivative of anything but the unknown, or of
! negative degree; a state whose extent is not the points times the
! components; a design that is not num_designs() values per point; a
! missing argument; a variation naming neither argument; a function
! index outside those defined; a stationarity in a multiplier the rule
! does not read. Each stops the program.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_expression

  use util_precision       , only : dp
  use operation_action     , only : operation, variation, contract
  use operation_action     , only : binding, seeded_argument, emit_real, bound_real_vector
  use field_calculus       , only : FIELD_REAL
  use view_directed        , only : directed_graph
  use field_calculus       , only : field
  use util_derivative_terms, only : derivative_terms, mixed_partial, max_subset_width
  use util_derivative_terms, only : partial, extended_terms, product_terms, quotient_terms
  use util_derivative_terms, only : integer_power_terms, real_power_terms, sine_terms, cosine_terms
  use util_derivative_terms, only : exponential_terms, logarithm_terms

  implicit none

  private
  public :: expression
  public :: unknown, design, constant, multiplier, coordinate, time_map, component, component_along
  public :: derivative, jet_form, represented_component, on_domain, paired
  public :: euler_lagrange, at_zero
  public :: partial_derivative, partial_in_design, total_derivative, stationarity, boundary_terms
  public :: flux_form, coefficient_form
  public :: operator(+), operator(-), operator(*), operator(/), operator(**)
  public :: sin, cos, exp, log, sqrt, heaviside

  ! The vertex kinds
  integer, parameter, public :: VERTEX_LEAF          = 1
  integer, parameter, public :: VERTEX_CONSTANT      = 2
  integer, parameter, public :: VERTEX_SUM           = 3
  integer, parameter, public :: VERTEX_DIFFERENCE    = 4
  integer, parameter, public :: VERTEX_PRODUCT       = 5
  integer, parameter, public :: VERTEX_QUOTIENT      = 6
  integer, parameter, public :: VERTEX_INTEGER_POWER = 7
  integer, parameter, public :: VERTEX_REAL_POWER    = 8
  integer, parameter, public :: VERTEX_FUNCTION      = 9
  ! the derivative of the vertex read, of an order along a coordinate:
  ! an operator, replaced by jet components before the law is evaluated
  integer, parameter, public :: VERTEX_DERIVATIVE    = 10

  ! The coordinate a derivative follows, named by its place among the
  ! coordinates the state is declared over. One is the first declared;
  ! nothing here specifies which that is, and nothing here limits how
  ! many there are.
  integer, parameter, public :: FIRST_COORDINATE = 1

  ! The arguments a leaf reads, in the operation's order; a multiplier
  ! is read from no argument, and a coordinate from the position of
  ! the point, and the position names each
  integer, parameter, public :: ARGUMENT_STATE      = 1
  integer, parameter, public :: ARGUMENT_DESIGN     = 2
  integer, parameter, public :: ARGUMENT_MULTIPLIER = 3
  integer, parameter, public :: ARGUMENT_COORDINATE = 4
  ! the time map t(tau), the position of the instant whose reference
  ! position is the coordinate of time: read by the pin of the
  ! instants t(tau) - tau alone, which the discrete image states as
  ! the condition fixing each instant at its reference position
  integer, parameter, public :: ARGUMENT_TIME_MAP   = 5

  ! The elementary functions. The step H(a), one where a >= 0 and zero
  ! elsewhere, is of a function of the coordinates alone: its
  ! derivative along a coordinate a reads is a measure, not a
  ! function, and is refused, and along any other zero.
  integer, parameter, public :: SINE        = 1
  integer, parameter, public :: COSINE      = 2
  integer, parameter, public :: EXPONENTIAL = 3
  integer, parameter, public :: LOGARITHM   = 4
  integer, parameter, public :: SQUARE_ROOT = 5
  integer, parameter, public :: UNIT_STEP   = 6

  !-------------------------------------------------------------------!
  ! Type that models one vertex as a value, for a view: the kind, the
  ! two vertices read (0 if none), a leaf's argument, field, order and
  ! coordinate, and a constant or an exponent
  !-------------------------------------------------------------------!

  type, public :: expression_vertex
     integer  :: kind = 0, first = 0, second = 0
     integer  :: position = 0, field = 0, order = 0, along = FIRST_COORDINATE, mixed_order = 0
     real(dp) :: coefficient = 0.0_dp
  end type expression_vertex

  !-------------------------------------------------------------------!
  ! Type that models an expression: its vertex table in evaluation
  ! order, and the domain read from its leaves
  !-------------------------------------------------------------------!

  type, extends(operation) :: expression

     ! The domain, read from the leaves at every construction and
     ! stored so that evaluation indexes the tuple without a pass over
     ! the vertices: one degree per coordinate (the first includes the
     ! component of order zero, the state itself, so it stores
     ! degrees(1) + 1 components; every later coordinate stores
     ! degrees(c), its orders running from one); the state fields,
     ! each with its degree along the first coordinate; the
     ! multipliers, absent from the tuple and zero in value; and
     ! varied, the multiplier whose stationarity the rule is, or zero.
     integer , allocatable, private :: degrees(:)
     integer              , private :: fields      = 0
     integer              , private :: multipliers = 0
     ! the designs the expression reads, one value per design at a
     ! point: the highest index named, one at least
     integer              , private :: designs     = 1
     integer              , private :: varied      = 0
     ! the coordinates the expression reads from the position of a
     ! point: the highest index named, zero when none is read
     integer              , private :: position_coordinates = 0
     integer , allocatable, private :: field_degree(:)
     ! per field, the highest order along the first coordinate of a
     ! mixed component, zero when the field has none
     integer , allocatable, private :: mixed_degree_of_field(:)
     ! per field, its first component in the tuple, and after the last
     ! field the count of components
     integer , allocatable, private :: field_offset(:)
     ! per vertex, the component of the tuple a state or multiplier
     ! leaf reads, formed with the domain; -1 for any other vertex
     integer , allocatable, private :: leaf_component(:)

     ! The vertex table, one entry per vertex in evaluation order
     integer , allocatable, private :: kind(:)
     integer , allocatable, private :: first(:), second(:)    ! the vertices read; 0 if none
     integer , allocatable, private :: position(:)            ! a leaf's argument
     integer , allocatable, private :: field(:)               ! a state leaf's field, a multiplier's index
     integer , allocatable, private :: order(:)               ! a leaf's component, or a function index
     integer , allocatable, private :: along(:)               ! the coordinate a leaf's derivative follows
     integer , allocatable, private :: mixed_order(:)         ! a leaf's order along the first coordinate, too
     real(dp), allocatable, private :: coefficient(:)         ! a constant, or an exponent

   contains

     ! The operation's symbols and the evaluation
     procedure :: apply           => expression_apply
     procedure :: defined_at_zero => expression_defined_at_zero
     procedure :: partial_action  => expression_partial_action
     procedure :: at_instant      => expression_at_instant
     procedure, private :: evaluate => expression_evaluate
     procedure :: value_at
     procedure :: gradient_at
     procedure :: value_on
     procedure, private :: values_over

     ! The queries of the leaves, the tuple and the domain
     procedure :: reads_state
     procedure :: reads_time_map
     procedure :: pins_instants
     procedure :: equation_degree
     procedure :: num_coordinates
     procedure :: num_position_coordinates
     procedure :: num_components
     procedure :: component_at
     procedure :: highest_degree_along
     procedure :: highest_order_read
     procedure :: read_components
     procedure :: root
     procedure :: vertices
     procedure :: num_vertices
     procedure :: num_fields
     procedure :: num_multipliers
     procedure :: num_designs
     procedure :: leaf_field
     procedure :: leaf_order
     procedure :: leaf_along
     procedure :: degree_along
     procedure :: degree_of_field
     procedure :: offset_of_field
     procedure :: mixed_degree

     ! The layout of the tuple, the evaluation over derivative terms,
     ! and the domain
     procedure, private :: leaf_at
     procedure, private :: components_per_field
     procedure, private :: stored_at
     procedure, private :: order_at
     procedure, private :: evaluated_over
     procedure, private :: derive_domain
     procedure, private :: declare_domain

  end type expression

  interface operator(+)
     module procedure plus, real_plus, plus_real
  end interface operator(+)

  interface operator(-)
     module procedure minus, real_minus, minus_real, negated
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

  interface heaviside
     module procedure step_of
  end interface heaviside

contains

  !===================================================================!
  ! Creates an expression of one vertex, a leaf or a constant, and
  ! reads its domain
  !===================================================================!

  function vertex(kind, position, order, coefficient, along, field, mixed_order) result(this)

    integer , intent(in)           :: kind, position, order
    real(dp), intent(in)           :: coefficient
    integer , intent(in), optional :: along, field, mixed_order
    type(expression)               :: this

    this % kind        = [kind]
    this % first       = [0]
    this % second      = [0]
    this % position    = [position]
    this % order       = [order]
    this % coefficient = [coefficient]
    this % along       = [FIRST_COORDINATE]
    this % field       = [0]
    this % mixed_order = [0]
    if (present(along)) this % along = [along]
    if (present(field)) this % field = [field]
    if (present(mixed_order)) this % mixed_order = [mixed_order]

    call this % derive_domain()

  end function vertex

  !===================================================================!
  ! Returns the i-th zeroth-order field, the first when none is named.
  ! An index below one stops the program.
  !===================================================================!

  function unknown(i) result(this)

    integer, intent(in), optional :: i
    type(expression)              :: this

    integer            :: field
    character(len=250) :: message

    field = 1
    if (present(i)) field = i
    if (field < 1) then
       write(message,'(a,i0)') 'operation_expression: a field must be named by a positive index; field = ', field
       error stop trim(message)
    end if

    this = vertex(VERTEX_LEAF, ARGUMENT_STATE, 0, 0.0_dp, field=field)

  end function unknown

  !===================================================================!
  ! Returns the i-th design: read from the design of the point, a
  ! vector of one value per design coordinate. An index below one
  ! stops the program; the first when none is given.
  !===================================================================!

  function design(field) result(this)

    integer, intent(in), optional :: field
    type(expression)              :: this

    integer            :: f
    character(len=250) :: message

    f = 1
    if (present(field)) f = field
    if (f < 1) then
       write(message,'(a,i0)') 'operation_expression: a design must be named by a positive index; field = ', f
       error stop trim(message)
    end if

    this = vertex(VERTEX_LEAF, ARGUMENT_DESIGN, 0, 0.0_dp, field=f)

  end function design

  !===================================================================!
  ! Returns the constant c
  !===================================================================!

  function constant(c) result(this)

    real(dp), intent(in) :: c
    type(expression)     :: this

    this = vertex(VERTEX_CONSTANT, 0, 0, c)

  end function constant

  !===================================================================!
  ! Returns the j-th multiplier: read from no argument, zero in value
  ! unless its stationarity is read. An index below one stops the
  ! program.
  !===================================================================!

  function multiplier(j) result(this)

    integer, intent(in) :: j
    type(expression)    :: this

    character(len=250) :: message

    if (j < 1) then
       write(message,'(a,i0)') 'operation_expression: a multiplier must be named by a positive index; j = ', j
       error stop trim(message)
    end if

    this = vertex(VERTEX_LEAF, ARGUMENT_MULTIPLIER, 0, 0.0_dp, field=j)

  end function multiplier

  !===================================================================!
  ! Returns the laws paired with their multipliers, multiplier(1) *
  ! laws(1) + multiplier(2) * laws(2) + ..., formed in one pass: the
  ! vertices the products and the sums formed one at a time would
  ! have, in their order, without the copy of the sum each would make.
  !
  !      lambda_1 | c_1 ... | x | lambda_2 | c_2 ... | x | + | ...
  !                           ^                        ^   ^
  !               lambda_1 c_1             lambda_2 c_2   their sum
  !===================================================================!

  function paired(laws) result(this)

    type(expression), intent(in) :: laws(:)
    type(expression)             :: this

    integer :: total, j, n, at, leaf, root

    total = sum([(laws(j) % root() + 2, j = 1, size(laws))]) + size(laws) - 1
    allocate(this % kind(total), this % first(total), this % second(total), this % position(total), &
         & this % field(total), this % order(total), this % along(total), this % mixed_order(total), &
         & this % coefficient(total))
    this % first       = 0
    this % second      = 0
    this % position    = 0
    this % field       = 0
    this % order       = 0
    this % along       = FIRST_COORDINATE
    this % mixed_order = 0
    this % coefficient = 0.0_dp
    at   = 0
    root = 0
    do j = 1, size(laws)
       ! the multiplier
       at = at + 1
       leaf = at
       this % kind(at)     = VERTEX_LEAF
       this % position(at) = ARGUMENT_MULTIPLIER
       this % field(at)    = j
       ! the law
       n = laws(j) % root()
       this % kind(at + 1:at + n)        = laws(j) % kind
       this % first(at + 1:at + n)       = shifted(laws(j) % first, at)
       this % second(at + 1:at + n)      = shifted(laws(j) % second, at)
       this % position(at + 1:at + n)    = laws(j) % position
       this % field(at + 1:at + n)       = laws(j) % field
       this % order(at + 1:at + n)       = laws(j) % order
       this % along(at + 1:at + n)       = laws(j) % along
       this % mixed_order(at + 1:at + n) = laws(j) % mixed_order
       this % coefficient(at + 1:at + n) = laws(j) % coefficient
       at = at + n
       ! their product, and its sum with the terms before
       at = at + 1
       this % kind(at)   = VERTEX_PRODUCT
       this % first(at)  = leaf
       this % second(at) = at - 1
       if (j > 1) then
          at = at + 1
          this % kind(at)   = VERTEX_SUM
          this % first(at)  = root
          this % second(at) = at - 1
       end if
       root = at
    end do
    call this % derive_domain()

  end function paired

  !===================================================================!
  ! Returns the c-th coordinate function: read from the position of the
  ! point the expression is evaluated at. An index below one stops the
  ! program.
  !===================================================================!

  function coordinate(c) result(this)

    integer, intent(in) :: c
    type(expression)    :: this

    character(len=250) :: message

    if (c < 1) then
       write(message,'(a,i0)') 'operation_expression: a coordinate must be named by a positive index; c = ', c
       error stop trim(message)
    end if

    this = vertex(VERTEX_LEAF, ARGUMENT_COORDINATE, 0, 0.0_dp, field=c)

  end function coordinate

  !===================================================================!
  ! Returns the time map t(tau), a function of the coordinate of time
  !===================================================================!

  function time_map() result(this)

    type(expression) :: this

    this = vertex(VERTEX_LEAF, ARGUMENT_TIME_MAP, 0, 0.0_dp, field=FIRST_COORDINATE)

  end function time_map

  !===================================================================!
  ! Returns the component of degree d of the unknown. Anything but a
  ! bare unknown, or a negative degree, stops the program: along the
  ! instants only the unknown has components.
  !===================================================================!

  function component(x, d) result(this)

    type(expression), intent(in) :: x
    integer         , intent(in) :: d
    type(expression)             :: this

    integer            :: n
    character(len=250) :: message

    n = x % root()
    if (n /= 1 .or. x % kind(1) /= VERTEX_LEAF .or. x % position(1) /= ARGUMENT_STATE) then
       write(message,'(a,i0,a,i0,a,i0)') 'operation_expression: component requires its &
            &argument to be the unknown; num_vertices = ', n, ', x % kind(1) = ', &
            & x % kind(1), ', x % position(1) = ', x % position(1)
       error stop trim(message)
    end if
    if (d < 0) then
       write(message,'(a,i0)') 'operation_expression: the degree of a derivative must not be &
            &negative; d = ', d
       error stop trim(message)
    end if

    this = vertex(VERTEX_LEAF, ARGUMENT_STATE, d, 0.0_dp, FIRST_COORDINATE, x % field(1))

  end function component

  !===================================================================!
  ! Returns the same component named along any coordinate: which one
  ! is given by its place among those declared, so a rule reads no
  ! differently whether the coordinate is a time, a length or a
  ! sample. With a mixed order the component is also of that order
  ! along the first coordinate.
  !
  !      component_along(u, x, 2)      u_xx
  !      component_along(u, x, 2, 1)   u_txx, with t the first
  !
  ! Invalid input: a mixed order below zero, or one along the first
  ! coordinate itself.
  !===================================================================!

  function component_along(x, coordinate, d, mixed_order) result(this)

    type(expression), intent(in)           :: x
    integer         , intent(in)           :: coordinate, d
    integer         , intent(in), optional :: mixed_order
    type(expression)                       :: this

    integer            :: n
    character(len=250) :: message

    n = x % root()
    if (n /= 1 .or. x % kind(1) /= VERTEX_LEAF .or. x % position(1) /= ARGUMENT_STATE) then
       write(message,'(a,i0,a,i0,a,i0)') 'operation_expression: component_along requires its &
            &argument to be the unknown; num_vertices = ', n, ', x % kind(1) = ', &
            & x % kind(1), ', x % position(1) = ', x % position(1)
       error stop trim(message)
    end if
    if (d < 0) then
       write(message,'(a,i0)') 'operation_expression: the degree of a derivative must not be &
            &negative; d = ', d
       error stop trim(message)
    end if
    if (coordinate < FIRST_COORDINATE) then
       write(message,'(a,i0,a,i0)') 'operation_expression: a coordinate must be one of those &
            &declared; coordinate = ', coordinate, ', FIRST_COORDINATE = ', FIRST_COORDINATE
       error stop trim(message)
    end if

    this = vertex(VERTEX_LEAF, ARGUMENT_STATE, d, 0.0_dp, coordinate, x % field(1))
    if (.not. present(mixed_order)) return
    if (mixed_order < 0 .or. (mixed_order > 0 .and. (coordinate == FIRST_COORDINATE .or. d < 1))) then
       error stop 'operation_expression: a mixed component is of positive order along a later &
            &coordinate and of nonnegative order along the first'
    end if
    this % mixed_order = [mixed_order]

  end function component_along

  !===================================================================!
  ! Returns the derivative operator of an order along a coordinate,
  ! applied to any expression and stored in the law as a vertex. The
  ! law is not evaluated until jet_form has replaced every such vertex
  ! by the jet components it denotes. Invalid input: an order below
  ! one, or a coordinate before the first.
  !===================================================================!

  function derivative(x, along, order) result(this)

    type(expression), intent(in) :: x
    integer         , intent(in) :: along, order
    type(expression)             :: this

    character(len=250) :: message

    if (order < 1 .or. along < FIRST_COORDINATE) then
       write(message,'(a,i0,a,i0)') 'operation_expression: a derivative operator is of order one or &
            &more along a declared coordinate; order = ', order, ', along = ', along
       error stop trim(message)
    end if
    this = applied(x, VERTEX_DERIVATIVE, order, 0.0_dp)
    this % along(size(this % along)) = along

  end function derivative

  !===================================================================!
  ! Returns the jet form of a law: every derivative operator replaced
  ! by what it denotes in the jet, the total derivative of the
  ! expression it reads by the chain rule, order times along its
  ! coordinate. A law with no operator is returned as it is. A mixed
  ! derivative of a component stops the program, as the total
  ! derivative does.
  !
  !      derivative(u_t * u, x, 1)   -->   u_tx u + u_t u_x
  !===================================================================!

  function jet_form(x) result(this)

    type(expression), intent(in) :: x
    type(expression)             :: this

    if (.not. any(x % kind == VERTEX_DERIVATIVE)) then
       this = x
       return
    end if
    this = subtree(x, x % root(), resolved=.true.)

  end function jet_form

  !===================================================================!
  ! Returns the component a representation law defines: for a law
  ! d^b(q) - p, an operator applied to a jet component less the jet
  ! component the operator denotes, the field, the order and the
  ! coordinate of p; zeros for any other law.
  !
  !      D_t u - u_(1)   -->   [field of u, 1, the first coordinate]
  !
  ! A law that is not the difference of operators on a leaf and a
  ! leaf, whose jet form is not the difference of two equal state
  ! components, or that has no operator, is not a representation law;
  ! the shape is read before the jet form, so that no other law is
  ! differentiated.
  !===================================================================!

  function represented_component(x) result(defined)

    type(expression), intent(in) :: x
    integer                      :: defined(3)

    type(expression) :: form
    integer          :: r, a, b

    defined = 0
    if (.not. any(x % kind == VERTEX_DERIVATIVE)) return
    r = x % root()
    if (x % kind(r) /= VERTEX_DIFFERENCE) return
    a = x % first(r)
    if (x % kind(x % second(r)) /= VERTEX_LEAF) return
    do while (x % kind(a) == VERTEX_DERIVATIVE)
       a = x % first(a)
    end do
    if (x % kind(a) /= VERTEX_LEAF) return
    form = jet_form(x)
    r = form % root()
    if (form % kind(r) /= VERTEX_DIFFERENCE) return
    a = form % first(r)
    b = form % second(r)
    if (form % kind(a) /= VERTEX_LEAF .or. form % kind(b) /= VERTEX_LEAF) return
    if (form % position(a) /= ARGUMENT_STATE .or. form % position(b) /= ARGUMENT_STATE) return
    if (form % field(a) /= form % field(b) .or. form % order(a) /= form % order(b)) return
    if (form % order(b) < 1 .or. form % along(a) /= form % along(b)) return
    if (form % mixed_order(a) /= form % mixed_order(b)) return
    if (form % mixed_order(b) > 0) then
       error stop 'operation_expression: a representation law of a mixed component is not stated; &
            &the concretion along space defines the mixed components'
    end if
    defined = [form % field(b), form % order(b), form % along(b)]

  end function represented_component

  !===================================================================!
  ! Returns the stationarity of a Lagrangian in its j-th multiplier:
  ! the partial of the rule along that multiplier, at zero. When the
  ! Lagrangian is a sum of terms, each a multiplier times a law free of
  ! multipliers, and the j-th multiplier is read in one term alone, the
  ! partial is that term's law: its vertices, on the domain and the
  ! arguments of the Lagrangian, so that it reads the Lagrangian's
  ! tuple. A multiplier the rule does not read stops the program.
  !
  !      L = lambda_1 c_1 + lambda_2 c_2
  !      euler_lagrange(L, 2)   c_2, its own vertices on L's domain
  !      any other form of L    L with lambda_2 varied, dL/dlambda_2
  !===================================================================!

  function euler_lagrange(lagrangian, multiplier, label) result(this)

    type(expression), intent(in)           :: lagrangian
    integer         , intent(in)           :: multiplier
    character(len=*), intent(in), optional :: label
    type(expression)                       :: this

    character(len=250) :: message
    integer            :: law, start

    if (multiplier < 1 .or. multiplier > lagrangian % multipliers) then
       write(message,'(a,i0,a,i0)') 'operation_expression: the stationarity must be in a &
            &multiplier the Lagrangian reads; multiplier = ', multiplier, &
            & ', lagrangian % multipliers = ', lagrangian % multipliers
       error stop trim(message)
    end if

    law = factor_of(lagrangian, multiplier)
    if (law > 0) then
       ! the vertices of a subtree are consecutive, its root last
       start = law
       do while (lagrangian % first(start) > 0)
          start = lagrangian % first(start)
       end do
       this % kind        = lagrangian % kind(start:law)
       this % first       = shifted(lagrangian % first(start:law), 1 - start)
       this % second      = shifted(lagrangian % second(start:law), 1 - start)
       this % position    = lagrangian % position(start:law)
       this % field       = lagrangian % field(start:law)
       this % order       = lagrangian % order(start:law)
       this % along       = lagrangian % along(start:law)
       this % mixed_order = lagrangian % mixed_order(start:law)
       this % coefficient = lagrangian % coefficient(start:law)
       this % fields                = lagrangian % fields
       this % multipliers           = lagrangian % multipliers
       this % designs               = lagrangian % designs
       this % position_coordinates  = lagrangian % position_coordinates
       this % field_degree          = lagrangian % field_degree
       this % mixed_degree_of_field = lagrangian % mixed_degree_of_field
       this % degrees               = lagrangian % degrees
       call this % share_arguments(lagrangian)
       call this % declare_domain(label)
       return
    end if
    this = lagrangian
    this % varied = multiplier
    call this % derive_domain(label)

  end function euler_lagrange

  !===================================================================!
  ! Returns the vertex of the law the j-th multiplier multiplies, in a
  ! Lagrangian whose terms are summed and the j-th multiplier read
  ! once, as a factor of a term whose other factor reads no
  ! multiplier; zero for a Lagrangian of any other form
  !===================================================================!

  function factor_of(lagrangian, j) result(law)

    type(expression), intent(in) :: lagrangian
    integer         , intent(in) :: j
    integer                      :: law

    integer, allocatable :: pending(:)
    integer              :: i, other, start

    law = 0
    if (count(lagrangian % kind == VERTEX_LEAF .and. lagrangian % position == ARGUMENT_MULTIPLIER &
         & .and. lagrangian % field == j) /= 1) return
    pending = [lagrangian % root()]
    do while (size(pending) > 0)
       i = pending(size(pending))
       pending = pending(1:size(pending) - 1)
       select case (lagrangian % kind(i))
       case (VERTEX_SUM)
          pending = [pending, lagrangian % first(i), lagrangian % second(i)]
       case (VERTEX_PRODUCT)
          if (read_multiplier(lagrangian % first(i))) then
             other = lagrangian % second(i)
          else if (read_multiplier(lagrangian % second(i))) then
             other = lagrangian % first(i)
          else
             cycle
          end if
          start = other
          do while (lagrangian % first(start) > 0)
             start = lagrangian % first(start)
          end do
          if (any(lagrangian % kind(start:other) == VERTEX_LEAF .and. &
               & lagrangian % position(start:other) == ARGUMENT_MULTIPLIER)) return
          law = other
          return
       end select
    end do

  contains

    !-----------------------------------------------------------------!
    ! Returns whether vertex k is the j-th multiplier
    !-----------------------------------------------------------------!

    logical function read_multiplier(k)

      integer, intent(in) :: k

      read_multiplier = lagrangian % kind(k) == VERTEX_LEAF .and. lagrangian % position(k) == ARGUMENT_MULTIPLIER &
           & .and. lagrangian % field(k) == j

    end function read_multiplier

  end function factor_of

  !===================================================================!
  ! Returns the Lagrangian with every multiplier at zero, which for a
  ! Lagrangian linear in them is the functional beside the constraints
  !===================================================================!

  function at_zero(lagrangian, label) result(this)

    type(expression), intent(in)           :: lagrangian
    character(len=*), intent(in), optional :: label
    type(expression)                       :: this

    this = lagrangian
    this % varied = 0
    call this % derive_domain(label)

  end function at_zero

  !===================================================================!
  ! Returns a binary vertex of kind reading two roots: the right
  ! operand's vertices are appended after the left's and renumbered,
  ! and the new root reads both
  !
  !      a: 1 .. na | b: na + 1 .. na + nb | root: reads na, na + nb
  !===================================================================!

  function joined(a, b, kind) result(this)

    type(expression), intent(in) :: a, b
    integer         , intent(in) :: kind
    type(expression)             :: this

    integer :: na, nb

    na = a % root()
    nb = b % root()

    this % kind        = [a % kind,        b % kind,        kind]
    this % first       = [a % first,       shifted(b % first,  na), na]
    this % second      = [a % second,      shifted(b % second, na), na + nb]
    this % position    = [a % position,    b % position,    0]
    this % field       = [a % field,       b % field,       0]
    this % order       = [a % order,       b % order,       0]
    this % along       = [a % along,       b % along,       FIRST_COORDINATE]
    this % mixed_order = [a % mixed_order, b % mixed_order, 0]
    this % coefficient = [a % coefficient, b % coefficient, 0.0_dp]

    call this % derive_domain()

  end function joined

  !===================================================================!
  ! Returns the copy of the subtree at i; resolved, every derivative
  ! operator replaced by the total derivatives it denotes, as jet_form
  ! requires. Symbolic differentiation forms new expressions over the
  ! same leaves by recursion over the tree with the operators of this
  ! module, a zero or a one met on the way absorbed.
  !===================================================================!

  recursive function subtree(x, i, resolved) result(this)

    type(expression), intent(in)           :: x
    integer         , intent(in)           :: i
    logical         , intent(in), optional :: resolved
    type(expression)                       :: this

    integer :: k, start
    logical :: jet

    jet = .false.
    if (present(resolved)) jet = resolved
    ! a subtree free of derivative operators is its consecutive
    ! vertices, its root last, renumbered: the vertices the joins below
    ! would form, in their order
    start = i
    do while (x % first(start) > 0)
       start = x % first(start)
    end do
    if (.not. any(x % kind(start:i) == VERTEX_DERIVATIVE)) then
       this % kind        = x % kind(start:i)
       this % first       = shifted(x % first(start:i), 1 - start)
       this % second      = shifted(x % second(start:i), 1 - start)
       this % position    = x % position(start:i)
       this % field       = x % field(start:i)
       this % order       = x % order(start:i)
       this % along       = x % along(start:i)
       this % mixed_order = x % mixed_order(start:i)
       this % coefficient = x % coefficient(start:i)
       call this % derive_domain()
       return
    end if
    select case (x % kind(i))
    case (VERTEX_LEAF)
       this = vertex(VERTEX_LEAF, x % position(i), x % order(i), 0.0_dp, x % along(i), x % field(i), &
            & x % mixed_order(i))
    case (VERTEX_CONSTANT)
       this = constant(x % coefficient(i))
    case (VERTEX_SUM, VERTEX_DIFFERENCE, VERTEX_PRODUCT, VERTEX_QUOTIENT)
       this = joined(subtree(x, x % first(i), resolved), subtree(x, x % second(i), resolved), x % kind(i))
    case (VERTEX_DERIVATIVE)
       this = subtree(x, x % first(i), resolved)
       if (.not. jet) then
          this = derivative(this, x % along(i), x % order(i))
          return
       end if
       do k = 1, x % order(i)
          this = total_derivative(this, x % along(i))
       end do
    case default
       this = applied(subtree(x, x % first(i), resolved), x % kind(i), x % order(i), x % coefficient(i))
    end select

  end function subtree

  !===================================================================!
  ! Returns whether x is the constant c alone
  !===================================================================!

  pure logical function is_constant(x, c)

    type(expression), intent(in) :: x
    real(dp)        , intent(in) :: c

    is_constant = .false.
    if (x % root() == 1) then
       if (x % kind(1) == VERTEX_CONSTANT) is_constant = x % coefficient(1) == c
    end if

  end function is_constant

  !===================================================================!
  ! Returns a b with a zero or a one absorbed
  !===================================================================!

  function product_of(a, b) result(this)

    type(expression), intent(in) :: a, b
    type(expression)             :: this

    if (is_constant(a, 0.0_dp) .or. is_constant(b, 0.0_dp)) then
       this = constant(0.0_dp)
    else if (is_constant(a, 1.0_dp)) then
       this = b
    else if (is_constant(b, 1.0_dp)) then
       this = a
    else
       this = a * b
    end if

  end function product_of

  !===================================================================!
  ! Returns a + b with a zero absorbed
  !===================================================================!

  function sum_of(a, b) result(this)

    type(expression), intent(in) :: a, b
    type(expression)             :: this

    if (is_constant(a, 0.0_dp)) then
       this = b
    else if (is_constant(b, 0.0_dp)) then
       this = a
    else
       this = a + b
    end if

  end function sum_of

  !===================================================================!
  ! Returns a - b with a zero absorbed
  !===================================================================!

  function difference_of(a, b) result(this)

    type(expression), intent(in) :: a, b
    type(expression)             :: this

    if (is_constant(b, 0.0_dp)) then
       this = a
    else if (is_constant(a, 0.0_dp)) then
       this = -b
    else
       this = a - b
    end if

  end function difference_of

  !===================================================================!
  ! Returns the derivative of the subtree at i, the leaves
  ! differentiated by the mode
  !
  !      mode   the derivative
  !      1      the partial in the state component (field, order,
  !             along)
  !      2      the partial in the design field
  !      3      the total derivative along the coordinate along
  !
  ! A subtree that reads no leaf the mode differentiates, reads(i)
  ! false, has the derivative zero, and an operand is copied only into
  ! a term that is not zero.
  !===================================================================!

  recursive function differentiated(x, i, mode, field, order, along, reads) result(this)

    type(expression), intent(in) :: x
    integer         , intent(in) :: i, mode, field, order, along
    logical         , intent(in) :: reads(:)
    type(expression)             :: this

    type(expression)   :: a, b, da, db
    integer            :: n
    character(len=250) :: message

    if (.not. reads(i)) then
       this = constant(0.0_dp)
       return
    end if
    select case (x % kind(i))
    case (VERTEX_DERIVATIVE)
       error stop 'operation_expression: a law with a derivative operator is differentiated only in its &
            &jet form; jet_form replaces the operator by the jet components it denotes'
    case (VERTEX_LEAF)
       this = constant(0.0_dp)
       select case (mode)
       case (1)
          if (x % position(i) == ARGUMENT_STATE) then
             if (x % field(i) == field .and. x % order(i) == order .and. x % mixed_order(i) == 0 .and. &
                  & (x % along(i) == along .or. order == 0)) this = constant(1.0_dp)
          end if
       case (2)
          if (x % position(i) == ARGUMENT_DESIGN) then
             if (x % field(i) == field) this = constant(1.0_dp)
          end if
       case (3)
          if (x % position(i) == ARGUMENT_STATE) then
             ! a mixed component is of order mixed_order along the
             ! first coordinate and order along a later one
             if (x % order(i) == 0) then
                this = vertex(VERTEX_LEAF, ARGUMENT_STATE, 1, 0.0_dp, along, x % field(i))
             else if (x % along(i) == along) then
                this = vertex(VERTEX_LEAF, ARGUMENT_STATE, x % order(i) + 1, 0.0_dp, along, x % field(i), &
                     & x % mixed_order(i))
             else if (x % along(i) == FIRST_COORDINATE) then
                this = vertex(VERTEX_LEAF, ARGUMENT_STATE, 1, 0.0_dp, along, x % field(i), x % order(i))
             else if (along == FIRST_COORDINATE) then
                this = vertex(VERTEX_LEAF, ARGUMENT_STATE, x % order(i), 0.0_dp, x % along(i), x % field(i), &
                     & x % mixed_order(i) + 1)
             else
                write(message,'(a,i0,a,i0)') 'operation_expression: the total derivative along coordinate ', &
                     & along, ' of a component along coordinate ', x % along(i)
                error stop trim(message) // ' is a mixed derivative of two later coordinates, which is not &
                     &a component of the jet'
             end if
          else if (x % position(i) == ARGUMENT_COORDINATE) then
             if (x % field(i) == along) this = constant(1.0_dp)
          end if
       end select
    case (VERTEX_CONSTANT)
       this = constant(0.0_dp)
    case (VERTEX_SUM)
       this = sum_of(differentiated(x, x % first(i), mode, field, order, along, reads), &
            &        differentiated(x, x % second(i), mode, field, order, along, reads))
    case (VERTEX_DIFFERENCE)
       this = difference_of(differentiated(x, x % first(i), mode, field, order, along, reads), &
            &               differentiated(x, x % second(i), mode, field, order, along, reads))
    case (VERTEX_PRODUCT)
       da = differentiated(x, x % first(i), mode, field, order, along, reads)
       db = differentiated(x, x % second(i), mode, field, order, along, reads)
       if (is_constant(da, 0.0_dp) .and. is_constant(db, 0.0_dp)) then
          this = constant(0.0_dp)
       else if (is_constant(da, 0.0_dp)) then
          this = product_of(subtree(x, x % first(i)), db)
       else if (is_constant(db, 0.0_dp)) then
          this = product_of(da, subtree(x, x % second(i)))
       else
          a  = subtree(x, x % first(i))
          b  = subtree(x, x % second(i))
          this = sum_of(product_of(da, b), product_of(a, db))
       end if
    case (VERTEX_QUOTIENT)
       da = differentiated(x, x % first(i), mode, field, order, along, reads)
       db = differentiated(x, x % second(i), mode, field, order, along, reads)
       if (is_constant(da, 0.0_dp) .and. is_constant(db, 0.0_dp)) then
          this = constant(0.0_dp)
          return
       end if
       a  = subtree(x, x % first(i))
       b  = subtree(x, x % second(i))
       this = difference_of(product_of(da, b), product_of(a, db))
       if (.not. is_constant(this, 0.0_dp)) this = this / (b * b)
    case (VERTEX_INTEGER_POWER)
       n  = nint(x % coefficient(i))
       da = differentiated(x, x % first(i), mode, field, order, along, reads)
       if (is_constant(da, 0.0_dp)) then
          this = constant(0.0_dp)
          return
       end if
       a  = subtree(x, x % first(i))
       if (n == 0) then
          this = constant(0.0_dp)
       else if (n == 1) then
          this = da
       else if (n == 2) then
          this = product_of(product_of(constant(2.0_dp), a), da)
       else
          this = product_of(product_of(constant(real(n, dp)), a ** (n - 1)), da)
       end if
    case (VERTEX_REAL_POWER)
       da = differentiated(x, x % first(i), mode, field, order, along, reads)
       if (is_constant(da, 0.0_dp)) then
          this = constant(0.0_dp)
          return
       end if
       a  = subtree(x, x % first(i))
       this = product_of(product_of(constant(x % coefficient(i)), a ** (x % coefficient(i) - 1.0_dp)), da)
    case (VERTEX_FUNCTION)
       da = differentiated(x, x % first(i), mode, field, order, along, reads)
       if (is_constant(da, 0.0_dp)) then
          this = constant(0.0_dp)
          return
       end if
       a  = subtree(x, x % first(i))
       select case (x % order(i))
       case (SINE)
          this = product_of(cos(a), da)
       case (COSINE)
          this = product_of(-sin(a), da)
       case (EXPONENTIAL)
          this = product_of(exp(a), da)
       case (LOGARITHM)
          if (is_constant(da, 0.0_dp)) then
             this = constant(0.0_dp)
          else
             this = da / a
          end if
       case (SQUARE_ROOT)
          if (is_constant(da, 0.0_dp)) then
             this = constant(0.0_dp)
          else
             this = da / (constant(2.0_dp) * sqrt(a))
          end if
       case (UNIT_STEP)
          error stop 'operation_expression: the derivative of a step along a variable its argument reads is a &
               &measure, not a function'
       case default
          error stop 'operation_expression: the function is not one this module differentiates'
       end select
    case default
       error stop 'operation_expression: the vertex kind is not one this module differentiates'
    end select

  end function differentiated

  !===================================================================!
  ! Returns the partial in the state component of the field of that
  ! order along the coordinate, every other leaf unchanged
  !===================================================================!

  function partial_derivative(x, field, order, along) result(this)

    type(expression), intent(in)           :: x
    integer         , intent(in)           :: field, order
    integer         , intent(in), optional :: along
    type(expression)                       :: this

    integer :: c

    c = FIRST_COORDINATE
    if (present(along)) c = along
    this = differentiated(x, x % root(), 1, field, order, c, leaves_read(x, 1, field, order, c))

  end function partial_derivative

  !===================================================================!
  ! Returns the partial in the i-th design
  !===================================================================!

  function partial_in_design(x, field) result(this)

    type(expression), intent(in) :: x
    integer         , intent(in) :: field
    type(expression)             :: this

    this = differentiated(x, x % root(), 2, field, 0, FIRST_COORDINATE, leaves_read(x, 2, field, 0, FIRST_COORDINATE))

  end function partial_in_design

  !===================================================================!
  ! Returns the total derivative along a coordinate c by the chain
  ! rule over every leaf
  !
  !      leaf                            its derivative along c
  !      u, order 0                      u_c
  !      order o along c                 order o + 1 along c
  !      order o along the first,        order 1 along c, mixed o
  !        c later
  !      order o along c' later,         order o along c', mixed
  !        c the first                   one more
  !      along c' later, c later         refused: a mixed derivative of
  !        and another                   two later coordinates is not a
  !                                      component of the jet
  !      the coordinate c                one
  !      any other leaf                  zero
  !===================================================================!

  function total_derivative(x, along) result(this)

    type(expression), intent(in) :: x
    integer         , intent(in) :: along
    type(expression)             :: this

    this = differentiated(x, x % root(), 3, 0, 0, along, leaves_read(x, 3, 0, 0, along))

  end function total_derivative

  !===================================================================!
  ! Returns whether the subtree at each vertex reads a leaf that the
  ! mode of differentiated differentiates to other than zero
  !===================================================================!

  pure function leaves_read(x, mode, field, order, along) result(reads)

    type(expression), intent(in) :: x
    integer         , intent(in) :: mode, field, order, along
    logical                      :: reads(x % num_vertices())

    integer :: i

    do i = 1, x % num_vertices()
       select case (x % kind(i))
       case (VERTEX_LEAF)
          select case (mode)
          case (1)
             reads(i) = x % position(i) == ARGUMENT_STATE .and. x % field(i) == field .and. &
                  & x % order(i) == order .and. x % mixed_order(i) == 0 .and. (x % along(i) == along .or. order == 0)
          case (2)
             reads(i) = x % position(i) == ARGUMENT_DESIGN .and. x % field(i) == field
          case default
             reads(i) = x % position(i) == ARGUMENT_STATE .or. &
                  & (x % position(i) == ARGUMENT_COORDINATE .and. x % field(i) == along)
          end select
       case (VERTEX_CONSTANT)
          reads(i) = .false.
       case default
          reads(i) = reads(x % first(i))
          if (x % second(i) > 0) reads(i) = reads(i) .or. reads(x % second(i))
       end select
    end do

  end function leaves_read

  !===================================================================!
  ! Returns the Euler-Lagrange derivative in the field: the sum over
  ! the orders o of (-1)^o times the o-th total derivative along the
  ! coordinate of the partial in the component of order o; the
  ! interior equation of the stationarity of the integral of x, the
  ! boundary terms of the integration by parts aside
  !
  !      E_u(x) = dx/du - D (dx/du_(1)) + D^2 (dx/du_(2)) - ...
  !===================================================================!

  function stationarity(x, field, along) result(this)

    type(expression), intent(in)           :: x
    integer         , intent(in)           :: field
    integer         , intent(in), optional :: along
    type(expression)                       :: this

    type(expression) :: term
    integer          :: c, o, k, top

    c = FIRST_COORDINATE
    if (present(along)) c = along
    if (field < 1 .or. field > x % fields) then
       error stop 'operation_expression: the stationarity is in a field the expression reads'
    end if
    if (any(x % mixed_order > 0)) then
       error stop 'operation_expression: the stationarity of an expression reading a mixed component &
            &is not formed; the partials in the mixed components are not taken'
    end if
    top = max(0, x % highest_degree_along(c, field))
    this = constant(0.0_dp)
    do o = 0, top
       term = partial_derivative(x, field, o, c)
       do k = 1, o
          term = total_derivative(term, c)
       end do
       if (mod(o, 2) == 0) then
          this = sum_of(this, term)
       else
          this = difference_of(this, term)
       end if
    end do

  end function stationarity

  !===================================================================!
  ! Returns the boundary terms of the stationarity of the integral of
  ! x in the field along a coordinate: the first variation, integrated
  ! by parts, leaves at the far end of the coordinate the sum over k of
  ! terms(k + 1) times the variation of the k-th derivative of the
  ! field, terms(k + 1) the sum over o > k of (-1)^(o - k - 1) times
  ! the (o - k - 1)-th total derivative of the partial in the component
  ! of order o. At the near end the same terms enter with the opposite
  ! sign. Where no datum fixes the field's k-th derivative,
  ! terms(k + 1) = 0 is the natural condition there.
  !
  !      top 2   terms(1) = dx/du_(1) - D (dx/du_(2))   of delta u
  !              terms(2) = dx/du_(2)                   of delta u_(1)
  !===================================================================!

  function boundary_terms(x, field, along) result(terms)

    type(expression), intent(in)           :: x
    integer         , intent(in)           :: field
    integer         , intent(in), optional :: along
    type(expression), allocatable          :: terms(:)

    type(expression) :: term
    integer          :: c, o, k, i, top

    c = FIRST_COORDINATE
    if (present(along)) c = along
    if (field < 1 .or. field > x % fields) then
       error stop 'operation_expression: the boundary terms are in a field the expression reads'
    end if
    if (any(x % mixed_order > 0)) then
       error stop 'operation_expression: the boundary terms of an expression reading a mixed component &
            &are not formed; the partials in the mixed components are not taken'
    end if
    top = max(0, x % highest_degree_along(c, field))
    allocate(terms(top))
    do k = 0, top - 1
       terms(k + 1) = constant(0.0_dp)
       do o = k + 1, top
          term = partial_derivative(x, field, o, c)
          do i = 1, o - k - 1
             term = total_derivative(term, c)
          end do
          if (mod(o - k - 1, 2) == 0) then
             terms(k + 1) = sum_of(terms(k + 1), term)
          else
             terms(k + 1) = difference_of(terms(k + 1), term)
          end if
       end do
    end do

  end function boundary_terms

  !===================================================================!
  ! Returns the flux form of a law tested against its multiplier. The
  ! law is a sum of terms, each times a constant; a term is the
  ! derivative of some order along one coordinate of space of a flux,
  ! any expression but a bare leaf, or any other expression.
  ! Integrated by parts, the pairing of the multiplier with the
  ! derivative of order o of a flux is (-1)^o times the pairing of the
  ! multiplier's derivative of that order with the flux. The terms are
  ! grouped by the derivative of the multiplier that tests them, of
  ! order(g) along coordinate(g), the order zero the multiplier itself;
  ! flux(g) is the jet form of the sum of the group's terms with their
  ! signs. The terms of the integration by parts on the faces are left
  ! to the conditions there. space is the first coordinate of space. A
  ! derivative of a bare unknown is not a flux and stays on the
  ! unknown.
  !
  !      law     q_tt - D_x (k q_x) + f
  !
  !      g   tested by       flux(g)
  !      1   lambda          q_tt + f
  !      2   D_x lambda      k q_x       (-1)^1 times the term's -1
  !
  !      lambda law = lambda (q_tt + f) + D_x lambda (k q_x),
  !                   less the terms on the faces
  !
  ! Invalid input: a flux differentiated along two coordinates.
  !===================================================================!

  subroutine flux_form(x, space, coordinate, order, flux)

    type(expression), intent(in)               :: x
    integer         , intent(in)               :: space
    integer         , allocatable, intent(out) :: coordinate(:), order(:)
    type(expression), allocatable, intent(out) :: flux(:)

    integer , allocatable :: term(:), along(:), times(:)
    real(dp), allocatable :: factor(:)
    type(expression)      :: group
    integer               :: k, g

    allocate(term(0), along(0), times(0), factor(0))
    call collect(x % root(), 1.0_dp)

    ! the groups in the order their first term appears, the
    ! multiplier itself first
    allocate(coordinate(0), order(0))
    if (any(times == 0)) then
       coordinate = [0]
       order      = [0]
    end if
    do k = 1, size(term)
       if (times(k) == 0) cycle
       if (any(coordinate == along(k) .and. order == times(k))) cycle
       coordinate = [coordinate, along(k)]
       order      = [order, times(k)]
    end do
    allocate(flux(size(order)))
    do g = 1, size(order)
       group = constant(0.0_dp)
       do k = 1, size(term)
          if (times(k) /= order(g) .or. (times(k) > 0 .and. along(k) /= coordinate(g))) cycle
          group = sum_of(group, product_of(constant(factor(k)), subtree(x, term(k), resolved=.true.)))
       end do
       flux(g) = group
    end do

  contains

    !-----------------------------------------------------------------!
    ! Collects the terms of the subtree at i, each times c: through
    ! sums, differences and products by a constant, and then each
    ! term's derivatives along one coordinate of space peeled off its
    ! flux
    !-----------------------------------------------------------------!

    recursive subroutine collect(i, c)

      integer , intent(in) :: i
      real(dp), intent(in) :: c

      integer :: j, a, o

      select case (x % kind(i))
      case (VERTEX_SUM)
         call collect(x % first(i), c)
         call collect(x % second(i), c)
         return
      case (VERTEX_DIFFERENCE)
         call collect(x % first(i), c)
         call collect(x % second(i), -c)
         return
      case (VERTEX_PRODUCT)
         if (x % kind(x % first(i)) == VERTEX_CONSTANT) then
            call collect(x % second(i), c * x % coefficient(x % first(i)))
            return
         end if
         if (x % kind(x % second(i)) == VERTEX_CONSTANT) then
            call collect(x % first(i), c * x % coefficient(x % second(i)))
            return
         end if
      end select
      ! the derivatives along one coordinate of space, peeled
      j = i
      a = 0
      o = 0
      do while (x % kind(j) == VERTEX_DERIVATIVE)
         if (x % along(j) < space) exit
         if (a /= 0 .and. x % along(j) /= a) then
            error stop 'operation_expression: a flux is differentiated along one coordinate of space; &
                 &a law in flux form is a sum of derivatives of fluxes, each along one coordinate'
         end if
         a = x % along(j)
         o = o + x % order(j)
         j = x % first(j)
      end do
      if (o > 0 .and. .not. (x % kind(j) == VERTEX_LEAF .and. x % position(j) == ARGUMENT_STATE)) then
         term   = [term, j]
         along  = [along, a]
         times  = [times, o]
         factor = [factor, c * real((-1)**o, dp)]
      else
         term   = [term, i]
         along  = [along, 0]
         times  = [times, 0]
         factor = [factor, c]
      end if

    end subroutine collect

  end subroutine flux_form

  !===================================================================!
  ! Returns the coefficients a node stores of every state field of x,
  ! as an expression that reads each once. With time the first
  ! coordinate: the orders along it that x stores, and along each
  ! coordinate of space from space to last the orders one to order,
  ! each with every order along time, the coefficient's own jet in
  ! time. Without time: the orders zero to order along the first
  ! coordinate and one to order along the others. Its layout is the
  ! tuple of a node.
  !
  !      time and x, order 2, u of degree 1 in t
  !      reads   u, u_t, u_x, u_tx, u_xx, u_txx
  !      tuple   [ u  u_t | u_x  u_xx | u_tx  u_txx ]
  !===================================================================!

  function coefficient_form(x, space, last, order) result(this)

    type(expression), intent(in) :: x
    integer         , intent(in) :: space, last, order
    type(expression)             :: this

    integer :: f, a, c, k, top

    this = constant(0.0_dp)
    do f = 1, x % fields
       top = 0
       if (space > FIRST_COORDINATE) top = x % field_degree(f)
       do a = 0, top
          this = this + vertex(VERTEX_LEAF, ARGUMENT_STATE, a, 0.0_dp, along=FIRST_COORDINATE, field=f)
       end do
       if (space == FIRST_COORDINATE) then
          do k = 1, order
             this = this + vertex(VERTEX_LEAF, ARGUMENT_STATE, k, 0.0_dp, along=FIRST_COORDINATE, field=f)
          end do
       end if
       do c = max(space, FIRST_COORDINATE + 1), last
          do k = 1, order
             do a = 0, top
                this = this + vertex(VERTEX_LEAF, ARGUMENT_STATE, k, 0.0_dp, along=c, field=f, mixed_order=a)
             end do
          end do
       end do
    end do

  end function coefficient_form

  !===================================================================!
  ! Returns a unary vertex of kind reading the root of a, with its
  ! order and coefficient: a power, a function or a derivative
  !===================================================================!

  function applied(a, kind, order, coefficient) result(this)

    type(expression), intent(in) :: a
    integer         , intent(in) :: kind, order
    real(dp)        , intent(in) :: coefficient
    type(expression)             :: this

    integer :: na

    na = a % root()

    this % kind        = [a % kind,        kind]
    this % first       = [a % first,       na]
    this % second      = [a % second,      0]
    this % position    = [a % position,    0]
    this % field       = [a % field,       0]
    this % order       = [a % order,       order]
    this % along       = [a % along,       FIRST_COORDINATE]
    this % mixed_order = [a % mixed_order, 0]
    this % coefficient = [a % coefficient, coefficient]

    call this % derive_domain()

  end function applied

  !===================================================================!
  ! Returns the reads of a vertex table moved by an offset, a read of
  ! none left zero
  !===================================================================!

  pure function shifted(reads, by) result(moved)

    integer, intent(in) :: reads(:), by
    integer             :: moved(size(reads))

    where (reads > 0)
       moved = reads + by
    elsewhere
       moved = 0
    end where

  end function shifted

  !===================================================================!
  ! Sets the domain, read from the leaves: the fields are those any
  ! unknown names, one at least; a field's degree along the first
  ! coordinate is the highest order named for it, zero when only its
  ! value is named; along a later coordinate the degree is the highest
  ! order named there by any field; the multipliers are those any
  ! multiplier leaf names. The argument contracts follow, with the
  ! highest exact degree, one less for a stationarity.
  !===================================================================!

  subroutine derive_domain(this, label)

    class(expression), intent(inout)        :: this
    character(len=*) , intent(in), optional :: label

    integer              :: i, f, c, coordinates
    integer, allocatable :: along_degree(:)

    this % fields      = 1
    this % multipliers = 0
    this % designs     = 1
    this % position_coordinates = 0
    coordinates        = FIRST_COORDINATE
    do i = 1, this % num_vertices()
       if (this % kind(i) /= VERTEX_LEAF) cycle
       select case (this % position(i))
       case (ARGUMENT_STATE)
          this % fields = max(this % fields, this % field(i))
          coordinates   = max(coordinates, this % along(i))
       case (ARGUMENT_MULTIPLIER)
          this % multipliers = max(this % multipliers, this % field(i))
       case (ARGUMENT_DESIGN)
          this % designs = max(this % designs, this % field(i))
       case (ARGUMENT_COORDINATE)
          this % position_coordinates = max(this % position_coordinates, this % field(i))
       end select
    end do

    ! in one pass over the leaves: the highest mixed order of each
    ! field, its highest order along the first coordinate, and the
    ! highest order along each later coordinate
    this % mixed_degree_of_field = [(0, f = 1, this % fields)]
    this % field_degree = [(-1, f = 1, this % fields)]
    allocate(along_degree(coordinates), source=-1)
    do i = 1, this % num_vertices()
       if (this % kind(i) /= VERTEX_LEAF .or. this % position(i) /= ARGUMENT_STATE) cycle
       f = this % field(i)
       this % mixed_degree_of_field(f) = max(this % mixed_degree_of_field(f), this % mixed_order(i))
       if (this % along(i) == FIRST_COORDINATE) this % field_degree(f) = max(this % field_degree(f), this % order(i))
       if (this % along(i) >= 1) along_degree(this % along(i)) = max(along_degree(this % along(i)), this % order(i))
    end do
    ! the orders along the first coordinate a mixed component is read
    ! from are stored too
    this % field_degree = [(max(0, this % field_degree(f), this % mixed_degree_of_field(f)), f = 1, this % fields)]
    this % degrees      = [this % field_degree(1), (max(0, along_degree(c)), c = FIRST_COORDINATE + 1, coordinates)]

    call this % declare_domain(label)

  end subroutine derive_domain

  !===================================================================!
  ! Sets the layout of the tuple and the argument contracts of the
  ! domain: the state, the design, and for a rule that reads a
  ! coordinate the positions of the points as a third argument, one
  ! value per coordinate read per point
  !===================================================================!

  subroutine declare_domain(this, label)

    class(expression), intent(inout)        :: this
    character(len=*) , intent(in), optional :: label

    type(contract), allocatable :: contracts(:)
    integer, allocatable        :: offset(:)
    integer                     :: f, i

    allocate(offset(this % fields + 1))
    offset(1) = 0
    do f = 1, this % fields
       offset(f + 1) = offset(f) + this % components_per_field(f)
    end do
    this % field_offset = offset
    this % leaf_component = [(-1, i = 1, this % num_vertices())]
    do i = 1, this % num_vertices()
       if (this % kind(i) /= VERTEX_LEAF) cycle
       if (this % position(i) /= ARGUMENT_STATE .and. this % position(i) /= ARGUMENT_MULTIPLIER) cycle
       this % leaf_component(i) = this % stored_at(this % position(i), this % field(i), this % along(i), &
            & this % order(i), this % mixed_order(i))
    end do
    contracts = [contract(FIELD_REAL, this % num_components()), contract(FIELD_REAL, this % designs)]
    if (this % position_coordinates > 0) contracts = [contracts, contract(FIELD_REAL, this % position_coordinates)]
    call this % declare_arguments(size(contracts), contracts, label=label, &
         & max_degree=max_subset_width() - merge(1, 0, this % varied > 0))

  end subroutine declare_domain

  !===================================================================!
  ! Returns the expression on the domain of another that contains its
  ! own: x read from the tuple, the designs and the positions of y, so
  ! that two rules evaluated on one block index one tuple. Invalid
  ! input: a domain that does not contain the expression's.
  !===================================================================!

  function on_domain(x, y) result(this)

    type(expression), intent(in) :: x, y
    type(expression)             :: this

    integer :: f, c

    if (x % fields > y % fields .or. x % multipliers > y % multipliers .or. x % designs > y % designs &
         & .or. x % position_coordinates > y % position_coordinates .or. size(x % degrees) > size(y % degrees)) then
       error stop 'operation_expression: the domain does not contain the expression''s own'
    end if
    do f = 1, x % fields
       if (x % field_degree(f) > y % field_degree(f) .or. x % mixed_degree_of_field(f) > y % mixed_degree_of_field(f)) &
            & error stop 'operation_expression: the domain does not contain the expression''s own'
    end do
    do c = 1, size(x % degrees)
       if (x % degrees(c) > y % degrees(c)) error stop 'operation_expression: the domain does not contain the &
            &expression''s own'
    end do
    this = x
    this % fields                = y % fields
    this % multipliers           = y % multipliers
    this % designs               = y % designs
    this % position_coordinates  = y % position_coordinates
    this % field_degree          = y % field_degree
    this % mixed_degree_of_field = y % mixed_degree_of_field
    this % degrees               = y % degrees
    call this % declare_domain()

  end function on_domain

  !===================================================================!
  ! Returns a + b
  !===================================================================!

  function plus(a, b) result(this)

    type(expression), intent(in) :: a, b
    type(expression)             :: this

    this = joined(a, b, VERTEX_SUM)

  end function plus

  !===================================================================!
  ! Returns x + b
  !===================================================================!

  function real_plus(x, b) result(this)

    real(dp)        , intent(in) :: x
    type(expression), intent(in) :: b
    type(expression)             :: this

    this = joined(constant(x), b, VERTEX_SUM)

  end function real_plus

  !===================================================================!
  ! Returns a + x
  !===================================================================!

  function plus_real(a, x) result(this)

    type(expression), intent(in) :: a
    real(dp)        , intent(in) :: x
    type(expression)             :: this

    this = joined(a, constant(x), VERTEX_SUM)

  end function plus_real

  !===================================================================!
  ! Returns a - b
  !===================================================================!

  function minus(a, b) result(this)

    type(expression), intent(in) :: a, b
    type(expression)             :: this

    this = joined(a, b, VERTEX_DIFFERENCE)

  end function minus

  !===================================================================!
  ! Returns x - b
  !===================================================================!

  function real_minus(x, b) result(this)

    real(dp)        , intent(in) :: x
    type(expression), intent(in) :: b
    type(expression)             :: this

    this = joined(constant(x), b, VERTEX_DIFFERENCE)

  end function real_minus

  !===================================================================!
  ! Returns a - x
  !===================================================================!

  function minus_real(a, x) result(this)

    type(expression), intent(in) :: a
    real(dp)        , intent(in) :: x
    type(expression)             :: this

    this = joined(a, constant(x), VERTEX_DIFFERENCE)

  end function minus_real

  !===================================================================!
  ! Returns -a, as (-1) a
  !===================================================================!

  function negated(a) result(this)

    type(expression), intent(in) :: a
    type(expression)             :: this

    this = joined(constant(-1.0_dp), a, VERTEX_PRODUCT)

  end function negated

  !===================================================================!
  ! Returns a * b
  !===================================================================!

  function times(a, b) result(this)

    type(expression), intent(in) :: a, b
    type(expression)             :: this

    this = joined(a, b, VERTEX_PRODUCT)

  end function times

  !===================================================================!
  ! Returns x * b
  !===================================================================!

  function real_times(x, b) result(this)

    real(dp)        , intent(in) :: x
    type(expression), intent(in) :: b
    type(expression)             :: this

    this = joined(constant(x), b, VERTEX_PRODUCT)

  end function real_times

  !===================================================================!
  ! Returns a * x
  !===================================================================!

  function times_real(a, x) result(this)

    type(expression), intent(in) :: a
    real(dp)        , intent(in) :: x
    type(expression)             :: this

    this = joined(a, constant(x), VERTEX_PRODUCT)

  end function times_real

  !===================================================================!
  ! Returns a / b
  !===================================================================!

  function over(a, b) result(this)

    type(expression), intent(in) :: a, b
    type(expression)             :: this

    this = joined(a, b, VERTEX_QUOTIENT)

  end function over

  !===================================================================!
  ! Returns x / b
  !===================================================================!

  function real_over(x, b) result(this)

    real(dp)        , intent(in) :: x
    type(expression), intent(in) :: b
    type(expression)             :: this

    this = joined(constant(x), b, VERTEX_QUOTIENT)

  end function real_over

  !===================================================================!
  ! Returns a / x
  !===================================================================!

  function over_real(a, x) result(this)

    type(expression), intent(in) :: a
    real(dp)        , intent(in) :: x
    type(expression)             :: this

    this = joined(a, constant(x), VERTEX_QUOTIENT)

  end function over_real

  !===================================================================!
  ! Returns a ** n
  !===================================================================!

  function to_integer(a, n) result(this)

    type(expression), intent(in) :: a
    integer         , intent(in) :: n
    type(expression)             :: this

    this = applied(a, VERTEX_INTEGER_POWER, 0, real(n, dp))

  end function to_integer

  !===================================================================!
  ! Returns a ** p
  !===================================================================!

  function to_real(a, p) result(this)

    type(expression), intent(in) :: a
    real(dp)        , intent(in) :: p
    type(expression)             :: this

    this = applied(a, VERTEX_REAL_POWER, 0, p)

  end function to_real

  !===================================================================!
  ! Returns sin(a)
  !===================================================================!

  function sine_of(a) result(this)

    type(expression), intent(in) :: a
    type(expression)             :: this

    this = applied(a, VERTEX_FUNCTION, SINE, 0.0_dp)

  end function sine_of

  !===================================================================!
  ! Returns cos(a)
  !===================================================================!

  function cosine_of(a) result(this)

    type(expression), intent(in) :: a
    type(expression)             :: this

    this = applied(a, VERTEX_FUNCTION, COSINE, 0.0_dp)

  end function cosine_of

  !===================================================================!
  ! Returns exp(a)
  !===================================================================!

  function exponential_of(a) result(this)

    type(expression), intent(in) :: a
    type(expression)             :: this

    this = applied(a, VERTEX_FUNCTION, EXPONENTIAL, 0.0_dp)

  end function exponential_of

  !===================================================================!
  ! Returns log(a)
  !===================================================================!

  function logarithm_of(a) result(this)

    type(expression), intent(in) :: a
    type(expression)             :: this

    this = applied(a, VERTEX_FUNCTION, LOGARITHM, 0.0_dp)

  end function logarithm_of

  !===================================================================!
  ! Returns sqrt(a)
  !===================================================================!

  function root_of(a) result(this)

    type(expression), intent(in) :: a
    type(expression)             :: this

    this = applied(a, VERTEX_FUNCTION, SQUARE_ROOT, 0.0_dp)

  end function root_of

  !===================================================================!
  ! Returns the step H(a); a step of anything but the coordinates
  ! stops the program
  !===================================================================!

  function step_of(a) result(this)

    type(expression), intent(in) :: a
    type(expression)             :: this

    integer :: i

    do i = 1, a % num_vertices()
       if (a % kind(i) == VERTEX_LEAF .and. a % position(i) /= ARGUMENT_COORDINATE) then
          error stop 'operation_expression: a step is a function of the coordinates alone'
       end if
    end do
    this = applied(a, VERTEX_FUNCTION, UNIT_STEP, 0.0_dp)

  end function step_of

  !===================================================================!
  ! Returns the rule at one point over derivative terms: the loop over
  ! the vertices. The tuple evaluated stores the state and then the
  ! multipliers; a multiplier, read from no input, is zero, and the
  ! varied one has one more direction, the highest, whose partial is
  ! returned over the caller's directions. A state component past
  ! those given, or a function index outside those defined, stops the
  ! program.
  !===================================================================!

  pure function expression_at_instant(this, q, nu, position) result(r)

    class(expression)     , intent(in)           :: this
    type(derivative_terms), intent(in)           :: q(0:)
    type(derivative_terms), intent(in)           :: nu(:)
    real(dp)              , intent(in), optional :: position(:)
    type(derivative_terms)                       :: r

    integer :: n

    if (size(q) /= this % num_components()) then
       error stop 'operation_expression: the point does not store one component per law component'
    end if
    if (size(nu) < this % designs) then
       error stop 'operation_expression: the point does not store one value per design the law reads'
    end if

    ! the tuple evaluated stores the state and then the multipliers;
    ! a multiplier, read from no input, is zero, and the varied one
    ! has one more direction, the highest, whose partial is returned
    ! over the caller's
    n = nu(1) % num_directions()
    if (size(q) > 0) n = max(n, q(0) % num_directions())
    if (this % multipliers > 0 .and. this % varied > 0) n = n + 1
    r = this % evaluated_over(q, nu, n, position)
    if (this % multipliers > 0 .and. this % varied > 0) r = partial(r, n)

  end function expression_at_instant

  !===================================================================!
  ! Returns the value at a position of an expression that reads no
  ! state and no multiplier: a function of the coordinates and of
  ! constants alone. Invalid input: an expression with a state or
  ! multiplier leaf, which has no value until the state is given.
  !===================================================================!

  pure real(dp) function value_at(this, position)

    class(expression), intent(in) :: this
    real(dp)         , intent(in) :: position(:)

    real(dp), allocatable :: q(:), v(:), dv(:)

    if (this % reads_state()) then
       error stop 'operation_expression: value_at requires an expression of the coordinates alone, &
            &but the expression reads the state or a multiplier'
    end if
    allocate(q(0:this % num_components() - 1), source=0.0_dp)
    call this % values_over(q, [0.0_dp], position, v, dv)
    value_at = v(size(v))

  end function value_at

  !===================================================================!
  ! Returns the values at every vertex, from the leaves to the root,
  ! for one state q, one design nu and, when a coordinate is read, one
  ! position: the evaluation of the graph in real arithmetic, with the
  ! derivative dv of every value in the varied multiplier (zero
  ! throughout when none is varied). The multipliers are zero, as in
  ! at_instant; the varied one has derivative one. The derivative
  ! arithmetic of evaluated_over reads the same leaves.
  !===================================================================!

  pure subroutine values_over(this, q, nu, position, v, dv)

    class(expression)    , intent(in)           :: this
    real(dp)             , intent(in)           :: q(0:)
    real(dp)             , intent(in)           :: nu(:)
    real(dp)             , intent(in), optional :: position(:)
    real(dp), allocatable, intent(out)          :: v(:), dv(:)

    integer  :: i, at, f, s, n
    real(dp) :: d

    if (size(nu) < this % designs) then
       error stop 'operation_expression: the point does not store one value per design the law reads'
    end if
    allocate(v(this % num_vertices()), dv(this % num_vertices()))
    dv = 0.0_dp

    do i = 1, this % num_vertices()
       f = this % first(i)
       s = this % second(i)
       select case (this % kind(i))
       case (VERTEX_LEAF)
          if (this % position(i) == ARGUMENT_DESIGN) then
             v(i) = nu(this % field(i))
          else if (this % position(i) == ARGUMENT_COORDINATE) then
             if (.not. present(position)) then
                error stop 'operation_expression: the expression reads a coordinate, but no position was given'
             end if
             if (this % field(i) > size(position)) then
                error stop 'operation_expression: the expression reads a coordinate the position does not have'
             end if
             v(i) = position(this % field(i))
          else if (this % position(i) == ARGUMENT_MULTIPLIER) then
             v(i) = 0.0_dp
             if (this % field(i) == this % varied) dv(i) = 1.0_dp
          else
             at = this % leaf_at(i)
             if (at > ubound(q, 1)) then
                error stop 'operation_expression: the state does not store the component read, &
                     &beyond the end of the stored tuple'
             end if
             v(i) = q(at)
          end if
       case (VERTEX_CONSTANT)
          v(i) = this % coefficient(i)
       case (VERTEX_SUM)
          v(i)  = v(f) + v(s)
          dv(i) = dv(f) + dv(s)
       case (VERTEX_DIFFERENCE)
          v(i)  = v(f) - v(s)
          dv(i) = dv(f) - dv(s)
       case (VERTEX_PRODUCT)
          v(i)  = v(f) * v(s)
          dv(i) = dv(f) * v(s) + v(f) * dv(s)
       case (VERTEX_QUOTIENT)
          v(i)  = v(f) / v(s)
          dv(i) = (dv(f) - v(i) * dv(s)) / v(s)
       case (VERTEX_INTEGER_POWER)
          n = nint(this % coefficient(i))
          v(i) = v(f) ** n
          d = 0.0_dp
          if (n /= 0) d = n * v(f) ** (n - 1)
          dv(i) = d * dv(f)
       case (VERTEX_REAL_POWER)
          v(i)  = v(f) ** this % coefficient(i)
          dv(i) = this % coefficient(i) * v(f) ** (this % coefficient(i) - 1.0_dp) * dv(f)
       case (VERTEX_FUNCTION)
          select case (this % order(i))
          case (SINE);        v(i) = sin(v(f));  dv(i) =  cos(v(f)) * dv(f)
          case (COSINE);      v(i) = cos(v(f));  dv(i) = -sin(v(f)) * dv(f)
          case (EXPONENTIAL); v(i) = exp(v(f));  dv(i) =  v(i) * dv(f)
          case (LOGARITHM);   v(i) = log(v(f));  dv(i) =  dv(f) / v(f)
          case (SQUARE_ROOT); v(i) = sqrt(v(f)); dv(i) =  dv(f) / (2.0_dp * v(i))
          case (UNIT_STEP)
             if (dv(f) /= 0.0_dp) error stop 'operation_expression: a step is differentiated along a direction &
                  &its argument reads, where its derivative is a measure'
             v(i)  = merge(1.0_dp, 0.0_dp, v(f) >= 0.0_dp)
             dv(i) = 0.0_dp
          case default
             error stop 'operation_expression: this vertex names a function that is not one of &
                  &those defined'
          end select
       case (VERTEX_DERIVATIVE)
          error stop 'operation_expression: a law with a derivative operator is evaluated only in its &
               &jet form; jet_form replaces the operator by the jet components it denotes'
       case default
          error stop 'operation_expression: this vertex names a kind that is not one of those &
               &defined'
       end select
    end do

  end subroutine values_over

  !===================================================================!
  ! Returns the rule's value on the real values of the state and the
  ! design, by the forward pass of gradient_at alone
  !===================================================================!

  pure real(dp) function value_on(this, q, nu, position) result(value)

    class(expression), intent(in)           :: this
    real(dp)         , intent(in)           :: q(0:)
    real(dp)         , intent(in)           :: nu(:)
    real(dp)         , intent(in), optional :: position(:)

    real(dp), allocatable :: v(:), dv(:)

    call this % values_over(q, nu, position, v, dv)
    if (this % varied > 0) then
       value = dv(size(v))
    else
       value = v(size(v))
    end if

  end function value_on

  !===================================================================!
  ! Returns the rule's value and its gradient in the state by one
  ! reverse pass: for a rule that is the stationarity in multiplier j,
  ! the value is d(root)/d(lambda_j) and g(k) = d^2(root)/d(lambda_j)
  ! d(q(k)); for a rule without a varied multiplier, the value is the
  ! root and g(k) = d(root)/d(q(k)). The gradient in the design is
  ! formed by the same pass when requested.
  !
  !      forward   v(i), dv(i) = dv(i)/dlambda_j      leaves to root
  !      reverse   a(i)  = d(root)/dv(i)              root to leaves
  !                da(i) = da(i)/dlambda_j
  !      leaf q(k) g(k) += da(i) for a stationarity, a(i) otherwise
  !
  ! The forward pass has every value with its derivative in lambda_j;
  ! the adjoint of each vertex is formed the same way from the root
  ! down, so that the derivative part of the adjoint at a state leaf is
  ! the mixed partial. The vertices are in evaluation order with the
  ! root last, so the pass from the last vertex to the first visits
  ! every vertex after all that read it. One pass gives every partial
  ! the jacobian of a rule needs, where the derivative arithmetic gives
  ! one direction per pass.
  !===================================================================!

  pure subroutine gradient_at(this, q, nu, value, g, position, g_design)

    class(expression), intent(in)            :: this
    real(dp)         , intent(in)            :: q(0:)
    real(dp)         , intent(in)            :: nu(:)
    real(dp)         , intent(out)           :: value
    real(dp)         , intent(out)           :: g(0:)
    real(dp)         , intent(in), optional  :: position(:)
    real(dp)         , intent(out), optional :: g_design(:)

    real(dp), allocatable :: v(:), dv(:), a(:), da(:)
    integer               :: i, f, s, at, n
    real(dp)              :: d, dd, c, dc

    call this % values_over(q, nu, position, v, dv)
    if (this % varied > 0) then
       value = dv(size(v))
    else
       value = v(size(v))
    end if

    allocate(a(size(v)), da(size(v)), source=0.0_dp)
    a(size(v)) = 1.0_dp
    g = 0.0_dp
    if (present(g_design)) g_design = 0.0_dp

    do i = this % num_vertices(), 1, -1
       if (a(i) == 0.0_dp .and. da(i) == 0.0_dp) cycle
       f = this % first(i)
       s = this % second(i)
       select case (this % kind(i))
       case (VERTEX_LEAF)
          if (this % position(i) == ARGUMENT_STATE) then
             at = this % leaf_at(i)
             if (this % varied > 0) then
                g(at) = g(at) + da(i)
             else
                g(at) = g(at) + a(i)
             end if
          else if (this % position(i) == ARGUMENT_DESIGN .and. present(g_design)) then
             ! the gradient in the design, by the same reverse pass
             if (this % varied > 0) then
                g_design(this % field(i)) = g_design(this % field(i)) + da(i)
             else
                g_design(this % field(i)) = g_design(this % field(i)) + a(i)
             end if
          end if
       case (VERTEX_CONSTANT)
          continue
       case (VERTEX_SUM)
          a(f)  = a(f)  + a(i)
          a(s)  = a(s)  + a(i)
          da(f) = da(f) + da(i)
          da(s) = da(s) + da(i)
       case (VERTEX_DIFFERENCE)
          a(f)  = a(f)  + a(i)
          a(s)  = a(s)  - a(i)
          da(f) = da(f) + da(i)
          da(s) = da(s) - da(i)
       case (VERTEX_PRODUCT)
          a(f)  = a(f)  + a(i) * v(s)
          da(f) = da(f) + da(i) * v(s) + a(i) * dv(s)
          a(s)  = a(s)  + a(i) * v(f)
          da(s) = da(s) + da(i) * v(f) + a(i) * dv(f)
       case (VERTEX_QUOTIENT)
          ! d(root)/d(first) = a / second, d(root)/d(second) = -a v / second
          c  = 1.0_dp / v(s)
          dc = -dv(s) * c * c
          a(f)  = a(f)  + a(i) * c
          da(f) = da(f) + da(i) * c + a(i) * dc
          a(s)  = a(s)  - a(i) * v(i) * c
          da(s) = da(s) - (da(i) * v(i) * c + a(i) * dv(i) * c + a(i) * v(i) * dc)
       case (VERTEX_INTEGER_POWER)
          n = nint(this % coefficient(i))
          d  = 0.0_dp
          dd = 0.0_dp
          if (n /= 0) d = n * v(f) ** (n - 1)
          if (n /= 0 .and. n /= 1) dd = n * (n - 1) * v(f) ** (n - 2) * dv(f)
          a(f)  = a(f)  + a(i) * d
          da(f) = da(f) + da(i) * d + a(i) * dd
       case (VERTEX_REAL_POWER)
          d  = this % coefficient(i) * v(f) ** (this % coefficient(i) - 1.0_dp)
          dd = this % coefficient(i) * (this % coefficient(i) - 1.0_dp) * v(f) ** (this % coefficient(i) - 2.0_dp) * dv(f)
          a(f)  = a(f)  + a(i) * d
          da(f) = da(f) + da(i) * d + a(i) * dd
       case (VERTEX_FUNCTION)
          select case (this % order(i))
          case (SINE);        d =  cos(v(f));          dd = -sin(v(f)) * dv(f)
          case (COSINE);      d = -sin(v(f));          dd = -cos(v(f)) * dv(f)
          case (EXPONENTIAL); d =  v(i);               dd =  dv(i)
          case (LOGARITHM);   d =  1.0_dp / v(f);      dd = -dv(f) / (v(f) * v(f))
          case (SQUARE_ROOT); d =  0.5_dp / v(i);      dd = -0.5_dp * dv(i) / (v(i) * v(i))
          case default
             ! the step, of the coordinates alone
             d  = 0.0_dp
             dd = 0.0_dp
          end select
          a(f)  = a(f)  + a(i) * d
          da(f) = da(f) + da(i) * d + a(i) * dd
       end select
    end do

  end subroutine gradient_at

  !===================================================================!
  ! Returns whether any leaf reads the state or a multiplier
  !===================================================================!

  pure logical function reads_state(this)

    class(expression), intent(in) :: this

    integer :: i

    reads_state = .false.
    do i = 1, this % num_vertices()
       if (this % kind(i) /= VERTEX_LEAF) cycle
       if (this % position(i) == ARGUMENT_STATE .or. this % position(i) == ARGUMENT_MULTIPLIER) then
          reads_state = .true.
          return
       end if
    end do

  end function reads_state

  !===================================================================!
  ! Returns whether any leaf reads the time map
  !===================================================================!

  pure logical function reads_time_map(this)

    class(expression), intent(in) :: this

    integer :: i

    reads_time_map = .false.
    do i = 1, this % num_vertices()
       if (this % kind(i) == VERTEX_LEAF .and. this % position(i) == ARGUMENT_TIME_MAP) reads_time_map = .true.
    end do

  end function reads_time_map

  !===================================================================!
  ! Returns whether the expression is t(tau) - tau, the time map less
  ! the coordinate of time: the pin of the instants
  !
  !      vertex   1 t(tau)   2 tau   3 difference(1, 2)
  !===================================================================!

  pure logical function pins_instants(this)

    class(expression), intent(in) :: this

    integer :: r, a, b

    pins_instants = .false.
    r = this % root()
    if (r /= 3) return
    if (this % kind(r) /= VERTEX_DIFFERENCE) return
    a = this % first(r)
    b = this % second(r)
    if (this % kind(a) /= VERTEX_LEAF .or. this % kind(b) /= VERTEX_LEAF) return
    pins_instants = this % position(a) == ARGUMENT_TIME_MAP .and. this % position(b) == ARGUMENT_COORDINATE &
         & .and. this % field(b) == FIRST_COORDINATE

  end function pins_instants

  !===================================================================!
  ! Returns whether the zero state lies in the expression's domain. A
  ! vertex reads the state when it is a state component or when either
  ! of its operands does. The expression has a value at the zero state
  ! unless such a vertex is
  !
  !      a divisor
  !      the base of a negative integer power
  !      the base of a real power of exponent below one
  !      the argument of a logarithm or a square root
  !
  ! each of which is singular, or has a singular derivative, at zero.
  ! The design is an input of its own and is not zeroed with the state.
  !===================================================================!

  pure logical function expression_defined_at_zero(this) result(defined)

    class(expression), intent(in) :: this

    logical, allocatable :: reads_state(:)
    integer              :: i

    allocate(reads_state(this % num_vertices()), source=.false.)
    defined = .true.

    do i = 1, this % num_vertices()
       select case (this % kind(i))
       case (VERTEX_LEAF)
          reads_state(i) = this % position(i) == ARGUMENT_STATE
       case (VERTEX_CONSTANT)
          reads_state(i) = .false.
       case (VERTEX_SUM, VERTEX_DIFFERENCE, VERTEX_PRODUCT)
          reads_state(i) = reads_state(this % first(i)) .or. reads_state(this % second(i))
       case (VERTEX_QUOTIENT)
          reads_state(i) = reads_state(this % first(i)) .or. reads_state(this % second(i))
          if (reads_state(this % second(i))) defined = .false.
       case (VERTEX_INTEGER_POWER)
          reads_state(i) = reads_state(this % first(i))
          if (reads_state(i) .and. nint(this % coefficient(i)) < 0) defined = .false.
       case (VERTEX_REAL_POWER)
          reads_state(i) = reads_state(this % first(i))
          if (reads_state(i) .and. this % coefficient(i) < 1.0_dp) defined = .false.
       case (VERTEX_FUNCTION)
          reads_state(i) = reads_state(this % first(i))
          if (reads_state(i) .and. (this % order(i) == LOGARITHM .or. this % order(i) == SQUARE_ROOT)) then
             defined = .false.
          end if
       case default
          defined = .false.
       end select
    end do

  end function expression_defined_at_zero

  !===================================================================!
  ! Returns the rule over n directions by the loop over the vertices,
  ! on the stored tuple: the coefficients of every vertex over the n
  ! directions in one array, each rule of util_derivative_terms
  ! applied to the coefficients of its operands. A state leaf reads its
  ! component of q and a design leaf its design, the directions either
  ! lacks zero; a multiplier is zero, the varied one with the
  ! derivative one along the direction n.
  !
  !      v(s, i)         the coefficient of vertex i on the subset s of
  !                      the directions, a set of bits
  !      v(0, i)         the value
  !      v(2**(k-1), i)  the partial along direction k alone
  !      v(2**n - 1, i)  the mixed partial along all n
  !===================================================================!

  pure function evaluated_over(this, q, nu, n, position) result(r)

    class(expression)     , intent(in)           :: this
    type(derivative_terms), intent(in)           :: q(0:)
    type(derivative_terms), intent(in)           :: nu(:)
    integer               , intent(in)           :: n
    real(dp)              , intent(in), optional :: position(:)
    type(derivative_terms)                       :: r

    real(dp), allocatable :: v(:,:)
    integer               :: i, at, f, s

    allocate(v(0:2**n - 1, this % num_vertices()))

    do i = 1, this % num_vertices()
       f = this % first(i)
       s = this % second(i)
       select case (this % kind(i))
       case (VERTEX_LEAF)
          select case (this % position(i))
          case (ARGUMENT_DESIGN)
             call extended_terms(nu(this % field(i)), v(:, i))
          case (ARGUMENT_COORDINATE)
             ! a coordinate is read from the position of the point,
             ! which an evaluation without a position does not have
             if (.not. present(position)) then
                error stop 'operation_expression: the expression reads a coordinate, but no position was given'
             end if
             if (this % field(i) > size(position)) then
                error stop 'operation_expression: the expression reads a coordinate the position does not have'
             end if
             v(:, i) = 0.0_dp
             v(0, i) = position(this % field(i))
          case (ARGUMENT_MULTIPLIER)
             v(:, i) = 0.0_dp
             if (this % field(i) == this % varied) v(2**(n - 1), i) = 1.0_dp
          case default
             ! the tuple lists every field's components together, the
             ! components along time first, then along space
             at = this % leaf_at(i)
             if (at > ubound(q, 1)) then
                error stop 'operation_expression: the state does not store the component read, &
                     &beyond the end of the stored tuple'
             end if
             call extended_terms(q(at), v(:, i))
          end select
       case (VERTEX_CONSTANT)
          v(:, i) = 0.0_dp
          v(0, i) = this % coefficient(i)
       case (VERTEX_SUM)
          v(:, i) = v(:, f) + v(:, s)
       case (VERTEX_DIFFERENCE)
          v(:, i) = v(:, f) - v(:, s)
       case (VERTEX_PRODUCT)
          call product_terms(v(:, f), v(:, s), v(:, i))
       case (VERTEX_QUOTIENT)
          call quotient_terms(v(:, f), v(:, s), v(:, i))
       case (VERTEX_INTEGER_POWER)
          call integer_power_terms(v(:, f), nint(this % coefficient(i)), v(:, i))
       case (VERTEX_REAL_POWER)
          call real_power_terms(v(:, f), this % coefficient(i), v(:, i))
       case (VERTEX_FUNCTION)
          select case (this % order(i))
          case (SINE);        call sine_terms(v(:, f), v(:, i))
          case (COSINE);      call cosine_terms(v(:, f), v(:, i))
          case (EXPONENTIAL); call exponential_terms(v(:, f), v(:, i))
          case (LOGARITHM);   call logarithm_terms(v(:, f), v(:, i))
          case (SQUARE_ROOT); call real_power_terms(v(:, f), 0.5_dp, v(:, i))
          case (UNIT_STEP)
             if (any(v(1:, f) /= 0.0_dp)) error stop 'operation_expression: a step is differentiated along a &
                  &direction its argument reads, where its derivative is a measure'
             v(:, i) = 0.0_dp
             v(0, i) = merge(1.0_dp, 0.0_dp, v(0, f) >= 0.0_dp)
          case default
             error stop 'operation_expression: this vertex names a function that is not one of &
                  &those defined'
          end select
       case (VERTEX_DERIVATIVE)
          error stop 'operation_expression: a law with a derivative operator is evaluated only in its &
               &jet form; jet_form replaces the operator by the jet components it denotes'
       case default
          error stop 'operation_expression: this vertex names a kind that is not one of those &
               &defined'
       end select
    end do

    r = derivative_terms(v(:, this % num_vertices()))

  end function evaluated_over

  !===================================================================!
  ! Returns the highest order this rule names along one coordinate,
  ! and minus one when it names none along it. One query serves every
  ! coordinate, so adding a coordinate requires no new procedure.
  !===================================================================!

  pure integer function highest_degree_along(this, coordinate, field)

    class(expression), intent(in)           :: this
    integer          , intent(in)           :: coordinate
    integer          , intent(in), optional :: field

    integer :: i

    highest_degree_along = -1
    do i = 1, this % num_vertices()
       if (this % kind(i) == VERTEX_LEAF .and. this % position(i) == ARGUMENT_STATE &
            & .and. this % along(i) == coordinate) then
          if (present(field)) then
             if (this % field(i) /= field) cycle
          end if
          highest_degree_along = max(highest_degree_along, this % order(i))
       end if
    end do

  end function highest_degree_along

  !===================================================================!
  ! Returns the highest order along the first coordinate at which this
  ! rule reads a field (a component along a later coordinate at the
  ! order along the first it is mixed with), and minus one when it
  ! reads none of the field
  !===================================================================!

  pure integer function highest_order_read(this, field)

    class(expression), intent(in) :: this
    integer          , intent(in) :: field

    integer :: i

    highest_order_read = -1
    do i = 1, this % num_vertices()
       if (this % kind(i) /= VERTEX_LEAF .or. this % position(i) /= ARGUMENT_STATE) cycle
       if (this % field(i) /= field) cycle
       if (this % along(i) == FIRST_COORDINATE) then
          highest_order_read = max(highest_order_read, this % order(i))
       else
          highest_order_read = max(highest_order_read, this % mixed_order(i))
       end if
    end do

  end function highest_order_read

  !===================================================================!
  ! Returns the components of the tuple the rule reads, each once, in
  ! increasing order: the state leaves' components, a multiplier's
  ! excluded since the tuple does not store it
  !===================================================================!

  subroutine read_components(this, components)

    class(expression)   , intent(in)  :: this
    integer, allocatable, intent(out) :: components(:)

    logical, allocatable :: read(:)
    integer              :: i, at

    allocate(read(0:this % num_components() - 1), source=.false.)
    do i = 1, this % num_vertices()
       if (this % kind(i) /= VERTEX_LEAF .or. this % position(i) /= ARGUMENT_STATE) cycle
       at = this % leaf_at(i)
       read(at) = .true.
    end do
    components = pack([(at, at = 0, size(read) - 1)], read)

  end subroutine read_components

  !===================================================================!
  ! Returns the count of vertices. Invalid input: a variable of type
  ! expression never assigned from unknown, design, constant or an
  ! operator has no vertex, so no root and no value. Every read of the
  ! vertex count passes through here, so this is the one refusal;
  ! without it the compiled code reads an array with no memory and
  ! faults in memmove with no message.
  !===================================================================!

  pure integer function num_vertices(this)

    class(expression), intent(in) :: this

    if (.not. allocated(this % kind)) then
       error stop 'operation_expression: the expression has no vertex, so no root and no value; &
            &a variable of type expression never assigned from unknown(), design(), constant() &
            &or an operator has no vertex'
    end if
    num_vertices = size(this % kind)

  end function num_vertices

  !===================================================================!
  ! Returns the root: the vertices are stored in evaluation order, so
  ! the vertex whose value is the value of the expression is the last
  !===================================================================!

  pure integer function root(this)

    class(expression), intent(in) :: this

    root = this % num_vertices()

  end function root

  !===================================================================!
  ! Returns the vertex table as values, in evaluation order, for a
  ! view
  !===================================================================!

  function vertices(this) result(v)

    class(expression), intent(in)        :: this
    type(expression_vertex), allocatable :: v(:)

    integer :: i

    allocate(v(this % num_vertices()))
    do i = 1, size(v)
       v(i) = expression_vertex(this % kind(i), this % first(i), this % second(i), this % position(i), &
            & this % field(i), this % order(i), this % along(i), this % mixed_order(i), this % coefficient(i))
    end do

  end function vertices

  !===================================================================!
  ! Returns the component a coordinate and an order name. The first
  ! coordinate's orders run from zero and every later one's from one,
  ! since order zero is the state itself and is named once. A mixed
  ! order names the component of that order along the first coordinate
  ! too.
  !===================================================================!

  pure integer function component_at(this, coordinate, order, field, mixed_order) result(at)

    class(expression), intent(in)           :: this
    integer          , intent(in)           :: coordinate, order
    integer          , intent(in), optional :: field, mixed_order

    integer :: f

    f = 1
    if (present(field)) f = field
    at = this % offset_of_field(f) + this % order_at(f, coordinate, order, mixed_order)

  end function component_at

  !===================================================================!
  ! Returns the component of the tuple a state or multiplier leaf
  ! reads: the one formed with the domain, where it was
  !===================================================================!

  pure integer function leaf_at(this, i) result(at)

    class(expression), intent(in) :: this
    integer          , intent(in) :: i

    if (allocated(this % leaf_component)) then
       if (i <= size(this % leaf_component)) then
          at = this % leaf_component(i)
          if (at >= 0) return
       end if
    end if
    at = this % stored_at(this % position(i), this % field(i), this % along(i), this % order(i), &
         & this % mixed_order(i))

  end function leaf_at

  !===================================================================!
  ! Returns the first component of a field in the tuple the rule
  ! reads, which stores the state fields alone. A multiplier stops the
  ! program.
  !===================================================================!

  pure integer function offset_of_field(this, field) result(at)

    class(expression), intent(in) :: this
    integer          , intent(in) :: field

    if (field < 1 .or. field > this % fields + this % multipliers) then
       error stop 'operation_expression: this component names a field the rule does not read'
    end if
    if (field > this % fields) then
       error stop 'operation_expression: this component names a multiplier, which is not stored &
            &in the tuple'
    end if
    at = this % field_offset(field)

  end function offset_of_field

  !===================================================================!
  ! Returns the same component in the tuple the rule evaluates, which
  ! stores every state field and then the multipliers, one value each
  !
  !      [ field 1 | field 2 | ... | lambda_1  lambda_2 ... ]
  !                                  ^ num_components()
  !===================================================================!

  pure integer function stored_at(this, position, field, coordinate, order, mixed_order) result(at)

    class(expression), intent(in)           :: this
    integer          , intent(in)           :: position, field, coordinate, order
    integer          , intent(in), optional :: mixed_order

    if (position == ARGUMENT_MULTIPLIER) then
       at = this % num_components() + field - 1
       return
    end if
    at = this % field_offset(field) + this % order_at(field, coordinate, order, mixed_order)

  end function stored_at

  !===================================================================!
  ! Returns the place of a component within its field's components:
  ! the orders along the first coordinate from zero, then one block of
  ! the later coordinates' orders for each mixed order from zero, the
  ! pure components first
  !
  !      u of degree 2 in t, degree 2 in x, mixed degree 1
  !
  !      place   0  1    2     3    4     5     6
  !              u  u_t  u_tt  u_x  u_xx  u_tx  u_txx
  !              along t       mixed 0    mixed 1
  !===================================================================!

  pure integer function order_at(this, field, coordinate, order, mixed_order) result(at)

    class(expression), intent(in)           :: this
    integer          , intent(in)           :: field, coordinate, order
    integer          , intent(in), optional :: mixed_order

    integer :: c, a

    a = 0
    if (present(mixed_order)) a = mixed_order

    if (coordinate < FIRST_COORDINATE .or. coordinate > this % num_coordinates()) then
       error stop 'operation_expression: this component names a coordinate the law does not read'
    end if
    if (coordinate > FIRST_COORDINATE .and. order < 1) then
       error stop 'operation_expression: this component names order zero away from the first &
            &coordinate, where only a positive order is named'
    end if
    if (order < 0 .or. order > merge(this % field_degree(field), this % degrees(coordinate), &
         & coordinate == FIRST_COORDINATE)) then
       error stop 'operation_expression: this component names an order the coordinate does not store'
    end if
    if (a < 0 .or. a > this % mixed_degree_of_field(field) .or. (a > 0 .and. coordinate == FIRST_COORDINATE)) then
       error stop 'operation_expression: this component names a mixed order the field does not store'
    end if

    at = order
    if (coordinate == FIRST_COORDINATE) return
    at = this % field_degree(field) + 1 + a * sum(this % degrees(FIRST_COORDINATE + 1:))
    do c = FIRST_COORDINATE + 1, coordinate - 1
       at = at + this % degrees(c)
    end do
    at = at + order - 1

  end function order_at

  !===================================================================!
  ! Returns how many components one field stores at a point: its own
  ! orders along the first coordinate including zero, and each later
  ! coordinate's orders from one, once for each mixed order from zero
  !===================================================================!

  pure integer function components_per_field(this, field)

    class(expression), intent(in) :: this
    integer          , intent(in) :: field

    components_per_field = this % field_degree(field) + 1 &
         & + (1 + this % mixed_degree_of_field(field)) * sum(this % degrees(FIRST_COORDINATE + 1:))

  end function components_per_field

  !===================================================================!
  ! Returns a field's highest order along the first coordinate among
  ! its mixed components; zero when it has none
  !===================================================================!

  pure integer function mixed_degree(this, field)

    class(expression), intent(in) :: this
    integer          , intent(in) :: field

    if (field < 1 .or. field > this % fields) then
       error stop 'operation_expression: the mixed degree is of a field the rule reads'
    end if
    mixed_degree = this % mixed_degree_of_field(field)

  end function mixed_degree

  !===================================================================!
  ! Returns how many components one point of the state stores: every
  ! state field, the multipliers excluded
  !===================================================================!

  pure integer function num_components(this)

    class(expression), intent(in) :: this

    integer :: f

    if (allocated(this % field_offset)) then
       num_components = this % field_offset(size(this % field_offset))
       return
    end if
    num_components = 0
    do f = 1, this % fields
       num_components = num_components + this % components_per_field(f)
    end do

  end function num_components

  !===================================================================!
  ! Returns the fields of the rule: the state fields and then the
  ! multipliers
  !===================================================================!

  pure integer function num_fields(this)

    class(expression), intent(in) :: this

    num_fields = this % fields + this % multipliers

  end function num_fields

  !===================================================================!
  ! Returns the degree along a coordinate: the first field's along the
  ! first coordinate, every field's along a later one
  !===================================================================!

  pure integer function degree_along(this, coordinate) result(degree)

    class(expression), intent(in) :: this
    integer          , intent(in) :: coordinate

    if (coordinate < FIRST_COORDINATE .or. coordinate > size(this % degrees)) then
       error stop 'operation_expression: this degree names a coordinate the law does not read'
    end if
    degree = this % degrees(coordinate)

  end function degree_along

  !===================================================================!
  ! Returns the count of multipliers
  !===================================================================!

  pure integer function num_multipliers(this)

    class(expression), intent(in) :: this

    num_multipliers = this % multipliers

  end function num_multipliers

  !===================================================================!
  ! Returns the count of designs a point stores for the rule
  !===================================================================!

  pure integer function num_designs(this)

    class(expression), intent(in) :: this

    num_designs = this % designs

  end function num_designs

  !===================================================================!
  ! Returns the field of the one leaf of an expression that is a leaf
  ! alone (u, or u % derivative([t])). An expression of several
  ! vertices, or one whose vertex is not a leaf, stops the program.
  !===================================================================!

  pure integer function leaf_field(this)

    class(expression), intent(in) :: this

    if (this % root() /= 1) error stop 'operation_expression: the leaf named is an expression of several vertices'
    if (this % kind(1) /= VERTEX_LEAF) error stop 'operation_expression: the vertex named is not a leaf'
    leaf_field = this % field(1)

  end function leaf_field

  !===================================================================!
  ! Returns the order of the one leaf, as leaf_field
  !===================================================================!

  pure integer function leaf_order(this)

    class(expression), intent(in) :: this

    if (this % root() /= 1) error stop 'operation_expression: the leaf named is an expression of several vertices'
    if (this % kind(1) /= VERTEX_LEAF) error stop 'operation_expression: the vertex named is not a leaf'
    leaf_order = this % order(1)

  end function leaf_order

  !===================================================================!
  ! Returns the coordinate the one leaf's derivative follows, as
  ! leaf_field
  !===================================================================!

  pure integer function leaf_along(this)

    class(expression), intent(in) :: this

    if (this % root() /= 1) error stop 'operation_expression: the leaf named is an expression of several vertices'
    if (this % kind(1) /= VERTEX_LEAF) error stop 'operation_expression: the vertex named is not a leaf'
    leaf_along = this % along(1)

  end function leaf_along

  !===================================================================!
  ! Returns a field's degree along the first coordinate; a
  ! multiplier's is zero
  !===================================================================!

  pure integer function degree_of_field(this, field) result(degree)

    class(expression), intent(in) :: this
    integer          , intent(in) :: field

    if (field < 1 .or. field > this % fields + this % multipliers) then
       error stop 'operation_expression: this degree names a field the rule does not read'
    end if
    degree = 0
    if (field <= this % fields) degree = this % field_degree(field)

  end function degree_of_field

  !===================================================================!
  ! Returns the coordinates a rule reads from the position of a point:
  ! the highest index named, zero when the rule reads none
  !===================================================================!

  pure integer function num_position_coordinates(this)

    class(expression), intent(in) :: this

    num_position_coordinates = this % position_coordinates

  end function num_position_coordinates

  !===================================================================!
  ! Returns the count of coordinates the state is declared over
  !===================================================================!

  pure integer function num_coordinates(this)

    class(expression), intent(in) :: this

    num_coordinates = size(this % degrees)

  end function num_coordinates

  !===================================================================!
  ! Returns the degree of the equation: the first field's along the
  ! first coordinate
  !===================================================================!

  pure integer function equation_degree(this)

    class(expression)     , intent(in) :: this

    equation_degree = this % degrees(FIRST_COORDINATE)

  end function equation_degree

  !===================================================================!
  ! Applies the rule at every point of the graph, emitting the
  ! coefficient of the full subset: the value with no directions, the
  ! mixed partial with n
  !===================================================================!

  subroutine expression_evaluate(this, input_graph, q, nu, output, positions)

    class(expression)        , intent(in)           :: this
    class(directed_graph)    , intent(in)           :: input_graph
    type(derivative_terms)   , intent(in)           :: q(:), nu(:)
    class(field), allocatable, intent(inout)        :: output
    real(dp)                 , intent(in), optional :: positions(:)

    real(dp), allocatable :: values(:)
    integer               :: k, nd, base, ndesign, coordinates
    character(len=250)    :: message

    nd      = this % num_components()
    ndesign = this % designs

    if (size(q) /= input_graph % num_vertices() * nd) then
       write(message,'(a,i0,a,i0,a,i0)') 'operation_expression: the state does not store one &
            &component per law component per point; size(q) = ', size(q), &
            & ', num_vertices = ', input_graph % num_vertices(), ', components per point = ', nd
       error stop trim(message)
    end if
    if (size(nu) /= input_graph % num_vertices() * ndesign) then
       write(message,'(a,i0,a,i0,a,i0)') 'operation_expression: the design does not store one value &
            &per design per point; size(nu) = ', size(nu), ', num_vertices = ', input_graph % num_vertices(), &
            & ', designs = ', ndesign
       error stop trim(message)
    end if

    coordinates = this % position_coordinates
    if (present(positions)) then
       if (size(positions) /= input_graph % num_vertices() * coordinates) then
          write(message,'(a,i0,a,i0,a,i0)') 'operation_expression: the positions do not store one &
               &value per coordinate read per point; size(positions) = ', size(positions), &
               & ', num_vertices = ', input_graph % num_vertices(), ', coordinates read = ', coordinates
          error stop trim(message)
       end if
    end if

    allocate(values(input_graph % num_vertices()))

    do k = 1, input_graph % num_vertices()
       base = (k - 1) * nd
       if (present(positions)) then
          values(k) = mixed_partial(this % at_instant(q(base + 1:base + nd), &
               & nu((k - 1) * ndesign + 1:k * ndesign), positions((k - 1) * coordinates + 1:k * coordinates)))
       else
          values(k) = mixed_partial(this % at_instant(q(base + 1:base + nd), &
               & nu((k - 1) * ndesign + 1:k * ndesign)))
       end if
    end do

    call emit_real(this % name(), input_graph % vertex_set(), &
         & input_graph % num_vertices(), values, output)

  end subroutine expression_evaluate

  !===================================================================!
  ! Applies the rule: its partial action of order zero
  !===================================================================!

  subroutine expression_apply(this, input_graph, inputs, output)

    class(expression)        , intent(in)           :: this
    class(directed_graph)    , intent(in)           :: input_graph
    type(binding)            , intent(in), optional :: inputs(:)
    class(field), allocatable, intent(inout)        :: output

    call this % value_by_partial_action(input_graph, inputs, output)

  end subroutine expression_apply

  !===================================================================!
  ! Applies the partial action: the state and the design as derivative
  ! terms, read from bindings by argument identity, each seeded by the
  ! variations naming it. A variation naming neither stops the program.
  !===================================================================!

  subroutine expression_partial_action(this, input_graph, inputs, variations, output)

    class(expression)        , intent(in)    :: this
    class(directed_graph)    , intent(in)    :: input_graph
    type(binding)            , intent(in)    :: inputs(:)
    type(variation)          , intent(in)    :: variations(:)
    class(field), allocatable, intent(inout) :: output

    type(derivative_terms), allocatable :: q(:), nu(:)
    real(dp), allocatable               :: positions(:)
    integer                             :: on_state, on_design
    character(len=250)                  :: message

    call this % require_variations(variations)
    call seeded_argument(this, inputs, variations, ARGUMENT_STATE , q , on_state)
    call seeded_argument(this, inputs, variations, ARGUMENT_DESIGN, nu, on_design)
    ! the positions, the third argument of a rule that reads a coordinate
    if (this % position_coordinates > 0) call bound_real_vector(inputs, this % argument(3), positions)
    if (on_state + on_design < size(variations)) then
       write(message,'(a,i0,a,i0,a,i0)') 'operation_expression: a variation names neither the &
            &state nor the design; on_state = ', on_state, ', on_design = ', on_design, &
            & ', size(variations) = ', size(variations)
       error stop trim(message)
    end if
    ! unallocated when the rule reads no coordinate, and then absent
    call this % evaluate(input_graph, q, nu, output, positions)

  end subroutine expression_partial_action

end module operation_expression
