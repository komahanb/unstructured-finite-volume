!=====================================================================!
! A continuous law in named coordinates, stored as data: a graph
! whose vertices are typed operators and whose edges are the reads
! between them.
!
! The law itself has no point graph. Its domain is read from its
! leaves: the fields, the highest derivative order of each field
! along each coordinate, and the multipliers. A discrete domain
! places the law on a directed graph when values are evaluated.
!
!      input graph     the points where the law is evaluated
!      state argument  the state, one field over those points with
!                      num_components() components per point
!      design argument the design, one value per point, so a design
!                      that varies over the points needs no other shape
!      output          one value per point
!
! A leaf reads a component of an argument - the state's component of
! degree d, or the design - or a multiplier, or stores a constant.
! Every other vertex is an arithmetic operation, a power, or an
! elementary function of the vertices it reads. The vertices are
! stored in evaluation order, every read before its reader, and the
! root is the last one.
!
! A rule is built by the intrinsic operators on expressions, so it
! reads as it is written and the compiler checks it:
!
!      q  = unknown()
!      nu = design()
!      r  = derivative(q, 2) + nu * derivative(q, 1) + sin(derivative(q, 0))
!
! derivative(q, d) is the component of degree d: along the instants
! the derivatives are unknowns the scheme relates, so the vertex
! selects and does not differentiate. The state r reads stores the
! components of degree 0, 1 and 2: the highest degree any leaf names.
!
!             SEVERAL FIELDS
!
! unknown(i) is the i-th zeroth-order field; unknown() is the first.
! A rule that reads fields 1..m is over m fields, the highest index
! read. Along the first coordinate each field stores its own jet, to
! the highest degree the rule names for it; along every later
! coordinate the fields share the highest degree named there. A
! point's tuple lists the fields in order, each field's components
! together.
!
!             THE LAGRANGIAN
!
! A rule that reads multiplier(j), j = 1..k, is a Lagrangian, its
! j-th multiplier paired with its j-th field, and the rules it
! generates are its partials. euler_lagrange(l, j) is the
! stationarity of l in its j-th multiplier, the coefficient of one
! direction seeded on that multiplier, and at_zero(l) is l with every
! multiplier at zero. A multiplier is not part of the tuple a rule
! reads, so the state contract excludes it; it is supplied inside the
! evaluation, zero in value, with its own direction when its
! stationarity is read. Nothing is rewritten: the same vertices are
! evaluated with one more direction.
!
!             THE PARTIALS
!
! The rule is evaluated over derivative_terms by one loop over the
! vertices, so the value and every mixed partial in the directions
! requested come from one pass and are exact; nothing is
! differentiated symbolically. A variation may name either argument,
! so the Q-partials, the X-partials and the mixed Q-X partials all
! come from one path and to any degree. A Newton block requires the
! Q-partials as a row of numbers, which is one call per degree with a
! direction that is one at that degree and zero elsewhere.
!
!             WHAT IS REJECTED
!
! A derivative of anything but the unknown, or of negative degree; a
! state whose extent is not the instants times the components; a
! design that is not one value per instant; a missing argument; a
! variation naming neither argument; a function index outside those
! defined; a stationarity in a multiplier the rule does not read.
! Each stops the program.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_expression

  use util_precision       , only : dp
  use operation_action     , only : operation, variation, contract
  use operation_action     , only : binding, seeded_argument, emit_real
  use field_calculus       , only : FIELD_REAL
  use view_directed        , only : directed_graph
  use field_calculus       , only : field
  use util_derivative_terms, only : derivative_terms, mixed_partial, max_subset_width, integer_power
  use util_derivative_terms, only : extend_directions, partial

  implicit none

  private
  public :: expression
  public :: unknown, design, constant, multiplier, derivative, derivative_along
  public :: euler_lagrange, at_zero
  public :: operator(+), operator(-), operator(*), operator(/), operator(**)
  public :: sin, cos, exp, log, sqrt

  ! the vertex kinds
  integer, parameter, public :: VERTEX_LEAF          = 1
  integer, parameter, public :: VERTEX_CONSTANT      = 2
  integer, parameter, public :: VERTEX_SUM           = 3
  integer, parameter, public :: VERTEX_DIFFERENCE    = 4
  integer, parameter, public :: VERTEX_PRODUCT       = 5
  integer, parameter, public :: VERTEX_QUOTIENT      = 6
  integer, parameter, public :: VERTEX_INTEGER_POWER = 7
  integer, parameter, public :: VERTEX_REAL_POWER    = 8
  integer, parameter, public :: VERTEX_FUNCTION      = 9

  ! THE COORDINATE A DERIVATIVE FOLLOWS, named by its place among the
  ! coordinates the state is declared over. One is the first declared;
  ! nothing here specifies which that is, and nothing here limits how
  ! many there are.
  integer, parameter, public :: FIRST_COORDINATE = 1

  ! the arguments a leaf reads, in the operation's order; a multiplier
  ! is read from no argument, and its position names it
  integer, parameter, public :: ARGUMENT_STATE      = 1
  integer, parameter, public :: ARGUMENT_DESIGN     = 2
  integer, parameter, public :: ARGUMENT_MULTIPLIER = 3

  ! the elementary functions
  integer, parameter, public :: SINE        = 1
  integer, parameter, public :: COSINE      = 2
  integer, parameter, public :: EXPONENTIAL = 3
  integer, parameter, public :: LOGARITHM   = 4
  integer, parameter, public :: SQUARE_ROOT = 5

  ! ONE VERTEX AS A VALUE, for a view: the kind, the two vertices
  ! read (0 if none), a leaf's argument, field, order and coordinate,
  ! and a constant or an exponent.
  type, public :: expression_vertex
     integer  :: kind = 0, first = 0, second = 0
     integer  :: position = 0, field = 0, order = 0, along = FIRST_COORDINATE
     real(dp) :: coefficient = 0.0_dp
  end type expression_vertex

  type, extends(operation) :: expression

     ! THE DOMAIN, READ FROM THE LEAVES at every construction and
     ! stored so that evaluation indexes the tuple without a pass over
     ! the vertices: one degree per coordinate (the first includes the
     ! component of order zero, the state itself, so it stores
     ! degrees(1) + 1 components; every later coordinate stores
     ! degrees(c), its orders running from one); the state fields,
     ! each with its degree along the first coordinate; the
     ! multipliers, absent from the tuple and zero in value; and
     ! varied, the multiplier whose stationarity the rule is, or zero.
     integer, allocatable, private :: degrees(:)
     integer, private :: fields      = 0
     integer, private :: multipliers = 0
     integer, private :: varied      = 0
     integer, allocatable, private :: field_degree(:)

     integer , allocatable, private :: kind(:)
     integer , allocatable, private :: first(:), second(:)    ! the vertices read; 0 if none
     integer , allocatable, private :: position(:)            ! a leaf's argument
     integer , allocatable, private :: field(:)               ! a state leaf's field, a multiplier's index
     integer , allocatable, private :: order(:)               ! a leaf's component, or a function index
     integer , allocatable, private :: along(:)               ! the coordinate a leaf's derivative follows
     real(dp), allocatable, private :: coefficient(:)         ! a constant, or an exponent

   contains

     procedure :: apply          => expression_apply
     procedure :: defined_at_zero => expression_defined_at_zero
     procedure :: partial_action => expression_partial_action
     procedure :: at_instant     => expression_at_instant
     procedure :: equation_degree
     procedure :: num_coordinates
     procedure :: num_components
     procedure :: component_at
     procedure :: highest_degree_along
     procedure :: read_components
     procedure :: root
     procedure :: vertices
     procedure :: num_vertices
     procedure :: num_fields
     procedure :: num_multipliers
     procedure :: degree_along
     procedure :: degree_of_field
     procedure :: offset_of_field
     procedure, private :: components_per_field
     procedure, private :: stored_at
     procedure, private :: order_at
     procedure, private :: evaluated_over
     procedure, private :: derive_domain

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

contains

  !===================================================================!
  ! THE LEAVES.
  !===================================================================!

  function vertex(kind, position, order, coefficient, along, field) result(this)

    integer , intent(in) :: kind, position, order
    real(dp), intent(in) :: coefficient
    integer , intent(in), optional :: along, field
    type(expression) :: this

    this % kind        = [kind]
    this % first       = [0]
    this % second      = [0]
    this % position    = [position]
    this % order       = [order]
    this % coefficient = [coefficient]
    this % along       = [FIRST_COORDINATE]
    this % field       = [0]
    if (present(along)) this % along = [along]
    if (present(field)) this % field = [field]

    call this % derive_domain()

  end function vertex

  !===================================================================!
  ! The i-th zeroth-order field, the first when none is named. An
  ! index below one stops the program.
  !===================================================================!

  function unknown(i) result(this)

    integer, intent(in), optional :: i
    type(expression) :: this

    integer :: field
    character(len=250) :: message

    field = 1
    if (present(i)) field = i
    if (field < 1) then
       write(message,'(a,i0)') 'operation_expression: a field must be named by a positive index; field = ', field
       error stop trim(message)
    end if

    this = vertex(VERTEX_LEAF, ARGUMENT_STATE, 0, 0.0_dp, field=field)

  end function unknown

  function design() result(this)

    type(expression) :: this

    this = vertex(VERTEX_LEAF, ARGUMENT_DESIGN, 0, 0.0_dp)

  end function design

  function constant(c) result(this)

    real(dp), intent(in) :: c
    type(expression) :: this

    this = vertex(VERTEX_CONSTANT, 0, 0, c)

  end function constant

  !===================================================================!
  ! The j-th multiplier: read from no argument, zero in value unless
  ! its stationarity is read. An index below one stops the program.
  !===================================================================!

  function multiplier(j) result(this)

    integer, intent(in) :: j
    type(expression) :: this

    character(len=250) :: message

    if (j < 1) then
       write(message,'(a,i0)') 'operation_expression: a multiplier must be named by a positive index; j = ', j
       error stop trim(message)
    end if

    this = vertex(VERTEX_LEAF, ARGUMENT_MULTIPLIER, 0, 0.0_dp, field=j)

  end function multiplier

  !===================================================================!
  ! The component of degree d of the unknown. Anything but a bare
  ! unknown, or a negative degree, stops the program: along the
  ! instants only the unknown has components.
  !===================================================================!

  function derivative(x, d) result(this)

    type(expression), intent(in) :: x
    integer         , intent(in) :: d
    type(expression) :: this

    integer :: n
    character(len=250) :: message

    n = x % root()
    if (n /= 1 .or. x % kind(1) /= VERTEX_LEAF .or. x % position(1) /= ARGUMENT_STATE) then
       write(message,'(a,i0,a,i0,a,i0)') 'operation_expression: derivative requires its &
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

  end function derivative

  !===================================================================!
  ! The same component named along any coordinate: which one is given
  ! by its place among those declared, so a rule reads no differently
  ! whether the coordinate is a time, a length or a sample.
  !===================================================================!

  function derivative_along(x, coordinate, d) result(this)

    type(expression), intent(in) :: x
    integer         , intent(in) :: coordinate, d
    type(expression) :: this

    integer :: n
    character(len=250) :: message

    n = x % root()
    if (n /= 1 .or. x % kind(1) /= VERTEX_LEAF .or. x % position(1) /= ARGUMENT_STATE) then
       write(message,'(a,i0,a,i0,a,i0)') 'operation_expression: derivative_along requires its &
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

  end function derivative_along

  !===================================================================!
  ! The stationarity of a Lagrangian in its j-th multiplier: the
  ! partial of the rule along that multiplier, at zero. A multiplier
  ! the rule does not read stops the program.
  !===================================================================!

  function euler_lagrange(lagrangian, multiplier, label) result(this)

    type(expression), intent(in) :: lagrangian
    integer         , intent(in) :: multiplier
    character(len=*), intent(in), optional :: label
    type(expression) :: this

    character(len=250) :: message

    if (multiplier < 1 .or. multiplier > lagrangian % multipliers) then
       write(message,'(a,i0,a,i0)') 'operation_expression: the stationarity must be in a &
            &multiplier the Lagrangian reads; multiplier = ', multiplier, &
            & ', lagrangian % multipliers = ', lagrangian % multipliers
       error stop trim(message)
    end if

    this = lagrangian
    this % varied = multiplier
    call this % derive_domain(label)

  end function euler_lagrange

  !===================================================================!
  ! The Lagrangian with every multiplier at zero, which for a
  ! Lagrangian linear in them is the functional beside the
  ! constraints.
  !===================================================================!

  function at_zero(lagrangian, label) result(this)

    type(expression), intent(in) :: lagrangian
    character(len=*), intent(in), optional :: label
    type(expression) :: this

    this = lagrangian
    this % varied = 0
    call this % derive_domain(label)

  end function at_zero

  !===================================================================!
  ! COMPOSITION. A binary vertex reads two roots: the right operand's
  ! vertices are appended after the left's and renumbered, and the
  ! new root reads both. A unary vertex reads one.
  !===================================================================!

  function joined(a, b, kind) result(this)

    type(expression), intent(in) :: a, b
    integer         , intent(in) :: kind
    type(expression) :: this

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
    this % coefficient = [a % coefficient, b % coefficient, 0.0_dp]

    call this % derive_domain()

  end function joined

  function applied(a, kind, order, coefficient) result(this)

    type(expression), intent(in) :: a
    integer         , intent(in) :: kind, order
    real(dp)        , intent(in) :: coefficient
    type(expression) :: this

    integer :: na

    na = a % root()

    this % kind        = [a % kind,        kind]
    this % first       = [a % first,       na]
    this % second      = [a % second,      0]
    this % position    = [a % position,    0]
    this % field       = [a % field,       0]
    this % order       = [a % order,       order]
    this % along       = [a % along,       FIRST_COORDINATE]
    this % coefficient = [a % coefficient, coefficient]

    call this % derive_domain()

  end function applied

  pure function shifted(reads, by) result(moved)

    integer, intent(in) :: reads(:), by
    integer :: moved(size(reads))

    where (reads > 0)
       moved = reads + by
    elsewhere
       moved = 0
    end where

  end function shifted

  !===================================================================!
  ! THE DOMAIN, READ FROM THE LEAVES: the fields are those any unknown
  ! names, one at least; a field's degree along the first coordinate
  ! is the highest order named for it, zero when only its value is
  ! named; along a later coordinate the degree is the highest order
  ! named there by any field; the multipliers are those any
  ! multiplier leaf names. The argument contracts follow, with the
  ! highest exact degree, one less for a stationarity.
  !===================================================================!

  subroutine derive_domain(this, label)

    class(expression), intent(inout) :: this
    character(len=*) , intent(in), optional :: label

    integer :: i, f, c, coordinates

    this % fields      = 1
    this % multipliers = 0
    coordinates        = FIRST_COORDINATE
    do i = 1, this % num_vertices()
       if (this % kind(i) /= VERTEX_LEAF) cycle
       select case (this % position(i))
       case (ARGUMENT_STATE)
          this % fields = max(this % fields, this % field(i))
          coordinates   = max(coordinates, this % along(i))
       case (ARGUMENT_MULTIPLIER)
          this % multipliers = max(this % multipliers, this % field(i))
       end select
    end do

    this % field_degree = [(max(0, this % highest_degree_along(FIRST_COORDINATE, f)), f = 1, this % fields)]
    this % degrees      = [this % field_degree(1), &
         & (max(0, this % highest_degree_along(c)), c = FIRST_COORDINATE + 1, coordinates)]

    call this % declare_arguments(2, [ &
         & contract(FIELD_REAL, this % num_components()), &
         & contract(FIELD_REAL, 1) ], label=label, &
         & max_degree=max_subset_width() - merge(1, 0, this % varied > 0))

  end subroutine derive_domain

  !===================================================================!
  ! THE OPERATORS.
  !===================================================================!

  function plus(a, b) result(this)
    type(expression), intent(in) :: a, b
    type(expression) :: this
    this = joined(a, b, VERTEX_SUM)
  end function plus

  function real_plus(x, b) result(this)
    real(dp)        , intent(in) :: x
    type(expression), intent(in) :: b
    type(expression) :: this
    this = joined(constant(x), b, VERTEX_SUM)
  end function real_plus

  function plus_real(a, x) result(this)
    type(expression), intent(in) :: a
    real(dp)        , intent(in) :: x
    type(expression) :: this
    this = joined(a, constant(x), VERTEX_SUM)
  end function plus_real

  function minus(a, b) result(this)
    type(expression), intent(in) :: a, b
    type(expression) :: this
    this = joined(a, b, VERTEX_DIFFERENCE)
  end function minus

  function real_minus(x, b) result(this)
    real(dp)        , intent(in) :: x
    type(expression), intent(in) :: b
    type(expression) :: this
    this = joined(constant(x), b, VERTEX_DIFFERENCE)
  end function real_minus

  function minus_real(a, x) result(this)
    type(expression), intent(in) :: a
    real(dp)        , intent(in) :: x
    type(expression) :: this
    this = joined(a, constant(x), VERTEX_DIFFERENCE)
  end function minus_real

  function negated(a) result(this)
    type(expression), intent(in) :: a
    type(expression) :: this
    this = joined(constant(-1.0_dp), a, VERTEX_PRODUCT)
  end function negated

  function times(a, b) result(this)
    type(expression), intent(in) :: a, b
    type(expression) :: this
    this = joined(a, b, VERTEX_PRODUCT)
  end function times

  function real_times(x, b) result(this)
    real(dp)        , intent(in) :: x
    type(expression), intent(in) :: b
    type(expression) :: this
    this = joined(constant(x), b, VERTEX_PRODUCT)
  end function real_times

  function times_real(a, x) result(this)
    type(expression), intent(in) :: a
    real(dp)        , intent(in) :: x
    type(expression) :: this
    this = joined(a, constant(x), VERTEX_PRODUCT)
  end function times_real

  function over(a, b) result(this)
    type(expression), intent(in) :: a, b
    type(expression) :: this
    this = joined(a, b, VERTEX_QUOTIENT)
  end function over

  function real_over(x, b) result(this)
    real(dp)        , intent(in) :: x
    type(expression), intent(in) :: b
    type(expression) :: this
    this = joined(constant(x), b, VERTEX_QUOTIENT)
  end function real_over

  function over_real(a, x) result(this)
    type(expression), intent(in) :: a
    real(dp)        , intent(in) :: x
    type(expression) :: this
    this = joined(a, constant(x), VERTEX_QUOTIENT)
  end function over_real

  function to_integer(a, n) result(this)
    type(expression), intent(in) :: a
    integer         , intent(in) :: n
    type(expression) :: this
    this = applied(a, VERTEX_INTEGER_POWER, 0, real(n, dp))
  end function to_integer

  function to_real(a, p) result(this)
    type(expression), intent(in) :: a
    real(dp)        , intent(in) :: p
    type(expression) :: this
    this = applied(a, VERTEX_REAL_POWER, 0, p)
  end function to_real

  function sine_of(a) result(this)
    type(expression), intent(in) :: a
    type(expression) :: this
    this = applied(a, VERTEX_FUNCTION, SINE, 0.0_dp)
  end function sine_of

  function cosine_of(a) result(this)
    type(expression), intent(in) :: a
    type(expression) :: this
    this = applied(a, VERTEX_FUNCTION, COSINE, 0.0_dp)
  end function cosine_of

  function exponential_of(a) result(this)
    type(expression), intent(in) :: a
    type(expression) :: this
    this = applied(a, VERTEX_FUNCTION, EXPONENTIAL, 0.0_dp)
  end function exponential_of

  function logarithm_of(a) result(this)
    type(expression), intent(in) :: a
    type(expression) :: this
    this = applied(a, VERTEX_FUNCTION, LOGARITHM, 0.0_dp)
  end function logarithm_of

  function root_of(a) result(this)
    type(expression), intent(in) :: a
    type(expression) :: this
    this = applied(a, VERTEX_FUNCTION, SQUARE_ROOT, 0.0_dp)
  end function root_of

  !===================================================================!
  ! EVALUATION: the loop over the vertices, over derivative_terms. A
  ! state component past those given, or a function index outside
  ! those defined, stops the program.
  !===================================================================!

  pure function expression_at_instant(this, q, nu) result(r)

    class(expression)     , intent(in) :: this
    type(derivative_terms), intent(in) :: q(0:)
    type(derivative_terms), intent(in) :: nu
    type(derivative_terms) :: r

    type(derivative_terms), allocatable :: stored(:)
    type(derivative_terms) :: design
    integer :: n, j, given

    if (size(q) /= this % num_components()) then
       error stop 'operation_expression: the point does not store one component per law component'
    end if
    if (this % multipliers == 0) then
       r = this % evaluated_over(q, nu)
       return
    end if

    ! the tuple evaluated stores the state and then the multipliers;
    ! a multiplier, read from no input, is zero, and the varied one
    ! has one more direction, the highest, whose partial is returned
    ! over the caller's
    n = nu % num_directions()
    if (this % varied > 0) n = n + 1
    design = extend_directions(nu, n)
    given  = this % num_components()
    allocate(stored(0:given + this % multipliers - 1))
    do j = 0, given - 1
       stored(j) = extend_directions(q(j), n)
    end do
    do j = 1, this % multipliers
       stored(given + j - 1) = derivative_terms(0.0_dp, design)
    end do
    if (this % varied > 0) call stored(given + this % varied - 1) % set_direction(n, 1.0_dp)

    r = this % evaluated_over(stored, design)
    if (this % varied > 0) r = partial(r, n)

  end function expression_at_instant

  !===================================================================!
  ! WHETHER THE ZERO STATE LIES IN THE EXPRESSION'S DOMAIN. A vertex
  ! reads the state when it is a state component or when either of
  ! its operands does. The expression has a value at the zero state
  ! unless such a vertex is a divisor, the base of a negative integer
  ! power, the base of a real power below the first, or the argument
  ! of a logarithm or a square root - each of which is singular, or
  ! has a singular derivative, at zero. The design is an input of its
  ! own and is not zeroed with the state.
  !===================================================================!

  pure logical function expression_defined_at_zero(this) result(defined)

    class(expression), intent(in) :: this

    logical, allocatable :: reads_state(:)
    integer :: i

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
  ! The loop over the vertices, on the stored tuple.
  !===================================================================!

  pure function evaluated_over(this, q, nu) result(r)

    ! the arithmetic over derivative_terms is read here alone, so the
    ! module's public operators are its own, on expressions
    use util_derivative_terms, only : operator(+), operator(-), operator(*), operator(/), &
         & operator(**), sin, cos, exp, log, sqrt

    class(expression)     , intent(in) :: this
    type(derivative_terms), intent(in) :: q(0:)
    type(derivative_terms), intent(in) :: nu
    type(derivative_terms) :: r

    type(derivative_terms), allocatable :: v(:)
    integer :: i, at

    allocate(v(this % num_vertices()))

    do i = 1, this % num_vertices()
       select case (this % kind(i))
       case (VERTEX_LEAF)
          if (this % position(i) == ARGUMENT_DESIGN) then
             v(i) = nu
          else
             ! the tuple lists every field's components together, the
             ! components along time first, then along space, and the
             ! multipliers after every field
             at = this % stored_at(this % position(i), this % field(i), this % along(i), this % order(i))
             if (at > ubound(q, 1)) then
                error stop 'operation_expression: the state does not store the component read, &
                     &beyond the end of the stored tuple'
             end if
             v(i) = q(at)
          end if
       case (VERTEX_CONSTANT)
          v(i) = derivative_terms(this % coefficient(i), nu)
       case (VERTEX_SUM)
          v(i) = v(this % first(i)) + v(this % second(i))
       case (VERTEX_DIFFERENCE)
          v(i) = v(this % first(i)) - v(this % second(i))
       case (VERTEX_PRODUCT)
          v(i) = v(this % first(i)) * v(this % second(i))
       case (VERTEX_QUOTIENT)
          v(i) = v(this % first(i)) / v(this % second(i))
       case (VERTEX_INTEGER_POWER)
          v(i) = integer_power(v(this % first(i)), nint(this % coefficient(i)))
       case (VERTEX_REAL_POWER)
          v(i) = v(this % first(i)) ** this % coefficient(i)
       case (VERTEX_FUNCTION)
          select case (this % order(i))
          case (SINE);        v(i) = sin(v(this % first(i)))
          case (COSINE);      v(i) = cos(v(this % first(i)))
          case (EXPONENTIAL); v(i) = exp(v(this % first(i)))
          case (LOGARITHM);   v(i) = log(v(this % first(i)))
          case (SQUARE_ROOT); v(i) = sqrt(v(this % first(i)))
          case default
             error stop 'operation_expression: this vertex names a function that is not one of &
                  &those defined'
          end select
       case default
          error stop 'operation_expression: this vertex names a kind that is not one of those &
               &defined'
       end select
    end do

    r = v(size(v))

  end function evaluated_over

  !===================================================================!
  ! The highest order this rule names along one coordinate, and minus
  ! one when it names none along it. One query serves every
  ! coordinate, so adding a coordinate requires no new procedure.
  !===================================================================!

  pure integer function highest_degree_along(this, coordinate, field)

    class(expression), intent(in) :: this
    integer          , intent(in) :: coordinate
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
  ! The components of the tuple the rule reads, each once, in
  ! increasing order: the state leaves' components, a multiplier's
  ! excluded since the tuple does not store it.
  !===================================================================!

  subroutine read_components(this, components)

    class(expression)   , intent(in)  :: this
    integer, allocatable, intent(out) :: components(:)

    logical, allocatable :: read(:)
    integer :: i, at

    allocate(read(0:this % num_components() - 1), source=.false.)
    do i = 1, this % num_vertices()
       if (this % kind(i) /= VERTEX_LEAF .or. this % position(i) /= ARGUMENT_STATE) cycle
       at = this % component_at(this % along(i), this % order(i), this % field(i))
       read(at) = .true.
    end do
    components = pack([(at, at = 0, size(read) - 1)], read)

  end subroutine read_components

  !===================================================================!
  ! Invalid input: a variable of type expression never assigned from
  ! unknown, design, constant or an operator has no vertex, so no
  ! root and no value. Every read of the vertex count passes through
  ! here, so this is the one refusal; without it the compiled code
  ! reads an array with no memory and faults in memmove with no
  ! message.
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
  ! The vertices are stored in evaluation order, so the root, the
  ! vertex whose value is the value of the expression, is the last.
  !===================================================================!

  pure integer function root(this)

    class(expression), intent(in) :: this

    root = this % num_vertices()

  end function root

  !===================================================================!
  ! The vertex table as values, in evaluation order, for a view.
  !===================================================================!

  function vertices(this) result(v)

    class(expression), intent(in) :: this
    type(expression_vertex), allocatable :: v(:)

    integer :: i

    allocate(v(this % num_vertices()))
    do i = 1, size(v)
       v(i) = expression_vertex(this % kind(i), this % first(i), this % second(i), this % position(i), &
            & this % field(i), this % order(i), this % along(i), this % coefficient(i))
    end do

  end function vertices

  !===================================================================!
  ! THE COMPONENT A COORDINATE AND AN ORDER NAME. The first
  ! coordinate's orders run from zero and every later one's from one,
  ! since order zero is the state itself and is named once.
  !===================================================================!

  pure integer function component_at(this, coordinate, order, field) result(at)

    class(expression), intent(in) :: this
    integer          , intent(in) :: coordinate, order
    integer          , intent(in), optional :: field

    integer :: f

    f = 1
    if (present(field)) f = field
    at = this % offset_of_field(f) + this % order_at(f, coordinate, order)

  end function component_at

  !===================================================================!
  ! The first component of a field in the tuple the rule reads, which
  ! stores the state fields alone. A multiplier stops the program.
  !===================================================================!

  pure integer function offset_of_field(this, field) result(at)

    class(expression), intent(in) :: this
    integer          , intent(in) :: field

    integer :: f

    if (field < 1 .or. field > this % fields + this % multipliers) then
       error stop 'operation_expression: this component names a field the rule does not read'
    end if
    if (field > this % fields) then
       error stop 'operation_expression: this component names a multiplier, which is not stored &
            &in the tuple'
    end if
    at = 0
    do f = 1, field - 1
       at = at + this % components_per_field(f)
    end do

  end function offset_of_field

  !===================================================================!
  ! The same component in the tuple the rule evaluates, which stores
  ! every state field and then the multipliers, one value each.
  !===================================================================!

  pure integer function stored_at(this, position, field, coordinate, order) result(at)

    class(expression), intent(in) :: this
    integer          , intent(in) :: position, field, coordinate, order

    integer :: f

    if (position == ARGUMENT_MULTIPLIER) then
       at = this % num_components() + field - 1
       return
    end if
    at = 0
    do f = 1, field - 1
       at = at + this % components_per_field(f)
    end do
    at = at + this % order_at(field, coordinate, order)

  end function stored_at

  pure integer function order_at(this, field, coordinate, order) result(at)

    class(expression), intent(in) :: this
    integer          , intent(in) :: field, coordinate, order

    integer :: c

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

    at = order
    if (coordinate == FIRST_COORDINATE) return
    at = this % field_degree(field) + 1
    do c = FIRST_COORDINATE + 1, coordinate - 1
       at = at + this % degrees(c)
    end do
    at = at + order - 1

  end function order_at

  !===================================================================!
  ! How many components one field stores at a point: its own orders
  ! along the first coordinate including zero, and each later
  ! coordinate's orders from one.
  !===================================================================!

  pure integer function components_per_field(this, field)

    class(expression), intent(in) :: this
    integer          , intent(in) :: field

    components_per_field = this % field_degree(field) + 1
    if (size(this % degrees) > 1) components_per_field = components_per_field + sum(this % degrees(2:))

  end function components_per_field

  !===================================================================!
  ! How many components one point of the state stores: every state
  ! field, the multipliers excluded.
  !===================================================================!

  pure integer function num_components(this)

    class(expression), intent(in) :: this

    integer :: f

    num_components = 0
    do f = 1, this % fields
       num_components = num_components + this % components_per_field(f)
    end do

  end function num_components

  !===================================================================!
  ! The fields of the rule: the state fields and then the multipliers.
  !===================================================================!

  pure integer function num_fields(this)

    class(expression), intent(in) :: this

    num_fields = this % fields + this % multipliers

  end function num_fields

  !===================================================================!
  ! The degree along a coordinate: the first field's along the first
  ! coordinate, every field's along a later one.
  !===================================================================!

  pure integer function degree_along(this, coordinate) result(degree)

    class(expression), intent(in) :: this
    integer          , intent(in) :: coordinate

    if (coordinate < FIRST_COORDINATE .or. coordinate > size(this % degrees)) then
       error stop 'operation_expression: this degree names a coordinate the law does not read'
    end if
    degree = this % degrees(coordinate)

  end function degree_along

  pure integer function num_multipliers(this)

    class(expression), intent(in) :: this

    num_multipliers = this % multipliers

  end function num_multipliers

  !===================================================================!
  ! A field's degree along the first coordinate; a multiplier's is
  ! zero.
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

  pure integer function num_coordinates(this)

    class(expression), intent(in) :: this

    num_coordinates = size(this % degrees)

  end function num_coordinates

  pure integer function equation_degree(this)

    class(expression)     , intent(in) :: this

    equation_degree = this % degrees(FIRST_COORDINATE)

  end function equation_degree

  !===================================================================!
  ! The rule at every instant, returning the coefficient of the full
  ! subset: the value with no directions, the mixed partial with n.
  !===================================================================!

  subroutine evaluated(this, input_graph, q, nu, output)

    class(expression)     , intent(in) :: this
    class(directed_graph) , intent(in) :: input_graph
    type(derivative_terms), intent(in) :: q(:), nu(:)
    class(field), allocatable, intent(inout) :: output

    real(dp), allocatable :: values(:)
    integer :: k, nd, base
    character(len=250) :: message

    nd = this % num_components()

    if (size(q) /= input_graph % num_vertices() * nd) then
       write(message,'(a,i0,a,i0,a,i0)') 'operation_expression: the state does not store one &
            &component per law component per point; size(q) = ', size(q), &
            & ', num_vertices = ', input_graph % num_vertices(), ', components per point = ', nd
       error stop trim(message)
    end if
    if (size(nu) /= input_graph % num_vertices()) then
       write(message,'(a,i0,a,i0)') 'operation_expression: the design does not store one value &
            &per point; size(nu) = ', size(nu), ', num_vertices = ', input_graph % num_vertices()
       error stop trim(message)
    end if

    allocate(values(input_graph % num_vertices()))

    do k = 1, input_graph % num_vertices()
       base = (k - 1) * nd
       values(k) = mixed_partial(this % at_instant(q(base + 1:base + nd), nu(k)))
    end do

    call emit_real(this % name(), input_graph % vertex_set(), &
         & input_graph % num_vertices(), values, output)

  end subroutine evaluated

  subroutine expression_apply(this, input_graph, inputs, output)

    class(expression)     , intent(in)       :: this
    class(directed_graph) , intent(in)       :: input_graph
    type(binding), intent(in), optional       :: inputs(:)
    class(field), allocatable, intent(inout) :: output

    call this % value_by_partial_action(input_graph, inputs, output)

  end subroutine expression_apply

  !===================================================================!
  ! The state and the design as derivative terms, read from bindings
  ! by argument identity, each seeded by the variations naming it. A
  ! variation naming neither stops the program.
  !===================================================================!

  subroutine expression_partial_action(this, input_graph, inputs, variations, output)

    class(expression)     , intent(in)       :: this
    class(directed_graph) , intent(in)       :: input_graph
    type(binding)          , intent(in)       :: inputs(:)
    type(variation)       , intent(in)       :: variations(:)
    class(field), allocatable, intent(inout) :: output

    type(derivative_terms), allocatable :: q(:), nu(:)
    integer :: on_state, on_design
    character(len=250) :: message

    call this % require_variations(variations)
    call seeded_argument(this, inputs, variations, ARGUMENT_STATE , q , on_state)
    call seeded_argument(this, inputs, variations, ARGUMENT_DESIGN, nu, on_design)
    if (on_state + on_design < size(variations)) then
       write(message,'(a,i0,a,i0,a,i0)') 'operation_expression: a variation names neither the &
            &state nor the design; on_state = ', on_state, ', on_design = ', on_design, &
            & ', size(variations) = ', size(variations)
       error stop trim(message)
    end if
    call evaluated(this, input_graph, q, nu, output)

  end subroutine expression_partial_action

end module operation_expression
