!=====================================================================!
! A continuous law in named coordinates, stored as data: a graph
! whose vertices are typed operators and whose edges are the reads
! between them.
!
! The law itself has no point graph. A caller first states the
! highest derivative degree for each coordinate, then a discrete
! domain places that law on a directed graph when values are
! evaluated.
!
!      input graph     the points where the law is evaluated
!      state argument  the state, one field over those points with
!                      num_components() components per point
!      design argument the design, one value per point, so a design
!                      that varies over the points needs no other shape
!      output          one value per point
!
! A leaf reads a component of an argument - the state's component of
! degree d, or the design - or stores a constant. Every other vertex
! is an arithmetic operation, a power, or an elementary function of
! the vertices it reads. The vertices are stored in evaluation order,
! every read before its reader, and the root is the last one.
!
! A rule is built by the intrinsic operators on expressions, so it
! reads as it is written and the compiler checks it:
!
!      q  = unknown()
!      nu = design()
!      r  = stated(derivative(q, 2) + nu * derivative(q, 1) &
!                + sin(derivative(q, 0)), 2, 'the residual')
!
! derivative(q, d) is the component of degree d: along the instants
! the derivatives are unknowns the scheme relates, so the vertex
! selects and does not differentiate.
!
!             SEVERAL FIELDS
!
! unknown(i) is the i-th zeroth-order field; unknown() is the first.
! A rule that reads fields 1..m is stated over m fields, the highest
! index read. Along the first coordinate each field stores its own
! jet, to the degree given for it when the rule is stated (the
! equation's degree for the first field, zero for a multiplier,
! unless listed); along every later coordinate the fields share the
! declared degrees. A point's tuple lists the fields in order, each
! field's components together.
!
!             THE LAGRANGIAN
!
! A rule stated with k multipliers is a Lagrangian: its last k fields
! are the multipliers, paired in order with its first k fields, and
! the rules it generates are its partials. euler_lagrange(l, j) is
! the stationarity of l in its j-th multiplier, the coefficient of
! one direction seeded on that field's component of order zero, and
! at_zero(l) is l with every multiplier at zero. A multiplier is not
! part of the tuple a rule reads, so the state contract excludes
! it; it is supplied inside the evaluation, zero in value, with its
! own direction when its stationarity is read. Nothing is rewritten:
! the same vertices are evaluated with one more direction.
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
! rule stated at a degree below one, or below the highest component
! it reads; a state whose extent is not the instants times N+1; a
! design that is not one value per instant; a missing argument; a
! variation naming neither argument; a function index outside those
! defined. Each stops the program.
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
  public :: unknown, design, constant, derivative, derivative_along, stated, stated_over
  public :: euler_lagrange, at_zero
  public :: operator(+), operator(-), operator(*), operator(/), operator(**)
  public :: sin, cos, exp, log, sqrt

  ! the vertex kinds
  integer, parameter :: VERTEX_LEAF          = 1
  integer, parameter :: VERTEX_CONSTANT      = 2
  integer, parameter :: VERTEX_SUM           = 3
  integer, parameter :: VERTEX_DIFFERENCE    = 4
  integer, parameter :: VERTEX_PRODUCT       = 5
  integer, parameter :: VERTEX_QUOTIENT      = 6
  integer, parameter :: VERTEX_INTEGER_POWER = 7
  integer, parameter :: VERTEX_REAL_POWER    = 8
  integer, parameter :: VERTEX_FUNCTION      = 9

  ! the arguments a leaf reads, in the operation's order
  ! THE COORDINATE A DERIVATIVE FOLLOWS, named by its place among the
  ! coordinates the state is declared over. One is the first declared;
  ! nothing here specifies which that is, and nothing here limits how
  ! many there are.
  integer, parameter, public :: FIRST_COORDINATE = 1

  integer, parameter :: ARGUMENT_STATE  = 1
  integer, parameter :: ARGUMENT_DESIGN = 2

  ! the elementary functions
  integer, parameter :: SINE        = 1
  integer, parameter :: COSINE      = 2
  integer, parameter :: EXPONENTIAL = 3
  integer, parameter :: LOGARITHM   = 4
  integer, parameter :: ROOT        = 5

  type, extends(operation) :: expression

     ! ONE DEGREE PER COORDINATE, in the order the coordinates are
     ! declared. The first includes the component of order zero, which
     ! is the state itself, so it stores degrees(1) + 1 components and
     ! every later coordinate stores degrees(c), its orders running
     ! from one.
     integer, allocatable, private :: degrees(:)

     ! the fields the rule is stated over, each with its degree along
     ! the first coordinate; the last `multipliers` of them are absent
     ! from the tuple and zero in value, and `varied` names the one
     ! whose stationarity the rule is, or is zero
     integer, private :: fields      = 0
     integer, private :: multipliers = 0
     integer, private :: varied      = 0
     integer, allocatable, private :: field_degree(:)

     integer , allocatable, private :: kind(:)
     integer , allocatable, private :: first(:), second(:)    ! the vertices read; 0 if none
     integer , allocatable, private :: position(:)            ! a leaf's argument
     integer , allocatable, private :: field(:)               ! a state leaf's field
     integer , allocatable, private :: order(:)               ! a leaf's component, or a function index
     integer , allocatable, private :: along(:)               ! the coordinate a leaf's derivative follows
     real(dp), allocatable, private :: coefficient(:)         ! a constant, or an exponent

   contains

     procedure :: apply          => expression_apply
     procedure :: partial_action => expression_partial_action
     procedure :: at_instant     => expression_at_instant
     procedure :: declared        => expression_declared
     procedure :: equation_degree
     procedure :: num_coordinates
     procedure :: num_components
     procedure :: component_at
     procedure :: declare_degree
     procedure :: highest_degree_along
     procedure :: read_components
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
     procedure, private :: restated

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

  end function vertex

  !===================================================================!
  ! The i-th zeroth-order field, the first when none is named. An
  ! index below one stops the program.
  !===================================================================!

  function unknown(i) result(this)

    integer, intent(in), optional :: i
    type(expression) :: this

    integer :: field

    field = 1
    if (present(i)) field = i
    if (field < 1) then
       error stop 'operation_expression: a field is named by a positive index'
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
  ! The component of degree d of the unknown. Anything but a bare
  ! unknown, or a negative degree, stops the program: along the
  ! instants only the unknown has components.
  !===================================================================!

  function derivative(x, d) result(this)

    type(expression), intent(in) :: x
    integer         , intent(in) :: d
    type(expression) :: this

    if (size(x % kind) /= 1 .or. x % kind(1) /= VERTEX_LEAF .or. x % position(1) /= ARGUMENT_STATE) then
       error stop 'operation_expression: a derivative is taken of the unknown'
    end if
    if (d < 0) then
       error stop 'operation_expression: the degree of a derivative is not negative'
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

    if (size(x % kind) /= 1 .or. x % kind(1) /= VERTEX_LEAF .or. x % position(1) /= ARGUMENT_STATE) then
       error stop 'operation_expression: a derivative is taken of the unknown'
    end if
    if (d < 0) then
       error stop 'operation_expression: the degree of a derivative is not negative'
    end if
    if (coordinate < FIRST_COORDINATE) then
       error stop 'operation_expression: a coordinate is one of those declared'
    end if

    this = vertex(VERTEX_LEAF, ARGUMENT_STATE, d, 0.0_dp, coordinate, x % field(1))

  end function derivative_along

  !===================================================================!
  ! A rule bound to the degree of the state it reads. A degree below
  ! the highest component read stops the program: the state stores no
  ! such component.
  !===================================================================!

  function stated(rule, degree, label, field_degrees, multipliers) result(this)

    type(expression), intent(in) :: rule
    integer         , intent(in) :: degree
    character(len=*), intent(in) :: label
    integer         , intent(in), optional :: field_degrees(:), multipliers
    type(expression) :: this

    this = stated_over(rule, [degree], label, field_degrees, multipliers)

  end function stated

  !===================================================================!
  ! A rule bound to a degree along each coordinate the state is
  ! declared over. The one-coordinate form is this with a list of
  ! one, so a rule over time alone reads no differently than before.
  !===================================================================!

  function stated_over(rule, degrees, label, field_degrees, multipliers) result(this)

    type(expression), intent(in) :: rule
    integer         , intent(in) :: degrees(:)
    character(len=*), intent(in) :: label
    integer         , intent(in), optional :: field_degrees(:), multipliers
    type(expression) :: this

    integer :: c, f

    ! a rule stated again keeps its fields' degrees, its multipliers
    ! and the stationarity it is, unless these are given again
    this = rule
    this % fields = max(1, maxval(this % field))
    if (present(multipliers)) then
       this % multipliers = multipliers
       this % varied      = 0
    end if
    if (this % multipliers < 0 .or. 2 * this % multipliers > this % fields) then
       error stop 'operation_expression: the multipliers are the last fields, each paired with a state field'
    end if

    if (size(degrees) < 1) then
       error stop 'operation_expression: a state is declared over one coordinate at least'
    end if
    if (.not. allocated(this % field_degree)) then
       allocate(this % field_degree(this % fields))
       this % field_degree = degrees(FIRST_COORDINATE)
       if (this % multipliers > 0) this % field_degree(this % fields - this % multipliers + 1:) = 0
    else if (size(this % field_degree) /= this % fields) then
       error stop 'operation_expression: a rule stated again reads the fields it was stated over'
    end if
    if (present(field_degrees)) then
       if (size(field_degrees) /= this % fields) then
          error stop 'operation_expression: one degree per field the rule is stated over'
       end if
       this % field_degree = field_degrees
    end if
    if (this % field_degree(1) /= degrees(FIRST_COORDINATE)) then
       error stop 'operation_expression: the first field stores the equation''s degree'
    end if

    do f = 1, this % fields
       if (this % highest_degree_along(FIRST_COORDINATE, f) > this % field_degree(f)) then
          error stop 'operation_expression: the rule reads a component the state stores'
       end if
    end do
    do c = FIRST_COORDINATE + 1, size(degrees)
       if (this % highest_degree_along(c) > degrees(c)) then
          error stop 'operation_expression: the rule reads a component the state stores'
       end if
    end do
    if (this % highest_degree_along(size(degrees) + 1) >= 0) then
       error stop 'operation_expression: the rule reads a coordinate the state is not declared over'
    end if

    call this % declare_degree(degrees, label)

  end function stated_over

  !===================================================================!
  ! The stationarity of a Lagrangian in its j-th multiplier: the
  ! partial of the rule along that field's component of order zero,
  ! at zero. A rule not stated, or a multiplier it does not declare,
  ! stops the program.
  !===================================================================!

  function euler_lagrange(lagrangian, multiplier, label) result(this)

    type(expression), intent(in) :: lagrangian
    integer         , intent(in) :: multiplier
    character(len=*), intent(in), optional :: label
    type(expression) :: this

    if (.not. lagrangian % declared()) then
       error stop 'operation_expression: a Lagrangian is stated before its stationarity is read'
    end if
    if (multiplier < 1 .or. multiplier > lagrangian % multipliers) then
       error stop 'operation_expression: the stationarity is in a multiplier the Lagrangian declares'
    end if

    this = lagrangian
    this % varied = lagrangian % fields - lagrangian % multipliers + multiplier
    call this % restated(label)

  end function euler_lagrange

  !===================================================================!
  ! The Lagrangian with every multiplier at zero, which for a
  ! Lagrangian linear in them is the functional beside the
  ! constraints. A rule not stated stops the program.
  !===================================================================!

  function at_zero(lagrangian, label) result(this)

    type(expression), intent(in) :: lagrangian
    character(len=*), intent(in), optional :: label
    type(expression) :: this

    if (.not. lagrangian % declared()) then
       error stop 'operation_expression: a Lagrangian is stated before its value at zero is read'
    end if

    this = lagrangian
    this % varied = 0
    call this % restated(label)

  end function at_zero

  subroutine restated(this, label)

    class(expression), intent(inout) :: this
    character(len=*) , intent(in), optional :: label

    integer, allocatable :: degrees(:)

    degrees = this % degrees
    if (present(label)) then
       call this % declare_degree(degrees, label)
    else
       call this % declare_degree(degrees, this % name())
    end if

  end subroutine restated

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

    na = size(a % kind)
    nb = size(b % kind)

    this % kind        = [a % kind,        b % kind,        kind]
    this % first       = [a % first,       shifted(b % first,  na), na]
    this % second      = [a % second,      shifted(b % second, na), na + nb]
    this % position    = [a % position,    b % position,    0]
    this % field       = [a % field,       b % field,       0]
    this % order       = [a % order,       b % order,       0]
    this % along       = [a % along,       b % along,       FIRST_COORDINATE]
    this % coefficient = [a % coefficient, b % coefficient, 0.0_dp]

  end function joined

  function applied(a, kind, order, coefficient) result(this)

    type(expression), intent(in) :: a
    integer         , intent(in) :: kind, order
    real(dp)        , intent(in) :: coefficient
    type(expression) :: this

    this % kind        = [a % kind,        kind]
    this % first       = [a % first,       size(a % kind)]
    this % second      = [a % second,      0]
    this % position    = [a % position,    0]
    this % field       = [a % field,       0]
    this % order       = [a % order,       order]
    this % along       = [a % along,       FIRST_COORDINATE]
    this % coefficient = [a % coefficient, coefficient]

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
    this = applied(a, VERTEX_FUNCTION, ROOT, 0.0_dp)
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
    integer :: per, f, n, j, at, given, total

    if (size(q) /= this % num_components()) then
       error stop 'operation_expression: a point stores one component per law component'
    end if
    if (this % multipliers == 0) then
       r = this % evaluated_over(q, nu)
       return
    end if

    ! the tuple evaluated stores every field; a multiplier, read from
    ! no input, is zero, and the varied one has one more direction,
    ! the highest, whose partial is returned over the caller's
    n = nu % num_directions()
    if (this % varied > 0) n = n + 1
    design = extend_directions(nu, n)
    total = 0
    do f = 1, this % fields
       total = total + this % components_per_field(f)
    end do
    allocate(stored(0:total - 1))
    at    = 0
    given = 0
    do f = 1, this % fields
       per = this % components_per_field(f)
       if (f > this % fields - this % multipliers) then
          stored(at:at + per - 1) = derivative_terms(0.0_dp, design)
          if (f == this % varied) call stored(at) % set_direction(n, 1.0_dp)
       else
          do j = 0, per - 1
             stored(at + j) = extend_directions(q(given + j), n)
          end do
          given = given + per
       end if
       at = at + per
    end do

    r = this % evaluated_over(stored, design)
    if (this % varied > 0) r = partial(r, n)

  end function expression_at_instant

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

    allocate(v(size(this % kind)))

    do i = 1, size(this % kind)
       select case (this % kind(i))
       case (VERTEX_LEAF)
          if (this % position(i) == ARGUMENT_STATE) then
             ! the tuple lists every field's components together; the
             ! components run along time first, then along space
             at = this % stored_at(this % field(i), this % along(i), this % order(i))
             if (at > ubound(q, 1)) then
                error stop 'operation_expression: the state stores the component read'
             end if
             v(i) = q(at)
          else
             v(i) = nu
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
          case (ROOT);        v(i) = sqrt(v(this % first(i)))
          case default
             error stop 'operation_expression: the function is one of those defined'
          end select
       case default
          error stop 'operation_expression: the vertex kind is one of those defined'
       end select
    end do

    r = v(size(v))

  end function evaluated_over

  !===================================================================!
  ! The highest state component read; minus one when none is.
  !===================================================================!

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
    do i = 1, size(this % kind)
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
    do i = 1, size(this % kind)
       if (this % kind(i) /= VERTEX_LEAF .or. this % position(i) /= ARGUMENT_STATE) cycle
       if (this % field(i) > this % fields - this % multipliers) cycle
       at = this % component_at(this % along(i), this % order(i), this % field(i))
       read(at) = .true.
    end do
    components = pack([(at, at = 0, size(read) - 1)], read)

  end subroutine read_components

  pure integer function num_vertices(this)

    class(expression), intent(in) :: this

    num_vertices = size(this % kind)

  end function num_vertices

  !===================================================================!
  ! The degree of the equation, declared when the rule is stated,
  ! with the name the rule reports and the highest exact degree: the
  ! width the subset masks can index. A degree below one stops the
  ! program: there is no highest derivative then.
  !===================================================================!

  subroutine declare_degree(this, degrees, label)

    class(expression)     , intent(inout) :: this
    integer               , intent(in)    :: degrees(:)
    character(len=*)      , intent(in), optional :: label

    if (size(degrees) < 1) then
       error stop 'operation_expression: a state is declared over one coordinate at least'
    end if
    if (degrees(FIRST_COORDINATE) < 1) then
       error stop 'operation_expression: the degree of the equation is positive'
    end if
    if (any(degrees < 0)) then
       error stop 'operation_expression: a degree along a coordinate is not negative'
    end if

    this % degrees = degrees
    if (this % fields < 1) then
       this % fields = 1
       this % field_degree = [degrees(FIRST_COORDINATE)]
    end if
    call this % declare_arguments(2, [ &
         & contract(FIELD_REAL, this % num_components()), &
         & contract(FIELD_REAL, 1) ], label=label, &
         & max_degree=max_subset_width() - merge(1, 0, this % varied > 0))

  end subroutine declare_degree

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

    if (field < 1 .or. field > this % fields) then
       error stop 'operation_expression: a component names a field the rule is stated over'
    end if
    if (field > this % fields - this % multipliers) then
       error stop 'operation_expression: a multiplier is not stored in the tuple'
    end if
    at = 0
    do f = 1, field - 1
       at = at + this % components_per_field(f)
    end do

  end function offset_of_field

  !===================================================================!
  ! The same component in the tuple the rule evaluates, which stores
  ! every field including the multipliers.
  !===================================================================!

  pure integer function stored_at(this, field, coordinate, order) result(at)

    class(expression), intent(in) :: this
    integer          , intent(in) :: field, coordinate, order

    integer :: f

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

    if (.not. this % declared()) then
       error stop 'operation_expression: the law is stated before its components are read'
    end if
    if (coordinate < FIRST_COORDINATE .or. coordinate > this % num_coordinates()) then
       error stop 'operation_expression: a component names a declared coordinate'
    end if
    if (coordinate > FIRST_COORDINATE .and. order < 1) then
       error stop 'operation_expression: a component away from the first coordinate names a positive order'
    end if
    if (order < 0 .or. order > merge(this % field_degree(field), this % degrees(coordinate), &
         & coordinate == FIRST_COORDINATE)) then
       error stop 'operation_expression: a component names an order declared on the coordinate'
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

    if (.not. this % declared()) then
       error stop 'operation_expression: the law is stated before its component count is read'
    end if

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
    do f = 1, this % fields - this % multipliers
       num_components = num_components + this % components_per_field(f)
    end do

  end function num_components

  pure integer function num_fields(this)

    class(expression), intent(in) :: this

    num_fields = this % fields

  end function num_fields

  !===================================================================!
  ! The declared degree along a coordinate: the first field's along
  ! the first coordinate, every field's along a later one.
  !===================================================================!

  pure integer function degree_along(this, coordinate) result(degree)

    class(expression), intent(in) :: this
    integer          , intent(in) :: coordinate

    if (.not. this % declared()) then
       error stop 'operation_expression: the law is stated before its degrees are read'
    end if
    if (coordinate < FIRST_COORDINATE .or. coordinate > size(this % degrees)) then
       error stop 'operation_expression: a degree names a declared coordinate'
    end if
    degree = this % degrees(coordinate)

  end function degree_along

  pure integer function num_multipliers(this)

    class(expression), intent(in) :: this

    num_multipliers = this % multipliers

  end function num_multipliers

  !===================================================================!
  ! A field's degree along the first coordinate.
  !===================================================================!

  pure integer function degree_of_field(this, field) result(degree)

    class(expression), intent(in) :: this
    integer          , intent(in) :: field

    if (.not. this % declared()) then
       error stop 'operation_expression: the law is stated before a field''s degree is read'
    end if
    if (field < 1 .or. field > this % fields) then
       error stop 'operation_expression: a degree names a field the rule is stated over'
    end if
    degree = this % field_degree(field)

  end function degree_of_field

  pure logical function expression_declared(this) result(is_declared)

    class(expression), intent(in) :: this

    is_declared = allocated(this % degrees)

  end function expression_declared

  pure integer function num_coordinates(this)

    class(expression), intent(in) :: this

    if (.not. this % declared()) then
       error stop 'operation_expression: the law is stated before its coordinate count is read'
    end if
    num_coordinates = size(this % degrees)

  end function num_coordinates

  pure integer function equation_degree(this)

    class(expression)     , intent(in) :: this

    if (.not. this % declared()) then
       error stop 'operation_expression: the law is stated before its equation degree is read'
    end if
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

    nd = this % num_components()

    if (size(q) /= input_graph % num_vertices() * nd) then
       error stop 'operation_expression: the state stores one component per law component per point'
    end if
    if (size(nu) /= input_graph % num_vertices()) then
       error stop 'operation_expression: the design stores one value per point'
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

    call this % require_variations(variations)
    call seeded_argument(this, inputs, variations, ARGUMENT_STATE , q , on_state)
    call seeded_argument(this, inputs, variations, ARGUMENT_DESIGN, nu, on_design)
    if (on_state + on_design < size(variations)) then
       error stop 'operation_expression: a variation names the state or the design'
    end if
    call evaluated(this, input_graph, q, nu, output)

  end subroutine expression_partial_action

end module operation_expression
