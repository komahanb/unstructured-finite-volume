!=====================================================================!
! A rule stated at one instant of a time hierarchy, held as data: a
! graph whose vertices are typed operators and whose edges are the
! reads between them.
!
!      input graph    the instants
!      input 1        the state, one flat field over instant and
!                     degree: the component of degree d at instant k
!                     is held at (k-1)(N+1) + d + 1
!      input 2        the design, one value per instant, so a design
!                     that varies in time needs no other shape
!      output         one value per instant
!
! A leaf reads a component of an argument - the state's component of
! degree d, or the design - or holds a constant. Every other vertex
! is an arithmetic operation, a power, or an elementary function of
! the vertices it reads. The vertices are stored in evaluation order,
! every read before its reader, and the root is the last one.
!
! A rule is built by the intrinsic operators on expressions, so it
! reads as it is written and the compiler checks it:
!
!      q  = unknown()
!      nu = design()
!      r  = stated(derivative(q, 2) - nu * (1.0_dp - derivative(q, 0)**2) * derivative(q, 1) &
!                + derivative(q, 0), 2, 'van der pol residual')
!
! derivative(q, d) is the component of degree d: along the instants
! the derivatives are unknowns the scheme relates, so the vertex
! selects and does not differentiate.
!
!             THE PARTIALS
!
! The rule is evaluated over derivative_terms by one loop over the
! vertices, so the value and every mixed partial in the directions
! asked for come from one pass and are exact; nothing is
! differentiated symbolically. A variation may name either argument,
! so the Q-partials, the X-partials and the mixed Q-X partials all
! come from one path and to any degree. A Newton block wants the
! Q-partials as a row of numbers, which is one call per degree with a
! direction that is one at that degree and zero elsewhere.
!
!             WHAT IS REFUSED
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
  use operation_action     , only : operation, variation
  use view_directed        , only : directed_graph
  use field_calculus       , only : field
  use graph_fractal        , only : graph
  use field_stored         , only : stored_field
  use util_derivative_terms, only : derivative_terms, mixed_partial, max_subset_width, integer_power

  implicit none

  private
  public :: expression
  public :: unknown, design, constant, derivative, stated
  public :: van_der_pol, van_der_pol_energy, van_der_pol_dissipation
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
  integer, parameter :: ARGUMENT_STATE  = 1
  integer, parameter :: ARGUMENT_DESIGN = 2

  ! the elementary functions
  integer, parameter :: SINE        = 1
  integer, parameter :: COSINE      = 2
  integer, parameter :: EXPONENTIAL = 3
  integer, parameter :: LOGARITHM   = 4
  integer, parameter :: ROOT        = 5

  type, extends(operation) :: expression

     integer, private :: degree = 2

     integer , allocatable, private :: kind(:)
     integer , allocatable, private :: first(:), second(:)    ! the vertices read; 0 if none
     integer , allocatable, private :: position(:)            ! a leaf's argument
     integer , allocatable, private :: order(:)               ! a leaf's component, or a function index
     real(dp), allocatable, private :: coefficient(:)         ! a constant, or an exponent
     character(len=:), allocatable, private :: label

   contains

     procedure :: name           => expression_name
     procedure :: domain         => expression_domain
     procedure :: apply          => expression_apply
     procedure :: max_degree     => expression_max_degree
     procedure :: partial_action => expression_partial_action
     procedure :: at_instant     => expression_at_instant
     procedure :: equation_degree
     procedure :: declare_degree
     procedure :: highest_degree
     procedure :: num_vertices

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

  function vertex(kind, position, order, coefficient) result(this)

    integer , intent(in) :: kind, position, order
    real(dp), intent(in) :: coefficient
    type(expression) :: this

    this % kind        = [kind]
    this % first       = [0]
    this % second      = [0]
    this % position    = [position]
    this % order       = [order]
    this % coefficient = [coefficient]

  end function vertex

  function unknown() result(this)

    type(expression) :: this

    this = vertex(VERTEX_LEAF, ARGUMENT_STATE, 0, 0.0_dp)

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

    this = vertex(VERTEX_LEAF, ARGUMENT_STATE, d, 0.0_dp)

  end function derivative

  !===================================================================!
  ! A rule bound to the degree of the state it reads. A degree below
  ! the highest component read stops the program: the state holds no
  ! such component.
  !===================================================================!

  function stated(rule, degree, label) result(this)

    type(expression), intent(in) :: rule
    integer         , intent(in) :: degree
    character(len=*), intent(in) :: label
    type(expression) :: this

    this = rule
    this % label = label

    if (this % highest_degree() > degree) then
       error stop 'operation_expression: the rule reads a component the state holds'
    end if

    call this % declare_degree(degree)

  end function stated

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
    this % order       = [a % order,       b % order,       0]
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
    this % order       = [a % order,       order]
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

    ! the arithmetic over derivative_terms is read here alone, so the
    ! module's public operators are its own, on expressions
    use util_derivative_terms, only : operator(+), operator(-), operator(*), operator(/), &
         & operator(**), sin, cos, exp, log, sqrt

    class(expression)     , intent(in) :: this
    type(derivative_terms), intent(in) :: q(0:)
    type(derivative_terms), intent(in) :: nu
    type(derivative_terms) :: r

    type(derivative_terms), allocatable :: v(:)
    integer :: i

    allocate(v(size(this % kind)))

    do i = 1, size(this % kind)
       select case (this % kind(i))
       case (VERTEX_LEAF)
          if (this % position(i) == ARGUMENT_STATE) then
             if (this % order(i) > ubound(q, 1)) then
                error stop 'operation_expression: the state holds the component read'
             end if
             v(i) = q(this % order(i))
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

  end function expression_at_instant

  !===================================================================!
  ! The highest state component read; minus one when none is.
  !===================================================================!

  pure integer function highest_degree(this)

    class(expression), intent(in) :: this

    integer :: i

    highest_degree = -1
    do i = 1, size(this % kind)
       if (this % kind(i) == VERTEX_LEAF .and. this % position(i) == ARGUMENT_STATE) then
          highest_degree = max(highest_degree, this % order(i))
       end if
    end do

  end function highest_degree

  pure integer function num_vertices(this)

    class(expression), intent(in) :: this

    num_vertices = size(this % kind)

  end function num_vertices

  pure function expression_name(this) result(name)

    class(expression), intent(in) :: this
    character(len=:), allocatable :: name

    if (allocated(this % label)) then
       name = this % label
    else
       name = 'expression'
    end if

  end function expression_name

  !===================================================================!
  ! The degree of the equation, declared when the rule is stated. A
  ! degree below one stops the program: there is no highest
  ! derivative then.
  !===================================================================!

  subroutine declare_degree(this, degree)

    class(expression)     , intent(inout) :: this
    integer               , intent(in)    :: degree

    if (degree < 1) then
       error stop 'operation_expression: the degree of the equation is positive'
    end if

    this % degree = degree
    call this % declare_arguments(2)

  end subroutine declare_degree

  pure integer function equation_degree(this)

    class(expression)     , intent(in) :: this

    equation_degree = this % degree

  end function equation_degree

  !===================================================================!
  ! THE MACHINERY.
  !
  ! One value per instant.
  !===================================================================!

  subroutine expression_domain(this, input_graph, domain, num_entries)

    class(expression)     , intent(in)  :: this
    class(directed_graph) , intent(in)  :: input_graph
    type(graph)           , intent(out) :: domain
    integer               , intent(out) :: num_entries

    associate (u1 => this); end associate
    domain      = input_graph % vertex_set()
    num_entries = input_graph % num_vertices()

  end subroutine expression_domain

  pure integer function expression_max_degree(this)

    class(expression)     , intent(in) :: this

    associate (u1 => this); end associate
    expression_max_degree = max_subset_width()

  end function expression_max_degree

  !===================================================================!
  ! The state and the design as derivative terms, each direction
  ! seeded on the argument its variation names.
  !===================================================================!

  subroutine seeded(this, input_data, variations, q, nu)

    class(expression)     , intent(in) :: this
    class(field)          , intent(in) :: input_data(:)
    type(variation)       , intent(in) :: variations(:)
    type(derivative_terms), allocatable, intent(out) :: q(:), nu(:)

    real(dp), allocatable :: state(:), design(:), v(:)
    integer :: n, i

    if (size(input_data) < 2) then
       error stop 'operation_expression: the state and the design are given'
    end if

    call input_data(1) % real_vector(state)
    call input_data(2) % real_vector(design)

    n = size(variations)
    call constants(state, n, q)
    call constants(design, n, nu)

    do i = 1, n
       call variations(i) % direction(v)
       if (variations(i) % argument_is(this % argument(1))) then
          call seed(q, i, v)
       else if (variations(i) % argument_is(this % argument(2))) then
          call seed(nu, i, v)
       else
          error stop 'operation_expression: a variation names the state or the design'
       end if
    end do

  end subroutine seeded

  !===================================================================!
  ! A real vector as terms that carry no derivative yet.
  !===================================================================!

  subroutine constants(x, n, terms)

    real(dp), intent(in) :: x(:)
    integer , intent(in) :: n
    type(derivative_terms), allocatable, intent(out) :: terms(:)

    integer :: j

    allocate(terms(size(x)))

    do j = 1, size(x)
       terms(j) = derivative_terms(x(j), n)
    end do

  end subroutine constants

  !===================================================================!
  ! One direction laid into every entry of a quantity.
  !===================================================================!

  subroutine seed(x, i, v)

    type(derivative_terms), intent(inout) :: x(:)
    integer               , intent(in)    :: i
    real(dp)              , intent(in)    :: v(:)

    integer :: j

    if (size(v) /= size(x)) then
       error stop 'operation_expression: a direction has one entry per unknown it varies'
    end if

    do j = 1, size(x)
       call x(j) % set_direction(i, v(j))
    end do

  end subroutine seed

  !===================================================================!
  ! The rule at every instant, answering the coefficient of the full
  ! subset: the value with no directions, the mixed partial with n.
  !===================================================================!

  subroutine evaluated(this, input_graph, q, nu, output)

    class(expression)     , intent(in) :: this
    class(directed_graph) , intent(in) :: input_graph
    type(derivative_terms), intent(in) :: q(:), nu(:)
    class(field), allocatable, intent(inout) :: output

    type(stored_field) :: out
    real(dp), allocatable :: values(:)
    integer :: k, nd, base

    nd = this % degree + 1

    if (size(q) /= input_graph % num_vertices() * nd) then
       error stop 'operation_expression: the state holds one component per degree per instant'
    end if
    if (size(nu) /= input_graph % num_vertices()) then
       error stop 'operation_expression: the design holds one value per instant'
    end if

    allocate(values(input_graph % num_vertices()))

    do k = 1, input_graph % num_vertices()
       base = (k - 1) * nd
       values(k) = mixed_partial(this % at_instant(q(base + 1:base + nd), nu(k)))
    end do

    out = stored_field(this % name(), input_graph % vertex_set(), &
         & input_graph % num_vertices())
    call out % set_real_vector(values)

    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine evaluated

  subroutine expression_apply(this, input_graph, input_data, output)

    class(expression)     , intent(in)       :: this
    class(directed_graph) , intent(in)       :: input_graph
    class(field), intent(in), optional       :: input_data(:)
    class(field), allocatable, intent(inout) :: output

    type(variation), allocatable :: none(:)
    type(derivative_terms), allocatable :: q(:), nu(:)

    if (.not. present(input_data)) then
       error stop 'operation_expression: the state and the design are given'
    end if

    allocate(none(0))
    call seeded(this, input_data, none, q, nu)
    call evaluated(this, input_graph, q, nu, output)

  end subroutine expression_apply

  subroutine expression_partial_action(this, input_graph, input_data, variations, output)

    class(expression)     , intent(in)       :: this
    class(directed_graph) , intent(in)       :: input_graph
    class(field)          , intent(in)       :: input_data(:)
    type(variation)       , intent(in)       :: variations(:)
    class(field), allocatable, intent(inout) :: output

    type(derivative_terms), allocatable :: q(:), nu(:)

    call this % require_owned(variations)

    if (size(variations) > this % max_degree()) then
       error stop 'operation_expression: the requested order is within max_degree'
    end if

    call seeded(this, input_data, variations, q, nu)
    call evaluated(this, input_graph, q, nu, output)

  end subroutine expression_partial_action

  function van_der_pol(degree) result(r)
    integer, intent(in) :: degree
    type(expression) :: r
    type(expression) :: q, nu
    q  = unknown()
    nu = design()
    r = stated(derivative(q, degree) - nu * (1.0_dp - derivative(q, 0)**2) * derivative(q, degree - 1) &
         & + derivative(q, 0), degree, 'van der pol residual')
  end function van_der_pol

  function van_der_pol_energy(degree) result(f)
    integer, intent(in) :: degree
    type(expression) :: f
    type(expression) :: q
    q = unknown()
    f = stated(0.5_dp * (derivative(q, 0)**2 + derivative(q, 1)**2), degree, 'van der pol energy')
  end function van_der_pol_energy

  function van_der_pol_dissipation(degree) result(f)
    integer, intent(in) :: degree
    type(expression) :: f
    type(expression) :: q, nu
    q  = unknown()
    nu = design()
    f = stated(nu * (1.0_dp - derivative(q, 0)**2) * derivative(q, 1) * derivative(q, 1), &
         & degree, 'van der pol dissipation')
  end function van_der_pol_dissipation

end module operation_expression
