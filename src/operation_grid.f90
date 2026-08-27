!=====================================================================!
! The discretisation of a coordinate: the operation that produces a
! finite measure on a duration - points of the axis and the weight
! each carries. A partition and a quadrature are the two readings of
! the one map: a partition's points are the instants and its weights
! the steps, exact on piecewise constants; the gauss kind's points
! are the Legendre nodes and its weights the rule's, exact on
! polynomials to degree 2n - 1.
!
!      input graph    the instants
!      input 1        the design, a real field
!      output         dt, one per instant: the step that ends there
!
! The first instant has no step ending at it and carries zero, so a
! grid over n instants is a partition of the duration into n-1
! steps. Every concretion supplies one unnormalised weight per step
! and nothing else; the weights are then scaled so that the steps sum
! to the duration exactly, which is written here once.
!
!      dt(k)  =  duration * w(k) / sum of the weights
!
! Because the normalisation is part of the operation rather than
! something a caller does afterwards, the partials carry it: a design
! that changes one weight changes every step, and the constraint that
! the steps sum to the duration holds along every direction. A grid
! whose weights are the design is therefore a design on the simplex,
! which is what grid design needs and what a caller normalising by
! hand would lose.
!
!             THE PARTIALS
!
! The weights are computed over derivative_terms, so the value and
! every mixed partial in the design are carried together and are
! exact. A grid whose weights do not read the design has zero
! partials and says so; it does not refuse the question.
!
!             WHAT IS REFUSED
!
! Fewer than two instants, since there is then no step; a weight that
! is not positive, which would give a step of zero or of negative
! length; and a variation on anything but the design.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_grid

  use iso_fortran_env, only : int64
  use util_precision  , only : dp
  use operation_action      , only : operation, variation
  use operation_action, only : emit
  use view_directed         , only : directed_graph
  use view_directed_stored  , only : stored_directed_graph
  use field_calculus        , only : field
  use graph_fractal         , only : graph
  use field_stored          , only : stored_field
  use util_derivative_terms , only : derivative_terms, value, mixed_partial, &
       & max_subset_width, operator(+), operator(*), operator(/)

  implicit none

  private
  public :: grid, uniform_grid, random_grid, designed_grid, fixed_grid, gauss_grid
  public :: partition, partitioned

  ! One representation, the weight rule a case: uniform weighs every
  ! step one, random a bounded deterministic draw, designed the design
  ! itself, fixed a stored partition carried as constants.
  integer, parameter :: GRID_UNIFORM  = 1
  integer, parameter :: GRID_RANDOM   = 2
  integer, parameter :: GRID_DESIGNED = 3
  integer, parameter :: GRID_FIXED    = 4
  integer, parameter :: GRID_GAUSS    = 5

  type, extends(operation) :: grid

     integer , private :: kind = GRID_UNIFORM
     real(dp), private :: span = 1.0_dp
     integer , private :: seed = 1
     real(dp), allocatable, private :: steps(:)

   contains

     procedure :: name           => grid_name
     procedure :: duration
     procedure :: apply          => grid_apply
     procedure :: max_degree     => grid_max_degree
     procedure :: partial_action => grid_partial_action
     procedure :: abscissae
     procedure, private :: weight_of

  end type grid

  interface uniform_grid
     module procedure create_uniform
  end interface uniform_grid

  interface random_grid
     module procedure create_random
  end interface random_grid

  interface designed_grid
     module procedure create_designed
  end interface designed_grid

  interface fixed_grid
     module procedure create_fixed
  end interface fixed_grid

  interface gauss_grid
     module procedure create_gauss
  end interface gauss_grid

  interface partition
     module procedure partition_duration, partition_values
  end interface partition

  interface partitioned
     module procedure partitioned_terms, partitioned_values
  end interface partitioned

contains

  !===================================================================!
  ! A duration that is not positive stops the program: there is no
  ! partition of it.
  !===================================================================!

  subroutine require_span(span)

    real(dp), intent(in) :: span

    if (span <= 0.0_dp) then
       error stop 'operation_grid: a duration is positive'
    end if

  end subroutine require_span

  function create_uniform(span) result(this)

    real(dp), intent(in) :: span
    type(grid) :: this

    call require_span(span)
    this % kind = GRID_UNIFORM
    this % span = span
    call this % declare_arguments(1)

  end function create_uniform

  function create_random(span, seed) result(this)

    real(dp), intent(in) :: span
    integer , intent(in) :: seed
    type(grid) :: this

    call require_span(span)
    this % kind = GRID_RANDOM
    this % span = span
    this % seed = seed
    call this % declare_arguments(1)

  end function create_random

  function create_designed(span) result(this)

    real(dp), intent(in) :: span
    type(grid) :: this

    call require_span(span)
    this % kind = GRID_DESIGNED
    this % span = span
    call this % declare_arguments(1)

  end function create_designed

  !===================================================================!
  ! A grid whose partition is given: the steps, all positive, summing
  ! to the duration they span. A nonpositive step stops the program.
  !===================================================================!

  function create_fixed(steps) result(this)

    real(dp), intent(in) :: steps(:)
    type(grid) :: this

    if (any(steps <= 0.0_dp)) then
       error stop 'operation_grid: every given step is positive'
    end if

    call require_span(sum(steps))
    this % kind  = GRID_FIXED
    this % span  = sum(steps)
    this % steps = steps
    call this % declare_arguments(1)

  end function create_fixed

  !===================================================================!
  ! The quadrature of a duration: every point of the graph applied on
  ! is a node of the Gauss-Legendre rule, and the weights sum to the
  ! duration.
  !===================================================================!

  function create_gauss(span) result(this)

    real(dp), intent(in) :: span
    type(grid) :: this

    call require_span(span)
    this % kind = GRID_GAUSS
    this % span = span
    call this % declare_arguments(1)

  end function create_gauss

  pure real(dp) function duration(this)

    class(grid), intent(in) :: this

    duration = this % span

  end function duration

  pure function grid_name(this) result(name)

    class(grid), intent(in) :: this
    character(len=:), allocatable :: name

    select case (this % kind)
    case (GRID_UNIFORM);  name = 'uniform grid'
    case (GRID_RANDOM);   name = 'random grid'
    case (GRID_DESIGNED); name = 'designed grid'
    case (GRID_GAUSS);    name = 'gauss grid'
    case default;         name = 'fixed grid'
    end select

  end function grid_name

  !===================================================================!
  ! The unnormalised weight of the step ending at instant k, of n:
  ! the rule the kind names.
  !===================================================================!

  pure function weight_of(this, design, k, n) result(w)

    class(grid)           , intent(in) :: this
    type(derivative_terms), intent(in) :: design(:)
    integer               , intent(in) :: k, n
    type(derivative_terms) :: w

    select case (this % kind)
    case (GRID_UNIFORM);  w = uniform_weight(this, design, k, n)
    case (GRID_RANDOM);   w = random_weight(this, design, k, n)
    case (GRID_DESIGNED); w = designed_weight(this, design, k, n)
    case (GRID_GAUSS);    w = gauss_weight(this, design, k, n)
    case default;         w = fixed_weight(this, design, k, n)
    end select

  end function weight_of

  !===================================================================!
  ! The given step ending at instant k: the (k-1)-th stored step, a
  ! constant in the design. An instant outside the partition stops the
  ! program.
  !===================================================================!

  pure function fixed_weight(this, design, k, n) result(w)

    class(grid)           , intent(in) :: this
    type(derivative_terms), intent(in) :: design(:)
    integer               , intent(in) :: k, n

    type(derivative_terms) :: w

    associate (u1 => n); end associate
    if (k - 1 < 1 .or. k - 1 > size(this % steps)) then
       error stop 'operation_grid: the instant is one of the given partition'
    end if
    w = derivative_terms(this % steps(k - 1), design(1))

  end function fixed_weight

  !===================================================================!
  ! The k-th weight of the n-point Gauss-Legendre rule on [0, span],
  ! a constant in the design. The node and weight on [-1, 1] come
  ! from Newton on the Legendre polynomial through its recurrence,
  ! seeded by the Chebyshev estimate; the iteration is stopped at the
  ! arithmetic's spacing of the node.
  !===================================================================!

  pure function gauss_weight(this, design, k, n) result(w)

    class(grid)           , intent(in) :: this
    type(derivative_terms), intent(in) :: design(:)
    integer               , intent(in) :: k, n

    type(derivative_terms) :: w
    real(dp) :: node, weight

    call gauss_node_weight(k, n, this % span, node, weight)
    w = derivative_terms(weight, design(1))

  end function gauss_weight

  pure subroutine gauss_node_weight(k, n, span, node, weight)

    integer , intent(in)  :: k, n
    real(dp), intent(in)  :: span
    real(dp), intent(out) :: node, weight

    real(dp) :: x, p, p_below, p_above, slope, pi
    integer  :: iteration, j

    if (k < 1 .or. k > n) then
       error stop 'operation_grid: the node is one of the rule'
    end if

    pi = acos(-1.0_dp)
    x  = -cos(pi * (real(k, dp) - 0.25_dp) / (real(n, dp) + 0.5_dp))

    do iteration = 1, 64
       ! P_n(x) and its slope by the three-term recurrence
       p_below = 1.0_dp
       p       = x
       do j = 2, n
          p_above = (real(2 * j - 1, dp) * x * p - real(j - 1, dp) * p_below) / real(j, dp)
          p_below = p
          p       = p_above
       end do
       slope = real(n, dp) * (x * p - p_below) / (x * x - 1.0_dp)
       if (abs(p / slope) <= spacing(abs(x) + 1.0_dp)) exit
       x = x - p / slope
    end do

    node   = span * (x + 1.0_dp) / 2.0_dp
    weight = span / ((1.0_dp - x * x) * slope * slope)

  end subroutine gauss_node_weight

  !===================================================================!
  ! THE RULES, one per kind.
  !===================================================================!

  pure function uniform_weight(this, design, k, n) result(w)

    class(grid)           , intent(in) :: this
    type(derivative_terms), intent(in) :: design(:)
    integer               , intent(in) :: k, n
    type(derivative_terms) :: w

    associate (u1 => this, u2 => k, u3 => n); end associate

    w = derivative_terms(1.0_dp, design(1))

  end function uniform_weight

  !===================================================================!
  ! A deterministic weight in one half to three halves, mixed from
  ! the seed and the instant. It is reproducible across runs, which
  ! a drawn number would not be, and it is bounded away from zero, so
  ! no step degenerates.
  !===================================================================!

  pure function random_weight(this, design, k, n) result(w)

    class(grid)           , intent(in) :: this
    type(derivative_terms), intent(in) :: design(:)
    integer               , intent(in) :: k, n
    type(derivative_terms) :: w

    integer(int64), parameter :: modulus    = 2147483648_int64
    integer(int64), parameter :: multiplier = 1103515245_int64
    integer(int64), parameter :: increment  = 12345_int64

    integer(int64) :: x

    associate (u1 => n); end associate

    ! modulo, not mod: mod carries the sign of its first argument, so
    ! a negative seed would leave x below zero and the weight at or
    ! below it. The two agree wherever the seed is positive.
    x = modulo(int(this % seed, int64) * 40503_int64 + int(k, int64) * 65537_int64, modulus)
    x = modulo(multiplier * x + increment, modulus)
    x = modulo(multiplier * x + increment, modulus)

    w = derivative_terms(0.5_dp + real(x, dp) / real(modulus, dp), design(1))

  end function random_weight

  !===================================================================!
  ! The design itself, one entry per step. The step ending at instant
  ! k is designed by entry k-1.
  !===================================================================!

  pure function designed_weight(this, design, k, n) result(w)

    class(grid)           , intent(in) :: this
    type(derivative_terms), intent(in) :: design(:)
    integer               , intent(in) :: k, n
    type(derivative_terms) :: w

    associate (u1 => this, u2 => n); end associate

    if (k - 1 < 1 .or. k - 1 > size(design)) then
       error stop 'operation_grid: the design holds one entry per step'
    end if

    w = design(k - 1)

  end function designed_weight

  pure integer function grid_max_degree(this)

    class(grid), intent(in) :: this

    associate (u1 => this); end associate
    grid_max_degree = max_subset_width()

  end function grid_max_degree

  !===================================================================!
  ! The design as terms, each direction seeded on it. A variation
  ! that names anything else stops the program.
  !===================================================================!

  subroutine seeded(this, input_data, variations, design)

    class(grid)    , intent(in) :: this
    class(field)   , intent(in) :: input_data(:)
    type(variation), intent(in) :: variations(:)
    type(derivative_terms), allocatable, intent(out) :: design(:)

    real(dp), allocatable :: x(:), v(:)
    integer :: n, i, j

    if (size(input_data) < 1) then
       error stop 'operation_grid: the design is given'
    end if

    call input_data(1) % real_vector(x)
    n = size(variations)

    allocate(design(max(size(x), 1)))
    design = derivative_terms(0.0_dp, n)
    do j = 1, size(x)
       design(j) = derivative_terms(x(j), n)
    end do

    do i = 1, n
       if (.not. variations(i) % argument_is(this % argument(1))) then
          error stop 'operation_grid: the steps vary with the design alone'
       end if
       call variations(i) % direction(v)
       do j = 1, size(x)
          call design(j) % set_direction(i, v(j))
       end do
    end do

  end subroutine seeded

  !===================================================================!
  ! The weights, scaled so that the steps sum to the duration. A
  ! weight that is not positive stops the program.
  !===================================================================!

  subroutine partitioned_terms(this, design, n, dt)

    class(grid)           , intent(in) :: this
    type(derivative_terms), intent(in) :: design(:)
    integer               , intent(in) :: n
    type(derivative_terms), allocatable, intent(out) :: dt(:)

    type(derivative_terms) :: total
    integer :: k, first

    ! a partition's first instant carries no step; a quadrature
    ! weighs every point
    first = merge(1, 2, this % kind == GRID_GAUSS)
    if (n < first + 1) then
       error stop 'operation_grid: a partition needs two instants'
    end if

    allocate(dt(n))
    total = derivative_terms(0.0_dp, design(1))
    if (first == 2) dt(1) = derivative_terms(0.0_dp, design(1))

    do k = first, n
       dt(k) = this % weight_of(design, k, n)
       if (value(dt(k)) <= 0.0_dp) then
          error stop 'operation_grid: every step weight is positive'
       end if
       total = total + dt(k)
    end do

    do k = first, n
       dt(k) = this % duration() * (dt(k) / total)
    end do

  end subroutine partitioned_terms

  subroutine partition_duration(duration, n, dt, t)
    real(dp), intent(in) :: duration
    integer , intent(in) :: n
    real(dp), allocatable, intent(out) :: dt(:), t(:)
    call partitioned(uniform_grid(duration), n, dt, t)
  end subroutine partition_duration

  subroutine partition_values(steps, num_instants, design, dt)
    class(grid), intent(in) :: steps
    integer    , intent(in) :: num_instants
    real(dp)   , intent(in) :: design(:)
    real(dp), allocatable, intent(out) :: dt(:)
    type(stored_directed_graph) :: instants
    type(stored_field) :: knobs
    class(field), allocatable :: out
    instants = stored_directed_graph(num_instants, tails=[integer ::], heads=[integer ::])
    knobs    = stored_field('design', instants % vertex_set(), max(size(design), 1))
    call knobs % set_real_vector(padded(design))
    call steps % apply(instants, [knobs], out)
    call out % real_vector(dt)
  end subroutine partition_values

  pure function padded(design) result(x)
    real(dp), intent(in) :: design(:)
    real(dp), allocatable :: x(:)
    if (size(design) == 0) then
       allocate(x(1), source=0.0_dp)
    else
       x = design
    end if
  end function padded

  subroutine partitioned_values(steps, n, dt, t, design)
    class(grid), intent(in) :: steps
    integer    , intent(in) :: n
    real(dp), allocatable, intent(out) :: dt(:), t(:)
    real(dp), intent(in), optional :: design(:)
    type(stored_directed_graph) :: instants
    type(stored_field) :: knobs
    class(field), allocatable :: out
    integer :: k
    instants = stored_directed_graph(n, tails=[integer ::], heads=[integer ::])
    if (present(design)) then
       knobs = stored_field('design', instants % vertex_set(), size(design))
       call knobs % set_real_vector(design)
    else
       knobs = stored_field('design', instants % vertex_set(), 1)
       call knobs % set_real_vector([0.0_dp])
    end if
    call steps % apply(instants, [knobs], out)
    call out % real_vector(dt)
    call steps % abscissae(n, dt, t)
  end subroutine partitioned_values

  !===================================================================!
  ! The points of the axis the weights sit at: a partition's are its
  ! instants, accumulated from the steps; the gauss kind's are its
  ! nodes. One weight per point either way.
  !===================================================================!

  pure subroutine abscissae(this, n, dt, t)

    class(grid), intent(in) :: this
    integer    , intent(in) :: n
    real(dp)   , intent(in) :: dt(:)
    real(dp), allocatable, intent(out) :: t(:)

    real(dp) :: weight
    integer  :: k

    allocate(t(n))
    if (this % kind == GRID_GAUSS) then
       do k = 1, n
          call gauss_node_weight(k, n, this % span, t(k), weight)
       end do
    else
       t(1) = 0.0_dp
       do k = 2, n
          t(k) = t(k - 1) + dt(k)
       end do
    end if

  end subroutine abscissae

  subroutine placed(this, input_graph, dt, output)

    class(grid)           , intent(in) :: this
    class(directed_graph) , intent(in) :: input_graph
    type(derivative_terms), intent(in) :: dt(:)
    class(field), allocatable, intent(inout) :: output

    type(stored_field) :: out
    real(dp), allocatable :: values(:)
    integer :: k

    allocate(values(size(dt)))
    do k = 1, size(dt)
       values(k) = mixed_partial(dt(k))
    end do

    out = stored_field(this % name(), input_graph % vertex_set(), size(dt))
    call out % set_real_vector(values)

    call emit(out, output)

  end subroutine placed

  subroutine grid_apply(this, input_graph, input_data, output)

    class(grid)          , intent(in)        :: this
    class(directed_graph), intent(in)        :: input_graph
    class(field), intent(in), optional       :: input_data(:)
    class(field), allocatable, intent(inout) :: output

    type(variation), allocatable :: none(:)
    type(derivative_terms), allocatable :: design(:), dt(:)

    if (.not. present(input_data)) then
       error stop 'operation_grid: the design is given'
    end if

    allocate(none(0))
    call seeded(this, input_data, none, design)
    call partitioned(this, design, input_graph % num_vertices(), dt)
    call placed(this, input_graph, dt, output)

  end subroutine grid_apply

  subroutine grid_partial_action(this, input_graph, input_data, variations, output)

    class(grid)          , intent(in)        :: this
    class(directed_graph), intent(in)        :: input_graph
    class(field)         , intent(in)        :: input_data(:)
    type(variation)      , intent(in)        :: variations(:)
    class(field), allocatable, intent(inout) :: output

    type(derivative_terms), allocatable :: design(:), dt(:)

    call this % require_owned(variations)

    if (size(variations) > this % max_degree()) then
       error stop 'operation_grid: the requested order is within max_degree'
    end if

    call seeded(this, input_data, variations, design)
    call partitioned(this, design, input_graph % num_vertices(), dt)
    call placed(this, input_graph, dt, output)

  end subroutine grid_partial_action

end module operation_grid
