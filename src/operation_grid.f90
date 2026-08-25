!=====================================================================!
! The time grid: the operation that produces a partition of a
! duration into steps.
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
  use field_calculus        , only : field
  use graph_fractal         , only : graph
  use field_stored          , only : stored_field
  use util_derivative_terms , only : derivative_terms, value, mixed_partial, &
       & max_subset_width, operator(+), operator(*), operator(/)

  implicit none

  private
  public :: grid, uniform_grid, random_grid, designed_grid

  type, abstract, extends(operation) :: grid

     real(dp), private :: span = 1.0_dp

   contains

     procedure(grid_weight_interface), deferred :: weight_of

     procedure :: duration
     procedure :: apply          => grid_apply
     procedure :: max_degree     => grid_max_degree
     procedure :: partial_action => grid_partial_action

  end type grid

  abstract interface

     !----------------------------------------------------------------!
     ! The unnormalised weight of the step ending at instant k, of n.
     !----------------------------------------------------------------!

     pure function grid_weight_interface(this, design, k, n) result(w)
       import :: grid, derivative_terms
       class(grid)           , intent(in) :: this
       type(derivative_terms), intent(in) :: design(:)
       integer               , intent(in) :: k, n
       type(derivative_terms) :: w
     end function grid_weight_interface

  end interface

  type, extends(grid) :: uniform_grid
   contains
     procedure :: name      => uniform_name
     procedure :: weight_of => uniform_weight
  end type uniform_grid

  type, extends(grid) :: random_grid
     integer, private :: seed = 1
   contains
     procedure :: name      => random_name
     procedure :: weight_of => random_weight
  end type random_grid

  type, extends(grid) :: designed_grid
   contains
     procedure :: name      => designed_name
     procedure :: weight_of => designed_weight
  end type designed_grid

  interface uniform_grid
     module procedure create_uniform
  end interface uniform_grid

  interface random_grid
     module procedure create_random
  end interface random_grid

  interface designed_grid
     module procedure create_designed
  end interface designed_grid

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
    type(uniform_grid) :: this

    call require_span(span)
    this % span = span
    call this % declare_arguments(1)

  end function create_uniform

  function create_random(span, seed) result(this)

    real(dp), intent(in) :: span
    integer , intent(in) :: seed
    type(random_grid) :: this

    call require_span(span)
    this % span = span
    this % seed = seed
    call this % declare_arguments(1)

  end function create_random

  function create_designed(span) result(this)

    real(dp), intent(in) :: span
    type(designed_grid) :: this

    call require_span(span)
    this % span = span
    call this % declare_arguments(1)

  end function create_designed

  pure real(dp) function duration(this)

    class(grid), intent(in) :: this

    duration = this % span

  end function duration

  pure function uniform_name(this) result(name)

    class(uniform_grid), intent(in) :: this
    character(len=:), allocatable :: name

    associate (u1 => this); end associate
    name = 'uniform grid'

  end function uniform_name

  pure function random_name(this) result(name)

    class(random_grid), intent(in) :: this
    character(len=:), allocatable :: name

    associate (u1 => this); end associate
    name = 'random grid'

  end function random_name

  pure function designed_name(this) result(name)

    class(designed_grid), intent(in) :: this
    character(len=:), allocatable :: name

    associate (u1 => this); end associate
    name = 'designed grid'

  end function designed_name

  !===================================================================!
  ! THE THREE RULES.
  !===================================================================!

  pure function uniform_weight(this, design, k, n) result(w)

    class(uniform_grid)   , intent(in) :: this
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

    class(random_grid)    , intent(in) :: this
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

    class(designed_grid)  , intent(in) :: this
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

  subroutine partitioned(this, design, n, dt)

    class(grid)           , intent(in) :: this
    type(derivative_terms), intent(in) :: design(:)
    integer               , intent(in) :: n
    type(derivative_terms), allocatable, intent(out) :: dt(:)

    type(derivative_terms) :: total
    integer :: k

    if (n < 2) then
       error stop 'operation_grid: a partition needs two instants'
    end if

    allocate(dt(n))
    dt(1) = derivative_terms(0.0_dp, design(1))
    total = derivative_terms(0.0_dp, design(1))

    do k = 2, n
       dt(k) = this % weight_of(design, k, n)
       if (value(dt(k)) <= 0.0_dp) then
          error stop 'operation_grid: every step weight is positive'
       end if
       total = total + dt(k)
    end do

    do k = 2, n
       dt(k) = this % duration() * (dt(k) / total)
    end do

  end subroutine partitioned

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
