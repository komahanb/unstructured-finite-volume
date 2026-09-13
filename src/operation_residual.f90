!=====================================================================!
! The residual operator: a physics law composed with one or two
! discretization stencils, summed, over points that each store a
! fixed number of components - and the fixed rows of that state,
! read off as an identity rather than solved for.
!
! LEVEL 3 OF THE STRATIFICATION. Two stencils apply to the whole
! state and are summed; the physics applies only at the points
! (`at`), one rule per row it governs (`primary`), and each rule's
! result is added into the summed stencils on its row. A Lagrangian
! with k multipliers governs k rows, its stationarity in the j-th
! multiplier on primary(j); a plain rule governs one, at every point.
! A row named fixed instead reads x(row) - fixed(row): the identity,
! not the physics.
!
! apply, explicit_tangent and partial_action are composed once here
! from the two stencils' own apply/explicit_tangent/partial_action
! and the physics expression's - a residual never differentiates
! anything itself, it only assembles what the two composed objects
! already differentiate.
!
! A concretion attaches a domain-specific meaning to the points and
! the stencils (a time march, a spatial mesh, or both together) by
! building them before construction; this type does not know which.
!
! THE DOMAINS. The residual owns its unknown domain U (`unknown_graph`,
! one vertex per unknown, no edges) and its design domain P
! (`point_domain`, one vertex per evaluation point). The state Q and
! every direction in the state are fields on U with one value per
! unknown; the design nu and every direction in the design are fields
! on P with one value per point; the residual, its tangents and its
! partials are fields on Y = U. A field of equal length on another
! domain, and a host graph of another identity, are refused: every
! consumer - value, explicit tangent, tangent and adjoint actions,
! higher partials, the frozen linearization and the constrained
! residual - reads one frozen tuple (Q, nu) on U x P.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_residual

  use util_precision   , only : dp
  use operation_action , only : operation, contract, variation
  use operation_action  , only : binding, bound_real_vector, bound_on
  use view_directed     , only : directed_graph
  use graph_fractal     , only : graph
  use view_directed_stored, only : stored_directed_graph
  use field_calculus    , only : field, FIELD_REAL
  use field_stored      , only : stored_field, typed_field_domain
  use operation_stencil  , only : stencil, combine_triples
  use operation_expression , only : expression, constant, stated, euler_lagrange
  use operation_domain     , only : continuous_domain, discrete_domain

  implicit none

  private
  public :: residual_operator

  type, extends(operation) :: residual_operator

     type(stencil)              , private :: primary_law
     type(stencil), allocatable , private :: connected_law
     type(expression)           , private :: physics
     type(expression), allocatable, private :: rules(:)
     type(stored_directed_graph), private :: points
     type(stored_directed_graph), private :: unknown_vertices
     integer , allocatable, private :: at(:)
     integer , allocatable, private :: fixed_rows(:)
     real(dp), allocatable, private :: fixed(:)
     integer                , private :: degrees  = 0
     integer                , private :: connected_degrees = 0
     integer                , private :: unknowns = 0
     integer , allocatable, private :: primary(:)

   contains

     procedure :: name           => residual_name
     procedure :: domain         => residual_domain
     procedure :: apply          => residual_apply
     procedure :: max_degree     => residual_max_degree
     procedure :: defined_at_zero => residual_defined_at_zero
     procedure :: partial_action => residual_partial_action
     procedure :: explicit_tangent => residual_explicit_tangent

     procedure :: num_unknowns
     procedure :: fixed_unknowns
     procedure :: fixed_indicator
     procedure :: fixed_values
     procedure :: num_fixed
     procedure :: first_fixed
     procedure :: stride
     procedure :: num_degrees
     procedure :: primary_degree
     procedure :: primary_row
     procedure :: num_rules
     procedure :: rule_of
     procedure :: num_points
     procedure :: points_at
     procedure :: rule
     procedure :: primary_stencil
     procedure :: has_connected_stencil
     procedure :: connected_stencil
     procedure :: attach_connected_stencil
     procedure :: point_domain
     procedure :: unknown_graph
     procedure :: unknown_domain
     procedure :: design_domain
     procedure :: state_fields
     procedure :: residual_fields
     procedure :: design_fields
     procedure :: frozen_tuple
     procedure :: selected_points
     procedure :: constrain
     procedure :: linearize

  end type residual_operator

  interface residual_operator
     module procedure create
  end interface residual_operator

contains

  !===================================================================!
  ! Build from the two stencils, the physics, and the points the
  ! physics reads: at(p) is the offset of point p's tuple, whose
  ! components follow in order, and primary(j) the row within it the
  ! j-th rule governs. Invalid input: a fixed row without a value to
  ! match, an evaluation point whose components run past the
  ! unknowns, a fixed row outside the unknowns, or a row count that
  ! is not the rule count.
  !===================================================================!

  function create(primary_law, rule, at, unknowns, degrees, primary, fixed_rows, fixed, &
       & connected_law) result(this)

    type(stencil)         , intent(in) :: primary_law
    type(expression)      , intent(in) :: rule
    integer               , intent(in) :: at(:), unknowns, degrees, primary(:)
    integer               , intent(in) :: fixed_rows(:)
    real(dp)              , intent(in) :: fixed(:)
    type(stencil)         , intent(in), optional :: connected_law
    type(residual_operator) :: this
    type(continuous_domain) :: domain
    integer :: j

    if (size(fixed_rows) /= size(fixed)) then
       error stop 'operation_residual: one value per fixed component'
    end if
    if (any(fixed_rows < 1) .or. any(fixed_rows > unknowns)) then
       error stop 'operation_residual: every fixed row names an unknown'
    end if
    if (any(at < 0) .or. any(at + degrees > unknowns)) then
       error stop 'operation_residual: an evaluation point''s degree components lie within the unknowns'
    end if
    domain = continuous_domain(rule)
    if (degrees < 1 .or. degrees > domain % num_components()) then
       error stop 'operation_residual: primary degree count is within the law''s component count'
    end if
    if (size(primary) /= max(1, rule % num_multipliers())) then
       error stop 'operation_residual: one governed row per rule'
    end if
    if (any(primary < 0) .or. any(primary >= degrees)) then
       error stop 'operation_residual: a governed row is one of the point''s degree components'
    end if

    this % primary_law = primary_law
    if (present(connected_law)) this % connected_law = connected_law
    this % physics    = rule
    ! a Lagrangian's rules are its stationarities, one per multiplier
    if (rule % num_multipliers() > 0) then
       allocate(this % rules(rule % num_multipliers()))
       do j = 1, rule % num_multipliers()
          this % rules(j) = euler_lagrange(rule, j)
       end do
    else
       this % rules = [rule]
    end if
    this % at      = at
    this % unknowns = unknowns
    this % degrees  = degrees
    ! the rule states how many components a point stores; the degrees
    ! given are the primary law's portion of them
    this % connected_degrees = domain % num_components() - degrees
    this % primary   = primary
    this % fixed_rows = fixed_rows
    this % fixed      = fixed
    this % points = stored_directed_graph(size(at), tails=[integer ::], heads=[integer ::])
    this % unknown_vertices = stored_directed_graph(unknowns, tails=[integer ::], heads=[integer ::])
    call this % declare_arguments(2, [contract(FIELD_REAL, 1), contract(FIELD_REAL, 1)])

  end function create

  pure integer function num_unknowns(this)
    class(residual_operator), intent(in) :: this
    num_unknowns = this % unknowns
  end function num_unknowns

  pure function fixed_unknowns(this) result(c)
    class(residual_operator), intent(in) :: this
    integer, allocatable :: c(:)
    c = this % fixed_rows
  end function fixed_unknowns

  pure function fixed_indicator(this) result(is_fixed)
    class(residual_operator), intent(in) :: this
    logical, allocatable :: is_fixed(:)
    allocate(is_fixed(this % unknowns), source=.false.)
    is_fixed(this % fixed_rows) = .true.
  end function fixed_indicator

  pure function fixed_values(this) result(h)
    class(residual_operator), intent(in) :: this
    real(dp), allocatable :: h(:)
    h = this % fixed
  end function fixed_values

  pure integer function stride(this)
    class(residual_operator), intent(in) :: this
    stride = this % degrees + this % connected_degrees
  end function stride

  pure integer function num_degrees(this)
    class(residual_operator), intent(in) :: this
    num_degrees = this % degrees
  end function num_degrees

  pure integer function primary_degree(this)
    class(residual_operator), intent(in) :: this
    primary_degree = this % primary(1)
  end function primary_degree

  pure integer function primary_row(this, j)
    class(residual_operator), intent(in) :: this
    integer                 , intent(in) :: j
    primary_row = this % primary(j)
  end function primary_row

  pure integer function num_rules(this)
    class(residual_operator), intent(in) :: this
    num_rules = size(this % rules)
  end function num_rules

  function rule_of(this, j) result(law)
    class(residual_operator), intent(in) :: this
    integer                 , intent(in) :: j
    type(expression) :: law
    law = this % rules(j)
  end function rule_of

  pure integer function num_points(this)
    class(residual_operator), intent(in) :: this
    num_points = size(this % at)
  end function num_points

  pure function points_at(this) result(at)
    class(residual_operator), intent(in) :: this
    integer, allocatable :: at(:)
    at = this % at
  end function points_at

  pure function first_fixed(this) result(x)
    class(residual_operator), intent(in) :: this
    real(dp), allocatable :: x(:)
    allocate(x(this % degrees), source=0.0_dp)
    if (size(this % fixed) >= this % degrees) x = this % fixed(1:this % degrees)
  end function first_fixed

  pure integer function num_fixed(this)
    class(residual_operator), intent(in) :: this
    num_fixed = size(this % fixed_rows)
  end function num_fixed

  pure function residual_name(this) result(name)
    class(residual_operator), intent(in) :: this
    character(len=:), allocatable :: name
    associate (u1 => this); end associate
    name = 'residual operator'
  end function residual_name

  pure integer function residual_max_degree(this)
    class(residual_operator), intent(in) :: this
    ! the discretization stencils are linear in the state, so every
    ! partial above the first is the physics expression's alone
    residual_max_degree = this % physics % max_degree()
  end function residual_max_degree

  !===================================================================!
  ! A LAW NEED NOT BE DEFINED AT THE ZERO STATE. The discretization
  ! stencils are linear in the state, so the zero state lies in the
  ! residual's domain exactly when it lies in the domain of every one
  ! of its rules: the radial oscillator's nu / q**3 excludes it, the
  ! frozen linearization, whose rule is the zero expression, does not.
  !===================================================================!

  pure logical function residual_defined_at_zero(this) result(defined)
    class(residual_operator), intent(in) :: this
    integer :: j
    defined = .true.
    do j = 1, size(this % rules)
       defined = defined .and. this % rules(j) % defined_at_zero()
    end do
  end function residual_defined_at_zero

  type(expression) function rule(this) result(law)
    class(residual_operator), intent(in) :: this
    law = this % physics
  end function rule

  type(stencil) function primary_stencil(this) result(law)
    class(residual_operator), intent(in) :: this
    law = this % primary_law
  end function primary_stencil

  pure logical function has_connected_stencil(this)
    class(residual_operator), intent(in) :: this
    has_connected_stencil = allocated(this % connected_law)
  end function has_connected_stencil

  type(stencil) function connected_stencil(this) result(law)
    class(residual_operator), intent(in) :: this
    law = this % connected_law
  end function connected_stencil

  !===================================================================!
  ! Attach the connected law after construction - a spatial
  ! discretization is laid over an already-built primary march, once
  ! the mesh side of a coupled problem is known.
  !===================================================================!

  subroutine attach_connected_stencil(this, law)
    class(residual_operator), intent(inout) :: this
    type(stencil)            , intent(in)    :: law
    this % connected_law = law
  end subroutine attach_connected_stencil

  type(stored_directed_graph) function point_domain(this) result(points)
    class(residual_operator), intent(in) :: this
    points = this % points
  end function point_domain

  !===================================================================!
  ! THE UNKNOWN DOMAIN U: the graph the residual is evaluated over
  ! (passed as the host of every apply, tangent and linearization),
  ! and its vertex set, the domain of the state, of every direction
  ! in the state and of the residual itself.
  !===================================================================!

  type(stored_directed_graph) function unknown_graph(this) result(unknowns)
    class(residual_operator), intent(in) :: this
    unknowns = this % unknown_vertices
  end function unknown_graph

  type(graph) function unknown_domain(this) result(domain)
    class(residual_operator), intent(in) :: this
    domain = this % unknown_vertices % vertex_set()
  end function unknown_domain

  !===================================================================!
  ! THE DESIGN DOMAIN P: the vertex set of the evaluation points, the
  ! domain of the design and of every direction in it, one value per
  ! point.
  !===================================================================!

  type(graph) function design_domain(this) result(domain)
    class(residual_operator), intent(in) :: this
    domain = this % points % vertex_set()
  end function design_domain

  !===================================================================!
  ! THE TYPED SUPPORTS. The state, a direction in the state and a
  ! tangent are fields on U with one value per unknown; the residual,
  ! a costate and a forcing are fields on Y = U; the design and a
  ! direction in the design are fields on P with one value per point.
  ! Each support reads its extent from the residual's own graph.
  !===================================================================!

  type(typed_field_domain) function state_fields(this) result(fields)
    class(residual_operator), intent(in) :: this
    fields = typed_field_domain(this % unknown_vertices)
  end function state_fields

  type(typed_field_domain) function residual_fields(this) result(fields)
    class(residual_operator), intent(in) :: this
    fields = typed_field_domain(this % unknown_vertices)
  end function residual_fields

  type(typed_field_domain) function design_fields(this) result(fields)
    class(residual_operator), intent(in) :: this
    fields = typed_field_domain(this % points)
  end function design_fields

  !===================================================================!
  ! THE FROZEN TUPLE (Q, nu) on U x P: the state x, one value per
  ! unknown, and the design nu, one value per point. Every consumer -
  ! the value, the explicit tangent, the tangent and adjoint actions,
  ! the higher partials - reads a tuple built here, so each
  ! linearizes the same function at the same point. Invalid input: a
  ! state or a design of another length.
  !===================================================================!

  function frozen_tuple(this, x, nu) result(inputs)
    class(residual_operator), intent(in) :: this
    real(dp)                , intent(in) :: x(:), nu(:)
    type(stored_field) :: inputs(2)
    type(typed_field_domain) :: states, designs
    if (size(x) /= this % unknowns) then
       error stop 'operation_residual: the state contains one component per degree per unknown point'
    end if
    if (size(nu) /= size(this % at)) then
       error stop 'operation_residual: the design contains one value per evaluation point'
    end if
    states  = this % state_fields()
    designs = this % design_fields()
    inputs(1) = states  % state(x)
    inputs(2) = designs % design(nu)
  end function frozen_tuple

  !===================================================================!
  ! The residual's own domain: Y = U, one row per unknown, whatever
  ! graph it is applied on.
  !===================================================================!

  subroutine residual_domain(this, input_graph, domain, num_entries)
    class(residual_operator), intent(in)  :: this
    class(directed_graph)   , intent(in)  :: input_graph
    type(graph)             , intent(out) :: domain
    integer                 , intent(out) :: num_entries
    associate (u1 => input_graph); end associate
    domain      = this % unknown_domain()
    num_entries = this % unknowns
  end subroutine residual_domain

  !===================================================================!
  ! Refuse a host graph other than U: the residual reads and writes
  ! fields on its own unknown domain, and a graph of another identity
  ! with an equal vertex count would type them on a domain they are
  ! not defined on.
  !===================================================================!

  subroutine require_host(this, input_graph)
    class(residual_operator), intent(in) :: this
    class(directed_graph)   , intent(in) :: input_graph
    type(graph) :: given, own
    given = input_graph % vertex_set()
    own   = this % unknown_domain()
    if (.not. own % same_as(given)) then
       error stop 'operation_residual: the residual is applied on its own unknown graph'
    end if
  end subroutine require_host

  !===================================================================!
  ! Refuse a field on another domain than the one stated, or with
  ! another number of values: equal length is not the claim.
  !===================================================================!

  subroutine require_field(given, num_given, domain, num_values, message)
    type(graph)     , intent(in) :: given, domain
    integer         , intent(in) :: num_given, num_values
    character(len=*), intent(in) :: message
    if (.not. domain % same_as(given) .or. num_given /= num_values) then
       error stop 'operation_residual: ' // message
    end if
  end subroutine require_field

  !===================================================================!
  ! The state gathered one point at a time, primary components then
  ! connected ones, in the order the physics reads a point's tuple.
  !===================================================================!

  pure function gathered(this, x) result(y)
    class(residual_operator), intent(in) :: this
    real(dp)                , intent(in) :: x(:)
    real(dp), allocatable :: y(:)
    integer :: p
    allocate(y(size(this % at) * this % stride()))
    do p = 1, size(this % at)
       y((p - 1) * this % stride() + 1:p * this % stride()) = &
            & x(this % at(p) + 1:this % at(p) + this % stride())
    end do
  end function gathered

  !===================================================================!
  ! The physics' inputs at the points: the state gathered per point
  ! and the design read from the frozen tuple. The design is a field
  ! on P with one value per point; another domain or extent is
  ! refused.
  !===================================================================!

  subroutine point_inputs(this, inputs, x, point_data)
    class(residual_operator), intent(in) :: this
    type(binding)            , intent(in) :: inputs(:)
    real(dp)                 , intent(in) :: x(:)
    type(stored_field), allocatable, intent(out) :: point_data(:)
    type(stored_field) :: state, design
    type(typed_field_domain) :: points, designs
    type(continuous_domain) :: continuous
    type(discrete_domain) :: domain
    real(dp), allocatable :: design_values(:)
    if (.not. bound_on(inputs, this % argument(2), this % design_domain(), size(this % at))) then
       error stop 'operation_residual: the design is defined on the point domain with one value per point'
    end if
    call bound_real_vector(inputs, this % argument(2), design_values)
    continuous = continuous_domain(this % physics)
    domain     = continuous % discrete(this % points)
    points  = domain % state_fields()
    designs = domain % design_fields()
    state  = points % state(gathered(this, x))
    design = designs % design(design_values)
    point_data = [state, design]
  end subroutine point_inputs

  pure subroutine placed(this, governing, r)
    class(residual_operator), intent(in)    :: this
    real(dp)                 , intent(in)    :: governing(:,:)
    real(dp)                 , intent(inout) :: r(:)
    integer :: p, j
    do j = 1, size(this % rules)
       do p = 1, size(this % at)
          r(this % at(p) + this % primary(j) + 1) = &
               & r(this % at(p) + this % primary(j) + 1) + governing(p, j)
       end do
    end do
  end subroutine placed

  !===================================================================!
  ! Every rule evaluated at the points, one column per rule: the
  ! value, or the partial along the variations given.
  !===================================================================!

  subroutine governed(this, point_data, governing, variations)
    class(residual_operator), intent(in) :: this
    type(stored_field)      , intent(in) :: point_data(:)
    real(dp), allocatable   , intent(out) :: governing(:,:)
    type(variation)         , intent(in), optional :: variations(:)
    class(field), allocatable :: half
    real(dp), allocatable :: column(:)
    integer :: j
    allocate(governing(size(this % at), size(this % rules)))
    do j = 1, size(this % rules)
       if (present(variations)) then
          call this % rules(j) % partial_action(this % points, this % rules(j) % bind(point_data), &
               & variations, half)
       else
          call this % rules(j) % apply(this % points, this % rules(j) % bind(point_data), half)
       end if
       call half % real_vector(column)
       governing(:, j) = column
    end do
  end subroutine governed

  pure subroutine accumulate_state(this, x, r)
    class(residual_operator), intent(in)    :: this
    real(dp)                 , intent(in)    :: x(:)
    real(dp)                 , intent(inout) :: r(:)
    integer :: i
    do i = 1, size(this % fixed_rows)
       r(this % fixed_rows(i)) = x(this % fixed_rows(i)) - this % fixed(i)
    end do
  end subroutine accumulate_state

  pure subroutine accumulate_direction(this, v, r)
    class(residual_operator), intent(in)    :: this
    real(dp)                 , intent(in)    :: v(:)
    real(dp)                 , intent(inout) :: r(:)
    integer :: i
    do i = 1, size(this % fixed_rows)
       r(this % fixed_rows(i)) = v(this % fixed_rows(i))
    end do
  end subroutine accumulate_direction

  pure subroutine zero_fixed_rows(this, r)
    class(residual_operator), intent(in)    :: this
    real(dp)                 , intent(inout) :: r(:)
    integer :: i
    do i = 1, size(this % fixed_rows)
       r(this % fixed_rows(i)) = 0.0_dp
    end do
  end subroutine zero_fixed_rows

  !===================================================================!
  ! The state read from the frozen tuple: a field on U with one value
  ! per unknown; another domain or extent is refused.
  !===================================================================!

  subroutine state_of(this, inputs, x, state)
    class(residual_operator), intent(in)  :: this
    type(binding)             , intent(in)  :: inputs(:)
    real(dp), allocatable, intent(out) :: x(:)
    type(stored_field)   , intent(out) :: state
    type(typed_field_domain) :: states
    if (.not. bound_on(inputs, this % argument(1), this % unknown_domain(), this % unknowns)) then
       error stop 'operation_residual: the state is defined on the unknown domain with one component per degree &
            &per unknown point'
    end if
    call bound_real_vector(inputs, this % argument(1), x)
    states = this % state_fields()
    state  = states % state(x)
  end subroutine state_of

  !===================================================================!
  ! A direction in the state is a field on U of the state's extent; a
  ! direction in the design is a field on P with one value per point.
  !===================================================================!

  subroutine require_direction(this, given)
    class(residual_operator), intent(in) :: this
    type(variation)         , intent(in) :: given
    type(graph) :: along
    real(dp), allocatable :: v(:)
    along = given % domain()
    call given % direction(v)
    if (given % argument_is(this % argument(1))) then
       call require_field(along, size(v), this % unknown_domain(), this % unknowns, &
            & 'a direction in the state is defined on the unknown domain with one value per unknown')
    else if (given % argument_is(this % argument(2))) then
       call require_field(along, size(v), this % design_domain(), size(this % at), &
            & 'a direction in the design is defined on the point domain with one value per point')
    else
       error stop 'operation_residual: a variation names the state or the design'
    end if
  end subroutine require_direction

  subroutine placed_output(this, r, output)
    class(residual_operator), intent(in) :: this
    real(dp)                 , intent(in) :: r(:)
    class(field), allocatable, intent(inout) :: output
    type(stored_field) :: out
    type(typed_field_domain) :: residuals
    residuals = this % residual_fields()
    out       = residuals % residual(r, this % name())
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)
  end subroutine placed_output

  subroutine residual_apply(this, input_graph, inputs, output)
    class(residual_operator), intent(in)      :: this
    class(directed_graph)    , intent(in)      :: input_graph
    type(binding), intent(in), optional       :: inputs(:)
    class(field), allocatable, intent(inout) :: output
    type(stored_field) :: state
    real(dp), allocatable :: r(:), governing(:,:), x(:)
    if (.not. present(inputs)) then
       error stop 'operation_residual: the state and the design are given'
    end if
    call require_host(this, input_graph)
    call state_of(this, inputs, x, state)
    call discretized(this, inputs, state, x, r, governing)
    call placed(this, governing, r)
    call accumulate_state(this, x, r)
    call placed_output(this, r, output)
  end subroutine residual_apply

  !===================================================================!
  ! THE TWO TERMS OF THE RESIDUAL: the two stencils applied to the
  ! state, summed into r, and the physics at the points, in
  ! governing; or, along a direction v in the state, the partial
  ! action of each in place of its value.
  !===================================================================!

  subroutine discretized(this, inputs, state, x, r, governing, v)
    class(residual_operator), intent(in) :: this
    type(binding)            , intent(in) :: inputs(:)
    type(stored_field)       , intent(in) :: state
    real(dp)                 , intent(in) :: x(:)
    real(dp), allocatable, intent(out) :: r(:), governing(:,:)
    real(dp), intent(in), optional    :: v(:)
    type(stored_field) :: direction
    type(stored_field), allocatable :: point_data(:)
    type(typed_field_domain) :: points
    class(field), allocatable :: half
    real(dp), allocatable :: coupled(:)
    call point_inputs(this, inputs, x, point_data)
    call stencil_term(this % primary_law, r)
    if (allocated(this % connected_law)) then
       call stencil_term(this % connected_law, coupled)
       r = r + coupled
    end if
    if (present(v)) then
       points    = typed_field_domain(this % points, this % stride())
       direction = points % direction(gathered(this, v))
       call governed(this, point_data, governing, [variation(this % physics % argument(1), direction)])
    else
       call governed(this, point_data, governing)
    end if
  contains
    subroutine stencil_term(op, y)
      type(stencil), intent(in) :: op
      real(dp), allocatable, intent(out) :: y(:)
      type(stored_field) :: along
      type(typed_field_domain) :: domain
      if (present(v)) then
         domain = this % state_fields()
         along  = domain % direction(v)
         call op % partial_action(this % unknown_vertices, op % bind([state]), &
              & [variation(op % argument(1), along)], half)
      else
         call op % apply(this % unknown_vertices, op % bind([state]), half)
      end if
      call half % real_vector(y)
    end subroutine stencil_term
  end subroutine discretized

  subroutine residual_explicit_tangent(this, input_graph, inputs, which, &
       & rows, columns, weights, tangent_defined)
    class(residual_operator), intent(in)  :: this
    class(directed_graph)    , intent(in)  :: input_graph
    type(binding)             , intent(in)  :: inputs(:)
    integer              , intent(in)  :: which
    integer , allocatable, intent(out) :: rows(:), columns(:)
    real(dp), allocatable, intent(out) :: weights(:)
    logical              , intent(out) :: tangent_defined
    type(stored_field) :: state, direction
    type(stored_field), allocatable :: point_data(:)
    type(typed_field_domain) :: points
    class(field), allocatable :: out
    real(dp), allocatable :: x(:), w(:), v(:), column(:)
    integer , allocatable :: r(:), c(:), reads(:)
    logical , allocatable :: is_fixed(:)
    integer :: e, d, p, npts, n, num_triples, count, j, k
    tangent_defined = which == 1
    if (.not. tangent_defined) return
    call require_host(this, input_graph)
    n    = this % unknowns
    npts = size(this % at)
    is_fixed = this % fixed_indicator()
    call state_of(this, inputs, x, state)
    call point_inputs(this, inputs, x, point_data)
    count = this % primary_law % pattern % num_edges() + npts * this % degrees * size(this % rules) &
         & + size(this % fixed_rows)
    if (allocated(this % connected_law)) count = count + this % connected_law % pattern % num_edges()
    allocate(r(count), c(count), w(count))
    num_triples = 0
    call stencil_triples(this % primary_law, is_fixed, r, c, w, num_triples)
    if (allocated(this % connected_law)) call stencil_triples(this % connected_law, is_fixed, r, c, w, num_triples)
    ! each rule's partials in the components it reads alone: a
    ! component no leaf of the rule names has a zero column
    allocate(v(npts * this % degrees))
    points = typed_field_domain(this % points, this % degrees)
    do j = 1, size(this % rules)
       call this % rules(j) % read_components(reads)
       do k = 1, size(reads)
          d = reads(k)
          v = 0.0_dp
          do p = 1, npts
             v((p - 1) * this % degrees + d + 1) = 1.0_dp
          end do
          direction = points % direction(v)
          call this % rules(j) % partial_action(this % points, this % rules(j) % bind(point_data), &
               & [variation(this % physics % argument(1), direction)], out)
          call out % real_vector(column)
          do p = 1, npts
             if (is_fixed(this % at(p) + this % primary(j) + 1)) cycle
             num_triples    = num_triples + 1
             r(num_triples) = this % at(p) + this % primary(j) + 1
             c(num_triples) = this % at(p) + d + 1
             w(num_triples) = column(p)
          end do
       end do
    end do
    do e = 1, size(this % fixed_rows)
       num_triples    = num_triples + 1
       r(num_triples) = this % fixed_rows(e)
       c(num_triples) = this % fixed_rows(e)
       w(num_triples) = 1.0_dp
    end do
    call combine_triples(n, n, r(1:num_triples), c(1:num_triples), w(1:num_triples), rows, columns, weights)
  end subroutine residual_explicit_tangent

  subroutine stencil_triples(op, is_fixed, r, c, w, num_triples)
    type(stencil), intent(in)    :: op
    logical      , intent(in)    :: is_fixed(:)
    integer      , intent(inout) :: r(:), c(:)
    real(dp)     , intent(inout) :: w(:)
    integer      , intent(inout) :: num_triples
    real(dp), allocatable :: weights(:)
    integer :: e, row
    call op % weights % real_vector(weights)
    do e = 1, op % pattern % num_edges()
       row = op % pattern % edge_head(e)
       if (is_fixed(row)) cycle
       num_triples    = num_triples + 1
       r(num_triples) = row
       c(num_triples) = op % pattern % edge_tail(e)
       w(num_triples) = weights(e)
    end do
  end subroutine stencil_triples

  subroutine residual_partial_action(this, input_graph, inputs, variations, output)
    class(residual_operator), intent(in)      :: this
    class(directed_graph)    , intent(in)      :: input_graph
    type(binding)             , intent(in)      :: inputs(:)
    type(variation)          , intent(in)      :: variations(:)
    class(field), allocatable, intent(inout) :: output
    type(stored_field) :: state
    real(dp), allocatable :: r(:), governing(:,:), v(:), x(:)
    integer :: i
    call this % require_owned(variations)
    if (size(variations) < 1 .or. size(variations) > this % max_degree()) then
       error stop 'operation_residual: the requested order is within max_degree'
    end if
    call require_host(this, input_graph)
    do i = 1, size(variations)
       call require_direction(this, variations(i))
    end do
    call state_of(this, inputs, x, state)
    if (size(variations) >= 2) then
       call second_tangent(this, inputs, variations, x, governing)
       allocate(r(this % num_unknowns()), source=0.0_dp)
       call placed(this, governing, r)
       call zero_fixed_rows(this, r)
       call placed_output(this, r, output)
       return
    end if
    call variations(1) % direction(v)
    if (variations(1) % argument_is(this % argument(1))) then
       call discretized(this, inputs, state, x, r, governing, v)
       call placed(this, governing, r)
       call accumulate_direction(this, v, r)
    else
       call design_tangent(this, inputs, variations, x, r, governing)
       call placed(this, governing, r)
       call zero_fixed_rows(this, r)
    end if
    call placed_output(this, r, output)
  end subroutine residual_partial_action

  subroutine second_tangent(this, inputs, variations, x, governing)
    class(residual_operator), intent(in) :: this
    type(binding)             , intent(in) :: inputs(:)
    type(variation)          , intent(in) :: variations(:)
    real(dp)                 , intent(in) :: x(:)
    real(dp), allocatable, intent(out) :: governing(:,:)
    type(stored_field), allocatable :: point_data(:)
    type(variation), allocatable :: at_points(:)
    integer :: i
    call point_inputs(this, inputs, x, point_data)
    allocate(at_points(size(variations)))
    do i = 1, size(variations)
       at_points(i) = physics_variation(this, variations(i))
    end do
    call governed(this, point_data, governing, at_points)
  end subroutine second_tangent

  function physics_variation(this, given) result(at_points)
    class(residual_operator), intent(in) :: this
    type(variation)          , intent(in) :: given
    type(variation) :: at_points
    type(stored_field) :: direction
    type(typed_field_domain) :: points
    real(dp), allocatable :: v(:)
    call given % direction(v)
    if (given % argument_is(this % argument(1))) then
       points    = typed_field_domain(this % points, this % stride())
       direction = points % direction(gathered(this, v))
       at_points = variation(this % physics % argument(1), direction)
    else if (given % argument_is(this % argument(2))) then
       at_points = given % with_argument(this % physics % argument(2))
    else
       error stop 'operation_residual: a variation names the state or the design'
    end if
  end function physics_variation

  subroutine design_tangent(this, inputs, variations, x, r, governing)
    class(residual_operator), intent(in) :: this
    type(binding)             , intent(in) :: inputs(:)
    type(variation)          , intent(in) :: variations(:)
    real(dp)                 , intent(in) :: x(:)
    real(dp), allocatable, intent(out) :: r(:), governing(:,:)
    type(stored_field), allocatable :: point_data(:)
    allocate(r(this % num_unknowns()), source=0.0_dp)
    call point_inputs(this, inputs, x, point_data)
    call governed(this, point_data, governing, [variations(1) % with_argument(this % physics % argument(2))])
  end subroutine design_tangent

  !===================================================================!
  ! THE RESIDUAL CONSTRAINED TO free UNKNOWNS OF ITS OWN NUMBERING,
  ! the rest set to values. Ch. 4.6.3 of the dissertation names this
  ! and linearize below the transpose-Jacobian-vector-product
  ! routines. The primary and, if present, connected stencil each
  ! restrict themselves the way a stencil already restricts for
  ! multigrid; a residual adds only what a stencil does not state -
  ! its own evaluation points and fixed rows, renumbered onto free. A
  ! point split across the boundary (some but not all of its degrees
  ! in free) is an invalid constraint.
  !===================================================================!

  function constrain(this, free, values) result(sub)

    class(residual_operator), intent(in) :: this
    integer                 , intent(in) :: free(:)
    real(dp)                , intent(in) :: values(:)
    type(residual_operator) :: sub

    type(stencil), allocatable :: secondary
    type(stencil) :: derived
    integer , allocatable :: sub_of(:), at(:), fixed_rows(:), points(:)
    real(dp), allocatable :: fixed(:)
    integer :: e, p, npts, ncar

    allocate(sub_of(this % unknowns), source=0)
    do e = 1, size(free)
       sub_of(free(e)) = e
    end do

    points = this % selected_points(free)
    npts   = size(points)
    allocate(at(npts))
    do p = 1, npts
       at(p) = sub_of(this % at(points(p)) + 1) - 1
    end do

    ncar = 0
    allocate(fixed_rows(size(this % fixed_rows)), fixed(size(this % fixed_rows)))
    do e = 1, size(this % fixed_rows)
       if (sub_of(this % fixed_rows(e)) == 0) cycle
       ncar             = ncar + 1
       fixed_rows(ncar) = sub_of(this % fixed_rows(e))
       fixed(ncar)      = this % fixed(e)
    end do

    derived = this % primary_law % restricted(free, values)
    if (allocated(this % connected_law)) then
       secondary = this % connected_law % restricted(free, values)
       sub = residual_operator(derived, this % physics, at, size(free), &
            & this % degrees, this % primary, fixed_rows(1:ncar), fixed(1:ncar), connected_law=secondary)
    else
       sub = residual_operator(derived, this % physics, at, size(free), &
            & this % degrees, this % primary, fixed_rows(1:ncar), fixed(1:ncar))
    end if

  end function constrain

  !===================================================================!
  ! THE POINTS A CONSTRAINT RETAINS: the evaluation points whose
  ! degree components all lie in free, in the order of this
  ! residual's points; the constrained residual's j-th point is the
  ! selected_points(j)-th point of this one, and a field on P
  ! restricts to the constrained residual's P' by these indices. A
  ! point split across the boundary (some but not all of its degrees
  ! in free) is an invalid constraint.
  !===================================================================!

  function selected_points(this, free) result(points)

    class(residual_operator), intent(in) :: this
    integer                 , intent(in) :: free(:)
    integer, allocatable :: points(:)

    logical, allocatable :: chosen(:)
    integer, allocatable :: selected(:)
    integer :: e, p, d, inside, npts

    if (any(free < 1) .or. any(free > this % unknowns)) then
       error stop 'operation_residual: a constraint selects unknowns of the residual'
    end if
    allocate(chosen(this % unknowns), source=.false.)
    do e = 1, size(free)
       chosen(free(e)) = .true.
    end do

    npts = 0
    allocate(selected(size(this % at)))
    do p = 1, size(this % at)
       inside = 0
       do d = 1, this % degrees
          if (chosen(this % at(p) + d)) inside = inside + 1
       end do
       if (inside == 0) cycle
       if (inside /= this % degrees) then
          error stop 'operation_residual: a member contains whole points'
       end if
       npts           = npts + 1
       selected(npts) = p
    end do
    points = selected(1:npts)

  end function selected_points

  !===================================================================!
  ! THE EXPLICIT TANGENT, FROZEN INTO A LINEAR RESIDUAL. Ch. 4.6.3 of
  ! the dissertation names this and constrain above the transpose-
  ! Jacobian-vector-product routines. The frozen matrix states the
  ! whole linear map, so the returned residual has no physics of its
  ! own (a zero rule) and no fixed rows. transposed states which of
  ! Jw = rhs or J^Tw = rhs the returned residual's own apply computes;
  ! version_number identifies the returned residual (versioned,
  ! this module's own procedure inherited from operation_action).
  !===================================================================!

  function linearize(this, input_graph, inputs, rhs, transposed, version_number) result(lin)

    class(residual_operator), intent(in) :: this
    class(directed_graph)   , intent(in) :: input_graph
    type(binding)           , intent(in) :: inputs(:)
    real(dp)                , intent(in) :: rhs(:)
    logical                 , intent(in) :: transposed
    integer                 , intent(in) :: version_number
    type(residual_operator) :: lin

    type(stencil) :: a
    integer , allocatable :: r(:), c(:)
    real(dp), allocatable :: w(:)
    logical :: tangent_defined

    if (size(rhs) /= this % unknowns) then
       error stop 'operation_residual: one right side per unknown'
    end if

    call this % explicit_tangent(input_graph, inputs, 1, r, c, w, tangent_defined)
    if (.not. tangent_defined) then
       error stop 'operation_residual: the tangent in the state is explicit'
    end if

    a = stencil(r, c, w, spread(0.0_dp, 1, this % unknowns), 'explicit tangent')
    if (transposed) call a % reverse()
    call a % constants % set_real_vector(-rhs)

    lin = residual_operator(a, stated(constant(0.0_dp), this % degrees - 1, 'zero'), this % at, &
         & this % unknowns, this % degrees, this % primary(1:1), [integer ::], [real(dp) ::])
    ! A = D_Q R maps U to Y = U: the frozen residual is on the same
    ! unknown and point domains as the residual it linearizes
    lin % unknown_vertices = this % unknown_vertices
    lin % points           = this % points
    call lin % versioned(version_number, transposed=transposed)

  end function linearize

end module operation_residual
