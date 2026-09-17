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
! not the physics. A fixed value may itself depend on the design -
! the law closes the highest component of the initial tuple - so each
! one states its design rate dh/dnu, and the row's design partial is
! -dh(row)/dnu, not zero. A rate left out is zero, the value being
! data of the problem rather than a function of the design.
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
  use operation_expression , only : expression, euler_lagrange, multiplier, constant
  use operation_expression , only : operator(+), operator(*)
  use operation_domain     , only : continuous_domain, discrete_domain
  use util_derivative_terms, only : max_subset_width, derivative_terms, mixed_partial
  use token_identity       , only : token
  use operation_field      , only : continuous_field, discrete_field
  use operation_manifold   , only : continuous_manifold, discrete_manifold
  use operation_family     , only : family, chain
  use operation_weight     , only : scheme_weight
  use operation_coupling   , only : weights_of
  use view_directed_connectivity, only : connectivity_graph
  use operation_scheme_stencil  , only : derived_constraints
  use operation_finite_difference, only : finite_difference
  use operation_exchange       , only : exchange
  use util_verbosity           , only : verbosity
  use transform_partitioner    , only : partitioner, PARTITION_BREADTH_FIRST
  use relation_partition       , only : partition_relation
  use operation_minimization, only : solve_result
  use operation_newton      , only : newton
  use operation_dense_direct, only : dense_direct
  use operation_gmres       , only : gmres
  use operation_gauss_seidel, only : gauss_seidel

  implicit none

  private
  public :: residual_operator
  public :: continuous_residual, discrete_residual
  public :: operator(+), operator(*)

  type, extends(operation) :: residual_operator

     type(stencil)              , private :: primary_law
     type(stencil), allocatable , private :: connected_law
     type(expression), allocatable, private :: physics
     type(expression), allocatable, private :: rules(:)
     type(stored_directed_graph), private :: points
     type(stored_directed_graph), private :: unknown_vertices
     integer , allocatable, private :: at(:)
     integer , allocatable, private :: fixed_rows(:)
     real(dp), allocatable, private :: fixed(:)
     real(dp), allocatable, private :: fixed_rate(:)
     ! rows stated as affine relations of the state in place of the
     ! residual's own: a fixed row is the case of one unknown
     type(stencil), allocatable, private :: prescribed
     logical      , allocatable, private :: prescribed_row(:)
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
     procedure :: fixed_rates
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

  !===================================================================!
  ! THE CONTINUOUS RESIDUAL: a sum of terms. A term is a list of
  ! equations, each a scalar field, on one manifold; paired with a
  ! multiplier, an unknown of that manifold with one component per
  ! equation, the term contributes to the Lagrangian the product of
  ! the two. The equations of a term on a part of a manifold - a face
  ! at an instant, the time factor - read the unknowns of the parent.
  !===================================================================!

  type :: residual_term

     type(continuous_manifold) :: manifold
     type(continuous_field), allocatable :: equation(:)
     type(continuous_field) :: multiplier
     logical :: paired = .false.

  end type residual_term

  type :: continuous_residual

     type(residual_term), allocatable :: term(:)

   contains

     procedure :: discretize => residual_discretize
     procedure :: num_terms

  end type continuous_residual

  interface continuous_residual
     module procedure create_continuous
  end interface continuous_residual

  interface operator(+)
     module procedure residual_sum
  end interface operator(+)

  interface operator(*)
     module procedure paired_term
  end interface operator(*)

  !===================================================================!
  ! THE DISCRETE RESIDUAL: the equations of the whole manifold as one
  ! Lagrangian rule, the conditions on its parts, the points, and the
  ! derivative approximations along time (a chain of families) and
  ! along space (finite differences of a degree). minimize assembles
  ! one residual_operator per block of the chain over the block's
  ! instants, its stages and the cells, and solves the blocks in
  ! order.
  !
  ! A condition on a face at an instant must be affine in one value
  ! of a parent unknown with coefficient one, u - g = 0: it becomes the
  ! fixed rows of that value at every cell of the instant, and its
  ! multiplier, the reaction, is not stored. A condition on the time
  ! factor must be the integral over the region of one parent unknown:
  ! it becomes, at every moment, the prescribed row that the sum of
  ! the cell volumes times the unknown vanishes, in place of that
  ! unknown's own row at the first cell. Any other condition is
  ! refused.
  !===================================================================!

  type :: discrete_residual

     type(discrete_manifold)   :: points
     type(continuous_manifold) :: manifold
     type(expression)          :: rule
     type(chain)               :: schemes
     logical                   :: with_chain = .false.
     type(finite_difference)   :: differences
     logical                   :: with_differences = .false.
     type(residual_term), allocatable :: condition(:)

   contains

     procedure :: minimize

  end type discrete_residual

contains

  !===================================================================!
  ! Build from the two stencils, the physics, and the points the
  ! physics reads: at(p) is the offset of point p's tuple, whose
  ! components follow in order, and primary(j) the row within it the
  ! j-th rule governs. A residual without a rule is the linear map of
  ! its stencils alone, the frozen linearization. Invalid input: a
  ! fixed row without a value to match, an evaluation point whose
  ! components run past the unknowns, a fixed row outside the
  ! unknowns, or a row count that is not the rule count.
  !===================================================================!

  function create(primary_law, rule, at, unknowns, degrees, primary, fixed_rows, fixed, &
       & connected_law, fixed_rate, prescribed) result(this)

    type(stencil)         , intent(in) :: primary_law
    type(expression)      , intent(in), optional :: rule
    integer               , intent(in) :: at(:), unknowns, degrees, primary(:)
    integer               , intent(in) :: fixed_rows(:)
    real(dp)              , intent(in) :: fixed(:)
    type(stencil)         , intent(in), optional :: connected_law
    real(dp)              , intent(in), optional :: fixed_rate(:)
    type(stencil)         , intent(in), optional :: prescribed
    type(residual_operator) :: this
    type(continuous_domain) :: domain
    integer :: j, components, governed, e
    character(len=250) :: message

    if (size(fixed_rows) /= size(fixed)) then
       write(message,'(a,i0,a,i0)') 'operation_residual: one value is required per fixed &
            &component; size(fixed_rows) = ', size(fixed_rows), ', size(fixed) = ', size(fixed)
       error stop trim(message)
    end if
    if (present(fixed_rate)) then
       if (size(fixed_rate) /= size(fixed)) then
          write(message,'(a,i0,a,i0)') 'operation_residual: one design rate is required per &
               &fixed value; size(fixed_rate) = ', size(fixed_rate), ', size(fixed) = ', size(fixed)
          error stop trim(message)
       end if
    end if
    if (any(fixed_rows < 1) .or. any(fixed_rows > unknowns)) then
       write(message,'(a,i0,a,i0,a,i0)') 'operation_residual: every fixed row must name an &
            &unknown 1..', unknowns, '; fixed_rows range from ', minval(fixed_rows), ' to ', &
            & maxval(fixed_rows)
       error stop trim(message)
    end if
    if (any(at < 0) .or. any(at + degrees > unknowns)) then
       write(message,'(a,i0,a,i0,a,i0)') 'operation_residual: an evaluation point''s degree &
            &components must lie within the unknowns; unknowns = ', unknowns, ', minval(at) = ', &
            & minval(at), ', maxval(at) + degrees = ', maxval(at) + degrees
       error stop trim(message)
    end if
    if (degrees < 1) then
       write(message,'(a,i0)') 'operation_residual: the primary degree count must be positive; &
            &degrees = ', degrees
       error stop trim(message)
    end if
    components = degrees
    governed   = 1
    if (present(rule)) then
       domain     = continuous_domain(rule)
       components = domain % num_components()
       governed   = max(1, rule % num_multipliers())
       if (degrees > components) then
          write(message,'(a,i0,a,i0)') 'operation_residual: the primary degree count must lie &
               &within the law''s component count; degrees = ', degrees, ', num_components = ', &
               & components
          error stop trim(message)
       end if
    end if
    if (size(primary) /= governed) then
       write(message,'(a,i0,a,i0)') 'operation_residual: one governed row is required per rule; &
            &size(primary) = ', size(primary), ', rules = ', governed
       error stop trim(message)
    end if
    if (any(primary < 0) .or. any(primary >= degrees)) then
       write(message,'(a,i0,a,i0,a,i0)') 'operation_residual: a governed row must be one of the &
            &point''s degree components 0..', degrees - 1, '; primary ranges from ', &
            & minval(primary), ' to ', maxval(primary)
       error stop trim(message)
    end if

    this % primary_law = primary_law
    if (present(connected_law)) this % connected_law = connected_law
    ! a Lagrangian's rules are its stationarities, one per multiplier;
    ! a residual without a rule governs no row
    if (present(rule)) then
       this % physics = rule
       if (rule % num_multipliers() > 0) then
          allocate(this % rules(rule % num_multipliers()))
          do j = 1, rule % num_multipliers()
             this % rules(j) = euler_lagrange(rule, j)
          end do
       else
          this % rules = [rule]
       end if
    else
       allocate(this % rules(0))
    end if
    this % at      = at
    this % unknowns = unknowns
    this % degrees  = degrees
    ! the rule states how many components a point stores; the degrees
    ! given are the primary law's portion of them
    this % connected_degrees = components - degrees
    this % primary   = primary
    this % fixed_rows = fixed_rows
    this % fixed      = fixed
    ! a fixed value that is data of the problem has a zero design rate
    allocate(this % fixed_rate(size(fixed)), source=0.0_dp)
    if (present(fixed_rate)) this % fixed_rate = fixed_rate
    this % points = stored_directed_graph(size(at), tails=[integer ::], heads=[integer ::])
    this % unknown_vertices = stored_directed_graph(unknowns, tails=[integer ::], heads=[integer ::])
    ! a prescribed row is one an edge of the prescribed stencil enters;
    ! a row both prescribed and fixed is invalid input
    if (present(prescribed)) then
       this % prescribed = prescribed
       allocate(this % prescribed_row(unknowns), source=.false.)
       do e = 1, prescribed % pattern % num_edges()
          if (prescribed % pattern % edge_head(e) < 1 .or. prescribed % pattern % edge_head(e) > unknowns) then
             error stop 'operation_residual: a prescribed row must name an unknown'
          end if
          this % prescribed_row(prescribed % pattern % edge_head(e)) = .true.
       end do
       do j = 1, size(fixed_rows)
          if (this % prescribed_row(fixed_rows(j))) then
             error stop 'operation_residual: a row is both fixed and prescribed'
          end if
       end do
    end if
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

  ! the rows the residual's own terms do not occupy: fixed and prescribed
  pure function fixed_indicator(this) result(is_fixed)
    class(residual_operator), intent(in) :: this
    logical, allocatable :: is_fixed(:)
    allocate(is_fixed(this % unknowns), source=.false.)
    is_fixed(this % fixed_rows) = .true.
    if (allocated(this % prescribed_row)) is_fixed = is_fixed .or. this % prescribed_row
  end function fixed_indicator

  !===================================================================!
  ! The prescribed rows written into r: the affine relation of the
  ! state, or its partial along a direction, in place of the residual's
  ! own terms on those rows.
  !===================================================================!

  subroutine write_prescribed(this, values, r, along)
    class(residual_operator), intent(in)    :: this
    real(dp)                , intent(in)    :: values(:)
    real(dp)                , intent(inout) :: r(:)
    logical                 , intent(in)    :: along
    type(typed_field_domain) :: states
    type(stored_field) :: given
    class(field), allocatable :: half
    real(dp), allocatable :: y(:)
    if (.not. allocated(this % prescribed)) return
    states = this % state_fields()
    if (along) then
       given = states % direction(values)
       call this % prescribed % partial_action(this % unknown_vertices, this % prescribed % bind([given]), &
            & [variation(this % prescribed % argument(1), given)], half)
    else
       given = states % state(values)
       call this % prescribed % apply(this % unknown_vertices, this % prescribed % bind([given]), half)
    end if
    call half % real_vector(y)
    where (this % prescribed_row) r = y
  end subroutine write_prescribed

  pure subroutine zero_prescribed_rows(this, r)
    class(residual_operator), intent(in)    :: this
    real(dp)                , intent(inout) :: r(:)
    if (.not. allocated(this % prescribed_row)) return
    where (this % prescribed_row) r = 0.0_dp
  end subroutine zero_prescribed_rows

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

  pure function fixed_rates(this) result(rate)
    class(residual_operator), intent(in) :: this
    real(dp), allocatable :: rate(:)
    rate = this % fixed_rate
  end function fixed_rates

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
    ! partial above the first is the physics expression's alone, and
    ! is exact to every order when there is no physics
    residual_max_degree = max_subset_width()
    if (allocated(this % physics)) residual_max_degree = this % physics % max_degree()
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
    if (.not. allocated(this % physics)) then
       error stop 'operation_residual: rule was read from a residual without a rule, the frozen linearization'
    end if
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
    character(len=250) :: message
    if (size(x) /= this % unknowns) then
       write(message,'(a,i0,a,i0)') 'operation_residual: the state must contain one component per &
            &degree per unknown point; size(x) = ', size(x), ', unknowns = ', this % unknowns
       error stop trim(message)
    end if
    if (size(nu) /= size(this % at)) then
       write(message,'(a,i0,a,i0)') 'operation_residual: the design must contain one value per &
            &evaluation point; size(nu) = ', size(nu), ', size(at) = ', size(this % at)
       error stop trim(message)
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
       error stop 'operation_residual: require_host received a graph that is not this residual''s own unknown graph'
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
    character(len=250) :: diagnostic
    if (.not. domain % same_as(given) .or. num_given /= num_values) then
       write(diagnostic,'(a,a,a,l1,a,i0,a,i0)') 'operation_residual: ', message, &
            & '; domain matches = ', domain % same_as(given), ', num_given = ', num_given, &
            & ', num_values = ', num_values
       error stop trim(diagnostic)
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
       error stop 'operation_residual: the design must be defined on the point domain with one value per point'
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

  !===================================================================!
  ! THE DESIGN PARTIAL OF A FIXED ROW. The row states x(row) - h(row)
  ! with h a function of the design, so its partial along a design
  ! direction v is -dh(row)/dnu times v at the point whose tuple contains
  ! the row. A value that is data of the problem has a zero rate and
  ! the row reads zero, as it did before rates were stated.
  !===================================================================!

  pure subroutine design_of_fixed_rows(this, v, r)
    class(residual_operator), intent(in)    :: this
    real(dp)                 , intent(in)    :: v(:)
    real(dp)                 , intent(inout) :: r(:)
    integer :: i
    do i = 1, size(this % fixed_rows)
       r(this % fixed_rows(i)) = -this % fixed_rate(i) * v(point_of(this, this % fixed_rows(i)))
    end do
  end subroutine design_of_fixed_rows

  ! the point whose tuple contains an unknown: the last point whose
  ! offset lies below it, the tuples being placed in order
  pure integer function point_of(this, row) result(p)
    class(residual_operator), intent(in) :: this
    integer                 , intent(in) :: row
    integer :: q
    p = 1
    do q = 1, size(this % at)
       if (this % at(q) < row .and. row <= this % at(q) + this % stride()) p = q
    end do
  end function point_of

  pure logical function any_fixed_rate(this) result(yes)
    class(residual_operator), intent(in) :: this
    yes = any(this % fixed_rate /= 0.0_dp)
  end function any_fixed_rate

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
    character(len=250) :: message
    if (.not. bound_on(inputs, this % argument(1), this % unknown_domain(), this % unknowns)) then
       write(message,'(a,i0)') 'operation_residual: the bound state must be defined on the unknown &
            &domain with one component per degree per unknown point; unknowns = ', this % unknowns
       error stop trim(message)
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
            & 'a direction in the state must be defined on the unknown domain with one value per unknown')
    else if (given % argument_is(this % argument(2))) then
       call require_field(along, size(v), this % design_domain(), size(this % at), &
            & 'a direction in the design must be defined on the point domain with one value per point')
    else
       error stop 'operation_residual: require_direction received a variation naming neither the state nor the design'
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
       error stop 'operation_residual: residual_apply requires inputs (the state and the design), but the optional inputs argument was not given'
    end if
    call require_host(this, input_graph)
    call state_of(this, inputs, x, state)
    call discretized(this, inputs, state, x, r, governing)
    call placed(this, governing, r)
    call accumulate_state(this, x, r)
    call write_prescribed(this, x, r, along=.false.)
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
    call stencil_term(this % primary_law, r)
    if (allocated(this % connected_law)) then
       call stencil_term(this % connected_law, coupled)
       r = r + coupled
    end if
    if (.not. allocated(this % physics)) then
       allocate(governing(size(this % at), 0))
       return
    end if
    call point_inputs(this, inputs, x, point_data)
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
    type(stored_field) :: state
    real(dp), allocatable :: x(:), w(:), column(:), xs(:), nu(:), g(:)
    integer , allocatable :: r(:), c(:), reads(:)
    logical , allocatable :: is_fixed(:)
    integer :: e, d, p, npts, n, num_triples, count, j, k, deg
    real(dp) :: value
    tangent_defined = which == 1
    if (.not. tangent_defined) return
    call require_host(this, input_graph)
    n    = this % unknowns
    npts = size(this % at)
    deg  = this % degrees
    is_fixed = this % fixed_indicator()
    call state_of(this, inputs, x, state)
    if (allocated(this % physics)) then
       if (.not. bound_on(inputs, this % argument(2), this % design_domain(), npts)) then
          error stop 'operation_residual: the design must be defined on the point domain with one value per point'
       end if
       call bound_real_vector(inputs, this % argument(2), nu)
       xs = gathered(this, x)
       allocate(g(0:deg - 1))
    end if
    count = this % primary_law % pattern % num_edges() + npts * this % degrees * size(this % rules) &
         & + size(this % fixed_rows)
    if (allocated(this % connected_law)) count = count + this % connected_law % pattern % num_edges()
    if (allocated(this % prescribed))    count = count + this % prescribed % pattern % num_edges()
    allocate(r(count), c(count), w(count))
    num_triples = 0
    call stencil_triples(this % primary_law, is_fixed, r, c, w, num_triples)
    if (allocated(this % connected_law)) call stencil_triples(this % connected_law, is_fixed, r, c, w, num_triples)
    ! each rule's partials in the components it reads, by one reverse
    ! pass per point: a component no leaf of the rule names has a zero
    ! column
    do j = 1, size(this % rules)
       call this % rules(j) % read_components(reads)
       do p = 1, npts
          if (is_fixed(this % at(p) + this % primary(j) + 1)) cycle
          call this % rules(j) % gradient_at(xs((p - 1) * deg + 1:p * deg), nu(p), value, g)
          do k = 1, size(reads)
             d = reads(k)
             num_triples    = num_triples + 1
             r(num_triples) = this % at(p) + this % primary(j) + 1
             c(num_triples) = this % at(p) + d + 1
             w(num_triples) = g(d)
          end do
       end do
    end do
    do e = 1, size(this % fixed_rows)
       num_triples    = num_triples + 1
       r(num_triples) = this % fixed_rows(e)
       c(num_triples) = this % fixed_rows(e)
       w(num_triples) = 1.0_dp
    end do
    if (allocated(this % prescribed)) then
       call this % prescribed % weights % real_vector(column)
       do e = 1, this % prescribed % pattern % num_edges()
          num_triples    = num_triples + 1
          r(num_triples) = this % prescribed % pattern % edge_head(e)
          c(num_triples) = this % prescribed % pattern % edge_tail(e)
          w(num_triples) = column(e)
       end do
    end if
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
    character(len=250) :: message
    call this % require_owned(variations)
    if (size(variations) < 1 .or. size(variations) > this % max_degree()) then
       write(message,'(a,i0,a,i0)') 'operation_residual: the requested order must lie within &
            &max_degree; size(variations) = ', size(variations), ', max_degree = ', this % max_degree()
       error stop trim(message)
    end if
    call require_host(this, input_graph)
    do i = 1, size(variations)
       call require_direction(this, variations(i))
    end do
    call state_of(this, inputs, x, state)
    if (size(variations) >= 2) then
       ! A mixed or repeated partial of x(row) - h(nu) in the state is
       ! zero at every order above the first. A repeated partial in the
       ! design alone is -d^m h/dnu^m, of which only the first rate is
       ! stated: refuse rather than report zero for a value whose rate
       ! is not zero.
       if (any_fixed_rate(this)) then
          if (all([(variations(i) % argument_is(this % argument(2)), i = 1, size(variations))])) then
             error stop 'operation_residual: a fixed value''s design rate is stated to first order &
                  &only; a repeated design partial of a fixed row is not stated'
          end if
       end if
       call second_tangent(this, inputs, variations, x, governing)
       allocate(r(this % num_unknowns()), source=0.0_dp)
       call placed(this, governing, r)
       call zero_fixed_rows(this, r)
       call zero_prescribed_rows(this, r)
       call placed_output(this, r, output)
       return
    end if
    call variations(1) % direction(v)
    if (variations(1) % argument_is(this % argument(1))) then
       call discretized(this, inputs, state, x, r, governing, v)
       call placed(this, governing, r)
       call accumulate_direction(this, v, r)
       call write_prescribed(this, v, r, along=.true.)
    else
       call design_tangent(this, inputs, variations, x, r, governing)
       call placed(this, governing, r)
       call design_of_fixed_rows(this, v, r)
       call zero_prescribed_rows(this, r)
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
    if (.not. allocated(this % physics)) then
       allocate(governing(size(this % at), 0))
       return
    end if
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
       error stop 'operation_residual: physics_variation received a variation naming neither the state nor the design'
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
    if (.not. allocated(this % physics)) then
       allocate(governing(size(this % at), 0))
       return
    end if
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
    real(dp), allocatable :: fixed(:), rate(:)
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
    allocate(fixed_rows(size(this % fixed_rows)), fixed(size(this % fixed_rows)), &
         &   rate(size(this % fixed_rows)))
    do e = 1, size(this % fixed_rows)
       if (sub_of(this % fixed_rows(e)) == 0) cycle
       ncar             = ncar + 1
       fixed_rows(ncar) = sub_of(this % fixed_rows(e))
       fixed(ncar)      = this % fixed(e)
       rate(ncar)       = this % fixed_rate(e)
    end do

    derived = this % primary_law % restricted(free, values)
    if (allocated(this % prescribed)) then
       error stop 'operation_residual: a residual with prescribed rows is not restricted; the &
            &prescribed relation would lose its unknowns outside the selection'
    end if
    if (allocated(this % connected_law)) then
       secondary = this % connected_law % restricted(free, values)
       sub = residual_operator(derived, this % physics, at, size(free), &
            & this % degrees, this % primary, fixed_rows(1:ncar), fixed(1:ncar), &
            & connected_law=secondary, fixed_rate=rate(1:ncar))
    else
       sub = residual_operator(derived, this % physics, at, size(free), &
            & this % degrees, this % primary, fixed_rows(1:ncar), fixed(1:ncar), &
            & fixed_rate=rate(1:ncar))
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
    character(len=250) :: message

    if (any(free < 1) .or. any(free > this % unknowns)) then
       write(message,'(a,i0,a,i0,a,i0)') 'operation_residual: a constraint must select unknowns &
            &of the residual 1..', this % unknowns, '; free ranges from ', minval(free), ' to ', &
            & maxval(free)
       error stop trim(message)
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
          write(message,'(a,i0,a,i0,a,i0)') 'operation_residual: point ', p, ' is split across &
               &the constraint boundary; ', inside, ' of its ', this % degrees, ' degree components are free'
          error stop trim(message)
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
  ! whole linear map, so the returned residual has no rule and no
  ! fixed rows. transposed states which of
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
    character(len=250) :: message

    if (size(rhs) /= this % unknowns) then
       write(message,'(a,i0,a,i0)') 'operation_residual: one right side value is required per &
            &unknown; size(rhs) = ', size(rhs), ', unknowns = ', this % unknowns
       error stop trim(message)
    end if

    call this % explicit_tangent(input_graph, inputs, 1, r, c, w, tangent_defined)
    if (.not. tangent_defined) then
       error stop 'operation_residual: linearize requires the tangent in the state to be &
            &explicit, but explicit_tangent() reported none'
    end if

    a = stencil(r, c, w, spread(0.0_dp, 1, this % unknowns), 'explicit tangent')
    if (transposed) call a % reverse()
    call a % constants % set_real_vector(-rhs)

    lin = residual_operator(a, at=this % at, unknowns=this % unknowns, degrees=this % degrees, &
         & primary=this % primary(1:1), fixed_rows=[integer ::], fixed=[real(dp) ::])
    ! A = D_Q R maps U to Y = U: the frozen residual is on the same
    ! unknown and point domains as the residual it linearizes
    lin % unknown_vertices = this % unknown_vertices
    lin % points           = this % points
    call lin % versioned(version_number, transposed=transposed)

  end function linearize

  !===================================================================!
  ! A term: the equations on a manifold. Invalid input: an equation
  ! that is not a scalar field, or one on neither the manifold nor
  ! its parent.
  !===================================================================!

  function create_continuous(manifold, equations) result(this)

    type(continuous_manifold), intent(in) :: manifold
    type(continuous_field)   , intent(in) :: equations(:)
    type(continuous_residual) :: this

    integer :: k
    character(len=250) :: message

    allocate(this % term(1))
    this % term(1) % manifold = manifold
    this % term(1) % equation = equations
    do k = 1, size(equations)
       if (equations(k) % num_components() /= 1) then
          write(message,'(a,i0,a,i0)') 'operation_residual: an equation is a scalar field; equation ', &
               & k, ' has components = ', equations(k) % num_components()
          error stop trim(message)
       end if
       if (.not. (manifold % same_support(equations(k) % on) .or. manifold % parent % matches(equations(k) % on))) then
          write(message,'(a,i0,a)') 'operation_residual: equation ', k, ' is a field on neither the &
               &manifold given nor its parent'
          error stop trim(message)
       end if
    end do

  end function create_continuous

  pure integer function num_terms(this)
    class(continuous_residual), intent(in) :: this
    num_terms = size(this % term)
  end function num_terms

  function residual_sum(a, b) result(this)

    type(continuous_residual), intent(in) :: a, b
    type(continuous_residual) :: this

    this % term = [a % term, b % term]

  end function residual_sum

  !===================================================================!
  ! The pairing lambda . g: the multiplier is an unknown of the term's
  ! manifold with one component per equation. Invalid input: a
  ! multiplier of another manifold, or one that is not an unknown, or
  ! a residual of several terms, or a component count other than the
  ! equation count.
  !===================================================================!

  function paired_term(lambda, g) result(this)

    type(continuous_field)   , intent(in) :: lambda
    type(continuous_residual), intent(in) :: g
    type(continuous_residual) :: this

    character(len=250) :: message

    if (size(g % term) /= 1) then
       error stop 'operation_residual: a multiplier pairs with one term; pair each term before summing'
    end if
    if (.not. lambda % is_unknown()) then
       error stop 'operation_residual: a multiplier is an unknown of the term''s manifold'
    end if
    if (.not. g % term(1) % manifold % same_support(lambda % on)) then
       error stop 'operation_residual: the multiplier is an unknown of another manifold than the term''s'
    end if
    if (lambda % num_components() /= size(g % term(1) % equation)) then
       write(message,'(a,i0,a,i0)') 'operation_residual: a multiplier has one component per equation; &
            &components = ', lambda % num_components(), ', equations = ', size(g % term(1) % equation)
       error stop trim(message)
    end if
    this = g
    this % term(1) % multiplier = lambda
    this % term(1) % paired     = .true.

  end function paired_term

  !===================================================================!
  ! THE DISCRETIZATION. Exactly one unpaired term on the whole
  ! manifold gives the rule, one equation per unknown in order; every
  ! other term is a paired condition on a part of that manifold. A
  ! coordinate present without its approximation, or an approximation
  ! given for a coordinate absent, is invalid input.
  !===================================================================!

  function residual_discretize(this, points, time, space) result(image)

    class(continuous_residual), intent(in) :: this
    type(discrete_manifold)   , intent(in) :: points
    type(chain)               , intent(in), optional :: time
    type(finite_difference)   , intent(in), optional :: space
    type(discrete_residual) :: image

    integer :: k, j, whole, num_conditions
    character(len=250) :: message

    whole = 0
    num_conditions = 0
    do k = 1, size(this % term)
       associate (t => this % term(k))
         if (t % manifold % part == 0) then
            if (t % paired) then
               error stop 'operation_residual: the equations on the whole manifold are not paired with a multiplier'
            end if
            if (whole /= 0) then
               error stop 'operation_residual: one term states the equations on the whole manifold'
            end if
            whole = k
         else
            if (.not. t % paired) then
               error stop 'operation_residual: a condition on a part of the manifold is paired with a multiplier'
            end if
            num_conditions = num_conditions + 1
         end if
       end associate
    end do
    if (whole == 0) then
       error stop 'operation_residual: no term states the equations on the whole manifold'
    end if

    associate (t => this % term(whole))
      if (.not. points % identity % matches(t % manifold % identity)) then
         error stop 'operation_residual: the points discretize another manifold than the equations'' own'
      end if
      if (size(t % equation) /= t % manifold % num_unknowns) then
         write(message,'(a,i0,a,i0)') 'operation_residual: one equation is required per unknown of the &
              &manifold; equations = ', size(t % equation), ', unknowns = ', t % manifold % num_unknowns
         error stop trim(message)
      end if
      image % manifold = t % manifold
      ! the Lagrangian of the equations: its stationarity in the j-th
      ! multiplier is the j-th equation, the row of the j-th unknown
      image % rule = multiplier(1) * t % equation(1) % graph(1)
      do j = 2, size(t % equation)
         image % rule = image % rule + multiplier(j) * t % equation(j) % graph(1)
      end do
    end associate

    allocate(image % condition(num_conditions))
    j = 0
    do k = 1, size(this % term)
       if (k == whole) cycle
       j = j + 1
       image % condition(j) = this % term(k)
       if (.not. image % condition(j) % manifold % parent % matches(image % manifold % identity)) then
          error stop 'operation_residual: a condition is stated on a part of another manifold'
       end if
    end do

    image % points = points
    if (points % with_time .neqv. present(time)) then
       error stop 'operation_residual: the derivatives along time are approximated by a chain of &
            &families, given exactly when the manifold has a time coordinate'
    end if
    if (points % with_space .neqv. present(space)) then
       error stop 'operation_residual: the derivatives along space are approximated by finite &
            &differences, given exactly when the manifold has a region'
    end if
    if (present(time)) then
       image % schemes    = time
       image % with_chain = .true.
       if (time % from(size(time % from)) > points % num_instants()) then
          error stop 'operation_residual: a block of the chain begins after the last instant'
       end if
    end if
    if (present(space)) then
       image % differences      = space
       image % with_differences = .true.
    end if

  end function residual_discretize

  !===================================================================!
  ! THE ZERO OF THE DISCRETE RESIDUAL from an estimate: one block of
  ! the chain after another, each a residual_operator over the block's
  ! moments (its instants, and the stages of a one-step family) and
  ! the cells. The solution is a discrete field of the manifold's
  ! unknowns at the instants.
  !===================================================================!

  subroutine minimize(this, estimate, solution)

    class(discrete_residual), intent(in)  :: this
    type(discrete_field)    , intent(in)  :: estimate
    type(discrete_field)    , intent(out) :: solution

    real(dp), allocatable :: tuple(:,:,:), value(:,:)
    type(stencil), allocatable :: rows_along(:,:)
    integer :: stride, nf, ncells, ninst, f, c, k, b, first, last, depth, nb, top, o
    character(len=250) :: message

    stride = this % rule % num_components()
    nf     = this % manifold % num_unknowns
    ncells = this % points % num_cells()
    ninst  = this % points % num_instants()

    if (.not. estimate % on % identity % matches(this % points % identity)) then
       error stop 'operation_residual: the estimate is a field on another manifold than the residual''s'
    end if
    if (estimate % num_components() /= nf) then
       write(message,'(a,i0,a,i0)') 'operation_residual: the estimate has one component per unknown; &
            &components = ', estimate % num_components(), ', unknowns = ', nf
       error stop trim(message)
    end if

    ! the tuples at every instant and cell: the values from the
    ! estimate, every derivative component zero
    allocate(tuple(stride, ncells, ninst), source=0.0_dp)
    do k = 1, ninst
       do c = 1, ncells
          do f = 1, nf
             tuple(this % rule % offset_of_field(f) + 1, c, k) = estimate % value(this % points % point_of(k, c), f)
          end do
       end do
    end do

    ! the finite-difference rows along each space axis and order, once
    if (this % with_differences) then
       top = 0
       do c = 2, this % rule % num_coordinates()
          top = max(top, this % rule % degree_along(c))
       end do
       allocate(rows_along(this % points % dimension, top))
       do c = 1, this % points % dimension
          do o = 1, top
             if (c + 1 > this % rule % num_coordinates()) cycle
             if (o > this % rule % degree_along(c + 1)) cycle
             rows_along(c, o) = this % differences % derivative_rows(this % points % cells, c, o)
          end do
       end do
    else
       allocate(rows_along(0, 0))
    end if

    if (this % with_chain) then
       nb = this % schemes % num_blocks()
       do b = 1, nb
          first = this % schemes % from(b)
          last  = ninst
          if (b < nb) last = this % schemes % from(b + 1) - 1
          depth = 0
          if (b > 1) depth = this % schemes % scheme(b) % history_depth(top_degree(this % rule, nf))
          if (first - depth < 1) then
             write(message,'(a,i0,a,i0,a,i0)') 'operation_residual: block ', b, ' reads ', depth, &
                  & ' instants of history before instant ', first
             error stop trim(message)
          end if
          call solve_block(this, this % schemes % scheme(b), first - depth, last, depth, rows_along, tuple)
       end do
    else
       call solve_point(this, tuple)
    end if

    allocate(value(this % points % num_points(), nf))
    do k = 1, ninst
       do c = 1, ncells
          do f = 1, nf
             value(this % points % point_of(k, c), f) = tuple(this % rule % offset_of_field(f) + 1, c, k)
          end do
       end do
    end do
    solution = discrete_field(this % points, value, this % manifold % unknown_name(1:nf))

  end subroutine minimize

  pure integer function top_degree(rule, nf)
    type(expression), intent(in) :: rule
    integer         , intent(in) :: nf
    integer :: f
    top_degree = 0
    do f = 1, nf
       top_degree = max(top_degree, rule % degree_of_field(f))
    end do
  end function top_degree

  !===================================================================!
  ! Whether a family marches by stages: it does when no derived row of
  ! its own determines a component below the top degree, so that the
  ! tableau's stages determine the step.
  !===================================================================!

  logical function marches_by_stages(scheme, nd)
    type(family), intent(in) :: scheme
    integer     , intent(in) :: nd
    integer, allocatable :: offset(:), degrees(:)
    integer :: d
    marches_by_stages = .true.
    do d = 0, nd - 1
       if (d == scheme % primary_degree(nd - 1)) cycle
       call scheme % row_pattern(d, nd - 1, offset, degrees)
       if (size(offset) > 0) marches_by_stages = .false.
    end do
  end function marches_by_stages

  !===================================================================!
  ! ONE BLOCK: the instants first..last, the first `depth` of them
  ! known from the block before, the stages between consecutive
  ! instants when the family marches by stages. The unknown vector
  ! lists the moments in time order, each moment's cells in order,
  ! each cell's tuple of stride components.
  !===================================================================!

  subroutine solve_block(this, scheme, first, last, depth, rows_along, tuple)

    class(discrete_residual), intent(in)    :: this
    type(family)            , intent(in)    :: scheme
    integer                 , intent(in)    :: first, last, depth
    type(stencil)           , intent(in)    :: rows_along(:,:)
    real(dp)                , intent(inout) :: tuple(:,:,:)

    type(residual_operator) :: rows
    type(stencil) :: relations, gauge
    type(connectivity_graph) :: connectivity
    integer , allocatable :: instant_moment(:), stage_moment(:,:), at(:), primary(:), fixed_rows(:)
    integer , allocatable :: determined(:), source(:), gauge_rows(:), gauge_columns(:)
    real(dp), allocatable :: weight(:), fixed(:), w(:), steps(:), q(:), gauge_weights(:), zeros(:)
    integer :: stride, nf, ncells, n, s, nm, m, k, i, f, c, e, count, nd, v, width, unknowns, npts, p
    integer :: comp, axis, o, j, ne, history_last
    logical :: staged

    stride = this % rule % num_components()
    nf     = this % manifold % num_unknowns
    ncells = this % points % num_cells()
    n      = last - first + 1
    if (verbosity >= 1 .and. this_image() == 1) then
       print '(a,a,a,i0,a,i0,a,es12.4,a,es12.4,a,i0,a)', 'block  ', trim(scheme % name()), '  instants ', first, &
            & '..', last, '  t = ', this % points % instant(first), ' .. ', this % points % instant(last), &
            & '  history ', depth, ' instants'
    end if
    nd     = top_degree(this % rule, nf) + 1
    staged = marches_by_stages(scheme, nd)
    s      = 0
    if (staged) s = scheme % num_stages()

    ! the moments: instant k of the block at instant_moment(k), and
    ! when staged, stage i of the step into instant k+1 at stage_moment(k, i)
    allocate(instant_moment(n), stage_moment(max(n - 1, 1), max(s, 1)), source=0)
    nm = 0
    do k = 1, n
       if (k > 1 .and. staged) then
          do i = 1, s
             nm = nm + 1
             stage_moment(k - 1, i) = nm
          end do
       end if
       nm = nm + 1
       instant_moment(k) = nm
    end do

    width    = stride * ncells
    unknowns = nm * width
    npts     = nm * ncells
    allocate(at(npts))
    do m = 1, nm
       do c = 1, ncells
          at((m - 1) * ncells + c) = (m - 1) * width + (c - 1) * stride
       end do
    end do

    ! THE ROWS ALONG TIME AND SPACE, as triples: the component
    ! determined, the component it reads, the weight
    count = 0
    do f = 1, nf
       if (this % rule % degree_of_field(f) < 1) cycle
       if (staged) then
          connectivity = scheme % stage_connectivity(this % rule % degree_of_field(f) + 1)
          count = count + (n - 1) * connectivity % num_edges() * ncells
       else
          connectivity = scheme % block_connectivity(this % rule % degree_of_field(f) + 1, n)
          count = count + connectivity % num_edges() * ncells
       end if
    end do
    if (this % with_differences) then
       do f = 1, nf
          do c = 2, this % rule % num_coordinates()
             do o = 1, this % rule % degree_along(c)
                count = count + rows_along(c - 1, o) % pattern % num_edges() * nm
             end do
          end do
       end do
    end if
    allocate(determined(count), source(count), weight(count))
    e = 0

    do f = 1, nf
       if (this % rule % degree_of_field(f) < 1) cycle
       if (staged) then
          connectivity = scheme % stage_connectivity(this % rule % degree_of_field(f) + 1)
          do k = 1, n - 1
             steps = spread(this % points % step(first + k), 1, s + 2)
             call weights_of(scheme_weight(scheme), connectivity, steps, w)
             do j = 1, connectivity % num_edges()
                do c = 1, ncells
                   e = e + 1
                   determined(e) = at_of(moment_of(connectivity % edge_head(j), k), c) &
                        & + this % rule % offset_of_field(f) + connectivity % head_degree(j) + 1
                   source(e)     = at_of(moment_of(connectivity % edge_tail(j), k), c) &
                        & + this % rule % offset_of_field(f) + connectivity % tail_degree(j) + 1
                   weight(e)     = w(j)
                end do
             end do
          end do
       else
          connectivity = scheme % block_connectivity(this % rule % degree_of_field(f) + 1, n)
          steps = [(this % points % step(first + k - 1), k = 1, n)]
          call weights_of(scheme_weight(scheme), connectivity, steps, w)
          do j = 1, connectivity % num_edges()
             do c = 1, ncells
                e = e + 1
                determined(e) = at_of(instant_moment(connectivity % edge_head(j)), c) &
                     & + this % rule % offset_of_field(f) + connectivity % head_degree(j) + 1
                source(e)     = at_of(instant_moment(connectivity % edge_tail(j)), c) &
                     & + this % rule % offset_of_field(f) + connectivity % tail_degree(j) + 1
                weight(e)     = w(j)
             end do
          end do
       end if
    end do

    if (this % with_differences) then
       do f = 1, nf
          do c = 2, this % rule % num_coordinates()
             axis = c - 1
             do o = 1, this % rule % degree_along(c)
                comp = this % rule % component_at(c, o, f)
                associate (op => rows_along(axis, o))
                  call op % weights % real_vector(w)
                  ne = op % pattern % num_edges()
                  do m = 1, nm
                     do j = 1, ne
                        e = e + 1
                        determined(e) = at_of(m, op % pattern % edge_head(j)) + comp + 1
                        source(e)     = at_of(m, op % pattern % edge_tail(j)) + this % rule % offset_of_field(f) + 1
                        weight(e)     = w(j)
                     end do
                  end do
                end associate
             end do
          end do
       end do
    end if

    relations = derived_constraints(determined(1:e), source(1:e), weight(1:e), unknowns, &
         & 'rows along time and space')

    ! the row of each equation within a point's tuple: the value row
    ! for a family that determines the derivatives, the top derivative
    ! row for one that determines the values
    allocate(primary(nf))
    do f = 1, nf
       primary(f) = this % rule % offset_of_field(f) &
            & + scheme % primary_degree(this % rule % degree_of_field(f))
    end do

    ! THE FIXED ROWS: every component of the history instants, and the
    ! values a face condition states at the first or the last instant
    call fixed_rows_of(this, first, last, depth, instant_moment, at, tuple, fixed_rows, fixed)

    ! THE GAUGE ROWS: at every moment after the history, the mean of
    ! the integrated unknown over the cells, in place of its own row at
    ! the first cell
    history_last = 0
    if (depth > 0) history_last = instant_moment(depth)
    call gauge_rows_of(this, nm, history_last, at, primary, gauge_rows, gauge_columns, gauge_weights)
    allocate(zeros(unknowns), source=0.0_dp)
    if (size(gauge_rows) > 0) then
       gauge = stencil(gauge_rows, gauge_columns, gauge_weights, zeros, 'gauge')
       rows  = residual_operator(relations, this % rule, at, unknowns, stride, primary, fixed_rows, fixed, &
            & prescribed=gauge)
    else
       rows  = residual_operator(relations, this % rule, at, unknowns, stride, primary, fixed_rows, fixed)
    end if

    ! the estimate: the tuple at each instant, a stage's the instant before it
    allocate(q(unknowns))
    do k = 1, n
       do c = 1, ncells
          q(at_of(instant_moment(k), c) + 1:at_of(instant_moment(k), c) + stride) = tuple(:, c, first + k - 1)
          if (k < n .and. staged) then
             do i = 1, s
                q(at_of(stage_moment(k, i), c) + 1:at_of(stage_moment(k, i), c) + stride) = tuple(:, c, first + k - 1)
             end do
          end if
       end do
    end do

    call solved(this, rows, unknowns, stride, npts, q)

    do k = depth + 1, n
       do c = 1, ncells
          tuple(:, c, first + k - 1) = q(at_of(instant_moment(k), c) + 1:at_of(instant_moment(k), c) + stride)
       end do
    end do

  contains

    pure integer function at_of(m, c)
      integer, intent(in) :: m, c
      at_of = (m - 1) * width + (c - 1) * stride
    end function at_of

    ! the moment of vertex v of a step's stage connectivity: 1 the
    ! instant behind, 2..s+1 the stages, s+2 the instant ahead
    pure integer function moment_of(v, k)
      integer, intent(in) :: v, k
      if (v == 1) then
         moment_of = instant_moment(k)
      else if (v == s + 2) then
         moment_of = instant_moment(k + 1)
      else
         moment_of = stage_moment(k, v - 1)
      end if
    end function moment_of

  end subroutine solve_block

  !===================================================================!
  ! The point manifold: one moment, one cell, no rows along any
  ! coordinate; the equations alone.
  !===================================================================!

  subroutine solve_point(this, tuple)

    class(discrete_residual), intent(in)    :: this
    real(dp)                , intent(inout) :: tuple(:,:,:)

    type(residual_operator) :: rows
    type(stencil) :: relations
    integer , allocatable :: primary(:)
    real(dp), allocatable :: q(:)
    integer :: stride, nf, f

    stride = this % rule % num_components()
    nf     = this % manifold % num_unknowns
    allocate(primary(nf))
    do f = 1, nf
       primary(f) = this % rule % offset_of_field(f)
    end do
    relations = stencil([integer ::], [integer ::], [real(dp) ::], spread(0.0_dp, 1, stride), 'no relation')
    rows = residual_operator(relations, this % rule, [0], stride, stride, primary, [integer ::], [real(dp) ::])
    q = tuple(:, 1, 1)
    call solved(this, rows, stride, stride, 1, q)
    tuple(:, 1, 1) = q

  end subroutine solve_point

  !===================================================================!
  ! Newton on one residual_operator: a direct inner solve on a point
  ! manifold, GMRES preconditioned by block Gauss-Seidel over the
  ! tuples elsewhere. A solve that does not converge stops the program
  ! with the solver's description.
  !===================================================================!

  subroutine solved(this, rows, unknowns, stride, npts, q)

    class(discrete_residual), intent(in)    :: this
    type(residual_operator) , intent(in)    :: rows
    integer                 , intent(in)    :: unknowns, stride, npts
    real(dp)                , intent(inout) :: q(:)

    type(newton)       :: solver
    type(dense_direct) :: direct
    type(gmres)        :: krylov
    type(gauss_seidel) :: sweeps
    type(solve_result) :: outcome
    type(typed_field_domain) :: designs
    type(stored_field) :: design
    real(dp) :: achieved

    if (this % points % with_space) then
       ! THE INNER SOLVE is inexact: each Newton step's linear system
       ! is solved to a relative residual of 1e-3, the constant forcing
       ! term of Dembo, Eisenstat and Steihaug, while Newton's own test
       ! decides the solution; the step count stays at five or six where
       ! a solve to 1e-8 needed the same, at a third of the linear
       ! work. The preconditioner is four block Gauss-Seidel sweeps:
       ! with two, restarted GMRES(60) stagnates on the 32 x 32 box
       ! (the residual unchanged over forty restarts) and never ends,
       ! where four sweeps end it in 26 s and a restart of 120 with two
       ! sweeps in 27 s at twice the basis storage.
       sweeps % block_width    = stride
       sweeps % max_iterations = 4
       krylov % restart        = 60
       krylov % tolerance      = 1.0e-3_dp
       allocate(krylov % preconditioner, source=sweeps)
       allocate(solver % inner, source=krylov)
       ! over several images, the cells are distributed and the
       ! linear solve with them
       if (num_images() > 1) solver % distribution = exchange(owners_of_unknowns(this, unknowns, stride, npts))
    else
       allocate(solver % inner, source=direct)
    end if

    designs = rows % design_fields()
    design  = designs % design(spread(0.0_dp, 1, npts))
    call solver % state(rows, rows % unknown_graph(), rows % unknown_domain(), unknowns, stored_inputs=[design])
    call solver % solve(spread(0.0_dp, 1, unknowns), q, achieved)
    outcome = solver % result()
    if (.not. outcome % converged()) then
       error stop 'operation_residual: minimize did not converge: ' // outcome % description()
    end if

  end subroutine solved

  !===================================================================!
  ! THE OWNER OF EVERY UNKNOWN of a block: the cells are partitioned
  ! over the images by the partitioner's breadth-first rule on the
  ! mesh, each part connected, and every unknown at a cell, at every
  ! moment, is owned by the cell's image.
  !===================================================================!

  function owners_of_unknowns(this, unknowns, stride, npts) result(owner)

    class(discrete_residual), intent(in) :: this
    integer                 , intent(in) :: unknowns, stride, npts
    integer, allocatable :: owner(:)

    type(partitioner) :: cut
    class(directed_graph), allocatable :: part
    type(partition_relation) :: relation
    integer, allocatable :: cell_owner(:)
    integer :: ncells, k, v, p, c, m, at, width

    ncells = this % points % num_cells()
    allocate(cell_owner(ncells), source=0)
    do k = 1, num_images()
       cut = partitioner(PARTITION_BREADTH_FIRST, num_images(), part=k)
       call cut % partition_graph(this % points % cells, part, relation)
       do v = 1, part % num_vertices()
          if (relation % vertex_owner_part(v) == k) cell_owner(relation % global_vertex_index(v)) = k
       end do
    end do
    if (any(cell_owner == 0)) then
       error stop 'operation_residual: the partition of the cells over the images leaves a cell unowned'
    end if

    width = stride * ncells
    allocate(owner(unknowns), source=0)
    do p = 1, npts
       m  = (p - 1) / ncells + 1
       c  = p - (m - 1) * ncells
       at = (m - 1) * width + (c - 1) * stride
       owner(at + 1:at + stride) = cell_owner(c)
    end do
    if (any(owner == 0)) then
       error stop 'operation_residual: an unknown of the block lies at no point of the manifold'
    end if

  end function owners_of_unknowns

  !===================================================================!
  ! THE FIXED ROWS of a block: every component at the history instants,
  ! from the tuples already solved; and, at the first or the last
  ! instant of the manifold, the value each face condition states.
  ! Invalid input: a face condition whose equation is not one parent
  ! value minus a function of the position, with coefficient one.
  !===================================================================!

  subroutine fixed_rows_of(this, first, last, depth, instant_moment, at, tuple, fixed_rows, fixed)

    class(discrete_residual), intent(in)  :: this
    integer                 , intent(in)  :: first, last, depth, instant_moment(:), at(:)
    real(dp)                , intent(in)  :: tuple(:,:,:)
    integer , allocatable   , intent(out) :: fixed_rows(:)
    real(dp), allocatable   , intent(out) :: fixed(:)

    integer :: stride, ncells, ninst, count, k, c, i, j, f, m, e, instant, comp
    real(dp) :: given, slope

    stride = this % rule % num_components()
    ncells = this % points % num_cells()
    ninst  = this % points % num_instants()

    count = depth * ncells * stride
    do j = 1, size(this % condition)
       if (.not. this % condition(j) % manifold % is_face()) cycle
       instant = face_instant(this, this % condition(j) % manifold % face_time)
       if (instant < first .or. instant > last) cycle
       count = count + size(this % condition(j) % equation) * ncells
    end do
    allocate(fixed_rows(count), fixed(count))
    e = 0

    do k = 1, depth
       m = instant_moment(k)
       do c = 1, ncells
          do i = 1, stride
             e = e + 1
             fixed_rows(e) = at((m - 1) * ncells + c) + i
             fixed(e)      = tuple(i, c, first + k - 1)
          end do
       end do
    end do

    do j = 1, size(this % condition)
       if (.not. this % condition(j) % manifold % is_face()) cycle
       instant = face_instant(this, this % condition(j) % manifold % face_time)
       if (instant < first .or. instant > last) cycle
       m = instant_moment(instant - first + 1)
       do i = 1, size(this % condition(j) % equation)
          call stated_value(this % condition(j) % equation(i), f, comp)
          do c = 1, ncells
             call affine_value(this % condition(j) % equation(i), comp, &
                  & this % points % position(this % points % point_of(instant, c)), given, slope)
             if (abs(slope - 1.0_dp) > 1.0e-12_dp) then
                error stop 'operation_residual: a condition on a face is stated as the value of an unknown &
                     &minus a function of the position'
             end if
             e = e + 1
             fixed_rows(e) = at((m - 1) * ncells + c) + this % rule % offset_of_field(f) + 1
             fixed(e)      = given
          end do
       end do
    end do

  end subroutine fixed_rows_of

  ! the instant of a face at a time of the interval
  integer function face_instant(this, time)
    class(discrete_residual), intent(in) :: this
    real(dp)                , intent(in) :: time
    integer :: k
    face_instant = 0
    do k = 1, this % points % num_instants()
       if (abs(this % points % instant(k) - time) <= 1.0e-12_dp * max(1.0_dp, abs(time))) face_instant = k
    end do
    if (face_instant == 0) then
       error stop 'operation_residual: a face condition is stated at a time that is not an instant'
    end if
  end function face_instant

  !===================================================================!
  ! The one parent unknown a face equation reads, at order zero: the
  ! field f and its component in the equation's own tuple. Invalid
  ! input: an equation reading no component, several, or a
  ! derivative.
  !===================================================================!

  subroutine stated_value(equation, f, comp)
    type(continuous_field), intent(in)  :: equation
    integer               , intent(out) :: f, comp
    type(expression) :: g
    integer, allocatable :: reads(:)
    integer :: k
    g = equation % graph(1)
    call g % read_components(reads)
    if (size(reads) /= 1) then
       error stop 'operation_residual: a condition on a face states the value of one unknown'
    end if
    comp = reads(1)
    f = 0
    do k = 1, g % num_fields() - g % num_multipliers()
       if (g % offset_of_field(k) == comp) f = k
    end do
    if (f == 0) then
       error stop 'operation_residual: a condition on a face states the value of an unknown, not a derivative'
    end if
  end subroutine stated_value

  !===================================================================!
  ! The equation g(u) at a position, evaluated at u = 0 and u = 1 on
  ! the component named: given = -g(0), slope = g(1) - g(0), so that a
  ! condition u - h(x) = 0 yields given = h(x), slope = 1.
  !===================================================================!

  subroutine affine_value(equation, comp, position, given, slope)
    type(continuous_field), intent(in)  :: equation
    integer               , intent(in)  :: comp
    real(dp)              , intent(in)  :: position(:)
    real(dp)              , intent(out) :: given, slope
    type(expression) :: g
    type(derivative_terms), allocatable :: q(:)
    type(derivative_terms) :: nu
    real(dp) :: at_zero_value, at_one_value
    g  = equation % graph(1)
    nu = derivative_terms(0.0_dp, 0)
    allocate(q(0:g % num_components() - 1), source=nu)
    at_zero_value = mixed_partial(g % at_instant(q, nu, position))
    q(comp) = derivative_terms(1.0_dp, 0)
    at_one_value  = mixed_partial(g % at_instant(q, nu, position))
    given = -at_zero_value
    slope = at_one_value - at_zero_value
  end subroutine affine_value

  !===================================================================!
  ! THE GAUGE ROWS: for each condition on the time factor, the
  ! integral of one parent unknown over the region, at every moment
  ! after the history: the sum over the cells of the volume times the
  ! unknown's value, in place of the unknown's own row at the first
  ! cell. Invalid input: a condition on the time factor that is not
  ! such an integral.
  !===================================================================!

  subroutine gauge_rows_of(this, nm, history_last, at, primary, rows, columns, weights)

    class(discrete_residual), intent(in)  :: this
    integer                 , intent(in)  :: nm, history_last, at(:), primary(:)
    integer , allocatable   , intent(out) :: rows(:), columns(:)
    real(dp), allocatable   , intent(out) :: weights(:)

    integer :: ncells, count, j, i, m, c, e, f

    ncells = this % points % num_cells()
    count  = 0
    do j = 1, size(this % condition)
       if (.not. this % condition(j) % manifold % is_time_factor()) cycle
       count = count + size(this % condition(j) % equation) * (nm - history_last) * ncells
    end do
    allocate(rows(count), columns(count), weights(count))
    e = 0
    do j = 1, size(this % condition)
       if (.not. this % condition(j) % manifold % is_time_factor()) cycle
       do i = 1, size(this % condition(j) % equation)
          associate (g => this % condition(j) % equation(i))
            if (.not. (g % is_integral() .and. g % is_unknown())) then
               error stop 'operation_residual: a condition on the time factor is the integral of one &
                    &unknown over the region'
            end if
            f = g % index(1)
          end associate
          do m = history_last + 1, nm
             do c = 1, ncells
                e = e + 1
                rows(e)    = at((m - 1) * ncells + 1) + primary(f) + 1
                columns(e) = at((m - 1) * ncells + c) + this % rule % offset_of_field(f) + 1
                weights(e) = this % points % volume(c)
             end do
          end do
       end do
    end do

  end subroutine gauge_rows_of

end module operation_residual
