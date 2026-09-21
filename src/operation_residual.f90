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
! apply, explicit_jacobian and partial_action are composed once here
! from the two stencils' own apply/explicit_jacobian/partial_action
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
! consumer - value, explicit jacobian, tangent and adjoint actions,
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
  use util_derivative_terms, only : max_subset_width, derivative_terms, mixed_partial, coefficient, symmetric_terms
  use util_derivative_terms, only : value
  use util_derivative_terms, only : operator(*)
  use token_identity       , only : token
  use operation_field      , only : continuous_field, discrete_field, WHOLE, TIME_FACE, TIME_FACTOR, SPACE_FACTOR
  use view_expression      , only : expression_view
  use operation_manifold   , only : continuous_manifold, discrete_manifold
  use operation_family     , only : family, chain
  use operation_weight     , only : scheme_weight
  use operation_coupling   , only : weights_of
  use view_directed_connectivity, only : connectivity_graph
  use operation_scheme_stencil  , only : derived_constraints
  use operation_derivative_approximation, only : derivative_approximation
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
     ! the designs at a point: the design is a vector of one value per
     ! design coordinate at every evaluation point
     integer                , private :: designs = 1

   contains

     procedure :: name           => residual_name
     procedure :: domain         => residual_domain
     procedure :: apply          => residual_apply
     procedure :: max_degree     => residual_max_degree
     procedure :: defined_at_zero => residual_defined_at_zero
     procedure :: partial_action => residual_partial_action
     procedure :: explicit_jacobian => residual_explicit_jacobian

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
  ! THE CONTINUOUS RESIDUAL: a sum of terms, and an objective. A term
  ! is a list of equations, each a scalar field, on one manifold;
  ! paired with a multiplier, an unknown of that manifold with one
  ! component per equation, declared by the term itself, the term
  ! contributes to the Lagrangian the product of the two. The
  ! equations of a term on a part of a manifold - a face at an
  ! instant, the time factor - read the unknowns of the parent. The
  ! objective J, a scalar integral over the manifold, enters the
  ! Lagrangian as a term of its own, L = J + lambda_r . r + lambda . g:
  ! its stationarity in the state is the adjoint equation, and the
  ! multipliers are the sensitivities of J to their equalities.
  !===================================================================!

  type :: residual_term

     type(continuous_manifold) :: manifold
     type(continuous_field), allocatable :: equation(:)
     type(continuous_field) :: multiplier
     logical :: paired = .false.

  end type residual_term

  type :: continuous_residual

     type(residual_term), allocatable :: term(:)
     type(continuous_field), allocatable :: objective

   contains

     procedure :: discretize => residual_discretize
     procedure :: multiplier => residual_multiplier
     procedure :: num_terms

  end type continuous_residual

  interface continuous_residual
     module procedure create_continuous
  end interface continuous_residual

  interface operator(+)
     module procedure residual_sum
     module procedure objective_sum
  end interface operator(+)

  interface operator(*)
     module procedure paired_term
  end interface operator(*)

  !===================================================================!
  ! THE DISCRETE RESIDUAL: the equations of the whole manifold as one
  ! Lagrangian rule, the conditions on its parts, the points, and the
  ! derivative approximations along time (a chain of families) and
  ! along space (a derivative_approximation of an order: finite
  ! differences or finite volumes). minimize assembles
  ! one residual_operator per block of the chain over the block's
  ! instants, its stages and the cells, and solves the blocks in
  ! order.
  !
  ! A condition on a face at an instant is affine in the components
  ! of the parent unknowns' jets, with coefficient one on exactly one
  ! of them, the component it determines: u - g = 0, or lambda_t -
  ! u_t = 0. It becomes, at every cell of the instant, the prescribed
  ! row of that component, the affine relation in place of the
  ! component's own row. A condition on the time factor must be the
  ! integral over the region of one parent unknown: it becomes, at
  ! every moment, the prescribed row that the sum of the cell volumes
  ! times the unknown vanishes, in place of that unknown's own row at
  ! the first cell. Any other condition is refused.
  !
  ! THE ADJOINT. With the equations paired with a multiplier, minimize
  ! follows the forward sweep over the blocks by a reverse sweep: each
  ! block's jacobian at its solution, transposed, determines the
  ! multipliers of the block's rows from the differential of the
  ! objective and the multipliers the later blocks pass back through
  ! the block's history instants. The multiplier of each equation at
  ! each instant is returned under the name the term declared - for
  ! a family marching by stages, the sum over the stages of the step
  ! into the instant of their multipliers, the sensitivity of the
  ! objective to a forcing of the equation over the step - the
  ! multiplier of each face condition, the reaction, under its own
  ! name at the face's points, the multiplier of each condition on
  ! the time factor under its name at every cell of the instant, and
  ! the multiplier of the design condition under its name at every
  ! point. Every multiplier returned is that of a discrete row: the
  ! sensitivity of the objective to a unit perturbation of that row.
  ! The row of an equation at an instant is the equation times the
  ! instant's step, so the multiplier is the continuous adjoint at
  ! the instant times the step; the multiplier of a datum's row is the
  ! sensitivity to the datum itself. A design field is an unknown
  ! with a paired equation
  ! f - f_0 = 0 among the equations of the manifold: the multiplier
  ! of that equation at each instant is the sensitivity of the
  ! objective to the datum there, the functional derivative.
  !===================================================================!

  type :: discrete_residual

     type(discrete_manifold)   :: points
     type(continuous_manifold) :: manifold
     type(expression)          :: rule
     type(chain)               :: schemes
     logical                   :: with_chain = .false.
     class(derivative_approximation), allocatable :: differences
     logical                   :: with_differences = .false.
     type(residual_term), allocatable :: condition(:)
     type(continuous_field), allocatable :: objective
     logical                   :: with_adjoint = .false.
     character(len=32)         :: adjoint_name = ''

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
    if (allocated(this % physics)) this % designs = this % physics % num_designs()
    call this % declare_arguments(2, [contract(FIELD_REAL, 1), contract(FIELD_REAL, this % designs)])

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
    type(stored_field) :: given(1)
    type(variation) :: varied(1)
    class(field), allocatable :: half
    real(dp), allocatable :: y(:)
    if (.not. allocated(this % prescribed)) return
    states = this % state_fields()
    ! the bound field and the variation are array variables, not
    ! array constructors: the compiler in use does not release the
    ! allocatable components of a constructor's temporary
    if (along) then
       given(1)  = states % direction(values)
       varied(1) = variation(this % prescribed % argument(1), given(1))
       call this % prescribed % partial_action(this % unknown_vertices, this % prescribed % bind(given), varied, half)
    else
       given(1) = states % state(values)
       call this % prescribed % apply(this % unknown_vertices, this % prescribed % bind(given), half)
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
    fields = typed_field_domain(this % points, this % designs)
  end function design_fields

  !===================================================================!
  ! THE FROZEN TUPLE (Q, nu) on U x P: the state x, one value per
  ! unknown, and the design nu, one value per point. Every consumer -
  ! the value, the explicit jacobian, the tangent and adjoint actions,
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
    if (size(nu) /= size(this % at) * this % designs) then
       write(message,'(a,i0,a,i0)') 'operation_residual: the design must contain one value per design &
            &per evaluation point; size(nu) = ', size(nu), ', size(at) * designs = ', &
            & size(this % at) * this % designs
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
    if (.not. bound_on(inputs, this % argument(2), this % design_domain(), size(this % at) * this % designs)) then
       error stop 'operation_residual: the design must be defined on the point domain with one value per &
            &design per point'
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
       call require_field(along, size(v), this % design_domain(), size(this % at) * this % designs, &
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
    type(variation) :: varied(1)
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
       varied(1) = variation(this % physics % argument(1), direction)
       call governed(this, point_data, governing, varied)
    else
       call governed(this, point_data, governing)
    end if
  contains
    subroutine stencil_term(op, y)
      type(stencil), intent(in) :: op
      real(dp), allocatable, intent(out) :: y(:)
      type(stored_field) :: along(1), bound(1)
      type(variation) :: varied(1)
      type(typed_field_domain) :: domain
      ! array variables, not array constructors: the compiler in use
      ! does not release the allocatable components of a constructor's
      ! temporary
      bound(1) = state
      if (present(v)) then
         domain    = this % state_fields()
         along(1)  = domain % direction(v)
         varied(1) = variation(op % argument(1), along(1))
         call op % partial_action(this % unknown_vertices, op % bind(bound), varied, half)
      else
         call op % apply(this % unknown_vertices, op % bind(bound), half)
      end if
      call half % real_vector(y)
    end subroutine stencil_term
  end subroutine discretized

  subroutine residual_explicit_jacobian(this, input_graph, inputs, which, &
       & rows, columns, weights, jacobian_defined)
    class(residual_operator), intent(in)  :: this
    class(directed_graph)    , intent(in)  :: input_graph
    type(binding)             , intent(in)  :: inputs(:)
    integer              , intent(in)  :: which
    integer , allocatable, intent(out) :: rows(:), columns(:)
    real(dp), allocatable, intent(out) :: weights(:)
    logical              , intent(out) :: jacobian_defined
    type(stored_field) :: state
    real(dp), allocatable :: x(:), w(:), column(:), xs(:), nu(:), g(:)
    integer , allocatable :: r(:), c(:), reads(:)
    logical , allocatable :: is_fixed(:)
    integer :: e, d, p, npts, n, num_triples, count, j, k, deg
    real(dp) :: value
    jacobian_defined = which == 1
    if (.not. jacobian_defined) return
    call require_host(this, input_graph)
    n    = this % unknowns
    npts = size(this % at)
    deg  = this % degrees
    is_fixed = this % fixed_indicator()
    call state_of(this, inputs, x, state)
    if (allocated(this % physics)) then
       if (.not. bound_on(inputs, this % argument(2), this % design_domain(), npts * this % designs)) then
          error stop 'operation_residual: the design must be defined on the point domain with one value per &
               &design per point'
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
          call this % rules(j) % gradient_at(xs((p - 1) * deg + 1:p * deg), &
               & nu((p - 1) * this % designs + 1:p * this % designs), value, g)
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
  end subroutine residual_explicit_jacobian

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
    logical :: jacobian_defined
    character(len=250) :: message

    if (size(rhs) /= this % unknowns) then
       write(message,'(a,i0,a,i0)') 'operation_residual: one right side value is required per &
            &unknown; size(rhs) = ', size(rhs), ', unknowns = ', this % unknowns
       error stop trim(message)
    end if

    call this % explicit_jacobian(input_graph, inputs, 1, r, c, w, jacobian_defined)
    if (.not. jacobian_defined) then
       error stop 'operation_residual: linearize requires the jacobian in the state to be &
            &explicit, but explicit_jacobian() reported none'
    end if

    a = stencil(r, c, w, spread(0.0_dp, 1, this % unknowns), 'explicit jacobian')
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
    if (allocated(a % objective) .and. allocated(b % objective)) then
       error stop 'operation_residual: one objective enters the Lagrangian'
    end if
    if (allocated(a % objective)) this % objective = a % objective
    if (allocated(b % objective)) this % objective = b % objective

  end function residual_sum

  !===================================================================!
  ! The objective added to a residual, J + lambda . r: a scalar
  ! field, the integral of a function of the state over the
  ! manifold. Invalid input: a field that is not an integral or not
  ! scalar, or a residual with an objective already.
  !===================================================================!

  function objective_sum(objective, residual) result(this)

    type(continuous_field)   , intent(in) :: objective
    type(continuous_residual), intent(in) :: residual
    type(continuous_residual) :: this

    if (.not. objective % is_integral()) then
       error stop 'operation_residual: the objective is the integral of a field over the manifold'
    end if
    if (objective % num_components() /= 1) then
       error stop 'operation_residual: the objective is a scalar field'
    end if
    if (allocated(residual % objective)) then
       error stop 'operation_residual: one objective enters the Lagrangian'
    end if
    this = residual
    this % objective = objective

  end function objective_sum

  !===================================================================!
  ! THE MULTIPLIER of a term: an unknown of the term's manifold with
  ! one component per equation, the count read from the term - the
  ! initial data of an equation of order n are n conditions, and
  ! their multiplier has n components - a function of the coordinates
  ! given as its arguments. Invalid input: a residual of several
  ! terms.
  !===================================================================!

  function residual_multiplier(this, name, arguments) result(lambda)

    class(continuous_residual), intent(inout) :: this
    character(len=*)          , intent(in)    :: name
    type(continuous_field)    , intent(in), optional :: arguments(:)
    type(continuous_field) :: lambda

    if (size(this % term) /= 1) then
       error stop 'operation_residual: a multiplier is declared by one term; declare it before summing'
    end if
    lambda = this % term(1) % manifold % unknown(name, arguments, components=size(this % term(1) % equation))

  end function residual_multiplier

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
    class(derivative_approximation), intent(in), optional :: space
    type(discrete_residual) :: image

    integer :: k, j, i, whole_term, num_conditions, nf
    logical :: covering
    real(dp) :: given
    real(dp), allocatable :: slope(:)
    character(len=250) :: message

    whole_term = 0
    num_conditions = 0
    do k = 1, size(this % term)
       associate (t => this % term(k))
         if (t % manifold % part == WHOLE) then
            if (whole_term /= 0) then
               error stop 'operation_residual: one term states the equations on the whole manifold'
            end if
            whole_term = k
         else
            if (.not. t % paired) then
               error stop 'operation_residual: a condition on a part of the manifold is paired with a multiplier'
            end if
            num_conditions = num_conditions + 1
         end if
       end associate
    end do
    if (whole_term == 0) then
       error stop 'operation_residual: no term states the equations on the whole manifold'
    end if

    associate (t => this % term(whole_term))
      if (.not. points % identity % matches(t % manifold % identity)) then
         error stop 'operation_residual: the points discretize another manifold than the equations'' own'
      end if
      ! the unknowns of the state: every unknown of the manifold but
      ! the multiplier of the equations, declared after them
      nf = t % manifold % num_unknowns
      if (t % paired) nf = nf - t % multiplier % num_components()
      if (size(t % equation) /= nf) then
         write(message,'(a,i0,a,i0)') 'operation_residual: one equation is required per unknown of the &
              &manifold, its multiplier aside; equations = ', size(t % equation), ', unknowns = ', nf
         error stop trim(message)
      end if
      image % manifold = t % manifold
      image % manifold % num_unknowns = nf
      if (t % paired) then
         image % with_adjoint = .true.
         image % adjoint_name = t % multiplier % name
      end if
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
       if (k == whole_term) cycle
       j = j + 1
       image % condition(j) = this % term(k)
       if (.not. image % condition(j) % manifold % parent % matches(image % manifold % identity)) then
          error stop 'operation_residual: a condition is stated on a part of another manifold'
       end if
       ! a condition on the design factor fixes the design coordinate
       ! at the value the manifold states: nu - nu_0 = 0, one equation
       if (image % condition(j) % manifold % is_design_factor()) then
          if (size(image % condition(j) % equation) /= points % num_designs()) then
             error stop 'operation_residual: the condition on the design factor fixes every design coordinate, &
                  &one equation per coordinate in the order of the coordinates'
          end if
          do i = 1, points % num_designs()
             call design_condition(image % condition(j) % equation(i), points % design_values(), given, slope)
             if (abs(given) > 1.0e-12_dp * max(1.0_dp, maxval(abs(points % design_values())))) then
                error stop 'operation_residual: the condition on the design factor is each design coordinate &
                     &minus the value the manifold fixes it at, nu - nu_0'
             end if
             slope(i) = slope(i) - 1.0_dp
             if (any(abs(slope) > 1.0e-12_dp)) then
                error stop 'operation_residual: the condition on the design factor is each design coordinate &
                     &minus its value, equation k the k-th coordinate'
             end if
          end do
       end if
    end do

    ! the objective covers the points: the integral over the whole
    ! manifold, over the one factor of a manifold of time or of space
    ! alone, or over the face at an instant
    if (allocated(this % objective)) then
       covering = .false.
       select case (this % objective % integrated_part)
       case (WHOLE)
          covering = this % objective % integrated_over % matches(image % manifold % identity)
       case (TIME_FACE)
          covering = this % objective % integrated_parent % matches(image % manifold % identity)
       case (TIME_FACTOR)
          covering = this % objective % integrated_parent % matches(image % manifold % identity) &
               & .and. .not. points % with_space
       case (SPACE_FACTOR)
          covering = this % objective % integrated_parent % matches(image % manifold % identity) &
               & .and. .not. points % with_time
       end select
       if (.not. covering) then
          error stop 'operation_residual: the objective is integrated over the whole manifold, over the &
               &one factor of a manifold of time or of space alone, or over the face at an instant; a &
               &factor of a product does not cover the points'
       end if
       if (.not. image % with_adjoint) then
          error stop 'operation_residual: an objective requires the equations paired with a multiplier, &
               &the adjoint, whose equation is the stationarity of the Lagrangian in the state'
       end if
       image % objective = this % objective
    end if

    image % points = points
    if (points % with_time .neqv. present(time)) then
       error stop 'operation_residual: the derivatives along time are approximated by a chain of &
            &families, given exactly when the manifold has a time coordinate'
    end if
    if (points % with_space .neqv. present(space)) then
       error stop 'operation_residual: the derivatives along space are approximated by a method of an &
            &order, finite differences or finite volumes, given exactly when the manifold has a region'
    end if
    if (present(time)) then
       image % schemes    = time
       image % with_chain = .true.
       if (time % from(size(time % from)) > points % num_instants()) then
          error stop 'operation_residual: a block of the chain begins after the last instant'
       end if
    end if
    if (present(space)) then
       allocate(image % differences, source=space)
       image % with_differences = .true.
    end if

  end function residual_discretize

  !===================================================================!
  ! THE ZERO OF THE DISCRETE RESIDUAL from an estimate: one block of
  ! the chain after another, each a residual_operator over the block's
  ! moments (its instants, and the stages of a one-step family) and
  ! the cells. The solution is a discrete field of the manifold's
  ! unknowns at the instants.
  !
  ! CHECKPOINTS: with checkpoint = k blocks, the stages of the steps
  ! and their jets along the design, which the reverse sweep reads at
  ! every block, are retained over one segment of k blocks at a time
  ! rather than over the chain: the forward sweep discards each
  ! segment's stages once the objective is summed over its quadrature
  ! nodes, and the reverse sweep forms each segment again from the
  ! state at its instants, which is retained, before the multipliers
  ! of its blocks. The storage of the stages falls from the chain to
  ! a segment at the cost of a second forward sweep, from the same
  ! estimate and through the same iterates. The solution then
  ! integrates the objective alone, whose value it stores; a chain
  ! pipelined over the images retains every stage.
  !===================================================================!

  subroutine minimize(this, estimate, solution, checkpoint)

    class(discrete_residual), intent(in)  :: this
    type(discrete_field)    , intent(in)  :: estimate
    type(discrete_field)    , intent(out) :: solution
    integer                 , intent(in), optional :: checkpoint

    real(dp), allocatable :: tuple(:,:,:), stages(:,:,:,:), tjet(:,:,:,:,:), sjet(:,:,:,:,:,:), value(:,:)
    real(dp), allocatable :: jet(:,:,:,:), adjoint(:,:,:,:,:), reaction(:,:,:,:), gauge(:,:,:,:)
    real(dp), allocatable :: design_multiplier(:,:,:)
    integer , allocatable :: slot(:), stage_node(:,:,:)
    type(residual_operator) :: rows
    real(dp), allocatable :: q(:)
    integer :: order, e, g, nm, col, nd, d
    character(len=32), allocatable :: names(:)
    type(stencil), allocatable :: rows_along(:,:)
    integer :: stride, nf, ncells, ninst, f, c, k, b, first, last, depth, nb, top, o, npts, nc, j
    integer :: segment, nseg, sg, b1, b2, lo, hi
    real(dp), allocatable :: objective_value(:), objective_partial(:,:,:), part_value(:), part_partial(:,:,:)
    real(dp), allocatable :: part_gradient(:,:,:,:)
    type(discrete_field) :: part
    type(expression_view) :: view
    integer(8) :: tick, tock, rate
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
    ! estimate, every derivative component zero; and their
    ! derivatives along the design to the order of the expansion
    order = this % points % design_order()
    nd    = max(1, this % points % num_expansions())
    segment = 0
    if (present(checkpoint)) segment = checkpoint
    if (segment < 0) error stop 'operation_residual: checkpoint is a number of blocks, at least one'
    if (segment > 0 .and. .not. this % with_chain) then
       error stop 'operation_residual: checkpoints divide a chain of blocks; the residual has none'
    end if
    if (segment > 0 .and. pipelined(this)) then
       error stop 'operation_residual: a chain pipelined over the images retains every stage; no checkpoints'
    end if
    allocate(tuple(stride, ncells, ninst), source=0.0_dp)
    allocate(tjet(stride, ncells, ninst, order, nd), source=0.0_dp)
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

    call system_clock(tick, rate)
    if (this % with_chain) then
       nb = this % schemes % num_blocks()
       ! the stages of the step into each instant, one set per family
       ! that marches by stages, retained for the adjoint
       top = 1
       do b = 1, nb
          top = max(top, this % schemes % scheme(b) % num_stages())
       end do
       do b = 1, nb
          call block_extent(this, b, first, last, depth)
          if (first < 1) then
             write(message,'(a,i0,a,i0,a,i0)') 'operation_residual: block ', b, ' reads ', depth, &
                  & ' instants of history before instant ', first + depth
             error stop trim(message)
          end if
       end do
       if (segment == 0) then
          allocate(stages(stride, ncells, top, ninst), sjet(stride, ncells, top, ninst, order, nd), source=0.0_dp)
          ! the wavefront over the blocks and the orders, when the images
          ! number the orders plus one on a manifold of time alone with
          ! one expansion; the block loop otherwise
          if (pipelined(this)) then
             call wavefront_expansion(this, rows_along, tuple, stages, tjet, sjet)
             nb = 0
          end if
          do b = 1, nb
             call block_extent(this, b, first, last, depth)
             call solve_block(this, this % schemes % scheme(b), first, last, depth, rows_along, tuple, &
                  & 1, stages, tjet, sjet)
          end do
       else
          ! segment by segment: the stages of a segment, the objective
          ! summed over its nodes, the stages released
          nseg = (nb + segment - 1) / segment
          allocate(objective_value(1 + nd * order), source=0.0_dp)
          allocate(objective_partial(size(this % points % design_values()), 0:order, nd), source=0.0_dp)
          do sg = 1, nseg
             call segment_extent(this, sg, segment, b1, b2, lo, hi)
             if (verbosity >= 1) print '(a,i0,a,i0,a,i0)', 'segment ', sg, '  blocks ', b1, '..', b2
             allocate(stages(stride, ncells, top, lo:hi), sjet(stride, ncells, top, lo:hi, order, nd), source=0.0_dp)
             call segment_solved(this, b1, b2, rows_along, tuple, lo, stages, tjet, sjet)
             if (allocated(this % objective)) then
                call segment_nodes(this, b1, b2, tuple, lo, stages, tjet, sjet, part, stage_node)
                call this % objective % functional(part, part_value, part_gradient, part_partial)
                objective_value   = objective_value + part_value
                objective_partial = objective_partial + part_partial
             end if
             deallocate(stages, sjet)
          end do
          allocate(stages(stride, ncells, top, 1:0), sjet(stride, ncells, top, 1:0, order, nd))
       end if
    else
       ! the point manifold: one system, its jets along the design
       allocate(stages(stride, 1, 1, 1), sjet(stride, 1, 1, 1, order, nd), source=0.0_dp)
       call solve_point(this, tuple, rows, q)
       if (order > 0) call design_jets(this, rows, q, [0], [1], 1, 1, 0, tjet, 1, sjet)
    end if
    ! the expansions shared among the images: their jets summed once;
    ! segment by segment above when the stages are not retained
    if (order > 0 .and. num_images() > 1 .and. .not. this % points % with_space .and. segment == 0) then
       call co_sum(tjet)
       call co_sum(sjet)
    end if

    ! THE SOLUTION: the state at every point, with the jet the rule
    ! lays out; with the equations paired, the multipliers of every
    ! equality of L_h by the reverse sweep - each equation's under
    ! the name its term declared, each face condition's under its
    ! multiplier's name at the face's points
    npts = this % points % num_points()
    nc = 0
    if (this % with_adjoint) then
       do j = 1, size(this % condition)
          if (this % condition(j) % manifold % is_face())          nc = nc + size(this % condition(j) % equation)
          if (this % condition(j) % manifold % is_time_factor())   nc = nc + size(this % condition(j) % equation)
          if (this % condition(j) % manifold % is_design_factor()) nc = nc + size(this % condition(j) % equation)
       end do
    end if
    ! the multiplier columns: the equations', then the conditions'
    nm = merge(nf + nc, 0, this % with_adjoint)
    allocate(value(npts, nf + nm), source=0.0_dp)
    allocate(names(size(value, 2)), jet(stride + nm, npts, order + 1, nd))
    allocate(slot(size(value, 2)), source=0)
    jet = 0.0_dp
    names(1:nf) = this % manifold % unknown_name(1:nf)
    do f = 1, nf
       slot(f) = this % rule % offset_of_field(f) + 1
    end do
    do k = 1, ninst
       do c = 1, ncells
          do f = 1, nf
             value(this % points % point_of(k, c), f) = tuple(this % rule % offset_of_field(f) + 1, c, k)
          end do
          do d = 1, nd
             jet(1:stride, this % points % point_of(k, c), 1, d) = tuple(:, c, k)
             do o = 1, order
                jet(1:stride, this % points % point_of(k, c), o + 1, d) = tjet(:, c, k, o, d)
             end do
          end do
       end do
    end do

    solution = discrete_field(this % points, value(:, 1:nf), names(1:nf))
    solution % jet  = jet
    solution % slot = slot(1:nf)
    solution % rule = this % rule
    if (segment == 0) then
       call quadrature_nodes(this, tuple, 1, stages, tjet, sjet, solution, stage_node)
    else
       ! the objective's value stored, the nodes not
       solution % objective_value   = objective_value
       solution % objective_formula = ''
       if (allocated(this % objective)) then
          view = expression_view(this % objective % component(1))
          solution % objective_formula = view % formula()
       end if
    end if
    call system_clock(tock)
    if (verbosity >= 1) print '(a,f10.3,a,i0,a)', 'the forward sweep with its expansions took ', &
         & real(tock - tick, dp) / real(rate, dp), ' s; resident memory ', resident_kilobytes() / 1024, ' MB'

    if (this % with_adjoint) then
       tick = tock
       call adjoint_sweep(this, rows_along, tuple, 1, stages, tjet, sjet, solution, stage_node, segment, &
            & objective_partial, estimate, adjoint, reaction, gauge, design_multiplier)
       call system_clock(tock)
       if (verbosity >= 1) print '(a,f10.3,a,i0,a)', 'the reverse sweep with its expansions took ', &
            & real(tock - tick, dp) / real(rate, dp), ' s; resident memory ', resident_kilobytes() / 1024, ' MB'
       ! every multiplier and its derivatives along the design as a
       ! column with a row of the jet of its own
       names(nf + 1:2 * nf) = this % adjoint_name
       do k = 1, ninst
          do c = 1, ncells
             do o = 0, order
                jet(stride + 1:stride + nf, this % points % point_of(k, c), o + 1, :) = adjoint(:, c, k, o, :)
             end do
          end do
       end do
       nc = 0
       e  = 0
       g  = 0
       do j = 1, size(this % condition)
          if (this % condition(j) % manifold % is_design_factor()) then
             ! the multipliers of the design conditions, one per design
             ! coordinate, functions on the design factor: one value
             ! each, at every point
             do f = 1, size(this % condition(j) % equation)
                nc = nc + 1
                names(2 * nf + nc) = this % condition(j) % multiplier % name
                do o = 0, order
                   do d = 1, nd
                      jet(stride + nf + nc, :, o + 1, d) = design_multiplier(f, o, d)
                   end do
                end do
             end do
          end if
          if (this % condition(j) % manifold % is_time_factor()) then
             ! the multiplier of a time-factor condition, a function
             ! of time alone: its value at every cell of the instant
             do f = 1, size(this % condition(j) % equation)
                nc = nc + 1
                g  = g + 1
                names(2 * nf + nc) = this % condition(j) % multiplier % name
                do k = 1, ninst
                   do c = 1, ncells
                      do o = 0, order
                         jet(stride + nf + nc, this % points % point_of(k, c), o + 1, :) = gauge(k, g, o, :)
                      end do
                   end do
                end do
             end do
          end if
          if (.not. this % condition(j) % manifold % is_face()) cycle
          k = face_instant(this, this % condition(j) % manifold % face_time)
          do f = 1, size(this % condition(j) % equation)
             nc = nc + 1
             e  = e + 1
             names(2 * nf + nc) = this % condition(j) % multiplier % name
             do c = 1, ncells
                do o = 0, order
                   jet(stride + nf + nc, this % points % point_of(k, c), o + 1, :) = reaction(c, e, o, :)
                end do
             end do
          end do
       end do
       do col = 1, nm
          slot(nf + col)     = stride + col
          value(:, nf + col) = jet(stride + col, :, 1, 1)
       end do
    end if

    estimate_nodes: block
      type(discrete_field) :: with_nodes
      with_nodes = solution
      solution = discrete_field(this % points, value, names)
      solution % jet  = jet
      solution % slot = slot
      solution % rule = this % rule
      if (segment == 0) then
         solution % node_jet      = with_nodes % node_jet
         solution % node_weight   = with_nodes % node_weight
         solution % node_position = with_nodes % node_position
         solution % node_point    = with_nodes % node_point
         solution % node_image    = with_nodes % node_image
      else
         solution % objective_value   = with_nodes % objective_value
         solution % objective_formula = with_nodes % objective_formula
      end if
    end block estimate_nodes

  end subroutine minimize

  !===================================================================!
  ! THE QUADRATURE NODES of a solution: the families' own rule for the
  ! integral of a functional over the manifold. Every point is a node
  ! - an instant at a cell - with the weight the multistep families'
  ! interpolatory rules over their steps give it, accumulated from
  ! every step whose rule reads the instant, the step times the
  ! weight over the unit step, and the cell's volume; every stage of a
  ! staged family's steps is a node with the tableau's weight times
  ! the step and the cell's volume, at the time
  ! the instant behind plus the abscissa times the step, its jet from
  ! the stages solved and their jets along the design; a point
  ! manifold's one point has weight one. stage_node(k, i, c) is the
  ! node of stage i of the step into instant k at cell c, zero where
  ! there is none.
  !===================================================================!

  subroutine quadrature_nodes(this, tuple, sbase, stages, tjet, sjet, solution, stage_node, blocks)

    class(discrete_residual), intent(in)    :: this
    real(dp)                , intent(in)    :: tuple(:,:,:)
    integer                 , intent(in)    :: sbase
    real(dp)                , intent(in)    :: stages(:,:,:,sbase:), tjet(:,:,:,:,:), sjet(:,:,:,sbase:,:,:)
    type(discrete_field)    , intent(inout) :: solution
    integer , allocatable   , intent(out)   :: stage_node(:,:,:)
    integer                 , intent(in), optional :: blocks(2)

    type(derivative_terms), allocatable :: dt(:), weight(:)
    integer , allocatable :: instant(:), cell_owner(:)
    real(dp), allocatable :: position(:)
    integer :: stride, ncells, ninst, npts, order, nd, nb, b, first, last, depth, n, s, smax, k, i, c, j
    integer :: nnodes, nn, nf, d, o, b1, b2, klo, khi
    logical :: staged
    real(dp) :: h

    stride = this % rule % num_components()
    nf     = this % manifold % num_unknowns
    ncells = this % points % num_cells()
    ninst  = this % points % num_instants()
    npts   = this % points % num_points()
    order  = this % points % design_order()
    nd     = max(1, this % points % num_expansions())

    ! the blocks whose nodes are formed: every block, or the range
    ! given; the points of every instant are nodes, those outside
    ! the range with no measure and no point
    smax   = 1
    nnodes = npts
    nb     = 0
    b1     = 1
    b2     = 0
    klo    = 1
    khi    = ninst
    if (this % with_chain) then
       nb = this % schemes % num_blocks()
       b2 = nb
       if (present(blocks)) then
          b1 = blocks(1)
          b2 = blocks(2)
          call block_extent(this, b1, first, last, depth)
          klo = first + depth
          call block_extent(this, b2, first, last, depth)
          khi = last
       end if
       do b = 1, nb
          smax = max(smax, this % schemes % scheme(b) % num_stages())
          if (b < b1 .or. b > b2) cycle
          call block_extent(this, b, first, last, depth)
          if (marches_by_stages(this % schemes % scheme(b), top_degree(this % rule, nf) + 1)) then
             nnodes = nnodes + (last - first) * this % schemes % scheme(b) % num_stages() * ncells
          end if
       end do
    end if
    allocate(stage_node(ninst, smax, ncells), source=0)
    allocate(solution % node_jet(stride, nnodes, order + 1, nd), source=0.0_dp)
    allocate(solution % node_weight(nnodes), source=0.0_dp)
    allocate(solution % node_position(this % points % num_coordinates(), nnodes), source=0.0_dp)
    allocate(solution % node_point(nnodes), source=0)
    allocate(solution % node_image(nnodes), source=this_image())
    cell_owner = cell_owners(this)

    ! the points as nodes, their jets those of the solution
    do k = klo, khi
       do c = 1, ncells
          nn = this % points % point_of(k, c)
          solution % node_point(nn)       = nn
          solution % node_image(nn)       = cell_owner(c)
          solution % node_position(:, nn) = this % points % position(nn)
          do d = 1, nd
             solution % node_jet(:, nn, 1, d) = tuple(:, c, k)
             do o = 1, order
                solution % node_jet(:, nn, o + 1, d) = tjet(:, c, k, o, d)
             end do
          end do
       end do
    end do

    if (.not. this % with_chain) then
       solution % node_weight = 1.0_dp
       return
    end if

    nn = npts
    do b = b1, b2
       call block_extent(this, b, first, last, depth)
       n      = last - first + 1
       staged = marches_by_stages(this % schemes % scheme(b), top_degree(this % rule, nf) + 1)
       associate (scheme => this % schemes % scheme(b))
         if (staged) then
            s = scheme % num_stages()
            do k = 1, n - 1
               h = this % points % step(first + k)
               do i = 1, s
                  do c = 1, ncells
                     nn = nn + 1
                     stage_node(first + k, i, c) = nn
                     solution % node_image(nn)   = cell_owner(c)
                     solution % node_weight(nn)  = h * scheme % stage_weight(i) * cell_volume(c)
                     position = this % points % position(this % points % point_of(first + k - 1, c))
                     position(1) = position(1) + scheme % stage_abscissa(i) * h
                     solution % node_position(:, nn) = position
                     do d = 1, nd
                        solution % node_jet(:, nn, 1, d) = stages(:, c, i, first + k)
                        do o = 1, order
                           solution % node_jet(:, nn, o + 1, d) = sjet(:, c, i, first + k, o, d)
                        end do
                     end do
                  end do
               end do
            end do
         else
            allocate(dt(n))
            do j = 1, n
               dt(j) = derivative_terms(this % points % step(first + j - 1), 0)
            end do
            do k = depth + 1, n
               call scheme % step_quadrature(dt, k, weight, instant=instant)
               h = this % points % step(first + k - 1)
               do j = 1, size(weight)
                  do c = 1, ncells
                     solution % node_weight(this % points % point_of(first + instant(j) - 1, c)) = &
                          & solution % node_weight(this % points % point_of(first + instant(j) - 1, c)) &
                          & + h * value(weight(j)) * cell_volume(c)
                  end do
               end do
            end do
            deallocate(dt)
         end if
       end associate
    end do

  contains

    ! the cell's volume, one without a region
    pure real(dp) function cell_volume(c)
      integer, intent(in) :: c
      cell_volume = 1.0_dp
      if (this % points % with_space) cell_volume = this % points % volume(c)
    end function cell_volume

  end subroutine quadrature_nodes

  ! the instants of block b: first - depth .. last, depth of history
  subroutine block_extent(this, b, first, last, depth)
    class(discrete_residual), intent(in)  :: this
    integer                 , intent(in)  :: b
    integer                 , intent(out) :: first, last, depth
    integer :: nb
    nb    = this % schemes % num_blocks()
    first = this % schemes % from(b)
    last  = this % points % num_instants()
    if (b < nb) last = this % schemes % from(b + 1) - 1
    depth = 0
    if (b > 1) depth = this % schemes % scheme(b) % history_depth(top_degree(this % rule, this % manifold % num_unknowns))
    first = first - depth
  end subroutine block_extent

  !===================================================================!
  ! The resident memory of the process in kilobytes, from the kernel's
  ! status of the process; zero where that is not readable.
  !===================================================================!

  integer function resident_kilobytes()
    character(len=128) :: line
    integer :: unit, status
    resident_kilobytes = 0
    open(newunit=unit, file='/proc/self/status', status='old', action='read', iostat=status)
    if (status /= 0) return
    do
       read(unit, '(a)', iostat=status) line
       if (status /= 0) exit
       if (line(1:6) == 'VmRSS:') then
          read(line(7:), *, iostat=status) resident_kilobytes
          exit
       end if
    end do
    close(unit)
  end function resident_kilobytes

  !===================================================================!
  ! The segment sg of k blocks: its blocks b1..b2, and the instants
  ! lo..hi of the stages its blocks determine, the step into each own
  ! instant of each block.
  !===================================================================!

  subroutine segment_extent(this, sg, k, b1, b2, lo, hi)
    class(discrete_residual), intent(in)  :: this
    integer                 , intent(in)  :: sg, k
    integer                 , intent(out) :: b1, b2, lo, hi
    integer :: first, last, depth
    b1 = (sg - 1) * k + 1
    b2 = min(this % schemes % num_blocks(), sg * k)
    call block_extent(this, b1, first, last, depth)
    lo = first + 1
    call block_extent(this, b2, first, last, depth)
    hi = last
  end subroutine segment_extent

  !===================================================================!
  ! The forward solve of the blocks b1..b2 into the stage arrays that
  ! begin at the instant lo; the jets of expansions the images share
  ! summed over the segment's instants once.
  !===================================================================!

  subroutine segment_solved(this, b1, b2, rows_along, tuple, lo, stages, tjet, sjet)

    class(discrete_residual), intent(in)    :: this
    integer                 , intent(in)    :: b1, b2, lo
    type(stencil)           , intent(in)    :: rows_along(:,:)
    real(dp)                , intent(inout) :: tuple(:,:,:), stages(:,:,:,lo:), tjet(:,:,:,:,:), sjet(:,:,:,lo:,:,:)

    real(dp), allocatable :: own(:,:,:,:,:)
    integer :: b, first, last, depth, hi

    do b = b1, b2
       call block_extent(this, b, first, last, depth)
       call solve_block(this, this % schemes % scheme(b), first, last, depth, rows_along, tuple, lo, stages, &
            & tjet, sjet)
    end do
    if (this % points % design_order() > 0 .and. num_images() > 1 .and. .not. this % points % with_space) then
       hi  = ubound(stages, 4)
       own = tjet(:, :, lo:hi, :, :)
       call co_sum(own)
       tjet(:, :, lo:hi, :, :) = own
       call co_sum(sjet)
    end if

  end subroutine segment_solved

  !===================================================================!
  ! The quadrature nodes of the blocks b1..b2 as a field: the
  ! manifold, the rule and the nodes, the jet at the points left out.
  !===================================================================!

  subroutine segment_nodes(this, b1, b2, tuple, lo, stages, tjet, sjet, part, stage_node)

    class(discrete_residual), intent(in)  :: this
    integer                 , intent(in)  :: b1, b2, lo
    real(dp)                , intent(in)  :: tuple(:,:,:), stages(:,:,:,lo:), tjet(:,:,:,:,:), sjet(:,:,:,lo:,:,:)
    type(discrete_field)    , intent(out) :: part
    integer , allocatable   , intent(out) :: stage_node(:,:,:)

    allocate(part % on, source=this % points)
    part % rule = this % rule
    call quadrature_nodes(this, tuple, lo, stages, tjet, sjet, part, stage_node, [b1, b2])

  end subroutine segment_nodes

  !===================================================================!
  ! THE REVERSE SWEEP: the multipliers of every row of L_h and their
  ! derivatives along each design coordinate to the order of the
  ! expansion, block by block from the last, each block's transposed
  ! jacobian at its solution once at order zero and once per
  ! coordinate and order above. The differential of the objective in
  ! the jet is the source at order zero; above it, the derivatives
  ! along the coordinate of that differential and of the product of
  ! the multipliers with the equations' jacobian, the multiplier's own
  ! derivative of that order left out. adjoint(f, c, k, o, d) is the
  ! o-th derivative along the coordinate d of the multiplier of
  ! equation f at cell c and instant k, reaction(c, e, o, d) of the
  ! e-th face condition equation at cell c, gauge(k, g, o, d) of the
  ! g-th time-factor condition equation at instant k, and
  ! design_multiplier(i, o, d) of the multiplier of the design
  ! condition on the coordinate i, kappa_i = -dJ/dnu_i, so that
  ! design_multiplier(i, o, d) = -d^(o+1) J / (dnu_i dnu_d^o). A
  ! residual without a chain has no blocks to sweep, and is refused.
  !===================================================================!

  subroutine adjoint_sweep(this, rows_along, tuple, sbase, stages, tjet, sjet, solution, stage_node, segment, &
       & objective_partial, estimate, adjoint, reaction, gauge, design_multiplier)

    class(discrete_residual), intent(in)  :: this
    type(stencil)           , intent(in)  :: rows_along(:,:)
    real(dp)                , intent(in)  :: tuple(:,:,:)
    integer                 , intent(in)  :: sbase
    real(dp)                , intent(in)  :: stages(:,:,:,sbase:), tjet(:,:,:,:,:), sjet(:,:,:,sbase:,:,:)
    type(discrete_field)    , intent(in)  :: solution
    integer , allocatable   , intent(in)  :: stage_node(:,:,:)
    integer                 , intent(in)  :: segment
    real(dp), allocatable   , intent(in)  :: objective_partial(:,:,:)
    type(discrete_field)    , intent(in)  :: estimate
    real(dp), allocatable   , intent(out) :: adjoint(:,:,:,:,:), reaction(:,:,:,:), gauge(:,:,:,:)
    real(dp), allocatable   , intent(out) :: design_multiplier(:,:,:)

    real(dp), allocatable :: gradient(:,:,:,:), partial(:,:,:), incoming(:,:,:,:,:), q(:), total(:)
    real(dp), allocatable :: tuple_again(:,:,:), tjet_again(:,:,:,:,:), stages_seg(:,:,:,:), sjet_seg(:,:,:,:,:,:)
    integer , allocatable :: stage_node_seg(:,:,:)
    type(discrete_field) :: part
    type(residual_operator) :: rows
    integer , allocatable :: primary(:)
    integer :: stride, nf, ncells, ninst, nc, ng, j, b, nb, first, last, depth, order, nd, ni, d
    integer :: nseg, sg, b1, b2, lo, hi, top, k, c, f

    stride = this % rule % num_components()
    nf     = this % manifold % num_unknowns
    ncells = this % points % num_cells()
    ninst  = this % points % num_instants()
    order  = this % points % design_order()
    nd     = max(1, this % points % num_expansions())
    ni     = size(this % points % design_values())

    if (segment > 0) then
       partial = objective_partial
    else if (allocated(this % objective)) then
       call this % objective % differential(solution, gradient, partial)
    else
       allocate(gradient(stride, size(solution % node_weight), 0:order, nd), partial(ni, 0:order, nd), source=0.0_dp)
    end if

    nc = 0
    ng = 0
    do j = 1, size(this % condition)
       if (this % condition(j) % manifold % is_face())        nc = nc + size(this % condition(j) % equation)
       if (this % condition(j) % manifold % is_time_factor()) ng = ng + size(this % condition(j) % equation)
    end do
    allocate(adjoint(nf, ncells, ninst, 0:order, nd), incoming(stride, ncells, ninst, 0:order, nd), source=0.0_dp)
    allocate(reaction(ncells, nc, 0:order, nd), gauge(ninst, ng, 0:order, nd), source=0.0_dp)
    allocate(design_multiplier(ni, 0:order, nd), source=0.0_dp)

    ! THE MULTIPLIERS OF THE DESIGN CONDITIONS nu_i - nu_i0 = 0 and
    ! their derivatives: the stationarity of L in nu_i, kappa_i =
    ! -(dJ/dnu_i at fixed state + the sum over every row of its
    ! multiplier times the row's partial in nu_i) = -dJ/dnu_i, the
    ! total derivative by the adjoint; the partial in nu_i here, the
    ! rows' sums by the blocks
    design_multiplier = -partial
    ! an expansion another image owns: its derivatives above order
    ! zero are left zero here, and summed over the images at the end
    do d = 1, nd
       if (.not. expansion_owned(this, d)) design_multiplier(:, 1:, d) = 0.0_dp
    end do

    if (.not. this % with_chain) then
       ! the point manifold: one system, its transposed solve
       call point_rows(this, tuple, rows, q, primary)
       call adjoint_rows(this, 1, 1, 0, rows, q, [0], [1], primary, tjet, sbase, sjet, gradient, stage_node, &
            & adjoint, incoming, reaction, gauge, design_multiplier)
       return
    end if

    if (pipelined(this)) then
       ! the wavefront in reverse: image o + 1 forms order o of every
       ! block, one block behind image o; each image's multipliers of
       ! its own order, summed over the images at the end
       design_multiplier = 0.0_dp
       design_multiplier(:, this_image() - 1, 1) = -partial(:, this_image() - 1, 1)
       call wavefront_adjoint(this, rows_along, tuple, stages, tjet, sjet, gradient, stage_node, &
            & adjoint, incoming, reaction, gauge, design_multiplier)
       call co_sum(adjoint)
       call co_sum(reaction)
       call co_sum(gauge)
       call co_sum(design_multiplier)
       return
    end if

    nb = this % schemes % num_blocks()
    if (segment == 0) then
       do b = nb, 1, -1
          call block_extent(this, b, first, last, depth)
          call adjoint_block(this, this % schemes % scheme(b), first, last, depth, rows_along, tuple, &
               & sbase, stages, tjet, sjet, gradient, stage_node, adjoint, incoming, reaction, gauge, design_multiplier)
       end do
    else
       ! segment by segment from the last: the segment's blocks solved
       ! again from the estimate at their own instants, the same
       ! iterates as the forward sweep's, the objective's differential
       ! at the segment's nodes, then its blocks in reverse
       top = 1
       do b = 1, nb
          top = max(top, this % schemes % scheme(b) % num_stages())
       end do
       nseg = (nb + segment - 1) / segment
       tuple_again = tuple
       tjet_again  = tjet
       do sg = nseg, 1, -1
          call segment_extent(this, sg, segment, b1, b2, lo, hi)
          if (verbosity >= 1) print '(a,i0,a,i0,a,i0,a,i0,a)', 'segment ', sg, ' again  blocks ', b1, '..', b2, &
               & '  resident memory ', resident_kilobytes() / 1024, ' MB'
          allocate(stages_seg(stride, ncells, top, lo:hi), sjet_seg(stride, ncells, top, lo:hi, order, nd), source=0.0_dp)
          do b = b1, b2
             call block_extent(this, b, first, last, depth)
             do k = first + depth, last
                do c = 1, ncells
                   tuple_again(:, c, k) = 0.0_dp
                   do f = 1, nf
                      tuple_again(this % rule % offset_of_field(f) + 1, c, k) = &
                           & estimate % value(this % points % point_of(k, c), f)
                   end do
                end do
             end do
          end do
          call segment_solved(this, b1, b2, rows_along, tuple_again, lo, stages_seg, tjet_again, sjet_seg)
          call segment_nodes(this, b1, b2, tuple, lo, stages_seg, tjet, sjet_seg, part, stage_node_seg)
          if (allocated(this % objective)) then
             call this % objective % functional(part, total, gradient, partial)
          else
             allocate(gradient(stride, size(part % node_weight), 0:order, nd), source=0.0_dp)
          end if
          do b = b2, b1, -1
             call block_extent(this, b, first, last, depth)
             call adjoint_block(this, this % schemes % scheme(b), first, last, depth, rows_along, tuple, &
                  & lo, stages_seg, tjet, sjet_seg, gradient, stage_node_seg, adjoint, incoming, reaction, gauge, &
                  & design_multiplier)
          end do
          deallocate(stages_seg, sjet_seg, gradient)
       end do
    end if
    call shared_expansions_summed(this, adjoint, reaction, gauge, design_multiplier)

  end subroutine adjoint_sweep

  ! the derivatives along the expansions the images shared, summed
  ! once: order zero is complete on every image and is left alone
  subroutine shared_expansions_summed(this, adjoint, reaction, gauge, design_multiplier)

    class(discrete_residual), intent(in)    :: this
    real(dp)                , intent(inout) :: adjoint(:,:,:,0:,:), reaction(:,:,0:,:), gauge(:,:,0:,:)
    real(dp)                , intent(inout) :: design_multiplier(:,0:,:)

    real(dp), allocatable :: above(:,:,:,:,:), above4(:,:,:,:), above3(:,:,:)
    integer :: order

    order = ubound(adjoint, 4)
    if (order < 1 .or. num_images() == 1 .or. this % points % with_space) return
    above = adjoint(:, :, :, 1:order, :)
    call co_sum(above)
    adjoint(:, :, :, 1:order, :) = above
    above4 = reaction(:, :, 1:order, :)
    call co_sum(above4)
    reaction(:, :, 1:order, :) = above4
    above4 = gauge(:, :, 1:order, :)
    call co_sum(above4)
    gauge(:, :, 1:order, :) = above4
    above3 = design_multiplier(:, 1:order, :)
    call co_sum(above3)
    design_multiplier(:, 1:order, :) = above3

  end subroutine shared_expansions_summed

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
  ! ONE BLOCK'S ROWS: the instants first..last, the first `depth` of
  ! them known from the block before, the stages between consecutive
  ! instants when the family marches by stages. The unknown vector
  ! lists the moments in time order, each moment's cells in order,
  ! each cell's tuple of stride components: at((m - 1) ncells + c) is
  ! the offset of the tuple of moment m and cell c, instant_moment(k)
  ! the moment of the block's k-th instant, primary(f) the row of
  ! equation f within a tuple, and q the estimate - each instant's
  ! tuple and, at a stage of the step into instant k + 1, the stage
  ! stored at stages(:, c, i, k + 1) when the block was solved before,
  ! the tuple of instant k otherwise.
  !===================================================================!

  subroutine block_rows(this, scheme, first, last, depth, rows_along, tuple, sbase, stages, solved_stages, &
       & rows, q, at, instant_moment, primary)

    class(discrete_residual), intent(in)  :: this
    type(family)            , intent(in)  :: scheme
    integer                 , intent(in)  :: first, last, depth
    type(stencil)           , intent(in)  :: rows_along(:,:)
    real(dp)                , intent(in)  :: tuple(:,:,:)
    integer                 , intent(in)  :: sbase
    real(dp)                , intent(in)  :: stages(:,:,:,sbase:)
    logical                 , intent(in)  :: solved_stages
    type(residual_operator) , intent(out) :: rows
    real(dp), allocatable   , intent(out) :: q(:)
    integer , allocatable   , intent(out) :: at(:), instant_moment(:), primary(:)

    type(stencil) :: relations, gauge
    type(connectivity_graph) :: connectivity
    integer , allocatable :: stage_moment(:,:), fixed_rows(:)
    integer , allocatable :: determined(:), source(:), gauge_rows(:), gauge_columns(:), face_rows(:), face_columns(:)
    real(dp), allocatable :: weight(:), fixed(:), w(:), steps(:), gauge_weights(:), zeros(:), face_weights(:)
    real(dp), allocatable :: face_constants(:)
    integer :: stride, nf, ncells, n, s, nm, m, k, i, f, c, e, count, nd, width, unknowns, npts
    integer :: comp, axis, o, j, ne, history_last
    logical :: staged

    stride = this % rule % num_components()
    nf     = this % manifold % num_unknowns
    ncells = this % points % num_cells()
    n      = last - first + 1
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

    ! THE PRESCRIBED ROWS: the gauge rows - at every moment after the
    ! history, the mean of the integrated unknown over the cells, in
    ! place of its own row at the first cell - and the face rows at
    ! the block's own instants, one stencil with the faces' constants
    history_last = 0
    if (depth > 0) history_last = instant_moment(depth)
    call gauge_rows_of(this, nm, history_last, at, primary, gauge_rows, gauge_columns, gauge_weights)
    call face_rows_of(this, first, last, depth, instant_moment, at, face_rows, face_columns, face_weights, face_constants)
    allocate(zeros(unknowns), source=0.0_dp)
    do e = 1, size(face_rows)
       zeros(face_rows(e)) = face_constants(e)
    end do
    if (size(gauge_rows) + size(face_rows) > 0) then
       gauge = stencil([gauge_rows, face_rows], [gauge_columns, face_columns], [gauge_weights, face_weights], &
            & zeros, 'prescribed')
       rows  = residual_operator(relations, this % rule, at, unknowns, stride, primary, fixed_rows, fixed, &
            & prescribed=gauge)
    else
       rows  = residual_operator(relations, this % rule, at, unknowns, stride, primary, fixed_rows, fixed)
    end if

    ! the estimate: the tuple at each instant; at a stage of the step
    ! into instant k + 1, the stage solved before when solved_stages,
    ! the tuple of instant k otherwise
    allocate(q(unknowns))
    do k = 1, n
       do c = 1, ncells
          q(at_of(instant_moment(k), c) + 1:at_of(instant_moment(k), c) + stride) = tuple(:, c, first + k - 1)
          if (k < n .and. staged) then
             do i = 1, s
                if (solved_stages) then
                   q(at_of(stage_moment(k, i), c) + 1:at_of(stage_moment(k, i), c) + stride) = &
                        & stages(:, c, i, first + k)
                else
                   q(at_of(stage_moment(k, i), c) + 1:at_of(stage_moment(k, i), c) + stride) = &
                        & tuple(:, c, first + k - 1)
                end if
             end do
          end if
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

  end subroutine block_rows

  !===================================================================!
  ! ONE BLOCK SOLVED: its rows formed at the tuples, Newton from the
  ! estimate, the block's own instants written back.
  !===================================================================!

  subroutine solve_block(this, scheme, first, last, depth, rows_along, tuple, sbase, stages, tjet, sjet, expand)

    class(discrete_residual), intent(in)    :: this
    type(family)            , intent(in)    :: scheme
    integer                 , intent(in)    :: first, last, depth
    type(stencil)           , intent(in)    :: rows_along(:,:)
    real(dp)                , intent(inout) :: tuple(:,:,:)
    integer                 , intent(in)    :: sbase
    real(dp)                , intent(inout) :: stages(:,:,:,sbase:), tjet(:,:,:,:,:), sjet(:,:,:,sbase:,:,:)
    logical                 , intent(in), optional :: expand

    type(residual_operator) :: rows
    real(dp), allocatable :: q(:)
    integer , allocatable :: at(:), instant_moment(:), primary(:)
    integer :: stride, ncells, n, k, c, base, i, s, width

    stride = this % rule % num_components()
    ncells = this % points % num_cells()
    n      = last - first + 1
    if (verbosity >= 1) then
       print '(a,a,a,i0,a,i0,a,es12.4,a,es12.4,a,i0,a)', 'block  ', trim(scheme % name()), '  instants ', first, &
            & '..', last, '  t = ', this % points % instant(first), ' .. ', this % points % instant(last), &
            & '  history ', depth, ' instants'
    end if

    call block_rows(this, scheme, first, last, depth, rows_along, tuple, sbase, stages, .false., &
         & rows, q, at, instant_moment, primary)
    call solved(this, rows, size(q), stride, size(at), q)

    ! the block's own instants, and the stages of every step, retained
    ! for the adjoint, which linearizes the block at its solution
    do k = depth + 1, n
       do c = 1, ncells
          base = at((instant_moment(k) - 1) * ncells + c)
          tuple(:, c, first + k - 1) = q(base + 1:base + stride)
       end do
    end do
    if (marches_by_stages(scheme, top_degree(this % rule, this % manifold % num_unknowns) + 1)) then
       s     = scheme % num_stages()
       width = stride * ncells
       do k = 1, n - 1
          do i = 1, s
             ! the moments run instant, then the stages of the step, then the instant ahead
             do c = 1, ncells
                base = (instant_moment(k) + i - 1) * width + (c - 1) * stride
                stages(:, c, i, first + k) = q(base + 1:base + stride)
             end do
          end do
       end do
    end if

    if (this % points % design_order() > 0) then
       if (present(expand)) then
          if (.not. expand) return
       end if
       call design_jets(this, rows, q, at, instant_moment, first, last, depth, tjet, sbase, sjet, scheme)
    end if

  end subroutine solve_block

  !===================================================================!
  ! THE JETS ALONG THE DESIGN of one block, coordinate by coordinate
  ! and order by order: the m-th derivative of the block's unknowns
  ! along the design coordinate d is the solution of A x = b, A the
  ! block's jacobian at its solution and b the m-th derivative along d
  ! of the rows with the m-th derivative of the state left out - at
  ! the row of an equation, by the jet arithmetic on the derivatives
  ! of lower order and the design; at the fixed row of a history
  ! instant, the derivative solved before; at the fixed row of a face
  ! condition, the m-th derivative of the value it states; zero at
  ! the rows along time and space, which are linear in the state.
  ! Every derivative of the block's own instants is stored at
  ! tjet(:, c, k, m, d), and of the stages of the step into instant k
  ! at sjet(:, c, i, k, m, d).
  !===================================================================!

  subroutine design_jets(this, rows, q, at, instant_moment, first, last, depth, tjet, sbase, sjet, scheme)

    class(discrete_residual), intent(in)    :: this
    type(residual_operator) , intent(in)    :: rows
    real(dp)                , intent(in)    :: q(:)
    integer                 , intent(in)    :: at(:), instant_moment(:), first, last, depth
    real(dp)                , intent(inout) :: tjet(:,:,:,:,:)
    integer                 , intent(in)    :: sbase
    real(dp)                , intent(inout) :: sjet(:,:,:,sbase:,:,:)
    type(family)            , intent(in), optional :: scheme

    type(newton) :: solver
    type(stored_field), allocatable :: inputs(:)
    real(dp), allocatable :: jets(:,:), x(:), y(:), values(:), seed(:)
    integer , allocatable :: owner(:)
    integer :: stride, ncells, n, order, unknowns, npts, m, k, c, base, i, d, nd
    logical :: staged

    stride   = this % rule % num_components()
    ncells   = this % points % num_cells()
    n        = last - first + 1
    order    = this % points % design_order()
    nd       = this % points % num_expansions()
    values   = this % points % design_values()
    unknowns = size(q)
    npts     = size(at)
    staged   = .false.
    if (present(scheme)) staged = marches_by_stages(scheme, top_degree(this % rule, this % manifold % num_unknowns) + 1)

    solver = solver_for(this, rows, unknowns, stride, npts)
    call solver % evaluate(q, y, inputs)
    owner = point_owners(this, npts)

    do d = 1, nd
       if (.not. expansion_owned(this, d)) cycle
       seed = this % points % expansion_seed(d)
       allocate(jets(unknowns, 0:order), source=0.0_dp)
       jets(:, 0) = q
       do k = 1, depth
          do c = 1, ncells
             base = at((instant_moment(k) - 1) * ncells + c)
             jets(base + 1:base + stride, 1:order) = tjet(:, c, first + k - 1, 1:order, d)
          end do
       end do

       do m = 1, order
          if (verbosity >= 1) print '(a,i0,a,i0)', 'expansion ', d, ' along the design, order ', m
          call expansion_order(this, rows, inputs, q, at, instant_moment, first, last, depth, m, values, seed, &
               & owner, jets(:, 0:m), x)
          jets(:, m) = x
       end do

       do k = depth + 1, n
          do c = 1, ncells
             base = at((instant_moment(k) - 1) * ncells + c)
             tjet(:, c, first + k - 1, 1:order, d) = jets(base + 1:base + stride, 1:order)
          end do
       end do
       if (staged) then
          do k = 1, n - 1
             do i = 1, scheme % num_stages()
                do c = 1, ncells
                   base = at((instant_moment(k) + i - 1) * ncells + c)
                   sjet(:, c, i, first + k, 1:order, d) = jets(base + 1:base + stride, 1:order)
                end do
             end do
          end do
       end if
       deallocate(jets)
    end do

  end subroutine design_jets

  !===================================================================!
  ! ONE ORDER of the expansion of a block: the m-th derivative of the
  ! block's unknowns along the seed, from jets(:, 0:m) - the state
  ! and the derivatives below m at every unknown, and the m-th at the
  ! history instants, solved before - as the solution x of A x = b
  ! described above, the jacobian at the block's solution bound in
  ! inputs.
  !===================================================================!

  subroutine expansion_order(this, rows, inputs, q, at, instant_moment, first, last, depth, m, values, seed, &
       & owner, jets, x)

    class(discrete_residual), intent(in)  :: this
    type(residual_operator) , intent(in)  :: rows
    type(stored_field)      , intent(in)  :: inputs(:)
    real(dp)                , intent(in)  :: q(:)
    integer                 , intent(in)  :: at(:), instant_moment(:), first, last, depth, m
    real(dp)                , intent(in)  :: values(:), seed(:)
    integer                 , intent(in)  :: owner(:)
    real(dp)                , intent(in)  :: jets(:,0:)
    real(dp), allocatable   , intent(out) :: x(:)

    type(residual_operator) :: lin
    real(dp), allocatable :: b(:), slopes(:)
    integer , allocatable :: fields(:), degrees(:)
    real(dp) :: given
    integer :: stride, ncells, unknowns, npts, k, c, base, j, i, instant, f, degree

    stride   = this % rule % num_components()
    ncells   = this % points % num_cells()
    unknowns = size(q)
    npts     = size(at)
    allocate(b(unknowns), source=0.0_dp)
    call design_source(rows, jets(:, 0:m - 1), m, values, seed, owner, b)
    do k = 1, depth
       do c = 1, ncells
          base = at((instant_moment(k) - 1) * ncells + c)
          b(base + 1:base + stride) = jets(base + 1:base + stride, m)
       end do
    end do
    do j = 1, size(this % condition)
       if (.not. this % condition(j) % manifold % is_face()) cycle
       instant = face_instant(this, this % condition(j) % manifold % face_time)
       if (instant < first + depth .or. instant > last) cycle
       do i = 1, size(this % condition(j) % equation)
          do c = 1, ncells
             call face_relation(this % condition(j) % equation(i), values, &
                  & this % points % position(this % points % point_of(instant, c)), fields, degrees, slopes, &
                  & given, f, degree)
             base = at((face_row_moment(this, instant, first, last, depth, instant_moment) - 1) * ncells + c)
             b(base + this % rule % offset_of_field(f) + degree + 1) = stated_derivative( &
                  & this % condition(j) % equation(i), m, values, seed, &
                  & this % points % position(this % points % point_of(instant, c)))
          end do
       end do
    end do
    lin = rows % linearize(rows % unknown_graph(), rows % bind(inputs), b, transposed=.false., &
         & version_number=rows % version())
    allocate(x(unknowns), source=0.0_dp)
    call solved(this, lin, unknowns, stride, npts, x)

  end subroutine expansion_order

  !===================================================================!
  ! THE WAVEFRONT over the blocks and the orders: on a manifold of
  ! time alone with one expansion, when the images number the orders
  ! plus one, image m + 1 forms order m of every block, one block
  ! behind image m. Order m of block b reads the state of block b from
  ! image 1, the orders below m at block b from the images that formed
  ! them, and its own order m at the block's history instants, formed
  ! at the block before; each block ends with a rendezvous with the
  ! next image, so that every image works while the others do. At the
  ! end every image keeps its own order's jets, zero elsewhere, for
  ! the sum over the images that minimize takes.
  !===================================================================!

  subroutine wavefront_expansion(this, rows_along, tuple, stages, tjet, sjet)

    class(discrete_residual), intent(in)    :: this
    type(stencil)           , intent(in)    :: rows_along(:,:)
    real(dp)                , intent(inout) :: tuple(:,:,:), stages(:,:,:,:), tjet(:,:,:,:,:), sjet(:,:,:,:,:,:)

    real(dp), allocatable :: tuple_co(:,:,:)[:], stages_co(:,:,:,:)[:], tjet_co(:,:,:,:)[:], sjet_co(:,:,:,:,:)[:]
    real(dp), allocatable :: lower(:,:,:,:), slower(:,:,:,:,:)
    type(residual_operator) :: rows
    type(newton) :: solver
    type(stored_field), allocatable :: inputs(:)
    real(dp), allocatable :: q(:), y(:), jets(:,:), x(:), values(:), seed(:)
    integer , allocatable :: at(:), instant_moment(:), primary(:), owner(:)
    integer :: m, me, b, nb, first, last, depth, n, k, c, i, o, base, stride, ncells, ninst, order, smax, s
    integer :: unknowns, npts
    logical :: staged

    stride = this % rule % num_components()
    ncells = this % points % num_cells()
    ninst  = this % points % num_instants()
    order  = this % points % design_order()
    smax   = size(stages, 3)
    me     = this_image()
    m      = me - 1
    values = this % points % design_values()
    seed   = this % points % expansion_seed(1)
    nb     = this % schemes % num_blocks()

    allocate(tuple_co(stride, ncells, ninst)[*], stages_co(stride, ncells, smax, ninst)[*])
    allocate(tjet_co(stride, ncells, ninst, order)[*], sjet_co(stride, ncells, smax, ninst, order)[*])
    tuple_co = 0.0_dp
    stages_co = 0.0_dp
    tjet_co = 0.0_dp
    sjet_co = 0.0_dp
    allocate(lower(stride, ncells, ninst, order), slower(stride, ncells, smax, ninst, order), source=0.0_dp)

    do b = 1, nb
       call block_extent(this, b, first, last, depth)
       n = last - first + 1
       associate (scheme => this % schemes % scheme(b))
         staged = marches_by_stages(scheme, top_degree(this % rule, this % manifold % num_unknowns) + 1)
         s = 0
         if (staged) s = scheme % num_stages()
         if (m == 0) then
            ! the state, published to the image of order one
            call solve_block(this, scheme, first, last, depth, rows_along, tuple, 1, stages, tjet, sjet, expand=.false.)
            tuple_co(:, :, first + depth:last) = tuple(:, :, first + depth:last)
            if (staged) stages_co(:, :, 1:s, first + 1:last) = stages(:, :, 1:s, first + 1:last)
            sync images(2)
         else
            sync images(me - 1)
            tuple(:, :, first:last) = tuple_co(:, :, first:last)[1]
            if (staged) stages(:, :, 1:s, first + 1:last) = stages_co(:, :, 1:s, first + 1:last)[1]
            do o = 1, m - 1
               lower(:, :, first:last, o) = tjet_co(:, :, first:last, o)[o + 1]
               if (staged) slower(:, :, 1:s, first + 1:last, o) = sjet_co(:, :, 1:s, first + 1:last, o)[o + 1]
            end do
            if (verbosity >= 1) print '(a,i0,a,i0,a,i0)', 'wavefront: image ', me, ' forms order ', m, ' of block ', b
            call block_rows(this, scheme, first, last, depth, rows_along, tuple, 1, stages, .true., &
                 & rows, q, at, instant_moment, primary)
            unknowns = size(q)
            npts     = size(at)
            solver = solver_for(this, rows, unknowns, stride, npts)
            call solver % evaluate(q, y, inputs)
            owner = point_owners(this, npts)
            allocate(jets(unknowns, 0:m), source=0.0_dp)
            jets(:, 0) = q
            do k = 1, n
               do c = 1, ncells
                  base = at((instant_moment(k) - 1) * ncells + c)
                  jets(base + 1:base + stride, 1:m - 1) = lower(:, c, first + k - 1, 1:m - 1)
                  if (k <= depth) jets(base + 1:base + stride, m) = tjet(:, c, first + k - 1, m, 1)
                  if (k < n .and. staged) then
                     do i = 1, s
                        base = at((instant_moment(k) + i - 1) * ncells + c)
                        jets(base + 1:base + stride, 1:m - 1) = slower(:, c, i, first + k, 1:m - 1)
                     end do
                  end if
               end do
            end do
            call expansion_order(this, rows, inputs, q, at, instant_moment, first, last, depth, m, values, seed, &
                 & owner, jets, x)
            do k = depth + 1, n
               do c = 1, ncells
                  base = at((instant_moment(k) - 1) * ncells + c)
                  tjet(:, c, first + k - 1, m, 1)  = x(base + 1:base + stride)
                  tjet_co(:, c, first + k - 1, m)  = x(base + 1:base + stride)
               end do
            end do
            if (staged) then
               do k = 1, n - 1
                  do i = 1, s
                     do c = 1, ncells
                        base = at((instant_moment(k) + i - 1) * ncells + c)
                        sjet(:, c, i, first + k, m, 1) = x(base + 1:base + stride)
                        sjet_co(:, c, i, first + k, m) = x(base + 1:base + stride)
                     end do
                  end do
               end do
            end if
            deallocate(jets)
            if (me < num_images()) sync images(me + 1)
         end if
       end associate
    end do
    sync all

  end subroutine wavefront_expansion

  !===================================================================!
  ! Whether the blocks and the orders are pipelined over the images:
  ! a chain on a manifold of time alone, one expansion, and as many
  ! images as orders plus one.
  !===================================================================!

  function pipelined(this)

    class(discrete_residual), intent(in) :: this
    logical :: pipelined
    integer :: order

    order = this % points % design_order()
    pipelined = this % with_chain .and. .not. this % points % with_space .and. order >= 1 &
         & .and. this % points % num_expansions() == 1 .and. num_images() == order + 1

  end function pipelined

  !===================================================================!
  ! THE WAVEFRONT IN REVERSE over the blocks and the orders of the
  ! adjoint: image o + 1 forms order o of the multipliers of every
  ! block from the last, one block behind image o. Order o of block b
  ! reads the multipliers of the orders below o at block b from the
  ! images that formed them, and its own order at the block's history
  ! rows from block b + 1, formed before; the reactions of the orders
  ! below enter the design conditions' multipliers. Every image leaves
  ! the multipliers of the other orders zero, for the sum over the
  ! images the sweep takes.
  !===================================================================!

  subroutine wavefront_adjoint(this, rows_along, tuple, stages, tjet, sjet, gradient, stage_node, &
       & adjoint, incoming, reaction, gauge, in_design)

    class(discrete_residual), intent(in)    :: this
    type(stencil)           , intent(in)    :: rows_along(:,:)
    real(dp)                , intent(in)    :: tuple(:,:,:), stages(:,:,:,:), tjet(:,:,:,:,:), sjet(:,:,:,:,:,:)
    real(dp)                , intent(in)    :: gradient(:,:,0:,:)
    integer                 , intent(in)    :: stage_node(:,:,:)
    real(dp)                , intent(inout) :: adjoint(:,:,:,0:,:), incoming(:,:,:,0:,:), reaction(:,:,0:,:)
    real(dp)                , intent(inout) :: gauge(:,:,0:,:), in_design(:,0:,:)

    real(dp), allocatable :: mu_co(:,:,:)[:], smu_co(:,:,:,:)[:], reaction_co(:,:,:)[:]
    type(residual_operator) :: rows
    type(newton) :: solver
    type(stored_field), allocatable :: inputs(:)
    real(dp), allocatable :: q(:), y(:), jets(:,:), mu(:,:), x(:), values(:), seed(:)
    integer , allocatable :: at(:), instant_moment(:), primary(:), owner(:)
    logical , allocatable :: gauged(:)
    integer :: o, me, b, nb, first, last, depth, n, k, c, i, j, base, stride, ncells, ninst, order, smax, s, nf
    integer :: unknowns, npts, nc
    logical :: staged

    stride = this % rule % num_components()
    nf     = this % manifold % num_unknowns
    ncells = this % points % num_cells()
    ninst  = this % points % num_instants()
    order  = this % points % design_order()
    smax   = size(stages, 3)
    nc     = size(reaction, 2)
    me     = this_image()
    o      = me - 1
    values = this % points % design_values()
    seed   = this % points % expansion_seed(1)
    gauged = gauged_unknowns(this)
    nb     = this % schemes % num_blocks()

    allocate(mu_co(stride, ncells, ninst)[*], smu_co(stride, ncells, smax, ninst)[*], reaction_co(ncells, nc, 0:order)[*])
    mu_co = 0.0_dp
    smu_co = 0.0_dp
    reaction_co = 0.0_dp

    do b = nb, 1, -1
       call block_extent(this, b, first, last, depth)
       n = last - first + 1
       associate (scheme => this % schemes % scheme(b))
         staged = marches_by_stages(scheme, top_degree(this % rule, nf) + 1)
         s = 0
         if (staged) s = scheme % num_stages()
         if (o > 0) sync images(me - 1)
         if (verbosity >= 1) print '(a,i0,a,i0,a,i0)', 'wavefront: image ', me, ' forms the adjoint of order ', o, &
              & ' of block ', b
         call block_rows(this, scheme, first, last, depth, rows_along, tuple, 1, stages, .true., &
              & rows, q, at, instant_moment, primary)
         unknowns = size(q)
         npts     = size(at)
         solver = solver_for(this, rows, unknowns, stride, npts)
         call solver % evaluate(q, y, inputs)
         owner = point_owners(this, npts)
         ! the jets of the state to order o, and the multipliers of the
         ! orders below at the block's own moments
         allocate(jets(unknowns, 0:o), mu(unknowns, 0:o), source=0.0_dp)
         jets(:, 0) = q
         do k = 1, n
            do c = 1, ncells
               base = at((instant_moment(k) - 1) * ncells + c)
               jets(base + 1:base + stride, 1:o) = tjet(:, c, first + k - 1, 1:o, 1)
               if (k > depth) then
                  do j = 1, o
                     mu(base + 1:base + stride, j - 1) = mu_co(:, c, first + k - 1)[j]
                  end do
               end if
               if (k < n .and. staged) then
                  do i = 1, s
                     base = at((instant_moment(k) + i - 1) * ncells + c)
                     jets(base + 1:base + stride, 1:o) = sjet(:, c, i, first + k, 1:o, 1)
                     do j = 1, o
                        mu(base + 1:base + stride, j - 1) = smu_co(:, c, i, first + k)[j]
                     end do
                  end do
               end if
            end do
         end do
         do j = 1, o
            reaction(:, :, j - 1, 1) = reaction_co(:, :, j - 1)[j]
         end do
         call adjoint_order(this, rows, inputs, at, instant_moment, first, last, depth, o, 1, values, seed, owner, &
              & jets, mu(:, 0:max(o - 1, 0)), gradient, stage_node, incoming, s, x)
         mu(:, o) = x
         call incoming_added(this, at, instant_moment, first, depth, o, 1, 1, x, incoming)
         do k = depth + 1, n
            do c = 1, ncells
               base = at((instant_moment(k) - 1) * ncells + c)
               mu_co(:, c, first + k - 1) = x(base + 1:base + stride)
            end do
         end do
         if (staged) then
            do k = 1, n - 1
               do i = 1, s
                  do c = 1, ncells
                     base = at((instant_moment(k) + i - 1) * ncells + c)
                     smu_co(:, c, i, first + k) = x(base + 1:base + stride)
                  end do
               end do
            end do
         end if
         call multipliers_of_order(this, rows, first, last, depth, at, instant_moment, primary, o, 1, values, seed, &
              & owner, gauged, s, jets, mu, adjoint, gauge, reaction, in_design)
         reaction_co(:, :, o) = reaction(:, :, o, 1)
         deallocate(jets, mu, x)
         if (me < num_images()) sync images(me + 1)
       end associate
    end do
    ! the reactions of the other orders, read from the other images,
    ! left to the sum
    do j = 0, order
       if (j /= o) reaction(:, :, j, 1) = 0.0_dp
    end do
    sync all

  end subroutine wavefront_adjoint

  !===================================================================!
  ! The m-th derivative along the design coordinate d of every
  ! equation row at the points not fixed, with the m-th derivative of
  ! the state left out, negated: each component of a point's tuple
  ! enters as the quantity whose k-th derivative is the k-th jet
  ! coefficient for k below m and zero at m, the coordinate d as the
  ! quantity with derivative one; the coefficient of the full subset
  ! of m directions is the m-th derivative by the product rule on
  ! subsets.
  !===================================================================!

  subroutine design_source(rows, jets, m, values, seed, owner, b)

    type(residual_operator), intent(in)    :: rows
    real(dp)               , intent(in)    :: jets(:,0:)
    integer                , intent(in)    :: m
    real(dp)               , intent(in)    :: values(:), seed(:)
    integer                , intent(in)    :: owner(:)
    real(dp)               , intent(inout) :: b(:)

    type(derivative_terms), allocatable :: q(:), design(:)
    type(derivative_terms) :: r
    logical, allocatable :: is_fixed(:)
    real(dp), allocatable :: source(:)
    integer :: deg, npts, j, p, c, k, row

    is_fixed = rows % fixed_indicator()
    deg      = rows % degrees
    npts     = size(rows % at)
    design   = design_terms(values, seed, m, m)
    allocate(q(0:deg - 1))
    allocate(source(size(b)), source=0.0_dp)
    do j = 1, size(rows % rules)
       do p = 1, npts
          if (owner(p) /= this_image()) cycle
          row = rows % at(p) + rows % primary(j) + 1
          if (is_fixed(row)) cycle
          do c = 0, deg - 1
             q(c) = derivative_terms(jets(rows % at(p) + c + 1, 0), m)
             do k = 1, m - 1
                call q(c) % set_symmetric(k, jets(rows % at(p) + c + 1, k))
             end do
          end do
          r = rows % rules(j) % at_instant(q, design)
          source(row) = -mixed_partial(r)
       end do
    end do
    if (any(owner /= this_image())) call co_sum(source)
    b = b + source

  end subroutine design_source

  !===================================================================!
  ! The m-th derivative along the design coordinate d of the value a
  ! face condition states: the condition g(u) = u - h(x, nu)
  ! evaluated on a state without derivatives, so that the m-th
  ! coefficient is -h^(m), negated.
  !===================================================================!

  real(dp) function stated_derivative(equation, m, values, seed, position)

    type(continuous_field), intent(in) :: equation
    integer               , intent(in) :: m
    real(dp)              , intent(in) :: values(:), seed(:)
    real(dp)              , intent(in) :: position(:)

    type(expression) :: g
    type(derivative_terms), allocatable :: q(:)

    g = equation % graph(1)
    allocate(q(0:g % num_components() - 1), source=derivative_terms(0.0_dp, m))
    stated_derivative = -mixed_partial(g % at_instant(q, design_terms(values, seed, m, m), position))

  end function stated_derivative

  !===================================================================!
  ! The derivatives along the design coordinate d, to order o, of the
  ! partial derivative in the coordinate i of the value a face
  ! condition states: slope(k) = d^k/dnu_d^k (dh/dnu_i), from one
  ! evaluation over o + 1 directions, the last seeded on i.
  !===================================================================!

  function stated_partial_jet(equation, i, seed, o, values, position) result(slope)

    type(continuous_field), intent(in) :: equation
    integer               , intent(in) :: i, o
    real(dp)              , intent(in) :: seed(:), values(:)
    real(dp)              , intent(in) :: position(:)
    real(dp), allocatable :: slope(:)

    type(expression) :: g
    type(derivative_terms), allocatable :: q(:), nu(:)
    type(derivative_terms) :: r
    integer :: k

    g  = equation % graph(1)
    nu = design_terms(values, seed, o, o + 1)
    call nu(i) % set_coefficient(2**o, 1.0_dp)
    allocate(q(0:g % num_components() - 1), source=derivative_terms(0.0_dp, o + 1))
    r = g % at_instant(q, nu, position)
    allocate(slope(0:o))
    do k = 0, o
       slope(k) = -coefficient(r, 2**k - 1 + 2**o)
    end do

  end function stated_partial_jet

  !===================================================================!
  ! ONE BLOCK'S ADJOINT, at order zero and then order by order along
  ! each design coordinate: the transposed linear system A^T mu = rhs
  ! at the block's solution, A the block's jacobian; rhs at every
  ! component of the block's own instants is the negative of the
  ! objective's differential there (its o-th derivative along the
  ! coordinate at order o) plus the multipliers the later blocks
  ! passed to the instant, less, above order zero, the o-th derivative
  ! of the product of the multipliers with the equations' jacobian
  ! with the o-th derivative of the multipliers left out; zero at the
  ! history instants. mu at the row of an equation is the multiplier
  ! of that equation at the instant; mu at the fixed row of a history
  ! instant is what this block passes to the block that determines
  ! the instant; mu at the fixed row of a face condition is the
  ! condition's multiplier, the reaction. The block's rows' partial
  ! derivatives in each design coordinate, each times its multiplier
  ! and differentiated along the coordinate of the expansion, are
  ! subtracted from the design conditions' multipliers.
  !===================================================================!

  subroutine adjoint_block(this, scheme, first, last, depth, rows_along, tuple, sbase, stages, tjet, sjet, gradient, &
       & stage_node, adjoint, incoming, reaction, gauge, in_design)

    class(discrete_residual), intent(in)    :: this
    type(family)            , intent(in)    :: scheme
    integer                 , intent(in)    :: first, last, depth
    type(stencil)           , intent(in)    :: rows_along(:,:)
    real(dp)                , intent(in)    :: tuple(:,:,:)
    integer                 , intent(in)    :: sbase
    real(dp)                , intent(in)    :: stages(:,:,:,sbase:), tjet(:,:,:,:,:), sjet(:,:,:,sbase:,:,:)
    real(dp)                , intent(in)    :: gradient(:,:,0:,:)
    integer                 , intent(in)    :: stage_node(:,:,:)
    real(dp)                , intent(inout) :: adjoint(:,:,:,0:,:), incoming(:,:,:,0:,:), reaction(:,:,0:,:)
    real(dp)                , intent(inout) :: gauge(:,:,0:,:), in_design(:,0:,:)

    type(residual_operator) :: rows
    real(dp), allocatable :: q(:)
    integer , allocatable :: at(:), instant_moment(:), primary(:)

    if (verbosity >= 1) then
       print '(a,a,a,i0,a,i0,a,es12.4,a,es12.4)', 'adjoint  ', trim(scheme % name()), '  instants ', &
            & first + depth, '..', last, '  t = ', this % points % instant(first + depth), ' .. ', &
            & this % points % instant(last)
    end if
    call block_rows(this, scheme, first, last, depth, rows_along, tuple, sbase, stages, .true., &
         & rows, q, at, instant_moment, primary)
    call adjoint_rows(this, first, last, depth, rows, q, at, instant_moment, primary, tjet, sbase, sjet, gradient, &
         & stage_node, adjoint, incoming, reaction, gauge, in_design, scheme)

  end subroutine adjoint_block

  !===================================================================!
  ! THE ADJOINT OF ONE SYSTEM OF ROWS - a block's, or the point
  ! manifold's, which has one moment, no history and no stages - as
  ! described above.
  !===================================================================!

  subroutine adjoint_rows(this, first, last, depth, rows, q, at, instant_moment, primary, tjet, sbase, sjet, gradient, &
       & stage_node, adjoint, incoming, reaction, gauge, in_design, scheme)

    class(discrete_residual), intent(in)    :: this
    integer                 , intent(in)    :: first, last, depth
    type(residual_operator) , intent(in)    :: rows
    real(dp)                , intent(in)    :: q(:)
    integer                 , intent(in)    :: at(:), instant_moment(:), primary(:)
    real(dp)                , intent(in)    :: tjet(:,:,:,:,:)
    integer                 , intent(in)    :: sbase
    real(dp)                , intent(in)    :: sjet(:,:,:,sbase:,:,:)
    real(dp)                , intent(in)    :: gradient(:,:,0:,:)
    integer                 , intent(in)    :: stage_node(:,:,:)
    real(dp)                , intent(inout) :: adjoint(:,:,:,0:,:), incoming(:,:,:,0:,:), reaction(:,:,0:,:)
    real(dp)                , intent(inout) :: gauge(:,:,0:,:), in_design(:,0:,:)
    type(family)            , intent(in), optional :: scheme

    type(newton) :: solver
    type(stored_field), allocatable :: inputs(:)
    real(dp), allocatable :: y(:), mu(:,:,:), x(:), jets(:,:,:), values(:)
    integer , allocatable :: owner(:)
    logical , allocatable :: gauged(:)
    integer :: stride, nf, ncells, n, k, c, base, unknowns, npts, i, s, o, order, d, nd, first_order
    logical :: staged

    stride = this % rule % num_components()
    nf     = this % manifold % num_unknowns
    ncells = this % points % num_cells()
    n      = last - first + 1
    order  = this % points % design_order()
    nd     = max(1, this % points % num_expansions())
    values = this % points % design_values()
    unknowns = size(q)
    npts     = size(at)
    staged   = .false.
    if (present(scheme)) staged = marches_by_stages(scheme, top_degree(this % rule, nf) + 1)
    s        = 0
    if (staged) s = scheme % num_stages()

    ! the jets of the state along each design coordinate at every
    ! unknown of the block: at the instants from tjet, at the stages
    ! from sjet
    allocate(jets(unknowns, 0:order, nd), source=0.0_dp)
    do d = 1, nd
       jets(:, 0, d) = q
       do k = 1, n
          do c = 1, ncells
             base = at((instant_moment(k) - 1) * ncells + c)
             jets(base + 1:base + stride, 1:order, d) = tjet(:, c, first + k - 1, 1:order, d)
             if (k < n .and. staged) then
                do i = 1, s
                   base = at((instant_moment(k) + i - 1) * ncells + c)
                   jets(base + 1:base + stride, 1:order, d) = sjet(:, c, i, first + k, 1:order, d)
                end do
             end if
          end do
       end do
    end do

    gauged = gauged_unknowns(this)
    solver = solver_for(this, rows, unknowns, stride, npts)
    call solver % evaluate(q, y, inputs)
    owner = point_owners(this, npts)
    allocate(mu(unknowns, 0:order, nd), source=0.0_dp)

    do d = 1, nd
       first_order = 0
       if (d > 1) first_order = 1
       do o = first_order, order
          if (o > 0 .and. .not. expansion_owned(this, d)) cycle
          call adjoint_order(this, rows, inputs, at, instant_moment, first, last, depth, o, d, values, &
               & this % points % expansion_seed(d), owner, jets(:, 0:o, d), mu(:, 0:max(o - 1, 0), d), gradient, &
               & stage_node, incoming, s, x)
          if (o == 0) then
             mu(:, 0, :) = spread(x, 2, nd)
          else
             mu(:, o, d) = x
          end if
          call incoming_added(this, at, instant_moment, first, depth, o, d, nd, x, incoming)
          deallocate(x)
       end do
    end do

    do d = 1, nd
       do o = 0, order
          if (o > 0 .and. .not. expansion_owned(this, d)) cycle
          call multipliers_of_order(this, rows, first, last, depth, at, instant_moment, primary, o, d, values, &
               & this % points % expansion_seed(d), owner, gauged, s, jets(:, 0:o, d), mu(:, 0:o, d), &
               & adjoint, gauge, reaction, in_design)
       end do
    end do

  end subroutine adjoint_rows

  !===================================================================!
  ! The unknowns a gauge replaces the first cell's row of.
  !===================================================================!

  function gauged_unknowns(this) result(gauged)

    class(discrete_residual), intent(in) :: this
    logical, allocatable :: gauged(:)
    integer :: i, j

    allocate(gauged(this % manifold % num_unknowns), source=.false.)
    do j = 1, size(this % condition)
       if (.not. this % condition(j) % manifold % is_time_factor()) cycle
       do i = 1, size(this % condition(j) % equation)
          gauged(this % condition(j) % equation(i) % index(1)) = .true.
       end do
    end do

  end function gauged_unknowns

  !===================================================================!
  ! ONE ORDER of the adjoint of a block: the o-th derivative along the
  ! expansion d of the multipliers of the block's rows, the solution x
  ! of the transposed system whose right side is the objective's
  ! differential of order o at the block's nodes, the multipliers
  ! passed in from the block ahead at the block's own instants, and,
  ! above order zero, the derivatives of the transposed products with
  ! the multipliers of the orders below.
  !===================================================================!

  subroutine adjoint_order(this, rows, inputs, at, instant_moment, first, last, depth, o, d, values, seed, owner, &
       & jets, mu, gradient, stage_node, incoming, s, x)

    class(discrete_residual), intent(in)  :: this
    type(residual_operator) , intent(in)  :: rows
    type(stored_field)      , intent(in)  :: inputs(:)
    integer                 , intent(in)  :: at(:), instant_moment(:), first, last, depth, o, d, s
    real(dp)                , intent(in)  :: values(:), seed(:)
    integer                 , intent(in)  :: owner(:)
    real(dp)                , intent(in)  :: jets(:,0:), mu(:,0:), gradient(:,:,0:,:)
    integer                 , intent(in)  :: stage_node(:,:,:)
    real(dp)                , intent(in)  :: incoming(:,:,:,0:,:)
    real(dp), allocatable   , intent(out) :: x(:)

    type(residual_operator) :: lin
    real(dp), allocatable :: rhs(:)
    integer :: stride, ncells, n, k, c, i, base, unknowns, npts

    stride   = this % rule % num_components()
    ncells   = this % points % num_cells()
    n        = last - first + 1
    unknowns = size(jets, 1)
    npts     = size(at)
    if (verbosity >= 1 .and. o > 0) print '(a,i0,a,i0)', 'adjoint expansion ', d, ' along the design, order ', o
    allocate(rhs(unknowns), source=0.0_dp)
    do k = depth + 1, n
       do c = 1, ncells
          base = at((instant_moment(k) - 1) * ncells + c)
          rhs(base + 1:base + stride) = -gradient(1:stride, this % points % point_of(first + k - 1, c), o, d) &
               & + incoming(:, c, first + k - 1, o, d)
       end do
    end do
    ! the objective's differential at the stages, the seeds of the
    ! stage quadrature
    if (s > 0) then
       do k = 1, n - 1
          do i = 1, s
             do c = 1, ncells
                base = at((instant_moment(k) + i - 1) * ncells + c)
                rhs(base + 1:base + stride) = -gradient(1:stride, stage_node(first + k, i, c), o, d)
             end do
          end do
       end do
    end if
    if (o > 0) call adjoint_source(rows, jets, mu(:, 0:o - 1), o, values, seed, owner, rhs)
    lin = rows % linearize(rows % unknown_graph(), rows % bind(inputs), rhs, transposed=.true., &
         & version_number=rows % version())
    allocate(x(unknowns), source=0.0_dp)
    call solved(this, lin, unknowns, stride, npts, x)

  end subroutine adjoint_order

  !===================================================================!
  ! The multipliers of order o at the block's history rows, passed to
  ! the block that determines those instants; at order zero to every
  ! expansion.
  !===================================================================!

  subroutine incoming_added(this, at, instant_moment, first, depth, o, d, nd, x, incoming)

    class(discrete_residual), intent(in)    :: this
    integer                 , intent(in)    :: at(:), instant_moment(:), first, depth, o, d, nd
    real(dp)                , intent(in)    :: x(:)
    real(dp)                , intent(inout) :: incoming(:,:,:,0:,:)

    integer :: stride, ncells, k, c, base

    stride = this % rule % num_components()
    ncells = this % points % num_cells()
    do k = 1, depth
       do c = 1, ncells
          base = at((instant_moment(k) - 1) * ncells + c)
          if (o == 0) then
             incoming(:, c, first + k - 1, 0, :) = incoming(:, c, first + k - 1, 0, :) &
                  & + spread(x(base + 1:base + stride), 2, nd)
          else
             incoming(:, c, first + k - 1, o, d) = incoming(:, c, first + k - 1, o, d) + x(base + 1:base + stride)
          end if
       end do
    end do

  end subroutine incoming_added

  !===================================================================!
  ! THE MULTIPLIERS OF THE EQUATIONS of order o along the expansion d
  ! at the block's own instants: the multiplier at the equation's row
  ! at the instant and, for a family marching by stages, at the
  ! stages of the step into it; a row the gauge of a time-factor
  ! condition replaces, the first cell's row of the integrated
  ! unknown, stores the gauge's multiplier instead, which is returned
  ! as the condition's; the multiplier at the fixed row of a face
  ! condition is the condition's, the reaction.
  !
  ! THE DESIGN CONDITIONS' MULTIPLIERS: minus the o-th derivative
  ! along the expansion of the sum over the block's own rows of the
  ! multiplier times the row's partial in the coordinate i - at the
  ! equation rows the partial of the equation, at the fixed rows of
  ! the face conditions minus the derivative of the datum, from the
  ! reactions of the orders to o; the rows along time and space, and
  ! the history rows, which the block does not own, add nothing.
  !===================================================================!

  subroutine multipliers_of_order(this, rows, first, last, depth, at, instant_moment, primary, o, d, values, seed, &
       & owner, gauged, s, jets, mu, adjoint, gauge, reaction, in_design)

    class(discrete_residual), intent(in)    :: this
    type(residual_operator) , intent(in)    :: rows
    integer                 , intent(in)    :: first, last, depth, at(:), instant_moment(:), primary(:), o, d, s
    real(dp)                , intent(in)    :: values(:), seed(:)
    integer                 , intent(in)    :: owner(:)
    logical                 , intent(in)    :: gauged(:)
    real(dp)                , intent(in)    :: jets(:,0:), mu(:,0:)
    real(dp)                , intent(inout) :: adjoint(:,:,:,0:,:), gauge(:,:,0:,:), reaction(:,:,0:,:)
    real(dp)                , intent(inout) :: in_design(:,0:,:)

    real(dp), allocatable :: slope(:), products(:), slopes(:)
    integer , allocatable :: fields(:), degrees(:)
    real(dp) :: given
    integer :: stride, nf, ncells, n, k, c, f, base, m, j, i, instant, degree, e, g, ni

    stride = this % rule % num_components()
    nf     = this % manifold % num_unknowns
    ncells = this % points % num_cells()
    n      = last - first + 1
    ni     = size(values)

    do k = depth + 1, n
       do c = 1, ncells
          do f = 1, nf
             adjoint(f, c, first + k - 1, o, d) = 0.0_dp
             do m = instant_moment(k) - merge(s, 0, k > 1), instant_moment(k)
                if (c == 1 .and. gauged(f)) cycle
                base = at((m - 1) * ncells + c)
                adjoint(f, c, first + k - 1, o, d) = adjoint(f, c, first + k - 1, o, d) + mu(base + primary(f) + 1, o)
             end do
          end do
       end do
    end do
    g = 0
    do j = 1, size(this % condition)
       if (.not. this % condition(j) % manifold % is_time_factor()) cycle
       do i = 1, size(this % condition(j) % equation)
          g = g + 1
          f = this % condition(j) % equation(i) % index(1)
          do k = depth + 1, n
             do m = instant_moment(k) - merge(s, 0, k > 1), instant_moment(k)
                base = at((m - 1) * ncells + 1)
                gauge(first + k - 1, g, o, d) = gauge(first + k - 1, g, o, d) + mu(base + primary(f) + 1, o)
             end do
          end do
       end do
    end do
    e = 0
    do j = 1, size(this % condition)
       if (.not. this % condition(j) % manifold % is_face()) cycle
       instant = face_instant(this, this % condition(j) % manifold % face_time)
       do i = 1, size(this % condition(j) % equation)
          e = e + 1
          if (instant < first + depth .or. instant > last) cycle
          m = face_row_moment(this, instant, first, last, depth, instant_moment)
          do c = 1, ncells
             call face_relation(this % condition(j) % equation(i), values, &
                  & this % points % position(this % points % point_of(instant, c)), fields, degrees, slopes, &
                  & given, f, degree)
             base = at((m - 1) * ncells + c)
             reaction(c, e, o, d) = mu(base + this % rule % offset_of_field(f) + degree + 1, o)
          end do
       end do
    end do

    if (.not. this % points % with_design) return
    products = design_products(rows, jets, mu, o, values, seed, owner)
    in_design(:, o, d) = in_design(:, o, d) - products
    e = 0
    do j = 1, size(this % condition)
       if (.not. this % condition(j) % manifold % is_face()) cycle
       instant = face_instant(this, this % condition(j) % manifold % face_time)
       do i = 1, size(this % condition(j) % equation)
          e = e + 1
          if (instant < first + depth .or. instant > last) cycle
          do c = 1, ncells
             do k = 1, ni
                slope = stated_partial_jet(this % condition(j) % equation(i), k, seed, o, values, &
                     & this % points % position(this % points % point_of(instant, c)))
                ! the row x - h(nu): its partial in nu_k is -dh/dnu_k
                in_design(k, o, d) = in_design(k, o, d) + mixed_partial( &
                     & symmetric_terms(reaction(c, e, 0, d), reaction(c, e, 1:o, d), o) &
                     & * symmetric_terms(-slope(0), -slope(1:o), o))
             end do
          end do
       end do
    end do

  end subroutine multipliers_of_order

  !===================================================================!
  ! THE SOURCE OF THE ADJOINT EXPANSION at order o along the design
  ! coordinate d: for every equation row not fixed and every component
  ! of its point, the o-th derivative along d of the multiplier times
  ! the partial of the equation in the component, the multiplier's
  ! o-th derivative left out, subtracted from rhs at the component.
  ! The partial's derivatives along d come from one evaluation of the
  ! equation on the jets over o + 1 directions, the last seeded on the
  ! component: its coefficients on the first k directions with the
  ! last are the k-th derivatives of the partial.
  !===================================================================!

  subroutine adjoint_source(rows, jets, mu, o, values, seed, owner, rhs)

    type(residual_operator), intent(in)    :: rows
    real(dp)               , intent(in)    :: jets(:,0:), mu(:,0:)
    integer                , intent(in)    :: o
    real(dp)               , intent(in)    :: values(:), seed(:)
    integer                , intent(in)    :: owner(:)
    real(dp)               , intent(inout) :: rhs(:)

    type(derivative_terms), allocatable :: q(:), design(:)
    type(derivative_terms) :: r, along
    real(dp), allocatable :: partial(:), lower(:), source(:)
    logical, allocatable :: is_fixed(:)
    integer :: deg, npts, j, p, c, k, row

    is_fixed = rows % fixed_indicator()
    deg      = rows % degrees
    npts     = size(rows % at)
    allocate(lower(o), partial(0:o), source=0.0_dp)
    allocate(source(size(rhs)), source=0.0_dp)
    design = design_terms(values, seed, o, o + 1)
    allocate(q(0:deg - 1))
    do j = 1, size(rows % rules)
       do p = 1, npts
          if (owner(p) /= this_image()) cycle
          row = rows % at(p) + rows % primary(j) + 1
          if (is_fixed(row)) cycle
          do c = 0, deg - 1
             q(c) = symmetric_terms(jets(rows % at(p) + c + 1, 0), jets(rows % at(p) + c + 1, 1:o), o + 1)
          end do
          lower(1:o - 1) = mu(row, 1:o - 1)
          lower(o)       = 0.0_dp
          along = symmetric_terms(mu(row, 0), lower, o)
          do c = 0, deg - 1
             call q(c) % set_coefficient(2**o, 1.0_dp)
             r = rows % rules(j) % at_instant(q, design)
             call q(c) % set_coefficient(2**o, 0.0_dp)
             do k = 0, o
                partial(k) = coefficient(r, 2**k - 1 + 2**o)
             end do
             source(rows % at(p) + c + 1) = source(rows % at(p) + c + 1) &
                  & - mixed_partial(along * symmetric_terms(partial(0), partial(1:o), o))
          end do
       end do
    end do
    if (any(owner /= this_image())) call co_sum(source)
    rhs = rhs + source

  end subroutine adjoint_source

  !===================================================================!
  ! The o-th derivative along the design coordinate d of the sum over
  ! the equation rows not fixed of the multiplier times the partial
  ! of the equation in the design coordinate i, every jet complete to
  ! order o: at order zero by the reverse pass of every rule, all
  ! coordinates at once; above, the partial's derivatives from one
  ! evaluation on the jets over o + 1 directions, the last seeded on
  ! the coordinate i.
  !===================================================================!

  function design_products(rows, jets, mu, o, values, seed, owner) result(total)

    type(residual_operator), intent(in) :: rows
    real(dp)               , intent(in) :: jets(:,0:), mu(:,0:)
    integer                , intent(in) :: o
    real(dp)               , intent(in) :: values(:), seed(:)
    integer                , intent(in) :: owner(:)
    real(dp), allocatable :: total(:)

    type(derivative_terms), allocatable :: q(:), design(:)
    type(derivative_terms) :: r
    real(dp), allocatable :: partial(:), g(:), gd(:), x(:)
    logical, allocatable :: is_fixed(:)
    real(dp) :: value
    integer :: deg, npts, j, p, c, k, i, row

    is_fixed = rows % fixed_indicator()
    deg      = rows % degrees
    npts     = size(rows % at)
    allocate(total(size(values)), source=0.0_dp)
    if (o == 0) then
       allocate(g(0:deg - 1), gd(size(values)), x(0:deg - 1))
       do j = 1, size(rows % rules)
          do p = 1, npts
             if (owner(p) /= this_image()) cycle
             row = rows % at(p) + rows % primary(j) + 1
             if (is_fixed(row)) cycle
             x = jets(rows % at(p) + 1:rows % at(p) + deg, 0)
             call rows % rules(j) % gradient_at(x, values, value, g, g_design=gd)
             total = total + mu(row, 0) * gd
          end do
       end do
       if (any(owner /= this_image())) call co_sum(total)
       return
    end if
    allocate(partial(0:o), q(0:deg - 1))
    do i = 1, size(values)
       design = design_terms(values, seed, o, o + 1)
       call design(i) % set_coefficient(2**o, 1.0_dp)
       do j = 1, size(rows % rules)
          do p = 1, npts
             if (owner(p) /= this_image()) cycle
             row = rows % at(p) + rows % primary(j) + 1
             if (is_fixed(row)) cycle
             do c = 0, deg - 1
                q(c) = symmetric_terms(jets(rows % at(p) + c + 1, 0), jets(rows % at(p) + c + 1, 1:o), o + 1)
             end do
             r = rows % rules(j) % at_instant(q, design)
             do k = 0, o
                partial(k) = coefficient(r, 2**k - 1 + 2**o)
             end do
             total(i) = total(i) + mixed_partial(symmetric_terms(mu(row, 0), mu(row, 1:o), o) &
                  & * symmetric_terms(partial(0), partial(1:o), o))
          end do
       end do
    end do
    if (any(owner /= this_image())) call co_sum(total)

  end function design_products

  !===================================================================!
  ! The condition on the design factor at the design values: given =
  ! -g(nu_0) and slope(i) = dg/dnu_i, so that nu_k - nu_k0 = 0 gives
  ! zero and the k-th unit vector.
  !===================================================================!

  subroutine design_condition(equation, values, given, slope)

    type(continuous_field), intent(in)  :: equation
    real(dp)              , intent(in)  :: values(:)
    real(dp)              , intent(out) :: given
    real(dp), allocatable , intent(out) :: slope(:)

    type(expression) :: g
    type(derivative_terms), allocatable :: q(:), nu(:)
    type(derivative_terms) :: r
    integer :: i, nd

    g  = equation % graph(1)
    nd = size(values)
    allocate(nu(nd))
    do i = 1, nd
       nu(i) = derivative_terms(values(i), nd)
       call nu(i) % set_direction(i, 1.0_dp)
    end do
    allocate(q(0:g % num_components() - 1), source=derivative_terms(0.0_dp, nd))
    r     = g % at_instant(q, nu)
    given = -coefficient(r, 0)
    allocate(slope(nd))
    do i = 1, nd
       slope(i) = coefficient(r, 2**(i - 1))
    end do

  end subroutine design_condition

  !===================================================================!
  ! The point manifold: one moment, one cell, no rows along any
  ! coordinate; the equations alone.
  !===================================================================!

  subroutine point_rows(this, tuple, rows, q, primary)

    class(discrete_residual), intent(in)  :: this
    real(dp)                , intent(in)  :: tuple(:,:,:)
    type(residual_operator) , intent(out) :: rows
    real(dp), allocatable   , intent(out) :: q(:)
    integer , allocatable   , intent(out) :: primary(:)

    type(stencil) :: relations
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

  end subroutine point_rows

  subroutine solve_point(this, tuple, rows, q)

    class(discrete_residual), intent(in)    :: this
    real(dp)                , intent(inout) :: tuple(:,:,:)
    type(residual_operator) , intent(out)   :: rows
    real(dp), allocatable   , intent(out)   :: q(:)

    integer , allocatable :: primary(:)
    integer :: stride

    stride = this % rule % num_components()
    call point_rows(this, tuple, rows, q, primary)
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
    type(solve_result) :: outcome
    real(dp) :: achieved

    solver = solver_for(this, rows, unknowns, stride, npts)
    call solver % solve(spread(0.0_dp, 1, unknowns), q, achieved)
    outcome = solver % result()
    if (.not. outcome % converged()) then
       error stop 'operation_residual: minimize did not converge: ' // outcome % description()
    end if

  end subroutine solved

  !===================================================================!
  ! The design vector at every point, in the point order: the
  ! manifold's values, as many as the rules read, repeated per point.
  ! A rule reading a design the manifold does not have is invalid
  ! input.
  !===================================================================!

  function design_at_points(values, designs, npts) result(v)

    real(dp), intent(in) :: values(:)
    integer , intent(in) :: designs, npts
    real(dp), allocatable :: v(:)

    integer :: p

    if (designs > size(values)) then
       error stop 'operation_residual: the equations read a design coordinate the manifold does not have'
    end if
    allocate(v(designs * npts))
    do p = 1, npts
       v((p - 1) * designs + 1:p * designs) = values(1:designs)
    end do

  end function design_at_points

  !===================================================================!
  ! The jets of the design vector over a width of directions: each
  ! coordinate the quantity whose derivative on the first direction
  ! (symmetric over the expansion's) is its entry of the seed, the
  ! seed being the unit vector of a coordinate or a direction in the
  ! design.
  !===================================================================!

  function design_terms(values, seed, order, width) result(nu)

    real(dp), intent(in) :: values(:), seed(:)
    integer , intent(in) :: order, width
    type(derivative_terms), allocatable :: nu(:)

    real(dp), allocatable :: unit(:)
    integer :: i

    if (size(seed) /= size(values)) then
       error stop 'operation_residual: a seed of the expansion has one entry per design coordinate'
    end if
    allocate(unit(order), source=0.0_dp)
    if (order > 0) unit(1) = 1.0_dp
    allocate(nu(size(values)))
    do i = 1, size(values)
       nu(i) = symmetric_terms(values(i), seed(i) * unit, width)
    end do

  end function design_terms

  !===================================================================!
  ! The Newton solver stated on one residual_operator, its inner
  ! minimizer chosen by the manifold and its state set, so that the
  ! operator's jacobian can be formed at any q through it.
  !===================================================================!

  function solver_for(this, rows, unknowns, stride, npts) result(solver)

    class(discrete_residual), intent(in) :: this
    type(residual_operator) , intent(in) :: rows
    integer                 , intent(in) :: unknowns, stride, npts
    type(newton) :: solver

    type(dense_direct) :: direct
    type(gmres)        :: krylov
    type(gauss_seidel) :: sweeps
    type(typed_field_domain) :: designs
    type(stored_field) :: design

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
       ! linear solve with them; a build of one image under a
       ! launcher of several would run the whole problem once per
       ! process, and is refused
       call require_images()
       if (num_images() > 1) solver % distribution = exchange(owners_of_unknowns(this, unknowns, stride, npts))
    else
       ! a direct inner solve converges Newton quadratically to
       ! round-off at the cost of a step, and the sensitivities read
       ! the converged state: the tolerance is tight here
       allocate(solver % inner, source=direct)
       solver % tolerance = 1.0e-13_dp
    end if

    designs = rows % design_fields()
    design  = designs % design(design_at_points(this % points % design_values(), rows % designs, npts))
    call solver % state(rows, rows % unknown_graph(), rows % unknown_domain(), unknowns, stored_inputs=[design])

  end function solver_for

  !===================================================================!
  ! A build of one image (-fcoarray=single) started by an MPI
  ! launcher with several processes is invalid: every process would
  ! solve the whole problem and print as image 1. The launcher's
  ! process count is read from the environment (OpenMPI, MPICH).
  !===================================================================!

  subroutine require_images()

    character(len=32) :: value
    integer :: status, ranks

    if (num_images() > 1) return
    ranks = 1
    call get_environment_variable('OMPI_COMM_WORLD_SIZE', value, status=status)
    if (status /= 0) call get_environment_variable('PMI_SIZE', value, status=status)
    if (status == 0) read(value, *, iostat=status) ranks
    if (ranks > 1) then
       error stop 'operation_residual: this program was built for one image, but the launcher &
            &started several processes; run it without mpirun, or build it with COARRAY=lib'
    end if

  end subroutine require_images

  !===================================================================!
  ! THE OWNER OF EVERY UNKNOWN of a block: the cells are partitioned
  ! over the images by the partitioner's breadth-first rule on the
  ! mesh, each part connected, and every unknown at a cell, at every
  ! moment, is owned by the cell's image.
  !===================================================================!

  !===================================================================!
  ! THE OWNER OF EVERY CELL: image 1 of everything on one image; over
  ! several, the partitioner's breadth-first rule on the mesh, each
  ! part connected. The points of a block, moment outer and cell
  ! inner, are owned by their cell's image, and every loop of the
  ! sensitivity work over the points - the sources of the expansions,
  ! the adjoint's, the design products, the functional over its nodes
  ! - runs over the owned points alone and sums its result over the
  ! images once.
  !===================================================================!

  function cell_owners(this) result(cell_owner)

    class(discrete_residual), intent(in) :: this
    integer, allocatable :: cell_owner(:)

    type(partitioner) :: cut
    class(directed_graph), allocatable :: part
    type(partition_relation) :: relation
    integer :: ncells, k, v

    ncells = this % points % num_cells()
    ! without a region nothing is shared: every image owns every cell
    ! of its own, and no sum over the images is taken
    allocate(cell_owner(ncells), source=this_image())
    if (num_images() == 1 .or. .not. this % points % with_space) return
    cell_owner = 0
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

  end function cell_owners

  !===================================================================!
  ! THE IMAGE OF AN EXPANSION. On a manifold with a region the images
  ! share the cells, and every expansion runs on every image over its
  ! owned points. On a manifold of time alone there are no cells to
  ! share, and the expansions are shared instead: expansion d runs on
  ! image mod(d - 1, images) + 1 alone, its jets and its multipliers'
  ! derivatives summed over the images once at the end, the other
  ! images having left them zero. The expansions never read each
  ! other, so nothing else is exchanged.
  !===================================================================!

  logical function expansion_owned(this, d)

    class(discrete_residual), intent(in) :: this
    integer                 , intent(in) :: d

    expansion_owned = .true.
    if (this % points % with_space .or. num_images() == 1) return
    expansion_owned = mod(d - 1, num_images()) + 1 == this_image()

  end function expansion_owned

  ! the owner of every point of a block, moment outer and cell inner
  function point_owners(this, npts) result(owner)

    class(discrete_residual), intent(in) :: this
    integer                 , intent(in) :: npts
    integer, allocatable :: owner(:)

    integer, allocatable :: cell_owner(:)
    integer :: ncells, p, m

    cell_owner = cell_owners(this)
    ncells     = this % points % num_cells()
    allocate(owner(npts))
    do p = 1, npts
       m = (p - 1) / ncells + 1
       owner(p) = cell_owner(p - (m - 1) * ncells)
    end do

  end function point_owners

  function owners_of_unknowns(this, unknowns, stride, npts) result(owner)

    class(discrete_residual), intent(in) :: this
    integer                 , intent(in) :: unknowns, stride, npts
    integer, allocatable :: owner(:)

    integer, allocatable :: cell_owner(:)
    integer :: ncells, p, c, m, at, width

    ncells     = this % points % num_cells()
    cell_owner = cell_owners(this)

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
  ! from the tuples already solved.
  !===================================================================!

  subroutine fixed_rows_of(this, first, last, depth, instant_moment, at, tuple, fixed_rows, fixed)

    class(discrete_residual), intent(in)  :: this
    integer                 , intent(in)  :: first, last, depth, instant_moment(:), at(:)
    real(dp)                , intent(in)  :: tuple(:,:,:)
    integer , allocatable   , intent(out) :: fixed_rows(:)
    real(dp), allocatable   , intent(out) :: fixed(:)

    integer :: stride, ncells, ninst, count, k, c, i, j, f, m, e, instant, comp, degree
    real(dp) :: given, slope

    stride = this % rule % num_components()
    ncells = this % points % num_cells()
    ninst  = this % points % num_instants()

    count = depth * ncells * stride
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

  end subroutine fixed_rows_of

  !===================================================================!
  ! THE FACE ROWS of a block: at the first or the last instant of the
  ! manifold when it is one of the block's own, each face condition
  ! at every cell as a prescribed row - the affine relation of the
  ! components it reads, in place of the row of the component it
  ! determines, the one with coefficient one - with the relation's
  ! value at zero state as the row's constant. A face at a history
  ! instant is fixed by the history already. Invalid input: a face
  ! condition with no component of coefficient one, or several.
  !===================================================================!

  subroutine face_rows_of(this, first, last, depth, instant_moment, at, rows, columns, weights, constants)

    class(discrete_residual), intent(in)  :: this
    integer                 , intent(in)  :: first, last, depth, instant_moment(:), at(:)
    integer , allocatable   , intent(out) :: rows(:), columns(:)
    real(dp), allocatable   , intent(out) :: weights(:), constants(:)

    integer , allocatable :: fields(:), degrees(:)
    real(dp), allocatable :: slopes(:)
    real(dp) :: given
    integer :: ncells, count, j, i, m, c, e, k, instant, base, row_field, row_degree

    ncells = this % points % num_cells()
    count  = 0
    do j = 1, size(this % condition)
       if (.not. this % condition(j) % manifold % is_face()) cycle
       instant = face_instant(this, this % condition(j) % manifold % face_time)
       if (instant < first + depth .or. instant > last) cycle
       do i = 1, size(this % condition(j) % equation)
          call face_relation(this % condition(j) % equation(i), this % points % design_values(), &
               & this % points % position(this % points % point_of(instant, 1)), fields, degrees, slopes, given, &
               & row_field, row_degree)
          count = count + size(fields) * ncells
       end do
    end do
    allocate(rows(count), columns(count), weights(count), constants(count))
    e = 0
    do j = 1, size(this % condition)
       if (.not. this % condition(j) % manifold % is_face()) cycle
       instant = face_instant(this, this % condition(j) % manifold % face_time)
       if (instant < first + depth .or. instant > last) cycle
       m = face_row_moment(this, instant, first, last, depth, instant_moment)
       do i = 1, size(this % condition(j) % equation)
          do c = 1, ncells
             call face_relation(this % condition(j) % equation(i), this % points % design_values(), &
                  & this % points % position(this % points % point_of(instant, c)), fields, degrees, slopes, given, &
                  & row_field, row_degree)
             ! the row at the moment the condition occupies, the columns at the face's
             base = at((instant_moment(instant - first + 1) - 1) * ncells + c)
             do k = 1, size(fields)
                e = e + 1
                rows(e)      = at((m - 1) * ncells + c) + this % rule % offset_of_field(row_field) + row_degree + 1
                columns(e)   = base + this % rule % offset_of_field(fields(k)) + degrees(k) + 1
                weights(e)   = slopes(k)
                constants(e) = -given
             end do
          end do
       end do
    end do

  end subroutine face_rows_of

  !===================================================================!
  ! THE MOMENT WHOSE ROWS A FACE CONDITION OCCUPIES. At the block's
  ! first instant the components of a field are determined by nothing
  ! - no family reads back past the first instant - so a condition at
  ! the manifold's first instant occupies the rows there. A condition
  ! at the manifold's last instant, whose components the family's
  ! relations determine from the instants before, occupies the free
  ! rows at the block's first instant instead: the two-point problem,
  ! the rows at one end and the columns at the other. It therefore
  ! requires a block without history whose instants reach the last
  ! instant; a chain of several blocks cannot satisfy it, and is
  ! refused.
  !===================================================================!

  integer function face_row_moment(this, instant, first, last, depth, instant_moment) result(m)

    class(discrete_residual), intent(in) :: this
    integer                 , intent(in) :: instant, first, last, depth, instant_moment(:)

    if (instant == 1) then
       m = instant_moment(instant - first + 1)
    else if (instant == this % points % num_instants()) then
       if (depth > 0 .or. first > 1 .or. last < instant) then
          error stop 'operation_residual: a condition at the last instant is satisfied by one block over every &
               &instant, its rows at the first instant; a chain of several blocks marches past it'
       end if
       m = instant_moment(1)
    else
       error stop 'operation_residual: a face is at the first or the last instant'
    end if

  end function face_row_moment

  !===================================================================!
  ! A face condition as an affine relation: the parent fields and the
  ! degrees along time of the components it reads, the coefficient of
  ! each, given = -g at zero state so that the relation is sum of
  ! slopes times components = given, and the field and degree of the
  ! component it determines: of the last field declared among those
  ! read with coefficient +1 or -1 - the adjoint's, declared after
  ! the state, at the far face - the component of highest order, the
  ! relation negated when its coefficient is -1. Invalid input: no
  ! component of coefficient +1 or -1.
  !===================================================================!

  subroutine face_relation(equation, values, position, fields, degrees, slopes, given, row_field, row_degree)

    type(continuous_field), intent(in)  :: equation
    real(dp)              , intent(in)  :: values(:), position(:)
    integer , allocatable , intent(out) :: fields(:), degrees(:)
    real(dp), allocatable , intent(out) :: slopes(:)
    real(dp)              , intent(out) :: given
    integer               , intent(out) :: row_field, row_degree

    type(expression) :: g
    integer , allocatable :: reads(:)
    real(dp), allocatable :: q(:), grad(:)
    real(dp) :: value, sign
    integer :: k, f, i, comp, units

    g = equation % graph(1)
    call g % read_components(reads)
    if (size(reads) < 1) then
       error stop 'operation_residual: a condition on a face reads a component of the parent unknowns'
    end if
    allocate(q(0:g % num_components() - 1), grad(0:g % num_components() - 1), source=0.0_dp)
    call g % gradient_at(q, values, value, grad, position)
    given = -value
    allocate(fields(size(reads)), degrees(size(reads)), slopes(size(reads)))
    units      = 0
    row_field  = 0
    row_degree = -1
    do k = 1, size(reads)
       comp = reads(k)
       f = 0
       do i = 1, g % num_fields() - g % num_multipliers()
          if (comp >= g % offset_of_field(i) .and. comp <= g % offset_of_field(i) + g % degree_of_field(i)) f = i
       end do
       if (f == 0) then
          error stop 'operation_residual: a condition on a face reads the values of the unknowns and their &
               &derivatives along time'
       end if
       fields(k)  = f
       degrees(k) = comp - g % offset_of_field(f)
       slopes(k)  = grad(comp)
       if (abs(abs(slopes(k)) - 1.0_dp) <= 1.0e-12_dp) then
          if (f > row_field .or. (f == row_field .and. degrees(k) > row_degree)) then
             row_field  = f
             row_degree = degrees(k)
             sign       = slopes(k)
             units      = 1
          end if
       end if
    end do
    if (row_field == 0) then
       error stop 'operation_residual: a condition on a face determines one component, read with coefficient &
            &+1 or -1: of the last field declared among those so read, the one of highest order'
    end if
    slopes = slopes / sign
    given  = given / sign

  end subroutine face_relation

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
  ! The one parent component a face equation reads: the field f, the
  ! component in the equation's own tuple, and the degree of the
  ! derivative along time it is - zero for the value of f, one for
  ! its first derivative - since the equation's tuple and the rule's
  ! lay the jets out differently. Invalid input: an equation reading
  ! no component, several, or one outside the jets of the unknowns.
  !===================================================================!

  subroutine stated_value(equation, f, comp, degree)
    type(continuous_field), intent(in)  :: equation
    integer               , intent(out) :: f, comp, degree
    type(expression) :: g
    integer, allocatable :: reads(:)
    integer :: k
    g = equation % graph(1)
    call g % read_components(reads)
    if (size(reads) /= 1) then
       error stop 'operation_residual: a condition on a face states the value of one component'
    end if
    comp = reads(1)
    f = 0
    do k = 1, g % num_fields() - g % num_multipliers()
       if (comp >= g % offset_of_field(k) .and. comp <= g % offset_of_field(k) + g % degree_of_field(k)) f = k
    end do
    if (f == 0) then
       error stop 'operation_residual: a condition on a face states the value of an unknown or of one of &
            &its derivatives along time'
    end if
    degree = comp - g % offset_of_field(f)
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
    at_zero_value = mixed_partial(g % at_instant(q, spread(nu, 1, g % num_designs()), position))
    q(comp) = derivative_terms(1.0_dp, 0)
    at_one_value  = mixed_partial(g % at_instant(q, spread(nu, 1, g % num_designs()), position))
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
