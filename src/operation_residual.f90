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
  use operation_expression , only : expression, euler_lagrange, multiplier, constant, jet_form, represented_component
  use operation_expression , only : flux_form, coefficient_form, total_derivative, on_domain, paired
  use operation_expression , only : operator(+), operator(*)
  use operation_expression , only : expression_vertex, VERTEX_LEAF, ARGUMENT_STATE, FIRST_COORDINATE
  use operation_domain     , only : continuous_domain, discrete_domain
  use util_derivative_terms, only : max_subset_width, derivative_terms, mixed_partial, coefficient, symmetric_terms
  use util_derivative_terms, only : value
  use util_derivative_terms, only : operator(*)
  use token_identity       , only : token
  use operation_field      , only : continuous_field, discrete_field, WHOLE, TIME_FACE, TIME_FACTOR, SPACE_FACTOR
  use view_expression      , only : expression_view
  use operation_manifold   , only : continuous_manifold, discrete_manifold
  use operation_expansion  , only : expansion
  use operation_family     , only : family, chain
  use operation_weight     , only : scheme_weight
  use operation_coupling   , only : weights_of
  use view_directed_connectivity, only : connectivity_graph
  use operation_scheme_stencil  , only : derived_constraints
  use operation_derivative_approximation, only : derivative_approximation, element_quadrature
  use operation_segment_trapezoid, only : segment_trapezoid
  use operation_exchange       , only : exchange, state_summed_rows
  use util_verbosity           , only : verbosity
  use transform_partitioner    , only : partitioner, PARTITION_BREADTH_FIRST
  use operation_minimization, only : solve_result
  use operation_newton      , only : newton
  use operation_dense_direct, only : dense_direct
  use operation_elimination , only : elimination
  use operation_gmres       , only : gmres
  use operation_gauss_seidel, only : gauss_seidel

  implicit none

  private
  public :: residual_operator
  public :: continuous_residual, discrete_residual, minimize_together
  public :: operator(+), operator(*)

  !===================================================================!
  ! THE PROJECTION OF A LAW ON THE CELLS, tested against a basis in
  ! place of placed at the points. An instance is an element at one
  ! moment: element_point(i, n) is the point of node i of instance n,
  ! negative where another part states the node's rows, and
  ! element_cell(n) the element of the quadrature. A quadrature point
  ! of an instance is a site, where the rules are evaluated: site
  ! (n - 1) nq + q is the point q of instance n. The tuple the rules
  ! read at a site is formed from the tuples of the instance's nodes:
  ! its component c is the sum over the nodes and their coefficients
  ! of the basis's derivative formed(0, c) times the node's component
  ! formed(k + 1, c), -1 for a coefficient that forms none. Rule j
  ! tests the equation tested(1, j) with the basis's derivative
  ! tested(2, j), and coefficient_row(k, f) is the row, within a node's
  ! tuple, of the coefficient k of the equation f: the rule's value at
  ! a site, times the weight and the test function there, is added to
  ! the row of every coefficient of every node the instance states.
  !===================================================================!

  type :: projection
     type(element_quadrature) :: elements
     integer , allocatable :: element_point(:,:)
     integer , allocatable :: element_cell(:)
     integer , allocatable :: tested(:,:)
     integer , allocatable :: coefficient_row(:,:)
     integer , allocatable :: formed(:,:)
  end type projection

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
     ! the position of every evaluation point, one column per point:
     ! the time, then the space coordinates; no rows when none is given
     real(dp), allocatable, private :: point_positions(:,:)
     ! the rows that read every point, the integral of an unknown over
     ! the region: over several images each part states the terms of
     ! its owned columns and summed_by adds the others'
     integer       , allocatable, private :: summed_rows(:)
     type(exchange), allocatable, private :: summed_by
     ! the law projected on the cells, in place of placed at the
     ! points; its sites as the vertices the rules are evaluated over,
     ! and the position of each
     type(projection), allocatable, private :: projected
     type(stored_directed_graph)  , private :: sites
     real(dp), allocatable        , private :: site_positions(:,:)
     ! the tying rows, each of one unknown on its own diagonal reading
     ! others before it: the history, the family's relations, and the
     ! copies and zero coefficients of the unknowns of time alone,
     ! eliminated before a direct solve of the rest
     logical, allocatable, private :: tying(:)
     ! THE RULES AT THE EXPLICIT POINTS, those of the moments at which
     ! the family's quadrature of the moments before determines every
     ! component below the top degree: there a law that reads none - a
     ! joint's condition - is stated by its derivative along time, as
     ! at the manifold's first instant; explicit(p) marks the points
     type(expression), allocatable, private :: explicit_physics
     type(expression), allocatable, private :: explicit_rules(:)
     logical         , allocatable, private :: explicit(:)

   contains

     procedure :: name           => residual_name
     procedure :: domain         => residual_domain
     procedure :: apply          => residual_apply
     procedure :: max_degree     => residual_max_degree
     procedure :: defined_at_zero => residual_defined_at_zero
     procedure :: partial_action => residual_partial_action
     procedure :: explicit_jacobian => residual_explicit_jacobian
     procedure, private :: positions => residual_positions

     procedure :: num_unknowns
     procedure :: fixed_unknowns
     procedure :: fixed_indicator
     procedure :: fixed_values
     procedure :: fixed_rates
     procedure :: num_fixed
     procedure :: first_fixed
     procedure :: stride
     procedure :: coupling
     procedure :: part_of
     procedure, private :: terms
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

     procedure, private :: any_fixed_rate
     procedure, private :: discretized
     procedure, private :: gathered
     procedure, private :: governed
     procedure, private :: governed_on_reals
     procedure, private :: physics_tangent
     procedure, private :: physics_variation
     procedure, private :: placed
     procedure, private :: placed_output
     procedure, private :: point_inputs
     procedure, private :: point_of
     procedure, private :: require_direction
     procedure, private :: require_host
     procedure, private :: state_of
     procedure, private :: write_prescribed
     procedure, private :: zero_prescribed_rows
     procedure, private :: formed_at_sites
     procedure, private :: design_at_sites
     procedure, private :: tested_into
     procedure, private :: projected_count
     procedure, private :: projected_columns
     procedure, private :: projected_triples
     procedure, private :: projected_edges
     procedure, private :: projected_part
     procedure, private :: explicit_marks
     procedure, private :: explicit_at
  end type residual_operator

  interface residual_operator
     module procedure create
  end interface residual_operator

  !===================================================================!
  ! THE CONTINUOUS RESIDUAL: a condition, a chain of conditions, a
  ! Lagrangian - one type. A condition is a list of equations, each a
  ! scalar field, on one manifold; paired with a multiplier, an
  ! unknown of that manifold with one component per equation,
  ! declared by the condition itself, it contributes to the Lagrangian
  ! the product of the two. The sum of two residuals is the residual
  ! of their conditions in the order summed, so that a Lagrangian is a
  ! chain of residuals and needs no type of its own. The equations of
  ! a condition on a part of a manifold - a face at an instant, the
  ! time factor - read the unknowns of the parent. The objective J, a
  ! scalar integral over the manifold, enters the Lagrangian as a
  ! term of its own, L = J + lambda_r . r + lambda . g: its
  ! stationarity in the state is the adjoint equation, and the
  ! multipliers are the sensitivities of J to their equalities.
  !
  ! A CONDITION NEEDS NO MULTIPLIER. The reverse sweep, which runs
  ! when the equations are paired, solves the transposed system for
  ! the multiplier of every row, whether or not its condition was
  ! paired: the multipliers mirror the conditions. Pairing a condition
  ! names its multiplier, and a named multiplier is returned with the
  ! solution; a condition summed bare, L = J + lambda_r . r + g, is
  ! satisfied the same way and its multiplier is not returned.
  !===================================================================!

  type :: residual_term

     type(continuous_manifold) :: manifold
     type(continuous_field), allocatable :: equation(:)
     type(continuous_field) :: multiplier
     logical :: paired = .false.
     ! an inequality on the design factor: each equation at most zero,
     ! its multiplier nonnegative and zero where the equation is below
     ! zero
     logical :: inequality = .false.

  end type residual_term

  ! the first instant, stage instant and node that the arrays of a
  ! sweep store: one for arrays over the chain, later for arrays over
  ! a window of it
  type :: lower_bounds
     integer :: instant = 1
     integer :: stage   = 1
     integer :: node    = 1
  end type lower_bounds

  type(lower_bounds), parameter :: origin = lower_bounds(1, 1, 1)

  !===================================================================!
  ! ONE BLOCK OF THE CHAIN: the instants first..last, the first
  ! `depth` of them the history, solved by the blocks before, and the
  ! rows stated on its moments - each instant and, for a family
  ! marching by stages, the stages of the step into it - with the
  ! unknowns q of those moments, the estimate and then the solution.
  ! The unknowns list the moments in the order of time, each moment's
  ! cells in order, each cell's tuple: instant_moment(k) is the moment
  ! of the k-th instant, and the stages of the step into the instant
  ! k + 1 are the moments after it. The point manifold is the block of
  ! one moment, without a history. Where the block is linearized,
  ! `inputs` are the rows' inputs at q, by which their jacobian is
  ! formed there, and `owner` the image of every point.
  !===================================================================!

  type :: chain_block

     type(residual_operator) :: rows
     real(dp), allocatable :: q(:)
     integer , allocatable :: instant_moment(:)
     integer :: first = 1, last = 1, depth = 0
     integer :: num_stages = 0
     integer :: num_cells  = 1
     type(stored_field), allocatable :: inputs(:)
     integer           , allocatable :: owner(:)

   contains

     procedure :: tuple_at
     procedure :: from_points
     procedure :: to_points

  end type chain_block

  ! the moments of a block at which a column of the manifold is read
  ! or written: the history instants; the block's own instants and the
  ! stages of its steps; every moment
  integer, parameter :: HISTORY_MOMENTS = 1, OWN_MOMENTS = 2, EVERY_MOMENT = 3

  ! the actions of a scheduled reverse sweep
  integer, parameter :: ACTION_RESTORE = 1, ACTION_ADVANCE = 2, ACTION_STORE = 3, ACTION_REVERSE = 4

  ! the iterations of Newton's method on a free design
  integer, parameter :: MAX_DESIGN_ITERATIONS = 50

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
     ! the representation laws paired with a multiplier, whose
     ! multipliers the reverse sweep returns under their names
     type(residual_term), allocatable :: representation(:)
     type(continuous_field), allocatable :: objective
     logical                   :: with_adjoint = .false.
     ! two images divide a checkpointed reverse sweep: image 2 solves
     ! the segments again, image 1 forms their multipliers; every
     ! image then owns every expansion
     logical                   :: overlapped = .false.
     character(len=32)         :: adjoint_name = ''
     ! THE LAW PROJECTED ON THE CELLS, for a method that tests it
     ! against a basis: its flux form as a Lagrangian, each multiplier
     ! of which tests the equation tested(1, m) with the basis's
     ! derivative tested(2, m), and the elements of the region; the
     ! rule is then the layout of the coefficients a node stores
     type(expression)        , allocatable :: flux
     integer                 , allocatable :: tested(:,:)
     type(element_quadrature), allocatable :: elements
     ! THE RULE AT THE MANIFOLD'S FIRST INSTANT, where a law the initial
     ! data satisfy is stated by its derivative along time; unallocated
     ! when no law is. Projected, the flux form there. The same rule
     ! states the laws at the explicit moments of a block, the instants
     ! the steps of a staged family reach
     type(expression)        , allocatable :: initial_rule
     type(expression)        , allocatable :: initial_flux
     ! each block after the first starts from the solution at the last
     ! instant before it (minimize, continued)
     logical                 :: continued = .false.
     ! THE MULTIPLIERS OF THE INEQUALITIES ON FUNCTIONALS of the state,
     ! one per inequality in the order stated, zero for one on the
     ! design alone: the reverse sweep adds each functional's
     ! differential times its multiplier to the objective's, the
     ! Lagrangian's, so that the multipliers of the design and their
     ! expansion are its gradient and its hessian in the design
     real(dp), allocatable   :: functional_weight(:)

   contains

     procedure :: minimize
     procedure :: disagreement
     procedure, private :: zero_at_design
     procedure, private :: inequalities_at
     procedure, private :: functional_inequalities
     procedure, private :: extrapolated
     procedure, private :: place_inequality_multipliers
     procedure, private :: check_representations
     procedure, private :: num_represented
     procedure, private :: block_positions

     procedure, private :: adjoint_block
     procedure, private :: adjoint_order
     procedure, private :: adjoint_rows
     procedure, private :: adjoint_sweep
     procedure, private :: allocate_multipliers
     procedure, private :: block_extent
     procedure, private :: block_inputs
     procedure, private :: block_rows
     procedure, private :: cell_owners
     procedure, private :: expansion_order
     procedure, private :: expansion_owned
     procedure, private :: face_instant
     procedure, private :: face_row_moment
     procedure, private :: face_rows_of
     procedure, private :: space_face_relation
     procedure, private :: face_relation
     procedure, private :: condition_relation
     procedure, private :: space_face_placement
     procedure, private :: segment_rows_of
     procedure, private :: open_line
     procedure, private :: whole_line
     procedure, private :: fixed_rows_of
     procedure, private :: gauge_rows_of
     procedure, private :: time_alone_rows_of
     procedure, private :: estimated_rates
     procedure, private :: gauged_unknowns
     procedure, private :: multipliers_of_order
     procedure, private :: owners_of_unknowns
     procedure, private :: pipelined
     procedure, private :: point_owners
     procedure, private :: point_rows
     procedure, private :: quadrature_nodes
     procedure, private :: scheduled_sweeps
     procedure, private :: shared_expansions_summed
     procedure, private :: solve_block
     procedure, private :: solve_rows
     procedure, private :: solved
     procedure, private :: solver_for
     procedure, private :: wavefront_adjoint
     procedure, private :: wavefront_expansion
     procedure, private :: coefficient_component
     procedure, private :: formed_components
     procedure, private :: points_on_face
     procedure, private :: stored_order
     procedure, private :: projected_nodes
     procedure, private :: site_jet
     procedure, private :: site_gradient
     procedure, private :: coefficient_derivative
     procedure, private :: fixed_coefficient
     procedure, private :: projected_datum
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
       & connected_law, fixed_rate, prescribed, positions, summed_rows, projected, tying, explicit_rule, explicit) &
       & result(this)

    type(stencil)         , intent(in) :: primary_law
    type(expression)      , intent(in), optional :: rule
    integer               , intent(in) :: at(:), unknowns, degrees, primary(:)
    integer               , intent(in) :: fixed_rows(:)
    real(dp)              , intent(in) :: fixed(:)
    type(stencil)         , intent(in), optional :: connected_law
    real(dp)              , intent(in), optional :: fixed_rate(:)
    type(stencil)         , intent(in), optional :: prescribed
    real(dp)              , intent(in), optional :: positions(:,:)
    integer               , intent(in), optional :: summed_rows(:)
    type(projection)      , intent(in), optional :: projected
    logical               , intent(in), optional :: tying(:)
    type(expression)      , intent(in), optional :: explicit_rule
    logical               , intent(in), optional :: explicit(:)
    type(residual_operator) :: this
    type(continuous_domain) :: domain
    integer :: j, components, governed, e, n, q, nq, s, dimension
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
       ! a projected law reads at its sites the tuple the basis forms,
       ! of its own layout, from the nodes' tuples of degrees
       ! components
       if (present(projected)) components = degrees
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
    ! the rules at the explicit points, of the rule's multipliers
    if (present(explicit_rule) .and. present(explicit)) then
       if (size(explicit) /= size(at)) then
          error stop 'operation_residual: one explicit mark is required per evaluation point'
       end if
       if (any(explicit)) then
          if (max(1, explicit_rule % num_multipliers()) /= size(this % rules)) then
             error stop 'operation_residual: the rule at the explicit points governs the rows the rule governs'
          end if
          this % explicit_physics = explicit_rule
          if (explicit_rule % num_multipliers() > 0) then
             allocate(this % explicit_rules(explicit_rule % num_multipliers()))
             do j = 1, explicit_rule % num_multipliers()
                this % explicit_rules(j) = euler_lagrange(explicit_rule, j)
             end do
          else
             this % explicit_rules = [explicit_rule]
          end if
          this % explicit = explicit
       end if
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
    ! the positions of the points, one column per point
    if (present(positions)) then
       if (size(positions, 2) /= size(at)) then
          write(message,'(a,i0,a,i0)') 'operation_residual: one position is required per evaluation &
               &point; size(positions, 2) = ', size(positions, 2), ', points = ', size(at)
          error stop trim(message)
       end if
       this % point_positions = positions
    else
       allocate(this % point_positions(0, size(at)))
    end if
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
    if (present(summed_rows)) this % summed_rows = summed_rows
    if (present(tying)) this % tying = tying
    ! the sites of a projection, at the time of their instance's moment
    ! and the position of their quadrature point
    if (present(projected)) then
       if (.not. present(positions)) then
          error stop 'operation_residual: a law projected on the cells is evaluated at the positions of its sites'
       end if
       this % projected = projected
       nq        = projected % elements % num_quadrature_points()
       dimension = projected % elements % dimension
       this % sites = stored_directed_graph(size(projected % element_cell) * nq, tails=[integer ::], heads=[integer ::])
       allocate(this % site_positions(size(positions, 1), size(projected % element_cell) * nq))
       do n = 1, size(projected % element_cell)
          do q = 1, nq
             s = (n - 1) * nq + q
             this % site_positions(:, s) = positions(:, abs(projected % element_point(1, n)))
             this % site_positions(size(positions, 1) - dimension + 1:, s) = &
                  & projected % elements % position(:, q, projected % element_cell(n))
          end do
       end do
    end if
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

  !===================================================================!
  ! The coupling of the unknowns: an edge from every column to every
  ! row of the stencils' entries. The physics placed at the points
  ! reads a point's own tuple alone, so it couples no two points; a
  ! projected law couples every node of an instance with every other.
  !===================================================================!

  function coupling(this) result(graph_of)
    class(residual_operator), intent(in) :: this
    type(stored_directed_graph) :: graph_of
    integer , allocatable :: heads(:), tails(:)
    real(dp), allocatable :: w(:)
    logical, allocatable :: is_summed(:)
    call this % terms(heads, tails, w)
    ! a projected law couples the nodes of each instance
    if (allocated(this % projected)) call this % projected_edges(heads, tails)
    ! a summed row couples no point to another: its sum is formed over
    ! the images
    allocate(is_summed(this % unknowns), source=.false.)
    if (allocated(this % summed_rows)) is_summed(this % summed_rows) = .true.
    graph_of = stored_directed_graph(this % unknowns, tails=pack(tails, .not. is_summed(heads)), &
         & heads=pack(heads, .not. is_summed(heads)))
  end function coupling

  ! every stencil's terms, as (row, column, weight)
  subroutine terms(this, heads, tails, w)
    class(residual_operator), intent(in)  :: this
    integer , allocatable   , intent(out) :: heads(:), tails(:)
    real(dp), allocatable   , intent(out) :: w(:)
    integer , allocatable :: h(:), t(:)
    real(dp), allocatable :: v(:)
    call this % primary_law % entries(heads, tails, w)
    if (allocated(this % connected_law)) then
       call this % connected_law % entries(h, t, v)
       heads = [heads, h]
       tails = [tails, t]
       w     = [w, v]
    end if
    if (allocated(this % prescribed)) then
       call this % prescribed % entries(h, t, v)
       heads = [heads, h]
       tails = [tails, t]
       w     = [w, v]
    end if
  end subroutine terms

  !===================================================================!
  ! THE PART of the residual an image forms: the rows of the points the
  ! part owns, over the unknowns of the part in its own numbering, the
  ! owned first and the halo after. The physics is stated at the owned
  ! points, each stencil restricted to the owned rows, and the fixed
  ! rows the part owns. Invalid input: the points' tuples not in the
  ! order of the unknowns.
  !===================================================================!

  function part_of(this, part) result(sub)
    class(residual_operator), intent(in)    :: this
    ! the part, given the terms of the summed rows it adds
    type(exchange)          , intent(inout) :: part
    type(residual_operator) :: sub
    type(stencil), allocatable :: connected, prescribed
    integer , allocatable :: points(:), at(:), place_of(:), fixed_rows(:), summed(:), heads(:), tails(:)
    integer , allocatable :: point_place(:)
    real(dp), allocatable :: fixed(:), rate(:), weights(:)
    logical , allocatable :: owned_fixed(:), owned_point(:)
    integer :: w, npts, i
    w = this % stride()
    if (any(this % at /= [((i - 1) * w, i = 1, size(this % at))])) then
       error stop 'operation_residual: a residual is cut into parts when its points'' tuples follow &
            &one another in the order of the unknowns'
    end if
    npts   = size(part % owned) / w
    ! a projected law reads the halo's tuples too: its points are the
    ! part's every point, the owned first
    if (allocated(this % projected)) npts = size(part % whole_index) / w
    points = [((part % whole_index((i - 1) * w + 1) - 1) / w + 1, i = 1, npts)]
    at     = [((i - 1) * w, i = 1, npts)]
    allocate(place_of(this % unknowns), source=0)
    place_of(part % whole_index) = [(i, i = 1, size(part % whole_index))]
    owned_fixed = this % fixed_rows > 0
    do i = 1, size(this % fixed_rows)
       owned_fixed(i) = place_of(this % fixed_rows(i)) > 0
       if (owned_fixed(i)) owned_fixed(i) = part % is_owned(place_of(this % fixed_rows(i)))
    end do
    fixed_rows = place_of(pack(this % fixed_rows, owned_fixed))
    fixed      = pack(this % fixed, owned_fixed)
    rate       = pack(this % fixed_rate, owned_fixed)
    allocate(summed(0))
    if (allocated(this % summed_rows)) summed = this % summed_rows
    call this % terms(heads, tails, weights)
    call state_summed_rows(part, summed, heads, tails, weights)
    if (allocated(this % connected_law)) connected = part_stencil(this % connected_law, part, summed)
    if (allocated(this % prescribed))    prescribed = part_stencil(this % prescribed, part, summed)
    if (allocated(this % projected)) then
       allocate(point_place(size(this % at)), source=0)
       point_place(points) = [(i, i = 1, npts)]
       owned_point = [(part % is_owned((i - 1) * w + 1), i = 1, npts)]
       sub = residual_operator(part_stencil(this % primary_law, part, summed), this % physics, at, &
            & size(part % whole_index), this % degrees, this % primary, fixed_rows, fixed, &
            & connected_law=connected, fixed_rate=rate, prescribed=prescribed, &
            & positions=this % point_positions(:, points), &
            & projected=this % projected_part(point_place, owned_point))
    else if (allocated(this % physics)) then
       sub = residual_operator(part_stencil(this % primary_law, part, summed), this % physics, at, &
            & size(part % whole_index), this % degrees, this % primary, fixed_rows, fixed, &
            & connected_law=connected, fixed_rate=rate, prescribed=prescribed, &
            & positions=this % point_positions(:, points))
    else
       sub = residual_operator(part_stencil(this % primary_law, part, summed), at=at, &
            & unknowns=size(part % whole_index), degrees=this % degrees, primary=this % primary, &
            & fixed_rows=fixed_rows, fixed=fixed, connected_law=connected, fixed_rate=rate, &
            & prescribed=prescribed, positions=this % point_positions(:, points))
    end if
    sub % summed_by = part
    call sub % explicit_at(this, points)
  end function part_of

  ! a stencil's rows the part owns, in the part's numbering, with their
  ! constants
  function part_stencil(op, part, summed) result(sub)
    type(stencil) , intent(in) :: op
    type(exchange), intent(in) :: part
    integer       , intent(in) :: summed(:)
    type(stencil) :: sub
    integer , allocatable :: heads(:), tails(:), rows(:), columns(:)
    real(dp), allocatable :: w(:), weights(:), constants(:)
    call op % entries(heads, tails, w)
    call op % constants % real_vector(constants)
    call part % restricted_triples(heads, tails, w, rows, columns, weights, summed)
    sub = stencil(rows, columns, weights, merge(constants(part % whole_index), 0.0_dp, part % is_owned), &
         & label=op % name())
  end function part_stencil

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
    if (allocated(this % explicit_physics)) residual_max_degree = max(residual_max_degree, &
         & this % explicit_physics % max_degree())
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
    if (allocated(this % explicit_rules)) then
       do j = 1, size(this % explicit_rules)
          defined = defined .and. this % explicit_rules(j) % defined_at_zero()
       end do
    end if
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
    if (allocated(this % projected)) then
       domain  = continuous % discrete(this % sites)
       points  = domain % state_fields()
       designs = domain % design_fields()
       state   = points % state(this % formed_at_sites(x))
       design  = designs % design(this % design_at_sites(design_values))
    else
       domain  = continuous % discrete(this % points)
       points  = domain % state_fields()
       designs = domain % design_fields()
       state   = points % state(this % gathered(x))
       design  = designs % design(design_values)
    end if
    point_data = [state, design]
  end subroutine point_inputs

  pure subroutine placed(this, governing, r)
    class(residual_operator), intent(in)    :: this
    real(dp)                 , intent(in)    :: governing(:,:)
    real(dp)                 , intent(inout) :: r(:)
    integer :: p, j
    if (allocated(this % projected)) then
       call this % tested_into(governing, r)
       return
    end if
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
    real(dp), allocatable :: column(:)
    logical , allocatable :: marked(:)
    integer :: j
    ! at the points, or at the sites of a projection
    if (allocated(this % projected)) then
       allocate(governing(this % sites % num_vertices(), size(this % rules)))
    else
       allocate(governing(size(this % at), size(this % rules)))
    end if
    marked = this % explicit_marks()
    do j = 1, size(this % rules)
       call evaluated(this % rules(j), column)
       governing(:, j) = column
       ! at the explicit sites, the rule there
       if (allocated(this % explicit_rules)) then
          call evaluated(this % explicit_rules(j), column)
          where (marked) governing(:, j) = column
       end if
    end do
  contains
    subroutine evaluated(law, column)
      type(expression)     , intent(in)  :: law
      real(dp), allocatable, intent(out) :: column(:)
      class(field), allocatable :: half
      type(stored_field), allocatable :: bound(:)
      integer :: coordinates
      ! a rule that reads a coordinate is given the positions of the
      ! points as its third argument
      bound       = point_data
      coordinates = law % num_position_coordinates()
      if (coordinates > 0) bound = [point_data, this % positions(coordinates)]
      if (allocated(this % projected)) then
         if (present(variations)) then
            call law % partial_action(this % sites, law % bind(bound), variations, half)
         else
            call law % apply(this % sites, law % bind(bound), half)
         end if
      else if (present(variations)) then
         call law % partial_action(this % points, law % bind(bound), variations, half)
      else
         call law % apply(this % points, law % bind(bound), half)
      end if
      call half % real_vector(column)
    end subroutine evaluated
  end subroutine governed

  !===================================================================!
  ! The explicit marks at the evaluation sites: at the points, or at
  ! each site of a projection its instance's, the mark of the instance's
  ! first point; none without rules at the explicit points.
  !===================================================================!

  function explicit_marks(this) result(marked)
    class(residual_operator), intent(in) :: this
    logical, allocatable :: marked(:)
    integer :: nq, s
    if (.not. allocated(this % explicit)) then
       allocate(marked(0))
       return
    end if
    if (allocated(this % projected)) then
       nq = this % projected % elements % num_quadrature_points()
       allocate(marked(size(this % projected % element_cell) * nq))
       do s = 1, size(marked)
          marked(s) = this % explicit(abs(this % projected % element_point(1, (s - 1) / nq + 1)))
       end do
    else
       marked = this % explicit
    end if
  end function explicit_marks

  ! every rule's value at every site of the projection, on the real
  ! values the basis forms there
  subroutine governed_on_reals(this, inputs, x, governing)
    class(residual_operator), intent(in)  :: this
    type(binding)           , intent(in)  :: inputs(:)
    real(dp)                , intent(in)  :: x(:)
    real(dp), allocatable   , intent(out) :: governing(:,:)
    real(dp), allocatable :: design_values(:), xs(:), ns(:)
    logical , allocatable :: marked(:)
    integer :: j, s, width, nd, nsites, coordinates, nq
    call bound_real_vector(inputs, this % argument(2), design_values)
    xs     = this % formed_at_sites(x)
    ns     = this % design_at_sites(design_values)
    width  = size(this % projected % formed, 2)
    nd     = this % designs
    nsites = this % sites % num_vertices()
    allocate(governing(nsites, size(this % rules)))
    nq = this % projected % elements % num_quadrature_points()
    marked = this % explicit_marks()
    do j = 1, size(this % rules)
       do s = 1, nsites
          ! no rule is evaluated where its test functions vanish
          if (.not. tested_at(this % projected % elements, this % projected % tested(2, j), s - ((s - 1) / nq) * nq, &
               & this % projected % element_cell((s - 1) / nq + 1))) then
             governing(s, j) = 0.0_dp
             cycle
          end if
          ! at an explicit site, the rule there
          if (allocated(this % explicit_rules)) then
             if (marked(s)) then
                governing(s, j) = value_of(this % explicit_rules(j))
                cycle
             end if
          end if
          governing(s, j) = value_of(this % rules(j))
       end do
    end do
  contains
    real(dp) function value_of(law)
      type(expression), intent(in) :: law
      coordinates = law % num_position_coordinates()
      if (coordinates > 0) then
         value_of = law % value_on(xs((s - 1) * width + 1:s * width), ns((s - 1) * nd + 1:s * nd), &
              & this % site_positions(1:coordinates, s))
      else
         value_of = law % value_on(xs((s - 1) * width + 1:s * width), ns((s - 1) * nd + 1:s * nd))
      end if
    end function value_of
  end subroutine governed_on_reals

  !===================================================================!
  ! The first `coordinates` coordinates of every point's position as a
  ! field on the points, one component per coordinate. A rule that
  ! reads a coordinate the positions do not store stops the program.
  !===================================================================!

  function residual_positions(this, coordinates) result(positions)
    class(residual_operator), intent(in) :: this
    integer                 , intent(in) :: coordinates
    type(stored_field) :: positions
    character(len=250) :: message
    if (coordinates > size(this % point_positions, 1)) then
       write(message,'(a,i0,a,i0)') 'operation_residual: a rule reads a coordinate the positions of the &
            &points do not store; coordinates read = ', coordinates, ', coordinates stored = ', &
            & size(this % point_positions, 1)
       error stop trim(message)
    end if
    if (allocated(this % projected)) then
       positions = stored_field('positions', this % sites % vertex_set(), this % sites % num_vertices(), &
            & num_components=coordinates)
       call positions % set_real_vector(reshape(this % site_positions(1:coordinates, :), &
            & [coordinates * this % sites % num_vertices()]))
       return
    end if
    positions = stored_field('positions', this % points % vertex_set(), size(this % at), num_components=coordinates)
    call positions % set_real_vector(reshape(this % point_positions(1:coordinates, :), [coordinates * size(this % at)]))
  end function residual_positions

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
    call this % require_host(input_graph)
    call this % state_of(inputs, x, state)
    call this % discretized(inputs, state, x, r, governing)
    call this % placed(governing, r)
    r(this % fixed_rows) = x(this % fixed_rows) - this % fixed
    call this % write_prescribed(x, r, along=.false.)
    if (allocated(this % summed_by)) call this % summed_by % add_summed_rows(x, r)
    call this % placed_output(r, output)
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
    ! with the tying rows eliminated, each rule's value at each site on
    ! the real values, by the pass its jacobian takes
    if (allocated(this % projected) .and. allocated(this % tying) .and. .not. present(v)) then
       call this % governed_on_reals(inputs, x, governing)
       return
    end if
    call this % point_inputs(inputs, x, point_data)
    if (present(v)) then
       if (allocated(this % projected)) then
          points    = typed_field_domain(this % sites, this % physics % num_components())
          direction = points % direction(this % formed_at_sites(v))
       else
          points    = typed_field_domain(this % points, this % stride())
          direction = points % direction(this % gathered(v))
       end if
       varied(1) = variation(this % physics % argument(1), direction)
       call this % governed(point_data, governing, varied)
    else
       call this % governed(point_data, governing)
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
    integer , allocatable :: r(:), c(:), reads(:), explicit_reads(:)
    logical , allocatable :: is_fixed(:), marked(:)
    integer :: e, d, p, npts, n, num_triples, count, j, k, deg
    real(dp) :: value
    jacobian_defined = which == 1
    if (.not. jacobian_defined) return
    call this % require_host(input_graph)
    n    = this % unknowns
    npts = size(this % at)
    deg  = this % degrees
    is_fixed = this % fixed_indicator()
    call this % state_of(inputs, x, state)
    if (allocated(this % physics)) then
       if (.not. bound_on(inputs, this % argument(2), this % design_domain(), npts * this % designs)) then
          error stop 'operation_residual: the design must be defined on the point domain with one value per &
               &design per point'
       end if
       call bound_real_vector(inputs, this % argument(2), nu)
       xs = this % gathered(x)
       allocate(g(0:deg - 1))
    end if
    count = this % primary_law % pattern % num_edges() + npts * this % degrees * size(this % rules) &
         & + size(this % fixed_rows)
    if (allocated(this % projected)) count = this % primary_law % pattern % num_edges() + size(this % fixed_rows) &
         & + this % projected_count()
    if (allocated(this % connected_law)) count = count + this % connected_law % pattern % num_edges()
    if (allocated(this % prescribed))    count = count + this % prescribed % pattern % num_edges()
    allocate(r(count), c(count), w(count))
    num_triples = 0
    call stencil_triples(this % primary_law, is_fixed, r, c, w, num_triples)
    if (allocated(this % connected_law)) call stencil_triples(this % connected_law, is_fixed, r, c, w, num_triples)
    ! each rule's partials in the components it reads, by one reverse
    ! pass per point: a component no leaf of the rule names has a zero
    ! column; projected, by one pass per site, through the basis
    if (allocated(this % projected)) then
       call this % projected_triples(x, nu, is_fixed, r, c, w, num_triples)
    end if
    marked = this % explicit_marks()
    do j = 1, size(this % rules)
       if (allocated(this % projected)) exit
       call this % rules(j) % read_components(reads)
       if (allocated(this % explicit_rules)) call this % explicit_rules(j) % read_components(explicit_reads)
       do p = 1, npts
          if (is_fixed(this % at(p) + this % primary(j) + 1)) cycle
          ! at an explicit point, the rule there
          if (allocated(this % explicit_rules)) then
             if (marked(p)) then
                call this % explicit_rules(j) % gradient_at(xs((p - 1) * deg + 1:p * deg), &
                     & nu((p - 1) * this % designs + 1:p * this % designs), value, g, this % point_positions(:, p))
                call triples_at(explicit_reads)
                cycle
             end if
          end if
          call this % rules(j) % gradient_at(xs((p - 1) * deg + 1:p * deg), &
               & nu((p - 1) * this % designs + 1:p * this % designs), value, g, this % point_positions(:, p))
          call triples_at(reads)
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
  contains
    ! the partials of rule j at point p in the components it reads
    subroutine triples_at(read)
      integer, intent(in) :: read(:)
      do k = 1, size(read)
         d = read(k)
         num_triples    = num_triples + 1
         r(num_triples) = this % at(p) + this % primary(j) + 1
         c(num_triples) = this % at(p) + d + 1
         w(num_triples) = g(d)
      end do
    end subroutine triples_at
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
    call this % require_host(input_graph)
    do i = 1, size(variations)
       call this % require_direction(variations(i))
    end do
    call this % state_of(inputs, x, state)
    if (size(variations) >= 2) then
       ! A mixed or repeated partial of x(row) - h(nu) in the state is
       ! zero at every order above the first. A repeated partial in the
       ! design alone is -d^m h/dnu^m, of which only the first rate is
       ! stated: refuse rather than report zero for a value whose rate
       ! is not zero.
       if (this % any_fixed_rate()) then
          if (all([(variations(i) % argument_is(this % argument(2)), i = 1, size(variations))])) then
             error stop 'operation_residual: a fixed value''s design rate is stated to first order &
                  &only; a repeated design partial of a fixed row is not stated'
          end if
       end if
       call this % physics_tangent(inputs, variations, x, governing)
       allocate(r(this % num_unknowns()), source=0.0_dp)
       call this % placed(governing, r)
       r(this % fixed_rows) = 0.0_dp
       call this % zero_prescribed_rows(r)
       call this % placed_output(r, output)
       return
    end if
    call variations(1) % direction(v)
    if (variations(1) % argument_is(this % argument(1))) then
       call this % discretized(inputs, state, x, r, governing, v)
       call this % placed(governing, r)
       r(this % fixed_rows) = v(this % fixed_rows)
       call this % write_prescribed(v, r, along=.true.)
       if (allocated(this % summed_by)) call this % summed_by % add_summed_rows(v, r)
    else
       allocate(r(this % num_unknowns()), source=0.0_dp)
       call this % physics_tangent(inputs, variations, x, governing)
       call this % placed(governing, r)
       ! a fixed row states x(row) - h(row); its partial along a design
       ! direction v is -dh(row)/dnu times v at the point of the row
       r(this % fixed_rows) = -this % fixed_rate * v([(this % point_of(this % fixed_rows(i)), i = 1, size(this % fixed_rows))])
       call this % zero_prescribed_rows(r)
    end if
    call this % placed_output(r, output)
  end subroutine residual_partial_action

  subroutine physics_tangent(this, inputs, variations, x, governing)
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
    call this % point_inputs(inputs, x, point_data)
    allocate(at_points(size(variations)))
    do i = 1, size(variations)
       at_points(i) = this % physics_variation(variations(i))
    end do
    call this % governed(point_data, governing, at_points)
  end subroutine physics_tangent

  function physics_variation(this, given) result(at_points)
    class(residual_operator), intent(in) :: this
    type(variation)          , intent(in) :: given
    type(variation) :: at_points
    type(stored_field) :: direction
    type(typed_field_domain) :: points
    real(dp), allocatable :: v(:)
    call given % direction(v)
    if (given % argument_is(this % argument(1)) .and. allocated(this % projected)) then
       points    = typed_field_domain(this % sites, this % physics % num_components())
       direction = points % direction(this % formed_at_sites(v))
       at_points = variation(this % physics % argument(1), direction)
    else if (given % argument_is(this % argument(1))) then
       points    = typed_field_domain(this % points, this % stride())
       direction = points % direction(this % gathered(v))
       at_points = variation(this % physics % argument(1), direction)
    else if (given % argument_is(this % argument(2)) .and. allocated(this % projected)) then
       points    = typed_field_domain(this % sites, this % designs)
       direction = points % direction(this % design_at_sites(v))
       at_points = variation(this % physics % argument(2), direction)
    else if (given % argument_is(this % argument(2))) then
       at_points = given % with_argument(this % physics % argument(2))
    else
       error stop 'operation_residual: physics_variation received a variation naming neither the state nor the design'
    end if
  end function physics_variation

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
    if (allocated(this % projected)) then
       error stop 'operation_residual: a residual projected on the cells is not restricted; its instances &
            &read the nodes outside the selection'
    end if
    if (allocated(this % prescribed)) then
       error stop 'operation_residual: a residual with prescribed rows is not restricted; the &
            &prescribed relation would lose its unknowns outside the selection'
    end if
    if (allocated(this % connected_law)) then
       secondary = this % connected_law % restricted(free, values)
       sub = residual_operator(derived, this % physics, at, size(free), &
            & this % degrees, this % primary, fixed_rows(1:ncar), fixed(1:ncar), &
            & connected_law=secondary, fixed_rate=rate(1:ncar), positions=this % point_positions(:, points))
    else
       sub = residual_operator(derived, this % physics, at, size(free), &
            & this % degrees, this % primary, fixed_rows(1:ncar), fixed(1:ncar), &
            & fixed_rate=rate(1:ncar), positions=this % point_positions(:, points))
    end if

    if (allocated(this % summed_rows)) sub % summed_rows = pack(sub_of(this % summed_rows), sub_of(this % summed_rows) > 0)
    call sub % explicit_at(this, points)

  end function constrain

  !===================================================================!
  ! The rules at the explicit points of a residual, at the points of it
  ! given: those of a part or of a constraint.
  !===================================================================!

  subroutine explicit_at(this, whole, points)
    class(residual_operator), intent(inout) :: this
    type(residual_operator) , intent(in)    :: whole
    integer                 , intent(in)    :: points(:)
    if (.not. allocated(whole % explicit)) return
    if (.not. any(whole % explicit(points))) return
    this % explicit_physics = whole % explicit_physics
    this % explicit_rules   = whole % explicit_rules
    this % explicit         = whole % explicit(points)
  end subroutine explicit_at

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
    if (allocated(this % summed_rows)) lin % summed_rows = this % summed_rows
    lin % points           = this % points
    lin % point_positions  = this % point_positions
    if (allocated(this % tying)) lin % tying = this % tying
    call lin % versioned(version_number, transposed=transposed)

  end function linearize

  !===================================================================!
  ! THE PROJECTION'S SITES. The tuples the rules read at every site,
  ! formed from the state by the basis, one tuple after another; the
  ! design at every site, that of its instance's first node.
  !===================================================================!

  pure function formed_at_sites(this, x) result(y)
    class(residual_operator), intent(in) :: this
    real(dp)                , intent(in) :: x(:)
    real(dp), allocatable :: y(:)
    integer :: width, nq, n, q, s, c, i, k, p, d, offset, e
    associate (pr => this % projected, el => this % projected % elements)
      width = size(pr % formed, 2)
      nq    = el % num_quadrature_points()
      allocate(y(width * size(pr % element_cell) * nq), source=0.0_dp)
      do n = 1, size(pr % element_cell)
         e = pr % element_cell(n)
         do q = 1, nq
            s = (n - 1) * nq + q
            do c = 1, width
               d = pr % formed(0, c)
               do i = 1, size(pr % element_point, 1)
                  p = abs(pr % element_point(i, n))
                  do k = 0, el % num_coefficients() - 1
                     offset = pr % formed(k + 1, c)
                     if (offset < 0) cycle
                     y((s - 1) * width + c) = y((s - 1) * width + c) + el % basis(d, k, i, q, e) &
                          & * x(this % at(p) + offset + 1)
                  end do
               end do
            end do
         end do
      end do
    end associate
  end function formed_at_sites

  pure function design_at_sites(this, v) result(w)
    class(residual_operator), intent(in) :: this
    real(dp)                , intent(in) :: v(:)
    real(dp), allocatable :: w(:)
    integer :: nq, n, q, s, p, nd
    nd = this % designs
    nq = this % projected % elements % num_quadrature_points()
    allocate(w(nd * size(this % projected % element_cell) * nq))
    do n = 1, size(this % projected % element_cell)
       p = abs(this % projected % element_point(1, n))
       do q = 1, nq
          s = (n - 1) * nq + q
          w((s - 1) * nd + 1:s * nd) = v((p - 1) * nd + 1:p * nd)
       end do
    end do
  end function design_at_sites

  !===================================================================!
  ! Each rule's value at every site, times the weight and the test
  ! function of every coefficient of every node the instance states,
  ! added to that coefficient's row.
  !===================================================================!

  pure subroutine tested_into(this, governing, r)
    class(residual_operator), intent(in)    :: this
    real(dp)                , intent(in)    :: governing(:,:)
    real(dp)                , intent(inout) :: r(:)
    integer :: j, f, d, n, e, q, s, i, p, k, nq
    real(dp) :: value
    associate (pr => this % projected, el => this % projected % elements)
      nq = el % num_quadrature_points()
      do j = 1, size(this % rules)
         f = pr % tested(1, j)
         d = pr % tested(2, j)
         do n = 1, size(pr % element_cell)
            e = pr % element_cell(n)
            do q = 1, nq
               s = (n - 1) * nq + q
               if (.not. tested_at(el, d, q, e)) cycle
               value = governing(s, j)
               do i = 1, size(pr % element_point, 1)
                  p = pr % element_point(i, n)
                  if (p < 1) cycle
                  do k = 0, el % num_coefficients() - 1
                     r(this % at(p) + pr % coefficient_row(k, f) + 1) = &
                          & r(this % at(p) + pr % coefficient_row(k, f) + 1) + value * el % test(d, k, i, q, e)
                  end do
               end do
            end do
         end do
      end do
    end associate
  end subroutine tested_into

  !===================================================================!
  ! THE JACOBIAN OF THE PROJECTION, element by element. At each site of
  ! an instance, the gradient of each rule in the tuple it reads, times
  ! the weight and the test function of a row and the basis's
  ! derivative that forms each component from each coefficient of each
  ! node, is summed into the instance's matrix: its rows the node, the
  ! coefficient and the equation; its columns the node and the
  ! component of a node's tuple the basis reads. The entries of the
  ! rows the instance states are then written once each, those of one
  ! row and column of two instances summed when combined.
  !===================================================================!

  integer function projected_count(this) result(count)
    class(residual_operator), intent(in) :: this
    integer , allocatable :: offsets(:), slot(:,:), columns_of(:), first(:), column(:)
    integer :: nodes, nk
    nodes = size(this % projected % element_point, 1)
    nk    = this % projected % elements % num_coefficients()
    call this % projected_columns(offsets, slot, columns_of, first, column)
    count = size(this % projected % element_cell) * nodes * nk * nodes * sum(columns_of)
  end function projected_count

  !===================================================================!
  ! THE COLUMNS EACH FIELD'S ROWS READ in an instance: of the
  ! components of a node's tuple the basis reads (offsets), those the
  ! rules testing the field read, in the order of offsets.
  ! columns_of(f) of them, column(first(f):first(f) + columns_of(f) -
  ! 1) their places in offsets, and slot(o, f) the rank of the o-th
  ! offset among the field's, zero where it reads none.
  !===================================================================!

  subroutine projected_columns(this, offsets, slot, columns_of, first, column)
    class(residual_operator), intent(in)  :: this
    integer , allocatable   , intent(out) :: offsets(:), slot(:,:), columns_of(:), first(:), column(:)
    integer , allocatable :: reads(:), place(:)
    integer :: no, nf, j, f, m, k2, o, r
    associate (pr => this % projected)
      offsets = read_offsets(pr % formed)
      no      = size(offsets)
      nf      = size(pr % coefficient_row, 2)
      allocate(place(0:maxval([offsets, 0])), source=0)
      place(offsets) = [(o, o = 1, no)]
      allocate(slot(no, nf), source=0)
      do j = 1, size(this % rules)
         f = pr % tested(1, j)
         call this % rules(j) % read_components(reads)
         call slotted()
         ! and those the rule at the explicit points reads
         if (allocated(this % explicit_rules)) then
            call this % explicit_rules(j) % read_components(reads)
            call slotted()
         end if
      end do
    end associate
    allocate(columns_of(nf), first(nf), column(count(slot > 0)))
    r = 0
    do f = 1, nf
       first(f) = r + 1
       do o = 1, no
          if (slot(o, f) == 0) cycle
          r = r + 1
          column(r)  = o
          slot(o, f) = r - first(f) + 1
       end do
       columns_of(f) = r - first(f) + 1
    end do
  contains
    subroutine slotted()
      do m = 1, size(reads)
         do k2 = 1, ubound(this % projected % formed, 1)
            if (this % projected % formed(k2, reads(m) + 1) < 0) cycle
            slot(place(this % projected % formed(k2, reads(m) + 1)), f) = 1
         end do
      end do
    end subroutine slotted
  end subroutine projected_columns

  ! the components of a node's tuple the basis reads, each once
  pure function read_offsets(formed) result(offsets)
    integer, intent(in) :: formed(0:,:)
    integer, allocatable :: offsets(:)
    integer :: c, k
    offsets = [integer ::]
    do c = 1, size(formed, 2)
       do k = 1, ubound(formed, 1)
          if (formed(k, c) < 0) cycle
          if (any(offsets == formed(k, c))) cycle
          offsets = [offsets, formed(k, c)]
       end do
    end do
  end function read_offsets

  subroutine projected_triples(this, x, nu, is_fixed, r, c, w, num_triples)
    class(residual_operator), intent(in)    :: this
    real(dp)                , intent(in)    :: x(:), nu(:)
    logical                 , intent(in)    :: is_fixed(:)
    integer                 , intent(inout) :: r(:), c(:)
    real(dp)                , intent(inout) :: w(:)
    integer                 , intent(inout) :: num_triples
    real(dp), allocatable :: xs(:), ns(:), g(:), local(:)
    integer , allocatable :: reads(:), offsets(:), place(:), slot(:,:), columns_of(:), first(:), column(:), base(:)
    integer :: width, nd, nq, nodes, nk, nf, no, j, f, d, n, e, q, s, i, p, k, m, comp, dd, i2, k2, offset, row, o
    integer :: at_row
    real(dp) :: value, tau
    logical , allocatable :: marked(:)
    logical :: explicit_site
    width   = size(this % projected % formed, 2)
    marked  = this % explicit_marks()
    nd      = this % designs
    xs      = this % formed_at_sites(x)
    ns      = this % design_at_sites(nu)
    call this % projected_columns(offsets, slot, columns_of, first, column)
    no      = size(offsets)
    allocate(place(0:maxval([offsets, 0])), source=0)
    place(offsets) = [(o, o = 1, no)]
    allocate(g(0:width - 1))
    associate (pr => this % projected, el => this % projected % elements)
      nq    = el % num_quadrature_points()
      nodes = size(pr % element_point, 1)
      nk    = el % num_coefficients()
      nf    = size(pr % coefficient_row, 2)
      ! the instance's matrix, a block for each field: its rows the
      ! node and the coefficient, its columns the node and the field's
      ! own offsets
      allocate(base(nf + 1))
      base(1) = 0
      do f = 1, nf
         base(f + 1) = base(f) + nodes * nk * nodes * columns_of(f)
      end do
      allocate(local(base(nf + 1)))
      do n = 1, size(pr % element_cell)
         e     = pr % element_cell(n)
         local = 0.0_dp
         do q = 1, nq
            s = (n - 1) * nq + q
            do j = 1, size(this % rules)
               f = pr % tested(1, j)
               d = pr % tested(2, j)
               if (.not. tested_at(el, d, q, e)) cycle
               ! at an explicit site, the rule there
               explicit_site = .false.
               if (allocated(this % explicit_rules)) explicit_site = marked(s)
               if (explicit_site) then
                  call this % explicit_rules(j) % read_components(reads)
                  call this % explicit_rules(j) % gradient_at(xs((s - 1) * width + 1:s * width), &
                       & ns((s - 1) * nd + 1:s * nd), value, g, this % site_positions(:, s))
               else
                  call this % rules(j) % read_components(reads)
                  call this % rules(j) % gradient_at(xs((s - 1) * width + 1:s * width), &
                       & ns((s - 1) * nd + 1:s * nd), value, g, this % site_positions(:, s))
               end if
               do i = 1, nodes
                  do k = 0, nk - 1
                     tau = el % test(d, k, i, q, e)
                     if (tau == 0.0_dp) cycle
                     at_row = base(f) + ((i - 1) * nk + k) * nodes * columns_of(f)
                     do m = 1, size(reads)
                        comp = reads(m)
                        dd   = pr % formed(0, comp + 1)
                        do i2 = 1, nodes
                           do k2 = 0, nk - 1
                              offset = pr % formed(k2 + 1, comp + 1)
                              if (offset < 0) cycle
                              o = at_row + (i2 - 1) * columns_of(f) + slot(place(offset), f)
                              local(o) = local(o) + tau * g(comp) * el % basis(dd, k2, i2, q, e)
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
         ! the rows the instance states, once each
         do f = 1, nf
            do i = 1, nodes
               p = pr % element_point(i, n)
               if (p < 1) cycle
               do k = 0, nk - 1
                  row = this % at(p) + pr % coefficient_row(k, f) + 1
                  if (is_fixed(row)) cycle
                  at_row = base(f) + ((i - 1) * nk + k) * nodes * columns_of(f)
                  do i2 = 1, nodes
                     do o = 1, columns_of(f)
                        if (local(at_row + (i2 - 1) * columns_of(f) + o) == 0.0_dp) cycle
                        num_triples    = num_triples + 1
                        r(num_triples) = row
                        c(num_triples) = this % at(abs(pr % element_point(i2, n))) + offsets(column(first(f) + o - 1)) + 1
                        w(num_triples) = local(at_row + (i2 - 1) * columns_of(f) + o)
                     end do
                  end do
               end do
            end do
         end do
      end do
    end associate
  end subroutine projected_triples

  !===================================================================!
  ! The edges a projected law adds to the coupling: from every
  ! component a node's tuple forms at a site to the row of every
  ! coefficient of every equation of every node of the instance.
  !===================================================================!

  subroutine projected_edges(this, heads, tails)
    class(residual_operator), intent(in)    :: this
    integer , allocatable   , intent(inout) :: heads(:), tails(:)
    integer , allocatable :: offsets(:), h(:), t(:)
    integer :: n, i, p, f, k, i2, p2, o, e
    associate (pr => this % projected, el => this % projected % elements)
      offsets = [integer ::]
      do o = 1, size(pr % formed, 2)
         do k = 1, size(pr % formed, 1) - 1
            if (pr % formed(k, o) < 0) cycle
            if (any(offsets == pr % formed(k, o))) cycle
            offsets = [offsets, pr % formed(k, o)]
         end do
      end do
      allocate(h(size(pr % element_cell) * size(pr % element_point, 1)**2 * size(pr % coefficient_row) &
           & * size(offsets)), t(size(pr % element_cell) * size(pr % element_point, 1)**2 &
           & * size(pr % coefficient_row) * size(offsets)))
      e = 0
      do n = 1, size(pr % element_cell)
         do i = 1, size(pr % element_point, 1)
            p = pr % element_point(i, n)
            if (p < 1) cycle
            do f = 1, size(pr % coefficient_row, 2)
               do k = 0, el % num_coefficients() - 1
                  do i2 = 1, size(pr % element_point, 1)
                     p2 = abs(pr % element_point(i2, n))
                     do o = 1, size(offsets)
                        e    = e + 1
                        h(e) = this % at(p) + pr % coefficient_row(k, f) + 1
                        t(e) = this % at(p2) + offsets(o) + 1
                     end do
                  end do
               end do
            end do
         end do
      end do
    end associate
    heads = [heads, h(1:e)]
    tails = [tails, t(1:e)]
  end subroutine projected_edges

  !===================================================================!
  ! THE PROJECTION OF A PART: the instances with a node the part owns,
  ! each node's point in the part's numbering, negative for a node of
  ! the halo, whose rows another part states.
  !===================================================================!

  function projected_part(this, place_of, owned) result(sub)
    class(residual_operator), intent(in) :: this
    ! the part's point of each point, zero outside the part, and
    ! whether the part owns it
    integer                 , intent(in) :: place_of(:)
    logical                 , intent(in) :: owned(:)
    type(projection) :: sub
    logical, allocatable :: stated(:)
    integer :: n, i, p, m
    associate (pr => this % projected)
      allocate(stated(size(pr % element_cell)), source=.false.)
      do n = 1, size(pr % element_cell)
         do i = 1, size(pr % element_point, 1)
            p = abs(pr % element_point(i, n))
            if (place_of(p) > 0) then
               if (owned(place_of(p)) .and. pr % element_point(i, n) > 0) stated(n) = .true.
            end if
         end do
      end do
      sub % elements        = pr % elements
      sub % tested          = pr % tested
      sub % coefficient_row = pr % coefficient_row
      sub % formed          = pr % formed
      sub % element_cell    = pack(pr % element_cell, stated)
      allocate(sub % element_point(size(pr % element_point, 1), count(stated)))
      m = 0
      do n = 1, size(pr % element_cell)
         if (.not. stated(n)) cycle
         m = m + 1
         do i = 1, size(pr % element_point, 1)
            p = abs(pr % element_point(i, n))
            if (place_of(p) == 0) then
               error stop 'operation_residual: a part reads every node of an instance it states'
            end if
            sub % element_point(i, m) = place_of(p)
            if (.not. owned(place_of(p)) .or. pr % element_point(i, n) < 1) sub % element_point(i, m) = -place_of(p)
         end do
      end do
    end associate
  end function projected_part

  !===================================================================!
  ! A term: the equations on a manifold; with inequality, on the
  ! design factor, each equation at most zero. Invalid input: an
  ! equation that is not a scalar field, or one on neither the
  ! manifold nor its parent; an inequality on another part.
  !===================================================================!

  function create_continuous(manifold, equations, inequality) result(this)

    type(continuous_manifold), intent(in) :: manifold
    type(continuous_field)   , intent(in) :: equations(:)
    logical                  , intent(in), optional :: inequality
    type(continuous_residual) :: this

    integer :: k
    character(len=250) :: message

    allocate(this % term(1))
    this % term(1) % manifold = manifold
    this % term(1) % equation = equations
    if (present(inequality)) this % term(1) % inequality = inequality
    if (this % term(1) % inequality .and. .not. manifold % is_design_factor()) then
       error stop 'operation_residual: an inequality is a condition on the design factor'
    end if
    do k = 1, size(equations)
       ! a smooth maximum of several components is a scalar
       if (equations(k) % num_components() /= 1 .and. .not. equations(k) % weight > 0.0_dp) then
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
  ! manifold, or the smooth maximum of several. Invalid input: a field
  ! that is neither, or a residual with an objective already.
  !===================================================================!

  function objective_sum(objective, residual) result(this)

    type(continuous_field)   , intent(in) :: objective
    type(continuous_residual), intent(in) :: residual
    type(continuous_residual) :: this

    if (.not. objective % is_integral()) then
       error stop 'operation_residual: the objective is the integral of a field over the manifold, or its &
            &smooth maximum'
    end if
    if (objective % num_components() /= 1 .and. .not. objective % weight > 0.0_dp) then
       error stop 'operation_residual: the objective is the integral of a scalar field'
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
  ! THE DISCRETIZATION. Exactly one term on the whole manifold that
  ! is not a representation law gives the rule, one equation per
  ! unknown in order; every other term is a condition on a part of
  ! that manifold, paired or bare. A coordinate present without its
  ! approximation, or an approximation given for a coordinate absent,
  ! is invalid input.
  !===================================================================!

  !===================================================================!
  ! THE REPRESENTATION LAWS STATED, each column of defined the field,
  ! the order and the coordinate of the jet component one law defines.
  ! Along a coordinate with laws stated, every component the rule
  ! stores along it, of order one to its degree, is defined exactly
  ! once, and the coordinate has its concretion: the chain of families
  ! along time, a derivative approximation along space. The rows stay
  ! those the concretion builds. Invalid input: a component defined
  ! twice or not at all, or a law along a coordinate without its
  ! concretion.
  !===================================================================!

  ! the equations of the representation laws paired with a multiplier
  pure integer function num_represented(this)

    class(discrete_residual), intent(in) :: this

    integer :: r

    num_represented = 0
    if (.not. allocated(this % representation)) return
    do r = 1, size(this % representation)
       num_represented = num_represented + size(this % representation(r) % equation)
    end do

  end function num_represented

  subroutine check_representations(this, defined, with_chain, with_differences)

    class(discrete_residual), intent(in) :: this
    integer                 , intent(in) :: defined(:,:)
    logical                 , intent(in) :: with_chain, with_differences

    integer :: c, f, o, top, times
    logical :: along_time
    character(len=250) :: message

    do c = 1, this % rule % num_coordinates()
       if (.not. any(defined(3, :) == c)) cycle
       along_time = this % manifold % with_time .and. c == 1
       if ((along_time .and. .not. with_chain) .or. (.not. along_time .and. .not. with_differences)) then
          write(message,'(a,i0)') 'operation_residual: a representation law is stated along a coordinate &
               &that has no concretion; coordinate = ', c
          error stop trim(message)
       end if
       do f = 1, this % manifold % num_unknowns
          top = merge(this % rule % degree_of_field(f), this % rule % degree_along(c), along_time)
          do o = 1, top
             times = count(defined(1, :) == f .and. defined(2, :) == o .and. defined(3, :) == c)
             if (times /= 1) then
                write(message,'(a,i0,a,i0,a,i0,a,i0,a)') 'operation_residual: along coordinate ', c, &
                     & ' the component of order ', o, ' of field ', f, ' is defined ', times, &
                     & ' times; a representation law defines each stored component once'
                error stop trim(message)
             end if
          end do
       end do
       if (count(defined(3, :) == c) /= sum([(merge(this % rule % degree_of_field(f), &
            & this % rule % degree_along(c), along_time), f = 1, this % manifold % num_unknowns)])) then
          write(message,'(a,i0)') 'operation_residual: a representation law defines a component the rule &
               &does not store; coordinate = ', c
          error stop trim(message)
       end if
    end do

  end subroutine check_representations

  function residual_discretize(this, points, time, space) result(image)

    class(continuous_residual), intent(in) :: this
    type(discrete_manifold)   , intent(in) :: points
    type(chain)               , intent(in), optional :: time
    class(derivative_approximation), intent(in), optional :: space
    type(discrete_residual) :: image

    integer :: k, j, i, whole_term, num_conditions, nf
    integer, allocatable :: defined(:,:)
    logical, allocatable :: representation(:)
    logical :: covering, projected
    type(expression) :: form
    real(dp) :: given
    real(dp), allocatable :: slope(:)
    character(len=250) :: message

    ! the terms on the whole manifold: the representation laws, every
    ! equation of which defines a jet component, and the one term of the
    ! physics
    allocate(representation(size(this % term)), source=.false.)
    allocate(defined(3, 0))
    whole_term = 0
    num_conditions = 0
    do k = 1, size(this % term)
       associate (t => this % term(k))
         if (t % manifold % part == WHOLE) then
            representation(k) = all([(represented_component(t % equation(i) % graph(1)) /= 0, &
                 & i = 1, size(t % equation))])
            if (representation(k)) then
               do i = 1, size(t % equation)
                  defined = reshape([defined, represented_component(t % equation(i) % graph(1))], &
                       & [3, size(defined, 2) + 1])
               end do
               cycle
            end if
            if (whole_term /= 0) then
               error stop 'operation_residual: one term states the physics on the whole manifold'
            end if
            whole_term = k
         else
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
      ! every derivative operator replaced by the jet components it
      ! denotes, here and nowhere else: the rule, the conditions and the
      ! objective are evaluated in their jet form. Under a projection on
      ! the cells the projection forms the rule from the equations as
      ! they are, the outer derivatives along space read as fluxes and
      ! never formed
      projected = .false.
      if (present(space)) projected = space % projects()
      if (.not. projected) image % rule = paired([(jet_form(t % equation(j) % graph(1)), j = 1, size(t % equation))])
    end associate

    if (projected) call projected_form(image, this % term(whole_term), space, points, defined)
    call image % check_representations(defined, present(time), present(space))
    if (any(image % manifold % of_time_alone(1:nf)) .and. .not. allocated(image % elements)) then
       error stop 'operation_residual: an unknown of time alone on a manifold of space is stated under a &
            &projection on the cells, which tests its equation with the weight alone'
    end if
    if (present(time) .and. .not. allocated(image % elements)) call initial_form(image, this % term(whole_term))
    image % representation = pack(this % term, representation .and. [(this % term(k) % paired, k = 1, size(this % term))])

    allocate(image % condition(num_conditions))
    j = 0
    do k = 1, size(this % term)
       if (k == whole_term .or. representation(k)) cycle
       j = j + 1
       image % condition(j) = this % term(k)
       do i = 1, size(image % condition(j) % equation)
          image % condition(j) % equation(i) = image % condition(j) % equation(i) % jet_form()
       end do
       if (.not. image % condition(j) % manifold % parent % matches(image % manifold % identity)) then
          error stop 'operation_residual: a condition is stated on a part of another manifold'
       end if
       ! a condition on a face of the region takes the rows of the
       ! physics at the points nearest the face: the derivatives along
       ! space come from a closed approximation at every point
       if (image % condition(j) % manifold % is_space_face() .and. .not. present(space)) then
          error stop 'operation_residual: a condition on a face of the region requires the derivatives along &
               &space approximated at every point, a derivative approximation'
       end if
       ! an inequality on the design factor is a function of the design
       ! coordinates alone, or a functional of the state over the whole
       ! manifold, an integral or a smooth maximum
       if (image % condition(j) % inequality) then
          do i = 1, size(image % condition(j) % equation)
             associate (equation => image % condition(j) % equation(i))
               form = equation % graph(1)
               if (equation % is_integral()) then
                  if (equation % integrated_part /= WHOLE .or. &
                       & .not. equation % integrated_over % matches(image % manifold % identity)) then
                     error stop 'operation_residual: an inequality on a functional of the state is over the whole &
                          &manifold'
                  end if
               else if (form % reads_state()) then
                  error stop 'operation_residual: an inequality on the design factor is a function of the design &
                       &coordinates alone, or a functional of the state over the manifold'
               end if
             end associate
          end do
       ! a condition on the design factor fixes the design coordinate
       ! at the value the manifold states: nu - nu_0 = 0, one equation
       else if (image % condition(j) % manifold % is_design_factor()) then
          if (size(image % condition(j) % equation) /= points % expansion % num_designs()) then
             error stop 'operation_residual: the condition on the design factor fixes every design coordinate, &
                  &one equation per coordinate in the order of the coordinates'
          end if
          do i = 1, points % expansion % num_designs()
             call design_condition(image % condition(j) % equation(i), points % expansion % values(), given, slope)
             if (abs(given) > 1.0e-12_dp * max(1.0_dp, maxval(abs(points % expansion % values())))) then
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
               & .and. .not. image % manifold % with_space
       case (SPACE_FACTOR)
          covering = this % objective % integrated_parent % matches(image % manifold % identity) &
               & .and. .not. image % manifold % with_time
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
       image % objective = this % objective % jet_form()
       if (allocated(image % elements) .and. this % objective % integrated_part /= WHOLE) then
          error stop 'operation_residual: under a projection on the cells the objective is integrated over the &
               &whole manifold, at the quadrature points of the cells'
       end if
    end if

    image % points = points
    if (image % manifold % with_time .neqv. present(time)) then
       error stop 'operation_residual: the derivatives along time are approximated by a chain of &
            &families, given exactly when the manifold has a time coordinate'
    end if
    if (image % manifold % with_space .neqv. present(space)) then
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
  ! THE ZERO OF THE DISCRETE RESIDUAL, at the design the image states,
  ! or, when no condition is stated on the design, at the design where
  ! dJ/ds = 0: the stationarity of L_h in the design once the
  ! condition that would hold it is removed. That design is found by
  ! Newton's iteration on s from the value the manifold states, each
  ! iterate the zero with the design prescribed there, the gradient -kappa
  ! and the hessian -dkappa/ds from the expansion of order one along
  ! each design coordinate: the reduced-space method. The iteration
  ! ends where |dJ/ds| is below 1e-10 times the larger of one and its
  ! first value, and the solution is the zero at that design. A free
  ! design requires an objective, the equations paired with their
  ! multiplier, and the expansion of order one along each design
  ! coordinate.
  !
  ! UNDER INEQUALITIES c(s) <= 0 on the design factor, each step is the
  ! minimizer of the quadratic model of J under the inequalities
  ! linearized at the iterate, c + A d <= 0, with their multipliers mu,
  ! nonnegative and zero where an inequality is not active; the
  ! iteration ends where |dJ/ds + A^T mu| is below the same tolerance
  ! and no inequality exceeds 1e-12 of its scale. The multipliers are
  ! returned under the inequality's name, those of the design returned:
  ! the stationarity of L = J + mu . c in the design there, by least
  ! squares on the active inequalities. Invalid input: an inequality
  ! and an equality on the design factor together.
  !
  ! An inequality on a FUNCTIONAL OF THE STATE - an integral or a
  ! smooth maximum over the manifold - has its value and its gradient
  ! from the expansion of the state, and enters the reverse sweep times
  ! its multiplier of the step before: the hessian of the model is that
  ! of the Lagrangian J + mu . c. The merit of an iterate is J plus
  ! each violation times its penalty, at least twice the inequality's
  ! own multiplier, zero while it has not been active - the next
  ! step's linearization then restores it: the inequalities are of
  ! their own units, and one penalty of the largest multiplier weighed
  ! the rotor's smooth maximum of the stress, of multiplier 0.04, by
  ! the 34.6 of its bounds on the thickness. A step that does not
  ! decrease the merit by a tenth of the decrease its model predicts
  ! is taken again from the design before it, at most half as long,
  ! with the hessian shifted by sigma times the identity, sigma raised
  ! fourfold from at least |dJ/ds| / (0.1 max(1, |s|)), and with each
  ! violated inequality linearized to remove a quarter of the fraction
  ! of its violation the step before removed, c + A d <= (1 - theta) c
  ! - the shift does not shorten a step its active inequalities fix; a
  ! hessian that is not positive definite is shifted from that value
  ! too, and the shift falls and theta rises fourfold, to one, after a
  ! step that decreases the merit by more than three quarters of its
  ! prediction. A positive definite hessian whose steps all decrease
  ! the merit is never shifted, and theta is one. Where no step meets
  ! the inequalities linearized, theta is quartered for that step, three
  ! times, and then zero: each violated inequality is then required
  ! only not to grow, which the null step meets.
  !
  ! CONTINUED, each block of the chain after the first starts Newton's
  ! iteration from the solution at the last instant before it, its
  ! unknowns of time alone extrapolated along time to the block's last
  ! instant, at its instants and the stages of its steps, rather than
  ! from the estimate, which then gives the first block alone, its
  ! derivatives along time those of the estimate: the march of a
  ! motion the estimate does not follow.
  !===================================================================!

  recursive subroutine minimize(this, estimate, solution, checkpoint, snapshots, continued)

    class(discrete_residual), intent(in)  :: this
    type(discrete_field)    , intent(in)  :: estimate
    type(discrete_field)    , intent(out) :: solution
    integer                 , intent(in), optional :: checkpoint, snapshots
    logical                 , intent(in), optional :: continued

    type(discrete_residual) :: iterate(1)
    type(discrete_field)    :: estimates(1), solutions(1)

    if (present(continued)) then
       if (continued .neqv. this % continued) then
          select type (this)
          type is (discrete_residual)
             iterate(1) = this
          end select
          iterate(1) % continued = continued
          call iterate(1) % minimize(estimate, solution, checkpoint, snapshots)
          return
       end if
    end if
    if (.not. free_design(this)) then
       if (any(this % condition % inequality)) then
          error stop 'operation_residual: an inequality bounds a free design; an equality on the design &
               &factor fixes it and leaves the inequality nothing to bound'
       end if
       call this % zero_at_design(estimate, solution, checkpoint, snapshots)
       return
    end if
    if (.not. (allocated(this % objective) .and. this % with_adjoint)) then
       error stop 'operation_residual: a design on which no condition is stated is free, and is found where &
            &dJ/ds = 0, which requires an objective and the equations paired with their multiplier'
    end if
    select type (this)
    type is (discrete_residual)
       iterate(1) = this
    end select
    ! the estimate in an array of its own: a constructor passed straight
    ! as an argument is a temporary the compiler in use mishandles
    estimates(1) = estimate
    call design_iteration(iterate, estimates, solutions, checkpoint, snapshots)
    solution = solutions(1)

  end subroutine minimize

  !===================================================================!
  ! SEVERAL LAGRANGIANS ON ONE DESIGN: residuals on manifolds of the
  ! same design coordinates, each with its own estimate and solution,
  ! the design free in every one - the scenarios of one design. The
  ! design is found where the sum of the Lagrangians is stationary in
  ! it: the sum of their objectives, any of which may have none, least
  ! under the inequalities of all, by the iteration of minimize, every
  ! residual's forward and reverse sweeps at each iterate and its
  ! multipliers returned with its solution. Invalid input: residuals
  ! of different design coordinates, a residual with an equality on
  ! the design, or none with an objective.
  !===================================================================!

  subroutine minimize_together(residuals, estimates, solutions, checkpoint, snapshots, continued)

    type(discrete_residual), intent(in)  :: residuals(:)
    type(discrete_field)   , intent(in)  :: estimates(:)
    type(discrete_field)   , intent(out) :: solutions(:)
    integer                , intent(in), optional :: checkpoint, snapshots
    logical                , intent(in), optional :: continued

    type(discrete_residual), allocatable :: iterates(:)
    integer :: k

    if (size(estimates) /= size(residuals) .or. size(solutions) /= size(residuals)) then
       error stop 'operation_residual: one estimate and one solution for each residual minimized together'
    end if
    iterates = residuals
    do k = 1, size(iterates)
       if (present(continued)) iterates(k) % continued = continued
       if (.not. free_design(iterates(k))) then
          error stop 'operation_residual: residuals minimized together have a free design, no equality on it'
       end if
       if (.not. iterates(k) % with_adjoint) then
          error stop 'operation_residual: residuals minimized together have their equations paired with their &
               &multiplier'
       end if
       if (iterates(k) % points % expansion % num_designs() /= iterates(1) % points % expansion % num_designs()) then
          error stop 'operation_residual: residuals minimized together have the same design coordinates'
       end if
       if (any(iterates(k) % points % expansion % values() /= iterates(1) % points % expansion % values())) then
          error stop 'operation_residual: residuals minimized together begin at the same design'
       end if
    end do
    if (.not. any([(allocated(iterates(k) % objective), k = 1, size(iterates))])) then
       error stop 'operation_residual: residuals minimized together have an objective, one at least'
    end if
    call design_iteration(iterates, estimates, solutions, checkpoint, snapshots)

  end subroutine minimize_together

  !===================================================================!
  ! THE ITERATION ON THE FREE DESIGN of one residual or of several on
  ! one design, as minimize states it: the gradient and the hessian
  ! the sums of every residual's, the inequalities those of all, the
  ! merit the sum of the objectives plus each violation times its
  ! penalty.
  !===================================================================!

  subroutine design_iteration(iterates, estimates, solutions, checkpoint, snapshots)

    type(discrete_residual), intent(inout) :: iterates(:)
    type(discrete_field)   , intent(in)    :: estimates(:)
    type(discrete_field)   , intent(out)   :: solutions(:)
    integer                , intent(in), optional :: checkpoint, snapshots

    type(discrete_field)    :: functional
    real(dp), allocatable :: kappa(:,:), gradient(:), hessian(:,:), c(:), a(:,:), step(:), mu(:), weight(:)
    real(dp), allocatable :: accepted(:), accepted_gradient(:), accepted_hessian(:,:), accepted_c(:), accepted_a(:,:)
    real(dp), allocatable :: ck(:), ak(:,:), functional_kappa(:,:), quasi_hessian(:,:), penalty(:)
    integer , allocatable :: active(:), first(:)
    logical , allocatable :: of_state(:)
    real(dp) :: tolerance, stationarity, objective, accepted_objective, shift, reach, predicted, decrease
    real(dp) :: reduction, limit, theta
    integer :: n, j, k, iteration, quartering
    logical :: with_accepted, refused, quasi, solvable
    character(len=250) :: message

    ! the hessian of the model: the expansion's of order one along each
    ! design coordinate, or, of order zero, the damped BFGS update's of
    ! the Lagrangian's gradient, each functional's by its own reverse
    ! sweep
    quasi = all([(iterates(k) % points % expansion % order == 0, k = 1, size(iterates))])
    do k = 1, size(iterates)
       if (.not. quasi .and. (iterates(k) % points % expansion % order < 1 .or. &
            & .not. along_each_coordinate(iterates(k) % points % expansion))) then
          error stop 'operation_residual: a free design is found by Newton''s iteration, whose hessian is the &
               &expansion of order one along each design coordinate, or, of order zero, its BFGS update'
       end if
    end do
    n = iterates(1) % points % expansion % num_designs()
    allocate(gradient(n), hessian(n, n), quasi_hessian(n, n), source=0.0_dp)
    tolerance = 0.0_dp
    ! the inequalities of every residual in turn, first(k) the first of
    ! the k-th residual's
    allocate(of_state(0), first(size(iterates) + 1))
    first(1) = 1
    do k = 1, size(iterates)
       of_state = [of_state, iterates(k) % functional_inequalities()]
       first(k + 1) = size(of_state) + 1
    end do
    allocate(weight(size(of_state)), penalty(size(of_state)), source=0.0_dp)
    shift     = 0.0_dp
    reduction = 1.0_dp
    predicted = 0.0_dp
    with_accepted = .false.
    do iteration = 1, MAX_DESIGN_ITERATIONS
       gradient  = 0.0_dp
       hessian   = 0.0_dp
       objective = 0.0_dp
       allocate(c(0), a(0, n))
       do k = 1, size(iterates)
          if (quasi) then
             call iterates(k) % zero_at_design(estimates(k), solutions(k), checkpoint, snapshots, kappa, &
                  & functional_kappa)
             call iterates(k) % inequalities_at(solutions(k), ck, ak, -functional_kappa)
          else
             if (any(of_state(first(k):first(k + 1) - 1))) &
                  & iterates(k) % functional_weight = weight(first(k):first(k + 1) - 1)
             call iterates(k) % zero_at_design(estimates(k), solutions(k), checkpoint, snapshots, kappa)
             call iterates(k) % inequalities_at(solutions(k), ck, ak)
          end if
          c = [c, ck]
          a = reshape([transpose(a), transpose(ak)], [size(c), n], order=[2, 1])
          ! the gradient of J: the Lagrangian's less the functionals' times
          ! the multipliers the reverse sweep read
          gradient = gradient - kappa(:, 0)
          if (.not. quasi) then
             do j = 1, n
                hessian(:, j) = hessian(:, j) - kappa(:, iterates(k) % points % expansion % column(j, 1))
             end do
          end if
          if (allocated(iterates(k) % objective) .and. (size(of_state) > 0 .or. quasi)) then
             functional = iterates(k) % objective % at(solutions(k))
             objective  = objective + functional % value(1, 1)
          end if
       end do
       if (verbosity >= 1) print '(a, i0, a, es24.16, a, *(es24.16))', 'design iterate ', iteration, '  J = ', &
            & objective, '  c = ', c
       if (any(of_state) .and. .not. quasi) gradient = gradient - matmul(weight, a)
       if (quasi) hessian = quasi_hessian
       if (iteration == 1) tolerance = 1.0e-10_dp * max(1.0_dp, norm2(gradient))
       if (size(c) == 0 .and. .not. quasi) then
          if (verbosity >= 1) print '(a, i0, a, es12.4, a, *(es24.16))', 'design iteration ', iteration, &
               & '  |dJ/ds| = ', norm2(gradient), '  s = ', iterates(1) % points % expansion % at
          if (norm2(gradient) < tolerance) return
          do k = 1, size(iterates)
             iterates(k) % points % expansion % at = iterates(k) % points % expansion % at &
                  & - linear_solution(hessian, gradient)
          end do
          deallocate(c, a)
          cycle
       end if

       ! the iterate against the design before it: refused, the step is
       ! taken again from there, shorter
       refused = .false.
       if (with_accepted) then
          decrease = accepted_objective + dot_product(penalty, max(accepted_c, 0.0_dp)) &
               & - (objective + dot_product(penalty, max(c, 0.0_dp)))
          ! a decrease predicted below the rounding of the merit is not
          ! tested
          refused  = decrease < 0.1_dp * predicted .and. predicted > 1.0e-12_dp * max(1.0_dp, &
               & abs(accepted_objective) + dot_product(penalty, max(accepted_c, 0.0_dp)))
          if (refused) then
             if (verbosity >= 1) print '(a, es12.4, a, es12.4)', 'refused: the merit decreased by ', decrease, &
                  & ' where the model predicted ', predicted
             limit = norm2(step) / 2.0_dp
             do k = 1, size(iterates)
                iterates(k) % points % expansion % at = accepted
             end do
             gradient  = accepted_gradient
             hessian   = accepted_hessian
             c         = accepted_c
             a         = accepted_a
             objective = accepted_objective
             shift     = max(4.0_dp * shift, reach)
             reduction = reduction / 4.0_dp
          else if (decrease > 0.75_dp * predicted) then
             shift = shift / 4.0_dp
             if (shift < 1.0e-3_dp * reach) shift = 0.0_dp
             reduction = min(1.0_dp, 4.0_dp * reduction)
          end if
       end if
       if (.not. refused) then
          ! the update of the Lagrangian's hessian by the step accepted and
          ! the change of its gradient, at the multipliers of that step
          if (quasi .and. with_accepted) then
             call bfgs_update(quasi_hessian, iterates(1) % points % expansion % at - accepted, &
                  & gradient + matmul(transpose(a), mu) - accepted_gradient - matmul(transpose(accepted_a), mu))
             hessian = quasi_hessian
          end if
          accepted           = iterates(1) % points % expansion % at
          accepted_gradient  = gradient
          accepted_hessian   = hessian
          accepted_c         = c
          accepted_a         = a
          accepted_objective = objective
          with_accepted      = .true.
       end if

       ! the step: the hessian shifted while it is not positive definite
       reach = norm2(gradient) / (0.1_dp * max(1.0_dp, maxval(abs(accepted))))
       if (.not. positive_definite(hessian, shift)) shift = max(shift, reach)
       do while (.not. positive_definite(hessian, shift))
          shift = 2.0_dp * shift
       end do
       ! where no step meets the inequalities linearized, each violated
       ! one is required to remove a quarter of its fraction, three
       ! times, and then only not to grow, which the null step meets
       theta = reduction
       do quartering = 0, 4
          call quadratic_step(shifted(hessian, shift), gradient, c - (1.0_dp - theta) * max(c, 0.0_dp), a, step, &
               & mu, solvable)
          if (solvable) exit
          theta = merge(0.0_dp, theta / 4.0_dp, quartering == 3)
       end do
       if (.not. solvable) then
          error stop 'operation_residual: no step meets the inequalities on the design, not even the null step'
       end if
       if (refused .and. norm2(step) > limit) step = step * (limit / norm2(step))
       stationarity = norm2(gradient + matmul(transpose(a), mu))
       if (verbosity >= 1) then
          write(message, '(a, i0, a)') '(a, i0, a, es12.4, a, es12.4, a, ', n, 'es24.16, a, *(es24.16))'
          print trim(message), 'design iteration ', iteration, '  |dJ/ds + A^T mu| = ', stationarity, &
               & '  shift = ', shift, '  s = ', iterates(1) % points % expansion % at, '  mu = ', mu
       end if
       if (.not. refused .and. stationarity < tolerance .and. maxval(c) <= 1.0e-12_dp * max(1.0_dp, maxval(abs(c)))) then
          ! the multipliers of the design returned, not of the step:
          ! dJ/ds + A^T mu = 0 by least squares on the active inequalities
          active = pack([(j, j = 1, size(c))], mu > 0.0_dp)
          if (size(active) > 0) mu(active) = linear_solution(matmul(a(active, :), transpose(a(active, :))), &
               & -matmul(a(active, :), gradient))
          do k = 1, size(iterates)
             call iterates(k) % place_inequality_multipliers(mu(first(k):first(k + 1) - 1), solutions(k))
          end do
          return
       end if
       ! each penalty at least twice its multiplier, and the decrease of
       ! the merit the model predicts for the step
       penalty   = max(penalty, 2.0_dp * mu)
       predicted = -(dot_product(gradient, step) + dot_product(step, matmul(hessian, step)) / 2.0_dp) &
            & + dot_product(penalty, max(c, 0.0_dp) - max(c + matmul(a, step), 0.0_dp))
       weight = merge(mu, 0.0_dp, of_state)
       do k = 1, size(iterates)
          iterates(k) % points % expansion % at = accepted + step
       end do
       deallocate(c, a)
    end do
    write(message, '(a,i0,a)') 'operation_residual: Newton''s iteration on the free design did not reach &
         &dJ/ds = 0 in ', MAX_DESIGN_ITERATIONS, ' iterations'
    error stop trim(message)

  end subroutine design_iteration

  !===================================================================!
  ! THE INEQUALITIES ON THE DESIGN FACTOR at the design the manifold
  ! states: c(i) the value of the i-th equation, a(i, :) its gradient
  ! in the design coordinates; none where no inequality is stated. A
  ! functional of the state is evaluated on the solution, its gradient
  ! the expansion's column of order one along each design coordinate,
  ! or, given, the m-th column of gradients, -kappa of its own reverse
  ! sweep.
  !===================================================================!

  subroutine inequalities_at(this, solution, c, a, gradients)

    class(discrete_residual), intent(in)  :: this
    type(discrete_field)    , intent(in)  :: solution
    real(dp), allocatable   , intent(out) :: c(:), a(:,:)
    real(dp)                , intent(in), optional :: gradients(:,:)

    type(discrete_field) :: functional
    real(dp), allocatable :: slope(:)
    real(dp) :: given
    integer :: j, i, m, d

    m = 0
    do j = 1, size(this % condition)
       if (this % condition(j) % inequality) m = m + size(this % condition(j) % equation)
    end do
    allocate(c(m), a(m, this % points % expansion % num_designs()))
    m = 0
    do j = 1, size(this % condition)
       if (.not. this % condition(j) % inequality) cycle
       do i = 1, size(this % condition(j) % equation)
          m = m + 1
          if (this % condition(j) % equation(i) % is_integral()) then
             functional = this % condition(j) % equation(i) % at(solution)
             c(m) = functional % value(1, 1)
             if (present(gradients)) then
                a(m, :) = gradients(:, m)
                cycle
             end if
             do d = 1, size(a, 2)
                a(m, d) = functional % value(1, 1 + this % points % expansion % column(d, 1))
             end do
             cycle
          end if
          call design_condition(this % condition(j) % equation(i), this % points % expansion % values(), given, slope)
          c(m)    = -given
          a(m, :) = slope
       end do
    end do

  end subroutine inequalities_at

  ! which of the inequalities, in the order stated, are on functionals
  ! of the state
  function functional_inequalities(this) result(of_state)

    class(discrete_residual), intent(in) :: this
    logical, allocatable :: of_state(:)

    integer :: j, i

    allocate(of_state(0))
    do j = 1, size(this % condition)
       if (.not. this % condition(j) % inequality) cycle
       of_state = [of_state, (this % condition(j) % equation(i) % is_integral(), &
            & i = 1, size(this % condition(j) % equation))]
    end do

  end function functional_inequalities

  ! the multipliers of the inequalities into the solution, each a
  ! function on the design factor: its value at every point
  subroutine place_inequality_multipliers(this, mu, solution)

    class(discrete_residual), intent(in)    :: this
    real(dp)                , intent(in)    :: mu(:)
    type(discrete_field)    , intent(inout) :: solution

    integer :: j, i, m, column, first

    m = 0
    do j = 1, size(this % condition)
       if (.not. this % condition(j) % inequality) cycle
       first = 1
       do i = 1, size(this % condition(j) % equation)
          m = m + 1
          if (.not. this % condition(j) % paired) cycle
          ! the i-th column of the multiplier's name
          do column = first, size(solution % name)
             if (trim(solution % name(column)) == trim(this % condition(j) % multiplier % name)) exit
          end do
          first = column + 1
          solution % value(:, column) = mu(m)
          if (allocated(solution % jet) .and. allocated(solution % slot)) then
             solution % jet(solution % slot(column), :, 1) = mu(m)
          end if
       end do
    end do

  end subroutine place_inequality_multipliers

  !===================================================================!
  ! THE STEP OF THE DESIGN UNDER ITS INEQUALITIES: the minimizer of
  ! g . d + d . H d / 2 over the steps with c + A d <= 0, H positive
  ! definite, and the multiplier of each inequality, nonnegative and
  ! zero where it is not active; solvable false where no step meets
  ! the inequalities.
  !
  ! By the dual method of Goldfarb and Idnani: from the minimizer
  ! without inequalities, d = -H^-1 g, the most violated inequality p
  ! joins the active set by raising its multiplier t from zero, the
  ! active ones remaining equalities - the step moves by -t z and their
  ! multipliers by -t r, H z + N r = a_p and N^T z = 0, N the active
  ! normals - until p is met; where an active multiplier vanishes
  ! first, that inequality leaves and p goes on joining. The active
  ! normals stay independent: a_p depends on them where a_p . z, which
  ! is z . H z, is below 1e-6 of a_p . H^-1 a_p, and only the
  ! multipliers then move; and no step exists where p is not met and no
  ! active multiplier vanishes, (1, -r) on p and the active set then
  ! the weights of the inequalities whose normals cancel and whose
  ! constants sum to more than zero. At every move the step and the
  ! multipliers are those of one solve of the active set's equalities,
  ! not the sum of the moves. Of 2000 random problems: summed, and the
  ! dependence at 1e-12, an inequality nearly dependent raised the
  ! multipliers to 1e15 and the active values away from zero; at 1e-8,
  ! two nearly dependent in turn made five normals active in four
  ! dimensions; at 1e-6 every step meets its conditions to 1.3e-11, and
  ! each of the 230 without one has the weights to 1.4e-5. Where the
  ! primal active set it replaces settled, 1738 of them, the step and
  ! the multipliers are its own to every digit; it took dependent
  ! normals in, and the rotor's three scenarios at 2 cm, its eight
  ! upper bounds active and three smooth maxima violated, made its
  ! system singular.
  !===================================================================!

  subroutine quadratic_step(hessian, gradient, c, a, step, mu, solvable)

    real(dp), intent(in)  :: hessian(:,:), gradient(:), c(:), a(:,:)
    real(dp), allocatable, intent(out) :: step(:), mu(:)
    logical , intent(out) :: solvable

    real(dp), allocatable :: u(:), r(:), z(:)
    integer , allocatable :: set(:)
    logical :: active(size(c))
    real(dp) :: tolerance, t1, t2, joining, met, along
    integer :: j, p, leaving, iteration

    tolerance = 1.0e-12_dp * max(1.0_dp, maxval(abs(c)))
    allocate(mu(size(c)), source=0.0_dp)
    active   = .false.
    solvable = .true.
    do iteration = 1, 10 * (size(c) + size(gradient) + 1)
       set = pack([(j, j = 1, size(c))], active)
       call active_solution(set, 0, 0.0_dp, step, u)
       if (all(active)) exit
       p = maxloc(c + matmul(a, step), dim=1, mask=.not. active)
       if (c(p) + dot_product(a(p, :), step) <= tolerance) exit
       ! p joins: its multiplier raised from zero until it is met
       joining = 0.0_dp
       along   = dot_product(a(p, :), linear_solution(hessian, a(p, :)))
       do
          call active_solution(set, p, joining, step, u, z, r)
          met = c(p) + dot_product(a(p, :), step)
          ! the largest rise before an active multiplier vanishes, and the
          ! rise that meets p, none where a_p depends on the active normals
          t1      = huge(1.0_dp)
          leaving = 0
          do j = 1, size(set)
             if (r(j) > 0.0_dp) then
                if (u(j) / r(j) < t1) then
                   t1      = u(j) / r(j)
                   leaving = j
                end if
             end if
          end do
          t2 = huge(1.0_dp)
          if (dot_product(a(p, :), z) > 1.0e-6_dp * along) t2 = met / dot_product(a(p, :), z)
          if (leaving == 0 .and. t2 == huge(1.0_dp)) then
             solvable = .false.
             return
          end if
          joining = joining + min(t1, t2)
          if (t2 <= t1) then
             active(p) = .true.
             exit
          end if
          active(set(leaving)) = .false.
          set = pack([(j, j = 1, size(c))], active)
       end do
    end do
    if (size(c) > 0) then
       if (maxval(c + matmul(a, step)) > tolerance) then
          error stop 'operation_residual: the dual method on the inequalities of the design does not settle'
       end if
    end if
    mu(set) = max(u, 0.0_dp)

  contains

    ! the step and the multipliers of the active set with the multiplier
    ! joining of p, H d + g + N u + joining a_p = 0 and N^T d = -c_N, by
    ! one solve of their equalities; and given z and r, H z + N r = a_p
    ! and N^T z = 0
    subroutine active_solution(set, p, joining, step, u, z, r)
      integer , intent(in)  :: set(:), p
      real(dp), intent(in)  :: joining
      real(dp), allocatable, intent(out) :: step(:), u(:)
      real(dp), allocatable, intent(out), optional :: z(:), r(:)
      real(dp), allocatable :: system(:,:), x(:), right(:)
      integer :: n, i
      n = size(gradient)
      allocate(system(n + size(set), n + size(set)), source=0.0_dp)
      system(1:n, 1:n)    = hessian
      system(1:n, n + 1:) = transpose(a(set, :))
      system(n + 1:, 1:n) = a(set, :)
      right = -gradient
      if (p > 0) right = right - joining * a(p, :)
      x    = linear_solution(system, [right, -c(set)])
      step = x(1:n)
      u    = x(n + 1:)
      if (present(z)) then
         x = linear_solution(system, [a(p, :), (0.0_dp, i = 1, size(set))])
         z = x(1:n)
         r = x(n + 1:)
      end if
    end subroutine active_solution

  end subroutine quadratic_step

  ! whether h + shift I is positive definite, the symmetric part of h
  ! read: its factor by Cholesky exists
  pure logical function positive_definite(h, shift)
    real(dp), intent(in) :: h(:,:), shift
    real(dp) :: l(size(h, 1), size(h, 1)), s
    integer :: i, j
    l = 0.0_dp
    positive_definite = .false.
    do j = 1, size(h, 1)
       s = h(j, j) + shift - sum(l(j, 1:j - 1)**2)
       if (.not. s > 0.0_dp) return
       l(j, j) = sqrt(s)
       do i = j + 1, size(h, 1)
          l(i, j) = ((h(i, j) + h(j, i)) / 2.0_dp - sum(l(i, 1:j - 1) * l(j, 1:j - 1))) / l(j, j)
       end do
    end do
    positive_definite = .true.
  end function positive_definite

  !===================================================================!
  ! THE DAMPED BFGS UPDATE of b by the step s and the change y of the
  ! gradient, Powell's: where s.y < s.b s / 5, y is replaced by
  ! theta y + (1 - theta) b s, theta = 0.8 s.b s / (s.b s - s.y), so
  ! that b stays positive definite. A zero b is first the identity
  ! times y.y / s.y, and is left zero where s.y is not positive.
  !===================================================================!

  pure subroutine bfgs_update(b, s, y)
    real(dp), intent(inout) :: b(:,:)
    real(dp), intent(in)    :: s(:), y(:)
    real(dp) :: bs(size(s)), r(size(s)), sbs, sy, theta
    integer :: i
    sy = dot_product(s, y)
    if (all(b == 0.0_dp)) then
       if (.not. sy > 0.0_dp) return
       do i = 1, size(s)
          b(i, i) = dot_product(y, y) / sy
       end do
    end if
    bs  = matmul(b, s)
    sbs = dot_product(s, bs)
    if (.not. sbs > 0.0_dp) return
    theta = 1.0_dp
    if (sy < 0.2_dp * sbs) theta = 0.8_dp * sbs / (sbs - sy)
    r = theta * y + (1.0_dp - theta) * bs
    b = b - spread(bs, 2, size(s)) * spread(bs, 1, size(s)) / sbs &
         & + spread(r, 2, size(s)) * spread(r, 1, size(s)) / dot_product(s, r)
  end subroutine bfgs_update

  ! h + shift I
  pure function shifted(h, shift) result(s)
    real(dp), intent(in) :: h(:,:), shift
    real(dp) :: s(size(h, 1), size(h, 2))
    integer :: i
    s = h
    if (shift == 0.0_dp) return
    do i = 1, size(h, 1)
       s(i, i) = s(i, i) + shift
    end do
  end function shifted

  ! a design on which no equality is stated
  logical function free_design(this)
    class(discrete_residual), intent(in) :: this
    integer :: j
    free_design = this % manifold % with_design
    do j = 1, size(this % condition)
       if (this % condition(j) % manifold % is_design_factor() .and. .not. this % condition(j) % inequality) &
            & free_design = .false.
    end do
  end function free_design

  ! whether the expansions are along the design coordinates, one each
  ! in their order
  logical function along_each_coordinate(expanded)
    type(expansion), intent(in) :: expanded
    integer :: e, i
    along_each_coordinate = expanded % num_expansions() == expanded % num_designs()
    if (.not. along_each_coordinate) return
    do e = 1, expanded % num_expansions()
       do i = 1, expanded % num_designs()
          if (expanded % along(i, e) /= merge(1.0_dp, 0.0_dp, i == e)) along_each_coordinate = .false.
       end do
    end do
  end function along_each_coordinate

  !===================================================================!
  ! The solution x of A x = b by elimination with partial pivoting. A
  ! singular A stops the program.
  !===================================================================!

  function linear_solution(A, b) result(x)

    real(dp), intent(in) :: A(:,:), b(:)
    real(dp) :: x(size(b))

    real(dp) :: M(size(b), size(b) + 1), row(size(b) + 1)
    integer :: n, k, i, p

    n = size(b)
    M(:, 1:n)   = A
    M(:, n + 1) = b
    do k = 1, n
       p = k - 1 + maxloc(abs(M(k:n, k)), 1)
       if (M(p, k) == 0.0_dp) error stop 'operation_residual: the hessian of the free design is singular'
       if (p /= k) then
          row     = M(k, :)
          M(k, :) = M(p, :)
          M(p, :) = row
       end if
       do i = k + 1, n
          M(i, k:) = M(i, k:) - M(i, k) / M(k, k) * M(k, k:)
       end do
    end do
    do k = n, 1, -1
       x(k) = (M(k, n + 1) - dot_product(M(k, k + 1:n), x(k + 1:n))) / M(k, k)
    end do

  end function linear_solution

  !===================================================================!
  ! THE ZERO AT THE DESIGN THE IMAGE STATES, from an estimate: one
  ! block of the chain after another, each a residual_operator over
  ! the block's moments (its instants, and the stages of a one-step
  ! family) and the cells. The solution is a discrete field of the
  ! manifold's unknowns at the instants.
  !
  ! SCHEDULES: with checkpoint = k or snapshots = s the state is not
  ! stored over the chain. Both sweeps run over a window of instants,
  ! the blocks being solved with the history they read; the state at
  ! the first instant of a block - the history instants before it -
  ! is stored as a snapshot where the schedule places one, and the
  ! reverse sweep solves blocks again from the snapshots, from the
  ! same estimate and through the same iterates, to form their
  ! multipliers: with checkpoint = k, a snapshot every k blocks and
  ! each segment of k blocks solved once more, its stages retained
  ! over the segment; with snapshots = s, the binomial schedule of
  ! Griewank, which reverses n blocks with s snapshots in the least
  ! number of block solves, t n - C(s + t + 1, t - 1) + n for the
  ! least t with C(s + 1 + t, t) >= n, the stages retained over one
  ! block. The solution stores the unknowns and the multipliers at
  ! every point, and the objective's value in place of the quadrature
  ! nodes, so that it integrates the objective alone. A chain
  ! pipelined over the images retains every stage and takes no
  ! schedule; two images on a manifold of time alone divide a
  ! scheduled reverse sweep instead. With kappa present, the
  ! multipliers of the design, -dJ/ds in every column of the
  ! expansion, are returned in it.
  !===================================================================!

  recursive subroutine zero_at_design(this, estimate, solution, checkpoint, snapshots, kappa, functional_kappa)

    class(discrete_residual), intent(in)  :: this
    type(discrete_field)    , intent(in)  :: estimate
    type(discrete_field)    , intent(out) :: solution
    integer                 , intent(in), optional :: checkpoint, snapshots
    real(dp), allocatable   , intent(out), optional :: kappa(:,:), functional_kappa(:,:)

    real(dp), allocatable :: tuple(:,:,:,:), stages(:,:,:,:,:), value(:,:)
    real(dp), allocatable :: jet(:,:,:), adjoint(:,:,:,:), reaction(:,:,:), gauge(:,:,:)
    real(dp), allocatable :: design_multiplier(:,:), above(:,:,:,:), stages_above(:,:,:,:,:)
    integer , allocatable :: slot(:), stage_node(:,:,:)
    type(chain_block) :: point
    integer :: order, e, g, nm, col, ncol
    character(len=32), allocatable :: names(:)
    type(stencil), allocatable :: rows_along(:,:)
    integer :: stride, nf, ncells, ninst, f, c, k, b, first, last, depth, nb, top, o, npts, nc, j
    integer :: segment, shots, nrow
    real(dp), allocatable :: objective_value(:), field_value(:,:,:,:)
    type(expression_view) :: view
    type(discrete_residual) :: halves
    logical :: wavefront, scheduled
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

    ! the tuples at every instant and cell, in every column of the
    ! expansion along the design: in column zero the values from the
    ! estimate, every derivative component zero; zero in the columns
    ! above
    order = this % points % expansion % order
    ncol  = this % points % expansion % num_columns()
    segment = 0
    shots   = 0
    if (present(checkpoint)) segment = checkpoint
    if (present(snapshots))  shots   = snapshots
    if (segment < 0) error stop 'operation_residual: checkpoint is a number of blocks, at least one'
    if (shots < 0)   error stop 'operation_residual: snapshots is a number of stored states, at least one'
    if (segment > 0 .and. shots > 0) then
       error stop 'operation_residual: the states are stored every checkpoint blocks, or at snapshots placed by the &
            &binomial schedule; one of the two is given'
    end if
    scheduled = segment > 0 .or. shots > 0
    if (scheduled .and. allocated(this % functional_weight)) then
       if (any(this % functional_weight /= 0.0_dp)) then
          error stop 'operation_residual: an inequality on a functional of the state is solved without a checkpoint &
               &or snapshot schedule, which integrates the objective alone'
       end if
    end if
    if (scheduled .and. allocated(this % elements)) then
       error stop 'operation_residual: a law projected on the cells is not solved under a checkpoint or snapshot &
            &schedule'
    end if
    if (scheduled .and. .not. this % with_chain) then
       error stop 'operation_residual: checkpoints and snapshots divide a chain of blocks; the residual has none'
    end if
    ! with a schedule the wavefront, which retains every stage, is not
    ! formed; two images on a manifold of time alone divide the reverse
    ! sweep instead, image 2 solving the blocks again while image 1
    ! forms the multipliers of the blocks after them
    wavefront = this % pipelined() .and. .not. scheduled
    if (scheduled .and. num_images() == 2 .and. .not. this % manifold % with_space .and. .not. this % overlapped) then
       select type (this)
       type is (discrete_residual)
          halves = this
       end select
       halves % overlapped = .true.
       call halves % zero_at_design(estimate, solution, checkpoint, snapshots, kappa)
       return
    end if
    if (.not. scheduled) then
       allocate(tuple(stride, ncells, ninst, 0:ncol - 1), source=0.0_dp)
       do k = 1, ninst
          do c = 1, ncells
             do f = 1, nf
                tuple(this % rule % offset_of_field(f) + 1, c, k, 0) = estimate % value(this % points % point_of(k, c), f)
             end do
          end do
       end do
       if (this % continued) call this % estimated_rates(tuple(:, :, :, 0), 1)
    end if

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
          call this % block_extent(b, first, last, depth)
          if (first < 1) then
             write(message,'(a,i0,a,i0,a,i0)') 'operation_residual: block ', b, ' reads ', depth, &
                  & ' instants of history before instant ', first + depth
             error stop trim(message)
          end if
       end do
       if (.not. scheduled) then
          allocate(stages(stride, ncells, top, ninst, 0:ncol - 1), source=0.0_dp)
          ! the wavefront over the blocks and the orders, when the images
          ! number the orders plus one on a manifold of time alone with
          ! one expansion; the block loop otherwise
          if (wavefront) then
             call this % wavefront_expansion(rows_along, tuple, stages)
             nb = 0
          end if
          do b = 1, nb
             call this % block_extent(b, first, last, depth)
             call this % solve_block(this % schemes % scheme(b), first, last, depth, rows_along, origin, tuple, stages)
          end do
       else
          ! both sweeps by the schedule, the state in a window
          call this % scheduled_sweeps(estimate, rows_along, segment, shots, field_value, objective_value, &
               & adjoint, reaction, gauge, design_multiplier)
       end if
    else
       ! the point manifold: one system, its columns along the design
       allocate(stages(stride, 1, 1, 1, 0:ncol - 1), source=0.0_dp)
       call this % point_rows(tuple(:, :, :, 0), point)
       call this % solve_rows(point, origin, tuple, stages)
    end if
    ! the expansions shared among the images: their columns summed
    ! once, column zero being complete on every image; block by block
    ! within a schedule
    if (order > 0 .and. num_images() > 1 .and. .not. this % manifold % with_space .and. .not. scheduled) then
       above = tuple(:, :, :, 1:)
       call co_sum(above)
       tuple(:, :, :, 1:) = above
       stages_above = stages(:, :, :, :, 1:)
       call co_sum(stages_above)
       stages(:, :, :, :, 1:) = stages_above
    end if

    ! THE SOLUTION: the state at every point, with the jet the rule
    ! lays out; with the equations paired, the multipliers of every
    ! equality of L_h by the reverse sweep, those named returned -
    ! each equation's under the name its term declared, each paired
    ! condition's under its multiplier's name at the points of its
    ! part
    npts = this % points % num_points()
    nc = 0
    if (this % with_adjoint) then
       do j = 1, size(this % condition)
          associate (term => this % condition(j))
            if (term % paired .and. (term % manifold % is_time_face() .or. term % manifold % is_time_factor() &
                 & .or. term % manifold % is_space_face() .or. term % manifold % is_design_factor())) &
                 & nc = nc + size(term % equation)
          end associate
       end do
    end if
    ! the multiplier columns: the equations', then the conditions'
    nm = merge(nf + nc + this % num_represented(), 0, this % with_adjoint)
    allocate(value(npts, nf + nm), source=0.0_dp)
    ! the rows of the jet: the rule's components over the chain; under
    ! a schedule the unknowns alone, the state having been a window
    nrow = stride
    if (scheduled) nrow = nf
    allocate(names(size(value, 2)), jet(nrow + nm, npts, ncol))
    allocate(slot(size(value, 2)), source=0)
    jet = 0.0_dp
    names(1:nf) = this % manifold % unknown_name(1:nf)
    do f = 1, nf
       slot(f) = this % rule % offset_of_field(f) + 1
       if (scheduled) slot(f) = f
    end do
    do k = 1, ninst
       do c = 1, ncells
          if (scheduled) then
             value(this % points % point_of(k, c), 1:nf) = field_value(:, c, k, 0)
             jet(1:nf, this % points % point_of(k, c), :) = field_value(:, c, k, :)
             cycle
          end if
          do f = 1, nf
             value(this % points % point_of(k, c), f) = tuple(this % rule % offset_of_field(f) + 1, c, k, 0)
          end do
          jet(1:stride, this % points % point_of(k, c), :) = tuple(:, c, k, :)
       end do
    end do

    solution = discrete_field(this % points, value(:, 1:nf), names(1:nf))
    solution % jet  = jet
    solution % slot = slot(1:nf)
    solution % rule = this % rule
    if (.not. scheduled) then
       call this % quadrature_nodes(origin, tuple, stages, solution, stage_node)
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
    if (verbosity >= 1 .and. .not. scheduled) print '(a,f10.3,a,i0,a)', 'the forward sweep with its expansions took ', &
         & real(tock - tick, dp) / real(rate, dp), ' s; resident memory ', resident_kilobytes() / 1024, ' MB'

    if (this % with_adjoint) then
       if (.not. scheduled) then
          tick = tock
          call this % adjoint_sweep(rows_along, tuple, stages, solution, stage_node, adjoint, reaction, gauge, &
               & design_multiplier)
          call system_clock(tock)
          if (verbosity >= 1) print '(a,f10.3,a,i0,a)', 'the reverse sweep with its expansions took ', &
               & real(tock - tick, dp) / real(rate, dp), ' s; resident memory ', resident_kilobytes() / 1024, ' MB'
       end if
       ! every multiplier, in every column of the expansion, as a
       ! component with a row of the jet of its own
       names(nf + 1:2 * nf) = this % adjoint_name
       do k = 1, ninst
          do c = 1, ncells
             jet(nrow + 1:nrow + nf, this % points % point_of(k, c), :) = adjoint(1:nf, c, k, :)
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
                if (.not. this % condition(j) % paired) cycle
                nc = nc + 1
                names(2 * nf + nc) = this % condition(j) % multiplier % name
                ! an inequality's multiplier is the design iteration's,
                ! placed by minimize
                if (this % condition(j) % inequality) cycle
                jet(nrow + nf + nc, :, :) = spread(design_multiplier(f, :), 1, npts)
             end do
          end if
          if (this % condition(j) % manifold % is_time_factor() .or. this % condition(j) % manifold % is_space_face()) then
             ! the multiplier of a condition on the time factor or on a
             ! face of the region, a function of time alone: its value at
             ! every cell of the instant
             do f = 1, size(this % condition(j) % equation)
                g  = g + 1
                if (.not. this % condition(j) % paired) cycle
                nc = nc + 1
                names(2 * nf + nc) = this % condition(j) % multiplier % name
                do k = 1, ninst
                   do c = 1, ncells
                      jet(nrow + nf + nc, this % points % point_of(k, c), :) = gauge(k, g, :)
                   end do
                end do
             end do
          end if
          if (.not. this % condition(j) % manifold % is_time_face()) cycle
          k = this % face_instant(this % condition(j) % manifold % face_value)
          do f = 1, size(this % condition(j) % equation)
             e  = e + 1
             if (.not. this % condition(j) % paired) cycle
             nc = nc + 1
             names(2 * nf + nc) = this % condition(j) % multiplier % name
             do c = 1, ncells
                jet(nrow + nf + nc, this % points % point_of(k, c), :) = reaction(c, e, :)
             end do
          end do
       end do
       ! the multipliers of the representation laws, after the
       ! conditions', one column per equation at every point; equation e
       ! of the laws is entry nf + e of the adjoint
       e = 0
       do j = 1, size(this % representation)
          do f = 1, size(this % representation(j) % equation)
             nc = nc + 1
             e  = e + 1
             names(2 * nf + nc) = this % representation(j) % multiplier % name
             do k = 1, ninst
                do c = 1, ncells
                   jet(nrow + nf + nc, this % points % point_of(k, c), :) = adjoint(nf + e, c, k, :)
                end do
             end do
          end do
       end do
       do col = 1, nm
          slot(nf + col)     = nrow + col
          value(:, nf + col) = jet(nrow + col, :, 1)
       end do
    end if

    if (present(kappa) .and. allocated(design_multiplier)) then
       allocate(kappa(size(design_multiplier, 1), 0:ubound(design_multiplier, 2)))
       kappa = design_multiplier
    end if

    ! the multipliers of the design of each functional of the state an
    ! inequality states, by a reverse sweep of its own on the same
    ! forward sweep, the functional alone the objective: -dc/ds in
    ! column m, the m-th inequality's, zero for one on the design alone
    if (present(functional_kappa)) then
       if (scheduled) then
          error stop 'operation_residual: the derivatives of the functionals of the inequalities are formed &
               &without a checkpoint or snapshot schedule'
       end if
       separate: block
         type(discrete_residual) :: alone
         real(dp), allocatable :: adjoint_m(:,:,:,:), reaction_m(:,:,:), gauge_m(:,:,:), multiplier_m(:,:)
         logical , allocatable :: of_state(:)
         integer :: m
         select type (this)
         type is (discrete_residual)
            alone = this
         end select
         if (allocated(alone % objective)) deallocate(alone % objective)
         of_state = this % functional_inequalities()
         allocate(functional_kappa(size(this % points % expansion % values()), size(of_state)), source=0.0_dp)
         do m = 1, size(of_state)
            if (.not. of_state(m)) cycle
            alone % functional_weight    = [(0.0_dp, k = 1, size(of_state))]
            alone % functional_weight(m) = 1.0_dp
            call alone % adjoint_sweep(rows_along, tuple, stages, solution, stage_node, adjoint_m, reaction_m, &
                 & gauge_m, multiplier_m)
            functional_kappa(:, m) = multiplier_m(:, 0)
         end do
       end block separate
    end if

    estimate_nodes: block
      type(discrete_field) :: with_nodes
      with_nodes = solution
      solution = discrete_field(this % points, value, names)
      solution % jet  = jet
      solution % slot = slot
      solution % rule = this % rule
      if (.not. scheduled) then
         solution % node_jet      = with_nodes % node_jet
         if (allocated(with_nodes % node_rule)) solution % node_rule = with_nodes % node_rule
         solution % node_weight   = with_nodes % node_weight
         solution % node_position = with_nodes % node_position
         solution % node_point    = with_nodes % node_point
         solution % node_image    = with_nodes % node_image
      else
         solution % objective_value   = with_nodes % objective_value
         solution % objective_formula = with_nodes % objective_formula
      end if
    end block estimate_nodes

  end subroutine zero_at_design

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

  subroutine quadrature_nodes(this, lb, tuple, stages, solution, stage_node, blocks)

    class(discrete_residual), intent(in)    :: this
    type(lower_bounds)      , intent(in)    :: lb
    real(dp)                , intent(in)    :: tuple(:,:,lb % instant:,0:), stages(:,:,:,lb % stage:,0:)
    type(discrete_field)    , intent(inout) :: solution
    integer , allocatable   , intent(out)   :: stage_node(:,:,:)
    integer                 , intent(in), optional :: blocks(2)

    type(derivative_terms), allocatable :: dt(:), weight(:)
    integer , allocatable :: instant(:), cell_owner(:)
    real(dp), allocatable :: position(:)
    integer :: stride, ncells, ninst, npts, nb, b, first, last, depth, n, s, smax, k, i, c, j
    integer :: nnodes, nn, nf, b1, b2, klo, khi, offset
    logical :: staged
    real(dp) :: h

    if (allocated(this % elements)) then
       call this % projected_nodes(tuple, stages, solution, stage_node)
       return
    end if

    stride = this % rule % num_components()
    nf     = this % manifold % num_unknowns
    ncells = this % points % num_cells()
    ninst  = this % points % num_instants()
    npts   = this % points % num_points()

    ! the blocks whose nodes are formed: every block, or the range
    ! given, whose nodes are the points from the first instant of the
    ! window, lb % instant, to the range's last, the history instants
    ! before the range's own among them: a multistep rule weighs them,
    ! and they are nodes without a point, outside every face
    smax   = 1
    nnodes = npts
    offset = lb % node - 1
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
          call this % block_extent(b1, first, last, depth)
          klo = first + depth
          call this % block_extent(b2, first, last, depth)
          khi = last
          nnodes = (khi - lb % instant + 1) * ncells
          if (lb % node /= this % points % point_of(lb % instant, 1)) then
             error stop 'operation_residual: the first node of a window is the first point of its first instant'
          end if
       end if
       do b = 1, nb
          smax = max(smax, this % schemes % scheme(b) % num_stages())
          if (b < b1 .or. b > b2) cycle
          call this % block_extent(b, first, last, depth)
          if (marches_by_stages(this % schemes % scheme(b), top_degree(this % rule, nf) + 1)) then
             nnodes = nnodes + (last - first) * this % schemes % scheme(b) % num_stages() * ncells
          end if
       end do
    end if
    allocate(stage_node(lb % stage:ubound(stages, 4), smax, ncells), source=0)
    if (allocated(solution % node_jet)) then
       deallocate(solution % node_jet, solution % node_weight, solution % node_position, solution % node_point, &
            & solution % node_image)
    end if
    allocate(solution % node_jet(stride, nnodes, this % points % expansion % num_columns()), source=0.0_dp)
    allocate(solution % node_weight(nnodes), source=0.0_dp)
    allocate(solution % node_position(this % points % num_coordinates(), nnodes), source=0.0_dp)
    allocate(solution % node_point(nnodes), source=0)
    allocate(solution % node_image(nnodes), source=this_image())
    cell_owner = this % cell_owners()

    ! the points as nodes, their jets those of the solution
    do k = lb % instant, khi
       do c = 1, ncells
          nn = this % points % point_of(k, c) - offset
          if (k >= klo) solution % node_point(nn) = nn + offset
          solution % node_image(nn)       = cell_owner(c)
          solution % node_position(:, nn) = this % points % position(nn)
          solution % node_jet(:, nn, :) = tuple(:, c, k, :)
       end do
    end do

    if (.not. this % with_chain) then
       solution % node_weight = 1.0_dp
       return
    end if

    nn = (khi - lb % instant + 1) * ncells
    do b = b1, b2
       call this % block_extent(b, first, last, depth)
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
                     stage_node(first + k, i, c) = nn + offset
                     solution % node_image(nn)   = cell_owner(c)
                     solution % node_weight(nn)  = h * scheme % stage_weight(i) * cell_volume(c)
                     position = this % points % position(this % points % point_of(first + k - 1, c))
                     position(1) = position(1) + scheme % stage_abscissa(i) * h
                     solution % node_position(:, nn) = position
                     solution % node_jet(:, nn, :) = stages(:, c, i, first + k, :)
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
                     solution % node_weight(this % points % point_of(first + instant(j) - 1, c) - offset) = &
                          & solution % node_weight(this % points % point_of(first + instant(j) - 1, c) - offset) &
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
      if (this % manifold % with_space) cell_volume = this % points % volume(c)
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
  ! BOTH SWEEPS BY A SCHEDULE, the state in a window of instants. The
  ! window of the blocks b1..b2 runs from reach instants before the
  ! first own instant of b1, reach the deepest history of the chain,
  ! to the last instant of b2; a snapshot of the boundary m is the
  ! state at the reach instants before the first own instant of the
  ! block m. The forward sweep solves block after block, stores the
  ! unknowns at every point, sums the objective over each block's
  ! quadrature nodes and stores the snapshots the schedule places.
  ! The reverse sweep performs the schedule's actions: a snapshot
  ! restored, blocks solved again from it, a snapshot stored, and the
  ! blocks b1..b2 reversed - solved again with their stages, the
  ! objective's differential formed at their nodes, their multipliers
  ! formed from the last block to the first. The multipliers passed
  ! to earlier instants are stored over the same window. A multistep
  ! rule of a block weighs history instants, whose differential is
  ! passed back with the multipliers to the blocks that own them.
  !
  ! With this % overlapped, image 2 performs the solves and image 1
  ! the reversals: image 2 writes the window of each reversal into the
  ! coarray half of its parity and meets image 1, which reads it, so
  ! that image 2 is one reversal ahead.
  !===================================================================!

  subroutine scheduled_sweeps(this, estimate, rows_along, segment, shots, field_value, objective_value, &
       & adjoint, reaction, gauge, design_multiplier)

    class(discrete_residual), intent(in)  :: this
    type(discrete_field)    , intent(in)  :: estimate
    type(stencil)           , intent(in)  :: rows_along(:,:)
    integer                 , intent(in)  :: segment, shots
    real(dp), allocatable   , intent(out) :: field_value(:,:,:,:), objective_value(:)
    real(dp), allocatable   , intent(out) :: adjoint(:,:,:,:), reaction(:,:,:), gauge(:,:,:)
    real(dp), allocatable   , intent(out) :: design_multiplier(:,:)

    real(dp), allocatable :: state(:,:,:,:), stages(:,:,:,:,:)
    real(dp), allocatable :: incoming(:,:,:,:), moved(:,:,:,:)
    real(dp), allocatable :: shot(:,:,:,:,:)
    real(dp), allocatable :: state_co(:,:,:,:,:)[:], stages_co(:,:,:,:,:,:)[:]
    real(dp), allocatable :: partial(:,:), part_value(:), part_partial(:,:), gradient(:,:,:)
    integer , allocatable :: taken(:), action(:,:), stage_node(:,:,:)
    type(discrete_field) :: part
    type(lower_bounds) :: lb
    integer :: stride, nf, ncells, ninst, order, ncol, ni, nb, top, reach, b, first, last, depth, nslots, again
    integer :: ia, k, c, f, b1, b2, wlo, whi, lo, width, half, reversals
    integer(8) :: optimum, tick, tock, rate
    logical :: solves, reverses

    stride = this % rule % num_components()
    nf     = this % manifold % num_unknowns
    ncells = this % points % num_cells()
    ninst  = this % points % num_instants()
    order  = this % points % expansion % order
    ncol   = this % points % expansion % num_columns()
    ni     = size(this % points % expansion % values())
    nb     = this % schemes % num_blocks()
    solves   = .not. this % overlapped .or. this_image() == 2
    reverses = .not. this % overlapped .or. this_image() == 1

    top   = 1
    reach = 0
    do b = 1, nb
       top = max(top, this % schemes % scheme(b) % num_stages())
       call this % block_extent(b, first, last, depth)
       reach = max(reach, depth)
    end do

    ! the field of the nodes of a window: the manifold and the rule
    allocate(part % on, source=this % points)
    part % rule = this % rule

    call reversal_schedule(nb, segment, shots, taken, action, again, optimum)
    nslots = max(maxval(taken), 0)
    do ia = 1, size(action, 2)
       if (action(1, ia) == ACTION_STORE) nslots = max(nslots, action(4, ia))
    end do
    if (verbosity >= 1) then
       print '(a,i0,a,i0,a,i0,a,i0)', 'schedule: blocks ', nb, '  snapshots ', nslots, '  block solves forward ', nb, &
            & '  again ', again
       if (shots > 0) print '(a,i0)', 'schedule: the binomial count of block solves, the forward sweep among them, ', optimum
    end if
    allocate(shot(stride, ncells, reach, 0:ncol - 1, nslots), source=0.0_dp)

    ! THE FORWARD SWEEP
    call system_clock(tick, rate)
    allocate(field_value(nf, ncells, ninst, 0:ncol - 1), source=0.0_dp)
    allocate(objective_value(ncol), partial(ni, 0:ncol - 1), source=0.0_dp)
    allocate(state(stride, ncells, 1:0, 0:ncol - 1))
    do b = 1, nb
       if (taken(b) > 0) call stored(b, taken(b))
       call window_moved(b, b)
       call blocks_solved(b, b)
       call this % block_extent(b, first, last, depth)
       do k = first + depth, last
          do c = 1, ncells
             do f = 1, nf
                field_value(f, c, k, :) = state(this % rule % offset_of_field(f) + 1, c, k, :)
             end do
          end do
       end do
       if (allocated(this % objective)) then
          call this % quadrature_nodes(lb, state, stages, part, stage_node, [b, b])
          call this % objective % functional(part, part_value, gradient, part_partial)
          objective_value = objective_value + part_value
          partial         = partial + part_partial
       end if
       deallocate(stages)
    end do
    call system_clock(tock)
    if (verbosity >= 1) print '(a,f10.3,a,i0,a)', 'the forward sweep with its expansions took ', &
         & real(tock - tick, dp) / real(rate, dp), ' s; resident memory ', resident_kilobytes() / 1024, ' MB'
    if (.not. this % with_adjoint) return

    ! THE REVERSE SWEEP
    tick = tock
    call this % allocate_multipliers(partial, 0, adjoint, incoming, reaction, gauge, design_multiplier)

    if (this % overlapped) then
       ! the widest window among the reversals, and the two halves
       width = 0
       do ia = 1, size(action, 2)
          if (action(1, ia) /= ACTION_REVERSE) cycle
          call this % block_extent(action(3, ia), first, last, depth)
          width = max(width, last - max(1, this % schemes % from(action(2, ia)) - reach) + 1)
       end do
       allocate(state_co(stride, ncells, width, 0:ncol - 1, 2)[*], stages_co(stride, ncells, top, width, 0:ncol - 1, 2)[*])
       state_co  = 0.0_dp
       stages_co = 0.0_dp
    end if

    reversals = 0
    do ia = 1, size(action, 2)
       select case (action(1, ia))
       case (ACTION_RESTORE)
          if (solves) call restored(action(2, ia), action(4, ia))
       case (ACTION_ADVANCE)
          if (solves) then
             do b = action(2, ia), action(3, ia) - 1
                call window_moved(b, b)
                call blocks_solved(b, b)
                deallocate(stages)
             end do
          end if
       case (ACTION_STORE)
          if (solves) call stored(action(2, ia), action(4, ia))
       case (ACTION_REVERSE)
          b1 = action(2, ia)
          b2 = action(3, ia)
          reversals = reversals + 1
          half = mod(reversals, 2) + 1
          if (solves) then
             call window_moved(b1, b2)
             call blocks_solved(b1, b2)
          end if
          if (this % overlapped) then
             call this % block_extent(b2, first, last, depth)
             whi = last
             wlo = max(1, this % schemes % from(b1) - reach)
             call this % block_extent(b1, first, last, depth)
             lo  = first + 1
             if (this_image() == 2) then
                state_co(:, :, 1:whi - wlo + 1, :, half) = state(:, :, wlo:whi, :)
                if (whi >= lo) stages_co(:, :, :, 1:whi - lo + 1, :, half) = stages(:, :, :, lo:whi, :)
                sync images(1)
                deallocate(stages)
                cycle
             end if
             sync images(2)
             if (allocated(state)) deallocate(state)
             allocate(state(stride, ncells, wlo:whi, 0:ncol - 1), source=0.0_dp)
             allocate(stages(stride, ncells, top, lo:whi, 0:ncol - 1), source=0.0_dp)
             ! sections on the left: an allocatable array assigned whole
             ! from another image may be allocated again
             state(:, :, wlo:whi, :) = state_co(:, :, 1:whi - wlo + 1, :, half)[2]
             if (whi >= lo) stages(:, :, :, lo:whi, :) = stages_co(:, :, :, 1:whi - lo + 1, :, half)[2]
             lb = lower_bounds(instant=wlo, stage=lo, node=this % points % point_of(wlo, 1))
          end if
          if (verbosity >= 1) print '(a,i0,a,i0,a,i0,a)', 'reversal of the blocks ', b1, '..', b2, &
               & '  resident memory ', resident_kilobytes() / 1024, ' MB'
          call reversed(b1, b2)
          deallocate(stages)
       end select
    end do
    if (this % overlapped) then
       ! the multipliers, formed by image 1, given to image 2
       sync all
       call co_broadcast(adjoint, source_image=1)
       call co_broadcast(reaction, source_image=1)
       call co_broadcast(gauge, source_image=1)
       call co_broadcast(design_multiplier, source_image=1)
    end if
    call this % shared_expansions_summed(adjoint, reaction, gauge, design_multiplier)
    call system_clock(tock)
    if (verbosity >= 1) print '(a,f10.3,a,i0,a)', 'the reverse sweep with its expansions took ', &
         & real(tock - tick, dp) / real(rate, dp), ' s; resident memory ', resident_kilobytes() / 1024, ' MB'

  contains

    ! the window moved to the blocks b1..b2: the history before b1 from
    ! the state, the own instants from the estimate, their columns
    ! above zero
    subroutine window_moved(b1, b2)
      integer, intent(in) :: b1, b2
      real(dp), allocatable :: window(:,:,:,:)
      integer :: own, first, last, depth, k, c, f
      own = this % schemes % from(b1)
      call this % block_extent(b2, first, last, depth)
      wlo = max(1, own - reach)
      whi = last
      allocate(window(stride, ncells, wlo:whi, 0:ncol - 1), source=0.0_dp)
      if (own > wlo) then
         if (lbound(state, 3) > wlo .or. ubound(state, 3) < own - 1) then
            error stop 'operation_residual: the state before a block of the schedule is not in the window'
         end if
         window(:, :, wlo:own - 1, :) = state(:, :, wlo:own - 1, :)
      end if
      do k = own, whi
         do c = 1, ncells
            do f = 1, nf
               window(this % rule % offset_of_field(f) + 1, c, k, 0) = estimate % value(this % points % point_of(k, c), f)
            end do
         end do
      end do
      if (this % continued .and. whi - own >= 2) call this % estimated_rates(window(:, :, own:whi, 0), own)
      call move_alloc(window, state)
    end subroutine window_moved

    ! the blocks b1..b2 solved in the window, their stages retained;
    ! the columns of expansions the images share summed once
    subroutine blocks_solved(b1, b2)
      integer, intent(in) :: b1, b2
      real(dp), allocatable :: own(:,:,:,:), stages_own(:,:,:,:,:)
      integer :: b, first, last, depth, from
      call this % block_extent(b1, first, last, depth)
      lo = first + 1
      lb = lower_bounds(instant=wlo, stage=lo, node=this % points % point_of(wlo, 1))
      allocate(stages(stride, ncells, top, lo:whi, 0:ncol - 1), source=0.0_dp)
      do b = b1, b2
         call this % block_extent(b, first, last, depth)
         call this % solve_block(this % schemes % scheme(b), first, last, depth, rows_along, lb, state, stages)
      end do
      if (order > 0 .and. num_images() > 1 .and. .not. this % manifold % with_space .and. .not. this % overlapped) then
         from = this % schemes % from(b1)
         own  = state(:, :, from:whi, 1:)
         call co_sum(own)
         state(:, :, from:whi, 1:) = own
         stages_own = stages(:, :, :, :, 1:)
         call co_sum(stages_own)
         stages(:, :, :, :, 1:) = stages_own
      end if
    end subroutine blocks_solved

    ! the snapshot of the boundary m: the state at the reach instants
    ! before the first own instant of the block m
    subroutine stored(m, place)
      integer, intent(in) :: m, place
      integer :: own, from
      own  = this % schemes % from(m)
      from = max(1, own - reach)
      if (own == from) return
      shot(:, :, 1:own - from, :, place) = state(:, :, from:own - 1, :)
    end subroutine stored

    subroutine restored(m, place)
      integer, intent(in) :: m, place
      integer :: own, from
      own  = this % schemes % from(m)
      from = max(1, own - reach)
      if (allocated(state)) deallocate(state)
      allocate(state(stride, ncells, from:own - 1, 0:ncol - 1))
      if (own == from) return
      state = shot(:, :, 1:own - from, :, place)
    end subroutine restored

    ! the multipliers of the blocks b1..b2, solved in the window with
    ! their stages: the objective's differential at their nodes, its
    ! part at the history instants passed back with the multipliers
    subroutine reversed(b1, b2)
      integer, intent(in) :: b1, b2
      integer :: b, first, last, depth, own, k, c, low, high, node
      own = this % schemes % from(b1)
      call this % quadrature_nodes(lb, state, stages, part, stage_node, [b1, b2])
      if (allocated(this % objective)) then
         call this % objective % functional(part, part_value, gradient, part_partial)
      else
         if (allocated(gradient)) deallocate(gradient)
         allocate(gradient(stride, size(part % node_weight), 0:ncol - 1), source=0.0_dp)
      end if
      ! the multipliers passed in, over the window
      allocate(moved(stride, ncells, lb % instant:ubound(state, 3), 0:ncol - 1), source=0.0_dp)
      low  = max(lbound(moved, 3), lbound(incoming, 3))
      high = min(ubound(moved, 3), ubound(incoming, 3))
      if (high >= low) moved(:, :, low:high, :) = incoming(:, :, low:high, :)
      call move_alloc(moved, incoming)
      do k = lb % instant, own - 1
         do c = 1, ncells
            node = this % points % point_of(k, c) - lb % node + 1
            incoming(:, c, k, :) = incoming(:, c, k, :) - gradient(:, node, :)
         end do
      end do
      do b = b2, b1, -1
         call this % block_extent(b, first, last, depth)
         call this % adjoint_block(this % schemes % scheme(b), first, last, depth, rows_along, lb, state, stages, &
              & gradient, stage_node, adjoint, incoming, reaction, gauge, design_multiplier)
      end do
    end subroutine reversed

  end subroutine scheduled_sweeps

  !===================================================================!
  ! THE SCHEDULE of a reverse sweep over nb blocks. taken(m) is the
  ! place of the snapshot the forward sweep stores at the boundary m,
  ! zero where it stores none; action(:, i) is the i-th action of the
  ! reverse sweep, (kind, a, b, place): the snapshot of the boundary a
  ! restored from the place (place zero the beginning, which stores
  ! nothing); the blocks a..b-1 solved again; the snapshot of the
  ! boundary a stored at the place; the blocks a..b reversed.
  !
  ! With segment = k: a snapshot at the first block of every segment
  ! of k blocks, each segment restored and reversed whole.
  !
  ! With shots = s: the binomial schedule. cost(s, n), the least
  ! number of block solves that reverse n blocks from a stored state
  ! with s free places, one block being reversed by one solve, is
  !
  !      cost(s, 1) = 1,     cost(0, n) = n (n + 1) / 2,
  !      cost(s, n) = min over l of  l + cost(s - 1, n - l) + cost(s, l),
  !
  ! l blocks solved to the place of a snapshot, the blocks after it
  ! reversed with one place fewer, the blocks before it with all; the
  ! least l that attains the minimum does not decrease with n and
  ! grows by one at most, so that the table is formed from two
  ! candidates per entry. cost(s, n) = t n - C(s + t + 1, t - 1) + n
  ! for the least t with C(s + 1 + t, t) >= n, Griewank's count. The
  ! forward sweep is the first descent: it solves every block once and
  ! stores the snapshots of that descent.
  !===================================================================!

  subroutine reversal_schedule(nb, segment, shots, taken, action, again, optimum)

    integer             , intent(in)  :: nb, segment, shots
    integer, allocatable, intent(out) :: taken(:), action(:,:)
    integer             , intent(out) :: again
    integer(8)          , intent(out) :: optimum

    integer(8), allocatable :: cost(:,:)
    integer   , allocatable :: split(:,:), grown(:,:)
    integer(8) :: candidate
    integer :: count, s, n, l, sg, nseg, b1, b2, i

    allocate(taken(nb + 1), source=0)
    allocate(action(4, 64), source=0)
    count   = 0
    optimum = 0

    if (segment > 0) then
       nseg = (nb + segment - 1) / segment
       do sg = 2, nseg
          taken((sg - 1) * segment + 1) = sg - 1
       end do
       do sg = nseg, 1, -1
          b1 = (sg - 1) * segment + 1
          b2 = min(nb, sg * segment)
          call emitted(ACTION_RESTORE, b1, 0, sg - 1)
          call emitted(ACTION_REVERSE, b1, b2, 0)
       end do
    else
       allocate(cost(0:shots, nb), split(0:shots, nb))
       split = 0
       do n = 1, nb
          cost(0, n) = int(n, 8) * int(n + 1, 8) / 2
       end do
       do s = 1, shots
          cost(s, 1) = 1
          do n = 2, nb
             cost(s, n) = huge(cost)
             do l = max(1, split(s, n - 1)), min(n - 1, split(s, n - 1) + 1)
                candidate = int(l, 8) + cost(s - 1, n - l) + cost(s, l)
                if (candidate < cost(s, n)) then
                   cost(s, n)  = candidate
                   split(s, n) = l
                end if
             end do
          end do
       end do
       optimum = cost(shots, nb)
       call reversal(1, nb + 1, shots, 0, .true.)
    end if

    grown  = action(:, 1:count)
    call move_alloc(grown, action)
    again = 0
    do i = 1, count
       if (action(1, i) == ACTION_ADVANCE) again = again + action(3, i) - action(2, i)
       if (action(1, i) == ACTION_REVERSE) again = again + action(3, i) - action(2, i) + 1
    end do

  contains

    subroutine emitted(kind, a, b, place)
      integer, intent(in) :: kind, a, b, place
      if (count == size(action, 2)) then
         allocate(grown(4, 2 * count), source=0)
         grown(:, 1:count) = action
         call move_alloc(grown, action)
      end if
      count = count + 1
      action(:, count) = [kind, a, b, place]
    end subroutine emitted

    ! the blocks a..b-1 reversed from the state of the boundary a,
    ! stored at place_a, with s free places; in the first descent
    ! the solves and the snapshots are the forward sweep's
    recursive subroutine reversal(a, b, s, place_a, first_descent)
      integer, intent(in) :: a, b, s, place_a
      logical, intent(in) :: first_descent
      integer :: n, m, j, place_m
      n = b - a
      if (n == 1) then
         call emitted(ACTION_RESTORE, a, 0, place_a)
         call emitted(ACTION_REVERSE, a, a, 0)
         return
      end if
      if (s == 0) then
         do j = b - 1, a, -1
            call emitted(ACTION_RESTORE, a, 0, place_a)
            if (j > a) call emitted(ACTION_ADVANCE, a, j, 0)
            call emitted(ACTION_REVERSE, j, j, 0)
         end do
         return
      end if
      m       = a + split(s, n)
      place_m = shots - s + 1
      if (first_descent) then
         taken(m) = place_m
      else
         call emitted(ACTION_RESTORE, a, 0, place_a)
         call emitted(ACTION_ADVANCE, a, m, 0)
         call emitted(ACTION_STORE, m, 0, place_m)
      end if
      call reversal(m, b, s - 1, place_m, first_descent)
      call reversal(a, m, s, place_a, .false.)
    end subroutine reversal

  end subroutine reversal_schedule

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
  ! THE REVERSE SWEEP: the multipliers of every row of L_h in every
  ! column of the expansion, block by block from the last, each
  ! block's transposed jacobian at its solution once per column. The
  ! differential of the objective in the jet is the source in column
  ! zero; in the column of order o along an expansion, the
  ! derivatives along it of that differential and of the product of
  ! the multipliers with the equations' jacobian, the multiplier's own
  ! derivative of that order left out. In the column of the order o
  ! along the expansion d, adjoint(f, c, k, :) is the o-th derivative
  ! along d of the multiplier of equation f at cell c and instant k,
  ! reaction(c, e, :) of the e-th face condition equation at cell c,
  ! gauge(k, g, :) of the g-th time-factor condition equation at
  ! instant k, and design_multiplier(i, :) of the multiplier of the
  ! design condition on the coordinate i, kappa_i = -dJ/dnu_i, which
  ! is -d^(o+1) J / (dnu_i dnu_d^o) there. A residual without a chain
  ! has no blocks to sweep, and is refused.
  !===================================================================!

  subroutine adjoint_sweep(this, rows_along, tuple, stages, solution, stage_node, adjoint, reaction, gauge, &
       & design_multiplier)

    class(discrete_residual), intent(in)  :: this
    type(stencil)           , intent(in)  :: rows_along(:,:)
    real(dp)                , intent(in)  :: tuple(:,:,:,0:), stages(:,:,:,:,0:)
    type(discrete_field)    , intent(in)  :: solution
    integer                 , intent(in)  :: stage_node(:,:,:)
    real(dp), allocatable   , intent(out) :: adjoint(:,:,:,:), reaction(:,:,:), gauge(:,:,:)
    real(dp), allocatable   , intent(out) :: design_multiplier(:,:)

    real(dp), allocatable :: gradient(:,:,:), partial(:,:), incoming(:,:,:,:), g(:,:,:), p(:,:)
    type(chain_block) :: point
    integer :: stride, ninst, b, nb, first, last, depth, ncol, ni, j, i, m

    stride = this % rule % num_components()
    ninst  = this % points % num_instants()
    ncol   = this % points % expansion % num_columns()
    ni     = size(this % points % expansion % values())

    if (allocated(this % objective)) then
       call this % objective % differential(solution, gradient, partial)
    else
       allocate(gradient(size(solution % node_jet, 1), size(solution % node_weight), 0:ncol - 1), &
            & partial(ni, 0:ncol - 1), source=0.0_dp)
    end if
    ! the functionals of the inequalities, each times its multiplier:
    ! the differential of the Lagrangian's functionals
    if (allocated(this % functional_weight)) then
       m = 0
       do j = 1, size(this % condition)
          if (.not. this % condition(j) % inequality) cycle
          do i = 1, size(this % condition(j) % equation)
             m = m + 1
             if (this % functional_weight(m) == 0.0_dp) cycle
             call this % condition(j) % equation(i) % differential(solution, g, p)
             gradient = gradient + this % functional_weight(m) * g
             partial  = partial  + this % functional_weight(m) * p
          end do
       end do
    end if

    call this % allocate_multipliers(partial, ninst, adjoint, incoming, reaction, gauge, design_multiplier)

    if (.not. this % with_chain) then
       ! the point manifold: one system, its transposed solve
       call this % point_rows(tuple(:, :, :, 0), point)
       call this % adjoint_rows(point, origin, tuple, stages, gradient, stage_node, adjoint, incoming, reaction, &
            & gauge, design_multiplier)
       return
    end if

    if (this % pipelined()) then
       ! the wavefront in reverse: image o + 1 forms order o of every
       ! block, one block behind image o; each image's multipliers of
       ! its own order, the column o of the one expansion, summed over
       ! the images at the end
       design_multiplier = 0.0_dp
       design_multiplier(:, this_image() - 1) = -partial(:, this_image() - 1)
       call this % wavefront_adjoint(rows_along, tuple, stages, gradient, stage_node, &
            & adjoint, incoming, reaction, gauge, design_multiplier)
       call co_sum(adjoint)
       call co_sum(reaction)
       call co_sum(gauge)
       call co_sum(design_multiplier)
       return
    end if

    nb = this % schemes % num_blocks()
    do b = nb, 1, -1
       call this % block_extent(b, first, last, depth)
       call this % adjoint_block(this % schemes % scheme(b), first, last, depth, rows_along, origin, tuple, &
            & stages, gradient, stage_node, adjoint, incoming, reaction, gauge, design_multiplier)
    end do
    call this % shared_expansions_summed(adjoint, reaction, gauge, design_multiplier)

  end subroutine adjoint_sweep

  !===================================================================!
  ! The multipliers of the reverse sweep, zero, in the columns of the
  ! partial given, with the incoming multipliers over ninst instants.
  ! THE MULTIPLIERS OF THE DESIGN CONDITIONS nu_i - nu_i0 = 0 and
  ! their derivatives: the stationarity of L in nu_i, kappa_i =
  ! -(dJ/dnu_i at fixed state + the sum over every row of its
  ! multiplier times the row's partial in nu_i) = -dJ/dnu_i, the total
  ! derivative by the adjoint; the partial in nu_i here, the rows'
  ! sums by the blocks. An expansion another image owns: its columns
  ! are left zero here, and summed over the images at the end.
  !===================================================================!

  subroutine allocate_multipliers(this, partial, ninst, adjoint, incoming, reaction, gauge, design_multiplier)

    class(discrete_residual), intent(in)  :: this
    real(dp)                , intent(in)  :: partial(:,0:)
    integer                 , intent(in)  :: ninst
    real(dp), allocatable   , intent(out) :: adjoint(:,:,:,:), incoming(:,:,:,:), reaction(:,:,:)
    real(dp), allocatable   , intent(out) :: gauge(:,:,:), design_multiplier(:,:)

    integer :: nf, ncells, top, nc, ng, j, d

    nf     = this % manifold % num_unknowns
    ncells = this % points % num_cells()
    top    = ubound(partial, 2)
    nc = 0
    ng = 0
    do j = 1, size(this % condition)
       if (this % condition(j) % manifold % is_time_face())        nc = nc + size(this % condition(j) % equation)
       ! the multipliers that are functions of time alone: of the
       ! conditions on the time factor and on the faces of the region
       if (this % condition(j) % manifold % is_time_factor() .or. this % condition(j) % manifold % is_space_face()) &
            & ng = ng + size(this % condition(j) % equation)
    end do
    allocate(adjoint(nf + this % num_represented(), ncells, this % points % num_instants(), 0:top), source=0.0_dp)
    allocate(incoming(this % rule % num_components(), ncells, ninst, 0:top), source=0.0_dp)
    allocate(reaction(ncells, nc, 0:top), gauge(this % points % num_instants(), ng, 0:top), source=0.0_dp)
    allocate(design_multiplier(size(partial, 1), 0:top), source=0.0_dp)
    design_multiplier = -partial
    associate (expanded => this % points % expansion)
      do d = 1, expanded % num_expansions()
         if (.not. this % expansion_owned(d)) then
            design_multiplier(:, expanded % column(d, 1):expanded % column(d, expanded % order)) = 0.0_dp
         end if
      end do
    end associate

  end subroutine allocate_multipliers

  ! the columns of the expansions the images shared, summed once:
  ! column zero is complete on every image and is left alone
  subroutine shared_expansions_summed(this, adjoint, reaction, gauge, design_multiplier)

    class(discrete_residual), intent(in)    :: this
    real(dp)                , intent(inout) :: adjoint(:,:,:,0:), reaction(:,:,0:), gauge(:,:,0:)
    real(dp)                , intent(inout) :: design_multiplier(:,0:)

    real(dp), allocatable :: above4(:,:,:,:), above3(:,:,:), above2(:,:)

    if (ubound(adjoint, 4) < 1 .or. num_images() == 1 .or. this % manifold % with_space .or. this % overlapped) return
    above4 = adjoint(:, :, :, 1:)
    call co_sum(above4)
    adjoint(:, :, :, 1:) = above4
    above3 = reaction(:, :, 1:)
    call co_sum(above3)
    reaction(:, :, 1:) = above3
    above3 = gauge(:, :, 1:)
    call co_sum(above3)
    gauge(:, :, 1:) = above3
    above2 = design_multiplier(:, 1:)
    call co_sum(above2)
    design_multiplier(:, 1:) = above2

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
  ! ONE BLOCK'S ROWS: the block of the instants first..last, the first
  ! `depth` of them solved by the blocks before, with the stages
  ! between consecutive instants when the family marches by stages.
  ! Its rows are stated with at((m - 1) ncells + c) the offset of the
  ! tuple of moment m and cell c and primary(f) the row of equation f
  ! within a tuple, and its q is the estimate - each instant's tuple
  ! and, at a stage of the step into instant k + 1, the stage stored
  ! at stages(:, c, i, k + 1) when the block was solved before, the
  ! tuple of instant k otherwise.
  !===================================================================!

  !===================================================================!
  ! The position of every point of a block, moment outer and cell
  ! inner: the time of the moment - an instant's own, or for stage i
  ! of the step into instant k the instant behind plus the stage's
  ! abscissa times the step - then the cell's centre.
  !===================================================================!

  function block_positions(this, scheme, first, n, instant_moment, stage_moment) result(positions)

    class(discrete_residual), intent(in) :: this
    type(family)            , intent(in) :: scheme
    integer                 , intent(in) :: first, n, instant_moment(:), stage_moment(:,:)
    real(dp), allocatable :: positions(:,:)

    integer :: ncells, k, i, c, m

    ncells = this % points % num_cells()
    allocate(positions(this % points % num_coordinates(), maxval(instant_moment) * ncells))
    do k = 1, n
       do c = 1, ncells
          positions(:, (instant_moment(k) - 1) * ncells + c) = &
               & this % points % position(this % points % point_of(first + k - 1, c))
       end do
       if (k == 1) cycle
       do i = 1, size(stage_moment, 2)
          m = stage_moment(k - 1, i)
          if (m == 0) cycle
          do c = 1, ncells
             positions(:, (m - 1) * ncells + c) = this % points % position(this % points % point_of(first + k - 2, c))
             positions(1, (m - 1) * ncells + c) = this % points % instant(first + k - 2) &
                  & + scheme % stage_abscissa(i) * this % points % step(first + k - 1)
          end do
       end do
    end do

  end function block_positions

  subroutine block_rows(this, scheme, first, last, depth, rows_along, lb, tuple, stages, solved_stages, blk)

    class(discrete_residual), intent(in)  :: this
    type(family)            , intent(in)  :: scheme
    integer                 , intent(in)  :: first, last, depth
    type(stencil)           , intent(in)  :: rows_along(:,:)
    type(lower_bounds)      , intent(in)  :: lb
    real(dp)                , intent(in)  :: tuple(:,:,lb % instant:), stages(:,:,:,lb % stage:)
    logical                 , intent(in)  :: solved_stages
    type(chain_block)       , intent(out) :: blk

    integer , allocatable :: at(:), instant_moment(:), primary(:)
    type(stencil) :: relations, gauge
    type(connectivity_graph) :: connectivity
    integer , allocatable :: stage_moment(:,:), fixed_rows(:)
    integer , allocatable :: determined(:), source(:), gauge_rows(:), gauge_columns(:), face_rows(:), face_columns(:)
    real(dp), allocatable :: weight(:), fixed(:), w(:), steps(:), gauge_weights(:), zeros(:), face_weights(:)
    real(dp), allocatable :: face_constants(:), point_positions(:,:), space_weights(:)
    integer , allocatable :: space_rows(:), space_columns(:)
    logical , allocatable :: implicit_moment(:)
    logical :: open_line
    integer :: stride, nf, ncells, n, s, nm, m, k, i, f, c, e, count, nd, width, unknowns, npts
    integer :: comp, axis, o, j, ne, history_last, time_relations, a, chain, chains, cells
    integer , allocatable :: rule_primary(:)
    type(projection) :: projected
    type(expression) :: rule
    logical :: staged
    integer :: nodes, factor
    integer , allocatable :: related(:), added(:), single(:)
    logical , allocatable :: tying(:), entries(:), explicit_point(:)

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
    ! determined, the component it reads, the weight. Projected on the
    ! cells, every coefficient a node stores has a chain along time of
    ! its own, related by the family as the value's
    chains = 1
    if (allocated(this % elements)) chains = this % elements % num_coefficients()
    count = 0
    do f = 1, nf
       if (this % rule % degree_of_field(f) < 1) cycle
       if (staged) then
          connectivity = scheme % stage_connectivity(this % rule % degree_of_field(f) + 1)
          count = count + (n - 1) * connectivity % num_edges() * ncells * chains
       else
          connectivity = scheme % block_connectivity(this % rule % degree_of_field(f) + 1, n)
          count = count + connectivity % num_edges() * ncells * chains
       end if
    end do
    if (this % with_differences) then
       do f = 1, nf
          do c = 2, this % rule % num_coordinates()
             do o = 1, this % rule % degree_along(c)
                count = count + rows_along(c - 1, o) % pattern % num_edges() * nm &
                     & * (1 + this % rule % mixed_degree(f))
             end do
          end do
       end do
    end if
    allocate(determined(count), source(count), weight(count))
    e = 0

    do f = 1, nf
       if (this % rule % degree_of_field(f) < 1) cycle
       do chain = 0, chains - 1
          if (staged) then
             connectivity = scheme % stage_connectivity(this % rule % degree_of_field(f) + 1)
             do k = 1, n - 1
                steps = spread(this % points % step(first + k), 1, s + 2)
                call weights_of(scheme_weight(scheme), connectivity, steps, w)
                do j = 1, connectivity % num_edges()
                   do c = 1, ncells
                      e = e + 1
                      determined(e) = at_of(moment_of(connectivity % edge_head(j), k), c) &
                           & + this % coefficient_component(f, chain, connectivity % head_degree(j)) + 1
                      source(e)     = at_of(moment_of(connectivity % edge_tail(j), k), c) &
                           & + this % coefficient_component(f, chain, connectivity % tail_degree(j)) + 1
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
                        & + this % coefficient_component(f, chain, connectivity % head_degree(j)) + 1
                   source(e)     = at_of(instant_moment(connectivity % edge_tail(j)), c) &
                        & + this % coefficient_component(f, chain, connectivity % tail_degree(j)) + 1
                   weight(e)     = w(j)
                end do
             end do
          end if
       end do
    end do

    time_relations = e

    ! the implicit moments: those at which a relation of the family
    ! determines a component from a component of the same moment; at
    ! any other moment the family's quadrature of the moments before
    ! determines the values
    allocate(implicit_moment(nm), source=.false.)
    do j = 1, time_relations
       if ((determined(j) - 1) / width == (source(j) - 1) / width) then
          implicit_moment((determined(j) - 1) / width + 1) = .true.
       end if
    end do
    open_line = this % open_line()

    ! the mixed component of order a along time and o along space is
    ! the space derivative of order o of the time component of order
    ! a, the path along space chosen by the concretion; a = 0 is the
    ! pure space derivative of the values
    if (this % with_differences) then
       do f = 1, nf
          do a = 0, this % rule % mixed_degree(f)
             do c = 2, this % rule % num_coordinates()
                axis = c - 1
                do o = 1, this % rule % degree_along(c)
                   comp = this % rule % component_at(c, o, f, a)
                   associate (op => rows_along(axis, o))
                     call op % weights % real_vector(w)
                     ne = op % pattern % num_edges()
                     do m = 1, nm
                        ! on an open line the segment relations state the
                        ! derivatives at the implicit moments
                        if (open_line .and. implicit_moment(m)) cycle
                        do j = 1, ne
                           e = e + 1
                           determined(e) = at_of(m, op % pattern % edge_head(j)) + comp + 1
                           source(e)     = at_of(m, op % pattern % edge_tail(j)) + this % rule % offset_of_field(f) + a + 1
                           weight(e)     = w(j)
                        end do
                     end do
                   end associate
                end do
             end do
          end do
       end do
    end if

    relations = derived_constraints(determined(1:e), source(1:e), weight(1:e), unknowns, &
         & 'rows along time and space')
    related = determined(1:e)

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
    call this % fixed_rows_of(first, last, depth, instant_moment, at, lb, tuple, fixed_rows, fixed)

    ! THE PRESCRIBED ROWS: the gauge rows - at every moment after the
    ! history, the mean of the integrated unknown over the cells, in
    ! place of its own row at the first cell - and the face rows at
    ! the block's own instants, one stencil with the faces' constants
    history_last = 0
    if (depth > 0) history_last = instant_moment(depth)
    call this % gauge_rows_of(nm, history_last, at, primary, gauge_rows, gauge_columns, gauge_weights)

    ! THE EXPLICIT POINTS: those of the moments after the history at
    ! which the family's quadrature of the moments before determines
    ! every component below the top degree - the instant a step of a
    ! staged family reaches. As at the manifold's first instant, a law
    ! there that reads no component of top degree - a joint's condition
    ! - is stated by its derivative along time of the least order that
    ! reads one, which determines its multiplier there; the block of
    ! the manifold's first instant takes that rule at every point
    allocate(explicit_point(npts), source=.false.)
    if (.not. (first == 1 .and. depth == 0)) then
       do m = history_last + 1, nm
          if (.not. implicit_moment(m)) explicit_point((m - 1) * ncells + 1:m * ncells) = .true.
       end do
    end if
    point_positions = this % block_positions(scheme, first, n, instant_moment, stage_moment)
    call this % face_rows_of(first, last, depth, nm, instant_moment, at, primary, point_positions, rows_along, &
         & implicit_moment, face_rows, face_columns, face_weights, face_constants)
    if (open_line) then
       call this % segment_rows_of(nm, at, implicit_moment, space_rows, space_columns, space_weights)
       face_rows      = [face_rows, space_rows]
       face_columns   = [face_columns, space_columns]
       face_weights   = [face_weights, space_weights]
       face_constants = [face_constants, spread(0.0_dp, 1, size(space_rows))]
    end if
    factor = merge(1, 0, any(this % manifold % of_time_alone(1:nf)))
    if (factor > 0) then
       call this % time_alone_rows_of(depth, nm, instant_moment, at, primary, fixed_rows, face_rows, &
            & face_columns, face_weights, face_constants, added)
       ! THE TYING ROWS: the history, the family's relations, the rows of
       ! the unknowns of time alone, and the prescribed rows of one entry
       ! on their own diagonal; a prescribed row of several entries is
       ! not, though it take the row of a relation
       allocate(tying(unknowns), entries(unknowns), source=.false.)
       tying(fixed_rows) = .true.
       tying(related)    = .true.
       allocate(single(unknowns), source=0)
       do j = 1, size(face_rows)
          single(face_rows(j)) = single(face_rows(j)) + 1
          if (face_columns(j) == face_rows(j)) entries(face_rows(j)) = .true.
       end do
       do j = 1, size(face_rows)
          tying(face_rows(j)) = single(face_rows(j)) == 1 .and. entries(face_rows(j))
       end do
       tying(added) = .true.
    end if
    allocate(zeros(unknowns), source=0.0_dp)
    do e = 1, size(face_rows)
       zeros(face_rows(e)) = face_constants(e)
    end do
    ! PROJECTED ON THE CELLS: each element at each moment an instance,
    ! the law in flux form tested at its quadrature points against the
    ! basis, its rules each in the row of the value's coefficient of
    ! the equation it tests
    if (allocated(this % elements)) then
       if (size(gauge_rows) > 0) then
          error stop 'operation_residual: under a projection on the cells the integral of an unknown over the &
               &region is not stated'
       end if
       cells = this % elements % num_cells()
       nodes = this % elements % num_nodes()
       allocate(projected % element_point(nodes + factor, nm * cells), projected % element_cell(nm * cells))
       do m = 1, nm
          do e = 1, cells
             projected % element_point(1:nodes, (m - 1) * cells + e) = (m - 1) * ncells + this % elements % node(:, e)
             ! the first point of the moment, where the equations of the
             ! unknowns of time alone are stated
             if (factor > 0) projected % element_point(nodes + 1, (m - 1) * cells + e) = (m - 1) * ncells + 1
             projected % element_cell((m - 1) * cells + e) = e
          end do
       end do
       projected % elements = this % elements
       if (factor > 0) call with_first_point(projected % elements)
       projected % tested   = this % tested
       allocate(projected % coefficient_row(0:this % elements % num_coefficients() - 1, nf))
       do f = 1, nf
          do k = 0, this % elements % num_coefficients() - 1
             projected % coefficient_row(k, f) = this % coefficient_component(f, k, &
                  & scheme % primary_degree(this % rule % degree_of_field(f)))
          end do
       end do
       call this % formed_components(projected % formed)
       rule_primary = [(primary(this % tested(1, j)), j = 1, size(this % tested, 2))]
       ! the block of the manifold's first instant states the laws there
       rule = this % flux
       if (allocated(this % initial_flux) .and. first == 1 .and. depth == 0) then
          if (n > 1) then
             error stop 'operation_residual: a law the initial data satisfy is stated at the first instant by its &
                  &derivative along time, in a block of that instant alone'
          end if
          rule = this % initial_flux
       end if
       if (size(face_rows) > 0) then
          gauge = stencil(face_rows, face_columns, face_weights, zeros, 'prescribed')
          if (factor > 0) then
             blk % rows = residual_operator(relations, rule, at, unknowns, stride, rule_primary, fixed_rows, &
                  & fixed, prescribed=gauge, positions=point_positions, projected=projected, tying=tying, &
                  & explicit_rule=this % initial_flux, explicit=explicit_point)
          else
             blk % rows = residual_operator(relations, rule, at, unknowns, stride, rule_primary, fixed_rows, &
                  & fixed, prescribed=gauge, positions=point_positions, projected=projected, &
                  & explicit_rule=this % initial_flux, explicit=explicit_point)
          end if
       else
          blk % rows = residual_operator(relations, rule, at, unknowns, stride, rule_primary, fixed_rows, &
               & fixed, positions=point_positions, projected=projected, explicit_rule=this % initial_flux, &
               & explicit=explicit_point)
       end if
    else
       ! the block of the manifold's first instant states the laws there
       rule = this % rule
       if (allocated(this % initial_rule) .and. first == 1 .and. depth == 0) then
          if (n > 1) then
             error stop 'operation_residual: a law the initial data satisfy is stated at the first instant by its &
                  &derivative along time, in a block of that instant alone'
          end if
          rule = this % initial_rule
       end if
       if (size(gauge_rows) + size(face_rows) > 0) then
          gauge = stencil([gauge_rows, face_rows], [gauge_columns, face_columns], [gauge_weights, face_weights], &
               & zeros, 'prescribed')
          ! the gauge rows, one per moment and gauged field, read every
          ! cell of the moment
          blk % rows = residual_operator(relations, rule, at, unknowns, stride, primary, fixed_rows, fixed, &
               & prescribed=gauge, positions=point_positions, summed_rows=gauge_rows(1:size(gauge_rows):ncells), &
               & explicit_rule=this % initial_rule, explicit=explicit_point)
       else
          blk % rows = residual_operator(relations, rule, at, unknowns, stride, primary, fixed_rows, fixed, &
               & positions=point_positions, explicit_rule=this % initial_rule, explicit=explicit_point)
       end if
    end if
    blk % first      = first
    blk % last       = last
    blk % depth      = depth
    blk % num_stages = s
    blk % num_cells  = ncells
    call move_alloc(instant_moment, blk % instant_moment)

    ! the estimate: the tuple at each instant; at a stage of the step
    ! into instant k + 1, the stage solved before when solved_stages,
    ! the tuple of instant k otherwise
    allocate(blk % q(unknowns))
    if (solved_stages) then
       call blk % from_points(EVERY_MOMENT, lb, tuple, blk % q, stages)
    else
       call blk % from_points(EVERY_MOMENT, lb, tuple, blk % q)
       do k = 1, n - 1
          do i = 1, s
             do c = 1, ncells
                blk % q(at_of(stage_moment(k, i), c) + 1:at_of(stage_moment(k, i), c) + stride) = &
                     & tuple(:, c, first + k - 1)
             end do
          end do
       end do
    end if

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
  ! A BLOCK AND THE POINTS OF THE MANIFOLD. tuple_at is the offset,
  ! among the block's unknowns, of the tuple of the moment m at the
  ! cell c. from_points reads one column of the tuples over the
  ! instants, and of the stages of the step into each, into a vector
  ! over the block's unknowns at the moments named; to_points writes
  ! the column from the vector. The stages are left out when their
  ! array is.
  !===================================================================!

  pure integer function tuple_at(this, m, c)
    class(chain_block), intent(in) :: this
    integer           , intent(in) :: m, c
    tuple_at = this % rows % at((m - 1) * this % num_cells + c)
  end function tuple_at

  pure subroutine from_points(this, moments, lb, tuple, x, stages)

    class(chain_block), intent(in)    :: this
    integer           , intent(in)    :: moments
    type(lower_bounds), intent(in)    :: lb
    real(dp)          , intent(in)    :: tuple(:,:,lb % instant:)
    real(dp)          , intent(inout) :: x(:)
    real(dp)          , intent(in), optional :: stages(:,:,:,lb % stage:)

    integer :: k, i, c, base, stride

    stride = this % rows % degrees
    do k = merge(this % depth + 1, 1, moments == OWN_MOMENTS), &
         & merge(this % depth, this % last - this % first + 1, moments == HISTORY_MOMENTS)
       do c = 1, this % num_cells
          base = this % tuple_at(this % instant_moment(k), c)
          x(base + 1:base + stride) = tuple(:, c, this % first + k - 1)
       end do
    end do
    if (moments == HISTORY_MOMENTS .or. .not. present(stages)) return
    do k = 1, this % last - this % first
       do i = 1, this % num_stages
          do c = 1, this % num_cells
             base = this % tuple_at(this % instant_moment(k) + i, c)
             x(base + 1:base + stride) = stages(:, c, i, this % first + k)
          end do
       end do
    end do

  end subroutine from_points

  pure subroutine to_points(this, moments, lb, x, tuple, stages)

    class(chain_block), intent(in)    :: this
    integer           , intent(in)    :: moments
    type(lower_bounds), intent(in)    :: lb
    real(dp)          , intent(in)    :: x(:)
    real(dp)          , intent(inout) :: tuple(:,:,lb % instant:)
    real(dp)          , intent(inout), optional :: stages(:,:,:,lb % stage:)

    integer :: k, i, c, base, stride

    stride = this % rows % degrees
    do k = merge(this % depth + 1, 1, moments == OWN_MOMENTS), &
         & merge(this % depth, this % last - this % first + 1, moments == HISTORY_MOMENTS)
       do c = 1, this % num_cells
          base = this % tuple_at(this % instant_moment(k), c)
          tuple(:, c, this % first + k - 1) = x(base + 1:base + stride)
       end do
    end do
    if (moments == HISTORY_MOMENTS .or. .not. present(stages)) return
    do k = 1, this % last - this % first
       do i = 1, this % num_stages
          do c = 1, this % num_cells
             base = this % tuple_at(this % instant_moment(k) + i, c)
             stages(:, c, i, this % first + k) = x(base + 1:base + stride)
          end do
       end do
    end do

  end subroutine to_points

  !===================================================================!
  ! The rows' inputs at the block's q, by which their jacobian is
  ! formed there, and the image of every point of the block.
  !===================================================================!

  subroutine block_inputs(this, blk)

    class(discrete_residual), intent(in)    :: this
    type(chain_block)       , intent(inout) :: blk

    type(newton) :: solver
    real(dp), allocatable :: y(:)
    integer :: npts

    npts   = size(blk % rows % at)
    solver = this % solver_for(blk % rows, size(blk % q), blk % rows % degrees, npts)
    call solver % evaluate(blk % q, y, blk % inputs)
    blk % owner = this % point_owners(npts)

  end subroutine block_inputs

  !===================================================================!
  ! ONE BLOCK SOLVED: its rows formed at the tuples of column zero,
  ! and solved as one system of rows.
  !===================================================================!

  subroutine solve_block(this, scheme, first, last, depth, rows_along, lb, tuple, stages, expand)

    class(discrete_residual), intent(in)    :: this
    type(family)            , intent(in)    :: scheme
    integer                 , intent(in)    :: first, last, depth
    type(stencil)           , intent(in)    :: rows_along(:,:)
    type(lower_bounds)      , intent(in)    :: lb
    real(dp)                , intent(inout) :: tuple(:,:,lb % instant:,0:), stages(:,:,:,lb % stage:,0:)
    logical                 , intent(in), optional :: expand

    type(chain_block) :: blk
    real(dp), allocatable :: predicted(:,:)
    integer :: m, c

    if (verbosity >= 1) then
       print '(a,a,a,i0,a,i0,a,es12.4,a,es12.4,a,i0,a)', 'block  ', trim(scheme % name()), '  instants ', first, &
            & '..', last, '  t = ', this % points % instant(first), ' .. ', this % points % instant(last), &
            & '  history ', depth, ' instants'
    end if
    call this % block_rows(scheme, first, last, depth, rows_along, lb, tuple(:, :, :, 0), stages(:, :, :, :, 0), &
         & .false., blk)
    ! continued, every moment after the history starts from the tuple
    ! of the last instant of the history, its unknowns of time alone
    ! extrapolated along time to the block's last instant
    if (this % continued .and. depth > 0) then
       associate (stride => this % rule % num_components())
         predicted = tuple(:, :, first + depth - 1, 0)
         call this % extrapolated(predicted, this % points % instant(last) - this % points % instant(first + depth - 1))
         do m = blk % instant_moment(depth) + 1, size(blk % rows % at) / blk % num_cells
            do c = 1, blk % num_cells
               blk % q(blk % tuple_at(m, c) + 1:blk % tuple_at(m, c) + stride) = predicted(:, c)
            end do
         end do
       end associate
    end if
    call this % solve_rows(blk, lb, tuple, stages, expand)

  end subroutine solve_block

  !===================================================================!
  ! ONE SYSTEM OF ROWS SOLVED - a block's, or the point manifold's:
  ! Newton from the estimate, the block's own instants and the stages
  ! of its steps written to column zero, the stages retained for the
  ! adjoint, which linearizes the block at its solution; then every
  ! column above, expansion by expansion and order by order.
  !===================================================================!

  subroutine solve_rows(this, blk, lb, tuple, stages, expand)

    class(discrete_residual), intent(in)    :: this
    type(chain_block)       , intent(inout) :: blk
    type(lower_bounds)      , intent(in)    :: lb
    real(dp)                , intent(inout) :: tuple(:,:,lb % instant:,0:), stages(:,:,:,lb % stage:,0:)
    logical                 , intent(in), optional :: expand

    integer :: d, m

    call this % solved(blk % rows, size(blk % q), blk % rows % degrees, size(blk % rows % at), blk % q)
    call blk % to_points(OWN_MOMENTS, lb, blk % q, tuple(:, :, :, 0), stages(:, :, :, :, 0))

    if (this % points % expansion % order == 0) return
    if (present(expand)) then
       if (.not. expand) return
    end if
    call this % block_inputs(blk)
    do d = 1, this % points % expansion % num_expansions()
       if (.not. this % expansion_owned(d)) cycle
       do m = 1, this % points % expansion % order
          if (verbosity >= 1) print '(a,i0,a,i0)', 'expansion ', d, ' along the design, order ', m
          call this % expansion_order(blk, lb, d, m, tuple, stages)
       end do
    end do

  end subroutine solve_rows

  !===================================================================!
  ! ONE ORDER OF THE EXPANSION of a block: the m-th derivative of the
  ! block's unknowns along the expansion d is the solution of A x = b,
  ! A the block's jacobian at its solution and b the m-th derivative
  ! along d of the rows with the m-th derivative of the state left out
  ! - at the row of an equation, by the jet arithmetic on the
  ! derivatives of lower order and the design; at the fixed row of a
  ! history instant, the derivative solved before; at the fixed row of
  ! a face condition, the m-th derivative of the value it states; zero
  ! at the rows along time and space, which are linear in the state.
  ! The derivatives below m are read from their columns of the tuples
  ! and of the stages, the m-th at the history instants from its own,
  ! and the solution is written to that column at the block's own
  ! instants and at the stages of its steps.
  !===================================================================!

  subroutine expansion_order(this, blk, lb, d, m, tuple, stages)

    class(discrete_residual), intent(in)    :: this
    type(chain_block)       , intent(in)    :: blk
    type(lower_bounds)      , intent(in)    :: lb
    integer                 , intent(in)    :: d, m
    real(dp)                , intent(inout) :: tuple(:,:,lb % instant:,0:), stages(:,:,:,lb % stage:,0:)

    type(residual_operator) :: lin
    real(dp), allocatable :: jets(:,:), rhs(:), x(:), slopes(:), values(:), seed(:)
    integer , allocatable :: columns(:)
    real(dp) :: given
    integer :: unknowns, column, o, c, base, j, i, instant, row

    unknowns = size(blk % q)
    values   = this % points % expansion % values()
    seed     = this % points % expansion % seed(d)

    allocate(jets(unknowns, 0:m - 1), rhs(unknowns), source=0.0_dp)
    jets(:, 0) = blk % q
    do o = 1, m - 1
       column = this % points % expansion % column(d, o)
       call blk % from_points(EVERY_MOMENT, lb, tuple(:, :, :, column), jets(:, o), stages(:, :, :, :, column))
    end do
    column = this % points % expansion % column(d, m)
    call design_source(blk % rows, jets, m, values, seed, blk % owner, rhs)
    call blk % from_points(HISTORY_MOMENTS, lb, tuple(:, :, :, column), rhs)
    do j = 1, size(this % condition)
       if (.not. this % condition(j) % manifold % is_time_face()) cycle
       instant = this % face_instant(this % condition(j) % manifold % face_value)
       if (instant < blk % first + blk % depth .or. instant > blk % last) cycle
       if (allocated(this % elements)) then
          ! projected: the derivative of each row's datum is the
          ! integral of the derivative of the condition's against the
          ! row's basis function
          do i = 1, size(this % condition(j) % equation)
             call this % projected_datum(this % condition(j) % equation(i), instant, m, values, seed, blk, rhs)
          end do
          cycle
       end if
       do i = 1, size(this % condition(j) % equation)
          do c = 1, blk % num_cells
             call this % face_relation(this % condition(j) % equation(i), values, &
                  & this % points % position(this % points % point_of(instant, c)), columns, slopes, given, row)
             base = blk % tuple_at(this % face_row_moment(instant, blk % first, blk % last, blk % depth, &
                  & blk % instant_moment), c)
             rhs(base + row + 1) = stated_derivative( &
                  & this % condition(j) % equation(i), m, values, seed, &
                  & this % points % position(this % points % point_of(instant, c)))
          end do
       end do
    end do
    lin = blk % rows % linearize(blk % rows % unknown_graph(), blk % rows % bind(blk % inputs), rhs, &
         & transposed=.false., version_number=blk % rows % version())
    allocate(x(unknowns), source=0.0_dp)
    call this % solved(lin, unknowns, blk % rows % degrees, size(blk % rows % at), x)
    call blk % to_points(OWN_MOMENTS, lb, x, tuple(:, :, :, column), stages(:, :, :, :, column))

  end subroutine expansion_order

  !===================================================================!
  ! THE WAVEFRONT over the blocks and the orders: on a manifold of
  ! time alone with one expansion, when the images number the orders
  ! plus one, image m + 1 forms order m of every block, one block
  ! behind image m. Order m of block b reads the state of block b and
  ! the orders below m, into their columns, from the images that
  ! formed them, and its own order m at the block's history instants,
  ! formed at the block before; each block ends with a rendezvous with
  ! the next image, so that every image works while the others do. The
  ! one expansion's column m is its order m. At the end every image
  ! has column zero and its own order's column, the others above zero
  ! set to zero, for the sum over the images that minimize takes.
  !===================================================================!

  subroutine wavefront_expansion(this, rows_along, tuple, stages)

    class(discrete_residual), intent(in)    :: this
    type(stencil)           , intent(in)    :: rows_along(:,:)
    real(dp)                , intent(inout) :: tuple(:,:,:,0:), stages(:,:,:,:,0:)

    real(dp), allocatable :: tuple_co(:,:,:,:)[:], stages_co(:,:,:,:,:)[:]
    type(chain_block) :: blk
    integer :: m, me, b, first, last, depth, o, s, order

    order = this % points % expansion % order
    me    = this_image()
    m     = me - 1

    allocate(tuple_co(size(tuple, 1), size(tuple, 2), size(tuple, 3), 0:order)[*])
    allocate(stages_co(size(stages, 1), size(stages, 2), size(stages, 3), size(stages, 4), 0:order)[*])
    tuple_co  = 0.0_dp
    stages_co = 0.0_dp

    do b = 1, this % schemes % num_blocks()
       call this % block_extent(b, first, last, depth)
       associate (scheme => this % schemes % scheme(b))
         s = 0
         if (marches_by_stages(scheme, top_degree(this % rule, this % manifold % num_unknowns) + 1)) then
            s = scheme % num_stages()
         end if
         if (m == 0) then
            ! the state, published to the image of order one
            call this % solve_block(scheme, first, last, depth, rows_along, origin, tuple, stages, expand=.false.)
         else
            sync images(me - 1)
            do o = 0, m - 1
               tuple(:, :, first:last, o) = tuple_co(:, :, first:last, o)[o + 1]
               if (s > 0) stages(:, :, 1:s, first + 1:last, o) = stages_co(:, :, 1:s, first + 1:last, o)[o + 1]
            end do
            if (verbosity >= 1) print '(a,i0,a,i0,a,i0)', 'wavefront: image ', me, ' forms order ', m, ' of block ', b
            call this % block_rows(scheme, first, last, depth, rows_along, origin, tuple(:, :, :, 0), &
                 & stages(:, :, :, :, 0), .true., blk)
            call this % block_inputs(blk)
            call this % expansion_order(blk, origin, 1, m, tuple, stages)
         end if
         tuple_co(:, :, first + depth:last, m) = tuple(:, :, first + depth:last, m)
         if (s > 0) stages_co(:, :, 1:s, first + 1:last, m) = stages(:, :, 1:s, first + 1:last, m)
         if (me < num_images()) sync images(me + 1)
       end associate
    end do
    ! the columns of the other orders, read from the other images, left
    ! to the sum
    do o = 1, order
       if (o == m) cycle
       tuple(:, :, :, o)     = 0.0_dp
       stages(:, :, :, :, o) = 0.0_dp
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

    order = this % points % expansion % order
    pipelined = this % with_chain .and. .not. this % manifold % with_space .and. order >= 1 &
         & .and. this % points % expansion % num_expansions() == 1 .and. num_images() == order + 1

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

  subroutine wavefront_adjoint(this, rows_along, tuple, stages, gradient, stage_node, &
       & adjoint, incoming, reaction, gauge, in_design)

    class(discrete_residual), intent(in)    :: this
    type(stencil)           , intent(in)    :: rows_along(:,:)
    real(dp)                , intent(in)    :: tuple(:,:,:,0:), stages(:,:,:,:,0:)
    real(dp)                , intent(in)    :: gradient(:,:,0:)
    integer                 , intent(in)    :: stage_node(:,:,:)
    real(dp)                , intent(inout) :: adjoint(:,:,:,0:), incoming(:,:,:,0:), reaction(:,:,0:)
    real(dp)                , intent(inout) :: gauge(:,:,0:), in_design(:,0:)

    real(dp), allocatable :: mu_co(:,:,:)[:], smu_co(:,:,:,:)[:], reaction_co(:,:,:)[:]
    real(dp), allocatable :: below(:,:,:), stages_below(:,:,:,:)
    type(chain_block) :: blk
    real(dp), allocatable :: jets(:,:), mu(:,:), x(:)
    logical , allocatable :: gauged(:)
    integer :: o, me, b, first, last, depth, j, s, stride, ncells, ninst, order, smax

    stride = this % rule % num_components()
    ncells = this % points % num_cells()
    ninst  = this % points % num_instants()
    order  = this % points % expansion % order
    smax   = size(stages, 3)
    me     = this_image()
    o      = me - 1
    gauged = this % gauged_unknowns()

    allocate(mu_co(stride, ncells, ninst)[*], smu_co(stride, ncells, smax, ninst)[*])
    allocate(reaction_co(ncells, size(reaction, 2), 0:order)[*])
    mu_co = 0.0_dp
    smu_co = 0.0_dp
    reaction_co = 0.0_dp
    allocate(below(stride, ncells, ninst), stages_below(stride, ncells, smax, ninst), source=0.0_dp)

    do b = this % schemes % num_blocks(), 1, -1
       call this % block_extent(b, first, last, depth)
       if (o > 0) sync images(me - 1)
       if (verbosity >= 1) print '(a,i0,a,i0,a,i0)', 'wavefront: image ', me, ' forms the adjoint of order ', o, &
            & ' of block ', b
       call this % block_rows(this % schemes % scheme(b), first, last, depth, rows_along, origin, tuple(:, :, :, 0), &
            & stages(:, :, :, :, 0), .true., blk)
       call this % block_inputs(blk)
       s = blk % num_stages
       ! the jets of the state to order o, and the multipliers of the
       ! orders below at the block's own moments, from the images that
       ! formed them
       allocate(jets(size(blk % q), 0:o), mu(size(blk % q), 0:o), source=0.0_dp)
       jets(:, 0) = blk % q
       do j = 1, o
          call blk % from_points(EVERY_MOMENT, origin, tuple(:, :, :, j), jets(:, j), stages(:, :, :, :, j))
          below(:, :, first + depth:last) = mu_co(:, :, first + depth:last)[j]
          if (s > 0) stages_below(:, :, 1:s, first + 1:last) = smu_co(:, :, 1:s, first + 1:last)[j]
          call blk % from_points(OWN_MOMENTS, origin, below, mu(:, j - 1), stages_below)
          reaction(:, :, j - 1) = reaction_co(:, :, j - 1)[j]
       end do
       call this % adjoint_order(origin, blk, o, 1, jets, mu(:, 0:max(o - 1, 0)), gradient, stage_node, incoming, x)
       mu(:, o) = x
       call blk % to_points(OWN_MOMENTS, origin, x, below, stages_below)
       mu_co(:, :, first + depth:last) = below(:, :, first + depth:last)
       if (s > 0) smu_co(:, :, 1:s, first + 1:last) = stages_below(:, :, 1:s, first + 1:last)
       call this % multipliers_of_order(blk, o, 1, gauged, jets, mu, adjoint, gauge, reaction, in_design)
       reaction_co(:, :, o) = reaction(:, :, o)
       deallocate(jets, mu, x)
       if (me < num_images()) sync images(me + 1)
    end do
    ! the reactions of the other orders, read from the other images,
    ! left to the sum
    do j = 0, order
       if (j /= o) reaction(:, :, j) = 0.0_dp
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

    if (allocated(rows % projected)) then
       call projected_design_source(rows, jets, m, values, seed, owner, b)
       return
    end if
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
          r = rows % rules(j) % at_instant(q, design, rows % point_positions(:, p))
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
  ! each expansion: the transposed linear system A^T mu = rhs
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

  subroutine adjoint_block(this, scheme, first, last, depth, rows_along, lb, tuple, stages, gradient, &
       & stage_node, adjoint, incoming, reaction, gauge, in_design)

    class(discrete_residual), intent(in)    :: this
    type(family)            , intent(in)    :: scheme
    integer                 , intent(in)    :: first, last, depth
    type(stencil)           , intent(in)    :: rows_along(:,:)
    type(lower_bounds)      , intent(in)    :: lb
    real(dp)                , intent(in)    :: tuple(:,:,lb % instant:,0:), stages(:,:,:,lb % stage:,0:)
    real(dp)                , intent(in)    :: gradient(:,lb % node:,0:)
    integer                 , intent(in)    :: stage_node(lb % stage:,:,:)
    real(dp)                , intent(inout) :: adjoint(:,:,:,0:), incoming(:,:,lb % instant:,0:), reaction(:,:,0:)
    real(dp)                , intent(inout) :: gauge(:,:,0:), in_design(:,0:)

    type(chain_block) :: blk

    if (verbosity >= 1) then
       print '(a,a,a,i0,a,i0,a,es12.4,a,es12.4)', 'adjoint  ', trim(scheme % name()), '  instants ', &
            & first + depth, '..', last, '  t = ', this % points % instant(first + depth), ' .. ', &
            & this % points % instant(last)
    end if
    call this % block_rows(scheme, first, last, depth, rows_along, lb, tuple(:, :, :, 0), stages(:, :, :, :, 0), &
         & .true., blk)
    call this % adjoint_rows(blk, lb, tuple, stages, gradient, stage_node, adjoint, incoming, reaction, gauge, &
         & in_design)

  end subroutine adjoint_block

  !===================================================================!
  ! THE ADJOINT OF ONE SYSTEM OF ROWS - a block's, or the point
  ! manifold's, which has one moment, no history and no stages - as
  ! described above.
  !===================================================================!

  subroutine adjoint_rows(this, blk, lb, tuple, stages, gradient, stage_node, adjoint, incoming, reaction, gauge, &
       & in_design)

    class(discrete_residual), intent(in)    :: this
    type(chain_block)       , intent(inout) :: blk
    type(lower_bounds)      , intent(in)    :: lb
    real(dp)                , intent(in)    :: tuple(:,:,lb % instant:,0:), stages(:,:,:,lb % stage:,0:)
    real(dp)                , intent(in)    :: gradient(:,lb % node:,0:)
    integer                 , intent(in)    :: stage_node(lb % stage:,:,:)
    real(dp)                , intent(inout) :: adjoint(:,:,:,0:), incoming(:,:,lb % instant:,0:), reaction(:,:,0:)
    real(dp)                , intent(inout) :: gauge(:,:,0:), in_design(:,0:)

    real(dp), allocatable :: mu(:,:,:), x(:), jets(:,:,:)
    logical , allocatable :: gauged(:)
    integer :: o, order, d, nd, column

    order = this % points % expansion % order
    nd    = max(1, this % points % expansion % num_expansions())

    ! the jets of the state along each expansion at every unknown of
    ! the block: at the instants and at the stages, from the columns
    ! of the expansion
    allocate(jets(size(blk % q), 0:order, nd), mu(size(blk % q), 0:order, nd), source=0.0_dp)
    do d = 1, nd
       jets(:, 0, d) = blk % q
       do o = 1, order
          column = this % points % expansion % column(d, o)
          call blk % from_points(EVERY_MOMENT, lb, tuple(:, :, :, column), jets(:, o, d), stages(:, :, :, :, column))
       end do
    end do
    gauged = this % gauged_unknowns()
    call this % block_inputs(blk)

    ! column zero once, with the first expansion; the columns above it
    ! expansion by expansion
    do d = 1, nd
       do o = merge(1, 0, d > 1), order
          if (o > 0 .and. .not. this % expansion_owned(d)) cycle
          call this % adjoint_order(lb, blk, o, d, jets(:, 0:o, d), mu(:, 0:max(o - 1, 0), d), gradient, stage_node, &
               & incoming, x)
          if (o == 0) then
             mu(:, 0, :) = spread(x, 2, nd)
          else
             mu(:, o, d) = x
          end if
          deallocate(x)
       end do
    end do

    do d = 1, nd
       do o = merge(1, 0, d > 1), order
          if (o > 0 .and. .not. this % expansion_owned(d)) cycle
          call this % multipliers_of_order(blk, o, d, gauged, jets(:, 0:o, d), mu(:, 0:o, d), adjoint, gauge, &
               & reaction, in_design)
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
  ! the multipliers of the orders below. The multipliers at the
  ! block's history rows are passed to the blocks that determine those
  ! instants.
  !===================================================================!

  subroutine adjoint_order(this, lb, blk, o, d, jets, mu, gradient, stage_node, incoming, x)

    class(discrete_residual), intent(in)    :: this
    type(lower_bounds)      , intent(in)    :: lb
    type(chain_block)       , intent(in)    :: blk
    integer                 , intent(in)    :: o, d
    real(dp)                , intent(in)    :: jets(:,0:), mu(:,0:), gradient(:,lb % node:,0:)
    integer                 , intent(in)    :: stage_node(lb % stage:,:,:)
    real(dp)                , intent(inout) :: incoming(:,:,lb % instant:,0:)
    real(dp), allocatable   , intent(out)   :: x(:)

    type(residual_operator) :: lin
    real(dp), allocatable :: rhs(:)
    integer :: stride, n, k, c, i, base, unknowns, column, nsites, site
    integer, allocatable :: formed(:,:)

    stride   = blk % rows % degrees
    n        = blk % last - blk % first + 1
    unknowns = size(blk % q)
    column   = this % points % expansion % column(d, o)
    if (verbosity >= 1 .and. o > 0) print '(a,i0,a,i0)', 'adjoint expansion ', d, ' along the design, order ', o
    allocate(rhs(unknowns), source=0.0_dp)
    if (allocated(this % elements)) then
       ! projected: the differential at the sites, returned to the
       ! nodes by the transpose of the basis that formed the sites' jets
       nsites = this % elements % num_cells() * this % elements % num_quadrature_points()
       call this % formed_components(formed)
       do k = blk % depth + 1, n
          do c = 1, blk % num_cells
             base = blk % tuple_at(blk % instant_moment(k), c)
             rhs(base + 1:base + stride) = incoming(:, c, blk % first + k - 1, column)
          end do
          do site = 1, nsites
             call this % site_gradient(formed, site, gradient(:, (blk % first + k - 2) * nsites + site, column), &
                  & blk, blk % instant_moment(k), rhs)
          end do
       end do
       do k = 1, n - 1
          do i = 1, blk % num_stages
             do site = 1, nsites
                call this % site_gradient(formed, site, gradient(:, stage_node(blk % first + k, i, site), column), &
                     & blk, blk % instant_moment(k) + i, rhs)
             end do
          end do
       end do
    else
       do k = blk % depth + 1, n
          do c = 1, blk % num_cells
             base = blk % tuple_at(blk % instant_moment(k), c)
             rhs(base + 1:base + stride) = -gradient(1:stride, this % points % point_of(blk % first + k - 1, c), &
                  & column) + incoming(:, c, blk % first + k - 1, column)
          end do
       end do
       ! the objective's differential at the stages, the seeds of the
       ! stage quadrature
       do k = 1, n - 1
          do i = 1, blk % num_stages
             do c = 1, blk % num_cells
                base = blk % tuple_at(blk % instant_moment(k) + i, c)
                rhs(base + 1:base + stride) = -gradient(1:stride, stage_node(blk % first + k, i, c), column)
             end do
          end do
       end do
    end if
    if (o > 0) call adjoint_source(blk % rows, jets, mu(:, 0:o - 1), o, this % points % expansion % values(), &
         & this % points % expansion % seed(d), blk % owner, rhs)
    lin = blk % rows % linearize(blk % rows % unknown_graph(), blk % rows % bind(blk % inputs), rhs, &
         & transposed=.true., version_number=blk % rows % version())
    allocate(x(unknowns), source=0.0_dp)
    call this % solved(lin, unknowns, stride, size(blk % rows % at), x)
    do k = 1, blk % depth
       do c = 1, blk % num_cells
          base = blk % tuple_at(blk % instant_moment(k), c)
          incoming(:, c, blk % first + k - 1, column) = incoming(:, c, blk % first + k - 1, column) &
               & + x(base + 1:base + stride)
       end do
    end do

  end subroutine adjoint_order

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

  subroutine multipliers_of_order(this, blk, o, d, gauged, jets, mu, adjoint, gauge, reaction, in_design)

    class(discrete_residual), intent(in)    :: this
    type(chain_block)       , intent(in)    :: blk
    integer                 , intent(in)    :: o, d
    logical                 , intent(in)    :: gauged(:)
    real(dp)                , intent(in)    :: jets(:,0:), mu(:,0:)
    real(dp)                , intent(inout) :: adjoint(:,:,:,0:), gauge(:,:,0:), reaction(:,:,0:)
    real(dp)                , intent(inout) :: in_design(:,0:)

    ! the orders of a slope from zero: an allocatable assigned the
    ! function's result would have them from one
    real(dp) :: slope(0:o)
    real(dp), allocatable :: products(:), slopes(:), values(:), seed(:)
    integer , allocatable :: columns(:), target(:), row_point(:), row_offset(:), field(:), determined(:,:)
    integer , allocatable :: moved(:), determined_degrees(:)
    logical , allocatable :: occupied(:)
    real(dp) :: given
    integer :: stride, nf, ncells, n, k, c, f, m, j, i, instant, degree, e, g, ni, face_point, row
    integer :: law(3), law_row, q, column, lo, first, last, depth, s
    logical :: along_time
    integer , allocatable :: equation_row(:), on_face(:)
    real(dp), allocatable :: position(:)
    logical , allocatable :: attributed(:)
    real(dp) :: tau
    integer :: node, el_cell, qp, i0, k0, row_field, row_order, kd

    stride = this % rule % num_components()
    nf     = this % manifold % num_unknowns
    ncells = this % points % num_cells()
    first  = blk % first
    last   = blk % last
    depth  = blk % depth
    s      = blk % num_stages
    n      = last - first + 1
    values = this % points % expansion % values()
    ! the row of each equation within a point's tuple: projected, the
    ! row of its value's coefficient
    allocate(equation_row(nf))
    do f = 1, nf
       if (allocated(this % elements)) then
          equation_row(f) = blk % rows % projected % coefficient_row(0, f)
       else
          equation_row(f) = blk % rows % primary(f)
       end if
    end do
    seed   = this % points % expansion % seed(d)
    ni     = size(values)
    ! the column of the order o along the expansion d, and the first
    ! of the expansion's own
    column = this % points % expansion % column(d, o)
    lo     = this % points % expansion % column(d, 1)

    ! the rows the conditions on the faces of the region occupy, and the
    ! physics rows their difference relations moved to: none of their
    ! multipliers is the adjoint of the physics
    allocate(occupied(size(mu, 1)), source=.false.)
    allocate(moved(size(mu, 1)), source=0)
    do j = 1, size(this % condition)
       if (.not. this % condition(j) % manifold % is_space_face()) cycle
       if (allocated(this % elements)) then
          ! projected: the rows of the coefficients the condition fixes at
          ! the nodes on the face
          on_face = this % points_on_face(j)
          do i = 1, size(this % condition(j) % equation)
             do m = 1, size(blk % rows % at) / ncells
                do node = 1, size(on_face)
                   row = blk % tuple_at(m, on_face(node)) + this % fixed_coefficient(j, i, f, &
                        & equation_row) + 1
                   if (blk % rows % prescribed_row(row)) occupied(row) = .true.
                end do
             end do
          end do
          cycle
       end if
       call this % space_face_placement(j, blk % rows % primary, face_point, target, row_point, row_offset, field, determined)
       do i = 1, size(this % condition(j) % equation)
          do m = 1, size(blk % rows % at) / ncells
             row = blk % tuple_at(m, row_point(i)) + row_offset(i) + 1
             if (.not. blk % rows % prescribed_row(row)) cycle
             occupied(row) = .true.
             if (target(i) > 0) occupied(blk % tuple_at(m, target(i)) + blk % rows % primary(field(i)) + 1) = .true.
             if (determined(2, i) > 0) moved(row) = blk % tuple_at(m, target(i)) + blk % rows % primary(field(i)) + 1
          end do
       end do
    end do

    do k = depth + 1, n
       do c = 1, ncells
          do f = 1, nf
             adjoint(f, c, first + k - 1, column) = 0.0_dp
             do m = blk % instant_moment(k) - merge(s, 0, k > 1), blk % instant_moment(k)
                if (c == 1 .and. gauged(f)) cycle
                row = blk % tuple_at(m, c) + equation_row(f) + 1
                if (occupied(row)) cycle
                adjoint(f, c, first + k - 1, column) = adjoint(f, c, first + k - 1, column) + mu(row, o)
             end do
          end do
       end do
    end do
    ! the multipliers of the representation laws, summed over the
    ! moments of the step into each instant: law k along time of a
    ! field, the family's row of the k-th component it determines, the
    ! initial instant aside, whose rows are the initial data's; law k
    ! along a space coordinate, the difference row of that component,
    ! or where a face condition took it, the row the difference moved to
    e = 0
    do j = 1, size(this % representation)
       do i = 1, size(this % representation(j) % equation)
          e   = e + 1
          law = represented_component(this % representation(j) % equation(i) % graph(1))
          f   = law(1)
          along_time = this % manifold % with_time .and. law(3) == FIRST_COORDINATE
          if (along_time) then
             determined_degrees = pack([(q, q = 0, this % rule % degree_of_field(f))], &
                  & [(q /= blk % rows % primary(f) - this % rule % offset_of_field(f), &
                  &   q = 0, this % rule % degree_of_field(f))])
             law_row = this % rule % offset_of_field(f) + determined_degrees(law(2))
          else
             law_row = this % rule % component_at(law(3), law(2), f)
          end if
          do k = depth + 1, n
             do c = 1, ncells
                adjoint(nf + e, c, first + k - 1, column) = 0.0_dp
                if (along_time .and. first + k - 1 == 1) cycle
                do m = blk % instant_moment(k) - merge(s, 0, k > 1), blk % instant_moment(k)
                   row = blk % tuple_at(m, c) + law_row + 1
                   if (moved(row) > 0) row = moved(row)
                   if (occupied(row)) cycle
                   adjoint(nf + e, c, first + k - 1, column) = adjoint(nf + e, c, first + k - 1, column) + mu(row, o)
                end do
             end do
          end do
       end do
    end do
    ! the multipliers that are functions of time alone, summed over the
    ! moments of the step into each instant: of a condition on the time
    ! factor, at its row at the first cell; of a condition on a face of
    ! the region, at the row it occupies where it acts
    g = 0
    do j = 1, size(this % condition)
       if (this % condition(j) % manifold % is_time_factor()) then
          do i = 1, size(this % condition(j) % equation)
             g = g + 1
             f = this % condition(j) % equation(i) % index(1)
             do k = depth + 1, n
                do m = blk % instant_moment(k) - merge(s, 0, k > 1), blk % instant_moment(k)
                   row = blk % tuple_at(m, 1) + blk % rows % primary(f) + 1
                   gauge(first + k - 1, g, column) = gauge(first + k - 1, g, column) + mu(row, o)
                end do
             end do
          end do
       else if (this % condition(j) % manifold % is_space_face() .and. allocated(this % elements)) then
          ! projected: the multiplier of the row the condition takes at
          ! the nodes on the face, summed over the nodes; a node on
          ! several faces belongs to the condition stated first
          on_face = this % points_on_face(j)
          if (.not. allocated(attributed)) allocate(attributed(size(mu, 1)), source=.false.)
          do i = 1, size(this % condition(j) % equation)
             g = g + 1
             do k = depth + 1, n
                do m = blk % instant_moment(k) - merge(s, 0, k > 1), blk % instant_moment(k)
                   do node = 1, size(on_face)
                      row = blk % tuple_at(m, on_face(node)) + this % fixed_coefficient(j, i, f, equation_row) + 1
                      if (.not. blk % rows % prescribed_row(row) .or. attributed(row)) cycle
                      attributed(row) = .true.
                      gauge(first + k - 1, g, column) = gauge(first + k - 1, g, column) + mu(row, o)
                   end do
                end do
             end do
          end do
       else if (this % condition(j) % manifold % is_space_face()) then
          call this % space_face_placement(j, blk % rows % primary, face_point, target, row_point, row_offset, field, determined)
          do i = 1, size(this % condition(j) % equation)
             g = g + 1
             do k = depth + 1, n
                do m = blk % instant_moment(k) - merge(s, 0, k > 1), blk % instant_moment(k)
                   row = blk % tuple_at(m, row_point(i)) + row_offset(i) + 1
                   if (.not. blk % rows % prescribed_row(row)) cycle
                   gauge(first + k - 1, g, column) = gauge(first + k - 1, g, column) + mu(row, o)
                end do
             end do
          end do
       end if
    end do
    e = 0
    do j = 1, size(this % condition)
       if (.not. this % condition(j) % manifold % is_time_face()) cycle
       instant = this % face_instant(this % condition(j) % manifold % face_value)
       do i = 1, size(this % condition(j) % equation)
          e = e + 1
          if (instant < first + depth .or. instant > last) cycle
          m = this % face_row_moment(instant, first, last, depth, blk % instant_moment)
          do c = 1, ncells
             call this % face_relation(this % condition(j) % equation(i), values, &
                  & this % points % position(this % points % point_of(instant, c)), columns, slopes, given, row)
             reaction(c, e, column) = mu(blk % tuple_at(m, c) + row + 1, o)
          end do
       end do
    end do

    if (.not. this % manifold % with_design) return
    products = design_products(blk % rows, jets, mu, o, values, seed, blk % owner)
    in_design(:, column) = in_design(:, column) - products
    e = 0
    do j = 1, size(this % condition)
       if (.not. this % condition(j) % manifold % is_time_face()) cycle
       instant = this % face_instant(this % condition(j) % manifold % face_value)
       do i = 1, size(this % condition(j) % equation)
          e = e + 1
          if (instant < first + depth .or. instant > last) cycle
          if (allocated(this % elements)) then
             ! projected: each row of the condition is its integral
             ! against a basis function, and so is the partial of its
             ! datum, weighted by the row's multiplier
             m = this % face_row_moment(instant, first, last, depth, blk % instant_moment)
             associate (el => this % elements)
               do el_cell = 1, el % num_cells()
                  do qp = 1, el % num_quadrature_points()
                     position = this % points % position(this % points % point_of(instant, el % node(1, el_cell)))
                     position(size(position) - el % dimension + 1:) = el % position(:, qp, el_cell)
                     call this % face_relation(this % condition(j) % equation(i), values, position, columns, slopes, &
                          & given, row)
                     call this % stored_order(row, row_field, row_order)
                     do kd = 1, ni
                        slope = stated_partial_jet(this % condition(j) % equation(i), kd, seed, o, values, position)
                        do i0 = 1, el % num_nodes()
                           do k0 = 0, el % num_coefficients() - 1
                              tau = el % test(0, k0, i0, qp, el_cell)
                              q   = blk % tuple_at(m, el % node(i0, el_cell)) &
                                   & + this % coefficient_component(row_field, k0, row_order) + 1
                              in_design(kd, column) = in_design(kd, column) + mixed_partial( &
                                   & symmetric_terms(mu(q, 0), mu(q, 1:o), o) &
                                   & * symmetric_terms(tau * slope(0), tau * slope(1:o), o))
                           end do
                        end do
                     end do
                  end do
               end do
             end associate
             cycle
          end if
          do c = 1, ncells
             do k = 1, ni
                slope = stated_partial_jet(this % condition(j) % equation(i), k, seed, o, values, &
                     & this % points % position(this % points % point_of(instant, c)))
                ! the row x - h(nu): its partial in nu_k is -dh/dnu_k, and
                ! the product of the row's multiplier with it is
                ! subtracted, as the equation rows' products are
                in_design(k, column) = in_design(k, column) + mixed_partial( &
                     & symmetric_terms(reaction(c, e, 0), reaction(c, e, lo:lo + o - 1), o) &
                     & * symmetric_terms(slope(0), slope(1:o), o))
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
    integer, allocatable :: reads(:)
    integer :: deg, npts, j, p, c, k, row, m

    if (allocated(rows % projected)) then
       call projected_adjoint_source(rows, jets, mu, o, values, seed, owner, rhs)
       return
    end if
    is_fixed = rows % fixed_indicator()
    deg      = rows % degrees
    npts     = size(rows % at)
    allocate(lower(o), partial(0:o), source=0.0_dp)
    allocate(source(size(rhs)), source=0.0_dp)
    design = design_terms(values, seed, o, o + 1)
    allocate(q(0:deg - 1))
    do j = 1, size(rows % rules)
       ! the components the rule reads: its partial in any other is zero
       call rows % rules(j) % read_components(reads)
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
          do m = 1, size(reads)
             c = reads(m)
             call q(c) % set_coefficient(2**o, 1.0_dp)
             r = rows % rules(j) % at_instant(q, design, rows % point_positions(:, p))
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

    if (allocated(rows % projected)) then
       total = projected_design_products(rows, jets, mu, o, values, seed, owner)
       return
    end if
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
             call rows % rules(j) % gradient_at(x, values, value, g, rows % point_positions(:, p), g_design=gd)
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
             r = rows % rules(j) % at_instant(q, design, rows % point_positions(:, p))
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
  ! THE PROJECTED SOURCES AND PRODUCTS, as the three above with the
  ! point replaced by the site: the tuple's jets along the design at a
  ! site formed by the basis from the nodes' columns; the m-th
  ! derivative of the rule's value tested into the rows of the
  ! instance's nodes; the multiplier at a site, of the rule's equation,
  ! the sum of the rows' multipliers times the weight and the test
  ! function; and the partial in a component of the site's tuple
  ! returned to the nodes' components by the basis that formed it. An
  ! instance belongs to the image of its first node, and the sums are
  ! taken over the images once.
  !===================================================================!

  subroutine projected_design_source(rows, jets, m, values, seed, owner, b)

    type(residual_operator), intent(in)    :: rows
    real(dp)               , intent(in)    :: jets(:,0:)
    integer                , intent(in)    :: m
    real(dp)               , intent(in)    :: values(:), seed(:)
    integer                , intent(in)    :: owner(:)
    real(dp)               , intent(inout) :: b(:)

    type(derivative_terms), allocatable :: q(:), design(:)
    type(derivative_terms) :: r
    logical , allocatable :: is_fixed(:)
    real(dp), allocatable :: source(:), formed(:,:)
    integer :: width, nq, j, f, d, n, e, qp, s, c, k, i, p, row, o
    real(dp) :: value

    is_fixed = rows % fixed_indicator()
    width    = size(rows % projected % formed, 2)
    design   = design_terms(values, seed, m, m)
    allocate(q(0:width - 1))
    allocate(source(size(b)), source=0.0_dp)
    allocate(formed(width * rows % sites % num_vertices(), 0:m - 1))
    do o = 0, m - 1
       formed(:, o) = rows % formed_at_sites(jets(:, o))
    end do
    associate (pr => rows % projected, el => rows % projected % elements)
      nq = el % num_quadrature_points()
      do j = 1, size(rows % rules)
         f = pr % tested(1, j)
         d = pr % tested(2, j)
         do n = 1, size(pr % element_cell)
            if (owner(abs(pr % element_point(1, n))) /= this_image()) cycle
            e = pr % element_cell(n)
            do qp = 1, nq
               s = (n - 1) * nq + qp
               if (.not. tested_at(el, d, qp, e)) cycle
               do c = 0, width - 1
                  q(c) = derivative_terms(formed((s - 1) * width + c + 1, 0), m)
                  do k = 1, m - 1
                     call q(c) % set_symmetric(k, formed((s - 1) * width + c + 1, k))
                  end do
               end do
               r     = rows % rules(j) % at_instant(q, design, rows % site_positions(:, s))
               value = -mixed_partial(r)
               do i = 1, size(pr % element_point, 1)
                  p = pr % element_point(i, n)
                  if (p < 1) cycle
                  do k = 0, el % num_coefficients() - 1
                     row = rows % at(p) + pr % coefficient_row(k, f) + 1
                     if (is_fixed(row)) cycle
                     source(row) = source(row) + value * el % test(d, k, i, qp, e)
                  end do
               end do
            end do
         end do
      end do
    end associate
    if (any(owner /= this_image())) call co_sum(source)
    b = b + source

  end subroutine projected_design_source

  subroutine projected_adjoint_source(rows, jets, mu, o, values, seed, owner, rhs)

    type(residual_operator), intent(in)    :: rows
    real(dp)               , intent(in)    :: jets(:,0:), mu(:,0:)
    integer                , intent(in)    :: o
    real(dp)               , intent(in)    :: values(:), seed(:)
    integer                , intent(in)    :: owner(:)
    real(dp)               , intent(inout) :: rhs(:)

    type(derivative_terms), allocatable :: q(:), design(:)
    type(derivative_terms) :: r, along
    logical , allocatable :: is_fixed(:)
    real(dp), allocatable :: source(:), formed(:,:), lambda(:), partial(:), lower(:)
    integer , allocatable :: reads(:)
    integer :: width, nq, j, f, d, n, e, qp, s, c, k, i, p, row, i2, p2, k2, offset, dd, m
    real(dp) :: value, tau

    is_fixed = rows % fixed_indicator()
    width    = size(rows % projected % formed, 2)
    design   = design_terms(values, seed, o, o + 1)
    allocate(q(0:width - 1), lambda(0:o - 1), partial(0:o), lower(o))
    allocate(source(size(rhs)), source=0.0_dp)
    allocate(formed(width * rows % sites % num_vertices(), 0:o))
    do k = 0, o
       formed(:, k) = rows % formed_at_sites(jets(:, k))
    end do
    associate (pr => rows % projected, el => rows % projected % elements)
      nq = el % num_quadrature_points()
      do j = 1, size(rows % rules)
         f = pr % tested(1, j)
         d = pr % tested(2, j)
         ! the components the rule reads: its partial in any other is zero
         call rows % rules(j) % read_components(reads)
         do n = 1, size(pr % element_cell)
            if (owner(abs(pr % element_point(1, n))) /= this_image()) cycle
            e = pr % element_cell(n)
            do qp = 1, nq
               s = (n - 1) * nq + qp
               if (.not. tested_at(el, d, qp, e)) cycle
               ! the multiplier of the rule at the site, to the order
               ! below o, the o-th left out
               lambda = 0.0_dp
               do i = 1, size(pr % element_point, 1)
                  p = pr % element_point(i, n)
                  if (p < 1) cycle
                  do k = 0, el % num_coefficients() - 1
                     row = rows % at(p) + pr % coefficient_row(k, f) + 1
                     if (is_fixed(row)) cycle
                     tau    = el % test(d, k, i, qp, e)
                     lambda = lambda + tau * mu(row, 0:o - 1)
                  end do
               end do
               lower(1:o - 1) = lambda(1:o - 1)
               lower(o)       = 0.0_dp
               along = symmetric_terms(lambda(0), lower, o)
               do c = 0, width - 1
                  q(c) = symmetric_terms(formed((s - 1) * width + c + 1, 0), &
                       & [(formed((s - 1) * width + c + 1, k), k = 1, o)], o + 1)
               end do
               do m = 1, size(reads)
                  c = reads(m)
                  call q(c) % set_coefficient(2**o, 1.0_dp)
                  r = rows % rules(j) % at_instant(q, design, rows % site_positions(:, s))
                  call q(c) % set_coefficient(2**o, 0.0_dp)
                  do k = 0, o
                     partial(k) = coefficient(r, 2**k - 1 + 2**o)
                  end do
                  value = -mixed_partial(along * symmetric_terms(partial(0), partial(1:o), o))
                  dd = pr % formed(0, c + 1)
                  do i2 = 1, size(pr % element_point, 1)
                     p2 = abs(pr % element_point(i2, n))
                     do k2 = 0, el % num_coefficients() - 1
                        offset = pr % formed(k2 + 1, c + 1)
                        if (offset < 0) cycle
                        source(rows % at(p2) + offset + 1) = source(rows % at(p2) + offset + 1) &
                             & + el % basis(dd, k2, i2, qp, e) * value
                     end do
                  end do
               end do
            end do
         end do
      end do
    end associate
    if (any(owner /= this_image())) call co_sum(source)
    rhs = rhs + source

  end subroutine projected_adjoint_source

  function projected_design_products(rows, jets, mu, o, values, seed, owner) result(total)

    type(residual_operator), intent(in) :: rows
    real(dp)               , intent(in) :: jets(:,0:), mu(:,0:)
    integer                , intent(in) :: o
    real(dp)               , intent(in) :: values(:), seed(:)
    integer                , intent(in) :: owner(:)
    real(dp), allocatable :: total(:)

    type(derivative_terms), allocatable :: q(:), design(:)
    type(derivative_terms) :: r
    logical , allocatable :: is_fixed(:)
    real(dp), allocatable :: formed(:,:), lambda(:), partial(:), g(:), gd(:)
    integer :: width, nq, j, f, d, n, e, qp, s, c, k, i, p, row, l
    real(dp) :: value, tau

    is_fixed = rows % fixed_indicator()
    width    = size(rows % projected % formed, 2)
    allocate(total(size(values)), source=0.0_dp)
    allocate(q(0:width - 1), lambda(0:o), partial(0:o), g(0:width - 1), gd(size(values)))
    allocate(formed(width * rows % sites % num_vertices(), 0:o))
    do k = 0, o
       formed(:, k) = rows % formed_at_sites(jets(:, k))
    end do
    associate (pr => rows % projected, el => rows % projected % elements)
      nq = el % num_quadrature_points()
      do l = 1, merge(1, size(values), o == 0)
         if (o > 0) then
            design = design_terms(values, seed, o, o + 1)
            call design(l) % set_coefficient(2**o, 1.0_dp)
         end if
         do j = 1, size(rows % rules)
            f = pr % tested(1, j)
            d = pr % tested(2, j)
            do n = 1, size(pr % element_cell)
               if (owner(abs(pr % element_point(1, n))) /= this_image()) cycle
               e = pr % element_cell(n)
               do qp = 1, nq
                  s = (n - 1) * nq + qp
                  if (.not. tested_at(el, d, qp, e)) cycle
                  lambda = 0.0_dp
                  do i = 1, size(pr % element_point, 1)
                     p = pr % element_point(i, n)
                     if (p < 1) cycle
                     do k = 0, el % num_coefficients() - 1
                        row = rows % at(p) + pr % coefficient_row(k, f) + 1
                        if (is_fixed(row)) cycle
                        tau    = el % test(d, k, i, qp, e)
                        lambda = lambda + tau * mu(row, 0:o)
                     end do
                  end do
                  if (o == 0) then
                     call rows % rules(j) % gradient_at(formed((s - 1) * width + 1:s * width, 0), values, value, g, &
                          & rows % site_positions(:, s), g_design=gd)
                     total = total + lambda(0) * gd
                     cycle
                  end if
                  do c = 0, width - 1
                     q(c) = symmetric_terms(formed((s - 1) * width + c + 1, 0), &
                          & [(formed((s - 1) * width + c + 1, k), k = 1, o)], o + 1)
                  end do
                  r = rows % rules(j) % at_instant(q, design, rows % site_positions(:, s))
                  do k = 0, o
                     partial(k) = coefficient(r, 2**k - 1 + 2**o)
                  end do
                  total(l) = total(l) + mixed_partial(symmetric_terms(lambda(0), lambda(1:o), o) &
                       & * symmetric_terms(partial(0), partial(1:o), o))
               end do
            end do
         end do
      end do
    end associate
    if (any(owner /= this_image())) call co_sum(total)

  end function projected_design_products

  !===================================================================!
  ! The condition on the design factor at the design values: given =
  ! -g(nu_0) and slope(i) = dg/dnu_i, so that nu_k - nu_k0 = 0 gives
  ! zero and the k-th unit vector. Each slope by a jet of one
  ! direction, seeded on nu_i: a jet of every coordinate's direction
  ! stores the mixed derivatives of every subset of them, 2^48 terms
  ! for the 48 strips of the rotor.
  !===================================================================!

  subroutine design_condition(equation, values, given, slope)

    type(continuous_field), intent(in)  :: equation
    real(dp)              , intent(in)  :: values(:)
    real(dp)              , intent(out) :: given
    real(dp), allocatable , intent(out) :: slope(:)

    type(expression) :: g
    type(derivative_terms), allocatable :: q(:), nu(:)
    type(derivative_terms) :: r
    integer :: i, j, nd

    g  = equation % graph(1)
    nd = size(values)
    allocate(nu(nd), slope(nd))
    allocate(q(0:g % num_components() - 1), source=derivative_terms(0.0_dp, 1))
    given = 0.0_dp
    do i = 1, max(nd, 1)
       do j = 1, nd
          nu(j) = derivative_terms(values(j), 1)
       end do
       if (nd > 0) call nu(i) % set_direction(1, 1.0_dp)
       r     = g % at_instant(q, nu)
       given = -coefficient(r, 0)
       if (nd > 0) slope(i) = coefficient(r, 1)
    end do

  end subroutine design_condition

  !===================================================================!
  ! The point manifold as a block: one moment, one cell, no rows along
  ! any coordinate; the equations alone.
  !===================================================================!

  subroutine point_rows(this, tuple, blk)

    class(discrete_residual), intent(in)  :: this
    real(dp)                , intent(in)  :: tuple(:,:,:)
    type(chain_block)       , intent(out) :: blk

    type(stencil) :: relations
    integer , allocatable :: primary(:)
    integer :: stride, nf, f

    stride = this % rule % num_components()
    nf     = this % manifold % num_unknowns
    allocate(primary(nf))
    do f = 1, nf
       primary(f) = this % rule % offset_of_field(f)
    end do
    relations = stencil([integer ::], [integer ::], [real(dp) ::], spread(0.0_dp, 1, stride), 'no relation')
    blk % rows = residual_operator(relations, this % rule, [0], stride, stride, primary, [integer ::], [real(dp) ::])
    blk % q    = tuple(:, 1, 1)
    blk % instant_moment = [1]

  end subroutine point_rows

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
    type(exchange)     :: whole, part
    type(residual_operator) :: sub
    real(dp), allocatable :: local(:)
    real(dp) :: achieved

    ! over several images each image solves its part of the rows, from
    ! the partitioner's cut of the residual's coupling, and the parts'
    ! owned entries are assembled whole; the line of an open line, and
    ! of a law projected on the cells of a line, is one block of the
    ! sweeps, which a cut would divide
    if (num_images() > 1 .and. this % manifold % with_space .and. .not. this % whole_line()) then
       call require_images()
       whole = exchange(this % owners_of_unknowns(unknowns, stride, npts), width=stride)
       part  = whole % part(rows % coupling())
       sub   = rows % part_of(part)
       solver = this % solver_for(sub, size(part % whole_index), stride, size(sub % at), part)
       local = part % restriction(q)
       call solver % solve(spread(0.0_dp, 1, size(local)), local, achieved)
       q = part % assembled(local)
    else
       solver = this % solver_for(rows, unknowns, stride, npts)
       call solver % solve(spread(0.0_dp, 1, unknowns), q, achieved)
    end if
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

  function solver_for(this, rows, unknowns, stride, npts, part) result(solver)

    class(discrete_residual), intent(in) :: this
    type(residual_operator) , intent(in) :: rows
    integer                 , intent(in) :: unknowns, stride, npts
    ! over several images, the part the rows are of
    type(exchange)          , intent(in), optional :: part
    type(newton) :: solver

    type(dense_direct) :: direct
    type(gmres)        :: krylov
    type(gauss_seidel) :: sweeps
    type(elimination)  :: reduced
    type(typed_field_domain) :: designs
    type(stored_field) :: design

    if (this % manifold % with_space) then
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
       ! on an open line the conditions of one end occupy rows at the
       ! other, and a law projected on the cells of a line couples every
       ! node of a moment through its elements: a block is the whole
       ! line at one moment
       if (this % whole_line()) sweeps % block_width = stride * this % points % num_cells()
       sweeps % max_iterations = 4
       krylov % restart        = 60
       krylov % tolerance      = 1.0e-3_dp
       allocate(krylov % preconditioner, source=sweeps)
       allocate(solver % inner, source=krylov)
       ! a build of one image under a launcher of several would run
       ! the whole problem once per process, and is refused
       call require_images()
       if (present(part)) solver % distribution = part
       ! a march begun from the instant before is solved to the rounding
       ! of its residual, as under a direct inner solve: its first
       ! residual, of components the instant before leaves inconsistent
       ! with the family, is not the law's own scale
       if (this % continued) solver % tolerance = 1.0e-13_dp
       ! WITH UNKNOWNS OF TIME ALONE, stored at every point, the tying
       ! rows are eliminated and the rest - the laws of the unknowns of
       ! time alone at the first point, those of the fields at every
       ! node - solved directly: Newton converges to round-off
       if (allocated(rows % tying)) then
          reduced % eliminated = rows % tying
          ! its rows of positions, rotations and reactions equilibrated
          direct % equilibrated = .true.
          allocate(reduced % inner, source=direct)
          deallocate(solver % inner)
          allocate(solver % inner, source=reduced)
          solver % tolerance = 1.0e-13_dp
       end if
    else
       ! a direct inner solve converges Newton quadratically to
       ! round-off at the cost of a step, and the sensitivities read
       ! the converged state: the tolerance is tight here
       allocate(solver % inner, source=direct)
       solver % tolerance = 1.0e-13_dp
    end if

    designs = rows % design_fields()
    design  = designs % design(design_at_points(this % points % expansion % values(), rows % designs, npts))
    ! the coupling of the unknowns, the summed rows excluded, is what
    ! the inner minimizer's sweep is coloured by; a direct solve after
    ! the elimination has no sweep
    if (allocated(rows % tying)) then
       call solver % state(rows, rows % unknown_graph(), rows % unknown_domain(), unknowns, stored_inputs=[design])
    else
       call solver % state(rows, rows % unknown_graph(), rows % unknown_domain(), unknowns, stored_inputs=[design], &
            & coupling=rows % coupling())
    end if

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

    ! without a region nothing is shared: every image owns every cell
    ! of its own, and no sum over the images is taken
    allocate(cell_owner(this % points % num_cells()), source=this_image())
    if (num_images() == 1 .or. .not. this % manifold % with_space) return
    cut = partitioner(PARTITION_BREADTH_FIRST, num_images())
    call cut % assign_owners(this % points % cells, cell_owner)

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
    if (this % manifold % with_space .or. num_images() == 1 .or. this % overlapped) return
    expansion_owned = mod(d - 1, num_images()) + 1 == this_image()

  end function expansion_owned

  ! the owner of every point of a block, moment outer and cell inner
  function point_owners(this, npts) result(owner)

    class(discrete_residual), intent(in) :: this
    integer                 , intent(in) :: npts
    integer, allocatable :: owner(:)

    integer, allocatable :: cell_owner(:)
    integer :: ncells, p, m

    cell_owner = this % cell_owners()
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
    cell_owner = this % cell_owners()

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
  ! THE PROJECTED FORM of the physics on the cells, for a method that
  ! tests the law against a basis. Each equation in flux form, each
  ! group of its terms paired with a multiplier of its own and tested
  ! with the derivative of the basis the group names: the Lagrangian of
  ! the projection. The coefficients a node stores reach one order
  ! below the highest derivative along space that the flux form puts
  ! on the unknown or on the multiplier, and their layout becomes the
  ! rule's; the basis's jet reaches that highest order. The equation
  ! of an unknown of time alone is one group, tested with the weight
  ! alone at the first point of each moment: its terms of time over the
  ! measure of the region, so that they sum to themselves, and the
  ! integrand of its integral over the region as it is. Invalid input:
  ! a representation law along space, which the basis states.
  !===================================================================!

  subroutine projected_form(image, physics, space, points, defined)

    type(discrete_residual)        , intent(inout) :: image
    type(residual_term)            , intent(in)    :: physics
    class(derivative_approximation), intent(in)    :: space
    type(discrete_manifold)        , intent(in)    :: points
    integer                        , intent(in)    :: defined(:,:)

    type(expression), allocatable :: flux(:), groups(:), over_region(:)
    type(expression) :: law, region_part
    integer, allocatable :: coordinate(:), order(:), group(:,:)
    logical, allocatable :: with_region(:)
    integer :: first, last, j, g, m, c, top, o
    real(dp) :: volume
    logical :: differentiated, reads_region

    first = merge(FIRST_COORDINATE + 1, FIRST_COORDINATE, image % manifold % with_time)
    last  = first + points % dimension - 1
    if (any(defined(3, :) >= first)) then
       error stop 'operation_residual: under a projection on the cells the basis states the representation &
            &laws along space; none is stated'
    end if
    allocate(group(3, 0), groups(0), over_region(0), with_region(0))
    m   = 0
    top = 0
    do j = 1, size(physics % equation)
       ! the equation of an unknown of time alone, one group without
       ! its flux form: its terms of time, and the integrand of its
       ! integral over the region
       associate (equation => physics % equation(j))
         region_part  = constant(0.0_dp)
         reads_region = .false.
         if (image % manifold % of_time_alone(j)) then
            if (allocated(equation % integrand)) then
               flux         = [jet_form(equation % graph(1))]
               region_part  = jet_form(equation % integrand(1))
               reads_region = .true.
            else if (equation % is_integral() .and. equation % integrated_part == SPACE_FACTOR) then
               flux         = [constant(0.0_dp)]
               region_part  = jet_form(equation % graph(1))
               reads_region = .true.
            else
               flux = [jet_form(equation % graph(1))]
            end if
            coordinate = [first]
            order      = [0]
         else
            if (allocated(equation % integrand)) then
               error stop 'operation_residual: an integral over the region is a term of the equation of an &
                    &unknown of time alone'
            end if
            call flux_form(equation % graph(1), first, coordinate, order, flux)
         end if
       end associate
       do g = 1, size(flux)
          m           = m + 1
          groups      = [groups, flux(g)]
          over_region = [over_region, region_part]
          with_region = [with_region, reads_region]
          group       = reshape([group, j, coordinate(g), order(g)], [3, m])
          top         = max(top, order(g))
          do c = first, last
             top = max(top, flux(g) % highest_degree_along(c), region_part % highest_degree_along(c))
          end do
       end do
    end do
    image % elements = space % elements(points % cells, max(top - 1, 0), top)
    volume = sum(image % elements % weight)
    law = constant(0.0_dp)
    allocate(image % tested(2, m))
    do g = 1, m
       image % tested(1, g) = group(1, g)
       image % tested(2, g) = 0
       if (group(3, g) > 0) image % tested(2, g) = (group(2, g) - first) * top + group(3, g)
       if (image % manifold % of_time_alone(group(1, g))) then
          ! tested with the weight alone, at the first point of each
          ! moment, over the measure of the region: its terms of time
          ! alone once, an integral over the region as it is; with no
          ! integral over the region, tested at one site of each moment
          ! with one, its terms evaluated there alone
          if (with_region(g)) then
             image % tested(2, g) = ubound(image % elements % test, 1) + 1
             groups(g) = groups(g) * constant(1.0_dp / volume) + over_region(g)
          else
             image % tested(2, g) = ubound(image % elements % test, 1) + 2
          end if
       end if
       law = law + multiplier(g) * groups(g)
    end do
    image % flux = law
    image % rule = coefficient_form(law, first, last, max(top - 1, 0))

    ! AT THE FIRST INSTANT an equation of time alone that reads no
    ! component of top degree - a joint's condition - is stated by its
    ! derivative along time of the least order that reads one, as on
    ! a manifold of time alone (initial_form); so too at the explicit
    ! moments of a block
    if (.not. image % manifold % with_time) return
    differentiated = .false.
    law = constant(0.0_dp)
    do g = 1, m
       if (image % manifold % of_time_alone(group(1, g))) then
          do o = 1, top_degree(image % rule, image % manifold % num_unknowns)
             if (reads_top_degree(groups(g))) exit
             groups(g) = total_derivative(groups(g), FIRST_COORDINATE)
             differentiated = .true.
          end do
       end if
       law = law + multiplier(g) * groups(g)
    end do
    if (differentiated) image % initial_flux = law

  contains

    logical function reads_top_degree(x)
      type(expression), intent(in) :: x
      integer :: f
      reads_top_degree = .false.
      do f = 1, image % manifold % num_unknowns
         if (x % highest_order_read(f) == image % rule % degree_of_field(f)) reads_top_degree = .true.
      end do
    end function reads_top_degree

  end subroutine projected_form

  !===================================================================!
  ! THE LAWS AT THE FIRST INSTANT. There the initial data fix every
  ! component of a field below its top degree along time, so a law
  ! that reads no component of top degree - a joint's condition on
  ! the positions - is satisfied by the data and determines nothing:
  ! the reaction paired with it has no row. Such a law is stated at
  ! the first instant by its derivative along time of the least order
  ! that reads a component of top degree, the condition hidden in its
  ! holding at every instant, which determines the reaction there. The
  ! rule at the first instant is the rule with those laws
  ! differentiated, on the rule's domain; it is not formed when no law
  ! is differentiated. It is the rule too at every explicit moment of
  ! a block - the instant a step of a staged family reaches - where the
  ! family's quadrature of the stages fixes every component below the
  ! top degree as the initial data fix them at the first instant.
  !===================================================================!

  subroutine initial_form(image, physics)

    type(discrete_residual), intent(inout) :: image
    type(residual_term)    , intent(in)    :: physics

    type(expression), allocatable :: laws(:)
    integer :: j, o, nf
    logical :: differentiated

    nf = image % manifold % num_unknowns
    differentiated = .false.
    allocate(laws(size(physics % equation)))
    do j = 1, size(physics % equation)
       laws(j) = jet_form(physics % equation(j) % graph(1))
       do o = 1, top_degree(image % rule, nf)
          if (reads_top_degree(laws(j))) exit
          laws(j) = total_derivative(laws(j), FIRST_COORDINATE)
          differentiated = .true.
       end do
    end do
    if (differentiated) image % initial_rule = on_domain(paired(laws), image % rule)

  contains

    logical function reads_top_degree(x)
      type(expression), intent(in) :: x
      integer :: f
      reads_top_degree = .false.
      do f = 1, nf
         if (x % highest_order_read(f) == image % rule % degree_of_field(f)) reads_top_degree = .true.
      end do
    end function reads_top_degree

  end subroutine initial_form

  !===================================================================!
  ! THE COMPONENT OF A NODE'S TUPLE storing the coefficient k of the
  ! field f at the order a along time: the value's own order for the
  ! coefficient zero, and otherwise the mixed component of the order
  ! of the coefficient's derivative along its axis. Without time the
  ! order a is zero.
  !===================================================================!

  pure integer function coefficient_component(this, f, k, a) result(at)

    class(discrete_residual), intent(in) :: this
    integer                 , intent(in) :: f, k, a

    integer :: first, axis, o

    if (k == 0) then
       at = this % rule % offset_of_field(f) + a
       return
    end if
    first = merge(FIRST_COORDINATE + 1, FIRST_COORDINATE, this % manifold % with_time)
    axis  = (k - 1) / this % elements % order + 1
    o     = k - (axis - 1) * this % elements % order
    if (this % manifold % with_time) then
       at = this % rule % component_at(first + axis - 1, o, f, a)
    else
       at = this % rule % component_at(first + axis - 1, o, f)
    end if

  end function coefficient_component

  !===================================================================!
  ! THE FIELD AND THE ORDER ALONG TIME of a component of a node's
  ! tuple that stores no derivative along space: the condition on a
  ! time face reads and determines such components alone.
  !===================================================================!

  pure subroutine stored_order(this, component, field, order)

    class(discrete_residual), intent(in)  :: this
    integer                 , intent(in)  :: component
    integer                 , intent(out) :: field, order

    integer :: f

    field = 0
    order = 0
    do f = 1, this % manifold % num_unknowns
       if (component >= this % rule % offset_of_field(f) .and. &
            & component <= this % rule % offset_of_field(f) + this % rule % degree_of_field(f)) then
          field = f
          order = component - this % rule % offset_of_field(f)
       end if
    end do
    if (field == 0) then
       error stop 'operation_residual: a condition on a time face reads the values of the unknowns and their &
            &derivatives along time'
    end if

  end subroutine stored_order

  !===================================================================!
  ! HOW THE TUPLE THE FLUX FORM READS IS FORMED at a quadrature point:
  ! for each of its components, formed(0, c) the derivative of the
  ! basis - zero the value, (a - 1) d + o the order o along the axis a,
  ! d the basis's highest order - and formed(k + 1, c) the component of
  ! a node's tuple that the coefficient k contributes, of the same
  ! field and the same order along time.
  !===================================================================!

  subroutine formed_components(this, formed)

    class(discrete_residual), intent(in)  :: this
    integer, allocatable    , intent(out) :: formed(:,:)

    integer :: f, a, c, cc, o, k, axis, first, top, nk

    nk    = this % elements % num_coefficients()
    top   = this % elements % derivatives
    first = merge(FIRST_COORDINATE + 1, FIRST_COORDINATE, this % manifold % with_time)
    allocate(formed(0:nk, this % flux % num_components()), source=-1)
    do f = 1, this % manifold % num_unknowns
       ! along the first coordinate: the orders along time, or without
       ! time the orders along the first axis of space
       do a = 0, this % flux % degree_of_field(f)
          c = this % flux % offset_of_field(f) + a + 1
          if (this % manifold % with_time) then
             formed(0, c) = 0
             do k = 0, nk - 1
                formed(k + 1, c) = this % coefficient_component(f, k, a)
             end do
          else
             formed(0, c) = a
             do k = 0, nk - 1
                formed(k + 1, c) = this % coefficient_component(f, k, 0)
             end do
          end if
       end do
       do cc = FIRST_COORDINATE + 1, this % flux % num_coordinates()
          axis = cc - first + 1
          do a = 0, merge(this % flux % mixed_degree(f), 0, this % manifold % with_time)
             do o = 1, this % flux % degree_along(cc)
                if (this % manifold % with_time) then
                   c = this % flux % component_at(cc, o, f, a) + 1
                else
                   c = this % flux % component_at(cc, o, f) + 1
                end if
                formed(0, c) = (axis - 1) * top + o
                do k = 0, nk - 1
                   formed(k + 1, c) = this % coefficient_component(f, k, a)
                end do
             end do
          end do
       end do
    end do

  end subroutine formed_components

  !===================================================================!
  ! THE QUADRATURE NODES OF A PROJECTED SOLUTION: every site of the
  ! cells at every instant and at every stage of a staged family's
  ! steps, the site (e - 1) nq + q at the quadrature point q of the
  ! element e; its weight the family's rule in time times the site's in
  ! space, its jet the one the tuples of the element's nodes form there
  ! by the basis, of the layout of a node's tuple. stage_node(k, i, s)
  ! is the node of stage i of the step into instant k at site s; a site
  ! is a node without a point.
  !===================================================================!

  subroutine projected_nodes(this, tuple, stages, solution, stage_node)

    class(discrete_residual), intent(in)    :: this
    real(dp)                , intent(in)    :: tuple(:,:,:,0:), stages(:,:,:,:,0:)
    type(discrete_field)    , intent(inout) :: solution
    integer , allocatable   , intent(out)   :: stage_node(:,:,:)

    type(derivative_terms), allocatable :: dt(:), weight(:)
    integer , allocatable :: instant(:), cell_owner(:), formed(:,:)
    real(dp), allocatable :: position(:)
    integer :: stride, ninst, nsites, nq, nnodes, nn, nb, b, first, last, depth, n, s, smax, k, i, j, e, q, site
    logical :: staged
    real(dp) :: h

    call this % formed_components(formed)
    stride = this % flux % num_components()
    ninst  = this % points % num_instants()
    nq     = this % elements % num_quadrature_points()
    nsites = this % elements % num_cells() * nq
    smax   = 1
    nnodes = ninst * nsites
    nb     = 0
    if (this % with_chain) nb = this % schemes % num_blocks()
    do b = 1, nb
       smax = max(smax, this % schemes % scheme(b) % num_stages())
       call this % block_extent(b, first, last, depth)
       if (marches_by_stages(this % schemes % scheme(b), top_degree(this % rule, this % manifold % num_unknowns) + 1)) &
            & nnodes = nnodes + (last - first) * this % schemes % scheme(b) % num_stages() * nsites
    end do
    allocate(stage_node(lbound(stages, 4):ubound(stages, 4), smax, nsites), source=0)
    if (allocated(solution % node_jet)) then
       deallocate(solution % node_jet, solution % node_weight, solution % node_position, solution % node_point, &
            & solution % node_image)
    end if
    allocate(solution % node_jet(stride, nnodes, this % points % expansion % num_columns()), source=0.0_dp)
    allocate(solution % node_weight(nnodes), source=0.0_dp)
    allocate(solution % node_position(this % points % num_coordinates(), nnodes), source=0.0_dp)
    allocate(solution % node_point(nnodes), source=0)
    allocate(solution % node_image(nnodes), source=this_image())
    solution % node_rule = this % flux
    cell_owner = this % cell_owners()

    ! the sites at every instant
    do k = 1, ninst
       do e = 1, this % elements % num_cells()
          do q = 1, nq
             site = (e - 1) * nq + q
             nn   = (k - 1) * nsites + site
             solution % node_image(nn) = cell_owner(this % elements % node(1, e))
             position = this % points % position(this % points % point_of(k, this % elements % node(1, e)))
             position(size(position) - this % elements % dimension + 1:) = this % elements % position(:, q, e)
             solution % node_position(:, nn) = position
             solution % node_jet(:, nn, :) = this % site_jet(formed, e, q, tuple(:, :, k, :))
          end do
       end do
    end do
    if (.not. this % with_chain) then
       do k = 1, ninst
          do e = 1, this % elements % num_cells()
             do q = 1, nq
                solution % node_weight((k - 1) * nsites + (e - 1) * nq + q) = this % elements % weight(q, e)
             end do
          end do
       end do
       return
    end if

    nn = ninst * nsites
    do b = 1, nb
       call this % block_extent(b, first, last, depth)
       n      = last - first + 1
       staged = marches_by_stages(this % schemes % scheme(b), top_degree(this % rule, this % manifold % num_unknowns) + 1)
       associate (scheme => this % schemes % scheme(b))
         if (staged) then
            s = scheme % num_stages()
            do k = 1, n - 1
               h = this % points % step(first + k)
               do i = 1, s
                  do e = 1, this % elements % num_cells()
                     do q = 1, nq
                        site = (e - 1) * nq + q
                        nn   = nn + 1
                        stage_node(first + k, i, site) = nn
                        solution % node_image(nn)  = cell_owner(this % elements % node(1, e))
                        solution % node_weight(nn) = h * scheme % stage_weight(i) * this % elements % weight(q, e)
                        position = this % points % position(this % points % point_of(first + k - 1, &
                             & this % elements % node(1, e)))
                        position(1) = position(1) + scheme % stage_abscissa(i) * h
                        position(size(position) - this % elements % dimension + 1:) = this % elements % position(:, q, e)
                        solution % node_position(:, nn) = position
                        solution % node_jet(:, nn, :) = this % site_jet(formed, e, q, stages(:, :, i, first + k, :))
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
                  do e = 1, this % elements % num_cells()
                     do q = 1, nq
                        site = (e - 1) * nq + q
                        nn   = (first + instant(j) - 2) * nsites + site
                        solution % node_weight(nn) = solution % node_weight(nn) &
                             & + h * value(weight(j)) * this % elements % weight(q, e)
                     end do
                  end do
               end do
            end do
            deallocate(dt)
         end if
       end associate
    end do

  end subroutine projected_nodes

  !===================================================================!
  ! THE JET AT A QUADRATURE POINT of the layout of the flux's tuple,
  ! in every column: each component - of a field, an order along time
  ! and an order along space - the sum over the element's nodes and
  ! their coefficients of the basis's derivative it names times the
  ! nodes' components of the same field and order along time, as
  ! formed_components lays them out. site_gradient is its transpose: a
  ! differential at the quadrature point returned to the nodes' rows
  ! at a moment of a block, subtracted from rhs.
  !===================================================================!

  pure function site_jet(this, formed, e, q, tuples) result(jet)

    class(discrete_residual), intent(in) :: this
    integer                 , intent(in) :: formed(0:,:), e, q
    real(dp)                , intent(in) :: tuples(:,:,0:)
    real(dp) :: jet(size(formed, 2), 0:ubound(tuples, 3))

    integer :: c, i, k

    jet = 0.0_dp
    associate (el => this % elements)
      do c = 1, size(formed, 2)
         if (formed(0, c) < 0) cycle
         do i = 1, el % num_nodes()
            do k = 0, el % num_coefficients() - 1
               jet(c, :) = jet(c, :) + el % basis(formed(0, c), k, i, q, e) * tuples(formed(k + 1, c) + 1, el % node(i, e), :)
            end do
         end do
      end do
    end associate

  end function site_jet

  subroutine site_gradient(this, formed, site, g, blk, m, rhs)

    class(discrete_residual), intent(in)    :: this
    integer                 , intent(in)    :: formed(0:,:), site, m
    real(dp)                , intent(in)    :: g(:)
    type(chain_block)       , intent(in)    :: blk
    real(dp)                , intent(inout) :: rhs(:)

    integer :: c, i, k, e, q, row

    associate (el => this % elements)
      e = (site - 1) / el % num_quadrature_points() + 1
      q = site - (e - 1) * el % num_quadrature_points()
      do c = 1, size(formed, 2)
         if (g(c) == 0.0_dp .or. formed(0, c) < 0) cycle
         do i = 1, el % num_nodes()
            do k = 0, el % num_coefficients() - 1
               row = blk % tuple_at(m, el % node(i, e)) + formed(k + 1, c) + 1
               rhs(row) = rhs(row) - el % basis(formed(0, c), k, i, q, e) * g(c)
            end do
         end do
      end do
    end associate

  end subroutine site_gradient

  !===================================================================!
  ! THE COMPONENT A CONDITION ON A FACE OF THE REGION FIXES under a
  ! projection, within a node's tuple: of the field of equation i of
  ! condition j, the coefficient it determines, at the order along time
  ! of the equations' rows (equation_row less the field's offset).
  !===================================================================!

  integer function fixed_coefficient(this, j, i, field, equation_row) result(component)

    class(discrete_residual), intent(in)  :: this
    integer                 , intent(in)  :: j, i, equation_row(:)
    integer                 , intent(out) :: field

    integer , allocatable :: column_of(:)
    real(dp), allocatable :: slopes(:), position(:)
    real(dp) :: given
    integer  :: along(2), coefficient, first

    first    = merge(FIRST_COORDINATE + 1, FIRST_COORDINATE, this % manifold % with_time)
    position = this % points % position(1)
    call this % space_face_relation(this % condition(j) % equation(i), position, column_of, slopes, given, field, along)
    coefficient = 0
    if (along(2) > 0) coefficient = (along(1) - first) * this % elements % order + along(2)
    component = this % coefficient_component(field, coefficient, equation_row(field) - this % rule % offset_of_field(field))

  end function fixed_coefficient

  !===================================================================!
  ! THE m-th DERIVATIVE ALONG THE DESIGN of the data of a condition on
  ! a time face, projected: at the row of each node's coefficient, the
  ! integral over the cells around the node of the coefficient's basis
  ! function times the derivative of the value the condition states,
  ! written into rhs at the row the condition occupies in the block.
  !===================================================================!

  subroutine projected_datum(this, equation, instant, m, values, seed, blk, rhs)

    class(discrete_residual), intent(in)    :: this
    type(continuous_field)  , intent(in)    :: equation
    integer                 , intent(in)    :: instant, m
    real(dp)                , intent(in)    :: values(:), seed(:)
    type(chain_block)       , intent(in)    :: blk
    real(dp)                , intent(inout) :: rhs(:)

    integer , allocatable :: columns(:)
    real(dp), allocatable :: slopes(:), position(:)
    real(dp) :: given, stated
    integer  :: e, q, i0, k0, row, row_field, row_order, moment, r
    logical , allocatable :: written(:)

    moment = this % face_row_moment(instant, blk % first, blk % last, blk % depth, blk % instant_moment)
    allocate(written(size(rhs)), source=.false.)
    associate (el => this % elements)
      do e = 1, el % num_cells()
         do q = 1, el % num_quadrature_points()
            position = this % points % position(this % points % point_of(instant, el % node(1, e)))
            position(size(position) - el % dimension + 1:) = el % position(:, q, e)
            call this % face_relation(equation, values, position, columns, slopes, given, row)
            call this % stored_order(row, row_field, row_order)
            stated = stated_derivative(equation, m, values, seed, position)
            do i0 = 1, el % num_nodes()
               do k0 = 0, el % num_coefficients() - 1
                  r = blk % tuple_at(moment, el % node(i0, e)) + this % coefficient_component(row_field, k0, row_order) + 1
                  if (.not. written(r)) rhs(r) = 0.0_dp
                  written(r) = .true.
                  rhs(r) = rhs(r) + el % test(0, k0, i0, q, e) * stated
               end do
            end do
         end do
      end do
    end associate

  end subroutine projected_datum

  ! the derivative of the basis a coefficient names: the value for the
  ! coefficient zero, the same order along the same axis otherwise
  pure integer function coefficient_derivative(this, k) result(d)
    class(discrete_residual), intent(in) :: this
    integer                 , intent(in) :: k
    integer :: axis, o
    d = 0
    if (k == 0) return
    axis = (k - 1) / this % elements % order + 1
    o    = k - (axis - 1) * this % elements % order
    d    = (axis - 1) * this % elements % derivatives + o
  end function coefficient_derivative

  !===================================================================!
  ! THE POINTS ON A FACE OF THE REGION: those whose coordinate along
  ! the face's axis is the face's value, to the rounding of the
  ! generated positions. A face with no point on it is invalid input.
  !===================================================================!

  function points_on_face(this, j) result(on_face)

    class(discrete_residual), intent(in) :: this
    integer                 , intent(in) :: j
    integer, allocatable :: on_face(:)

    real(dp), allocatable :: distance(:)
    real(dp) :: tolerance
    integer :: ncells, d, axis, c

    ncells = this % points % num_cells()
    d      = this % points % dimension
    associate (face => this % condition(j) % manifold)
      axis      = face % face_axis - merge(1, 0, this % manifold % with_time)
      distance  = abs(this % points % centre(axis:d * ncells:d) - face % face_value)
      tolerance = 4.0_dp * epsilon(1.0_dp) * max(abs(face % face_value), maxval(distance))
      on_face   = pack([(c, c = 1, ncells)], distance <= tolerance)
    end associate
    if (size(on_face) == 0) then
       error stop 'operation_residual: a condition on a face of the region requires a point on the face'
    end if

  end function points_on_face

  !===================================================================!
  ! THE ESTIMATE'S DERIVATIVES ALONG TIME at its instants, those of the
  ! parabola through its values at three instants about each: the
  ! first and the second derivative of every field, of its value's
  ! coefficient; above the second, none. A march begun from the
  ! instant before (continued) begins its first instant from them,
  ! where the initial data state the rates, rather than from zero
  ! rates, which a law quadratic in them reaches by a damped
  ! iteration. tuple(:, :, k) is the instant first + k - 1.
  !===================================================================!

  subroutine estimated_rates(this, tuple, first)

    class(discrete_residual), intent(in)    :: this
    real(dp)                , intent(inout) :: tuple(:,:,:)
    integer                 , intent(in)    :: first

    real(dp) :: s(3), w1(3), w2(3), tk
    integer :: n, k, j, f, o, i

    n = size(tuple, 3)
    if (n < 3) return
    do k = 1, n
       j  = min(max(k - 1, 1), n - 2)
       s  = [(this % points % instant(first + j + i - 1), i = 0, 2)]
       tk = this % points % instant(first + k - 1)
       w1 = [((tk - s(2)) + (tk - s(3))) / ((s(1) - s(2)) * (s(1) - s(3))), &
            & ((tk - s(1)) + (tk - s(3))) / ((s(2) - s(1)) * (s(2) - s(3))), &
            & ((tk - s(1)) + (tk - s(2))) / ((s(3) - s(1)) * (s(3) - s(2)))]
       w2 = [2.0_dp / ((s(1) - s(2)) * (s(1) - s(3))), 2.0_dp / ((s(2) - s(1)) * (s(2) - s(3))), &
            & 2.0_dp / ((s(3) - s(1)) * (s(3) - s(2)))]
       do f = 1, this % manifold % num_unknowns
          o = this % rule % offset_of_field(f)
          if (this % rule % degree_of_field(f) >= 1) tuple(o + 2, :, k) = matmul(tuple(o + 1, :, j:j + 2), w1)
          if (this % rule % degree_of_field(f) >= 2) tuple(o + 3, :, k) = matmul(tuple(o + 1, :, j:j + 2), w2)
       end do
    end do

  end subroutine estimated_rates

  !===================================================================!
  ! A TUPLE EXTRAPOLATED ALONG TIME by h: of every unknown that is a
  ! function of no coordinate of space and of the second degree in
  ! time at least, the value and each derivative along time below the
  ! top one replaced by the taylor polynomial, at h, of the
  ! derivatives above it and below the top - of the field, or of each
  ! coefficient of the basis of a projection.
  !
  ! The top derivative is not read: at the instant ahead of a step of
  ! a tableau the law alone determines it, from the value and the
  ! rates the weights b form, the error of their stiff modes
  ! multiplied by those modes' stiffness - 2e2 where the acceleration
  ! is 8 on the beam of case 8 on its hub - and the polynomial that
  ! read it began Newton 1e4 times further from the solution than the
  ! tuple itself. Nor is a field of space extrapolated: its stiff
  ! modes' rates times h and their stiffness begin Newton at six times
  ! the residual of the tuple itself, its test relative to that
  ! residual then ends at the edge of the rounding, and the misfit of
  ! case 8 varies by 5e-15 between designs 1e-9 apart where it varies
  ! by 8e-16 from the tuple.
  !===================================================================!

  subroutine extrapolated(this, tuple, h)

    class(discrete_residual), intent(in)    :: this
    real(dp)                , intent(inout) :: tuple(:,:)
    real(dp)                , intent(in)    :: h

    integer, allocatable :: chain(:)
    integer :: f, k, a, top, nk

    nk = 1
    if (allocated(this % elements)) nk = this % elements % num_coefficients()
    do f = 1, this % manifold % num_unknowns
       top = this % rule % degree_of_field(f)
       if (top < 2) cycle
       if (this % manifold % with_space .and. .not. this % manifold % of_time_alone(f)) cycle
       do k = 0, nk - 1
          if (allocated(this % elements)) then
             chain = [(this % coefficient_component(f, k, a) + 1, a = 0, top - 1)]
          else
             chain = [(this % rule % offset_of_field(f) + a + 1, a = 0, top - 1)]
          end if
          call taylor(chain)
       end do
    end do

  contains

    ! the components of the chain, from the value up, each the sum of
    ! those above it times h^i / i!
    subroutine taylor(chain)
      integer, intent(in) :: chain(0:)
      real(dp) :: term(size(tuple, 2)), factor
      integer :: j, i
      do j = 0, ubound(chain, 1) - 1
         term   = tuple(chain(j), :)
         factor = 1.0_dp
         do i = j + 1, ubound(chain, 1)
            factor = factor * h / real(i - j, dp)
            term   = term + factor * tuple(chain(i), :)
         end do
         tuple(chain(j), :) = term
      end do
    end subroutine taylor

  end subroutine extrapolated

  !===================================================================!
  ! THE ROWS OF THE UNKNOWNS OF TIME ALONE in a block. Such an unknown
  ! is stored at every point, its equation stated at the first point
  ! of each moment; at every other point the component the equation
  ! determines equals the one at the point before, and every
  ! coefficient of the basis but the value is zero at every order, so
  ! that the basis forms the unknown itself at every quadrature point.
  ! Appended to the prescribed rows, none where the history or a face
  ! condition takes the row; added lists them.
  !===================================================================!

  subroutine time_alone_rows_of(this, depth, nm, instant_moment, at, primary, fixed_rows, rows, columns, weights, &
       & constants, added)

    class(discrete_residual), intent(in)    :: this
    integer                 , intent(in)    :: depth, nm, instant_moment(:), at(:), primary(:), fixed_rows(:)
    integer , allocatable   , intent(inout) :: rows(:), columns(:)
    real(dp), allocatable   , intent(inout) :: weights(:), constants(:)
    integer , allocatable   , intent(out)   :: added(:)

    logical , allocatable :: taken(:)
    integer , allocatable :: r(:), c(:)
    real(dp), allocatable :: w(:)
    integer :: ncells, nk, f, a, m, p, k, row, component, e, o

    ncells = this % points % num_cells()
    nk     = this % elements % num_coefficients()
    allocate(taken(maxval(at) + this % rule % num_components()), source=.false.)
    taken(fixed_rows) = .true.
    taken(rows)       = .true.
    allocate(r(2 * nm * ncells * nk * (top_degree(this % rule, this % manifold % num_unknowns) + 1) &
         & * count(this % manifold % of_time_alone(1:this % manifold % num_unknowns))))
    allocate(c(size(r)), w(size(r)), added(size(r)))
    e = 0
    o = 0
    do f = 1, this % manifold % num_unknowns
       if (.not. this % manifold % of_time_alone(f)) cycle
       a = primary(f) - this % rule % offset_of_field(f)
       do m = 1, nm
          if (any(instant_moment(1:depth) == m)) cycle
          do p = 1, ncells
             ! the value's coefficient at the order the equation
             ! determines: the one at the point before
             component = this % coefficient_component(f, 0, a)
             row = at((m - 1) * ncells + p) + component + 1
             if (p > 1 .and. .not. taken(row)) then
                call add(row, row, 1.0_dp)
                call add(row, at((m - 1) * ncells + p - 1) + component + 1, -1.0_dp)
                o = o + 1
                added(o) = row
             end if
             ! every other coefficient zero, at every order
             do k = 1, nk - 1
                do a = 0, this % rule % degree_of_field(f)
                   row = at((m - 1) * ncells + p) + this % coefficient_component(f, k, a) + 1
                   if (taken(row)) cycle
                   call add(row, row, 1.0_dp)
                   o = o + 1
                   added(o) = row
                end do
             end do
             a = primary(f) - this % rule % offset_of_field(f)
          end do
       end do
    end do
    rows      = [rows, r(1:e)]
    columns   = [columns, c(1:e)]
    weights   = [weights, w(1:e)]
    constants = [constants, spread(0.0_dp, 1, e)]
    added     = added(1:o)

  contains

    subroutine add(row, column, weight)
      integer , intent(in) :: row, column
      real(dp), intent(in) :: weight
      e    = e + 1
      r(e) = row
      c(e) = column
      w(e) = weight
    end subroutine add

  end subroutine time_alone_rows_of

  !===================================================================!
  ! THE FIRST POINT OF A MOMENT AS A NODE OF EVERY ELEMENT: no basis
  ! function, so that the tuple at a quadrature point is formed as
  ! before, and, at the derivative index past the basis's last, a test
  ! function that is the weight of each quadrature point. The equation
  ! of an unknown of time alone, tested so, is its integral over the
  ! cells, stated at the first point of the moment. At the index after
  ! it, a test function that is one at the first quadrature point of
  ! the first cell and zero elsewhere: an equation tested so is stated
  ! by its value there.
  !===================================================================!

  pure subroutine with_first_point(el)

    type(element_quadrature), intent(inout) :: el

    real(dp), allocatable :: basis(:,:,:,:,:), test(:,:,:,:,:)
    integer , allocatable :: node(:,:)
    integer :: d, nn

    d  = ubound(el % test, 1) + 1
    nn = el % num_nodes()
    allocate(basis(0:d + 1, 0:ubound(el % basis, 2), nn + 1, size(el % basis, 4), size(el % basis, 5)), &
         & source=0.0_dp)
    allocate(test(0:d + 1, 0:ubound(el % basis, 2), nn + 1, size(el % basis, 4), size(el % basis, 5)), &
         & source=0.0_dp)
    basis(0:d - 1, :, 1:nn, :, :) = el % basis
    test(0:d - 1, :, 1:nn, :, :)  = el % test
    test(d, 0, nn + 1, :, :)      = el % weight
    test(d + 1, 0, nn + 1, 1, 1)  = 1.0_dp
    allocate(node(nn + 1, size(el % node, 2)))
    node(1:nn, :)   = el % node
    node(nn + 1, :) = 1
    call move_alloc(basis, el % basis)
    call move_alloc(test, el % test)
    call move_alloc(node, el % node)

  end subroutine with_first_point

  ! whether a test function of the derivative index d is other than
  ! zero at the quadrature point q of the cell e
  pure logical function tested_at(el, d, q, e)
    type(element_quadrature), intent(in) :: el
    integer                 , intent(in) :: d, q, e
    tested_at = any(el % test(d, :, :, q, e) /= 0.0_dp)
  end function tested_at

  !===================================================================!
  ! THE FIXED ROWS of a block: every component at the history instants,
  ! from the tuples already solved.
  !===================================================================!

  subroutine fixed_rows_of(this, first, last, depth, instant_moment, at, lb, tuple, fixed_rows, fixed)

    class(discrete_residual), intent(in)  :: this
    integer                 , intent(in)  :: first, last, depth, instant_moment(:), at(:)
    type(lower_bounds)      , intent(in)  :: lb
    real(dp)                , intent(in)  :: tuple(:,:,lb % instant:)
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

  !===================================================================!
  ! THE FACE ROWS of a block, each condition's affine relation in
  ! place of the row of the component it determines, with its value at
  ! zero state as the row's constant.
  !
  ! On a face of the time factor, the first or the last instant of the
  ! manifold when it is one of the block's own, at every cell: the row
  ! of the component read with coefficient one. A face at a history
  ! instant is fixed by the history already.
  !
  ! On a face of the region, at every moment of the block after its history but the
  ! manifold's first instant, where the initial data and the physics
  ! at every point determine the jet, and but the explicit moments,
  ! where the family's quadrature of the moments before, which satisfy
  ! the condition, determines the values. Each condition takes the row of
  ! the physics at one of the points nearest the face, so the physics
  ! is stated on the points away from it. A condition on the value takes
  ! the physics row at the face's own point. A condition on a space
  ! derivative takes the row of the difference relation of that
  ! derivative at the face's point, and the relation moves to the
  ! physics row of the next point: the same equations, each point's
  ! block keeping an entry on its own components. The condition is its
  ! affine relation in the components at the face's point, with its
  ! value at zero state as the constant. Invalid input: a face with no
  ! point on it, or more conditions on a face than the region has
  ! points.
  !===================================================================!

  subroutine face_rows_of(this, first, last, depth, nm, instant_moment, at, primary, positions, rows_along, &
       & implicit_moment, rows, columns, weights, constants)

    class(discrete_residual), intent(in)  :: this
    integer                 , intent(in)  :: first, last, depth, nm, instant_moment(:), at(:), primary(:)
    real(dp)                , intent(in)  :: positions(:,:)
    type(stencil)           , intent(in)  :: rows_along(:,:)
    logical                 , intent(in)  :: implicit_moment(:)
    integer , allocatable   , intent(out) :: rows(:), columns(:)
    real(dp), allocatable   , intent(out) :: weights(:), constants(:)

    integer , allocatable :: column_of(:), target(:), row_point(:), row_offset(:), field(:), determined(:,:)
    real(dp), allocatable :: slopes(:), w(:)
    real(dp) :: given
    integer :: ncells, j, i, m, k, e, face_point, row, space, f, along(2), bound, entries, instant, base
    logical :: acting
    logical, allocatable :: taken(:), claimed(:)

    ! an equation on a time face enters a column per component at each
    ! cell; one on a face of the region enters at each moment at most
    ! one moved row, the edges into the face's point, and a column per
    ! component. Projected, an equation on a time face enters at each
    ! quadrature point of each element a column per component read,
    ! node and coefficient, for each node and coefficient tested
    ncells = this % points % num_cells()
    bound  = 0
    do j = 1, size(this % condition)
       associate (face => this % condition(j) % manifold, count => size(this % condition(j) % equation))
         if (face % is_time_face()) bound = bound + count * ncells * this % rule % num_components()
         if (face % is_space_face()) bound = bound + count * nm * (1 + ncells + this % rule % num_components())
         if (face % is_time_face() .and. allocated(this % elements)) bound = bound + count &
              & * this % elements % num_cells() * this % elements % num_quadrature_points() &
              & * (this % elements % num_nodes() * this % elements % num_coefficients())**2 &
              & * this % rule % num_components()
       end associate
    end do
    allocate(rows(bound), columns(bound), weights(bound), constants(bound))
    entries = 0
    ! projected, the initial data are projected onto the functions that
    ! satisfy the conditions on the faces of the region: at the
    ! manifold's first instant each such condition takes the row of
    ! the coefficient it fixes in place of the projection's
    allocate(taken(maxval(at) + this % rule % num_components()), source=.false.)
    ! a node on several faces: the condition stated first takes its row
    allocate(claimed(size(taken)), source=.false.)
    if (allocated(this % elements) .and. first == 1 .and. depth == 0) call taken_at_first_instant()
    do j = 1, size(this % condition)
       if (.not. this % condition(j) % manifold % is_time_face()) cycle
       instant = this % face_instant(this % condition(j) % manifold % face_value)
       if (instant < first + depth .or. instant > last) cycle
       m = this % face_row_moment(instant, first, last, depth, instant_moment)
       if (allocated(this % elements)) then
          do i = 1, size(this % condition(j) % equation)
             ! a condition on unknowns of time alone, at every point
             if (of_time_alone(this % condition(j) % equation(i))) then
                do k = 1, ncells
                   call this % face_relation(this % condition(j) % equation(i), this % points % expansion % values(), &
                        & this % points % position(this % points % point_of(instant, k)), column_of, slopes, given, row)
                   base = at((instant_moment(instant - first + 1) - 1) * ncells + k)
                   do e = 1, size(column_of)
                      call add(at((m - 1) * ncells + k) + row + 1, base + column_of(e) + 1, slopes(e), -given)
                   end do
                end do
                cycle
             end if
             call projected_condition(this % condition(j) % equation(i), instant, &
                  & (instant_moment(instant - first + 1) - 1) * ncells, (m - 1) * ncells)
          end do
          cycle
       end if
       do i = 1, size(this % condition(j) % equation)
          do k = 1, ncells
             call this % face_relation(this % condition(j) % equation(i), this % points % expansion % values(), &
                  & this % points % position(this % points % point_of(instant, k)), column_of, slopes, given, row)
             ! the row at the moment the condition occupies, the columns at the face's
             base = at((instant_moment(instant - first + 1) - 1) * ncells + k)
             do e = 1, size(column_of)
                call add(at((m - 1) * ncells + k) + row + 1, base + column_of(e) + 1, slopes(e), -given)
             end do
          end do
       end do
    end do
    do j = 1, size(this % condition)
       if (.not. this % condition(j) % manifold % is_space_face()) cycle
       if (allocated(this % elements)) then
          call coefficient_condition(j)
          cycle
       end if
       call this % space_face_placement(j, primary, face_point, target, row_point, row_offset, field, determined)
       do i = 1, size(this % condition(j) % equation)
          do m = 1, nm
             ! the implicit moments of the block after its history, the
             ! first instant aside: at an explicit moment the family's
             ! quadrature of the moments before, which satisfy the
             ! condition, satisfies it too
             acting = implicit_moment(m)
             do k = 1, depth
                if (instant_moment(k) == m) acting = .false.
             end do
             if (first == 1 .and. instant_moment(1) == m) acting = .false.
             if (.not. acting) cycle
             call this % space_face_relation(this % condition(j) % equation(i), &
                  & positions(:, (m - 1) * ncells + face_point), column_of, slopes, given, f, along)
             row = at((m - 1) * ncells + row_point(i)) + row_offset(i) + 1
             if (determined(2, i) > 0) then
                ! the difference relation of the derivative at the face's
                ! point, moved to the physics row of the target point
                space = determined(1, i) - merge(1, 0, this % manifold % with_time)
                associate (op => rows_along(space, determined(2, i)))
                  call op % weights % real_vector(w)
                  call add(at((m - 1) * ncells + target(i)) + primary(field(i)) + 1, row, 1.0_dp, 0.0_dp)
                  do e = 1, op % pattern % num_edges()
                     if (op % pattern % edge_head(e) /= face_point) cycle
                     call add(at((m - 1) * ncells + target(i)) + primary(field(i)) + 1, &
                          & at((m - 1) * ncells + op % pattern % edge_tail(e)) + this % rule % offset_of_field(field(i)) &
                          & + 1, -w(e), 0.0_dp)
                  end do
                end associate
             end if
             do e = 1, size(column_of)
                call add(row, at((m - 1) * ncells + face_point) + column_of(e) + 1, slopes(e), -given)
             end do
          end do
       end do
    end do
    rows      = rows(1:entries)
    columns   = columns(1:entries)
    weights   = weights(1:entries)
    constants = constants(1:entries)

  contains

    subroutine add(row, column, weight, constant)
      integer , intent(in) :: row, column
      real(dp), intent(in) :: weight, constant
      entries            = entries + 1
      rows(entries)      = row
      columns(entries)   = column
      weights(entries)   = weight
      constants(entries) = constant
    end subroutine add

    ! a condition that reads the unknowns of time alone and no other
    logical function of_time_alone(equation)
      type(continuous_field), intent(in) :: equation
      type(expression) :: form
      integer :: f
      form = equation % graph(1)
      of_time_alone = any(this % manifold % of_time_alone(1:this % manifold % num_unknowns))
      do f = 1, this % manifold % num_unknowns
         if (.not. this % manifold % of_time_alone(f) .and. form % highest_order_read(f) >= 0) of_time_alone = .false.
      end do
    end function of_time_alone

    !-----------------------------------------------------------------!
    ! A CONDITION ON A TIME FACE, PROJECTED: for every node and every
    ! coefficient, the integral over the cells around the node of the
    ! coefficient's basis function times the condition, in place of
    ! the row of that coefficient at the order along time the
    ! condition determines; the condition's jet at each quadrature
    ! point formed by the basis from the coefficients at the face's
    ! instant. The constant of a row is the integral of the basis
    ! function times the condition at zero state, written on each of
    ! the row's entries. columns_at and rows_at offset the points of
    ! the moment of the columns and of the rows, the face's instant
    ! being the manifold's instant.
    !-----------------------------------------------------------------!

    subroutine projected_condition(equation, instant, columns_at, rows_at)

      type(continuous_field), intent(in) :: equation
      integer               , intent(in) :: instant, columns_at, rows_at

      real(dp), allocatable :: position(:), constant_of(:)
      integer  :: e, q, i0, k0, i, k, r, cr, row_field, row_order, column_field, column_order, start, row
      real(dp) :: tau

      allocate(constant_of(maxval(at) + this % rule % num_components()), source=0.0_dp)
      start = entries
      associate (el => this % elements)
        do e = 1, el % num_cells()
           do q = 1, el % num_quadrature_points()
              position = this % points % position(this % points % point_of(instant, el % node(1, e)))
              position(size(position) - el % dimension + 1:) = el % position(:, q, e)
              call this % face_relation(equation, this % points % expansion % values(), position, column_of, &
                   & slopes, given, row)
              call this % stored_order(row, row_field, row_order)
              do i0 = 1, el % num_nodes()
                 do k0 = 0, el % num_coefficients() - 1
                    tau = el % test(0, k0, i0, q, e)
                    r   = at(rows_at + el % node(i0, e)) + this % coefficient_component(row_field, k0, row_order) + 1
                    if (taken(r)) cycle
                    constant_of(r) = constant_of(r) - tau * given
                    do cr = 1, size(column_of)
                       call this % stored_order(column_of(cr), column_field, column_order)
                       do i = 1, el % num_nodes()
                          do k = 0, el % num_coefficients() - 1
                             call add(r, at(columns_at + el % node(i, e)) &
                                  & + this % coefficient_component(column_field, k, column_order) + 1, &
                                  & tau * slopes(cr) * el % basis(0, k, i, q, e), 0.0_dp)
                          end do
                       end do
                    end do
                 end do
              end do
           end do
        end do
      end associate
      do r = start + 1, entries
         constants(r) = constant_of(rows(r))
      end do

    end subroutine projected_condition

    !-----------------------------------------------------------------!
    ! A CONDITION ON A FACE OF THE REGION under a projection fixes, at
    ! every node on the face, the coefficient it determines - the value,
    ! or a derivative along space to the element's order - in place of
    ! that coefficient's row, at the moments the conditions of the
    ! points act at. A condition on a derivative above the element's
    ! order is natural: the projection satisfies it through the terms
    ! of the integration by parts, and it is not stated as a row.
    !-----------------------------------------------------------------!

    subroutine coefficient_condition(j)

      integer, intent(in) :: j

      integer, allocatable :: on_face(:)
      integer :: i, m, k, n, c, coefficient, start_space
      logical :: initial

      start_space = merge(FIRST_COORDINATE + 1, FIRST_COORDINATE, this % manifold % with_time)
      on_face = this % points_on_face(j)
      do i = 1, size(this % condition(j) % equation)
         do m = 1, nm
            acting = implicit_moment(m)
            do k = 1, depth
               if (instant_moment(k) == m) acting = .false.
            end do
            ! the manifold's first instant, the block's own
            initial = first == 1 .and. depth == 0 .and. instant_moment(1) == m
            if (initial) acting = .true.
            if (.not. acting) cycle
            do n = 1, size(on_face)
               call this % space_face_relation(this % condition(j) % equation(i), &
                    & positions(:, (m - 1) * ncells + on_face(n)), column_of, slopes, given, f, along)
               if (along(2) > this % elements % order) then
                  error stop 'operation_residual: under a projection on the cells a condition on a face fixes a &
                       &coefficient of the element, the value or a derivative to its order; a condition on a &
                       &higher derivative is natural and is not stated'
               end if
               coefficient = 0
               if (along(2) > 0) coefficient = (along(1) - start_space) * this % elements % order + along(2)
               ! at the first instant the row of the coefficient's value
               ! the initial data project, elsewhere the row of the law
               if (initial) then
                  row = at((m - 1) * ncells + on_face(n)) + this % coefficient_component(f, coefficient, 0) + 1
               else
                  row = at((m - 1) * ncells + on_face(n)) + this % coefficient_component(f, coefficient, &
                       & primary(f) - this % rule % offset_of_field(f)) + 1
               end if
               if (claimed(row)) cycle
               claimed(row) = .true.
               do c = 1, size(column_of)
                  call add(row, at((m - 1) * ncells + on_face(n)) + column_of(c) + 1, slopes(c), -given)
               end do
            end do
         end do
      end do

    end subroutine coefficient_condition

    !-----------------------------------------------------------------!
    ! The rows the conditions on the faces of the region take at the
    ! manifold's first instant: the value's row of each coefficient
    ! they fix at the nodes on the faces.
    !-----------------------------------------------------------------!

    subroutine taken_at_first_instant()

      integer, allocatable :: on_face(:)
      integer :: jj, ii, n, coefficient, start_space, m1

      start_space = merge(FIRST_COORDINATE + 1, FIRST_COORDINATE, this % manifold % with_time)
      m1 = instant_moment(1)
      do jj = 1, size(this % condition)
         if (.not. this % condition(jj) % manifold % is_space_face()) cycle
         on_face = this % points_on_face(jj)
         do ii = 1, size(this % condition(jj) % equation)
            do n = 1, size(on_face)
               call this % space_face_relation(this % condition(jj) % equation(ii), &
                    & positions(:, (m1 - 1) * ncells + on_face(n)), column_of, slopes, given, f, along)
               coefficient = 0
               if (along(2) > 0) coefficient = (along(1) - start_space) * this % elements % order + along(2)
               taken(at((m1 - 1) * ncells + on_face(n)) + this % coefficient_component(f, coefficient, 0) + 1) = .true.
            end do
         end do
      end do

    end subroutine taken_at_first_instant

  end subroutine face_rows_of

  !===================================================================!
  ! THE DISAGREEMENT OF ANOTHER REPRESENTATION with a solution: for
  ! each law of a stated representation along space, d^k w - w_(k) = 0
  ! in the row of the component w_(k), the other concretion's
  ! derivative of the solution's values less the solution's component
  ! at every point of every instant, times the law's multiplier there,
  ! negated. The solution satisfies its own representation; were the
  ! other one exact, the solution's rows would be perturbed by these
  ! mismatches, and the multiplier of a row is -dJ/d(its datum), so the
  ! contributions, one per point, sum to an estimate of
  ! J(other) - J(solution). On an open line the row of a law is the
  ! trapezoid relation of a segment, and its mismatch is that relation
  ! evaluated on the other concretion's jet of the solution's values,
  ! at the segment's upper point; at the first instant the rows are the
  ! closed ones. Invalid input: a solution without the multipliers of a
  ! representation along space.
  !===================================================================!

  function disagreement(this, solution, other) result(contribution)

    class(discrete_residual)       , intent(in) :: this
    type(discrete_field)           , intent(in) :: solution
    class(derivative_approximation), intent(in) :: other
    real(dp), allocatable :: contribution(:)

    type(stencil) :: rows, lower_rows
    type(discrete_field) :: multipliers
    real(dp), allocatable :: values(:), derivative(:), lower(:)
    real(dp) :: h
    integer :: ncells, ninst, j, i, k, c, f, law(3), point, axis
    logical :: stated, open_line

    ncells = this % points % num_cells()
    ninst  = this % points % num_instants()
    allocate(contribution(this % points % num_points()), source=0.0_dp)
    open_line = this % open_line()
    stated = .false.
    if (allocated(this % representation)) then
       do j = 1, size(this % representation)
          law = represented_component(this % representation(j) % equation(1) % graph(1))
          if (this % manifold % with_time .and. law(3) == FIRST_COORDINATE) cycle
          stated = .true.
          multipliers = solution % fields([this % representation(j) % multiplier % name])
          do i = 1, size(this % representation(j) % equation)
             law  = represented_component(this % representation(j) % equation(i) % graph(1))
             f    = law(1)
             axis = law(3) - merge(1, 0, this % manifold % with_time)
             rows = other % derivative_rows(this % points % cells, axis, law(2))
             if (law(2) > 1) lower_rows = other % derivative_rows(this % points % cells, axis, law(2) - 1)
             allocate(values(ncells), derivative(ncells), lower(ncells))
             do k = 1, ninst
                values = solution % jet(this % rule % offset_of_field(f) + 1, (k - 1) * ncells + 1:k * ncells, 1)
                call rows % row_product([(c, c = 1, ncells)], values, derivative)
                if (open_line .and. k > 1) then
                   ! the other jet's order law(2) - 1, the values at order zero
                   lower = values
                   if (law(2) > 1) call lower_rows % row_product([(c, c = 1, ncells)], values, lower)
                   do c = 1, ncells - 1
                      point = (k - 1) * ncells + c + 1
                      h = this % points % centre(c + 1) - this % points % centre(c)
                      contribution(point) = contribution(point) - multipliers % value(point, i) &
                           & * (lower(c + 1) - lower(c) - h / 2.0_dp * (derivative(c + 1) + derivative(c)))
                   end do
                   cycle
                end if
                do c = 1, ncells
                   point = (k - 1) * ncells + c
                   contribution(point) = contribution(point) - multipliers % value(point, i) &
                        & * (derivative(c) - solution % jet(this % rule % component_at(law(3), law(2), f) + 1, point, 1))
                end do
             end do
             deallocate(values, derivative, lower)
          end do
       end do
    end if
    if (.not. stated) then
       error stop 'operation_residual: the disagreement reads the multipliers of a representation along space, &
            &stated as a residual and paired'
    end if

  end function disagreement

  !===================================================================!
  ! Whether the derivatives along space are relations on the segments
  ! of a line, which leave rows open for the conditions on its faces.
  !===================================================================!

  logical function open_line(this)

    class(discrete_residual), intent(in) :: this

    open_line = .false.
    if (.not. this % with_differences) return
    select type (differences => this % differences)
    type is (segment_trapezoid)
       open_line = .true.
    end select

  end function open_line

  ! whether a block of the sweeps is the whole line at one moment: on
  ! an open line, and for a law projected on the cells of a line
  logical function whole_line(this)
    class(discrete_residual), intent(in) :: this
    whole_line = this % open_line()
    if (allocated(this % elements)) whole_line = whole_line .or. this % points % dimension == 1
  end function whole_line

  !===================================================================!
  ! THE SEGMENT ROWS of a block, at its implicit moments: for each
  ! field, each space coordinate and each order k of its chain, the
  ! trapezoid relation between the components of orders k - 1 and k at
  ! the two ends of every segment, in the row of the component of order
  ! k at the segment's upper point. The points of the line lie in
  ! ascending order. Invalid input: a region of several coordinates,
  ! points that do not ascend, or a mixed component, whose chain along
  ! the line would need conditions of its own at the ends.
  !===================================================================!

  subroutine segment_rows_of(this, nm, at, implicit_moment, rows, columns, weights)

    class(discrete_residual), intent(in)  :: this
    integer                 , intent(in)  :: nm, at(:)
    logical                 , intent(in)  :: implicit_moment(:)
    integer , allocatable   , intent(out) :: rows(:), columns(:)
    real(dp), allocatable   , intent(out) :: weights(:)

    integer :: ncells, m, f, c, k, i, lower, upper, row, here, next, e, entries
    real(dp) :: h

    ncells = this % points % num_cells()
    if (this % points % dimension /= 1) then
       error stop 'operation_residual: the segment relations act on a line; the region has several coordinates'
    end if
    if (any(this % points % centre(2:) <= this % points % centre(:ncells - 1))) then
       error stop 'operation_residual: the segment relations read the points of the line in ascending order'
    end if
    if (any([(this % rule % mixed_degree(f) > 0, f = 1, this % manifold % num_unknowns)])) then
       error stop 'operation_residual: the segment relations do not state a mixed component; its chain &
            &along the line has no conditions at the ends'
    end if
    ! four entries per segment, order and coordinate at each implicit moment
    entries = 0
    do c = merge(2, 1, this % manifold % with_time), this % rule % num_coordinates()
       entries = entries + this % rule % degree_along(c)
    end do
    entries = 4 * (ncells - 1) * entries * this % manifold % num_unknowns * count(implicit_moment)
    allocate(rows(entries), columns(entries), weights(entries))
    e = 0
    do m = 1, nm
       if (.not. implicit_moment(m)) cycle
       do f = 1, this % manifold % num_unknowns
          do c = merge(2, 1, this % manifold % with_time), this % rule % num_coordinates()
             do k = 1, this % rule % degree_along(c)
                lower = this % rule % offset_of_field(f)
                if (k > 1) lower = this % rule % component_at(c, k - 1, f)
                upper = this % rule % component_at(c, k, f)
                do i = 1, ncells - 1
                   h    = this % points % centre(i + 1) - this % points % centre(i)
                   here = at((m - 1) * ncells + i)
                   next = at((m - 1) * ncells + i + 1)
                   row  = next + upper + 1
                   rows(e + 1:e + 4)    = row
                   columns(e + 1:e + 4) = [next + lower + 1, here + lower + 1, next + upper + 1, here + upper + 1]
                   weights(e + 1:e + 4) = [1.0_dp, -1.0_dp, -h / 2.0_dp, -h / 2.0_dp]
                   e = e + 4
                end do
             end do
          end do
       end do
    end do

  end subroutine segment_rows_of

  !===================================================================!
  ! WHERE THE CONDITIONS ON ONE FACE OF THE REGION ARE STATED: the
  ! face's point, and for equation i of condition j the point whose
  ! physics row it displaces, the point and the offset within its
  ! tuple of the row the condition occupies, the field, and the space
  ! derivative it determines (the order zero for a value). Conditions
  ! on the value come first, on the face's own point; the conditions
  ! stated before on the same face take the nearer points. Invalid
  ! input: a face with no point on it, or more conditions on a face
  ! than the region has points.
  !===================================================================!

  subroutine space_face_placement(this, j, primary, face_point, target, row_point, row_offset, field, determined)

    class(discrete_residual), intent(in)  :: this
    integer                 , intent(in)  :: j, primary(:)
    integer                 , intent(out) :: face_point
    integer , allocatable   , intent(out) :: target(:), row_point(:), row_offset(:), field(:), determined(:,:)

    integer , allocatable :: nearest(:), column_of(:)
    real(dp), allocatable :: distance(:), slopes(:)
    logical , allocatable :: taken(:)
    real(dp) :: given
    integer :: ncells, d, n, i, k, used, axis, pass, f, along(2), top

    ncells = this % points % num_cells()
    d      = this % points % dimension
    n      = size(this % condition(j) % equation)
    associate (face => this % condition(j) % manifold)
      ! the points by their distance to the face along its axis, ties
      ! by index
      axis     = face % face_axis - merge(1, 0, this % manifold % with_time)
      distance = abs(this % points % centre(axis:d * ncells:d) - face % face_value)
      allocate(nearest(ncells))
      allocate(taken(ncells), source=.false.)
      do k = 1, ncells
         nearest(k) = minloc(distance, dim=1, mask=.not. taken)
         taken(nearest(k)) = .true.
      end do
      face_point = nearest(1)
      ! a point on the face, up to the rounding of the generated
      ! positions: the vertex-centred mesh of an interval has one at
      ! each end, a cell-centred mesh has none
      if (distance(face_point) > 4.0_dp * epsilon(1.0_dp) * max(abs(face % face_value), distance(nearest(2)))) then
         error stop 'operation_residual: a condition on a face of the region requires a point on the face; &
              &the region of an interval discretized by points has one at each end'
      end if
      used = 0
      do k = 1, j - 1
         if (.not. this % condition(k) % manifold % is_space_face()) cycle
         if (this % condition(k) % manifold % face_axis == face % face_axis .and. &
              & this % condition(k) % manifold % face_value == face % face_value) then
            used = used + size(this % condition(k) % equation)
         end if
      end do
      if (used + n > ncells) then
         error stop 'operation_residual: a face of the region has more conditions than the region has points'
      end if
    end associate
    allocate(target(n), row_point(n), row_offset(n), field(n), determined(2, n))
    if (this % open_line()) then
       ! the segment relations leave the rows of the chain's components
       ! open at the line's lower point: the conditions of both faces
       ! take them in the order they are stated, the rows at one end and
       ! the columns at the other, no physics row displaced
       used = 0
       top  = 0
       do k = 1, size(this % condition)
          if (.not. this % condition(k) % manifold % is_space_face()) cycle
          if (k < j) used = used + size(this % condition(k) % equation)
          top = top + size(this % condition(k) % equation)
       end do
       axis = this % condition(j) % manifold % face_axis
       if (top /= this % rule % degree_along(axis)) then
          error stop 'operation_residual: the faces of an open line state as many conditions as its chain has &
               &orders; each takes one open row'
       end if
       do i = 1, n
          call this % space_face_relation(this % condition(j) % equation(i), &
               & this % points % position(face_point), column_of, slopes, given, f, along)
          target(i)        = 0
          field(i)         = f
          determined(:, i) = [axis, 0]
          row_point(i)     = minloc(this % points % centre(1:ncells), dim=1)
          row_offset(i)    = this % rule % component_at(axis, used + i, f)
       end do
       return
    end if
    do pass = 1, 2
       do i = 1, n
          call this % space_face_relation(this % condition(j) % equation(i), &
               & this % points % position(face_point), column_of, slopes, given, f, along)
          if ((pass == 1) .neqv. (along(2) == 0)) cycle
          used            = used + 1
          target(i)       = nearest(used)
          field(i)        = f
          determined(:, i) = along
          if (along(2) == 0) then
             row_point(i)  = target(i)
             row_offset(i) = primary(f)
          else
             row_point(i)  = face_point
             row_offset(i) = this % rule % component_at(along(1), along(2), f)
          end if
       end do
    end do

  end subroutine space_face_placement

  !===================================================================!
  ! A CONDITION AS AN AFFINE RELATION of the components of the rule's
  ! tuple at one point: the sum over k of slopes(k) times the component
  ! columns(k) equals given, given = -g at zero state. The components
  ! are in the order of the condition's own tuple; leaves(:, k) is
  ! the field, the coordinate and the order of each. Invalid input: a
  ! condition that reads no component of the unknowns.
  !===================================================================!

  subroutine condition_relation(this, equation, values, position, columns, slopes, given, leaves)

    class(discrete_residual), intent(in)  :: this
    type(continuous_field)  , intent(in)  :: equation
    real(dp)                , intent(in)  :: values(:), position(:)
    integer , allocatable   , intent(out) :: columns(:), leaves(:,:)
    real(dp), allocatable   , intent(out) :: slopes(:)
    real(dp)                , intent(out) :: given

    type(expression) :: g
    type(expression_vertex), allocatable :: v(:)
    integer , allocatable :: reads(:)
    real(dp), allocatable :: q(:), grad(:)
    real(dp) :: value
    integer :: k, i, along

    g = equation % graph(1)
    call g % read_components(reads)
    if (size(reads) < 1) then
       error stop 'operation_residual: a condition on a face reads a component of the unknowns'
    end if
    allocate(q(0:g % num_components() - 1), grad(0:g % num_components() - 1), source=0.0_dp)
    call g % gradient_at(q, values, value, grad, position)
    given  = -value
    slopes = grad(reads)
    v = g % vertices()
    allocate(columns(size(reads)), leaves(3, size(reads)))
    do k = 1, size(reads)
       do i = 1, size(v)
          if (v(i) % kind /= VERTEX_LEAF .or. v(i) % position /= ARGUMENT_STATE) cycle
          along = merge(FIRST_COORDINATE, v(i) % along, v(i) % order == 0)
          if (g % component_at(along, v(i) % order, v(i) % field, v(i) % mixed_order) /= reads(k)) cycle
          leaves(:, k) = [v(i) % field, along, v(i) % order]
          columns(k)   = this % rule % component_at(along, v(i) % order, v(i) % field, v(i) % mixed_order)
          exit
       end do
    end do

  end subroutine condition_relation

  !===================================================================!
  ! A condition on a face of the region as that relation, with the
  ! field of its last component, and the coordinate and the order of
  ! the space derivative of highest order it reads (the order zero
  ! when it reads none).
  !===================================================================!

  subroutine space_face_relation(this, equation, position, columns, slopes, given, field, determined)

    class(discrete_residual), intent(in)  :: this
    type(continuous_field)  , intent(in)  :: equation
    real(dp)                , intent(in)  :: position(:)
    integer , allocatable   , intent(out) :: columns(:)
    real(dp), allocatable   , intent(out) :: slopes(:)
    real(dp)                , intent(out) :: given
    integer                 , intent(out) :: field, determined(2)

    integer, allocatable :: leaves(:,:)
    integer :: k

    call this % condition_relation(equation, this % points % expansion % values(), position, columns, slopes, given, &
         & leaves)
    field      = leaves(1, size(columns))
    determined = [FIRST_COORDINATE, 0]
    do k = 1, size(columns)
       if (leaves(2, k) == FIRST_COORDINATE .and. this % manifold % with_time) cycle
       if (leaves(3, k) > determined(2)) determined = leaves(2:3, k)
    end do

  end subroutine space_face_relation

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
  ! A condition on a face of the time factor as that relation, and the
  ! component it determines: of the last field declared among those
  ! read with coefficient +1 or -1 - the adjoint's, declared after the
  ! state, at the far face - the component of highest order, the
  ! relation negated when its coefficient is -1. Invalid input: a
  ! derivative along space, or no component of coefficient +1 or -1.
  !===================================================================!

  subroutine face_relation(this, equation, values, position, columns, slopes, given, row)

    class(discrete_residual), intent(in)  :: this
    type(continuous_field)  , intent(in)  :: equation
    real(dp)                , intent(in)  :: values(:), position(:)
    integer , allocatable   , intent(out) :: columns(:)
    real(dp), allocatable   , intent(out) :: slopes(:)
    real(dp)                , intent(out) :: given
    integer                 , intent(out) :: row

    integer, allocatable :: leaves(:,:)
    real(dp) :: sign
    integer :: k, row_field, row_degree

    call this % condition_relation(equation, values, position, columns, slopes, given, leaves)
    row_field  = 0
    row_degree = -1
    do k = 1, size(columns)
       if (leaves(2, k) /= FIRST_COORDINATE .and. leaves(3, k) > 0) then
          error stop 'operation_residual: a condition on a face reads the values of the unknowns and their &
               &derivatives along time'
       end if
       if (abs(abs(slopes(k)) - 1.0_dp) > 1.0e-12_dp) cycle
       if (leaves(1, k) > row_field .or. (leaves(1, k) == row_field .and. leaves(3, k) > row_degree)) then
          row_field  = leaves(1, k)
          row_degree = leaves(3, k)
          row        = columns(k)
          sign       = slopes(k)
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
