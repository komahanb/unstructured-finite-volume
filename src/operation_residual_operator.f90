!=====================================================================!
! Module that states the residual operator: a physics law composed
! with one or two discretization stencils, summed, over points that
! each store a fixed number of components, and the fixed rows of that
! state, read off as an identity rather than solved for.
!
! Stratification. This is level 3. Two stencils apply to the whole
! state and are summed; the physics applies only at the points (at),
! one rule per row it governs (primary), and each rule's result is
! added into the summed stencils on its row. A Lagrangian with k
! multipliers governs k rows, its stationarity in the j-th multiplier
! on primary(j); a plain rule governs one, at every point.
!
!      x = [ tuple of point 1 | tuple of point 2 | ... ]
!            ^ at(1) + 1        ^ at(2) + 1
!
!      r  = S_1(x) + S_2(x)                  the two stencils, affine
!      r(at(p) + primary(j) + 1) += rule_j   the physics at point p
!      r(row) = x(row) - h(row)              a fixed row
!
! Fixed rows. A row named fixed reads x(row) - h(row): the identity,
! not the physics. A fixed value may itself depend on the design (the
! law closes the highest component of the initial tuple), so each one
! states its design rate dh/dnu, and the row's design partial is
! -dh(row)/dnu, not zero. A rate left out is zero, the value being
! data of the problem rather than a function of the design. A
! prescribed row is stated by a third stencil, an affine relation of
! the state, in place of the residual's own terms; a fixed row is the
! case of one unknown.
!
! Derivatives. apply, explicit_jacobian and partial_action are
! composed once here from the two stencils' own apply,
! explicit_jacobian and partial_action and the physics expression's.
! A residual never differentiates anything itself; it only assembles
! what the two composed objects already differentiate.
!
! Concretions. A concretion attaches a domain-specific meaning to the
! points and the stencils (a time march, a spatial mesh, or both
! together) by building them before construction; this type does not
! know which.
!
! Domains. The residual owns its unknown domain U (unknown_graph, one
! vertex per unknown, no edges) and its design domain P (point_domain,
! one vertex per evaluation point). A field of equal length on another
! domain, and a host graph of another identity, are refused: every
! consumer (value, explicit jacobian, tangent and adjoint actions,
! higher partials, the frozen linearization and the constrained
! residual) reads one frozen tuple (Q, nu) on U x P.
!
!      field                             domain   values
!      state Q and its directions        U        one per unknown
!      design nu and its directions      P        designs per point
!      residual, tangents, partials      Y = U    one per unknown
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_residual_operator

  use util_precision                    , only : dp
  use operation_action                  , only : operation, contract, variation
  use operation_action                  , only : binding, bound_real_vector, bound_on
  use view_directed                     , only : directed_graph
  use graph_fractal                     , only : graph
  use view_directed_stored              , only : stored_directed_graph
  use field_calculus                    , only : field, FIELD_REAL
  use field_stored                      , only : stored_field, typed_field_domain
  use operation_stencil                 , only : stencil, combine_triples
  use operation_expression              , only : expression, euler_lagrange, multiplier
  use operation_expression              , only : operator(+), operator(-), operator(*), unknown
  use operation_domain                  , only : continuous_domain, discrete_domain
  use util_derivative_terms             , only : max_subset_width, coefficient
  use util_derivative_terms             , only : value
  use util_derivative_terms             , only : operator(*)
  use operation_derivative_approximation, only : element_quadrature, tested_at
  use operation_exchange                , only : exchange, state_summed_rows

  implicit none

  private
  public :: residual_operator, projection, stationarities, chain_block, lower_bounds
  public :: origin, HISTORY_MOMENTS, OWN_MOMENTS, EVERY_MOMENT

  !-------------------------------------------------------------------!
  ! Type that models the projection of a law on the cells: the law
  ! tested against a basis in place of placed at the points.
  !
  ! An instance is an element at one moment: element_point(i, n) is
  ! the point of node i of instance n, negative where another part
  ! states the node's rows, and element_cell(n) the element e of the
  ! quadrature. A quadrature point of an instance is a site, where the
  ! rules are evaluated.
  !
  !      instance n, quadrature points q = 1..nq
  !      site s = (n - 1) nq + q
  !
  ! The tuple the rules read at a site is formed from the tuples of the
  ! instance's nodes: its component c is the sum over the nodes i and
  ! their coefficients k of
  !
  !      b(formed(0, c), k, i, q, e) x(at(p_i) + formed(k + 1, c) + 1)
  !
  ! where b is the basis, p_i the point of node i, and formed(k + 1, c)
  ! is -1 for a coefficient that forms none. Rule j tests the equation
  ! tested(1, j) with the basis's derivative tested(2, j), and
  ! coefficient_row(k, f) is the row, within a node's tuple, of the
  ! coefficient k of the equation f: the rule's value at a site, times
  ! the weight and the test function there, is added to the row of
  ! every coefficient of every node the instance states.
  !
  !      node tuples --basis--> tuple at site s --rule j--> value
  !      value x test --> row at(p_i) + coefficient_row(k, f) + 1
  !-------------------------------------------------------------------!

  type :: projection
     type(element_quadrature) :: elements
     integer , allocatable    :: element_point(:,:)
     integer , allocatable    :: element_cell(:)
     integer , allocatable    :: tested(:,:)
     integer , allocatable    :: coefficient_row(:,:)
     integer , allocatable    :: formed(:,:)
  end type projection

  !-------------------------------------------------------------------!
  ! Type that models the residual operator: its stencils and rules,
  ! its points and unknowns, its fixed and prescribed rows, and, where
  ! the law is projected on the cells, its sites
  !-------------------------------------------------------------------!

  type, extends(operation) :: residual_operator

     ! The laws: the two stencils, the physics and its rules
     type(stencil)                            , private :: primary_law
     type(stencil)              , allocatable , private :: connected_law
     type(expression)           , allocatable , private :: physics
     type(expression)           , allocatable           :: rules(:)
     ! The points and the unknowns, each a graph without edges
     type(stored_directed_graph)                        :: points
     type(stored_directed_graph)              , private :: unknown_vertices
     ! The offset of each point's tuple, and the fixed rows with their
     ! values and design rates
     integer                    , allocatable           :: at(:)
     integer                    , allocatable , private :: fixed_rows(:)
     real(dp)                   , allocatable , private :: fixed(:)
     real(dp)                   , allocatable , private :: fixed_rate(:)
     ! The rows stated as affine relations of the state in place of the
     ! residual's own: a fixed row is the case of one unknown
     type(stencil)              , allocatable , private :: prescribed
     logical                    , allocatable           :: prescribed_row(:)
     ! The components per point, the primary law's and the connected
     ! law's; the count of unknowns; the row within a tuple each rule
     ! governs
     integer                                            :: degrees           = 0
     integer                                  , private :: connected_degrees = 0
     integer                                  , private :: unknowns          = 0
     integer                    , allocatable           :: primary(:)
     ! The designs at a point: the design is a vector of one value per
     ! design coordinate at every evaluation point
     integer                                            :: designs = 1
     ! The position of every evaluation point, one column per point:
     ! the time, then the space coordinates; no rows when none is given
     real(dp)                   , allocatable           :: point_positions(:,:)
     ! The rows that read every point, the integral of an unknown over
     ! the region: over several images each part states the terms of
     ! its owned columns and summed_by adds the others'
     integer                    , allocatable , private :: summed_rows(:)
     type(exchange)             , allocatable , private :: summed_by
     ! The law projected on the cells, in place of placed at the
     ! points; its sites as the vertices the rules are evaluated over,
     ! and the position of each
     type(projection)           , allocatable           :: projected
     type(stored_directed_graph)                        :: sites
     real(dp)                   , allocatable           :: site_positions(:,:)
     ! The instances of the projection divided over the images, where
     ! every image solves the same rows: each image forms the values
     ! and the matrices of the instances it owns (instance_owned), the
     ! others zero, and their sum over the images is the whole on every
     ! image, each entry formed by one image alone
     logical                                            :: divided = .false.
     ! The tying rows, each of one unknown on its own diagonal reading
     ! others before it: the history, the family's relations, and the
     ! copies and zero coefficients of the unknowns of time alone,
     ! eliminated before a direct solve of the rest
     logical                    , allocatable           :: tying(:)
     ! At a reduced point, a moment at which the law is stated in its
     ! form at the first instant (block_rows), rule j is
     ! rules(variant(j)): the reduced rule's, a law that reads no
     ! component of top degree stated by its derivative along time
     ! (initial_form), appended to the rules, and rule j itself where
     ! the law reads one. Under a projection on the cells the flags are
     ! of the sites, a site reduced where its instance's moment is.
     integer                    , allocatable           :: variant(:)
     logical                    , allocatable , private :: reduced_point(:)
     ! The count of the rules, one per multiplier, their variants after
     integer                                            :: laws = 0

   contains

     ! The operation's symbols and partials
     procedure :: name              => residual_name
     procedure :: domain            => residual_domain
     procedure :: apply             => residual_apply
     procedure :: max_degree        => residual_max_degree
     procedure :: defined_at_zero   => residual_defined_at_zero
     procedure :: partial_action    => residual_partial_action
     procedure :: explicit_jacobian => residual_explicit_jacobian
     procedure, private :: positions => residual_positions

     ! The unknowns, the fixed rows, the rules, the points, the
     ! domains, and the residuals formed from this one
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

     ! The evaluation at the points and at the sites of a projection
     procedure, private :: any_fixed_rate
     procedure, private :: discretized
     procedure, private :: gathered
     procedure, private :: governed
     procedure, private :: governed_on_reals
     procedure, private :: governed_at_points
     procedure, private :: physics_tangent
     procedure, private :: physics_variation
     procedure, private :: placed
     procedure, private :: placed_output
     procedure, private :: point_inputs
     procedure          :: point_of
     procedure, private :: require_direction
     procedure, private :: require_host
     procedure, private :: state_of
     procedure, private :: write_prescribed
     procedure, private :: zero_prescribed_rows
     procedure          :: formed_at_sites
     procedure, private :: design_at_sites
     procedure, private :: tested_into
     procedure, private :: projected_count
     procedure, private :: projected_columns
     procedure, private :: projected_triples
     procedure, private :: instance_owned
     procedure, private :: projected_edges
     procedure, private :: projected_part
     procedure          :: law_at
     procedure, private :: restrict_reduced_points
  end type residual_operator

  interface residual_operator
     module procedure create
  end interface residual_operator

  !-------------------------------------------------------------------!
  ! Type that models the first instant, stage instant and node that
  ! the arrays of a sweep store: one for arrays over the chain, later
  ! for arrays over a window of it
  !-------------------------------------------------------------------!

  type :: lower_bounds
     integer :: instant = 1
     integer :: stage   = 1
     integer :: node    = 1
  end type lower_bounds

  type(lower_bounds), parameter :: origin = lower_bounds(1, 1, 1)

  !-------------------------------------------------------------------!
  ! Type that models one block of the chain: the instants first..last,
  ! the first depth of them the history, solved by the blocks before,
  ! and the rows stated on its moments (each instant and, for a family
  ! marching by stages, the stages of the step into it), with the
  ! unknowns q of those moments, the estimate and then the solution.
  !
  !      first ... first + depth - 1 | first + depth ... last
  !      the history, solved before   | solved here
  !
  ! The unknowns list the moments in the order of time, each moment's
  ! cells in order, each cell's tuple: instant_moment(k) is the moment
  ! of the k-th instant, and the stages of the step into the instant
  ! k + 1 are the moments after it.
  !
  !      moment     m_k | m_k + 1 ... m_k + num_stages | m_(k+1)
  !                 instant k    the stages of the step   instant k + 1
  !      unknowns   moment 1 [cell 1 | cell 2 ...] | moment 2 ...
  !
  !      where m_k = instant_moment(k)
  !
  ! The point manifold is the block of one moment, without a history.
  ! Where the block is linearized, inputs are the rows' inputs at q, by
  ! which their jacobian is formed there, and owner the image of every
  ! point.
  !-------------------------------------------------------------------!

  type :: chain_block

     type(residual_operator)         :: rows
     real(dp), allocatable           :: q(:)
     integer , allocatable           :: instant_moment(:)
     integer                         :: first = 1, last = 1, depth = 0
     integer                         :: num_stages = 0
     integer                         :: num_cells  = 1
     type(stored_field), allocatable :: inputs(:)
     integer           , allocatable :: owner(:)
     ! The design at every point of the block, the rows' stored input
     type(stored_field)              :: design
     ! Where an expansion moves the instants: the rows along time, each
     ! determining component time_row(e) from time_column(e), and
     ! weight_derivative(e, o, d) the o-th derivative of the edge's
     ! weight along the expansion d
     integer , allocatable           :: time_row(:), time_column(:)
     real(dp), allocatable           :: weight_derivative(:,:,:)
     ! Where the instants are pinned: instant_weight(e, a, k, d) the
     ! derivative in the instant instant_first + a - 1 of the k-th
     ! derivative along the expansion d of the edge's weight
     integer                         :: instant_first = 1
     real(dp), allocatable           :: instant_weight(:,:,:,:)

   contains

     procedure :: tuple_at
     procedure :: from_points
     procedure :: to_points

  end type chain_block

  ! The moments of a block at which a column of the manifold is read
  ! or written
  !
  !      HISTORY_MOMENTS   the history instants, no stages
  !      OWN_MOMENTS       the block's own instants, the stages
  !      EVERY_MOMENT      every instant, the stages
  integer, parameter :: HISTORY_MOMENTS = 1, OWN_MOMENTS = 2, EVERY_MOMENT = 3

contains

  !===================================================================!
  ! Returns the laws of a rule: the stationarities of its Lagrangian in
  ! each multiplier, the rule itself where it pairs none
  !
  !      L = lambda_1 c_1 + lambda_2 c_2   -->   [c_1, c_2]
  !      f, no multiplier                  -->   [f]
  !===================================================================!

  function stationarities(rule) result(laws)

    type(expression), intent(in)  :: rule
    type(expression), allocatable :: laws(:)

    integer :: j

    if (rule % num_multipliers() == 0) then
       laws = [rule]
       return
    end if
    allocate(laws(rule % num_multipliers()))
    do j = 1, rule % num_multipliers()
       laws(j) = euler_lagrange(rule, j)
    end do

  end function stationarities

  !===================================================================!
  ! Creates the residual from the two stencils, the physics, and the
  ! points the physics reads. A residual without a rule is the linear
  ! map of its stencils alone, the frozen linearization.
  !
  !      primary_law       the primary stencil, on every row
  !      rule              the physics; absent, the stencils alone
  !      at                the offset of each point's tuple, whose
  !                        components follow in order
  !      unknowns          the count of unknowns
  !      degrees           the primary law's components per point
  !      primary           the row within a tuple the j-th rule governs
  !      fixed_rows, fixed the fixed rows and their values
  !      connected_law     the second stencil
  !      fixed_rate        the design rate of each fixed value
  !      prescribed        rows stated as affine relations of x
  !      positions         the position of every point, a column each
  !      summed_rows       the rows that read every point
  !      projected         the projection of the law on the cells
  !      tying             the tying rows
  !      reduced           the reduced points
  !      reduced_laws      the laws stated there, where varying(j)
  !      laws              the rule's stationarities, formed elsewhere
  !
  ! Invalid input: a fixed row without a value to match, an evaluation
  ! point whose components run past the unknowns, a fixed row outside
  ! the unknowns, or a row count that is not the rule count.
  !===================================================================!

  function create(primary_law, rule, at, unknowns, degrees, primary, fixed_rows, fixed, &
       & connected_law, fixed_rate, prescribed, positions, summed_rows, projected, tying, reduced, &
       & reduced_laws, varying, laws) result(this)

    type(stencil)   , intent(in)           :: primary_law
    type(expression), intent(in), optional :: rule
    integer         , intent(in)           :: at(:), unknowns, degrees, primary(:)
    integer         , intent(in)           :: fixed_rows(:)
    real(dp)        , intent(in)           :: fixed(:)
    type(stencil)   , intent(in), optional :: connected_law
    real(dp)        , intent(in), optional :: fixed_rate(:)
    type(stencil)   , intent(in), optional :: prescribed
    real(dp)        , intent(in), optional :: positions(:,:)
    integer         , intent(in), optional :: summed_rows(:)
    type(projection), intent(in), optional :: projected
    logical         , intent(in), optional :: tying(:)
    logical         , intent(in), optional :: reduced(:)
    type(expression), intent(in), optional :: reduced_laws(:)
    logical         , intent(in), optional :: varying(:)
    type(expression), intent(in), optional :: laws(:)
    type(residual_operator)                :: this

    type(continuous_domain) :: domain
    integer                 :: j, k, variants, components, governed, e, n, q, nq, s, dimension
    character(len=250)      :: message

    if (size(fixed_rows) /= size(fixed)) then
       write(message,'(a,i0,a,i0)') 'operation_lagrangian: one value is required per fixed &
            &component; size(fixed_rows) = ', size(fixed_rows), ', size(fixed) = ', size(fixed)
       error stop trim(message)
    end if
    if (present(fixed_rate)) then
       if (size(fixed_rate) /= size(fixed)) then
          write(message,'(a,i0,a,i0)') 'operation_lagrangian: one design rate is required per &
               &fixed value; size(fixed_rate) = ', size(fixed_rate), ', size(fixed) = ', size(fixed)
          error stop trim(message)
       end if
    end if
    if (any(fixed_rows < 1) .or. any(fixed_rows > unknowns)) then
       write(message,'(a,i0,a,i0,a,i0)') 'operation_lagrangian: every fixed row must name an &
            &unknown 1..', unknowns, '; fixed_rows range from ', minval(fixed_rows), ' to ', &
            & maxval(fixed_rows)
       error stop trim(message)
    end if
    if (any(at < 0) .or. any(at + degrees > unknowns)) then
       write(message,'(a,i0,a,i0,a,i0)') 'operation_lagrangian: an evaluation point''s degree &
            &components must lie within the unknowns; unknowns = ', unknowns, ', minval(at) = ', &
            & minval(at), ', maxval(at) + degrees = ', maxval(at) + degrees
       error stop trim(message)
    end if
    if (degrees < 1) then
       write(message,'(a,i0)') 'operation_lagrangian: the primary degree count must be positive; &
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
          write(message,'(a,i0,a,i0)') 'operation_lagrangian: the primary degree count must lie &
               &within the law''s component count; degrees = ', degrees, ', num_components = ', &
               & components
          error stop trim(message)
       end if
    end if
    if (size(primary) /= governed) then
       write(message,'(a,i0,a,i0)') 'operation_lagrangian: one governed row is required per rule; &
            &size(primary) = ', size(primary), ', rules = ', governed
       error stop trim(message)
    end if
    if (any(primary < 0) .or. any(primary >= degrees)) then
       write(message,'(a,i0,a,i0,a,i0)') 'operation_lagrangian: a governed row must be one of the &
            &point''s degree components 0..', degrees - 1, '; primary ranges from ', &
            & minval(primary), ' to ', maxval(primary)
       error stop trim(message)
    end if

    this % primary_law = primary_law
    if (present(connected_law)) this % connected_law = connected_law
    ! the rules at the reduced points: law j of reduced_laws in place
    ! of its own where varying(j), stored after the laws' own
    variants = 0
    if (present(reduced) .and. present(reduced_laws) .and. present(varying) .and. present(rule)) then
       if (any(reduced)) variants = count(varying)
    end if
    ! a Lagrangian's rules are its stationarities, one per multiplier,
    ! given as laws where they are formed once for many operators; a
    ! residual without a rule governs no row
    if (present(rule)) then
       this % physics = rule
       this % laws    = max(1, rule % num_multipliers())
       allocate(this % rules(this % laws + variants))
       if (present(laws)) then
          if (size(laws) /= this % laws) then
             write(message,'(a,i0,a,i0)') 'operation_lagrangian: one law is given per stationarity of the rule; &
                  &laws = ', size(laws), ', stationarities = ', this % laws
             error stop trim(message)
          end if
          this % rules(1:this % laws) = laws
       else
          this % rules(1:this % laws) = stationarities(rule)
       end if
    else
       allocate(this % rules(0))
       this % laws = 0
    end if
    if (variants > 0) then
       allocate(this % variant(this % laws))
       k = this % laws
       do j = 1, this % laws
          this % variant(j) = j
          if (.not. varying(j)) cycle
          k = k + 1
          this % rules(k) = reduced_laws(j)
          call this % rules(k) % share_arguments(rule)
          this % variant(j) = k
       end do
       this % reduced_point = reduced
    end if
    this % at       = at
    this % unknowns = unknowns
    this % degrees  = degrees
    ! the rule states how many components a point stores; the degrees
    ! given are the primary law's portion of them
    this % connected_degrees = components - degrees
    this % primary    = primary
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
          write(message,'(a,i0,a,i0)') 'operation_lagrangian: one position is required per evaluation &
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
             error stop 'operation_lagrangian: a prescribed row must name an unknown'
          end if
          this % prescribed_row(prescribed % pattern % edge_head(e)) = .true.
       end do
       do j = 1, size(fixed_rows)
          if (this % prescribed_row(fixed_rows(j))) then
             error stop 'operation_lagrangian: a row is both fixed and prescribed'
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
          error stop 'operation_lagrangian: a law projected on the cells is evaluated at the positions of its sites'
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
       ! a site is reduced where the moment of its instance is: every
       ! node of an instance is at that moment
       if (allocated(this % reduced_point)) then
          this % reduced_point = [((reduced(abs(projected % element_point(1, n))), q = 1, nq), &
               & n = 1, size(projected % element_cell))]
       end if
    end if
    call this % declare_arguments(2, [contract(FIELD_REAL, 1), contract(FIELD_REAL, this % designs)])

  end function create

  !===================================================================!
  ! Returns the count of unknowns
  !===================================================================!

  pure integer function num_unknowns(this)

    class(residual_operator), intent(in) :: this

    num_unknowns = this % unknowns

  end function num_unknowns

  !===================================================================!
  ! Returns the fixed rows
  !===================================================================!

  pure function fixed_unknowns(this) result(c)

    class(residual_operator), intent(in) :: this
    integer, allocatable                 :: c(:)

    c = this % fixed_rows

  end function fixed_unknowns

  !===================================================================!
  ! Returns whether each row is one the residual's own terms do not
  ! occupy: fixed or prescribed
  !===================================================================!

  pure function fixed_indicator(this) result(is_fixed)

    class(residual_operator), intent(in) :: this
    logical, allocatable                 :: is_fixed(:)

    allocate(is_fixed(this % unknowns), source=.false.)
    is_fixed(this % fixed_rows) = .true.
    if (allocated(this % prescribed_row)) is_fixed = is_fixed .or. this % prescribed_row

  end function fixed_indicator

  !===================================================================!
  ! Writes the prescribed rows into r: the affine relation of the
  ! state, or its partial along a direction, in place of the residual's
  ! own terms on those rows
  !===================================================================!

  subroutine write_prescribed(this, values, r, along)

    class(residual_operator), intent(in)    :: this
    real(dp)                , intent(in)    :: values(:)
    real(dp)                , intent(inout) :: r(:)
    logical                 , intent(in)    :: along

    type(typed_field_domain)  :: states
    type(stored_field)        :: given(1)
    type(variation)           :: varied(1)
    class(field), allocatable :: half
    real(dp), allocatable     :: y(:)

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

  !===================================================================!
  ! Sets the prescribed rows of r to zero
  !===================================================================!

  pure subroutine zero_prescribed_rows(this, r)

    class(residual_operator), intent(in)    :: this
    real(dp)                , intent(inout) :: r(:)

    if (.not. allocated(this % prescribed_row)) return
    where (this % prescribed_row) r = 0.0_dp

  end subroutine zero_prescribed_rows

  !===================================================================!
  ! Returns the values of the fixed rows
  !===================================================================!

  pure function fixed_values(this) result(h)

    class(residual_operator), intent(in) :: this
    real(dp), allocatable                :: h(:)

    h = this % fixed

  end function fixed_values

  !===================================================================!
  ! Returns the components a point stores: the primary law's and the
  ! connected law's
  !===================================================================!

  pure integer function stride(this)

    class(residual_operator), intent(in) :: this

    stride = this % degrees + this % connected_degrees

  end function stride

  !===================================================================!
  ! Returns the coupling of the unknowns: an edge from every column to
  ! every row of the stencils' entries. The physics placed at the
  ! points reads a point's own tuple alone, so it couples no two
  ! points; a projected law couples every node of an instance with
  ! every other.
  !
  !      entry (row, column) of a stencil   -->   edge column -> row
  !      nodes i, i2 of an instance         -->   edges from the
  !                                               components of i2
  !                                               to the rows of i
  !      a summed row                       -->   no edge into it
  !===================================================================!

  function coupling(this) result(graph_of)

    class(residual_operator), intent(in) :: this
    type(stored_directed_graph)          :: graph_of

    integer , allocatable :: heads(:), tails(:)
    real(dp), allocatable :: w(:)
    logical , allocatable :: is_summed(:)

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

  !===================================================================!
  ! Returns every stencil's terms, as (row, column, weight): the
  ! primary, the connected and the prescribed, in that order
  !===================================================================!

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
  ! Returns the part of the residual an image forms: the rows of the
  ! points the part owns, over the unknowns of the part in its own
  ! numbering, the owned first and the halo after. The physics is
  ! stated at the owned points, each stencil restricted to the owned
  ! rows, and the fixed rows the part owns.
  !
  !      whole     1 2 3 4 5 6 7 8
  !      part      owned [3 4 5 6] | halo [2 7]
  !      numbering        1 2 3 4  |      5 6
  !
  ! Invalid input: the points' tuples not in the order of the unknowns.
  !===================================================================!

  function part_of(this, part) result(sub)

    class(residual_operator), intent(in)    :: this
    ! the part, given the terms of the summed rows it adds
    type(exchange)          , intent(inout) :: part
    type(residual_operator)                 :: sub

    type(stencil), allocatable :: connected, prescribed
    integer      , allocatable :: points(:), at(:), place_of(:), fixed_rows(:), summed(:), heads(:), tails(:)
    integer      , allocatable :: point_place(:)
    real(dp)     , allocatable :: fixed(:), rate(:), weights(:)
    logical      , allocatable :: owned_fixed(:), owned_point(:)
    integer                    :: w, npts, i

    w = this % stride()
    if (any(this % at /= [((i - 1) * w, i = 1, size(this % at))])) then
       error stop 'operation_lagrangian: a residual is cut into parts when its points'' tuples follow &
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
       if (allocated(this % variant)) then
          error stop 'operation_lagrangian: a law projected on the cells and stated by its variant at the &
               &explicit moments is solved on one image; its rows are not cut into parts'
       end if
       allocate(point_place(size(this % at)), source=0)
       point_place(points) = [(i, i = 1, npts)]
       owned_point = [(part % is_owned((i - 1) * w + 1), i = 1, npts)]
       sub = residual_operator(part_stencil(this % primary_law, part, summed), this % physics, at, &
            & size(part % whole_index), this % degrees, this % primary, fixed_rows, fixed, &
            & connected_law=connected, fixed_rate=rate, prescribed=prescribed, &
            & positions=this % point_positions(:, points), &
            & projected=this % projected_part(point_place, owned_point), laws=this % rules(1:this % laws))
    else if (allocated(this % physics)) then
       sub = residual_operator(part_stencil(this % primary_law, part, summed), this % physics, at, &
            & size(part % whole_index), this % degrees, this % primary, fixed_rows, fixed, &
            & connected_law=connected, fixed_rate=rate, prescribed=prescribed, &
            & positions=this % point_positions(:, points), laws=this % rules(1:this % laws))
       call this % restrict_reduced_points(sub, points)
    else
       sub = residual_operator(part_stencil(this % primary_law, part, summed), at=at, &
            & unknowns=size(part % whole_index), degrees=this % degrees, primary=this % primary, &
            & fixed_rows=fixed_rows, fixed=fixed, connected_law=connected, fixed_rate=rate, &
            & prescribed=prescribed, positions=this % point_positions(:, points))
    end if
    sub % summed_by = part

  end function part_of

  !===================================================================!
  ! Returns a stencil's rows the part owns, in the part's numbering,
  ! with their constants
  !===================================================================!

  function part_stencil(op, part, summed) result(sub)

    type(stencil) , intent(in) :: op
    type(exchange), intent(in) :: part
    integer       , intent(in) :: summed(:)
    type(stencil)              :: sub

    integer , allocatable :: heads(:), tails(:), rows(:), columns(:)
    real(dp), allocatable :: w(:), weights(:), constants(:)

    call op % entries(heads, tails, w)
    call op % constants % real_vector(constants)
    call part % restricted_triples(heads, tails, w, rows, columns, weights, summed)
    sub = stencil(rows, columns, weights, merge(constants(part % whole_index), 0.0_dp, part % is_owned), &
         & label=op % name())

  end function part_stencil

  !===================================================================!
  ! Returns the primary law's components per point
  !===================================================================!

  pure integer function num_degrees(this)

    class(residual_operator), intent(in) :: this

    num_degrees = this % degrees

  end function num_degrees

  !===================================================================!
  ! Returns the row within a tuple the first rule governs
  !===================================================================!

  pure integer function primary_degree(this)

    class(residual_operator), intent(in) :: this

    primary_degree = this % primary(1)

  end function primary_degree

  !===================================================================!
  ! Returns the row within a tuple the j-th rule governs
  !===================================================================!

  pure integer function primary_row(this, j)

    class(residual_operator), intent(in) :: this
    integer                 , intent(in) :: j

    primary_row = this % primary(j)

  end function primary_row

  !===================================================================!
  ! Returns the count of the rules, one per multiplier, without their
  ! variants
  !===================================================================!

  pure integer function num_rules(this)

    class(residual_operator), intent(in) :: this

    num_rules = this % laws

  end function num_rules

  !===================================================================!
  ! Returns the rule j at point p: at a reduced point its variant there
  !===================================================================!

  pure integer function law_at(this, j, p)

    class(residual_operator), intent(in) :: this
    integer                 , intent(in) :: j, p

    law_at = j
    if (.not. allocated(this % reduced_point)) return
    if (this % reduced_point(p)) law_at = this % variant(j)

  end function law_at

  !===================================================================!
  ! Sets the rules at the reduced points of a residual formed of the
  ! points given of this one
  !===================================================================!

  subroutine restrict_reduced_points(this, sub, points)

    class(residual_operator), intent(in)    :: this
    type(residual_operator) , intent(inout) :: sub
    integer                 , intent(in)    :: points(:)

    if (.not. allocated(this % reduced_point)) return
    sub % rules         = this % rules
    sub % variant       = this % variant
    sub % reduced_point = this % reduced_point(points)

  end subroutine restrict_reduced_points

  !===================================================================!
  ! Returns the j-th rule
  !===================================================================!

  function rule_of(this, j) result(law)

    class(residual_operator), intent(in) :: this
    integer                 , intent(in) :: j
    type(expression)                     :: law

    law = this % rules(j)

  end function rule_of

  !===================================================================!
  ! Returns the count of evaluation points
  !===================================================================!

  pure integer function num_points(this)

    class(residual_operator), intent(in) :: this

    num_points = size(this % at)

  end function num_points

  !===================================================================!
  ! Returns the offset of each point's tuple
  !===================================================================!

  pure function points_at(this) result(at)

    class(residual_operator), intent(in) :: this
    integer, allocatable                 :: at(:)

    at = this % at

  end function points_at

  !===================================================================!
  ! Returns the first degrees fixed values, or zeros where fewer are
  ! fixed
  !===================================================================!

  pure function first_fixed(this) result(x)

    class(residual_operator), intent(in) :: this
    real(dp), allocatable                :: x(:)

    allocate(x(this % degrees), source=0.0_dp)
    if (size(this % fixed) >= this % degrees) x = this % fixed(1:this % degrees)

  end function first_fixed

  !===================================================================!
  ! Returns the design rate of each fixed value
  !===================================================================!

  pure function fixed_rates(this) result(rate)

    class(residual_operator), intent(in) :: this
    real(dp), allocatable                :: rate(:)

    rate = this % fixed_rate

  end function fixed_rates

  !===================================================================!
  ! Returns the count of fixed rows
  !===================================================================!

  pure integer function num_fixed(this)

    class(residual_operator), intent(in) :: this

    num_fixed = size(this % fixed_rows)

  end function num_fixed

  !===================================================================!
  ! Returns the name of the operation
  !===================================================================!

  pure function residual_name(this) result(name)

    class(residual_operator), intent(in) :: this
    character(len=:), allocatable        :: name

    associate (u1 => this); end associate
    name = 'residual operator'

  end function residual_name

  !===================================================================!
  ! Returns the highest exact degree. The discretization stencils are
  ! linear in the state, so every partial above the first is the
  ! physics expression's alone, or its variants' at the reduced points,
  ! and is exact to every order when there is no physics.
  !===================================================================!

  pure integer function residual_max_degree(this)

    class(residual_operator), intent(in) :: this

    integer :: j

    residual_max_degree = max_subset_width()
    if (allocated(this % physics)) residual_max_degree = this % physics % max_degree()
    if (allocated(this % variant)) then
       residual_max_degree = max(residual_max_degree, &
            & maxval([(this % rules(j) % max_degree(), j = this % laws + 1, size(this % rules))]))
    end if

  end function residual_max_degree

  !===================================================================!
  ! Returns whether the zero state lies in the residual's domain. A law
  ! need not be defined there. The discretization stencils are linear
  ! in the state, so the zero state lies in the residual's domain
  ! exactly when it lies in the domain of every one of its rules: the
  ! radial oscillator's nu / q**3 excludes it; the frozen
  ! linearization, whose rule is the zero expression, does not.
  !===================================================================!

  pure logical function residual_defined_at_zero(this) result(defined)

    class(residual_operator), intent(in) :: this

    integer :: j

    defined = .true.
    do j = 1, size(this % rules)
       defined = defined .and. this % rules(j) % defined_at_zero()
    end do

  end function residual_defined_at_zero

  !===================================================================!
  ! Returns the physics; a residual without a rule, the frozen
  ! linearization, stops the program
  !===================================================================!

  type(expression) function rule(this) result(law)

    class(residual_operator), intent(in) :: this

    if (.not. allocated(this % physics)) then
       error stop 'operation_lagrangian: rule was read from a residual without a rule, the frozen linearization'
    end if
    law = this % physics

  end function rule

  !===================================================================!
  ! Returns the primary stencil
  !===================================================================!

  type(stencil) function primary_stencil(this) result(law)

    class(residual_operator), intent(in) :: this

    law = this % primary_law

  end function primary_stencil

  !===================================================================!
  ! Returns whether the residual has a connected stencil
  !===================================================================!

  pure logical function has_connected_stencil(this)

    class(residual_operator), intent(in) :: this

    has_connected_stencil = allocated(this % connected_law)

  end function has_connected_stencil

  !===================================================================!
  ! Returns the connected stencil
  !===================================================================!

  type(stencil) function connected_stencil(this) result(law)

    class(residual_operator), intent(in) :: this

    law = this % connected_law

  end function connected_stencil

  !===================================================================!
  ! Sets the connected law after construction: a spatial
  ! discretization is laid over an already built primary march, once
  ! the mesh side of a coupled problem is known
  !===================================================================!

  subroutine attach_connected_stencil(this, law)

    class(residual_operator), intent(inout) :: this
    type(stencil)           , intent(in)    :: law

    this % connected_law = law

  end subroutine attach_connected_stencil

  !===================================================================!
  ! Returns the point domain P as a graph, one vertex per evaluation
  ! point, no edges
  !===================================================================!

  type(stored_directed_graph) function point_domain(this) result(points)

    class(residual_operator), intent(in) :: this

    points = this % points

  end function point_domain

  !===================================================================!
  ! Returns the unknown domain U as a graph: the graph the residual is
  ! evaluated over, passed as the host of every apply, tangent and
  ! linearization
  !===================================================================!

  type(stored_directed_graph) function unknown_graph(this) result(unknowns)

    class(residual_operator), intent(in) :: this

    unknowns = this % unknown_vertices

  end function unknown_graph

  !===================================================================!
  ! Returns the vertex set of U: the domain of the state, of every
  ! direction in the state and of the residual itself
  !===================================================================!

  type(graph) function unknown_domain(this) result(domain)

    class(residual_operator), intent(in) :: this

    domain = this % unknown_vertices % vertex_set()

  end function unknown_domain

  !===================================================================!
  ! Returns the design domain P: the vertex set of the evaluation
  ! points, the domain of the design and of every direction in it
  !===================================================================!

  type(graph) function design_domain(this) result(domain)

    class(residual_operator), intent(in) :: this

    domain = this % points % vertex_set()

  end function design_domain

  !===================================================================!
  ! Returns the typed support of the state, of a direction in the
  ! state and of a tangent: fields on U with one value per unknown.
  ! Each support reads its extent from the residual's own graph.
  !===================================================================!

  type(typed_field_domain) function state_fields(this) result(fields)

    class(residual_operator), intent(in) :: this

    fields = typed_field_domain(this % unknown_vertices)

  end function state_fields

  !===================================================================!
  ! Returns the typed support of the residual, a costate and a
  ! forcing: fields on Y = U
  !===================================================================!

  type(typed_field_domain) function residual_fields(this) result(fields)

    class(residual_operator), intent(in) :: this

    fields = typed_field_domain(this % unknown_vertices)

  end function residual_fields

  !===================================================================!
  ! Returns the typed support of the design and of a direction in the
  ! design: fields on P with designs values per point
  !===================================================================!

  type(typed_field_domain) function design_fields(this) result(fields)

    class(residual_operator), intent(in) :: this

    fields = typed_field_domain(this % points, this % designs)

  end function design_fields

  !===================================================================!
  ! Returns the frozen tuple (Q, nu) on U x P: the state x, one value
  ! per unknown, and the design nu, designs values per point. Every
  ! consumer (the value, the explicit jacobian, the tangent and adjoint
  ! actions, the higher partials) reads a tuple built here, so each
  ! linearizes the same function at the same point. Invalid input: a
  ! state or a design of another length.
  !===================================================================!

  function frozen_tuple(this, x, nu) result(inputs)

    class(residual_operator), intent(in) :: this
    real(dp)                , intent(in) :: x(:), nu(:)
    type(stored_field)                   :: inputs(2)

    type(typed_field_domain) :: states, designs
    character(len=250)       :: message

    if (size(x) /= this % unknowns) then
       write(message,'(a,i0,a,i0)') 'operation_lagrangian: the state must contain one component per &
            &degree per unknown point; size(x) = ', size(x), ', unknowns = ', this % unknowns
       error stop trim(message)
    end if
    if (size(nu) /= size(this % at) * this % designs) then
       write(message,'(a,i0,a,i0)') 'operation_lagrangian: the design must contain one value per design &
            &per evaluation point; size(nu) = ', size(nu), ', size(at) * designs = ', &
            & size(this % at) * this % designs
       error stop trim(message)
    end if
    states    = this % state_fields()
    designs   = this % design_fields()
    inputs(1) = states  % state(x)
    inputs(2) = designs % design(nu)

  end function frozen_tuple

  !===================================================================!
  ! Returns the residual's own domain: Y = U, one row per unknown,
  ! whatever graph it is applied on
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
  ! Refuses a host graph other than U: the residual reads and writes
  ! fields on its own unknown domain, and a graph of another identity
  ! with an equal vertex count would type them on a domain they are
  ! not defined on
  !===================================================================!

  subroutine require_host(this, input_graph)

    class(residual_operator), intent(in) :: this
    class(directed_graph)   , intent(in) :: input_graph

    type(graph) :: given, own

    given = input_graph % vertex_set()
    own   = this % unknown_domain()
    if (.not. own % same_as(given)) then
       error stop 'operation_lagrangian: require_host received a graph that is not this residual''s own unknown graph'
    end if

  end subroutine require_host

  !===================================================================!
  ! Refuses a field on another domain than the one stated, or with
  ! another number of values: equal length is not the claim
  !===================================================================!

  subroutine require_field(given, num_given, domain, num_values, message)

    type(graph)     , intent(in) :: given, domain
    integer         , intent(in) :: num_given, num_values
    character(len=*), intent(in) :: message

    character(len=250) :: diagnostic

    if (.not. domain % same_as(given) .or. num_given /= num_values) then
       write(diagnostic,'(a,a,a,l1,a,i0,a,i0)') 'operation_lagrangian: ', message, &
            & '; domain matches = ', domain % same_as(given), ', num_given = ', num_given, &
            & ', num_values = ', num_values
       error stop trim(diagnostic)
    end if

  end subroutine require_field

  !===================================================================!
  ! Returns the state gathered one point at a time, primary components
  ! then connected ones, in the order the physics reads a point's
  ! tuple
  !
  !      x   [ ... | x(at(p) + 1 : at(p) + stride) | ... ]
  !      y   [ point 1's tuple | point 2's tuple | ... ]
  !===================================================================!

  pure function gathered(this, x) result(y)

    class(residual_operator), intent(in) :: this
    real(dp)                , intent(in) :: x(:)
    real(dp), allocatable                :: y(:)

    integer :: p

    allocate(y(size(this % at) * this % stride()))
    do p = 1, size(this % at)
       y((p - 1) * this % stride() + 1:p * this % stride()) = &
            & x(this % at(p) + 1:this % at(p) + this % stride())
    end do

  end function gathered

  !===================================================================!
  ! Returns the physics' inputs at the points: the state gathered per
  ! point and the design read from the frozen tuple; at the sites of a
  ! projection, the tuples and the designs the sites read. The design
  ! is a field on P with designs values per point; another domain or
  ! extent is refused.
  !===================================================================!

  subroutine point_inputs(this, inputs, x, point_data)

    class(residual_operator)       , intent(in)  :: this
    type(binding)                  , intent(in)  :: inputs(:)
    real(dp)                       , intent(in)  :: x(:)
    type(stored_field), allocatable, intent(out) :: point_data(:)

    type(stored_field)       :: state, design
    type(typed_field_domain) :: points, designs
    type(continuous_domain)  :: continuous
    type(discrete_domain)    :: domain
    real(dp), allocatable    :: design_values(:)

    if (.not. bound_on(inputs, this % argument(2), this % design_domain(), size(this % at) * this % designs)) then
       error stop 'operation_lagrangian: the design must be defined on the point domain with one value per &
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

  !===================================================================!
  ! Adds every rule's value at every point into r on the row it
  ! governs; under a projection, tested into the rows of the nodes
  !
  !      r(at(p) + primary(j) + 1) += governing(p, j)
  !===================================================================!

  pure subroutine placed(this, governing, r)

    class(residual_operator), intent(in)    :: this
    real(dp)                , intent(in)    :: governing(:,:)
    real(dp)                , intent(inout) :: r(:)

    integer :: p, j

    if (allocated(this % projected)) then
       call this % tested_into(governing, r)
       return
    end if
    do j = 1, this % laws
       do p = 1, size(this % at)
          r(this % at(p) + this % primary(j) + 1) = &
               & r(this % at(p) + this % primary(j) + 1) + governing(p, j)
       end do
    end do

  end subroutine placed

  !===================================================================!
  ! Returns every rule evaluated at the points, or at the sites of a
  ! projection, one column per rule: the value, or the partial along
  ! the variations given; at a reduced point the rule's variant
  !===================================================================!

  subroutine governed(this, point_data, governing, variations)

    class(residual_operator), intent(in)           :: this
    type(stored_field)      , intent(in)           :: point_data(:)
    real(dp), allocatable   , intent(out)          :: governing(:,:)
    type(variation)         , intent(in), optional :: variations(:)

    class(field), allocatable       :: half
    real(dp), allocatable           :: column(:)
    type(stored_field), allocatable :: bound(:)
    integer                         :: j, coordinates

    ! at the points, or at the sites of a projection
    if (allocated(this % projected)) then
       allocate(governing(this % sites % num_vertices(), this % laws))
    else
       allocate(governing(size(this % at), this % laws))
    end if
    do j = 1, this % laws
       governing(:, j) = evaluated(j)
       ! at the reduced points the rule's variant, evaluated at every
       ! point and taken at those
       if (allocated(this % variant)) then
          if (this % variant(j) /= j) governing(:, j) = merge(evaluated(this % variant(j)), governing(:, j), &
               & this % reduced_point)
       end if
    end do

  contains

    !-----------------------------------------------------------------!
    ! Returns rule k at every point or site: its value, or its partial
    ! along the variations
    !-----------------------------------------------------------------!

    function evaluated(k) result(column)

      integer, intent(in)   :: k
      real(dp), allocatable :: column(:)

      ! a rule that reads a coordinate is given the positions of the
      ! points as its third argument
      bound       = point_data
      coordinates = this % rules(k) % num_position_coordinates()
      if (coordinates > 0) bound = [point_data, this % positions(coordinates)]
      if (allocated(this % projected)) then
         if (present(variations)) then
            call this % rules(k) % partial_action(this % sites, this % rules(k) % bind(bound), variations, half)
         else
            call this % rules(k) % apply(this % sites, this % rules(k) % bind(bound), half)
         end if
      else if (present(variations)) then
         call this % rules(k) % partial_action(this % points, this % rules(k) % bind(bound), &
              & variations, half)
      else
         call this % rules(k) % apply(this % points, this % rules(k) % bind(bound), half)
      end if
      call half % real_vector(column)

    end function evaluated

  end subroutine governed

  !===================================================================!
  ! Returns every rule's value at every point, on the real values of
  ! its tuple, at a reduced point its variant
  !===================================================================!

  subroutine governed_at_points(this, inputs, x, governing)

    class(residual_operator), intent(in)  :: this
    type(binding)           , intent(in)  :: inputs(:)
    real(dp)                , intent(in)  :: x(:)
    real(dp), allocatable   , intent(out) :: governing(:,:)

    real(dp), allocatable :: nu(:), xs(:)
    integer               :: j, p, deg, nd

    call bound_real_vector(inputs, this % argument(2), nu)
    xs  = this % gathered(x)
    deg = this % degrees
    nd  = this % designs
    allocate(governing(size(this % at), this % laws))
    do j = 1, this % laws
       do p = 1, size(this % at)
          governing(p, j) = this % rules(this % law_at(j, p)) % value_on(xs((p - 1) * deg + 1:p * deg), &
               & nu((p - 1) * nd + 1:p * nd), this % point_positions(:, p))
       end do
    end do

  end subroutine governed_at_points

  !===================================================================!
  ! Returns every rule's value at every site of the projection, on the
  ! real values the basis forms there, at a reduced site its variant.
  ! A rule is zero where its test functions vanish and at the sites of
  ! another image's instance; divided over the images, the values are
  ! then summed over them.
  !===================================================================!

  subroutine governed_on_reals(this, inputs, x, governing)

    class(residual_operator), intent(in)  :: this
    type(binding)           , intent(in)  :: inputs(:)
    real(dp)                , intent(in)  :: x(:)
    real(dp), allocatable   , intent(out) :: governing(:,:)

    real(dp), allocatable :: design_values(:), xs(:), ns(:)
    integer               :: j, s, width, nd, nsites, coordinates, nq, k

    call bound_real_vector(inputs, this % argument(2), design_values)
    xs     = this % formed_at_sites(x)
    ns     = this % design_at_sites(design_values)
    width  = size(this % projected % formed, 2)
    nd     = this % designs
    nsites = this % sites % num_vertices()
    allocate(governing(nsites, this % laws))
    nq = this % projected % elements % num_quadrature_points()
    do j = 1, this % laws
       do s = 1, nsites
          ! no rule is evaluated where its test functions vanish, nor at
          ! the sites of another image's instance
          if (.not. this % instance_owned((s - 1) / nq + 1)) then
             governing(s, j) = 0.0_dp
             cycle
          end if
          if (.not. tested_at(this % projected % elements, this % projected % tested(2, j), s - ((s - 1) / nq) * nq, &
               & this % projected % element_cell((s - 1) / nq + 1))) then
             governing(s, j) = 0.0_dp
             cycle
          end if
          k = this % law_at(j, s)
          coordinates = this % rules(k) % num_position_coordinates()
          if (coordinates > 0) then
             governing(s, j) = this % rules(k) % value_on(xs((s - 1) * width + 1:s * width), &
                  & ns((s - 1) * nd + 1:s * nd), this % site_positions(1:coordinates, s))
          else
             governing(s, j) = this % rules(k) % value_on(xs((s - 1) * width + 1:s * width), ns((s - 1) * nd + 1:s * nd))
          end if
       end do
    end do
    if (this % divided) call co_sum(governing)

  end subroutine governed_on_reals

  !===================================================================!
  ! Returns whether this image forms instance n of the projection:
  ! every instance where the instances are not divided; divided, the
  ! instances dealt to the images in turn
  !
  !      instance   1  2  3  4  5  6  7
  !      image      1  2  3  1  2  3  1      on three images
  !===================================================================!

  logical function instance_owned(this, n)

    class(residual_operator), intent(in) :: this
    integer                 , intent(in) :: n

    instance_owned = .not. this % divided
    if (this % divided) instance_owned = mod(n - 1, num_images()) + 1 == this_image()

  end function instance_owned

  !===================================================================!
  ! Returns the first coordinates coordinates of every point's
  ! position, or of every site's under a projection, as a field on the
  ! points, one component per coordinate. A rule that reads a
  ! coordinate the positions do not store stops the program.
  !===================================================================!

  function residual_positions(this, coordinates) result(positions)

    class(residual_operator), intent(in) :: this
    integer                 , intent(in) :: coordinates
    type(stored_field)                   :: positions

    character(len=250) :: message

    if (coordinates > size(this % point_positions, 1)) then
       write(message,'(a,i0,a,i0)') 'operation_lagrangian: a rule reads a coordinate the positions of the &
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

  !===================================================================!
  ! Returns the point whose tuple contains an unknown: the last point
  ! whose offset lies below it, the tuples being placed in order
  !===================================================================!

  pure integer function point_of(this, row) result(p)

    class(residual_operator), intent(in) :: this
    integer                 , intent(in) :: row

    integer :: q

    p = 1
    do q = 1, size(this % at)
       if (this % at(q) < row .and. row <= this % at(q) + this % stride()) p = q
    end do

  end function point_of

  !===================================================================!
  ! Returns whether any fixed value has a nonzero design rate
  !===================================================================!

  pure logical function any_fixed_rate(this) result(yes)

    class(residual_operator), intent(in) :: this

    yes = any(this % fixed_rate /= 0.0_dp)

  end function any_fixed_rate

  !===================================================================!
  ! Returns the state read from the frozen tuple: a field on U with one
  ! value per unknown; another domain or extent is refused
  !===================================================================!

  subroutine state_of(this, inputs, x, state)

    class(residual_operator), intent(in)  :: this
    type(binding)           , intent(in)  :: inputs(:)
    real(dp), allocatable   , intent(out) :: x(:)
    type(stored_field)      , intent(out) :: state

    type(typed_field_domain) :: states
    character(len=250)       :: message

    if (.not. bound_on(inputs, this % argument(1), this % unknown_domain(), this % unknowns)) then
       write(message,'(a,i0)') 'operation_lagrangian: the bound state must be defined on the unknown &
            &domain with one component per degree per unknown point; unknowns = ', this % unknowns
       error stop trim(message)
    end if
    call bound_real_vector(inputs, this % argument(1), x)
    states = this % state_fields()
    state  = states % state(x)

  end subroutine state_of

  !===================================================================!
  ! Refuses a direction off its domain: a direction in the state is a
  ! field on U of the state's extent; a direction in the design is a
  ! field on P with designs values per point
  !===================================================================!

  subroutine require_direction(this, given)

    class(residual_operator), intent(in) :: this
    type(variation)         , intent(in) :: given

    type(graph)           :: along
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
       error stop 'operation_lagrangian: require_direction received a variation naming neither the state nor the design'
    end if

  end subroutine require_direction

  !===================================================================!
  ! Emits r as the residual field on Y = U
  !===================================================================!

  subroutine placed_output(this, r, output)

    class(residual_operator) , intent(in)    :: this
    real(dp)                 , intent(in)    :: r(:)
    class(field), allocatable, intent(inout) :: output

    type(stored_field)       :: out
    type(typed_field_domain) :: residuals

    residuals = this % residual_fields()
    out       = residuals % residual(r, this % name())
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine placed_output

  !===================================================================!
  ! Applies the residual at the frozen tuple: the stencils and the
  ! rules placed on their rows, then the fixed, prescribed and summed
  ! rows
  !
  !      r  = S_1(x) + S_2(x)                  the stencils
  !      r(at(p) + primary(j) + 1) += rule_j   the physics at point p
  !      r(row) = x(row) - h(row)              a fixed row
  !      r(row) = the prescribed relation      a prescribed row
  !      r(row) += the other parts' terms      a summed row
  !===================================================================!

  subroutine residual_apply(this, input_graph, inputs, output)

    class(residual_operator) , intent(in)           :: this
    class(directed_graph)    , intent(in)           :: input_graph
    type(binding)            , intent(in), optional :: inputs(:)
    class(field), allocatable, intent(inout)        :: output

    type(stored_field)    :: state
    real(dp), allocatable :: r(:), governing(:,:), x(:)

    if (.not. present(inputs)) then
       error stop 'operation_lagrangian: residual_apply requires inputs (the state and the design), but the optional inputs argument was not given'
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
  ! Returns the two terms of the residual: the two stencils applied to
  ! the state, summed into r, and the physics at the points, in
  ! governing; or, along a direction v in the state, the partial
  ! action of each in place of its value
  !===================================================================!

  subroutine discretized(this, inputs, state, x, r, governing, v)

    class(residual_operator), intent(in)           :: this
    type(binding)           , intent(in)           :: inputs(:)
    type(stored_field)      , intent(in)           :: state
    real(dp)                , intent(in)           :: x(:)
    real(dp), allocatable   , intent(out)          :: r(:), governing(:,:)
    real(dp)                , intent(in), optional :: v(:)

    type(stored_field)              :: direction
    type(stored_field), allocatable :: point_data(:)
    type(typed_field_domain)        :: points
    type(variation)                 :: varied(1)
    class(field), allocatable       :: half
    real(dp), allocatable           :: coupled(:)

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
    if (allocated(this % tying) .and. .not. present(v)) then
       if (allocated(this % projected)) then
          call this % governed_on_reals(inputs, x, governing)
       else
          call this % governed_at_points(inputs, x, governing)
       end if
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

    !-----------------------------------------------------------------!
    ! Returns a stencil applied to the state, or its partial along v
    !-----------------------------------------------------------------!

    subroutine stencil_term(op, y)

      type(stencil)        , intent(in)  :: op
      real(dp), allocatable, intent(out) :: y(:)

      type(stored_field)       :: along(1), bound(1)
      type(variation)          :: varied(1)
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

  !===================================================================!
  ! Returns the explicit tangent in the state (which = 1) as triples;
  ! in the design none is stated. Duplicate pairs are combined.
  !
  !      row                      entries
  !      a row of the stencils    S_1 + S_2, unless fixed or prescribed
  !      at(p) + primary(j) + 1   rule j's gradient at p, in the column
  !                               at(p) + d + 1 of each component d it
  !                               reads
  !      a fixed row              1 on the diagonal
  !      a prescribed row         the prescribed stencil's entries
  !
  ! Each rule's partials come from one reverse pass per point; a
  ! component no leaf of the rule names has a zero column. Under a
  ! projection, one pass per site, through the basis.
  !===================================================================!

  subroutine residual_explicit_jacobian(this, input_graph, inputs, which, &
       & rows, columns, weights, jacobian_defined)

    class(residual_operator), intent(in)  :: this
    class(directed_graph)   , intent(in)  :: input_graph
    type(binding)           , intent(in)  :: inputs(:)
    integer                 , intent(in)  :: which
    integer , allocatable   , intent(out) :: rows(:), columns(:)
    real(dp), allocatable   , intent(out) :: weights(:)
    logical                 , intent(out) :: jacobian_defined

    type(stored_field)    :: state
    real(dp), allocatable :: x(:), w(:), column(:), xs(:), nu(:), g(:)
    integer , allocatable :: r(:), c(:), reads(:), variant_reads(:)
    logical , allocatable :: is_fixed(:)
    integer               :: e, d, p, npts, n, num_triples, count, j, k, deg
    real(dp)              :: value

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
          error stop 'operation_lagrangian: the design must be defined on the point domain with one value per &
               &design per point'
       end if
       call bound_real_vector(inputs, this % argument(2), nu)
       xs = this % gathered(x)
       allocate(g(0:deg - 1))
    end if
    count = this % primary_law % pattern % num_edges() + npts * this % degrees * this % laws &
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
    do j = 1, this % laws
       if (allocated(this % projected)) exit
       call this % rules(j) % read_components(reads)
       variant_reads = reads
       if (allocated(this % variant)) call this % rules(this % variant(j)) % read_components(variant_reads)
       do p = 1, npts
          if (is_fixed(this % at(p) + this % primary(j) + 1)) cycle
          call this % rules(this % law_at(j, p)) % gradient_at(xs((p - 1) * deg + 1:p * deg), &
               & nu((p - 1) * this % designs + 1:p * this % designs), value, g, this % point_positions(:, p))
          if (this % law_at(j, p) == j) then
             call place(reads)
          else
             call place(variant_reads)
          end if
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

    !-----------------------------------------------------------------!
    ! Appends the rule's partials in the components it reads, in its
    ! row at p
    !-----------------------------------------------------------------!

    subroutine place(components)

      integer, intent(in) :: components(:)

      do k = 1, size(components)
         d = components(k)
         num_triples    = num_triples + 1
         r(num_triples) = this % at(p) + this % primary(j) + 1
         c(num_triples) = this % at(p) + d + 1
         w(num_triples) = g(d)
      end do

    end subroutine place

  end subroutine residual_explicit_jacobian

  !===================================================================!
  ! Appends a stencil's entries as triples, off the fixed rows
  !===================================================================!

  subroutine stencil_triples(op, is_fixed, r, c, w, num_triples)

    type(stencil), intent(in)    :: op
    logical      , intent(in)    :: is_fixed(:)
    integer      , intent(inout) :: r(:), c(:)
    real(dp)     , intent(inout) :: w(:)
    integer      , intent(inout) :: num_triples

    real(dp), allocatable :: weights(:)
    integer               :: e, row

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

  !===================================================================!
  ! Applies the partial action along the variations
  !
  !      first order in the state    the stencils' and the rules'
  !                                  partials; a fixed row v(row), a
  !                                  prescribed row its partial, and
  !                                  the summed rows' terms added
  !      first order in the design   the rules' partials; a fixed row
  !                                  -dh/dnu v(point), a prescribed
  !                                  row zero
  !      second order and above      the rules' partials; the fixed
  !                                  and prescribed rows zero
  !
  ! A repeated partial of a fixed row in the design alone is refused
  ! where any fixed value has a design rate: only the first rate is
  ! stated.
  !===================================================================!

  subroutine residual_partial_action(this, input_graph, inputs, variations, output)

    class(residual_operator) , intent(in)    :: this
    class(directed_graph)    , intent(in)    :: input_graph
    type(binding)            , intent(in)    :: inputs(:)
    type(variation)          , intent(in)    :: variations(:)
    class(field), allocatable, intent(inout) :: output

    type(stored_field)    :: state
    real(dp), allocatable :: r(:), governing(:,:), v(:), x(:)
    integer               :: i
    character(len=250)    :: message

    call this % require_owned(variations)
    if (size(variations) < 1 .or. size(variations) > this % max_degree()) then
       write(message,'(a,i0,a,i0)') 'operation_lagrangian: the requested order must lie within &
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
             error stop 'operation_lagrangian: a fixed value''s design rate is stated to first order &
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

  !===================================================================!
  ! Returns every rule's partial at the points along the variations,
  ! each restated on the rules' own arguments; without physics, no
  ! column
  !===================================================================!

  subroutine physics_tangent(this, inputs, variations, x, governing)

    class(residual_operator), intent(in)  :: this
    type(binding)           , intent(in)  :: inputs(:)
    type(variation)         , intent(in)  :: variations(:)
    real(dp)                , intent(in)  :: x(:)
    real(dp), allocatable   , intent(out) :: governing(:,:)

    type(stored_field), allocatable :: point_data(:)
    type(variation)   , allocatable :: at_points(:)
    integer                         :: i

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

  !===================================================================!
  ! Returns a variation of the residual restated as one of the
  ! physics
  !
  !      along        at the points              at the sites
  !      the state    gathered per point         formed by the basis
  !      the design   as given, on the rules'    that of each site's
  !                   own argument               instance
  !===================================================================!

  function physics_variation(this, given) result(at_points)

    class(residual_operator), intent(in) :: this
    type(variation)         , intent(in) :: given
    type(variation)                      :: at_points

    type(stored_field)       :: direction
    type(typed_field_domain) :: points
    real(dp), allocatable    :: v(:)

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
       error stop 'operation_lagrangian: physics_variation received a variation naming neither the state nor the design'
    end if

  end function physics_variation

  !===================================================================!
  ! Returns the residual constrained to the unknowns free, numbered in
  ! the order free lists them, the rest set to values. Section 4.6.3 of
  ! the dissertation names this and linearize below the transpose-
  ! Jacobian-vector-product routines. The primary and, if present,
  ! the connected stencil each restrict themselves, the way a stencil
  ! already restricts for multigrid; a residual adds only what a
  ! stencil does not state: its own evaluation points and fixed rows,
  ! renumbered onto free.
  !
  !      unknowns   1  2  3  4  5  6
  !      free          2     4  5          sub_of = [0 1 0 2 3 0]
  !
  ! A point split across the boundary (some but not all of its degrees
  ! in free), a law projected on the cells and prescribed rows are
  ! refused.
  !===================================================================!

  function constrain(this, free, values) result(sub)

    class(residual_operator), intent(in) :: this
    integer                 , intent(in) :: free(:)
    real(dp)                , intent(in) :: values(:)
    type(residual_operator)              :: sub

    type(stencil), allocatable :: secondary
    type(stencil)              :: derived
    integer , allocatable      :: sub_of(:), at(:), fixed_rows(:), points(:)
    real(dp), allocatable      :: fixed(:), rate(:)
    integer                    :: e, p, npts, ncar

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
       error stop 'operation_lagrangian: a residual projected on the cells is not restricted; its instances &
            &read the nodes outside the selection'
    end if
    if (allocated(this % prescribed)) then
       error stop 'operation_lagrangian: a residual with prescribed rows is not restricted; the &
            &prescribed relation would lose its unknowns outside the selection'
    end if
    if (allocated(this % connected_law)) then
       secondary = this % connected_law % restricted(free, values)
       sub = residual_operator(derived, this % physics, at, size(free), &
            & this % degrees, this % primary, fixed_rows(1:ncar), fixed(1:ncar), &
            & connected_law=secondary, fixed_rate=rate(1:ncar), positions=this % point_positions(:, points), &
            & laws=this % rules(1:this % laws))
    else
       sub = residual_operator(derived, this % physics, at, size(free), &
            & this % degrees, this % primary, fixed_rows(1:ncar), fixed(1:ncar), &
            & fixed_rate=rate(1:ncar), positions=this % point_positions(:, points), laws=this % rules(1:this % laws))
    end if

    if (allocated(this % summed_rows)) sub % summed_rows = pack(sub_of(this % summed_rows), sub_of(this % summed_rows) > 0)
    call this % restrict_reduced_points(sub, points)

  end function constrain

  !===================================================================!
  ! Returns the points a constraint retains: the evaluation points
  ! whose degree components all lie in free, in the order of this
  ! residual's points. The constrained residual's j-th point is the
  ! selected_points(j)-th point of this one, and a field on P
  ! restricts to the constrained residual's P' by these indices.
  !
  !      degrees = 2, at = [0, 2, 4]
  !      free = [3, 4, 5, 6]   -->   points = [2, 3]
  !      free = [2, 3]         -->   point 1 split: refused
  !
  ! A point split across the boundary (some but not all of its degrees
  ! in free) is an invalid constraint.
  !===================================================================!

  function selected_points(this, free) result(points)

    class(residual_operator), intent(in) :: this
    integer                 , intent(in) :: free(:)
    integer, allocatable                 :: points(:)

    logical, allocatable :: chosen(:)
    integer, allocatable :: selected(:)
    integer              :: e, p, d, inside, npts
    character(len=250)   :: message

    if (any(free < 1) .or. any(free > this % unknowns)) then
       write(message,'(a,i0,a,i0,a,i0)') 'operation_lagrangian: a constraint must select unknowns &
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
          write(message,'(a,i0,a,i0,a,i0)') 'operation_lagrangian: point ', p, ' is split across &
               &the constraint boundary; ', inside, ' of its ', this % degrees, ' degree components are free'
          error stop trim(message)
       end if
       npts           = npts + 1
       selected(npts) = p
    end do
    points = selected(1:npts)

  end function selected_points

  !===================================================================!
  ! Returns the explicit tangent, frozen into a linear residual.
  ! Section 4.6.3 of the dissertation names this and constrain above
  ! the transpose-Jacobian-vector-product routines. The frozen matrix
  ! states the whole linear map, so the returned residual has no rule
  ! and no fixed rows. Its value is J w - rhs, or J^T w - rhs where
  ! transposed, so that its zero solves J w = rhs or J^T w = rhs; the
  ! returned residual is versioned by version_number.
  !
  !      R  --explicit_jacobian-->  J = D_Q R  -->  w |-> J w - rhs
  !===================================================================!

  function linearize(this, input_graph, inputs, rhs, transposed, version_number) result(lin)

    class(residual_operator), intent(in) :: this
    class(directed_graph)   , intent(in) :: input_graph
    type(binding)           , intent(in) :: inputs(:)
    real(dp)                , intent(in) :: rhs(:)
    logical                 , intent(in) :: transposed
    integer                 , intent(in) :: version_number
    type(residual_operator)              :: lin

    type(stencil)         :: a
    integer , allocatable :: r(:), c(:)
    real(dp), allocatable :: w(:)
    logical               :: jacobian_defined
    character(len=250)    :: message

    if (size(rhs) /= this % unknowns) then
       write(message,'(a,i0,a,i0)') 'operation_lagrangian: one right side value is required per &
            &unknown; size(rhs) = ', size(rhs), ', unknowns = ', this % unknowns
       error stop trim(message)
    end if

    call this % explicit_jacobian(input_graph, inputs, 1, r, c, w, jacobian_defined)
    if (.not. jacobian_defined) then
       error stop 'operation_lagrangian: linearize requires the jacobian in the state to be &
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
  ! Returns the tuples the rules read at every site of the projection,
  ! formed from the state by the basis, one tuple after another
  !
  !      y((s - 1) width + c) = sum over the nodes i and coefficients k
  !                             of b(formed(0, c), k, i, q, e)
  !                                x(at(p_i) + formed(k + 1, c) + 1)
  !===================================================================!

  pure function formed_at_sites(this, x) result(y)

    class(residual_operator), intent(in) :: this
    real(dp)                , intent(in) :: x(:)
    real(dp), allocatable                :: y(:)

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

  !===================================================================!
  ! Returns the design at every site of the projection: that of its
  ! instance's first node
  !===================================================================!

  pure function design_at_sites(this, v) result(w)

    class(residual_operator), intent(in) :: this
    real(dp)                , intent(in) :: v(:)
    real(dp), allocatable                :: w(:)

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
  ! Adds each rule's value at every site, times the weight and the
  ! test function of every coefficient of every node the instance
  ! states, to that coefficient's row
  !
  !      r(at(p_i) + coefficient_row(k, f) + 1) +=
  !                     governing(s, j) test(d, k, i, q, e)
  !===================================================================!

  pure subroutine tested_into(this, governing, r)

    class(residual_operator), intent(in)    :: this
    real(dp)                , intent(in)    :: governing(:,:)
    real(dp)                , intent(inout) :: r(:)

    integer  :: j, f, d, n, e, q, s, i, p, k, nq
    real(dp) :: value

    associate (pr => this % projected, el => this % projected % elements)
      nq = el % num_quadrature_points()
      do j = 1, this % laws
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
  ! Returns the count of triples the projection's jacobian may state:
  ! for every instance and every field, nodes x coefficients rows,
  ! each reading nodes x the field's columns
  !===================================================================!

  integer function projected_count(this) result(count)

    class(residual_operator), intent(in) :: this

    integer , allocatable :: offsets(:), slot(:,:), columns_of(:), first(:), column(:)
    integer               :: nodes, nk

    nodes = size(this % projected % element_point, 1)
    nk    = this % projected % elements % num_coefficients()
    call this % projected_columns(offsets, slot, columns_of, first, column)
    count = size(this % projected % element_cell) * nodes * nk * nodes * sum(columns_of)

  end function projected_count

  !===================================================================!
  ! Returns the columns each field's rows read in an instance: of the
  ! components of a node's tuple the basis reads (offsets), those the
  ! rules testing the field read, in the order of offsets.
  !
  !      offsets(o)        o = 1..no, the components the basis reads
  !      columns_of(f)     the count of field f's columns
  !      column(first(f) : first(f) + columns_of(f) - 1)
  !                        their places in offsets
  !      slot(o, f)        the rank of the o-th offset among the
  !                        field's columns, zero where it reads none
  !===================================================================!

  subroutine projected_columns(this, offsets, slot, columns_of, first, column)

    class(residual_operator), intent(in)  :: this
    integer , allocatable   , intent(out) :: offsets(:), slot(:,:), columns_of(:), first(:), column(:)

    integer , allocatable :: reads(:), place(:)
    integer               :: no, nf, j, f, m, k2, o, r

    associate (pr => this % projected)
      offsets = read_offsets(pr % formed)
      no      = size(offsets)
      nf      = size(pr % coefficient_row, 2)
      allocate(place(0:maxval([offsets, 0])), source=0)
      place(offsets) = [(o, o = 1, no)]
      allocate(slot(no, nf), source=0)
      do j = 1, this % laws
         f = pr % tested(1, j)
         call this % rules(j) % read_components(reads)
         call slots_of(reads)
         ! a variant at the reduced sites reads its own components
         if (allocated(this % variant)) then
            call this % rules(this % variant(j)) % read_components(reads)
            call slots_of(reads)
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

    !-----------------------------------------------------------------!
    ! Marks, for the field f, every offset that forms one of the
    ! components given
    !-----------------------------------------------------------------!

    subroutine slots_of(components)

      integer, intent(in) :: components(:)

      associate (pr => this % projected)
        do m = 1, size(components)
           do k2 = 1, ubound(pr % formed, 1)
              if (pr % formed(k2, components(m) + 1) < 0) cycle
              slot(place(pr % formed(k2, components(m) + 1)), f) = 1
           end do
        end do
      end associate

    end subroutine slots_of

  end subroutine projected_columns

  !===================================================================!
  ! Returns the components of a node's tuple the basis reads, each
  ! once, in the order they are first read
  !===================================================================!

  pure function read_offsets(formed) result(offsets)

    integer, intent(in)  :: formed(0:,:)
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

  !===================================================================!
  ! Appends the jacobian of the projection as triples, element by
  ! element. At each site of an instance, the gradient of each rule in
  ! the tuple it reads, times the weight and the test function of a
  ! row and the basis's derivative that forms each component from each
  ! coefficient of each node, is summed into the instance's matrix: its
  ! rows the node, the coefficient and the equation; its columns the
  ! node and the component of a node's tuple the basis reads. The
  ! entries of the rows the instance states are then written once
  ! each, those of one row and column of two instances summed when
  ! combined.
  !
  !      the block of field f in the instance's matrix, from base(f)
  !
  !                           node i2, column o = 1..columns_of(f)
  !                         +----------------------------------+
  !      node i,            |  tau g(comp) b(dd, k2, i2, q, e) |
  !      coefficient k      |  summed over the sites           |
  !                         +----------------------------------+
  !
  !      entry  base(f) + ((i - 1) nk + k) nodes columns_of(f)
  !                     + (i2 - 1) columns_of(f) + o
  !
  ! Divided over the images, each image forms the matrices of the
  ! instances it owns, and the matrices are summed over the images
  ! before their rows are placed.
  !===================================================================!

  subroutine projected_triples(this, x, nu, is_fixed, r, c, w, num_triples)

    class(residual_operator), intent(in)    :: this
    real(dp)                , intent(in)    :: x(:), nu(:)
    logical                 , intent(in)    :: is_fixed(:)
    integer                 , intent(inout) :: r(:), c(:)
    real(dp)                , intent(inout) :: w(:)
    integer                 , intent(inout) :: num_triples

    real(dp), allocatable :: xs(:), ns(:), g(:), local(:), matrices(:,:)
    integer , allocatable :: reads(:), offsets(:), place(:), slot(:,:), columns_of(:), first(:), column(:), base(:)
    integer               :: width, nd, nq, nodes, nk, nf, no, f, n, i, p, k, row, o, i2
    integer               :: at_row

    width   = size(this % projected % formed, 2)
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
      ! divided over the images, the matrices of the instances each
      ! image owns, summed over the images before their rows are placed
      if (this % divided) then
         allocate(matrices(base(nf + 1), size(pr % element_cell)), source=0.0_dp)
         do n = 1, size(pr % element_cell)
            if (this % instance_owned(n)) call instance_matrix(n, matrices(:, n))
         end do
         call co_sum(matrices)
      end if
      allocate(local(base(nf + 1)))
      do n = 1, size(pr % element_cell)
         if (this % divided) then
            local = matrices(:, n)
         else
            call instance_matrix(n, local)
         end if
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

  contains

    !-----------------------------------------------------------------!
    ! Returns the matrix of instance n: the partials of every rule at
    ! its sites, through the basis, times the test functions, in its
    ! block
    !-----------------------------------------------------------------!

    subroutine instance_matrix(n, local)

      integer , intent(in)  :: n
      real(dp), intent(out) :: local(:)

      integer  :: e, q, s, j, f, d, i, k, m, comp, dd, i2, k2, offset, o, at_row
      real(dp) :: value, tau

      associate (pr => this % projected, el => this % projected % elements)
        e     = pr % element_cell(n)
        local = 0.0_dp
        do q = 1, nq
           s = (n - 1) * nq + q
           do j = 1, this % laws
              f = pr % tested(1, j)
              d = pr % tested(2, j)
              if (.not. tested_at(el, d, q, e)) cycle
              call this % rules(this % law_at(j, s)) % read_components(reads)
              call this % rules(this % law_at(j, s)) % gradient_at(xs((s - 1) * width + 1:s * width), &
                   & ns((s - 1) * nd + 1:s * nd), value, g, this % site_positions(:, s))
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
      end associate

    end subroutine instance_matrix

  end subroutine projected_triples

  !===================================================================!
  ! Appends the edges a projected law adds to the coupling: from every
  ! component a node's tuple forms at a site to the row of every
  ! coefficient of every equation of every node the instance states
  !===================================================================!

  subroutine projected_edges(this, heads, tails)

    class(residual_operator), intent(in)    :: this
    integer , allocatable   , intent(inout) :: heads(:), tails(:)

    integer , allocatable :: offsets(:), h(:), t(:)
    integer               :: n, i, p, f, k, i2, p2, o, e

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
  ! Returns the projection of a part: the instances with a node the
  ! part owns, each node's point in the part's numbering, negative for
  ! a node of the halo, whose rows another part states
  !===================================================================!

  function projected_part(this, place_of, owned) result(sub)

    class(residual_operator), intent(in) :: this
    ! the part's point of each point, zero outside the part, and
    ! whether the part owns it
    integer                 , intent(in) :: place_of(:)
    logical                 , intent(in) :: owned(:)
    type(projection)                     :: sub

    logical, allocatable :: stated(:)
    integer              :: n, i, p, m

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
               error stop 'operation_lagrangian: a part reads every node of an instance it states'
            end if
            sub % element_point(i, m) = place_of(p)
            if (.not. owned(place_of(p)) .or. pr % element_point(i, n) < 1) sub % element_point(i, m) = -place_of(p)
         end do
      end do
    end associate

  end function projected_part

  !===================================================================!
  ! Returns the offset, among the block's unknowns, of the tuple of
  ! the moment m at the cell c
  !===================================================================!

  pure integer function tuple_at(this, m, c)

    class(chain_block), intent(in) :: this
    integer           , intent(in) :: m, c

    tuple_at = this % rows % at((m - 1) * this % num_cells + c)

  end function tuple_at

  !===================================================================!
  ! Reads one column of the tuples over the instants, and of the
  ! stages of the step into each, into a vector over the block's
  ! unknowns at the moments named. The stages are left out when their
  ! array is.
  !
  !      instant k of the block, cell c
  !        tuple(:, c, first + k - 1)  -->  x at tuple_at(m_k, c)
  !      stage i of the step into instant k + 1, cell c
  !        stages(:, c, i, first + k)  -->  x at tuple_at(m_k + i, c)
  !
  !      where m_k = instant_moment(k)
  !===================================================================!

  pure subroutine from_points(this, moments, lb, tuple, x, stages)

    class(chain_block), intent(in)           :: this
    integer           , intent(in)           :: moments
    type(lower_bounds), intent(in)           :: lb
    real(dp)          , intent(in)           :: tuple(:,:,lb % instant:)
    real(dp)          , intent(inout)        :: x(:)
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

  !===================================================================!
  ! Writes the column of the tuples, and of the stages, from a vector
  ! over the block's unknowns at the moments named: the inverse of
  ! from_points
  !===================================================================!

  pure subroutine to_points(this, moments, lb, x, tuple, stages)

    class(chain_block), intent(in)              :: this
    integer           , intent(in)              :: moments
    type(lower_bounds), intent(in)              :: lb
    real(dp)          , intent(in)              :: x(:)
    real(dp)          , intent(inout)           :: tuple(:,:,lb % instant:)
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

end module operation_residual_operator
