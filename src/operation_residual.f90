!=====================================================================!
! The residual operator: a physics law composed with one or two
! discretization stencils, summed, over points that each store a
! fixed number of components - and the fixed rows of that state,
! read off as an identity rather than solved for.
!
! LEVEL 3 OF THE STRATIFICATION. Two stencils apply to the whole
! state and are summed; the physics applies only at the points
! (`at`), one component read at a time (`primary`), and its result
! is added into the summed stencils there. A row named fixed instead
! reads x(row) - fixed(row): the identity, not the physics.
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
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_residual

  use util_precision   , only : dp
  use operation_action , only : operation, contract, variation
  use operation_action  , only : binding, bound_real_vector
  use view_directed     , only : directed_graph
  use view_directed_stored, only : stored_directed_graph
  use field_calculus    , only : field, FIELD_REAL
  use field_stored      , only : stored_field, typed_field_domain
  use operation_stencil  , only : stencil, combine_triples
  use operation_expression , only : expression, constant, stated

  implicit none

  private
  public :: residual_operator

  type, extends(operation) :: residual_operator

     type(stencil)              , private :: primary_law
     type(stencil), allocatable , private :: connected_law
     type(expression)           , private :: physics
     type(stored_directed_graph), private :: points
     integer , allocatable, private :: at(:)
     integer , allocatable, private :: fixed_rows(:)
     real(dp), allocatable, private :: fixed(:)
     integer                , private :: degrees  = 0
     integer                , private :: connected_degrees = 0
     integer                , private :: unknowns = 0
     integer                , private :: primary  = 0

   contains

     procedure :: name           => residual_name
     procedure :: apply          => residual_apply
     procedure :: max_degree     => residual_max_degree
     procedure :: partial_action => residual_partial_action
     procedure :: explicit_tangent => residual_explicit_tangent

     procedure :: num_unknowns
     procedure :: fixed_unknowns
     procedure :: fixed_mask
     procedure :: fixed_values
     procedure :: num_fixed
     procedure :: first_fixed
     procedure :: stride
     procedure :: num_degrees
     procedure :: primary_degree
     procedure :: num_points
     procedure :: points_at
     procedure :: rule
     procedure :: primary_stencil
     procedure :: has_connected_stencil
     procedure :: connected_stencil
     procedure :: attach_connected_stencil
     procedure :: point_domain
     procedure :: constrain
     procedure :: linearize

  end type residual_operator

  interface residual_operator
     module procedure create
  end interface residual_operator

contains

  !===================================================================!
  ! Build from the two stencils, the physics, and the points the
  ! physics reads: at(p) is the offset of point p's primary
  ! component, always followed by the rest of that point's
  ! components in order. Invalid input: a fixed row without a value
  ! to match, an evaluation point whose components run past the
  ! unknowns, or a fixed row outside the unknowns.
  !===================================================================!

  function create(primary_law, rule, at, unknowns, degrees, primary, fixed_rows, fixed, &
       & connected_law) result(this)

    type(stencil)         , intent(in) :: primary_law
    type(expression)      , intent(in) :: rule
    integer               , intent(in) :: at(:), unknowns, degrees, primary
    integer               , intent(in) :: fixed_rows(:)
    real(dp)              , intent(in) :: fixed(:)
    type(stencil)         , intent(in), optional :: connected_law
    type(residual_operator) :: this

    if (size(fixed_rows) /= size(fixed)) then
       error stop 'operation_residual: one value per fixed component'
    end if
    if (any(fixed_rows < 1) .or. any(fixed_rows > unknowns)) then
       error stop 'operation_residual: every fixed row names an unknown'
    end if
    if (any(at < 0) .or. any(at + degrees > unknowns)) then
       error stop 'operation_residual: an evaluation point''s degree components lie within the unknowns'
    end if

    this % primary_law = primary_law
    if (present(connected_law)) this % connected_law = connected_law
    this % physics    = rule
    this % at      = at
    this % unknowns = unknowns
    this % degrees  = degrees
    ! the rule states how many components a point stores; the degrees
    ! given are the primary law's portion of them
    this % connected_degrees = rule % num_components() - degrees
    this % primary   = primary
    this % fixed_rows = fixed_rows
    this % fixed      = fixed
    this % points = stored_directed_graph(size(at), tails=[integer ::], heads=[integer ::])
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

  pure function fixed_mask(this) result(is_fixed)
    class(residual_operator), intent(in) :: this
    logical, allocatable :: is_fixed(:)
    allocate(is_fixed(this % unknowns), source=.false.)
    is_fixed(this % fixed_rows) = .true.
  end function fixed_mask

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
    primary_degree = this % primary
  end function primary_degree

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

  subroutine point_inputs(this, inputs, x, point_data)
    class(residual_operator), intent(in) :: this
    type(binding)            , intent(in) :: inputs(:)
    real(dp)                 , intent(in) :: x(:)
    type(stored_field), allocatable, intent(out) :: point_data(:)
    type(stored_field) :: state, design
    type(typed_field_domain) :: points, designs
    real(dp), allocatable :: design_values(:)
    call bound_real_vector(inputs, this % argument(2), design_values)
    points = typed_field_domain(this % points % vertex_set(), size(this % at), this % stride())
    designs = typed_field_domain(this % points % vertex_set(), size(design_values))
    state  = points % state(gathered(this, x))
    design = designs % design(design_values)
    point_data = [state, design]
  end subroutine point_inputs

  pure subroutine placed(this, governing, r)
    class(residual_operator), intent(in)    :: this
    real(dp)                 , intent(in)    :: governing(:)
    real(dp)                 , intent(inout) :: r(:)
    integer :: p
    do p = 1, size(this % at)
       r(this % at(p) + this % primary + 1) = &
            & r(this % at(p) + this % primary + 1) + governing(p)
    end do
  end subroutine placed

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

  subroutine state_of(this, inputs, input_graph, x, state)
    class(residual_operator), intent(in)  :: this
    type(binding)             , intent(in)  :: inputs(:)
    class(directed_graph)    , intent(in)  :: input_graph
    real(dp), allocatable, intent(out) :: x(:)
    type(stored_field)   , intent(out) :: state
    type(typed_field_domain) :: states
    call bound_real_vector(inputs, this % argument(1), x)
    if (size(x) /= this % num_unknowns()) then
       error stop 'operation_residual: the state contains one component per degree per unknown point'
    end if
    states = typed_field_domain(input_graph % vertex_set(), size(x))
    state  = states % state(x)
  end subroutine state_of

  subroutine placed_output(this, input_graph, r, output)
    class(residual_operator), intent(in) :: this
    class(directed_graph)    , intent(in) :: input_graph
    real(dp)                 , intent(in) :: r(:)
    class(field), allocatable, intent(inout) :: output
    type(stored_field) :: out
    type(typed_field_domain) :: residuals
    residuals = typed_field_domain(input_graph % vertex_set(), size(r))
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
    real(dp), allocatable :: r(:), governing(:), x(:)
    if (.not. present(inputs)) then
       error stop 'operation_residual: the state and the design are given'
    end if
    call state_of(this, inputs, input_graph, x, state)
    call discretized(this, input_graph, inputs, state, x, r, governing)
    call placed(this, governing, r)
    call accumulate_state(this, x, r)
    call placed_output(this, input_graph, r, output)
  end subroutine residual_apply

  !===================================================================!
  ! THE TWO TERMS OF THE RESIDUAL: the two stencils applied to the
  ! state, summed into r, and the physics at the points, in
  ! governing; or, along a direction v in the state, the partial
  ! action of each in place of its value.
  !===================================================================!

  subroutine discretized(this, input_graph, inputs, state, x, r, governing, v)
    class(residual_operator), intent(in) :: this
    class(directed_graph)    , intent(in) :: input_graph
    type(binding)            , intent(in) :: inputs(:)
    type(stored_field)       , intent(in) :: state
    real(dp)                 , intent(in) :: x(:)
    real(dp), allocatable, intent(out) :: r(:), governing(:)
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
       points    = typed_field_domain(this % points % vertex_set(), size(this % at), this % stride())
       direction = points % direction(gathered(this, v))
       call this % physics % partial_action(this % points, this % physics % bind(point_data), &
            & [variation(this % physics % argument(1), direction)], half)
    else
       call this % physics % apply(this % points, this % physics % bind(point_data), half)
    end if
    call half % real_vector(governing)
  contains
    subroutine stencil_term(op, y)
      type(stencil), intent(in) :: op
      real(dp), allocatable, intent(out) :: y(:)
      type(stored_field) :: along
      type(typed_field_domain) :: domain
      if (present(v)) then
         domain = typed_field_domain(input_graph % vertex_set(), size(v))
         along  = domain % direction(v)
         call op % partial_action(input_graph, op % bind([state]), &
              & [variation(op % argument(1), along)], half)
      else
         call op % apply(input_graph, op % bind([state]), half)
      end if
      call half % real_vector(y)
    end subroutine stencil_term
  end subroutine discretized

  subroutine residual_explicit_tangent(this, input_graph, inputs, which, &
       & rows, columns, weights, available)
    class(residual_operator), intent(in)  :: this
    class(directed_graph)    , intent(in)  :: input_graph
    type(binding)             , intent(in)  :: inputs(:)
    integer              , intent(in)  :: which
    integer , allocatable, intent(out) :: rows(:), columns(:)
    real(dp), allocatable, intent(out) :: weights(:)
    logical              , intent(out) :: available
    type(stored_field) :: state, direction
    type(stored_field), allocatable :: point_data(:)
    type(typed_field_domain) :: points
    class(field), allocatable :: out
    real(dp), allocatable :: x(:), w(:), governing(:), v(:)
    integer , allocatable :: r(:), c(:)
    logical , allocatable :: is_fixed(:)
    integer :: e, d, p, npts, n, kept, count
    available = which == 1
    if (.not. available) return
    n    = this % unknowns
    npts = size(this % at)
    is_fixed = this % fixed_mask()
    call state_of(this, inputs, input_graph, x, state)
    call point_inputs(this, inputs, x, point_data)
    count = this % primary_law % pattern % num_edges() + npts * this % degrees + size(this % fixed_rows)
    if (allocated(this % connected_law)) count = count + this % connected_law % pattern % num_edges()
    allocate(r(count), c(count), w(count))
    kept = 0
    call stencil_triples(this % primary_law, is_fixed, r, c, w, kept)
    if (allocated(this % connected_law)) call stencil_triples(this % connected_law, is_fixed, r, c, w, kept)
    allocate(v(npts * this % degrees))
    do d = 0, this % degrees - 1
       v = 0.0_dp
       do p = 1, npts
          v((p - 1) * this % degrees + d + 1) = 1.0_dp
       end do
       points    = typed_field_domain(this % points % vertex_set(), npts, this % degrees)
       direction = points % direction(v)
       call this % physics % partial_action(this % points, this % physics % bind(point_data), &
            & [variation(this % physics % argument(1), direction)], out)
       call out % real_vector(governing)
       do p = 1, npts
          if (is_fixed(this % at(p) + this % primary + 1)) cycle
          kept    = kept + 1
          r(kept) = this % at(p) + this % primary + 1
          c(kept) = this % at(p) + d + 1
          w(kept) = governing(p)
       end do
    end do
    do e = 1, size(this % fixed_rows)
       kept    = kept + 1
       r(kept) = this % fixed_rows(e)
       c(kept) = this % fixed_rows(e)
       w(kept) = 1.0_dp
    end do
    call combine_triples(n, n, r(1:kept), c(1:kept), w(1:kept), rows, columns, weights)
  end subroutine residual_explicit_tangent

  subroutine stencil_triples(op, is_fixed, r, c, w, kept)
    type(stencil), intent(in)    :: op
    logical      , intent(in)    :: is_fixed(:)
    integer      , intent(inout) :: r(:), c(:)
    real(dp)     , intent(inout) :: w(:)
    integer      , intent(inout) :: kept
    real(dp), allocatable :: weights(:)
    integer :: e, row
    call op % weights % real_vector(weights)
    do e = 1, op % pattern % num_edges()
       row = op % pattern % edge_head(e)
       if (is_fixed(row)) cycle
       kept    = kept + 1
       r(kept) = row
       c(kept) = op % pattern % edge_tail(e)
       w(kept) = weights(e)
    end do
  end subroutine stencil_triples

  subroutine residual_partial_action(this, input_graph, inputs, variations, output)
    class(residual_operator), intent(in)      :: this
    class(directed_graph)    , intent(in)      :: input_graph
    type(binding)             , intent(in)      :: inputs(:)
    type(variation)          , intent(in)      :: variations(:)
    class(field), allocatable, intent(inout) :: output
    type(stored_field) :: state
    real(dp), allocatable :: r(:), governing(:), v(:), x(:)
    call this % require_owned(variations)
    if (size(variations) < 1 .or. size(variations) > this % max_degree()) then
       error stop 'operation_residual: the requested order is within max_degree'
    end if
    call state_of(this, inputs, input_graph, x, state)
    if (size(variations) >= 2) then
       call second_tangent(this, inputs, variations, x, governing)
       allocate(r(this % num_unknowns()), source=0.0_dp)
       call placed(this, governing, r)
       call zero_fixed_rows(this, r)
       call placed_output(this, input_graph, r, output)
       return
    end if
    call variations(1) % direction(v)
    if (variations(1) % argument_is(this % argument(1))) then
       call discretized(this, input_graph, inputs, state, x, r, governing, v)
       call placed(this, governing, r)
       call accumulate_direction(this, v, r)
    else if (variations(1) % argument_is(this % argument(2))) then
       call design_tangent(this, inputs, variations, x, r, governing)
       call placed(this, governing, r)
       call zero_fixed_rows(this, r)
    else
       error stop 'operation_residual: a variation names the state or the design'
    end if
    call placed_output(this, input_graph, r, output)
  end subroutine residual_partial_action

  subroutine second_tangent(this, inputs, variations, x, governing)
    class(residual_operator), intent(in) :: this
    type(binding)             , intent(in) :: inputs(:)
    type(variation)          , intent(in) :: variations(:)
    real(dp)                 , intent(in) :: x(:)
    real(dp), allocatable, intent(out) :: governing(:)
    type(stored_field), allocatable :: point_data(:)
    type(variation), allocatable :: at_points(:)
    class(field), allocatable :: half
    integer :: i
    call point_inputs(this, inputs, x, point_data)
    allocate(at_points(size(variations)))
    do i = 1, size(variations)
       at_points(i) = physics_variation(this, variations(i))
    end do
    call this % physics % partial_action(this % points, this % physics % bind(point_data), at_points, half)
    call half % real_vector(governing)
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
       points    = typed_field_domain(this % points % vertex_set(), size(this % at), this % stride())
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
    real(dp), allocatable, intent(out) :: r(:), governing(:)
    type(stored_field), allocatable :: point_data(:)
    class(field), allocatable :: half
    allocate(r(this % num_unknowns()), source=0.0_dp)
    call point_inputs(this, inputs, x, point_data)
    call this % physics % partial_action(this % points, this % physics % bind(point_data), &
         & [variations(1) % with_argument(this % physics % argument(2))], half)
    call half % real_vector(governing)
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
    integer , allocatable :: sub_of(:), at(:), fixed_rows(:)
    real(dp), allocatable :: fixed(:)
    integer :: e, p, d, inside, npts, ncar

    allocate(sub_of(this % unknowns), source=0)
    do e = 1, size(free)
       sub_of(free(e)) = e
    end do

    npts = 0
    allocate(at(size(this % at)))
    do p = 1, size(this % at)
       inside = 0
       do d = 1, this % degrees
          if (sub_of(this % at(p) + d) > 0) inside = inside + 1
       end do
       if (inside == 0) cycle
       if (inside /= this % degrees) then
          error stop 'operation_residual: a member contains whole points'
       end if
       npts     = npts + 1
       at(npts) = sub_of(this % at(p) + 1) - 1
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
       sub = residual_operator(derived, this % physics, at(1:npts), size(free), &
            & this % degrees, this % primary, fixed_rows(1:ncar), fixed(1:ncar), connected_law=secondary)
    else
       sub = residual_operator(derived, this % physics, at(1:npts), size(free), &
            & this % degrees, this % primary, fixed_rows(1:ncar), fixed(1:ncar))
    end if

  end function constrain

  !===================================================================!
  ! THE EXPLICIT TANGENT, FROZEN INTO A LINEAR RESIDUAL. Ch. 4.6.3 of
  ! the dissertation names this and constrain above the transpose-
  ! Jacobian-vector-product routines. The frozen matrix states the
  ! whole linear map, so the returned residual has no physics of its
  ! own (a zero rule) and no fixed rows. transposed states which of
  ! Jw = rhs or J^Tw = rhs the returned residual's own apply computes;
  ! mark is the tag the returned residual is recorded under (versioned,
  ! this module's own procedure inherited from operation_action).
  !===================================================================!

  function linearize(this, input_graph, inputs, rhs, transposed, mark) result(lin)

    class(residual_operator), intent(in) :: this
    class(directed_graph)   , intent(in) :: input_graph
    type(binding)           , intent(in) :: inputs(:)
    real(dp)                , intent(in) :: rhs(:)
    logical                 , intent(in) :: transposed
    integer                 , intent(in) :: mark
    type(residual_operator) :: lin

    type(stencil) :: a
    integer , allocatable :: r(:), c(:)
    real(dp), allocatable :: w(:)
    logical :: available

    if (size(rhs) /= this % unknowns) then
       error stop 'operation_residual: one right side per unknown'
    end if

    call this % explicit_tangent(input_graph, inputs, 1, r, c, w, available)
    if (.not. available) then
       error stop 'operation_residual: the tangent in the state is explicit'
    end if

    a = stencil(r, c, w, spread(0.0_dp, 1, this % unknowns), 'explicit tangent')
    if (transposed) a = a % transpose()
    call a % constants % set_real_vector(-rhs)

    lin = residual_operator(a, stated(constant(0.0_dp), this % degrees - 1, 'zero'), this % at, &
         & this % unknowns, this % degrees, this % primary, [integer ::], [real(dp) ::])
    call lin % versioned(mark, transposed=transposed)

  end function linearize

end module operation_residual
