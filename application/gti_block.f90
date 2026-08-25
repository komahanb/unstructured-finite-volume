!=====================================================================!
! The residual of one block, as one operation.
!
! A block's rows are of three kinds and they are added here into a
! single statement, because a minimizer drives one statement to zero:
!
!      derived     the scheme's own rows, already assembled as a
!                  stencil; linear in the state, so that stencil is
!                  also their jacobian
!      governing   the physics, at each evaluation point's primary
!                  degree - the one degree no derived row determines
!      carried     the instants a block reaches back over, whose
!                  components are known before it starts; their rows
!                  are the identity less what they hold, so the block
!                  is square and nonsingular
!
! Its two arguments are the state and the design, in that order,
! which is what a minimizer supplies when the design is handed to it
! as a held input.
!
!             WHERE THE PHYSICS IS EVALUATED
!
! At the points given, and nowhere else. A multistep block evaluates
! at its instants, and its points are the instants in order, so the
! components it hands the physics are the state unchanged. A stage
! block evaluates at its stages and recovers its instants from them,
! so its points are the stages and the components are gathered out
! from between them. The physics is nodal either way and never learns
! which it is being asked about.
!
!             THE JACOBIAN
!
! Both halves carry exact partials - the stencil by being linear, the
! physics by differentiating its own rule - so the tangent is exact
! and nothing is differenced. A variation arrives named for this
! statement's argument and is renamed for each half before it is
! passed on, since each half checks the variation against its own,
! and a variation in the state is gathered to the points along with
! the state itself.
!
! A variation in the design is answered too, and it is a different
! statement: the scheme's rows are frozen at the steps they were
! built from and the carried rows hold given numbers, so neither
! varies with the design and only the governing rows do. That partial
! is what a sensitivity reads, by either the tangent or the adjoint.
!
!             WHAT IS REFUSED
!
! A state that is not one component per degree per unknown point; a
! missing argument; a carried row outside the unknowns; an evaluation
! point whose degrees run past them.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_block

  use util_precision  , only : dp
  use operation_action     , only : operation, variation
  use view_directed        , only : directed_graph
  use view_directed_stored , only : stored_directed_graph
  use field_calculus       , only : field
  use field_stored         , only : stored_field
  use graph_fractal        , only : graph
  use operation_stencil    , only : combine_triples, stencil
  use physics_integrand    , only : nodal_integrand

  implicit none

  private
  public :: block_residual

  type, extends(operation) :: block_residual

     type(stencil)                       , private :: derived
     class(nodal_integrand), allocatable  , private :: physics

     ! THE LEVEL BELOW. A stencil over the same unknowns coupling the
     ! components of one instant across the nodes of a spatial mesh:
     ! the spatial operator, linear in the state and independent of
     ! the design. It adds to the derived rows in the apply and in
     ! the tangent, and nowhere else, having no design partial and no
     ! partial above the first. Absent, the block is one node's.
     type(stencil), allocatable, private :: spatial
     type(stored_directed_graph)         , private :: points
     integer , allocatable               , private :: at(:)
     integer , allocatable               , private :: carried(:)
     real(dp), allocatable               , private :: held(:)
     integer                             , private :: degrees  = 0
     integer                             , private :: unknowns = 0
     integer                             , private :: primary  = 0

   contains

     procedure :: name           => block_name
     procedure :: domain         => block_domain
     procedure :: apply          => block_apply
     procedure :: max_degree     => block_max_degree
     procedure :: partial_action => block_partial_action
     procedure :: compiled_tangent => block_compiled_tangent
     procedure :: restricted => block_restricted
     procedure :: num_unknowns
     procedure :: num_degrees
     procedure :: num_points
     procedure :: points_at
     procedure :: num_carried
     procedure :: carried_unknowns
     procedure :: held_values
     procedure :: first_held

  end type block_residual

  interface block_residual
     module procedure create
  end interface block_residual

contains

  function create(derived, physics, at, unknowns, degrees, primary, carried, held, &
       & spatial) result(this)

    type(stencil)         , intent(in) :: derived
    class(nodal_integrand), intent(in) :: physics
    integer               , intent(in) :: at(:), unknowns, degrees, primary
    integer               , intent(in) :: carried(:)
    real(dp)              , intent(in) :: held(:)
    type(stencil)         , intent(in), optional :: spatial
    type(block_residual) :: this

    if (size(carried) /= size(held)) then
       error stop 'gti_block: one value per carried component'
    end if
    if (any(carried < 1) .or. any(carried > unknowns)) then
       error stop 'gti_block: every carried row names an unknown'
    end if
    if (any(at < 0) .or. any(at + degrees > unknowns)) then
       error stop 'gti_block: an evaluation point holds its degrees within the unknowns'
    end if

    this % derived  = derived
    if (present(spatial)) this % spatial = spatial
    allocate(this % physics, source=physics)
    this % at       = at
    this % unknowns = unknowns
    this % degrees  = degrees
    this % primary  = primary
    this % carried  = carried
    this % held     = held

    this % points = stored_directed_graph(size(at), tails=[integer ::], heads=[integer ::])
    call this % declare_arguments(2)

  end function create

  pure integer function num_unknowns(this)

    class(block_residual), intent(in) :: this

    num_unknowns = this % unknowns

  end function num_unknowns

  pure function carried_unknowns(this) result(c)

    class(block_residual), intent(in) :: this
    integer, allocatable :: c(:)

    c = this % carried

  end function carried_unknowns

  pure function held_values(this) result(h)

    class(block_residual), intent(in) :: this
    real(dp), allocatable :: h(:)

    h = this % held

  end function held_values

  pure integer function num_degrees(this)

    class(block_residual), intent(in) :: this

    num_degrees = this % degrees

  end function num_degrees

  pure integer function num_points(this)

    class(block_residual), intent(in) :: this

    num_points = size(this % at)

  end function num_points

  pure function points_at(this) result(at)

    class(block_residual), intent(in) :: this
    integer, allocatable :: at(:)

    at = this % at

  end function points_at

  !===================================================================!
  ! What the first instant a block was given holds. A solver starting
  ! from it begins near the trajectory rather than at nothing, which
  ! for a state of any size is much the same thing as starting at the
  ! wrong end of it.
  !===================================================================!

  pure function first_held(this) result(x)

    class(block_residual), intent(in) :: this
    real(dp), allocatable :: x(:)

    allocate(x(this % degrees), source=0.0_dp)
    if (size(this % held) >= this % degrees) x = this % held(1:this % degrees)

  end function first_held

  pure integer function num_carried(this)

    class(block_residual), intent(in) :: this

    num_carried = size(this % carried)

  end function num_carried

  pure function block_name(this) result(name)

    class(block_residual), intent(in) :: this
    character(len=:), allocatable :: name

    associate (u1 => this); end associate
    name = 'block residual'

  end function block_name

  subroutine block_domain(this, input_graph, domain, num_entries)

    class(block_residual), intent(in)  :: this
    class(directed_graph), intent(in)  :: input_graph
    type(graph)          , intent(out) :: domain
    integer              , intent(out) :: num_entries

    associate (u1 => this); end associate
    domain      = input_graph % vertex_set()
    num_entries = input_graph % num_vertices()

  end subroutine block_domain

  !===================================================================!
  ! The stencil supplies one order and the physics many, so the
  ! statement supplies one.
  !===================================================================!

  pure integer function block_max_degree(this)

    class(block_residual), intent(in) :: this

    associate (u1 => this); end associate
    block_max_degree = 1

  end function block_max_degree

  !===================================================================!
  ! The components held at the evaluation points, taken out from
  ! among the unknowns so that a nodal rule reads them one point at a
  ! time.
  !===================================================================!

  pure function gathered(this, x) result(y)

    class(block_residual), intent(in) :: this
    real(dp)             , intent(in) :: x(:)
    real(dp), allocatable :: y(:)

    integer :: p

    allocate(y(size(this % at) * this % degrees))

    do p = 1, size(this % at)
       y((p - 1) * this % degrees + 1:p * this % degrees) = &
            & x(this % at(p) + 1:this % at(p) + this % degrees)
    end do

  end function gathered

  !===================================================================!
  ! What the physics is handed: the gathered components and the
  ! design, both over the points.
  !===================================================================!

  subroutine point_inputs(this, input_data, x, inputs)

    class(block_residual), intent(in) :: this
    class(field)         , intent(in) :: input_data(:)
    real(dp)             , intent(in) :: x(:)
    type(stored_field), allocatable, intent(out) :: inputs(:)

    type(stored_field) :: state, design
    real(dp), allocatable :: knob(:)

    call input_data(2) % real_vector(knob)

    state = stored_field('state', this % points % vertex_set(), &
         & size(this % at) * this % degrees)
    call state % set_real_vector(gathered(this, x))

    design = stored_field('design', this % points % vertex_set(), size(knob))
    call design % set_real_vector(knob)

    inputs = [state, design]

  end subroutine point_inputs

  !===================================================================!
  ! The governing value at each point, on the row its primary degree
  ! holds.
  !===================================================================!

  pure subroutine placed(this, governing, r)

    class(block_residual), intent(in)    :: this
    real(dp)             , intent(in)    :: governing(:)
    real(dp)             , intent(inout) :: r(:)

    integer :: p

    do p = 1, size(this % at)
       r(this % at(p) + this % primary + 1) = &
            & r(this % at(p) + this % primary + 1) + governing(p)
    end do

  end subroutine placed

  !===================================================================!
  ! The rows of the instants a block reaches back over: what they
  ! hold, less what is proposed for them.
  !===================================================================!

  pure subroutine carry(this, x, r)

    class(block_residual), intent(in)    :: this
    real(dp)             , intent(in)    :: x(:)
    real(dp)             , intent(inout) :: r(:)

    integer :: i

    do i = 1, size(this % carried)
       r(this % carried(i)) = x(this % carried(i)) - this % held(i)
    end do

  end subroutine carry

  pure subroutine carry_direction(this, v, r)

    class(block_residual), intent(in)    :: this
    real(dp)             , intent(in)    :: v(:)
    real(dp)             , intent(inout) :: r(:)

    integer :: i

    do i = 1, size(this % carried)
       r(this % carried(i)) = v(this % carried(i))
    end do

  end subroutine carry_direction

  !===================================================================!
  ! A carried row holds a given number, which varies with nothing.
  !===================================================================!

  pure subroutine carry_held(this, r)

    class(block_residual), intent(in)    :: this
    real(dp)             , intent(inout) :: r(:)

    integer :: i

    do i = 1, size(this % carried)
       r(this % carried(i)) = 0.0_dp
    end do

  end subroutine carry_held

  subroutine state_of(this, input_data, input_graph, x, state)

    class(block_residual), intent(in)  :: this
    class(field)         , intent(in)  :: input_data(:)
    class(directed_graph), intent(in)  :: input_graph
    real(dp), allocatable, intent(out) :: x(:)
    type(stored_field)   , intent(out) :: state

    if (size(input_data) < 2) then
       error stop 'gti_block: the state and the design are given'
    end if

    call input_data(1) % real_vector(x)
    if (size(x) /= this % num_unknowns()) then
       error stop 'gti_block: the state holds one component per degree per unknown point'
    end if

    state = stored_field('state', input_graph % vertex_set(), size(x))
    call state % set_real_vector(x)

  end subroutine state_of

  subroutine placed_output(this, input_graph, r, output)

    class(block_residual), intent(in) :: this
    class(directed_graph), intent(in) :: input_graph
    real(dp)             , intent(in) :: r(:)
    class(field), allocatable, intent(inout) :: output

    type(stored_field) :: out

    out = stored_field(this % name(), input_graph % vertex_set(), size(r))
    call out % set_real_vector(r)

    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine placed_output

  subroutine block_apply(this, input_graph, input_data, output)

    class(block_residual), intent(in)        :: this
    class(directed_graph), intent(in)        :: input_graph
    class(field), intent(in), optional       :: input_data(:)
    class(field), allocatable, intent(inout) :: output

    type(stored_field) :: state
    type(stored_field), allocatable :: inputs(:)
    class(field), allocatable :: half
    real(dp), allocatable :: r(:), governing(:), x(:), coupled(:)

    if (.not. present(input_data)) then
       error stop 'gti_block: the state and the design are given'
    end if

    call state_of(this, input_data, input_graph, x, state)
    call point_inputs(this, input_data, x, inputs)

    call this % derived % apply(input_graph, [state], half)
    call half % real_vector(r)

    if (allocated(this % spatial)) then
       call this % spatial % apply(input_graph, [state], half)
       call half % real_vector(coupled)
       r = r + coupled
    end if

    call this % physics % apply(this % points, inputs, half)
    call half % real_vector(governing)

    call placed(this, governing, r)
    call carry(this, x, r)
    call placed_output(this, input_graph, r, output)

  end subroutine block_apply

  !===================================================================!
  ! The tangent: each half differentiated in its own argument, the
  ! variation renamed for it and, in the state, gathered to the
  ! points the physics reads.
  !===================================================================!

  !===================================================================!
  ! THE BLOCK RESTRICTED to a member of a level - some of its unknowns
  ! - with the rest held at the values given. The derived and the
  ! spatial rows restrict as stencils do, the outside taken into
  ! their constants; the points whose components all lie inside stay
  ! points; the carried rows inside stay carried. A point half inside
  ! stops the program, a member being whole points or nothing.
  !===================================================================!

  function block_restricted(this, kept, values) result(sub)

    class(block_residual), intent(in) :: this
    integer              , intent(in) :: kept(:)
    real(dp)             , intent(in) :: values(:)
    type(block_residual) :: sub

    type(stencil) :: derived, spatial
    integer , allocatable :: sub_of(:), at(:), carried(:)
    real(dp), allocatable :: held(:)
    integer :: e, p, d, inside, npts, ncar

    allocate(sub_of(this % unknowns), source=0)
    do e = 1, size(kept)
       sub_of(kept(e)) = e
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
          error stop 'gti_block: a member holds whole points'
       end if
       npts     = npts + 1
       at(npts) = sub_of(this % at(p) + 1) - 1
    end do

    ncar = 0
    allocate(carried(size(this % carried)), held(size(this % carried)))
    do e = 1, size(this % carried)
       if (sub_of(this % carried(e)) == 0) cycle
       ncar          = ncar + 1
       carried(ncar) = sub_of(this % carried(e))
       held(ncar)    = this % held(e)
    end do

    derived = this % derived % restricted(kept, values)

    if (allocated(this % spatial)) then
       spatial = this % spatial % restricted(kept, values)
       sub = block_residual(derived, this % physics, at(1:npts), size(kept), &
            & this % degrees, this % primary, carried(1:ncar), held(1:ncar), spatial=spatial)
    else
       sub = block_residual(derived, this % physics, at(1:npts), size(kept), &
            & this % degrees, this % primary, carried(1:ncar), held(1:ncar))
    end if

  end function block_restricted

  !===================================================================!
  ! THE COMPILED TANGENT in the state. The block knows its own
  ! structure: the derived rows and the spatial rows are stencils
  ! already, and the physics is nodal, so its tangent at every point
  ! comes from one partial action per degree - a direction of one on
  ! that degree at every point at once, the points being independent.
  ! A carried row is an identity. Triples landing on one entry are
  ! combined. Only the state's tangent is compiled; any other
  ! argument is not available.
  !===================================================================!

  subroutine block_compiled_tangent(this, input_graph, input_data, which, &
       & rows, columns, weights, available)

    class(block_residual), intent(in)  :: this
    class(directed_graph), intent(in)  :: input_graph
    class(field)         , intent(in)  :: input_data(:)
    integer              , intent(in)  :: which
    integer , allocatable, intent(out) :: rows(:), columns(:)
    real(dp), allocatable, intent(out) :: weights(:)
    logical              , intent(out) :: available

    type(stored_field) :: state, direction
    type(stored_field), allocatable :: inputs(:)
    class(field), allocatable :: out
    real(dp), allocatable :: x(:), w(:), governing(:), v(:)
    integer , allocatable :: r(:), c(:)
    logical , allocatable :: is_carried(:)
    integer :: e, d, p, ne, npts, n, kept, count

    available = which == 1
    if (.not. available) return

    n    = this % unknowns
    npts = size(this % at)
    allocate(is_carried(n), source=.false.)
    is_carried(this % carried) = .true.

    call state_of(this, input_data, input_graph, x, state)
    call point_inputs(this, input_data, x, inputs)

    ! room for the derived and spatial triples, the physics's degrees
    ! per point, and the carried identities
    count = this % derived % pattern % num_edges() + npts * this % degrees + size(this % carried)
    if (allocated(this % spatial)) count = count + this % spatial % pattern % num_edges()
    allocate(r(count), c(count), w(count))
    kept = 0

    call stencil_triples(this % derived, is_carried, r, c, w, kept)
    if (allocated(this % spatial)) call stencil_triples(this % spatial, is_carried, r, c, w, kept)

    allocate(v(npts * this % degrees))
    do d = 0, this % degrees - 1
       v = 0.0_dp
       do p = 1, npts
          v((p - 1) * this % degrees + d + 1) = 1.0_dp
       end do
       direction = stored_field('direction', this % points % vertex_set(), size(v))
       call direction % set_real_vector(v)
       call this % physics % partial_action(this % points, inputs, &
            & [variation(this % physics % argument(1), direction)], out)
       call out % real_vector(governing)
       do p = 1, npts
          if (is_carried(this % at(p) + this % primary + 1)) cycle
          kept    = kept + 1
          r(kept) = this % at(p) + this % primary + 1
          c(kept) = this % at(p) + d + 1
          w(kept) = governing(p)
       end do
    end do

    do e = 1, size(this % carried)
       kept    = kept + 1
       r(kept) = this % carried(e)
       c(kept) = this % carried(e)
       w(kept) = 1.0_dp
    end do

    call combine_triples(n, n, r(1:kept), c(1:kept), w(1:kept), rows, columns, weights)

    associate (u1 => ne); end associate

  end subroutine block_compiled_tangent

  !-------------------------------------------------------------------!
  ! A stencil's triples, less those on carried rows, appended.
  !-------------------------------------------------------------------!

  subroutine stencil_triples(op, is_carried, r, c, w, kept)

    type(stencil), intent(in)    :: op
    logical      , intent(in)    :: is_carried(:)
    integer      , intent(inout) :: r(:), c(:)
    real(dp)     , intent(inout) :: w(:)
    integer      , intent(inout) :: kept

    real(dp), allocatable :: weights(:)
    integer :: e, row

    call op % weights % real_vector(weights)
    do e = 1, op % pattern % num_edges()
       row = op % pattern % edge_head(e)
       if (is_carried(row)) cycle
       kept    = kept + 1
       r(kept) = row
       c(kept) = op % pattern % edge_tail(e)
       w(kept) = weights(e)
    end do

  end subroutine stencil_triples

  subroutine block_partial_action(this, input_graph, input_data, variations, output)

    class(block_residual), intent(in)        :: this
    class(directed_graph), intent(in)        :: input_graph
    class(field)         , intent(in)        :: input_data(:)
    type(variation)      , intent(in)        :: variations(:)
    class(field), allocatable, intent(inout) :: output

    type(stored_field) :: state
    real(dp), allocatable :: r(:), governing(:), v(:), x(:)

    call this % require_owned(variations)

    if (size(variations) /= 1) then
       error stop 'gti_block: the requested order is within max_degree'
    end if

    call state_of(this, input_data, input_graph, x, state)
    call variations(1) % direction(v)

    if (variations(1) % argument_is(this % argument(1))) then
       call state_tangent(this, input_graph, input_data, variations, state, x, v, &
            & r, governing)
       call placed(this, governing, r)
       call carry_direction(this, v, r)
    else if (variations(1) % argument_is(this % argument(2))) then
       call design_tangent(this, input_data, variations, x, r, governing)
       call placed(this, governing, r)
       call carry_held(this, r)
    else
       error stop 'gti_block: a variation names the state or the design'
    end if

    call placed_output(this, input_graph, r, output)

  end subroutine block_partial_action

  subroutine state_tangent(this, input_graph, input_data, variations, state, x, v, &
       & r, governing)

    class(block_residual), intent(in) :: this
    class(directed_graph), intent(in) :: input_graph
    class(field)         , intent(in) :: input_data(:)
    type(variation)      , intent(in) :: variations(:)
    type(stored_field)   , intent(in) :: state
    real(dp)             , intent(in) :: x(:), v(:)
    real(dp), allocatable, intent(out) :: r(:), governing(:)

    type(stored_field) :: direction
    type(stored_field), allocatable :: inputs(:)
    class(field), allocatable :: half
    real(dp), allocatable :: coupled(:)

    call this % derived % partial_action(input_graph, [state], &
         & [variations(1) % with_argument(this % derived % argument(1))], half)
    call half % real_vector(r)

    if (allocated(this % spatial)) then
       call this % spatial % partial_action(input_graph, [state], &
            & [variations(1) % with_argument(this % spatial % argument(1))], half)
       call half % real_vector(coupled)
       r = r + coupled
    end if

    call point_inputs(this, input_data, x, inputs)

    direction = stored_field('direction', this % points % vertex_set(), &
         & size(this % at) * this % degrees)
    call direction % set_real_vector(gathered(this, v))

    call this % physics % partial_action(this % points, inputs, &
         & [variation(this % physics % argument(1), direction)], half)
    call half % real_vector(governing)

  end subroutine state_tangent

  !===================================================================!
  ! The design half: the scheme's rows are frozen and the carried
  ! rows hold given numbers, so only the governing rows vary.
  !===================================================================!

  subroutine design_tangent(this, input_data, variations, x, r, governing)

    class(block_residual), intent(in) :: this
    class(field)         , intent(in) :: input_data(:)
    type(variation)      , intent(in) :: variations(:)
    real(dp)             , intent(in) :: x(:)
    real(dp), allocatable, intent(out) :: r(:), governing(:)

    type(stored_field), allocatable :: inputs(:)
    class(field), allocatable :: half

    allocate(r(this % num_unknowns()), source=0.0_dp)
    call point_inputs(this, input_data, x, inputs)

    call this % physics % partial_action(this % points, inputs, &
         & [variations(1) % with_argument(this % physics % argument(2))], half)
    call half % real_vector(governing)

  end subroutine design_tangent

end module gti_block
