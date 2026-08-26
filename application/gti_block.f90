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
!                  degree - the one degree no time discretization stencil row determines
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
  use gti_expansion        , only : expansion
  use view_level           , only : level_member, level_num_members, level_is_leaf
  use operation_stencil    , only : combine_triples, stencil
  use operation_family     , only : family
  use operation_weight     , only : scheme_weight
  use operation_coupling   , only : weights_varied
  use operation_expression    , only : expression, constant, stated
  use view_directed        , only : forward

  implicit none

  private
  public :: block_residual, coupling_reach

  !-------------------------------------------------------------------!
  ! THE REACH a time discretization stencil row was built from, in the family's own
  ! numbering of vertices, kept so that the rows can be weighted again
  ! along a direction in the steps: which of the block's steps each
  ! vertex reads, each edge's tail and head vertex and degrees, and
  ! the block unknown each edge determines and reads at the first
  ! node - every other node lies a degrees' width further on.
  !-------------------------------------------------------------------!
  type :: coupling_reach
     integer :: vertices = 0
     integer, allocatable :: step_of(:)
     integer, allocatable :: tails(:), heads(:), source_degree(:), determines(:)
     integer, allocatable :: row(:), column(:)
  end type coupling_reach

  type, extends(operation) :: block_residual

     type(stencil)                       , private :: time_discretization_stencil
     type(expression)     , private :: physics

     ! THE SPATIAL DISCRETIZATION STENCIL. A stencil over the same unknowns coupling the
     ! components of one moment across the nodes of a spatial mesh:
     ! the spatial operator, linear in the state and independent of
     ! the design, laid on the block by spatial_discretization_laid. It adds to the
     ! time discretization stencil rows in the apply and in the tangent, and nowhere else,
     ! having no design partial and no partial above the first.
     ! Absent, the block is one node's.
     type(stencil), allocatable, private :: spatial_discretization_stencil
     type(stored_directed_graph)         , private :: points
     integer , allocatable               , private :: at(:)
     integer , allocatable               , private :: carried(:)
     real(dp), allocatable               , private :: held(:)
     integer                             , private :: degrees  = 0
     integer                             , private :: unknowns = 0
     integer                             , private :: primary  = 0

     ! WHERE THE BLOCK LIES in the graph: its node of the expansion,
     ! and the expansion that node belongs to. Where each unknown lies
     ! - the member of the time level, an instant or a step; its node
     ! of the space level; its moment, the instant or stage whose
     ! values it is among - is read from the graph whenever asked and
     ! held nowhere else: a sweep reads its members and their coupling
     ! from it, the spatial discretization stencil is laid on the moments, the aggregates
     ! a multigrid coarsens by are read off it. A member of a block
     ! keeps the block's node and the unknowns it was restricted to.
     type(graph)    , pointer, private :: node  => null()
     type(expansion), pointer, private :: tower => null()
     integer, allocatable    , private :: kept(:)
     type(coupling_reach), allocatable, private :: reach(:)

   contains

     procedure :: name           => block_name
     procedure :: domain         => block_domain
     procedure :: apply          => block_apply
     procedure :: max_degree     => block_max_degree
     procedure :: partial_action => block_partial_action
     procedure :: compiled_tangent => block_compiled_tangent
     procedure :: restricted => block_restricted
     procedure :: placed_on
     procedure :: slice_of
     procedure :: node_of
     procedure :: moment_of
     procedure, private :: labels_of
     procedure :: num_nodes
     procedure :: spatial_discretization_laid
     procedure :: aggregates
     procedure :: with_reach
     procedure :: rows_varied
     procedure :: linear_block
     procedure :: member_order
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
       & spatial_discretization_stencil) result(this)

    type(stencil)         , intent(in) :: derived
    type(expression)      , intent(in) :: physics
    integer               , intent(in) :: at(:), unknowns, degrees, primary
    integer               , intent(in) :: carried(:)
    real(dp)              , intent(in) :: held(:)
    type(stencil)         , intent(in), optional :: spatial_discretization_stencil
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

    this % time_discretization_stencil  = derived
    if (present(spatial_discretization_stencil)) this % spatial_discretization_stencil = spatial_discretization_stencil
    this % physics = physics
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
    block_max_degree = 2

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

    call this % time_discretization_stencil % apply(input_graph, [state], half)
    call half % real_vector(r)

    if (allocated(this % spatial_discretization_stencil)) then
       call this % spatial_discretization_stencil % apply(input_graph, [state], half)
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
  ! Where the unknowns lie: one time member and one space member per
  ! unknown. Labels of the wrong extent, or one below one, stop the
  ! program.
  !===================================================================!

  subroutine placed_on(this, tower, node)

    class(block_residual), intent(inout)     :: this
    type(expansion)      , intent(in), target :: tower
    type(graph)          , intent(in), target :: node

    this % tower => tower
    this % node  => node

  end subroutine placed_on

  !-------------------------------------------------------------------!
  ! Where every unknown lies, read from the block's node: the slices
  ! are its members; a slice whose first member is a leaf is one
  ! moment, any other slice's members are its moments; a moment holds
  ! one freedom per node of its first component's extent at every
  ! degree, node by node, degrees within a node; the moments lie one
  ! after another. A member of a block reads the block's and keeps
  ! the unknowns it was restricted to. Invalid input: a block placed
  ! nowhere.
  !-------------------------------------------------------------------!

  subroutine labels_of(this, slice, node, moment)

    class(block_residual), intent(in) :: this
    integer, allocatable , intent(out) :: slice(:), node(:), moment(:)

    type(graph), pointer :: one_slice, first
    integer, allocatable :: whole_slice(:), whole_node(:), whole_moment(:)
    integer :: n, k, j, members, moments, g, m, count, u, i, d

    if (.not. associated(this % node)) then
       error stop 'gti_block: the block has not been placed in the graph'
    end if

    n = level_num_members(this % node)
    moments = 0
    do k = 1, n
       moments = moments + members_of(level_member(this % node, k))
    end do
    first => level_member(this % node, 1)
    if (.not. level_is_leaf(level_member(first, 1))) first => level_member(first, 1)
    m     = this % tower % extent_of(level_member(first, 1))
    count = moments * m * this % degrees
    allocate(whole_slice(count), whole_node(count), whole_moment(count))

    g = 0
    do k = 1, n
       one_slice => level_member(this % node, k)
       members = members_of(one_slice)
       do j = 1, members
          g = g + 1
          do i = 1, m
             do d = 0, this % degrees - 1
                u = ((g - 1) * m + (i - 1)) * this % degrees + d + 1
                whole_slice(u)  = k
                whole_node(u)   = i
                whole_moment(u) = g
             end do
          end do
       end do
    end do

    if (allocated(this % kept)) then
       slice  = whole_slice(this % kept)
       node   = whole_node(this % kept)
       moment = whole_moment(this % kept)
    else
       call move_alloc(whole_slice , slice)
       call move_alloc(whole_node  , node)
       call move_alloc(whole_moment, moment)
    end if

  contains

    integer function members_of(one_slice)

      type(graph), intent(in) :: one_slice

      if (level_is_leaf(level_member(one_slice, 1))) then
         members_of = 1
      else
         members_of = level_num_members(one_slice)
      end if

    end function members_of

  end subroutine labels_of

  function slice_of(this) result(slice)

    class(block_residual), intent(in) :: this
    integer, allocatable :: slice(:)

    integer, allocatable :: node(:), moment(:)

    call this % labels_of(slice, node, moment)

  end function slice_of

  function node_of(this) result(node)

    class(block_residual), intent(in) :: this
    integer, allocatable :: node(:)

    integer, allocatable :: slice(:), moment(:)

    call this % labels_of(slice, node, moment)

  end function node_of

  function moment_of(this) result(moment)

    class(block_residual), intent(in) :: this
    integer, allocatable :: moment(:)

    integer, allocatable :: slice(:), node(:)

    call this % labels_of(slice, node, moment)

  end function moment_of

  ! the largest node label: a member of a block keeps the block's
  ! numbering, so this is the extent a map over the nodes must reach
  integer function num_nodes(this)

    class(block_residual), intent(in) :: this

    integer, allocatable :: slice(:), node(:), moment(:)

    num_nodes = 1
    if (.not. associated(this % node)) return
    call this % labels_of(slice, node, moment)
    num_nodes = maxval(node)

  end function num_nodes

  !-------------------------------------------------------------------!
  ! THE SPATIAL DISCRETIZATION STENCIL, laid on this block: a stencil over the nodes is
  ! placed at every moment the block evaluates its physics at, on the
  ! row the physics sits on, and reads the values of that moment. A
  ! moment with no evaluation point - an instant a stage block
  ! recovers - takes no spatial rows, since no physics is stated
  ! there. Invalid input: a stencil over other than the nodes; a
  ! stencil carrying a constant, which would be a source the block
  ! has no place for; a moment holding some nodes and not others.
  !-------------------------------------------------------------------!

  subroutine spatial_discretization_laid(this, spatial_discretization_stencil)

    class(block_residual), intent(inout) :: this
    type(stencil)        , intent(in)    :: spatial_discretization_stencil

    integer , allocatable :: base(:,:), r(:), c(:), slice(:), node(:), moment(:)
    real(dp), allocatable :: lw(:), held(:), w(:)
    integer :: nodes, moments, p, u, e, g, ne, n, rc, cc

    call this % labels_of(slice, node, moment)
    nodes   = maxval(node)
    moments = maxval(moment)
    if (spatial_discretization_stencil % pattern % num_vertices() /= nodes) then
       error stop 'gti_block: the spatial discretization stencil is a stencil over the nodes'
    end if
    call spatial_discretization_stencil % constants % real_vector(held)
    if (any(abs(held) > 0.0_dp)) then
       error stop 'gti_block: the spatial discretization stencil carries no constant'
    end if
    call spatial_discretization_stencil % weights % real_vector(lw)

    ! where each node's components lie at each moment with a point
    allocate(base(nodes, moments), source=-1)
    do p = 1, size(this % at)
       u = this % at(p) + 1
       base(node(u), moment(u)) = this % at(p)
    end do

    ne = spatial_discretization_stencil % pattern % num_edges()
    allocate(r(ne * moments), c(ne * moments), w(ne * moments))
    n = 0
    do g = 1, moments
       do e = 1, ne
          rc = spatial_discretization_stencil % pattern % edge_head(e)
          cc = spatial_discretization_stencil % pattern % edge_tail(e)
          if (base(rc, g) < 0) cycle
          if (base(cc, g) < 0) then
             error stop 'gti_block: a moment holds every node or none'
          end if
          n    = n + 1
          r(n) = base(rc, g) + this % primary + 1
          c(n) = base(cc, g) + 1
          w(n) = lw(e)
       end do
    end do

    this % spatial_discretization_stencil = stencil(r(1:n), c(1:n), w(1:n), spread(0.0_dp, 1, this % unknowns), &
         & 'spatial discretization stencil')

  end subroutine spatial_discretization_laid

  !-------------------------------------------------------------------!
  ! The reach the time discretization stencil rows were built from, given to the block.
  !-------------------------------------------------------------------!

  subroutine with_reach(this, reach)

    class(block_residual), intent(inout) :: this
    type(coupling_reach) , intent(in)    :: reach(:)

    this % reach = reach

  end subroutine with_reach

  !-------------------------------------------------------------------!
  ! The time discretization stencil rows weighted again along a direction in the block's
  ! steps - and a second, given - as the partial of the rows: the
  ! family's weight action carries its partials in the steps, and the
  ! determined component, entering with one, takes no part. Invalid
  ! input: a block built without its reach.
  !-------------------------------------------------------------------!

  function rows_varied(this, scheme, dt, along, along2) result(varied)

    class(block_residual), intent(in) :: this
    class(family)        , intent(in) :: scheme
    real(dp)             , intent(in) :: dt(:), along(:)
    real(dp)             , intent(in), optional :: along2(:)
    type(stencil) :: varied

    integer , allocatable :: r(:), c(:)
    real(dp), allocatable :: w(:), dw(:)
    integer :: k, e, i, nodes, count, n

    if (.not. allocated(this % reach)) then
       error stop 'gti_block: the block was built without its reach'
    end if
    nodes = this % num_nodes()
    count = 0
    do k = 1, size(this % reach)
       count = count + size(this % reach(k) % tails) * nodes
    end do
    allocate(r(count), c(count), w(count))

    n = 0
    do k = 1, size(this % reach)
       associate (reach => this % reach(k))
         if (present(along2)) then
            call weights_varied(scheme_weight(scheme), reach % vertices, reach % tails, &
                 & reach % heads, dt(reach % step_of), along(reach % step_of), &
                 & reach % source_degree, reach % determines, dw, along2(reach % step_of))
         else
            call weights_varied(scheme_weight(scheme), reach % vertices, reach % tails, &
                 & reach % heads, dt(reach % step_of), along(reach % step_of), &
                 & reach % source_degree, reach % determines, dw)
         end if
         do i = 1, nodes
            do e = 1, size(reach % tails)
               n    = n + 1
               r(n) = reach % row(e)    + (i - 1) * this % degrees
               c(n) = reach % column(e) + (i - 1) * this % degrees
               w(n) = -dw(e)
            end do
         end do
       end associate
    end do

    varied = stencil(r, c, w, spread(0.0_dp, 1, this % unknowns), 'varied time discretization stencil')

  end function rows_varied

  !-------------------------------------------------------------------!
  ! The aggregates a multigrid coarsens this block by: the coarse
  ! cell of each unknown's node, at its own moment and degree, and
  ! numbered from one in the order met. Invalid input: a map that
  ! does not reach every node label.
  !-------------------------------------------------------------------!

  function aggregates(this, cell) result(aggregate)

    class(block_residual), intent(in) :: this
    integer              , intent(in) :: cell(:)
    integer, allocatable :: aggregate(:)

    integer, allocatable :: numbered(:), slice(:), node(:), moment(:)
    integer :: u, coarse, key, count

    call this % labels_of(slice, node, moment)
    if (size(cell) < maxval(node)) then
       error stop 'gti_block: a coarse cell for every node'
    end if

    coarse = maxval(cell)
    allocate(aggregate(this % unknowns))
    allocate(numbered(maxval(moment) * coarse * this % degrees), source=0)
    count = 0
    do u = 1, this % unknowns
       key = ((moment(u) - 1) * coarse + cell(node(u)) - 1) * this % degrees &
            & + mod(u - 1, this % degrees) + 1
       if (numbered(key) == 0) then
          count         = count + 1
          numbered(key) = count
       end if
       aggregate(u) = numbered(key)
    end do

  end function aggregates

  !===================================================================!
  ! THE LINEAR BLOCK: the tangent in the state at the inputs given,
  ! frozen, as a block of its own, so that a linear statement A w = b
  ! - or A^T w = b - goes through the same solve as the block it came
  ! from and converges in one newton step. Its derived stencil is A
  ! with -b as its constant, its physics is zero, its points and
  ! degrees are this block's, and it carries no rows: the identities
  ! on the carried components are already in A. A negative stamp is
  ! given to the transpose, which a direct solver reads as the same
  ! factors the other way round.
  !===================================================================!

  function linear_block(this, input_graph, input_data, rhs, transposed, mark) result(lin)

    class(block_residual), intent(in) :: this
    class(directed_graph), intent(in) :: input_graph
    class(field)         , intent(in) :: input_data(:)
    real(dp)             , intent(in) :: rhs(:)
    logical              , intent(in) :: transposed
    integer              , intent(in) :: mark
    type(block_residual) :: lin

    type(stencil) :: a
    integer , allocatable :: r(:), c(:)
    real(dp), allocatable :: w(:)
    logical :: available

    if (size(rhs) /= this % unknowns) then
       error stop 'gti_block: one right side per unknown'
    end if

    call this % compiled_tangent(input_graph, input_data, 1, r, c, w, available)
    if (.not. available) then
       error stop 'gti_block: the tangent in the state compiles'
    end if

    a = stencil(r, c, w, spread(0.0_dp, 1, this % unknowns), 'frozen tangent')
    if (transposed) a = a % transpose()
    call a % constants % set_real_vector(-rhs)

    lin = block_residual(a, stated(constant(0.0_dp), this % degrees - 1, 'zero'), this % at, this % unknowns, &
         & this % degrees, this % primary, [integer ::], [real(dp) ::])
    call lin % stamped(mark, transposed=a % pattern % transposed())
    lin % tower => this % tower
    lin % node  => this % node
    if (allocated(this % kept)) lin % kept = this % kept

  end function linear_block

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

    type(stencil) :: derived, spatial_discretization_stencil
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

    derived = this % time_discretization_stencil % restricted(kept, values)

    if (allocated(this % spatial_discretization_stencil)) then
       spatial_discretization_stencil = this % spatial_discretization_stencil % restricted(kept, values)
       sub = block_residual(derived, this % physics, at(1:npts), size(kept), &
            & this % degrees, this % primary, carried(1:ncar), held(1:ncar), spatial_discretization_stencil=spatial_discretization_stencil)
    else
       sub = block_residual(derived, this % physics, at(1:npts), size(kept), &
            & this % degrees, this % primary, carried(1:ncar), held(1:ncar))
    end if
    ! the member lies where the block lies, on the unknowns kept
    sub % tower => this % tower
    sub % node  => this % node
    if (allocated(this % kept)) then
       sub % kept = this % kept(kept)
    else
       sub % kept = kept
    end if

  end function block_restricted

  !===================================================================!
  ! THE COMPILED TANGENT in the state. The block knows its own
  ! structure: the time discretization stencil rows and the spatial rows are stencils
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
    count = this % time_discretization_stencil % pattern % num_edges() + npts * this % degrees + size(this % carried)
    if (allocated(this % spatial_discretization_stencil)) count = count + this % spatial_discretization_stencil % pattern % num_edges()
    allocate(r(count), c(count), w(count))
    kept = 0

    call stencil_triples(this % time_discretization_stencil, is_carried, r, c, w, kept)
    if (allocated(this % spatial_discretization_stencil)) call stencil_triples(this % spatial_discretization_stencil, is_carried, r, c, w, kept)

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
    if (size(variations) < 1 .or. size(variations) > this % max_degree()) then
       error stop 'gti_block: the requested order is within max_degree'
    end if
    call state_of(this, input_data, input_graph, x, state)

    ! THE SECOND PARTIAL. The time discretization stencil rows, the spatial discretization stencil and the
    ! carried rows are linear in the state and read no design, so
    ! only the physics has one: its second partial at every point,
    ! along both directions, on the row its primary degree holds.
    if (size(variations) == 2) then
       call second_tangent(this, input_data, variations, x, governing)
       allocate(r(this % num_unknowns()), source=0.0_dp)
       call placed(this, governing, r)
       call carry_held(this, r)
       call placed_output(this, input_graph, r, output)
       return
    end if

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

    call this % time_discretization_stencil % partial_action(input_graph, [state], &
         & [variations(1) % with_argument(this % time_discretization_stencil % argument(1))], half)
    call half % real_vector(r)

    if (allocated(this % spatial_discretization_stencil)) then
       call this % spatial_discretization_stencil % partial_action(input_graph, [state], &
            & [variations(1) % with_argument(this % spatial_discretization_stencil % argument(1))], half)
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
  ! The physics' second partial along two variations, each on the
  ! argument its variation names: a state direction is gathered over
  ! the points, a design direction is read as given.
  !===================================================================!

  subroutine second_tangent(this, input_data, variations, x, governing)

    class(block_residual), intent(in) :: this
    class(field)         , intent(in) :: input_data(:)
    type(variation)      , intent(in) :: variations(:)
    real(dp)             , intent(in) :: x(:)
    real(dp), allocatable, intent(out) :: governing(:)

    type(stored_field), allocatable :: inputs(:)
    type(variation) :: at_points(2)
    class(field), allocatable :: half
    integer :: i

    call point_inputs(this, input_data, x, inputs)
    do i = 1, 2
       at_points(i) = physics_variation(this, variations(i))
    end do
    call this % physics % partial_action(this % points, inputs, at_points, half)
    call half % real_vector(governing)

  end subroutine second_tangent

  function physics_variation(this, given) result(at_points)

    class(block_residual), intent(in) :: this
    type(variation)      , intent(in) :: given
    type(variation) :: at_points

    type(stored_field) :: direction
    real(dp), allocatable :: v(:)

    call given % direction(v)
    if (given % argument_is(this % argument(1))) then
       direction = stored_field('direction', this % points % vertex_set(), &
            & size(this % at) * this % degrees)
       call direction % set_real_vector(gathered(this, v))
       at_points = variation(this % physics % argument(1), direction)
    else if (given % argument_is(this % argument(2))) then
       at_points = given % with_argument(this % physics % argument(2))
    else
       error stop 'gti_block: a variation names the state or the design'
    end if

  end function physics_variation

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

  !===================================================================!
  ! THE ORDER A LEVEL IS SWEPT IN, derived and not declared. Swept by
  ! instants, the members are coupled by the time discretization stencil rows - a row at
  ! one instant reading a point at another - and that coupling is
  ! acyclic for any march, since every scheme reads backward; its
  ! loop is the sweep. The transposed block's pattern is the same
  ! graph read the other way, so it sweeps from the last instant by
  ! the same rule and no one says so. Swept by nodes, the
  ! coupling is the spatial discretization stencil's and symmetric, so no node is before
  ! another and they are swept as they lie.
  !===================================================================!

  subroutine member_order(this, by_instants, order)

    class(block_residual), intent(in)  :: this
    logical              , intent(in)  :: by_instants
    integer, allocatable , intent(out) :: order(:)

    integer, allocatable :: table(:,:), label(:), slice(:), node(:), moment(:)
    type(stored_directed_graph) :: coupling
    integer :: ne, e, n, t, h, k, members

    call this % labels_of(slice, node, moment)

    ! the space level's members couple both ways through the mesh,
    ! which has no loop, and are swept as numbered
    if (.not. by_instants) then
       order = [(k, k = 1, maxval(node))]
       return
    end if

    ! the time level's coupling: a time discretization stencil row at one member that
    ! reads an unknown at another, which for every family looks one way
    label   = slice
    members = maxval(label)
    ne      = this % time_discretization_stencil % pattern % num_edges()
    allocate(table(2, ne))
    n = 0
    do e = 1, ne
       t = label(this % time_discretization_stencil % pattern % edge_tail(e))
       h = label(this % time_discretization_stencil % pattern % edge_head(e))
       if (t == h) cycle
       n = n + 1
       table(:, n) = [t, h]
    end do
    coupling = stored_directed_graph(members, tails=table(1, 1:n), heads=table(2, 1:n))
    order    = coupling % loop(forward)

  end subroutine member_order

end module gti_block
