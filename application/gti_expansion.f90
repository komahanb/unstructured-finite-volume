!=====================================================================!
! The expansion: the graph a time integrator is, with what hangs on
! it.
!
! One object, one identity. Its branches are the six levels -
! expansion, sweep, horizon, block, slice, component - and everything
! that is not structure is kept in a map keyed on a node's identity:
! what it is called, what it holds, and how many members the set it
! denotes has. That is the arrangement view_mesh uses, one step
! further: a mesh attaches its measurements as components because a
! mesh has a fixed set of them, and the levels here do not, so they
! are attached by identity instead.
!
!             THE SLOTS
!
! Handed a physics, one family per block, the instants each block
! covers, a grid and a design, this builds a hierarchy in which every
! level is consistent and every coupling is relationally valid. The
! slots are the only things it is told; nothing else about the
! problem is written here.
!
!             WHAT IS NOT ASSIGNABLE
!
! The storage lends pointers into its own nodes, so an expansion is
! refused assignment: a copy would share them and either release
! would strand the other. That is why it is built into a variable
! rather than returned from one - a constructor's result would be
! assigned, and the assignment is what is refused.
!
!             THE ROWS A BLOCK CARRIES
!
! A block holds only the rows that fit inside it. A row on the d-th
! derivative at instant k reads a fixed number of instants back, and
! at the first instants of a block there are not that many, so those
! rows are absent and their components are carried in instead: they
! are marked as holding a value from the start, and every component
! after them waits on a march. Joining one block's last instants to
! the next block's first is the junction constraint, which is not
! built here.
!
!             WHAT IS REFUSED
!
! A block with fewer instants than its family reaches; a design the
! grid cannot read; an expansion built twice. Every level is checked
! as it is assembled, so a coupling that names another level's
! members stops the build where it is made.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_expansion

  use util_precision  , only : dp
  use graph_fractal         , only : graph, branch
  use view_level            , only : level_storage, level_consistent, &
       & level_is_leaf, level_members, level_couples, level_coupling
  use view_sequence         , only : sequence_empty, sequence_first, sequence_rest
  use view_relational       , only : relational_binding, relational_valid, &
       & num_relations, relation_at
  use relation_finitary     , only : relation
  use relation_binary       , only : csr_relation, binary_relation
  use map_value             , only : value_map, VALUE_UNKNOWN, VALUE_KNOWN
  use map_label             , only : label_map
  use map_set               , only : set_map
  use map_set_representation, only : counted_set_representation
  use view_directed_stored  , only : stored_directed_graph
  use field_calculus        , only : field
  use field_stored          , only : stored_field
  use operation_family      , only : family
  use operation_grid        , only : grid
  use operation_coupling    , only : weights_of
  use operation_weight      , only : scheme_weight
  use physics_integrand     , only : nodal_integrand

  implicit none

  private
  public :: expansion, family_holder
  public :: marches_by_stages
  public :: block_reach

  !===================================================================!
  ! One family per block. Families of different kinds cannot share an
  ! array, so each is held in its own allocatable slot.
  !===================================================================!

  type :: family_holder
     class(family), allocatable :: scheme
  end type family_holder

  type :: expansion

     type(level_storage)     , private :: nodes
     type(label_map)         , private :: labels
     type(value_map)         , private :: values
     type(set_map)           , private :: extents
     type(relational_binding), private :: bindings
     integer                 , private :: root_at = 0
     integer                 , private :: degrees = 0

   contains

     procedure :: build
     procedure :: root
     procedure :: node
     procedure :: num_nodes
     procedure :: label_of
     procedure :: status_of
     procedure :: value_of
     procedure :: extent_of
     procedure :: consistent
     procedure :: tuples_of
     procedure, private :: refuse_assignment
     generic :: assignment(=) => refuse_assignment

  end type expansion

contains

  !===================================================================!
  ! An expansion lends pointers into its own storage, so a copy would
  ! share them.
  !===================================================================!

  subroutine refuse_assignment(lhs, rhs)

    class(expansion), intent(out) :: lhs
    class(expansion), intent(in)  :: rhs

    associate (u1 => lhs, u2 => rhs); end associate

    error stop 'gti_expansion: an expansion is not assignable'

  end subroutine refuse_assignment

  pure integer function root(this)

    class(expansion), intent(in) :: this

    root = this % root_at

  end function root

  function node(this, at) result(g)

    class(expansion), intent(in) :: this
    integer         , intent(in) :: at
    type(graph), pointer :: g

    g => this % nodes % node(at)

  end function node

  pure integer function num_nodes(this)

    class(expansion), intent(in) :: this

    num_nodes = this % nodes % num_nodes()

  end function num_nodes

  function label_of(this, g) result(text)

    class(expansion), intent(in) :: this
    type(graph)     , intent(in) :: g
    character(len=:), allocatable :: text

    text = ''
    if (this % labels % labelled(g)) text = this % labels % label_of(g)

  end function label_of

  pure integer function status_of(this, g)

    class(expansion), intent(in) :: this
    type(graph)     , intent(in) :: g

    status_of = this % values % status_of(g)

  end function status_of

  subroutine value_of(this, g, x)

    class(expansion)     , intent(in)  :: this
    type(graph)          , intent(in)  :: g
    real(dp), allocatable, intent(out) :: x(:)

    call this % values % value_of(g, x)

  end subroutine value_of

  !===================================================================!
  ! How many members the set a node denotes has, or zero where no
  ! extent was recorded.
  !===================================================================!

  integer function extent_of(this, g) result(n)

    class(expansion), intent(in) :: this
    type(graph)     , intent(in) :: g

    n = 0
    if (this % extents % describes(g)) n = this % extents % num_members_of(g)

  end function extent_of

  !===================================================================!
  ! The tuples of a coupling's relation, in the relation's own order,
  ! which is the order the coupling's value holds its weights in.
  ! Invalid input: a node that is not a coupling of one relation.
  !===================================================================!

  subroutine tuples_of(this, coupling, table)

    class(expansion), intent(in) :: this
    type(graph)     , intent(in) :: coupling
    integer, allocatable, intent(out) :: table(:,:)

    class(relation), pointer :: r

    if (num_relations(coupling) /= 1) then
       error stop 'gti_expansion: a coupling holds one relation'
    end if
    r => relation_at(coupling, this % bindings, 1)
    select type (r)
    class is (binary_relation)
       call r % tuples(table)
    class default
       error stop 'gti_expansion: a coupling''s relation is binary'
    end select

  end subroutine tuples_of

  !===================================================================!
  ! The weights of a coupling in the order its relation holds the
  ! tuples: the relation groups them by source and keeps each once,
  ! so a weight computed per tuple as given is placed where the
  ! relation put its tuple. Invalid input: a tuple given twice, which
  ! would leave one weight with no place.
  !===================================================================!

  function in_relation_order(this, coupling, table, w) result(placed)

    class(expansion), intent(in) :: this
    type(graph)     , intent(in) :: coupling
    integer         , intent(in) :: table(:,:)
    real(dp)        , intent(in) :: w(:)
    real(dp), allocatable :: placed(:)

    integer, allocatable :: kept(:,:), at(:,:)
    integer :: e, n

    call this % tuples_of(coupling, kept)
    if (size(kept, 2) /= size(table, 2)) then
       error stop 'gti_expansion: a coupling names each tuple once'
    end if
    n = max(maxval(table(1, :)), maxval(table(2, :)))
    allocate(at(n, n), source=0)
    do e = 1, size(table, 2)
       at(table(1, e), table(2, e)) = e
    end do
    allocate(placed(size(kept, 2)))
    do e = 1, size(kept, 2)
       placed(e) = w(at(kept(1, e), kept(2, e)))
    end do

  end function in_relation_order

  !===================================================================!
  ! THE BUILD.
  !===================================================================!

  subroutine build(this, physics, schemes, instants, steps, &
       & max_derivative_degree, design)

    class(expansion)      , intent(inout) :: this
    class(nodal_integrand), intent(in)    :: physics
    type(family_holder)   , intent(in)    :: schemes(:)
    integer               , intent(in)    :: instants(:)
    class(grid)           , intent(in)    :: steps
    integer               , intent(in)    :: max_derivative_degree
    real(dp)              , intent(in)    :: design(:)

    real(dp), allocatable :: dt(:)
    integer , allocatable :: sweeps(:)
    integer :: s

    if (this % root_at /= 0) then
       error stop 'gti_expansion: an expansion is built once'
    end if
    if (size(schemes) /= size(instants)) then
       error stop 'gti_expansion: one family and one instant count per block'
    end if

    this % degrees = physics % equation_degree() + 1
    call partition(steps, sum(instants), design, dt)

    allocate(sweeps(max_derivative_degree + 1))
    do s = 0, max_derivative_degree
       sweeps(s + 1) = one_sweep(this, physics, schemes, instants, dt, s)
    end do

    this % root_at = this % nodes % assemble(sweeps, 0)
    call this % labels % bind(this % node(this % root_at), &
         & 'expansion of ' // physics % name() // ' in the design')
    call this % extents % bind(this % node(this % root_at), &
         & counted_set_representation(max(size(design), 1)))
    call attach_known(this, this % root_at, design)

  end subroutine build

  !===================================================================!
  ! The steps, from the grid, over every instant of the horizon.
  !===================================================================!

  subroutine partition(steps, num_instants, design, dt)

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

  end subroutine partition

  pure function padded(design) result(x)

    real(dp), intent(in) :: design(:)
    real(dp), allocatable :: x(:)

    if (size(design) == 0) then
       allocate(x(1), source=0.0_dp)
    else
       x = design
    end if

  end function padded

  !===================================================================!
  ! A row attached and marked in one step.
  !===================================================================!

  subroutine attach_known(this, at, x)

    class(expansion), intent(inout) :: this
    integer         , intent(in)    :: at
    real(dp)        , intent(in)    :: x(:)

    call this % values % attach_unknown(this % node(at))
    call this % values % mark_known(this % node(at), padded(x))

  end subroutine attach_known

  !===================================================================!
  ! One sweep, over a horizon of its own. Sweep zero is the primal
  ! and the sweeps above it are the tangents; each owns its horizon,
  ! because the value rows are keyed on identity and one sweep's
  ! state must not stand for another's.
  !===================================================================!

  integer function one_sweep(this, physics, schemes, instants, dt, sensitivity) result(at)

    class(expansion)      , intent(inout) :: this
    class(nodal_integrand), intent(in)    :: physics
    type(family_holder)   , intent(in)    :: schemes(:)
    integer               , intent(in)    :: instants(:)
    real(dp)              , intent(in)    :: dt(:)
    integer               , intent(in)    :: sensitivity

    at = this % nodes % assemble([one_horizon(this, physics, schemes, instants, dt)], 0)

    if (sensitivity == 0) then
       call this % labels % bind(this % node(at), 'sweep 0, the functional itself')
    else
       call this % labels % bind(this % node(at), 'sweep ' // written(sensitivity) // &
            & ', derivative ' // written(sensitivity) // ' in the design')
    end if
    call this % values % attach_unknown(this % node(at))

  end function one_sweep

  !===================================================================!
  ! One horizon: the blocks that partition the instants, in order.
  !===================================================================!

  integer function one_horizon(this, physics, schemes, instants, dt) result(at)

    class(expansion)      , intent(inout) :: this
    class(nodal_integrand), intent(in)    :: physics
    type(family_holder)   , intent(in)    :: schemes(:)
    integer               , intent(in)    :: instants(:)
    real(dp)              , intent(in)    :: dt(:)

    integer, allocatable :: blocks(:)
    integer :: b, first

    allocate(blocks(size(instants)))
    first = 1

    do b = 1, size(instants)
       blocks(b) = one_block(this, physics, schemes(b) % scheme, &
            & first, first + instants(b) - 1, dt)
       first = first + instants(b)
    end do

    at = this % nodes % assemble(blocks, 0)
    call this % labels % bind(this % node(at), 'horizon of duration ' // &
         & written(sum(dt)))

  end function one_horizon

  !===================================================================!
  ! One block: its slices, the family that marches them, its steps,
  ! and the coupling that carries the scheme reach.
  !===================================================================!

  integer function one_block(this, physics, scheme, first, last, dt) result(at)

    class(expansion)      , intent(inout) :: this
    class(nodal_integrand), intent(in)    :: physics
    class(family)         , intent(in)    :: scheme
    integer               , intent(in)    :: first, last
    real(dp)              , intent(in)    :: dt(:)

    integer, allocatable :: slices(:)
    integer :: k, coupling

    if (last - first + 1 <= scheme % history_depth(this % degrees - 1)) then
       error stop 'gti_expansion: a block holds more instants than its family reaches'
    end if

    allocate(slices(last - first + 1))
    do k = first, last
       slices(k - first + 1) = one_slice(this, physics, scheme, k, first, dt(k))
    end do

    if (marches_by_stages(scheme, this % degrees)) then
       coupling = carry_coupling(this, scheme, slices, first, last, dt)
    else
       coupling = block_coupling(this, physics, scheme, slices, first, last, dt)
    end if
    at       = this % nodes % assemble(slices, coupling)

    call this % labels % bind(this % node(at), scheme % name() // ' block')
    call attach_known(this, at, dt(first:last))

  end function one_block

  !===================================================================!
  ! One slice: the components of every degree at one instant.
  !===================================================================!

  integer function one_slice(this, physics, scheme, instant, first, step) result(at)

    class(expansion)      , intent(inout) :: this
    class(nodal_integrand), intent(in)    :: physics
    class(family)         , intent(in)    :: scheme
    integer               , intent(in)    :: instant, first
    real(dp)              , intent(in)    :: step

    associate (u1 => physics); end associate

    if (marches_by_stages(scheme, this % degrees)) then
       at = staged_slice(this, scheme, instant, first, step)
    else
       at = plain_slice(this, scheme, instant, first)
    end if

    call this % labels % bind(this % node(at), 'slice at instant ' // written(instant))

  end function one_slice

  !===================================================================!
  ! A slice of a family whose rows run between instants: the
  ! components of every degree, and no stage level at all.
  !===================================================================!

  integer function plain_slice(this, scheme, instant, first) result(at)

    class(expansion), intent(inout) :: this
    class(family)   , intent(in)    :: scheme
    integer         , intent(in)    :: instant, first

    integer, allocatable :: components(:)
    integer :: d

    allocate(components(this % degrees))
    do d = 0, this % degrees - 1
       components(d + 1) = one_component(this, scheme, instant, first, d)
    end do

    at = this % nodes % assemble(components, 0)

  end function plain_slice

  !===================================================================!
  ! Whether a family's derived rows run between the stages of one
  ! step rather than between instants. A stage family gives an empty
  ! pattern at every degree, which is how it says its rows are not
  ! offsets back through the instants; that question is asked here,
  ! so no family declares its kind twice.
  !===================================================================!

  logical function marches_by_stages(scheme, nd) result(staged)

    class(family), intent(in) :: scheme
    integer      , intent(in) :: nd

    integer, allocatable :: offset(:), degrees_of(:)
    integer :: d, primary

    primary = scheme % primary_degree(nd - 1)
    staged  = .true.

    do d = 0, nd - 1
       if (d == primary) cycle
       call scheme % row_pattern(d, nd - 1, offset, degrees_of)
       if (size(offset) > 0) then
          staged = .false.
          return
       end if
    end do

  end function marches_by_stages

  !===================================================================!
  ! A slice of a stage family: the stages of the step arriving at
  ! this instant, then the instant itself, which is the step's
  ! closing evaluation and holds components like any stage. The
  ! block's first instant has no step arriving at it and so holds no
  ! stages; its components are the initial data.
  !===================================================================!

  integer function staged_slice(this, scheme, instant, first, step) result(at)

    class(expansion), intent(inout) :: this
    class(family)   , intent(in)    :: scheme
    integer         , intent(in)    :: instant, first
    real(dp)        , intent(in)    :: step

    integer, allocatable :: members(:)
    integer :: s, i

    s = scheme % num_stages()

    if (instant == first) then
       at = this % nodes % assemble([stage_node(this, scheme, instant, first, 0)], 0)
       return
    end if

    allocate(members(s + 1))
    do i = 1, s
       members(i) = stage_node(this, scheme, instant, first, i)
    end do
    members(s + 1) = stage_node(this, scheme, instant, first, 0)

    at = this % nodes % assemble(members, &
         & slice_coupling(this, scheme, members, s, step))

  end function staged_slice

  !===================================================================!
  ! One stage, or the arriving instant when the index is zero: the
  ! components of every degree held at one evaluation point.
  !===================================================================!

  integer function stage_node(this, scheme, instant, first, index) result(at)

    class(expansion), intent(inout) :: this
    class(family)   , intent(in)    :: scheme
    integer         , intent(in)    :: instant, first, index

    integer, allocatable :: components(:)
    integer :: d

    allocate(components(this % degrees))
    do d = 0, this % degrees - 1
       components(d + 1) = one_component(this, scheme, instant, first, d)
    end do

    at = this % nodes % assemble(components, 0)

    if (index == 0) then
       call this % labels % bind(this % node(at), 'the arriving instant')
    else
       call this % labels % bind(this % node(at), 'stage ' // written(index))
    end if

  end function stage_node

  !===================================================================!
  ! One component: a leaf holding its freedoms. The instants a block
  ! reaches back over carry their values from the start; every
  ! component after them waits on a march.
  !===================================================================!

  integer function one_component(this, scheme, instant, first, degree) result(at)

    class(expansion), intent(inout) :: this
    class(family)   , intent(in)    :: scheme
    integer         , intent(in)    :: instant, first, degree

    at = this % nodes % assemble([integer ::], 0)

    call this % labels % bind(this % node(at), 'component of degree ' // written(degree))
    call this % extents % bind(this % node(at), counted_set_representation(1))

    if (instant - first < scheme % history_depth(this % degrees - 1)) then
       call attach_known(this, at, [0.0_dp])
    else
       call this % values % attach_unknown(this % node(at))
    end if

  end function one_component

  !===================================================================!
  ! An integer and a real as text, for the labels.
  !===================================================================!

  function written(n) result(text)

    class(*), intent(in) :: n
    character(len=:), allocatable :: text

    character(len=24) :: buffer

    select type (n)
    type is (integer)
       write(buffer,'(i0)') n
    type is (real(dp))
       write(buffer,'(f0.4)') n
    class default
       error stop 'gti_expansion: a label is written from a number'
    end select

    text = trim(buffer)

  end function written

  !===================================================================!
  ! THE COUPLING OF A BLOCK.
  !
  ! How many edges the rows that fit inside this block make, and, on
  ! the second pass, what they are. A row on degree d at local
  ! instant kk fits when every source it reads lies at or after the
  ! block's first instant.
  !===================================================================!

  subroutine block_reach(scheme, nd, n, tails, heads, source_degree, determines)

    class(family), intent(in) :: scheme
    integer      , intent(in) :: nd, n
    integer, allocatable, intent(out) :: tails(:), heads(:)
    integer, allocatable, intent(out) :: source_degree(:), determines(:)

    integer, allocatable :: offset(:), degrees_of(:)
    integer :: primary, kk, d, e, counted, at, pass

    primary = scheme % primary_degree(nd - 1)

    do pass = 1, 2
       counted = 0
       do kk = 1, n
          do d = 0, nd - 1
             if (d == primary) cycle
             call scheme % row_pattern(d, nd - 1, offset, degrees_of)
             if (size(offset) == 0) cycle
             if (kk - maxval(offset) < 1) cycle
             do e = 1, size(offset)
                counted = counted + 1
                if (pass == 2) then
                   at = counted
                   tails(at)         = kk - offset(e)
                   heads(at)         = kk
                   source_degree(at) = degrees_of(e)
                   determines(at)    = d
                end if
             end do
          end do
       end do
       if (pass == 1) allocate(tails(counted), heads(counted), &
            & source_degree(counted), determines(counted))
    end do

  end subroutine block_reach

  !===================================================================!
  ! The weights the reach carries, from the family and the steps of
  ! this block alone: a row that fits reads no instant before the
  ! block's first, so the block's own steps are all it needs.
  !===================================================================!

  subroutine reach_weights(scheme, n, tails, heads, source_degree, determines, dt, w)

    class(family), intent(in) :: scheme
    integer      , intent(in) :: n, tails(:), heads(:), source_degree(:), determines(:)
    real(dp)     , intent(in) :: dt(:)
    real(dp), allocatable, intent(out) :: w(:)

    call weights_of(scheme_weight(scheme), n, tails, heads, dt, source_degree, determines, w)

  end subroutine reach_weights

  !===================================================================!
  ! The coupling itself: this block's slices as carriers, then the
  ! two sets the relation runs between, then the relation. The
  ! weights are the coupling's own value, so the sparsity and the
  ! numbers hang on one identity.
  !===================================================================!

  integer function block_coupling(this, physics, scheme, slices, first, last, dt) result(at)

    class(expansion)      , intent(inout) :: this
    class(nodal_integrand), intent(in)    :: physics
    class(family)         , intent(in)    :: scheme
    integer               , intent(in)    :: slices(:), first, last
    real(dp)              , intent(in)    :: dt(:)

    integer, allocatable :: tails(:), heads(:), source_degree(:), determines(:)
    integer, allocatable :: table(:,:)
    real(dp), allocatable :: w(:)
    integer :: n, nd, components, constraints, holder

    associate (u1 => physics); end associate

    n  = last - first + 1
    nd = this % degrees

    call block_reach(scheme, nd, n, tails, heads, source_degree, determines)
    call reach_weights(scheme, n, tails, heads, source_degree, determines, &
         & dt(first:last), w)

    components  = named_set(this, n * nd, 'the components of this block')
    constraints = named_set(this, n * nd, 'the constraint instances of this block')
    table       = tuples(nd, tails, heads, source_degree, determines)

    holder = this % nodes % assemble([integer ::], 0)
    call this % labels % bind(this % node(holder), 'the scheme reach')

    at = this % nodes % couple([slices, components, constraints], [holder])

    call bind_carriers(this, [slices, components, constraints])
    call bind_reach(this, holder, components, constraints, table)

    call this % labels % bind(this % node(at), scheme % name() // ' coupling')
    call attach_known(this, at, in_relation_order(this, this % node(at), table, w))

  end function block_coupling

  !===================================================================!
  ! Each edge as a tuple: the unknown its source is, and the unknown
  ! its constraint determines. Rows and columns share one index
  ! space, so the block is square.
  !===================================================================!

  pure function tuples(nd, tails, heads, source_degree, determines) result(table)

    integer, intent(in) :: nd, tails(:), heads(:), source_degree(:), determines(:)
    integer, allocatable :: table(:,:)

    integer :: e

    allocate(table(2, size(tails)))

    table(1,:) = [((tails(e) - 1) * nd + source_degree(e) + 1, e = 1, size(tails))]
    table(2,:) = [((heads(e) - 1) * nd + determines(e) + 1, e = 1, size(heads))]

  end function tuples

  !===================================================================!
  ! A carrier denoting a set of the given extent.
  !===================================================================!

  integer function named_set(this, n, text) result(at)

    class(expansion), intent(inout) :: this
    integer         , intent(in)    :: n
    character(len=*), intent(in)    :: text

    at = this % nodes % assemble([integer ::], 0)
    call this % labels % bind(this % node(at), text)
    call this % extents % bind(this % node(at), counted_set_representation(n))

  end function named_set

  subroutine bind_carriers(this, carriers)

    class(expansion), intent(inout) :: this
    integer         , intent(in)    :: carriers(:)

    type(graph), pointer :: g
    integer :: i

    do i = 1, size(carriers)
       g => this % nodes % node(carriers(i))
       call this % bindings % bind_set(g, g)
    end do

  end subroutine bind_carriers

  subroutine bind_reach(this, holder, components, constraints, table)

    class(expansion), intent(inout) :: this
    integer         , intent(in)    :: holder, components, constraints, table(:,:)

    type(graph), pointer :: g, from, into

    from => this % nodes % node(components)
    into => this % nodes % node(constraints)
    g    => this % nodes % node(holder)

    call this % bindings % bind_relation(g, &
         & csr_relation('scheme reach', from, into, table, this % extents))

  end subroutine bind_reach

  !===================================================================!
  ! Every level consistent and every coupling relationally valid.
  ! Each level is already checked as it is assembled; this checks the
  ! built hierarchy once more, from the outside.
  !===================================================================!

  recursive logical function consistent(this, g) result(ok)

    class(expansion), intent(in) :: this
    type(graph)     , intent(in) :: g

    type(graph), pointer :: coupling

    ok = level_consistent(g)
    if (.not. ok) return

    if (level_couples(g)) then
       coupling => level_coupling(g)
       ok = relational_valid(coupling, this % bindings)
       if (.not. ok) return
    end if

    if (level_is_leaf(g)) return
    ok = every_member(this, level_members(g))

  end function consistent

  recursive logical function every_member(this, members) result(ok)

    class(expansion), intent(in) :: this
    type(branch)    , intent(in) :: members

    type(graph), pointer :: first

    ok = .true.
    if (sequence_empty(members)) return

    first => sequence_first(members)
    ok = this % consistent(first)
    if (ok) ok = every_member(this, sequence_rest(members))

  end function every_member

  !===================================================================!
  ! THE ROWS OF ONE STEP.
  !
  ! In the numbering a stage family uses, vertex one is the instant
  ! the step leaves from, vertices two to one plus s are the stages,
  ! and the last is the instant it arrives at. The rows made here are
  ! the ones inside the step, so no edge leaves vertex one: those
  ! carry the previous instant across and belong to the block.
  !
  !      degree d below the highest
  !          stage i reads stage j at degree d+1, for j at or before i
  !          the arriving instant reads every stage at degree d+1
  !      the highest degree
  !          the arriving instant reads every stage at that degree,
  !          which is the recovery
  !===================================================================!

  subroutine stage_reach(nd, s, tails, heads, source_degree, determines)

    integer, intent(in) :: nd, s
    integer, allocatable, intent(out) :: tails(:), heads(:)
    integer, allocatable, intent(out) :: source_degree(:), determines(:)

    integer :: d, i, j, at

    at = (nd - 1) * (s * (s + 1) / 2 + s) + s
    allocate(tails(at), heads(at), source_degree(at), determines(at))
    at = 0

    do d = 0, nd - 1
       do i = 1, s
          if (d == nd - 1) cycle
          do j = 1, i
             at = at + 1
             tails(at) = 1 + j
             heads(at) = 1 + i
             source_degree(at) = d + 1
             determines(at) = d
          end do
       end do
       do j = 1, s
          at = at + 1
          tails(at) = 1 + j
          heads(at) = 2 + s
          source_degree(at) = min(d + 1, nd - 1)
          determines(at) = d
       end do
    end do

  end subroutine stage_reach

  !===================================================================!
  ! A vertex of the family's numbering as an unknown of this step:
  ! stage i is member i, the arriving instant is member s+1.
  !===================================================================!

  pure integer function stage_unknown(vertex, degree, s, nd) result(at)

    integer, intent(in) :: vertex, degree, s, nd

    integer :: member

    if (vertex == 2 + s) then
       member = s + 1
    else
       member = vertex - 1
    end if

    at = (member - 1) * nd + degree + 1

  end function stage_unknown

  !===================================================================!
  ! The weights of one step. The step is read at every vertex of the
  ! family's numbering, the stages included, because that is where
  ! the scaling reads it.
  !===================================================================!

  subroutine stage_weights(scheme, s, tails, heads, source_degree, determines, step, w)

    class(family), intent(in) :: scheme
    integer      , intent(in) :: s, tails(:), heads(:), source_degree(:), determines(:)
    real(dp)     , intent(in) :: step
    real(dp), allocatable, intent(out) :: w(:)

    call weights_of(scheme_weight(scheme), s + 2, tails, heads, spread(step, 1, s + 2), &
         & source_degree, determines, w)

  end subroutine stage_weights

  !===================================================================!
  ! The coupling of one step: its stages and its arriving instant as
  ! carriers, then the two sets the tableau runs between, then the
  ! relation. The weights are the coupling's own value.
  !===================================================================!

  integer function slice_coupling(this, scheme, members, s, step) result(at)

    class(expansion), intent(inout) :: this
    class(family)   , intent(in)    :: scheme
    integer         , intent(in)    :: members(:), s
    real(dp)        , intent(in)    :: step

    integer, allocatable :: tails(:), heads(:), source_degree(:), determines(:)
    integer, allocatable :: table(:,:)
    real(dp), allocatable :: w(:)
    integer :: nd, components, constraints, holder, e

    nd = this % degrees

    call stage_reach(nd, s, tails, heads, source_degree, determines)
    call stage_weights(scheme, s, tails, heads, source_degree, determines, step, w)

    components  = named_set(this, (s + 1) * nd, 'the components of this step')
    constraints = named_set(this, (s + 1) * nd, 'the constraint instances of this step')

    allocate(table(2, size(tails)))
    table(1,:) = [(stage_unknown(tails(e), source_degree(e), s, nd), e = 1, size(tails))]
    table(2,:) = [(stage_unknown(heads(e), determines(e), s, nd), e = 1, size(heads))]

    holder = this % nodes % assemble([integer ::], 0)
    call this % labels % bind(this % node(holder), 'the butcher reach')

    at = this % nodes % couple([members, components, constraints], [holder])

    call bind_carriers(this, [members, components, constraints])
    call bind_reach(this, holder, components, constraints, table)

    call this % labels % bind(this % node(at), scheme % name() // ' stage coupling')
    call attach_known(this, at, in_relation_order(this, this % node(at), table, w))

  end function slice_coupling

  !===================================================================!
  ! THE CARRY BETWEEN STEPS.
  !
  ! Where a block's rows run inside its steps, what joins one step to
  ! the next is the instant they share: every row below the highest
  ! degree reads the previous instant's own degree, carried across
  ! unchanged. Those are the edges this coupling holds, and the
  ! family gives their weight rather than this module assuming it.
  !===================================================================!

  pure integer function slice_base(kk, s, nd) result(at)

    integer, intent(in) :: kk, s, nd

    if (kk == 1) then
       at = 0
    else
       at = (1 + (kk - 2) * (s + 1)) * nd
    end if

  end function slice_base

  pure integer function closing_instant(kk, s, nd) result(at)

    integer, intent(in) :: kk, s, nd

    if (kk == 1) then
       at = slice_base(kk, s, nd)
    else
       at = slice_base(kk, s, nd) + s * nd
    end if

  end function closing_instant

  subroutine carry_reach(this, s, n, table, sources)

    class(expansion), intent(in) :: this
    integer         , intent(in) :: s, n
    integer, allocatable, intent(out) :: table(:,:)
    integer, allocatable, intent(out) :: sources(:)

    integer :: nd, kk, d, m, counted, pass, from, into

    nd = this % degrees

    do pass = 1, 2
       counted = 0
       do kk = 2, n
          from = closing_instant(kk - 1, s, nd)
          do d = 0, nd - 2
             do m = 1, s + 1
                counted = counted + 1
                into = slice_base(kk, s, nd) + (m - 1) * nd + d + 1
                if (pass == 2) then
                   table(1, counted) = from + d + 1
                   table(2, counted) = into
                   sources(counted)  = merge(2 + s, 1 + m, m == s + 1)
                end if
             end do
          end do
       end do
       if (pass == 1) then
          allocate(table(2, counted), sources(counted))
       end if
    end do

  end subroutine carry_reach

  integer function carry_coupling(this, scheme, slices, first, last, dt) result(at)

    class(expansion), intent(inout) :: this
    class(family)   , intent(in)    :: scheme
    integer         , intent(in)    :: slices(:), first, last
    real(dp)        , intent(in)    :: dt(:)

    integer, allocatable :: table(:,:), sources(:)
    real(dp), allocatable :: w(:)
    integer :: n, nd, s, unknowns, components, constraints, holder, e

    n  = last - first + 1
    nd = this % degrees
    s  = scheme % num_stages()

    call carry_reach(this, s, n, table, sources)

    unknowns = (1 + (n - 1) * (s + 1)) * nd
    call stage_weights(scheme, s, [(1, e = 1, size(sources))], sources, &
         & [((mod(table(1, e) - 1, nd)), e = 1, size(sources))], &
         & [((mod(table(2, e) - 1, nd)), e = 1, size(sources))], dt(first), w)

    components  = named_set(this, unknowns, 'the components of this block')
    constraints = named_set(this, unknowns, 'the constraint instances of this block')

    holder = this % nodes % assemble([integer ::], 0)
    call this % labels % bind(this % node(holder), 'the carry between steps')

    at = this % nodes % couple([slices, components, constraints], [holder])

    call bind_carriers(this, [slices, components, constraints])
    call bind_reach(this, holder, components, constraints, table)

    call this % labels % bind(this % node(at), scheme % name() // ' carry coupling')
    call attach_known(this, at, in_relation_order(this, this % node(at), table, w))

  end function carry_coupling

end module gti_expansion
