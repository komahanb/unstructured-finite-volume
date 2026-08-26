!=====================================================================!
! Building one block and solving it.
!
! The pieces are the same ones the assembly uses - the rows a family
! reaches over, the weights on them, the stencil they make, and the
! statement that adds the governing and carried rows to it - gathered
! here so that a caller marching a block and a caller differentiating
! one write them once.
!
!             WHERE THE BLOCKS OF A HORIZON SIT
!
! horizon_bounds says which instants each block of a chain spans. A
! block reaches back over instants that begin before it does, so
! every block after the first overlaps what came before it by
! exactly what its family reaches. What is done with that overlap -
! the junction, and the layouts either side of it - belongs to
! gti_chain, which marches them.
!
! A block must add more instants than its family reaches back over,
! or it would consist of nothing but what it was given.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_march

  use util_precision  , only : dp, least_kind_for
  use iso_fortran_env , only : real128
  use operation_coupling      , only : weights_of, weights_varied
  use gti_configuration       , only : refuse_unknown
  use operation_weight        , only : scheme_weight
  use view_directed_stored    , only : stored_directed_graph
  use view_directed           , only : directed_graph
  use field_calculus          , only : field
  use field_stored            , only : stored_field
  use operation_action      , only : variation
  use operation_stencil       , only : stencil
  use operation_newton        , only : newton
  use operation_minimization  , only : minimizer, relative, absolute, &
       & by_count, by_rate
  use operation_dense_direct  , only : dense_direct
  use operation_gmres         , only : gmres
  use operation_family        , only : family
  use operation_grid          , only : grid, uniform_grid
  use operation_weight        , only : scheme_weight
  use operation_scheme_stencil, only : derived_constraints
  use physics_integrand       , only : nodal_integrand
  use gti_expansion           , only : block_reach, family_holder
  use gti_block               , only : block_residual
  use gti_sweeps              , only : jacobian_of, assembly_present, multigrid_on, &
       & set_aggregates, coarse_nodes, take_inner, keep_inner, forget_inner
  use util_tally              , only : tally_record, tangent_loops, adjoint_loops

  implicit none

  !-------------------------------------------------------------------!
  ! WHAT AN UNCONVERGED MARCH LEFT, by aspect. The imbalance is a
  ! vector with one entry per unknown, and the unknowns lie in slots
  ! of one instant's - or one stage's - components each, so its norm
  ! splits exactly, ||r||^2 = sum over slots and degrees of r^2, and
  ! its steepest direction in the state is the gradient of the norm,
  ! d||r||/dq = A^T r / ||r||, one transposed matvec. The largest entry
  ! of each names where the imbalance sits and which state drives it.
  !-------------------------------------------------------------------!

  type :: imbalance

     logical  :: converged = .true.
     logical  :: diverging = .false.
     real(dp) :: norm      = 0.0_dp
     real(dp) :: began     = 0.0_dp
     real(dp), allocatable :: by_degree(:)
     integer  :: worst_slot = 0, worst_degree = 0
     integer  :: steepest_slot = 0, steepest_degree = 0
     real(dp) :: steepest = 0.0_dp

  end type imbalance

  !-------------------------------------------------------------------!
  ! WHICH LEVEL IS SWEPT. The block is one nonlinear statement over
  ! every instant, node and component; solving it is a choice of
  ! which level's members are solved exactly inside and which level
  ! is swept over them with the rest held:
  !
  !      space-time    no level: the whole block at once
  !      time          the instants, in order: each instant's nodes
  !                    and components solved with the instants before
  !                    it held - the classical step, exact in one pass
  !                    since every scheme looks backward
  !      space         the nodes: each node's whole history solved
  !                    with its neighbours' histories held, the sweep
  !                    repeated until the coupling agrees
  !
  ! The same fixed point in all three. Nothing below this loop knows
  ! which was chosen.
  !-------------------------------------------------------------------!

  character(len=16), save :: sweep_level = 'space-time'

  ! Stamps handed out to statements, so that a direct solver can tell
  ! a statement it has factorised from a new one.
  integer, save :: stamps_given = 0

  !-------------------------------------------------------------------!
  ! HOW A MARCH STOPS. A caller that sets nothing gets a tolerance
  ! measured against the imbalance the march began at, and a budget
  ! taken from the rate the march itself shows. The count is a
  ! backstop and not the operative limit.
  !-------------------------------------------------------------------!

  real(dp), save :: stopping_tolerance  = 1.0e-12_dp
  integer , save :: stopping_criterion  = relative
  integer , save :: stopping_budget     = by_rate
  integer , save :: stopping_iterations = 100

  private
  public :: partition, partitioned, scheme_rows, block_of, solved, unknowns_graph, step_partials
  public :: step_second_partials
  public :: unknown, consistent_states, frozen_inputs
  public :: set_stopping
  public :: consistent_state
  public :: imbalance
  public :: swept, set_sweep, sweep_named
  public :: solved_linear, by_tangent, by_adjoint, fresh_stamp
  public :: weight_of, precision_needed
  public :: horizon_bounds

contains

  !===================================================================!
  ! THE WEIGHT A BLOCK CARRIES: ||A||_inf from the family and the step,
  ! without the matrix. A row determining degree d reads the sources
  ! its pattern names, each weighted alpha dt^(sigma - d) by the same
  ! scheme_weight that builds the block, so the row's absolute sum is
  ! one apply on a coupling of that pattern at the step given, plus
  ! the one the row carries on the column it determines. The largest
  ! over the degrees is the norm. A stage family has no pattern in
  ! instants and its rows lie within one step: the incoming instant
  ! and the stages at or before, read the same way.
  !===================================================================!

  real(dp) function weight_of(scheme, degrees, step) result(w)

    class(family), intent(in) :: scheme
    integer      , intent(in) :: degrees
    real(dp)     , intent(in) :: step

    integer, allocatable :: offset(:), source_degree(:)
    integer :: d, reach, s, i, k
    logical :: any_pattern

    w = 1.0_dp
    any_pattern = .false.

    do d = 0, degrees - 1
       call scheme % row_pattern(d, degrees - 1, offset, source_degree)
       if (size(offset) == 0) cycle
       any_pattern = .true.
       reach = maxval(offset)
       w = max(w, 1.0_dp + row_weight(scheme, reach + 1, &
            & [(reach + 1 - offset(k), k = 1, size(offset))], reach + 1, &
            & source_degree, d, step))
    end do

    if (any_pattern) return

    s = scheme % num_stages()
    do d = 0, degrees - 2
       do i = 1, s
          w = max(w, 1.0_dp + row_weight(scheme, s + 2, &
               & [1, (1 + k, k = 1, i)], 1 + i, &
               & [d, (d + 1, k = 1, i)], d, step))
       end do
       w = max(w, 1.0_dp + row_weight(scheme, s + 2, &
            & [1, (1 + k, k = 1, s)], s + 2, &
            & [d, (d + 1, k = 1, s)], d, step))
    end do

  end function weight_of

  real(dp) function row_weight(scheme, num_vertices, tails, head, source_degree, &
       & determines, step) result(total)

    class(family), intent(in) :: scheme
    integer      , intent(in) :: num_vertices, tails(:), head, source_degree(:), determines
    real(dp)     , intent(in) :: step

    real(dp), allocatable :: c(:)
    integer :: k

    call weights_of(scheme_weight(scheme), num_vertices, tails, [(head, k = 1, size(tails))], &
         & [(step, k = 1, num_vertices)], source_degree, [(determines, k = 1, size(tails))], c)

    total = sum(abs(c))

  end function row_weight

  !===================================================================!
  ! THE PRECISION A TARGET NEEDS. The floor a march reaches is
  ! eps ||A|| ||q||, so a target is reachable at a kind whose spacing
  ! is under target / (||A|| ||q||). The target is the tolerance times
  ! the starting imbalance where the criterion is relative, and the
  ! tolerance itself where it is absolute.
  !===================================================================!

  subroutine precision_needed(weight, state_size, began, spacing_needed, least_kind)

    real(dp)        , intent(in)  :: weight, state_size, began
    real(real128)   , intent(out) :: spacing_needed
    character(len=:), allocatable, intent(out) :: least_kind

    real(dp) :: target

    select case (stopping_criterion)
    case (relative)
       target = stopping_tolerance * began
    case default
       target = stopping_tolerance
    end select

    spacing_needed = real(target, real128) / real(max(weight * state_size, tiny(1.0_dp)), real128)
    least_kind     = least_kind_for(spacing_needed)

  end subroutine precision_needed

  !===================================================================!
  ! THE CONSISTENT INITIAL STATE. Given the components below the
  ! highest at one instant, the highest is what the physics says it
  ! is there: q^(N) with R(q, q', ..., q^(N)) = 0, solved at that one
  ! instant with everything below it held.
  !
  ! This is the smallest block there is - one evaluation point, no
  ! scheme rows, the lower components carried and the highest the one
  ! unknown - and it is solved by the same newton as every other
  ! block. Nothing about the physics is assumed: whatever R is, its
  ! zero at the instant is what comes back. A lower vector of the
  ! wrong extent stops the program.
  !===================================================================!

  function consistent_state(physics, degrees, lower, design_value) result(q)

    class(nodal_integrand), intent(in) :: physics
    integer               , intent(in) :: degrees
    real(dp)              , intent(in) :: lower(:), design_value
    real(dp), allocatable :: q(:)

    if (size(lower) /= degrees - 1) then
       error stop 'gti_march: the components below the highest are given, and no others'
    end if

    q = consistent_states(physics, degrees, reshape(lower, [degrees - 1, 1]), design_value)

  end function consistent_state

  !===================================================================!
  ! The same at every node of a mesh: the components below the
  ! highest given at each node, lower(d + 1, i), and the highest at
  ! every node what the physics - and the level below, where one is
  ! attached - say it is there. One node is the ordinary state.
  !===================================================================!

  function consistent_states(physics, degrees, lower, design_value, spatial) result(q)

    class(nodal_integrand), intent(in)           :: physics
    integer               , intent(in)           :: degrees
    real(dp)              , intent(in)           :: lower(:,:), design_value
    type(stencil)         , intent(in), optional :: spatial
    real(dp), allocatable :: q(:)

    type(block_residual) :: rows
    type(stencil) :: none
    integer , allocatable :: carried(:), at(:)
    real(dp), allocatable :: held(:)
    real(dp) :: achieved
    integer :: i, d, nodes

    nodes = size(lower, 2)
    if (size(lower, 1) /= degrees - 1) then
       error stop 'gti_march: the components below the highest are given at every node'
    end if

    none = stencil([integer ::], [integer ::], [real(dp) ::], &
         & spread(0.0_dp, 1, nodes * degrees), 'none')

    carried = [((unknown(1, d, degrees, i, nodes), d = 0, degrees - 2), i = 1, nodes)]
    held    = [((lower(d + 1, i), d = 0, degrees - 2), i = 1, nodes)]
    at      = [(unknown(1, 0, degrees, i, nodes) - 1, i = 1, nodes)]

    rows = block_residual(none, physics, at, nodes * degrees, degrees, degrees - 1, &
         & carried, held)
    call rows % placed_in(spread(1, 1, nodes * degrees), &
         & [((i, d = 0, degrees - 1), i = 1, nodes)], spread(1, 1, nodes * degrees))
    if (present(spatial)) call rows % spatial_laid(spatial)

    call solved(rows, design_value, q, achieved)

    if (.not. achieved <= stopping_tolerance * max(1.0_dp, norm2(lower))) then
       write(*,'(a,es12.3)') ' the physics at the initial instant left a residual of ', achieved
       error stop 'gti_march: the initial state is consistent with the physics'
    end if

  end function consistent_states

  !===================================================================!
  ! A state and a design as the two inputs a block's rows read: the
  ! state on its unknowns, one design value per point.
  !===================================================================!

  subroutine frozen_inputs(q, design, num_points, unknowns, inputs)

    real(dp), intent(in) :: q(:), design
    integer , intent(in) :: num_points
    type(stored_directed_graph)    , intent(out) :: unknowns
    type(stored_field), allocatable, intent(out) :: inputs(:)

    unknowns = stored_directed_graph(size(q), tails=[integer ::], heads=[integer ::])

    allocate(inputs(2))
    inputs(1) = stored_field('state' , unknowns % vertex_set(), size(q))
    inputs(2) = stored_field('design', unknowns % vertex_set(), num_points)
    call inputs(1) % set_real_vector(q)
    call inputs(2) % set_real_vector(spread(design, 1, num_points))

  end subroutine frozen_inputs

  !===================================================================!
  ! How every march that follows stops. A criterion or a budget that
  ! is neither of its two stops the program.
  !===================================================================!

  subroutine set_stopping(tolerance, criterion, budget, iterations)

    real(dp), intent(in) :: tolerance
    integer , intent(in) :: criterion, budget, iterations

    if (tolerance <= 0.0_dp) then
       error stop 'gti_march: a tolerance is positive'
    end if
    if (criterion /= relative .and. criterion /= absolute) then
       error stop 'gti_march: a tolerance is measured relative or absolute'
    end if
    if (budget /= by_count .and. budget /= by_rate) then
       error stop 'gti_march: a budget is counted or taken from the rate'
    end if
    if (iterations < 1) then
       error stop 'gti_march: an iteration budget is positive'
    end if

    stopping_tolerance  = tolerance
    stopping_criterion  = criterion
    stopping_budget     = budget
    stopping_iterations = iterations

  end subroutine set_stopping

  !===================================================================!
  ! A uniform partition of the duration, and the instants it makes.
  !===================================================================!

  subroutine partition(duration, n, dt, t)

    real(dp), intent(in) :: duration
    integer , intent(in) :: n
    real(dp), allocatable, intent(out) :: dt(:), t(:)

    call partitioned(uniform_grid(duration), n, dt, t)

  end subroutine partition

  !===================================================================!
  ! The instants a grid makes over the duration it was given.
  !===================================================================!

  subroutine partitioned(steps, n, dt, t, design)

    class(grid), intent(in) :: steps
    integer    , intent(in) :: n
    real(dp), allocatable, intent(out) :: dt(:), t(:)
    real(dp), intent(in), optional :: design(:)

    type(stored_directed_graph) :: instants
    type(stored_field) :: knobs
    class(field), allocatable :: out
    integer :: k

    instants = stored_directed_graph(n, tails=[integer ::], heads=[integer ::])

    if (present(design)) then
       knobs = stored_field('design', instants % vertex_set(), size(design))
       call knobs % set_real_vector(design)
    else
       knobs = stored_field('design', instants % vertex_set(), 1)
       call knobs % set_real_vector([0.0_dp])
    end if

    call steps % apply(instants, [knobs], out)
    call out % real_vector(dt)

    allocate(t(n))
    t(1) = 0.0_dp
    do k = 2, n
       t(k) = t(k - 1) + dt(k)
    end do

  end subroutine partitioned

  !===================================================================!
  ! The partial of every step in every entry of a grid's design, one
  ! column per entry, read from the grid's own partial action: exact,
  ! and carrying the normalisation that keeps the steps summing to
  ! the duration.
  !===================================================================!

  subroutine step_partials(steps, n, design, v)

    class(grid), intent(in) :: steps
    integer    , intent(in) :: n
    real(dp)   , intent(in) :: design(:)
    real(dp), allocatable, intent(out) :: v(:,:)

    type(stored_directed_graph) :: instants
    type(stored_field) :: knobs, direction
    class(field), allocatable :: out
    real(dp), allocatable :: e(:), column(:)
    integer :: j

    instants = stored_directed_graph(n, tails=[integer ::], heads=[integer ::])
    knobs    = stored_field('design', instants % vertex_set(), size(design))
    call knobs % set_real_vector(design)

    allocate(v(n, size(design)), e(size(design)))
    do j = 1, size(design)
       e    = 0.0_dp
       e(j) = 1.0_dp
       direction = stored_field('direction', instants % vertex_set(), size(design))
       call direction % set_real_vector(e)
       call steps % partial_action(instants, [knobs], &
            & [variation(steps % argument(1), direction)], out)
       call out % real_vector(column)
       v(:, j) = column
    end do

  end subroutine step_partials

  !===================================================================!
  ! The mixed second partial of every step in two entries of a grid's
  ! design, from the grid's own partial action along both.
  !===================================================================!

  subroutine step_second_partials(steps, n, design, j, k, u)

    class(grid), intent(in) :: steps
    integer    , intent(in) :: n, j, k
    real(dp)   , intent(in) :: design(:)
    real(dp), allocatable, intent(out) :: u(:)

    type(stored_directed_graph) :: instants
    type(stored_field) :: knobs, first, second
    class(field), allocatable :: out
    real(dp), allocatable :: e(:)

    instants = stored_directed_graph(n, tails=[integer ::], heads=[integer ::])
    knobs    = stored_field('design', instants % vertex_set(), size(design))
    call knobs % set_real_vector(design)

    allocate(e(size(design)))
    e = 0.0_dp
    e(j) = 1.0_dp
    first = stored_field('direction', instants % vertex_set(), size(design))
    call first % set_real_vector(e)
    e = 0.0_dp
    e(k) = 1.0_dp
    second = stored_field('direction', instants % vertex_set(), size(design))
    call second % set_real_vector(e)
    call steps % partial_action(instants, [knobs], &
         & [variation(steps % argument(1), first), variation(steps % argument(1), second)], out)
    call out % real_vector(u)

  end subroutine step_second_partials

  !===================================================================!
  ! Where a component lies: instants follow one another, nodes lie
  ! within an instant, and the components of one point stay together,
  !
  !      ((instant - 1) nodes + (node - 1)) degrees + degree + 1
  !
  ! which at one node is (instant - 1) degrees + degree + 1, the
  ! ordinary block. A field over a mesh is this at nodes > 1 and
  ! nothing else.
  !===================================================================!

  pure integer function unknown(instant, degree, degrees, node, nodes) result(at)

    integer, intent(in)           :: instant, degree, degrees
    integer, intent(in), optional :: node, nodes

    integer :: i, m

    i = 1
    m = 1
    if (present(node))  i = node
    if (present(nodes)) m = nodes

    at = ((instant - 1) * m + (i - 1)) * degrees + degree + 1

  end function unknown

  function unknowns_graph(n, degrees) result(g)

    integer, intent(in) :: n, degrees
    type(stored_directed_graph) :: g

    g = stored_directed_graph(n * degrees, tails=[integer ::], heads=[integer ::])

  end function unknowns_graph

  !===================================================================!
  ! The derived rows of a block, as a stencil: the rows that fit, the
  ! weights on them, and the sign convention operation_scheme_stencil
  ! owns.
  !===================================================================!

  function scheme_rows(scheme, degrees, n, dt, nodes, along, along2) result(rows)

    class(family), intent(in)           :: scheme
    integer      , intent(in)           :: degrees, n
    real(dp)     , intent(in)           :: dt(:)
    integer      , intent(in), optional :: nodes
    real(dp)     , intent(in), optional :: along(:), along2(:)
    type(stencil) :: rows

    integer , allocatable :: tails(:), heads(:), source_degree(:), determines(:)
    real(dp), allocatable :: w(:)
    integer :: e, i, m

    m = 1
    if (present(nodes)) m = nodes

    call block_reach(scheme, degrees, n, tails, heads, source_degree, determines)

    ! the family's rows once per node: each node's history is its own.
    ! Along a direction in the steps the rows are the partial of the
    ! weights, which the determined component, entering with one,
    ! takes no part in.
    if (present(along)) then
       call weights_varied(scheme_weight(scheme), n, tails, heads, dt, along, &
            & source_degree, determines, w, along2)
       rows = stencil( &
            & [((unknown(heads(e), determines(e), degrees, i, m), e = 1, size(heads)), i = 1, m)], &
            & [((unknown(tails(e), source_degree(e), degrees, i, m), e = 1, size(tails)), i = 1, m)], &
            & [(-w, i = 1, m)], spread(0.0_dp, 1, n * m * degrees), 'varied rows')
    else
       call weights_of(scheme_weight(scheme), n, tails, heads, dt, source_degree, determines, w)
       rows = derived_constraints( &
            & [((unknown(heads(e), determines(e), degrees, i, m), e = 1, size(heads)), i = 1, m)], &
            & [((unknown(tails(e), source_degree(e), degrees, i, m), e = 1, size(tails)), i = 1, m)], &
            & [(w, i = 1, m)], n * m * degrees, 'derived rows')
    end if

  end function scheme_rows

  !===================================================================!
  ! The whole statement of one block. The instants the family reaches
  ! back over are carried, and the values given for them are what
  ! their rows hold.
  !===================================================================!

  function block_of(scheme, physics, degrees, n, dt, held, nodes, spatial) result(rows)

    class(family)         , intent(in)           :: scheme
    class(nodal_integrand), intent(in)           :: physics
    integer               , intent(in)           :: degrees, n
    real(dp)              , intent(in)           :: dt(:), held(:)
    integer               , intent(in), optional :: nodes
    type(stencil)         , intent(in), optional :: spatial
    type(block_residual) :: rows

    integer, allocatable :: carried(:), at(:)
    integer :: h, k, i, d, m

    m = 1
    if (present(nodes)) m = nodes

    h       = scheme % history_depth(degrees - 1)
    carried = [(((unknown(k, d, degrees, i, m), d = 0, degrees - 1), i = 1, m), k = 1, h)]
    at      = [((unknown(k, 0, degrees, i, m) - 1, i = 1, m), k = 1, n)]

    if (size(held) /= size(carried)) then
       error stop 'gti_march: one value per carried component'
    end if

    rows = block_residual(scheme_rows(scheme, degrees, n, dt, m), physics, at, &
         & n * m * degrees, degrees, scheme % primary_degree(degrees - 1), carried, held)

    ! where every unknown lies: instant k, node i, and instant k is its
    ! moment; the level below, a stencil over the nodes, is laid on
    ! every instant
    call rows % placed_in([(((k, d = 0, degrees - 1), i = 1, m), k = 1, n)], &
         &                [(((i, d = 0, degrees - 1), i = 1, m), k = 1, n)], &
         &                [(((k, d = 0, degrees - 1), i = 1, m), k = 1, n)])
    if (present(spatial)) call rows % spatial_laid(spatial)

  end function block_of

  !===================================================================!
  ! Newton over the whole block. The design is held while the state
  ! varies, which is what a minimizer supplies as an extra input.
  !
  !             WHAT COUNTS AS SOLVED
  !
  ! A scheme's rows carry a power of the step, so a difference on the
  ! second derivative weighs its sources by the inverse square of it.
  ! Refining the grid therefore raises the size of a residual for the
  ! same trajectory, and the smallest one reachable in the arithmetic
  ! rises with it: at a hundredth of a unit it is near ten to the
  ! minus thirteen, and finer than that it passes any fixed target.
  !
  ! Asked for a fixed one, newton reaches the trajectory in two steps
  ! and then spends its whole budget failing to better it. Measured on a
  ! degree-two problem over three units: a hundred and twenty instants
  ! took a hundred and sixty seconds to produce what forty iterations
  ! produce in a sixth of one, to the same six digits.
  !
  ! So the target is set against the residual the first guess gives,
  ! which is the only scale in the problem that is known before it is
  ! solved, and the budget is a backstop rather than a cost.
  !===================================================================!

  subroutine solved(rows, design_value, q, achieved, left, seed)

    type(block_residual), intent(in)  :: rows
    real(dp)            , intent(in)  :: design_value
    real(dp), allocatable, intent(out) :: q(:)
    real(dp)            , intent(out) :: achieved
    type(imbalance), intent(out), optional :: left
    real(dp)       , intent(in) , optional :: seed(:)

    type(newton) :: solver
    type(stored_directed_graph) :: unknowns
    type(stored_field) :: design
    integer :: count, width

    count    = rows % num_unknowns()
    unknowns = stored_directed_graph(count, tails=[integer ::], heads=[integer ::])
    design   = stored_field('nu', unknowns % vertex_set(), rows % num_points())
    call design % set_real_vector(spread(design_value, 1, rows % num_points()))

    ! every unknown lies in a point of degrees consecutive components,
    ! a stage's as much as an instant's, and a point is smoothed whole
    width = rows % num_degrees()
    ! multigrid coarsens by aggregates read off the block: the coarse
    ! cell of each unknown's node, at its own moment and degree
    if (multigrid_on()) call set_aggregates(rows % aggregates(coarse_nodes(rows % num_nodes())))
    call take_inner(solver % inner, count, width)
    call solver % attach(rows, unknowns, unknowns % vertex_set(), count, &
         & held_inputs = [design])

    if (present(seed)) then
       q = seed
    else
       q = at_first_instant(rows, count)
    end if

    solver % compiled       = assembly_present()
    solver % max_iterations = stopping_iterations
    solver % tolerance      = stopping_tolerance
    solver % criterion      = stopping_criterion
    solver % budget         = stopping_budget

    call solver % solve(spread(0.0_dp, 1, count), q, achieved)
    call keep_inner(solver % inner)

    if (present(left)) then
       left % converged = solver % converged(achieved)
       left % diverging = solver % diverging(achieved)
       left % norm      = achieved
       left % began     = solver % began()
       if (.not. left % converged) call by_aspect(rows, unknowns, q, design, left)
    end if

  end subroutine solved

  !===================================================================!
  ! A stamp no statement has had before.
  !===================================================================!

  integer function fresh_stamp() result(mark)

    stamps_given = stamps_given + 1
    mark = stamps_given

  end function fresh_stamp

  !===================================================================!
  ! A LINEAR SYSTEM IN THE TANGENT, A w = rhs or A^T w = rhs, solved
  ! as the block it came from is solved: as a linear block through the
  ! sweep, where newton stops after one step. The stamp given is the
  ! tangent's; every right side against the same tangent gives the
  ! same stamp, and a direct solver then factorises once.
  !===================================================================!

  subroutine solved_linear(rows, unknowns, inputs, rhs, transposed, mark, w)

    type(block_residual)       , intent(in)  :: rows
    class(directed_graph)      , intent(in)  :: unknowns
    type(stored_field)         , intent(in)  :: inputs(:)
    real(dp)                   , intent(in)  :: rhs(:)
    logical                    , intent(in)  :: transposed
    integer                    , intent(in)  :: mark
    real(dp), allocatable      , intent(out) :: w(:)

    type(block_residual) :: lin
    real(dp) :: achieved

    if (transposed) then
       call tally_record(adjoint_loops)
    else
       call tally_record(tangent_loops)
    end if

    lin = rows % linear_block(unknowns, inputs, rhs, transposed, mark)
    call swept(lin, 0.0_dp, w, achieved)

  end subroutine solved_linear

  !===================================================================!
  ! The gradient in the design by the tangent - one solve in the
  ! state, the gradient read along it - and by the adjoint - one solve
  ! against the transpose, the design partial read along it. Both
  ! through the sweep, against one tangent, stamped once.
  !===================================================================!

  real(dp) function by_tangent(rows, unknowns, inputs, g, design_rate, explicit, mark) &
       & result(df)

    type(block_residual) , intent(in) :: rows
    class(directed_graph), intent(in) :: unknowns
    type(stored_field)   , intent(in) :: inputs(:)
    real(dp)             , intent(in) :: g(:), design_rate(:), explicit
    integer              , intent(in) :: mark

    real(dp), allocatable :: w(:)

    call solved_linear(rows, unknowns, inputs, -design_rate, .false., mark, w)
    df = explicit + dot_product(g, w)

  end function by_tangent

  real(dp) function by_adjoint(rows, unknowns, inputs, g, design_rate, explicit, mark) &
       & result(df)

    type(block_residual) , intent(in) :: rows
    class(directed_graph), intent(in) :: unknowns
    type(stored_field)   , intent(in) :: inputs(:)
    real(dp)             , intent(in) :: g(:), design_rate(:), explicit
    integer              , intent(in) :: mark

    real(dp), allocatable :: lambda(:)

    call solved_linear(rows, unknowns, inputs, g, .true., mark, lambda)
    df = explicit - dot_product(lambda, design_rate)

  end function by_adjoint

  !===================================================================!
  ! Which level the marches sweep. A name that is none of the three
  ! stops the program.
  !===================================================================!

  subroutine set_sweep(name)

    character(len=*), intent(in) :: name

    call refuse_unknown(name, ['space-time', 'time      ', 'space     '], 'sweep')
    sweep_level = name

  end subroutine set_sweep

  pure function sweep_named() result(name)

    character(len=:), allocatable :: name

    name = trim(sweep_level)

  end function sweep_named

  !===================================================================!
  ! THE SWEEP. The block's points lie instant by instant, nodes within
  ! an instant, so a member of the time level is one instant's points
  ! and a member of the space level is one node's points across the
  ! instants. Each member is solved as a block of its own, restricted
  ! from the whole with the rest held at the current state, and its
  ! solution written back; a member with nothing to solve - every
  ! component carried - is passed over. A pass is judged on the whole
  ! block's residual by the same criteria as any iteration, so the
  ! time sweep, exact after one pass, stops on the second, and the
  ! space sweep stops where the coupling has settled.
  !===================================================================!

  subroutine swept(rows, design_value, q, achieved, left)

    type(block_residual), intent(in)  :: rows
    real(dp)            , intent(in)  :: design_value
    real(dp), allocatable, intent(out) :: q(:)
    real(dp)            , intent(out) :: achieved
    type(imbalance), intent(out), optional :: left

    type(block_residual) :: sub
    type(newton) :: judge
    type(stored_directed_graph) :: unknowns
    type(stored_field) :: design
    integer , allocatable :: member(:), order(:), label(:)
    real(dp), allocatable :: piece(:)
    logical , allocatable :: is_carried(:)
    real(dp) :: sub_achieved, before
    integer :: count, npts, members, m, mm, pass, k

    if (trim(sweep_level) == 'space-time') then
       call solved(rows, design_value, q, achieved, left)
       return
    end if

    count   = rows % num_unknowns()
    npts    = rows % num_points()

    ! the level's members, read off the block's own labels: instants
    ! or steps for the time level, nodes for the space level
    if (trim(sweep_level) == 'time') then
       label = rows % slice_of()
    else
       label = rows % node_of()
    end if
    members = maxval(label)



    allocate(is_carried(count), source=.false.)
    is_carried(rows % carried_unknowns()) = .true.

    ! the seed, and the carried components at what they are held at:
    ! a member that is all carried is then already solved
    q = at_first_instant(rows, count)
    q(rows % carried_unknowns()) = rows % held_values()

    unknowns = stored_directed_graph(count, tails=[integer ::], heads=[integer ::])
    design   = stored_field('nu', unknowns % vertex_set(), npts)
    call design % set_real_vector(spread(design_value, 1, npts))

    judge % max_iterations = stopping_iterations
    judge % tolerance      = stopping_tolerance
    judge % criterion      = stopping_criterion
    judge % budget         = stopping_budget
    call judge % begin_imbalance()


    ! the order the members are swept in is the coupling's own: a
    ! transposed statement, upper triangular in time, sweeps from the
    ! last instant because its pattern says so
    call rows % member_order(trim(sweep_level) == 'time', order)


    ! The residual where the sweep begins is what a relative target
    ! is measured against, as the first residual is for any march.
    achieved = whole_residual(rows, unknowns, design, q)
    call judge % note_imbalance(achieved)

    do pass = 1, stopping_iterations

       before = achieved

       do mm = 1, members
          m = order(mm)
          member = pack([(k, k = 1, count)], label == m)
          if (all(is_carried(member))) cycle

          ! a member not yet solved is seeded from the one before it in
          ! the order swept, which is continuation, the seed every step
          ! of a march has: a member of the same extent is copied, and
          ! a step's stages and arriving instant each take the instant
          ! before them
          if (pass == 1 .and. trim(sweep_level) == 'time' .and. mm > 1) then
             call continued(q, member, pack([(k, k = 1, count)], label == order(mm - 1)))
          end if

          sub = rows % restricted(member, q)
          if (rows % stamp() /= 0) then
             call sub % stamped(abs(rows % stamp()) * members + m, rows % stamp_transposed())
          end if
          call solved(sub, design_value, piece, sub_achieved, seed=q(member))
          q(member) = piece
       end do

       achieved = whole_residual(rows, unknowns, design, q)
       call judge % note_imbalance(achieved)

       if (judge % converged(achieved)) exit
       if (judge % exhausted(pass)) exit

       ! A pass that left the residual exactly where it was has
       ! reached the sweep's fixed point; another would do the same.
       if (achieved == before) exit

    end do


    if (present(left)) then
       left % converged = judge % converged(achieved)
       left % diverging = judge % diverging(achieved)
       left % norm      = achieved
       left % began     = judge % began()
       if (.not. left % converged) call by_aspect(rows, unknowns, q, design, left)
    end if

  end subroutine swept

  real(dp) function whole_residual(rows, unknowns, design, q) result(norm)

    type(block_residual)       , intent(in) :: rows
    type(stored_directed_graph), intent(in) :: unknowns
    type(stored_field)         , intent(in) :: design
    real(dp)                   , intent(in) :: q(:)

    type(stored_field) :: state
    class(field), allocatable :: out
    real(dp), allocatable :: r(:)

    state = stored_field('state', unknowns % vertex_set(), size(q))
    call state % set_real_vector(q)
    call rows % apply(unknowns, [state, design], out)
    call out % real_vector(r)
    norm = norm2(r)

  end function whole_residual

  !-------------------------------------------------------------------!
  ! The unknowns of the m-th member: the points of instant m, or the
  ! points of node m across the instants, each point's components.
  !-------------------------------------------------------------------!

  !-------------------------------------------------------------------!
  ! The seed of a member from the member solved before it. Of the
  ! same extent, the values are copied; otherwise the last point of
  ! the earlier member - the instant a step arrives at - is laid on
  ! every point of the later one, which is where a step's stages and
  ! its own arriving instant begin.
  !-------------------------------------------------------------------!

  subroutine continued(q, member, earlier)

    real(dp), intent(inout) :: q(:)
    integer , intent(in)    :: member(:), earlier(:)

    integer :: pieces, i, width

    ! a member of the same extent is copied; a step's stages and
    ! arriving instant each take the instant before them; an instant
    ! after a step takes the step's last piece
    if (size(member) == size(earlier)) then
       q(member) = q(earlier)
    else if (mod(size(member), size(earlier)) == 0) then
       width  = size(earlier)
       pieces = size(member) / width
       do i = 1, pieces
          q(member((i - 1) * width + 1:i * width)) = q(earlier)
       end do
    else if (mod(size(earlier), size(member)) == 0) then
       width = size(member)
       q(member) = q(earlier(size(earlier) - width + 1:))
    else
       error stop 'gti_march: a member is seeded from one of its extent, a multiple of it, or a divisor'
    end if

  end subroutine continued

  !===================================================================!
  ! The aspects of what was left: the norm split by degree, the
  ! largest entry, and the largest entry of A^T r / ||r||. A is formed
  ! here in full, which is O(n^2) and is paid only on a march that
  ! did not converge.
  !===================================================================!

  subroutine by_aspect(rows, unknowns, q, design, left)

    type(block_residual)       , intent(in)    :: rows
    type(stored_directed_graph), intent(in)    :: unknowns
    real(dp)                   , intent(in)    :: q(:)
    type(stored_field)         , intent(in)    :: design
    type(imbalance)            , intent(inout) :: left

    type(stored_field) :: state
    class(field), allocatable :: out
    real(dp), allocatable :: r(:), a(:,:), slope(:)
    integer :: i, d, nd, n

    n  = size(q)
    nd = rows % num_degrees()

    state = stored_field('state', unknowns % vertex_set(), n)
    call state % set_real_vector(q)
    call rows % apply(unknowns, [state, design], out)
    call out % real_vector(r)

    allocate(left % by_degree(0:nd - 1), source=0.0_dp)
    do i = 1, n
       d = mod(i - 1, nd)
       left % by_degree(d) = left % by_degree(d) + r(i) ** 2
    end do
    left % by_degree = sqrt(left % by_degree)

    i = maxloc(abs(r), dim=1)
    left % worst_slot   = (i - 1) / nd + 1
    left % worst_degree = mod(i - 1, nd)

    if (left % norm <= 0.0_dp) return

    call jacobian_of(rows, unknowns, [state, design], n, unknowns % vertex_set(), a)
    slope = matmul(r, a) / left % norm

    i = maxloc(abs(slope), dim=1)
    left % steepest_slot   = (i - 1) / nd + 1
    left % steepest_degree = mod(i - 1, nd)
    left % steepest        = slope(i)

  end subroutine by_aspect

  !===================================================================!
  ! How large a residual the first guess gives, which is the scale
  ! the target is set against.
  !===================================================================!


  !===================================================================!
  ! A first guess: every point of the block holding what its first
  ! instant was given. It costs nothing to form and it starts newton
  ! near the trajectory rather than at zero, which for a state of any
  ! size is far away and is where a jacobian is most likely to be
  ! singular.
  !===================================================================!

  function at_first_instant(rows, count) result(q)

    type(block_residual), intent(in) :: rows
    integer             , intent(in) :: count
    real(dp), allocatable :: q(:)

    real(dp), allocatable :: one(:)
    integer , allocatable :: at(:)
    integer :: p, nd

    one = rows % first_held()
    nd  = size(one)
    at  = rows % points_at()

    allocate(q(count), source=0.0_dp)

    do p = 1, size(at)
       q(at(p) + 1:at(p) + nd) = one
    end do

  end function at_first_instant

  !===================================================================!
  ! The linear solver inside newton. A dense factorisation forms the
  ! jacobian column by column - one application of the statement per
  ! unknown - and then costs the cube of the count to factor, so it
  ! wins while the count is small and loses badly once it is not. The
  ! statement supplies a matvec through its partial action, so a
  ! krylov solver forms no matrix at all.
  !
  ! Where the crossing sits, and why it is settable, is stated in
  ! gti_sweeps, which owns it.
  !===================================================================!


  !===================================================================!
  ! Where each block begins and ends. A block adds the instants given
  ! for it and reaches back over its predecessor's last, so the
  ! blocks overlap by exactly what each family reaches.
  !===================================================================!

  subroutine horizon_bounds(schemes, added, equation_degree, first, last)

    type(family_holder), intent(in) :: schemes(:)
    integer            , intent(in) :: added(:), equation_degree
    integer, allocatable, intent(out) :: first(:), last(:)

    integer :: b

    if (size(schemes) /= size(added)) then
       error stop 'gti_march: one family and one instant count per block'
    end if

    allocate(first(size(added)), last(size(added)))

    do b = 1, size(added)
       if (b == 1) then
          first(b) = 1
          last(b)  = added(b)
       else
          first(b) = last(b - 1) - schemes(b) % scheme % history_depth(equation_degree) + 1
          last(b)  = last(b - 1) + added(b)
       end if

       if (added(b) <= schemes(b) % scheme % history_depth(equation_degree)) then
          error stop 'gti_march: a block adds more instants than its family reaches'
       end if
       if (first(b) < 1) then
          error stop 'gti_march: the horizon holds every instant its blocks reach back over'
       end if
    end do

  end subroutine horizon_bounds

end module gti_march
