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
  use operation_coupling      , only : weights_of
  use gti_configuration       , only : refuse_unknown
  use operation_weight        , only : scheme_weight
  use view_directed_stored    , only : stored_directed_graph
  use view_directed           , only : directed_graph
  use field_calculus          , only : field
  use field_stored            , only : stored_field
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
       & aggregates_given, aggregates_of, set_aggregates, take_inner, keep_inner, forget_inner
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
  public :: partition, partitioned, scheme_rows, block_of, solved, unknowns_graph
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

    type(block_residual) :: rows
    type(stencil) :: none
    real(dp) :: achieved
    integer :: k

    if (size(lower) /= degrees - 1) then
       error stop 'gti_march: the components below the highest are given, and no others'
    end if

    none = stencil([integer ::], [integer ::], [real(dp) ::], &
         & spread(0.0_dp, 1, degrees), 'none')

    rows = block_residual(none, physics, at=[0], unknowns=degrees, degrees=degrees, &
         & primary=degrees - 1, carried=[(k, k = 1, degrees - 1)], held=lower)

    call solved(rows, design_value, q, achieved)

    if (.not. achieved <= stopping_tolerance * max(1.0_dp, norm2(lower))) then
       write(*,'(a,es12.3)') ' the physics at the initial instant left a residual of ', achieved
       error stop 'gti_march: the initial state is consistent with the physics'
    end if

  end function consistent_state

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

  pure integer function unknown(instant, degree, degrees) result(at)

    integer, intent(in) :: instant, degree, degrees

    at = (instant - 1) * degrees + degree + 1

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

  function scheme_rows(scheme, degrees, n, dt) result(rows)

    class(family), intent(in) :: scheme
    integer      , intent(in) :: degrees, n
    real(dp)     , intent(in) :: dt(:)
    type(stencil) :: rows

    integer , allocatable :: tails(:), heads(:), source_degree(:), determines(:)
    real(dp), allocatable :: w(:)
    integer :: e

    call block_reach(scheme, degrees, n, tails, heads, source_degree, determines)
    call weights_of(scheme_weight(scheme), n, tails, heads, dt, source_degree, determines, w)

    rows = derived_constraints( &
         & [(unknown(heads(e), determines(e), degrees), e = 1, size(heads))], &
         & [(unknown(tails(e), source_degree(e), degrees), e = 1, size(tails))], &
         & w, n * degrees, 'derived rows')

  end function scheme_rows

  !===================================================================!
  ! The whole statement of one block. The instants the family reaches
  ! back over are carried, and the values given for them are what
  ! their rows hold.
  !===================================================================!

  function block_of(scheme, physics, degrees, n, dt, held) result(rows)

    class(family)         , intent(in) :: scheme
    class(nodal_integrand), intent(in) :: physics
    integer               , intent(in) :: degrees, n
    real(dp)              , intent(in) :: dt(:), held(:)
    type(block_residual) :: rows

    integer, allocatable :: carried(:)
    integer :: h, k, d

    h = scheme % history_depth(degrees - 1)
    carried = [((unknown(k, d, degrees), d = 0, degrees - 1), k = 1, h)]

    if (size(held) /= size(carried)) then
       error stop 'gti_march: one value per carried component'
    end if

    rows = block_residual(scheme_rows(scheme, degrees, n, dt), physics, &
         & [((k - 1) * degrees, k = 1, n)], n * degrees, degrees, &
         & scheme % primary_degree(degrees - 1), carried, held)

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
    integer :: count

    count    = rows % num_unknowns()
    unknowns = stored_directed_graph(count, tails=[integer ::], heads=[integer ::])
    design   = stored_field('nu', unknowns % vertex_set(), rows % num_points())
    call design % set_real_vector(spread(design_value, 1, rows % num_points()))

    call take_inner(solver % inner, count)
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

  subroutine solved_linear(rows, unknowns, inputs, rhs, transposed, mark, nodes, w)

    type(block_residual)       , intent(in)  :: rows
    class(directed_graph)      , intent(in)  :: unknowns
    type(stored_field)         , intent(in)  :: inputs(:)
    real(dp)                   , intent(in)  :: rhs(:)
    logical                    , intent(in)  :: transposed
    integer                    , intent(in)  :: mark, nodes
    real(dp), allocatable      , intent(out) :: w(:)

    type(block_residual) :: lin
    real(dp) :: achieved

    if (transposed) then
       call tally_record(adjoint_loops)
    else
       call tally_record(tangent_loops)
    end if

    lin = rows % linear_block(unknowns, inputs, rhs, transposed, mark)
    call swept(lin, 0.0_dp, nodes, w, achieved, backward=transposed)

  end subroutine solved_linear

  !===================================================================!
  ! The gradient in the design by the tangent - one solve in the
  ! state, the gradient read along it - and by the adjoint - one solve
  ! against the transpose, the design partial read along it. Both
  ! through the sweep, against one tangent, stamped once.
  !===================================================================!

  real(dp) function by_tangent(rows, unknowns, inputs, g, design_rate, explicit, nodes, &
       & mark) result(df)

    type(block_residual) , intent(in) :: rows
    class(directed_graph), intent(in) :: unknowns
    type(stored_field)   , intent(in) :: inputs(:)
    real(dp)             , intent(in) :: g(:), design_rate(:), explicit
    integer              , intent(in) :: nodes, mark

    real(dp), allocatable :: w(:)

    call solved_linear(rows, unknowns, inputs, -design_rate, .false., mark, nodes, w)
    df = explicit + dot_product(g, w)

  end function by_tangent

  real(dp) function by_adjoint(rows, unknowns, inputs, g, design_rate, explicit, nodes, &
       & mark) result(df)

    type(block_residual) , intent(in) :: rows
    class(directed_graph), intent(in) :: unknowns
    type(stored_field)   , intent(in) :: inputs(:)
    real(dp)             , intent(in) :: g(:), design_rate(:), explicit
    integer              , intent(in) :: nodes, mark

    real(dp), allocatable :: lambda(:)

    call solved_linear(rows, unknowns, inputs, g, .true., mark, nodes, lambda)
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

  subroutine swept(rows, design_value, nodes, q, achieved, left, backward)

    type(block_residual), intent(in)  :: rows
    real(dp)            , intent(in)  :: design_value
    integer             , intent(in)  :: nodes
    real(dp), allocatable, intent(out) :: q(:)
    real(dp)            , intent(out) :: achieved
    type(imbalance), intent(out), optional :: left
    logical        , intent(in) , optional :: backward

    type(block_residual) :: sub
    type(newton) :: judge
    type(stored_directed_graph) :: unknowns
    type(stored_field) :: design
    integer , allocatable :: at(:), member(:), whole(:)
    real(dp), allocatable :: piece(:)
    logical , allocatable :: is_carried(:)
    real(dp) :: sub_achieved, before
    integer :: count, degrees, npts, instants, members, m, mm, pass, neighbour
    logical :: reversed

    if (trim(sweep_level) == 'space-time') then
       call solved(rows, design_value, q, achieved, left)
       return
    end if

    count    = rows % num_unknowns()
    degrees  = rows % num_degrees()
    at       = rows % points_at()
    npts     = size(at)
    instants = npts / nodes

    ! A stage block keeps its instants between its stages, and its
    ! points - the stages - do not tile its unknowns, so a sweep by
    ! instants is not defined on it and it is solved whole. A member
    ! of one step, its stages and the instant it arrives at, is the
    ! sweep such a block would take, and is not built.
    if (count /= npts * degrees) then
       call solved(rows, design_value, q, achieved, left)
       return
    end if

    if (instants * nodes /= npts) then
       error stop 'gti_march: the points lie instant by instant, the nodes within'
    end if

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

    if (trim(sweep_level) == 'time') then
       members = instants
    else
       members = nodes
    end if

    ! a transposed statement is upper triangular in time, and its
    ! instants are swept from the last
    reversed = .false.
    if (present(backward)) reversed = backward

    ! the block's aggregates, kept aside while the members set their own
    if (multigrid_on()) call aggregates_of(whole)

    ! The residual where the sweep begins is what a relative target
    ! is measured against, as the first residual is for any march.
    achieved = whole_residual(rows, unknowns, design, q)
    call judge % note_imbalance(achieved)

    do pass = 1, stopping_iterations

       before = achieved

       do mm = 1, members
          m = mm
          if (reversed) m = members - mm + 1
          member = member_unknowns(at, degrees, nodes, instants, m)
          if (all(is_carried(member))) cycle

          ! an instant not yet solved is seeded from the one before it
          ! in the order swept, which is continuation: the seed every
          ! step of a march has
          neighbour = m - 1
          if (reversed) neighbour = m + 1
          if (pass == 1 .and. trim(sweep_level) == 'time' .and. mm > 1) then
             q(member) = q(member_unknowns(at, degrees, nodes, instants, neighbour))
          end if

          sub = rows % restricted(member, q)
          if (rows % stamp() /= 0) then
             call sub % stamped(sign(abs(rows % stamp()) * members + m, rows % stamp()))
          end if
          if (multigrid_on()) call set_aggregates(member_aggregates(member, whole))
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

    if (multigrid_on()) call set_aggregates(whole)

    if (present(left)) then
       left % converged = judge % converged(achieved)
       left % diverging = judge % diverging(achieved)
       left % norm      = achieved
       left % began     = judge % began()
       if (.not. left % converged) call by_aspect(rows, unknowns, q, design, left)
    end if

  end subroutine swept

  !-------------------------------------------------------------------!
  ! The aggregates of a member, renumbered from one in the order they
  ! first appear, so multigrid on the member coarsens as the whole
  ! block would.
  !-------------------------------------------------------------------!

  function member_aggregates(member, whole) result(agg)

    integer, intent(in) :: member(:)
    integer, intent(in), allocatable :: whole(:)
    integer, allocatable :: agg(:)

    integer, allocatable :: renumbered(:)
    integer :: i, next

    if (.not. allocated(whole)) then
       error stop 'gti_march: multigrid coarsens by aggregates, and none were given'
    end if

    allocate(agg(size(member)), renumbered(maxval(whole)), source=0)
    next = 0
    do i = 1, size(member)
       if (renumbered(whole(member(i))) == 0) then
          next = next + 1
          renumbered(whole(member(i))) = next
       end if
       agg(i) = renumbered(whole(member(i)))
    end do

  end function member_aggregates

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

  function member_unknowns(at, degrees, nodes, instants, m) result(member)

    integer, intent(in) :: at(:), degrees, nodes, instants, m
    integer, allocatable :: member(:)

    integer :: k, i, d, p, e

    if (trim(sweep_level) == 'time') then
       allocate(member(nodes * degrees))
       e = 0
       do i = 1, nodes
          p = (m - 1) * nodes + i
          do d = 1, degrees
             e = e + 1
             member(e) = at(p) + d
          end do
       end do
    else
       allocate(member(instants * degrees))
       e = 0
       do k = 1, instants
          p = (k - 1) * nodes + m
          do d = 1, degrees
             e = e + 1
             member(e) = at(p) + d
          end do
       end do
    end if

  end function member_unknowns

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
