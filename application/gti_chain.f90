!=====================================================================!
! A horizon of blocks whose layouts need not agree.
!
! A multistep block holds one set of components per instant. A stage
! block holds, for each step, its stages and then the instant it
! arrives at. So a horizon that changes from one family to the other
! cannot keep its state in one array indexed the same way throughout,
! and the junction between two blocks stops being a contiguous copy.
!
! What it becomes is an index map, and the only thing it needs is a
! question each block already resolves for itself: where among its
! unknowns its k-th instant sits. A block hands its successor
! components, never rows, so that question is the whole of the
! interface between them.
!
! A block may reach back further than the block before it is long - a
! stage family spans one instant and a backward difference of order
! three on a degree-three equation looks back over nine - so what a
! block is given is gathered from whichever earlier block computed
! each instant, not from the one immediately before it.
!
!             WHOSE INSTANT IS IT
!
! An instant shared by two blocks is computed by the earlier and
! carried by the later, so the functional counts it once, under the
! block that computed it. The first block additionally owns the
! instants it was given, whose values are initial conditions: they
! contribute to the functional and not to any derivative of it,
! because they do not move.
!
!             THE EXPANSION ALONG A CHAIN
!
! Every order travels the junction the way the trajectory does. At
! order m each block solves against its own jacobian for a right side
! its own physics determines, with its carried rows set to what its
! predecessor found at those instants at that same order. Order zero
! is the march itself.
!
!             THE ORDER A COEFFICIENT SITS AT
!
! A series is indexed from zero, the order being the index, and every
! array that holds one is allocated that way on purpose. A section of
! such an array is indexed from one, so copying one into a fresh
! array and then striking out an order by its number strikes out the
! order below it. Where that happens the trajectory itself is wiped
! and every derivative comes out exactly zero, which is what it did.
!
!             WHAT IS REFUSED
!
! A chain of no blocks; a block that adds no more instants than its
! family reaches back over; and initial conditions that are not one
! value per degree over the instants the first block was given.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_chain

  use util_precision  , only : dp
  use operation_family , only : family
  use operation_grid   , only : grid
  use operation_expression, only : expression
  use gti_expansion    , only : family_holder, marches_by_stages, expansion, &
       & design_of_physics, design_of_steps
  use gti_block        , only : block_residual
  use gti_march        , only : imbalance, swept, solved_linear, fresh_stamp, partitioned, horizon_bounds, solved, &
       & frozen_inputs
  use gti_march        , only : block_from
  use gti_sweeps       , only : functional_design_partial, varied_by
  use operation_stencil, only : stencil
  use operation_family_dirk, only : crouzeix_three_stage
  use view_directed_stored, only : stored_directed_graph
  use field_calculus   , only : field
  use field_stored     , only : stored_field
  use operation_action , only : variation
  use gti_sweeps       , only : design_partial, route_of, functional_gradient, &
       & route_substitutions, forward_route, reverse_route
  use util_tally            , only : tally_order, tally_enter, tally_leave, &
       & at_horizon, at_block, at_stage
  use gti_taylor       , only : nodal_coefficient, order_of_series

  implicit none

  private
  public :: chain_block, march_chain, chain_expansion, instant_components
  public :: functional_holder, one_functional, first_of, chain_hessian
  public :: chain_system, chain_systems, chain_by_tangent, chain_by_adjoint
  public :: expansion_substitutions


  !===================================================================!
  ! One block of a chain: its statement, where its instants sit among
  ! its unknowns, which global instants it spans, how many it was
  ! given, and what it computed.
  !===================================================================!

  type :: chain_block
     type(block_residual)  :: rows
     integer , allocatable :: instants_at(:)
     real(dp), allocatable :: state(:)
     ! WHERE THE BLOCK LIES on the horizon: its first and last instant
     ! as fine indices, and its stride, the fine indices between its
     ! own instants - one for a startup block over refined steps, the
     ! refinement for every block over the march's own steps, and one
     ! for both when there is no startup
     integer               :: first  = 0
     integer               :: last   = 0
     integer               :: stride = 1
     integer               :: given = 0
     integer               :: primary = 0
     ! the components one instant holds: the degrees at every node
     integer :: width = 0
     integer :: nodes = 1
     ! the family and the steps the block was built from, so that its
     ! rows can be differentiated along a direction in the steps; each
     ! step as a fraction of the march's step it lies in, and which
     ! one, so that a direction in the march's steps reads on its own
     integer , allocatable :: coarse_step(:)
     class(family), allocatable :: scheme
     real(dp)     , allocatable :: dt(:)
     real(dp)              :: fraction = 1.0_dp
     ! whether its instants count in the functional: a startup's do
     ! not, its instants being what the first block is given
     logical               :: counted = .true.
     ! the imbalance its solve began at, against which its relative
     ! tolerance was measured
     real(dp) :: began = 0.0_dp
  end type chain_block

  !===================================================================!
  ! What one block contributes to a sensitivity: its jacobian in the
  ! state, its partial in the design, and the part of the
  ! functional's gradient it owns.
  !===================================================================!

  !===================================================================!
  ! One functional, held so that several can be handed over at once.
  !===================================================================!

  type :: functional_holder
     type(expression) :: rule
  end type functional_holder

  type :: chain_system
     ! the block's partial in every design, one column per design:
     ! the physics' first, then one per entry of the grid's
     real(dp), allocatable :: rate(:,:)
     ! the part of every functional's gradient the block owns, one
     ! column per functional
     real(dp), allocatable :: g(:,:)
     ! what each functional itself adds through the designs, one
     ! row per functional and one column per design: the integrand's
     ! partial in the physics' design, the measure's in the grid's
     real(dp), allocatable :: explicit(:,:)

     ! The stamp of the tangent at the frozen state. Every order of
     ! the expansion, the tangent and the adjoint solve against that
     ! one tangent, and a direct solver factorises it once.
     integer :: mark = 0

  end type chain_system

contains

  !===================================================================!
  ! The block holding a fine instant, the latest that does - a block
  ! recomputes the instants it was given, so the latest is what the
  ! next block reads - and where the instant lies in it. None holds
  ! it, and the block is zero.
  !===================================================================!

  pure subroutine locate(chain, fine, held_by, local)

    type(chain_block), intent(in)  :: chain(:)
    integer          , intent(in)  :: fine
    integer          , intent(out) :: held_by, local

    integer :: b

    held_by = 0
    local   = 0
    do b = size(chain), 1, -1
       if (fine < chain(b) % first .or. fine > chain(b) % last) cycle
       if (mod(fine - chain(b) % first, chain(b) % stride) /= 0) cycle
       held_by = b
       local   = (fine - chain(b) % first) / chain(b) % stride + 1
       return
    end do

  end subroutine locate

  !===================================================================!
  ! The components a chain holds at one fine instant, and at one of
  ! the march's own instants, which is the fine one at the march's
  ! stride.
  !===================================================================!

  pure function fine_components(chain, fine) result(x)

    type(chain_block), intent(in) :: chain(:)
    integer          , intent(in) :: fine
    real(dp), allocatable :: x(:)

    integer :: b, local, at

    call locate(chain, fine, b, local)
    if (b == 0) error stop 'gti_chain: that instant lies outside the chain'
    at = chain(b) % instants_at(local)
    x  = chain(b) % state(at + 1:at + chain(b) % width)

  end function fine_components

  pure function instant_components(chain, instant) result(x)

    type(chain_block), intent(in) :: chain(:)
    integer          , intent(in) :: instant
    real(dp), allocatable :: x(:)

    x = fine_components(chain, 1 + (instant - 1) * chain(size(chain)) % stride)

  end function instant_components


  !===================================================================!
  ! What a block is given at the instants it shares with the ones
  ! before it, laid out the way its own carried rows expect: instant
  ! by instant at its own stride, the components within an instant.
  !===================================================================!

  pure function handed_over(earlier, first, stride, given) result(held)

    type(chain_block), intent(in) :: earlier(:)
    integer          , intent(in) :: first, stride, given
    real(dp), allocatable :: held(:)

    integer :: i, width

    width = earlier(1) % width
    allocate(held(given * width))

    do i = 1, given
       held((i - 1) * width + 1:i * width) = fine_components(earlier, first + (i - 1) * stride)
    end do

  end function handed_over

  !===================================================================!
  ! One block's statement and where its instants lie, read from its
  ! node of the expansion graph.
  !===================================================================!

  subroutine built(tower, b, scheme, physics, held, rows, instants_at)

    type(expansion)       , intent(in), target :: tower
    integer               , intent(in)  :: b
    class(family)         , intent(in)  :: scheme
    type(expression)      , intent(in)  :: physics
    real(dp)              , intent(in)  :: held(:)
    type(block_residual)  , intent(out) :: rows
    integer, allocatable  , intent(out) :: instants_at(:)

    call block_from(tower, b, scheme, physics, held, rows, instants_at)

  end subroutine built

  !===================================================================!
  ! The whole chain built and marched, block after block, each given
  ! what its predecessor computed at the instants they share.
  !===================================================================!

  subroutine march_chain(schemes, added, physics, degrees, steps, &
       & design, initial, chain, tower, dt, t, achieved, grid_design, left, nodes, spatial_discretization_stencil, &
       & startup)

    type(family_holder)   , intent(in) :: schemes(:)
    integer               , intent(in) :: added(:), degrees
    type(expression)      , intent(in) :: physics
    real(dp)              , intent(in) :: design, initial(:)
    class(grid)           , intent(in) :: steps
    type(chain_block), allocatable, intent(out) :: chain(:)
    ! THE GRAPH the chain is read from, the caller's, built here and
    ! outliving the march: every block lies at its node of it
    type(expansion), allocatable, intent(inout), target :: tower
    real(dp)         , allocatable, intent(out) :: dt(:), t(:)
    real(dp)              , intent(out) :: achieved
    real(dp), intent(in), optional     :: grid_design(:)
    type(imbalance), intent(out), optional :: left
    integer        , intent(in) , optional :: nodes
    type(stencil)  , intent(in) , optional :: spatial_discretization_stencil
    ! given, the first block's given instants are marched first by a
    ! stage family of order four, every step split this many ways,
    ! as block one of the chain: a startup that is part of the chain
    ! and so of every derivative, and reads the initial state at the
    ! first instant alone
    integer        , intent(in) , optional :: startup

    type(family_holder), allocatable :: every(:)
    type(imbalance) :: one_left
    integer , allocatable :: first(:), last(:), spans(:)
    real(dp), allocatable :: fine(:), knobs(:)
    real(dp) :: one_achieved
    integer :: b, k, r, given, before
    logical :: with_startup

    if (size(added) < 1) then
       error stop 'gti_chain: a chain holds at least one block'
    end if

    call horizon_bounds(schemes, added, degrees - 1, first, last)
    call partitioned(steps, last(size(added)), dt, t, grid_design)

    given        = schemes(1) % scheme % history_depth(degrees - 1)
    with_startup = .false.
    r            = 1
    if (present(startup)) then
       if (given > 1) then
          with_startup = .true.
          r            = max(startup, 1)
       end if
    end if
    before = merge(1, 0, with_startup)
    allocate(chain(size(added) + before))

    ! THE GRAPH the blocks are read from: one node per block, slice
    ! and component, with the couplings' relations. The expansion lays
    ! its blocks end to end, where the chain's blocks share the
    ! instants one hands the next, so the tower is built over each
    ! block's own span with each block's own steps laid end to end -
    ! a block reads only its own steps, and the sharing is the
    ! chain's junction. The startup, over its refined steps, is the
    ! first block of the same tower.
    knobs = [real(dp) ::]
    allocate(every(size(added) + before))
    if (with_startup) then
       fine  = [0.0_dp, (dt(1 + (k - 1) / r + 1) / real(r, dp), k = 1, (given - 1) * r)]
       knobs = fine(2:)
       allocate(every(1) % scheme, source=crouzeix_three_stage())
    end if
    do b = 1, size(added)
       ! the step ending at a block's first instant: none at the
       ! horizon's first, and after a startup the two share an instant,
       ! so a positive placeholder no row reads takes the place of the zero
       ! a partition refuses
       if (b == 1 .and. with_startup) then
          knobs = [knobs, fine(size(fine)), dt(2:last(1))]
       else
          knobs = [knobs, dt(first(b) + merge(1, 0, b == 1):last(b))]
       end if
       allocate(every(before + b) % scheme, source=schemes(b) % scheme)
    end do
    allocate(spans(size(every)))
    if (with_startup) spans(1) = (given - 1) * r + 1
    do b = 1, size(added)
       spans(before + b) = last(b) - first(b) + 1
    end do
    if (allocated(tower)) deallocate(tower)
    allocate(tower)
    call tower % build(physics, every, spans, steps, 0, design, nodes, spatial_discretization_stencil, &
         & weights=grid_design, block_steps=knobs)

    achieved = 0.0_dp
    call tally_enter(at_horizon)
    if (with_startup) then
       call one_block(chain, 1, tower, 1, every(1) % scheme, physics, degrees, 1, &
            & (given - 1) * r + 1, 1, fine, [0, (1 + (k - 1) / r + 1, k = 1, (given - 1) * r)], &
            & 1.0_dp / real(r, dp), .false., design, initial, one_achieved, one_left, nodes, &
            & spatial_discretization_stencil)
       achieved = one_achieved
       if (present(left)) left = one_left
    end if
    do b = 1, size(added)
       call one_block(chain, before + b, tower, before + b, schemes(b) % scheme, physics, degrees, &
            & 1 + (first(b) - 1) * r, 1 + (last(b) - 1) * r, r, dt(first(b):last(b)), &
            & [(k, k = first(b), last(b))], 1.0_dp, .true., design, initial, &
            & one_achieved, one_left, nodes, spatial_discretization_stencil)
       achieved = max(achieved, one_achieved)
       ! The report kept is the first block's that did not converge:
       ! every block after it reads a state it never reached.
       if (present(left)) then
          if (before + b == 1) left = one_left
          if (left % converged .and. .not. one_left % converged) left = one_left
       end if
    end do
    call tally_leave()

  end subroutine march_chain

  !===================================================================!
  ! One block of a chain: given what its predecessor computed at the
  ! instants they share, or the initial conditions if it is first,
  ! then built and solved.
  !===================================================================!

  subroutine one_block(chain, b, tower, in_tower, scheme, physics, degrees, first, last, &
       & stride, dt, coarse_step, fraction, counted, design, initial, achieved, left, nodes, &
       & spatial_discretization_stencil)

    type(chain_block)     , intent(inout) :: chain(:)
    type(expansion)       , intent(in), target :: tower
    integer               , intent(in)    :: b, in_tower, degrees, first, last, stride
    integer               , intent(in)    :: coarse_step(:)
    class(family)         , intent(in)    :: scheme
    type(expression)      , intent(in)    :: physics
    real(dp)              , intent(in)    :: dt(:), fraction, design, initial(:)
    logical               , intent(in)    :: counted
    real(dp)              , intent(out)   :: achieved
    type(imbalance)       , intent(out)   :: left
    integer      , intent(in), optional   :: nodes
    type(stencil), intent(in), optional   :: spatial_discretization_stencil

    real(dp), allocatable :: held(:)

    chain(b) % first    = first
    chain(b) % last     = last
    chain(b) % stride   = stride
    chain(b) % given    = scheme % history_depth(degrees - 1)
    chain(b) % primary  = scheme % primary_degree(degrees - 1)
    chain(b) % width    = degrees
    if (present(nodes)) then
       chain(b) % width = degrees * nodes
       chain(b) % nodes = nodes
    end if
    allocate(chain(b) % scheme, source=scheme)
    chain(b) % dt          = dt
    chain(b) % coarse_step = coarse_step
    chain(b) % fraction    = fraction
    chain(b) % counted     = counted

    if (b == 1) then
       if (size(initial) /= chain(b) % given * chain(b) % width) then
          error stop 'gti_chain: the initial state holds the first block''s given instants'
       end if
       held = initial
    else
       held = handed_over(chain(1:b - 1), first, stride, chain(b) % given)
    end if

    ! A block whose scheme keeps stages within a step is filed under
    ! the stage level, every other under the block level.
    if (scheme % num_stages() > 1) then
       call tally_enter(at_stage)
    else
       call tally_enter(at_block)
    end if
    call built(tower, in_tower, scheme, physics, held, chain(b) % rows, &
         & chain(b) % instants_at)
    associate (u1 => nodes, u2 => spatial_discretization_stencil); end associate
    call swept(chain(b) % rows, design, chain(b) % state, achieved, left)
    chain(b) % began = left % began
    call tally_leave()

  end subroutine one_block

  !===================================================================!
  ! The instants one block owns, as its own local indices: the ones
  ! it computed, and the ones it was given as well when no counted
  ! block before it holds them. A startup block owns none.
  !===================================================================!

  pure subroutine owned(chain, b, from, to)

    type(chain_block), intent(in)  :: chain(:)
    integer          , intent(in)  :: b
    integer          , intent(out) :: from, to

    to = size(chain(b) % instants_at)
    if (.not. chain(b) % counted) then
       from = to + 1
       return
    end if
    from = 1
    if (b > 1) then
       if (chain(b - 1) % counted) from = 1 + chain(b) % given
    end if

  end subroutine owned

  !===================================================================!
  ! The functional and every derivative of it in the design, along a
  ! chain. Each order sweeps the blocks forward, handing its
  ! coefficient over at every junction just as the trajectory does.
  !===================================================================!

  subroutine chain_expansion(chain, tower, functionals, degrees, max_order, f, node_measure)

    type(chain_block)      , intent(in) :: chain(:)
    type(expansion)        , intent(in) :: tower
    type(functional_holder), intent(in) :: functionals(:)
    integer                , intent(in) :: degrees, max_order
    real(dp), allocatable  , intent(out) :: f(:,:)
    real(dp), intent(in), optional      :: node_measure(:)

    ! One design: the expansion is in the physics' design alone, and
    ! that is what the gate is asked about at every order.
    integer, parameter :: num_designs = 1

    type(chain_system), allocatable :: systems(:)
    real(dp), allocatable :: series(:,:,:)
    type(expression) :: physics
    real(dp) :: design
    integer :: b, m, widest, route

    physics = tower % rule()
    design  = tower % parameter()

    widest = 0
    do b = 1, size(chain)
       widest = max(widest, chain(b) % rows % num_unknowns())
    end do

    allocate(series(0:max_order, widest, size(chain)), source=0.0_dp)

    do b = 1, size(chain)
       series(0, 1:size(chain(b) % state), b) = chain(b) % state
    end do

    ! The jacobian of every block, factorised once. Every order below
    ! substitutes against it.
    if (max_order >= 1) then
       call chain_systems(chain, tower, functionals, degrees, systems, node_measure)
    end if

    do m = 1, max_order
       call tally_order(m)
       call tally_enter(at_horizon)

       ! THE GATE. The route is chosen from the counts and not assumed.
       ! Only the forward route is built for an expansion, so a choice
       ! of the other stops the program and says so rather than taking
       ! the dearer one silently.
       route = route_of(num_designs, size(functionals), m)
       if (route /= forward_route) then
          write(*,'(a,i0,a)') ' the reverse route is the cheaper at order ', m, &
               & ' and an expansion by it is not built.'
          error stop 'gti_chain: an expansion by the reverse route is not built'
       end if

       do b = 1, size(chain)
          call tally_enter(at_block)
          call one_order(chain, systems, b, physics, degrees, design, m, series)
          call tally_leave()
       end do
       call tally_leave()
    end do
    call tally_order(0)

    call chain_functional(chain, functionals, degrees, design, max_order, series, f, &
         & node_measure)

  end subroutine chain_expansion

  !===================================================================!
  ! What one block's tangent is frozen at: its own trajectory and the
  ! design, over its own unknowns.
  !===================================================================!

  subroutine frozen_at(b, design, unknowns, inputs)

    type(chain_block), intent(in) :: b
    real(dp)         , intent(in) :: design
    type(stored_directed_graph), intent(out) :: unknowns
    type(stored_field), allocatable, intent(out) :: inputs(:)

    call frozen_inputs(b % state, design, b % rows % num_points(), unknowns, inputs)

  end subroutine frozen_at

  !===================================================================!
  ! One block at one order: its own physics on the rows it governs,
  ! and its predecessor's coefficient on the rows it was given.
  !===================================================================!

  subroutine one_order(chain, systems, b, physics, degrees, design, m, series)

    type(chain_block)     , intent(in)    :: chain(:)
    type(chain_system)    , intent(in)    :: systems(:)
    integer               , intent(in)    :: b, degrees, m
    type(expression)      , intent(in)    :: physics
    real(dp)              , intent(in)    :: design
    real(dp)              , intent(inout) :: series(0:, :, :)

    type(stored_directed_graph) :: unknowns
    type(stored_field), allocatable :: inputs(:)
    real(dp), allocatable :: w(:), held(:)
    integer :: count, mark

    count = chain(b) % rows % num_unknowns()
    mark  = systems(b) % mark

    call frozen_at(chain(b), design, unknowns, inputs)

    if (b == 1) then
       call order_of_series(chain(b) % rows, physics, unknowns, inputs, degrees, &
            & chain(b) % rows % points_at(), chain(b) % primary, chain(b) % rows % num_carried(), &
            & design, m, series(:, 1:count, b), mark, w)
    else
       held = coefficients_handed(chain(1:b - 1), chain(b) % first, chain(b) % stride, &
            & chain(b) % given, series, m)
       call order_of_series(chain(b) % rows, physics, unknowns, inputs, degrees, &
            & chain(b) % rows % points_at(), chain(b) % primary, chain(b) % rows % num_carried(), &
            & design, m, series(:, 1:count, b), mark, w, handed=held)
    end if
    series(m, 1:count, b) = w

  end subroutine one_order

  !===================================================================!
  ! What the earlier blocks found at one order, at the instants a
  ! later one carries. A coefficient travels the junction the way the
  ! trajectory does, and comes from whichever block computed that
  ! instant.
  !===================================================================!

  pure function coefficients_handed(earlier, first, stride, given, series, order) &
       & result(held)

    type(chain_block), intent(in) :: earlier(:)
    integer          , intent(in) :: first, stride, given, order
    real(dp)         , intent(in) :: series(0:, :, :)
    real(dp), allocatable :: held(:)

    integer :: i, b, local, at, width

    width = earlier(1) % width
    allocate(held(given * width), source=0.0_dp)

    do i = 1, given
       call locate(earlier, first + (i - 1) * stride, b, local)
       if (b == 0) cycle
       at = earlier(b) % instants_at(local)
       held((i - 1) * width + 1:i * width) = series(order, at + 1:at + width, b)
    end do

  end function coefficients_handed

  !===================================================================!
  ! The functional at every order: each block's own points, each
  ! counted once, weighted by the step that ends at its instant times
  ! the measure of its node.
  !===================================================================!

  subroutine chain_functional(chain, functionals, degrees, design, max_order, &
       & series, f, node_measure)

    type(chain_block)      , intent(in) :: chain(:)
    type(functional_holder), intent(in) :: functionals(:)
    integer                , intent(in) :: degrees, max_order
    real(dp)               , intent(in) :: design, series(0:, :, :)
    real(dp), allocatable  , intent(out) :: f(:,:)
    real(dp), intent(in), optional      :: node_measure(:)

    real(dp), allocatable :: values(:), weight(:)
    integer , allocatable :: at(:)
    integer :: b, m, i, from, to, count

    allocate(f(0:max_order, size(functionals)), source=0.0_dp)

    do b = 1, size(chain)
       call owned(chain, b, from, to)
       if (from > to) cycle
       count = chain(b) % rows % num_unknowns()
       call owned_points(chain(b), from, to, node_measure, at, weight)
       do i = 1, size(functionals)
          do m = 0, max_order
             call nodal_coefficient(functionals(i) % rule, degrees, at, series(:, 1:count, b), &
                  & design, m, values)
             f(m, i) = f(m, i) + sum(weight * values)
          end do
       end do
    end do

  end subroutine chain_functional

  !===================================================================!
  ! Every block's system at the trajectory already marched. The
  ! gradient is shared out by ownership, each block taking the
  ! instants it computed and none of the instants it was given, so
  ! that a shared instant counts once.
  !===================================================================!

  subroutine chain_systems(chain, tower, functionals, degrees, systems, node_measure)

    type(chain_block)      , intent(in) :: chain(:)
    type(expansion)        , intent(in) :: tower
    type(functional_holder), intent(in) :: functionals(:)
    integer                , intent(in) :: degrees
    type(chain_system), allocatable, intent(out) :: systems(:)
    real(dp), intent(in), optional      :: node_measure(:)

    type(stored_directed_graph) :: unknowns, points
    type(stored_field), allocatable :: inputs(:)
    type(stored_field) :: state, knobs
    real(dp), allocatable :: rate(:), g(:), weight(:), varied_weight(:), values(:), series(:,:)
    real(dp), allocatable :: step_partials(:,:)
    integer , allocatable :: at(:), varied_at(:)
    real(dp) :: design
    integer :: b, i, j, count, num_designs, num_functionals, from, to

    ! THE DESIGNS, read from the tower: the physics' parameter first,
    ! then, when the steps are designs, one column per weight, the
    ! steps' partials in them from the tower's own grid
    call designs_of(tower, design, step_partials)
    num_functionals = size(functionals)
    num_designs     = 1
    if (allocated(step_partials)) num_designs = 1 + size(step_partials, 2)

    allocate(systems(size(chain)))
    do b = 1, size(chain)
       count = chain(b) % rows % num_unknowns()
       call frozen_at(chain(b), design, unknowns, inputs)
       systems(b) % mark = fresh_stamp()
       allocate(systems(b) % rate(count, num_designs), source=0.0_dp)
       allocate(systems(b) % g(count, num_functionals), source=0.0_dp)
       allocate(systems(b) % explicit(num_functionals, num_designs), source=0.0_dp)

       ! the block's partial in the physics' design, then in each
       ! entry of the grid's: its rows differentiated along the steps'
       ! partial in that entry, applied to its state
       call design_partial(chain(b) % rows, unknowns, inputs, &
            & chain(b) % rows % num_points(), unknowns % vertex_set(), rate)
       systems(b) % rate(:, 1) = rate
       do j = 2, num_designs
          systems(b) % rate(:, j) = varied_rate(chain(b), degrees, &
               & along_of(chain(b), step_partials(:, j - 1)))
       end do

       ! every functional: its gradient over the owned points, its own
       ! partial in the physics' design, and the measure's partial in
       ! each entry of the grid's
       ! a block with no owned points - a startup - has no gradient
       ! and nothing to add
       call owned(chain, b, from, to)
       if (from > to) cycle
       call owned_points(chain(b), from, to, node_measure, at, weight)
       points = stored_directed_graph(size(at), tails=[integer ::], heads=[integer ::])
       call at_owned_points(inputs, degrees, design, at, points, state, knobs)
       series = reshape(chain(b) % state, [1, count])
       do i = 1, num_functionals
          call owned_gradient(chain, b, functionals(i) % rule, degrees, design, inputs, &
               & unknowns, g, node_measure)
          systems(b) % g(:, i) = g
          call functional_design_partial(functionals(i) % rule, points, [state, knobs], &
               & weight, size(at), points % vertex_set(), systems(b) % explicit(i, 1))
          if (num_designs > 1) then
             call nodal_coefficient(functionals(i) % rule, degrees, at, series, design, 0, values)
          end if
          do j = 2, num_designs
             ! the measure's partial: the points again, weighted by the
             ! steps' partial in this entry in place of the steps
             call owned_points(chain(b), from, to, node_measure, varied_at, varied_weight, &
                  & along=along_of(chain(b), step_partials(:, j - 1)))
             systems(b) % explicit(i, j) = sum(varied_weight * values)
          end do
       end do
    end do

  end subroutine chain_systems

  !===================================================================!
  ! What the tower says the designs are: the physics' parameter, and
  ! the steps' partials in the weights when the weights are designs.
  ! Invalid input: a tower whose first design is not the parameter.
  !===================================================================!

  subroutine designs_of(tower, design, step_partials)

    type(expansion), intent(in) :: tower
    real(dp)       , intent(out) :: design
    real(dp), allocatable, intent(out) :: step_partials(:,:)

    integer :: k

    if (tower % design_kind_of(1) /= design_of_physics) then
       error stop 'gti_chain: the physics'' parameter is the first design'
    end if
    design = tower % parameter()
    do k = 2, tower % num_designs()
       if (tower % design_kind_of(k) == design_of_steps) call tower % step_partials(step_partials)
    end do

  end subroutine designs_of

  !===================================================================!
  ! A block's rows differentiated along a direction in its steps,
  ! applied to its state: the block's partial in whatever the steps
  ! read, by the chain rule through the steps. The carried rows, the
  ! physics and the spatial discretization stencil read no step and take no part.
  !===================================================================!

  function varied_rate(b, degrees, along) result(r)

    type(chain_block), intent(in) :: b
    integer          , intent(in) :: degrees
    real(dp)         , intent(in) :: along(:)
    real(dp), allocatable :: r(:)

    r = varied_applied(b, degrees, along, b % state)

  end function varied_rate

  !===================================================================!
  ! One functional, held; and the first entry of a table of
  ! derivatives, for a caller with one functional and one design.
  !===================================================================!

  function one_functional(rule) result(holder)

    type(expression)      , intent(in) :: rule
    type(functional_holder) :: holder

    holder % rule = rule

  end function one_functional

  pure real(dp) function first_of(table)

    real(dp), intent(in) :: table(:,:)

    first_of = table(1, 1)

  end function first_of

  !===================================================================!
  ! The points one block owns and the measure of each: the owned
  ! instants, node by node, weighted by the block's own step ending
  ! at the instant times the node's measure - one when no measure is
  ! given, as for one node's equation. Invalid input: a measure for
  ! other than every node.
  !===================================================================!

  subroutine owned_points(b, from, to, node_measure, at, weight, along)

    type(chain_block), intent(in) :: b
    integer          , intent(in) :: from, to
    real(dp), intent(in), optional :: node_measure(:)
    integer , allocatable, intent(out) :: at(:)
    real(dp), allocatable, intent(out) :: weight(:)
    ! given, a direction in the block's steps in place of the steps
    real(dp), intent(in), optional :: along(:)

    real(dp), allocatable :: measure(:), step(:)
    integer :: degrees, k, i

    degrees = b % width / b % nodes
    if (present(node_measure)) then
       if (size(node_measure) /= b % nodes) then
          error stop 'gti_chain: one measure per node'
       end if
       measure = node_measure
    else
       measure = spread(1.0_dp, 1, b % nodes)
    end if
    if (present(along)) then
       step = along
    else
       step = b % dt
    end if

    at     = [((b % instants_at(k) + (i - 1) * degrees, i = 1, b % nodes), k = from, to)]
    weight = [((step(k) * measure(i), i = 1, b % nodes), k = from, to)]

  end subroutine owned_points

  !===================================================================!
  ! A direction in the march's steps read on a block's own: each of
  ! its steps takes the direction of the march's step it lies in,
  ! times its fraction of that step.
  !===================================================================!

  pure function along_of(b, v) result(along)

    type(chain_block), intent(in) :: b
    real(dp)         , intent(in) :: v(:)
    real(dp) :: along(size(b % dt))

    integer :: k

    along(1) = 0.0_dp
    do k = 2, size(b % dt)
       along(k) = v(b % coarse_step(k)) * b % fraction
    end do

  end function along_of

  !===================================================================!
  ! The functional's gradient over the points one block owns - its
  ! instants, node by node over a field - each weighted by the step
  ! that ends at its instant times the measure of its node. The
  ! integrand reads one point at a time, so one partial action per
  ! degree gives the whole of it rather than one per unknown.
  !===================================================================!

  subroutine owned_gradient(chain, b, integrand, degrees, design, inputs, &
       & unknowns, g, node_measure)

    type(chain_block)          , intent(in) :: chain(:)
    integer                    , intent(in) :: b, degrees
    type(expression)           , intent(in) :: integrand
    real(dp)                   , intent(in) :: design
    type(stored_field)         , intent(in) :: inputs(:)
    type(stored_directed_graph), intent(in) :: unknowns
    real(dp), allocatable      , intent(out) :: g(:)
    real(dp), intent(in), optional          :: node_measure(:)

    type(stored_directed_graph) :: points
    type(stored_field) :: state, knobs
    real(dp), allocatable :: owned_g(:), weight(:)
    integer , allocatable :: at(:)
    integer :: count, from, to, p, d

    count = chain(b) % rows % num_unknowns()
    allocate(g(count), source=0.0_dp)

    call owned(chain, b, from, to)
    if (from > to) return
    call owned_points(chain(b), from, to, node_measure, at, weight)
    points = stored_directed_graph(size(at), tails=[integer ::], heads=[integer ::])
    call at_owned_points(inputs, degrees, design, at, points, state, knobs)

    ! the gradient over the owned points, then scattered to where
    ! those points lie in the block
    call functional_gradient(integrand, points, [state, knobs], weight, size(at), degrees, &
         & points % vertex_set(), owned_g)
    do p = 1, size(at)
       do d = 0, degrees - 1
          g(at(p) + d + 1) = owned_g((p - 1) * degrees + d + 1)
       end do
    end do
    associate (u1 => unknowns); end associate

  end subroutine owned_gradient

  !===================================================================!
  ! The trajectory at the points one block owns, laid out one point
  ! at a time so that a nodal rule reads them.
  !===================================================================!

  subroutine at_owned_points(inputs, degrees, design, at, points, state, knobs)

    type(stored_field)         , intent(in)  :: inputs(:)
    integer                    , intent(in)  :: degrees, at(:)
    real(dp)                   , intent(in)  :: design
    type(stored_directed_graph), intent(in)  :: points
    type(stored_field)         , intent(out) :: state, knobs

    real(dp), allocatable :: whole(:), v(:)
    integer :: p

    call inputs(1) % real_vector(whole)
    allocate(v(size(at) * degrees))
    do p = 1, size(at)
       v((p - 1) * degrees + 1:p * degrees) = whole(at(p) + 1:at(p) + degrees)
    end do

    state = stored_field('state', points % vertex_set(), size(at) * degrees)
    knobs = stored_field('design', points % vertex_set(), size(at))
    call state % set_real_vector(v)
    call knobs % set_real_vector(spread(design, 1, size(at)))

  end subroutine at_owned_points

  !===================================================================!
  ! Which block holds one fine instant, the latest that does, and
  ! where in it.
  !===================================================================!

  pure subroutine holder_of(chain, fine, held_by, at)

    type(chain_block), intent(in)  :: chain(:)
    integer          , intent(in)  :: fine
    integer          , intent(out) :: held_by, at

    integer :: local

    call locate(chain, fine, held_by, local)
    at = 0
    if (held_by > 0) at = chain(held_by) % instants_at(local)

  end subroutine holder_of

  !===================================================================!
  ! THE SECOND DERIVATIVES BY THE REVERSE ROUTE: one row of every
  ! functional's hessian per design, from the direction's tangent and
  ! a second-order costate.
  !
  ! With R(q, p) = 0 and f(q, p), the first derivative along the
  ! costate lambda (A^T lambda = f_q) reads f_p - lambda^T R_p. Its
  ! derivative in design j, with the state moving by the tangent w_j
  ! (A w_j = -R_pj), is
  !
  !    H_jk = f_pk_q w_j + f_pk_pj - lambda_j'^T R_pk
  !           - lambda^T (R_pk_q w_j + R_pk_pj)
  !
  ! where the costate's own derivative solves
  !
  !    A^T lambda_j' = f_qq w_j + f_q_pj - (R_qq w_j + R_q_pj)^T lambda.
  !
  ! Every second partial is exact: the physics' from its rule over
  ! two directions, the rows' from the weight action along two step
  ! directions, and the steps' from the grid. The rows are linear in
  ! the state, so R_qq is the physics' alone and R_q_pj is the varied
  ! rows for a grid design; through a nodal physics, a transposed
  ! contraction is local to the point. Along a chain the tangents are
  ! handed forward and the costates back exactly as at first order.
  ! The cost is D tangent solves, F costates and F D second-order
  ! costates, which is what the gate counts for the reverse route at
  ! order two. The table is symmetric in theory, and is not made so.
  !===================================================================!

  subroutine chain_hessian(chain, tower, systems, functionals, degrees, hessian, node_measure)

    type(chain_block)      , intent(in) :: chain(:)
    type(expansion)        , intent(in) :: tower
    type(chain_system)     , intent(in) :: systems(:)
    type(functional_holder), intent(in) :: functionals(:)
    integer                , intent(in) :: degrees
    real(dp), allocatable  , intent(out) :: hessian(:,:,:)
    real(dp), intent(in), optional      :: node_measure(:)

    type(stored_directed_graph) :: unknowns, points
    type(stored_field), allocatable :: inputs(:)
    type(stored_field) :: state, knobs
    real(dp), allocatable :: w(:,:,:), lambda(:,:,:), rhs(:,:), mu(:), explicit(:,:)
    real(dp), allocatable :: r2(:), g2(:), weight(:), values(:), series(:,:), step_partials(:,:)
    integer , allocatable :: at(:)
    real(dp) :: design
    integer :: nf, nd, b, i, j, k, count, widest, from, to, p, d, instant, held_by, where

    call designs_of(tower, design, step_partials)
    nf = size(functionals)
    nd = size(systems(1) % rate, 2)
    if (nd > 1 .and. .not. allocated(step_partials)) then
       error stop 'gti_chain: the systems name more designs than the tower holds'
    end if
    widest = 0
    do b = 1, size(chain)
       widest = max(widest, chain(b) % rows % num_unknowns())
    end do

    ! every direction's tangent and every functional's costate
    call tangents(chain, systems, degrees, design, widest, w)
    call costates(chain, systems, degrees, design, widest, lambda)

    allocate(hessian(nf, nd, nd), source=0.0_dp)
    allocate(rhs(widest, size(chain)))
    do i = 1, nf
       do j = 1, nd
          ! the second-order costate's right side, block by block:
          ! the functional's second partials along the tangent and in
          ! the design, less the transposed contraction of the rows'
          rhs = 0.0_dp
          do b = 1, size(chain)
             count = chain(b) % rows % num_unknowns()
             call frozen_at(chain(b), design, unknowns, inputs)
             call rows_second_transposed(chain(b), degrees, unknowns, inputs, &
                  & w(1:count, b, j), lambda(1:count, b, i), j, step_partials, mu)
             rhs(1:count, b) = -mu
             call owned(chain, b, from, to)
             if (from > to) cycle
             call owned_points(chain(b), from, to, node_measure, at, weight)
             points = stored_directed_graph(size(at), tails=[integer ::], heads=[integer ::])
             call at_owned_points(inputs, degrees, design, at, points, state, knobs)
             call functional_second(functionals(i) % rule, points, state, knobs, weight, at, &
                  & degrees, count, w(1:count, b, j), j, step_partials, chain(b), from, to, &
                  & node_measure, g2)
             rhs(1:count, b) = rhs(1:count, b) + g2
          end do
          ! the costate's derivative, back along the chain
          do b = size(chain), 1, -1
             count = chain(b) % rows % num_unknowns()
             call frozen_at(chain(b), design, unknowns, inputs)
             call solved_linear(chain(b) % rows, unknowns, inputs, rhs(1:count, b), .true., &
                  & systems(b) % mark, mu)
             do k = 1, nd
                hessian(i, j, k) = hessian(i, j, k) - dot_product(mu, systems(b) % rate(:, k))
             end do
             do p = 1, chain(b) % given * chain(b) % width
                instant = chain(b) % first + ((p - 1) / chain(b) % width) * chain(b) % stride
                d       = mod(p - 1, chain(b) % width)
                call holder_of(chain(1:b - 1), instant, held_by, where)
                if (held_by > 0) rhs(where + d + 1, held_by) = rhs(where + d + 1, held_by) + mu(p)
             end do
          end do
          ! the rest of the row: the costate against the rows' mixed
          ! partials, and the functional's own mixed partials over the
          ! points the block owns
          do b = 1, size(chain)
             count = chain(b) % rows % num_unknowns()
             call frozen_at(chain(b), design, unknowns, inputs)
             do k = 1, nd
                call rows_mixed(chain(b), degrees, unknowns, inputs, w(1:count, b, j), j, k, &
                     & step_partials, tower, r2)
                hessian(i, j, k) = hessian(i, j, k) - dot_product(lambda(1:count, b, i), r2)
             end do
             call owned(chain, b, from, to)
             if (from > to) cycle
             call owned_points(chain(b), from, to, node_measure, at, weight)
             points = stored_directed_graph(size(at), tails=[integer ::], heads=[integer ::])
             call at_owned_points(inputs, degrees, design, at, points, state, knobs)
             series = reshape(chain(b) % state, [1, count])
             call nodal_coefficient(functionals(i) % rule, degrees, at, series, design, 0, values)
             do k = 1, nd
                call functional_mixed(functionals(i) % rule, points, state, knobs, weight, at, &
                     & degrees, count, w(1:count, b, j), j, k, step_partials, tower, &
                     & chain(b), from, to, node_measure, values, explicit)
                hessian(i, j, k) = hessian(i, j, k) + explicit(1, 1)
             end do
          end do
       end do
    end do

  end subroutine chain_hessian

  !===================================================================!
  ! Every design's tangent over every block, handed forward.
  !===================================================================!

  subroutine tangents(chain, systems, degrees, design, widest, w)

    type(chain_block) , intent(in) :: chain(:)
    type(chain_system), intent(in) :: systems(:)
    integer           , intent(in) :: degrees, widest
    real(dp)          , intent(in) :: design
    real(dp), allocatable, intent(out) :: w(:,:,:)

    type(stored_directed_graph) :: unknowns
    type(stored_field), allocatable :: inputs(:)
    real(dp), allocatable :: rhs(:), one(:)
    integer :: b, i, j, k, d, held_by, at, nd

    nd = size(systems(1) % rate, 2)
    allocate(w(widest, size(chain), nd), source=0.0_dp)
    do j = 1, nd
       do b = 1, size(chain)
          rhs = -systems(b) % rate(:, j)
          do i = 1, chain(b) % given * chain(b) % width
             k = chain(b) % first + ((i - 1) / chain(b) % width) * chain(b) % stride
             d = mod(i - 1, chain(b) % width)
             call holder_of(chain(1:b - 1), k, held_by, at)
             if (held_by > 0) rhs(i) = w(at + d + 1, held_by, j)
          end do
          call frozen_at(chain(b), design, unknowns, inputs)
          call solved_linear(chain(b) % rows, unknowns, inputs, rhs, .false., systems(b) % mark, one)
          w(1:size(one), b, j) = one
       end do
    end do

  end subroutine tangents

  !===================================================================!
  ! Every functional's costate over every block, handed back.
  !===================================================================!

  subroutine costates(chain, systems, degrees, design, widest, lambda)

    type(chain_block) , intent(in) :: chain(:)
    type(chain_system), intent(in) :: systems(:)
    integer           , intent(in) :: degrees, widest
    real(dp)          , intent(in) :: design
    real(dp), allocatable, intent(out) :: lambda(:,:,:)

    type(stored_directed_graph) :: unknowns
    type(stored_field), allocatable :: inputs(:)
    real(dp), allocatable :: rhs(:,:), one(:)
    integer :: b, i, k, d, held_by, at, instant, nf, count

    nf = size(systems(1) % g, 2)
    allocate(lambda(widest, size(chain), nf), source=0.0_dp)
    allocate(rhs(widest, size(chain)))
    do i = 1, nf
       rhs = 0.0_dp
       do b = 1, size(chain)
          rhs(1:size(systems(b) % g, 1), b) = systems(b) % g(:, i)
       end do
       do b = size(chain), 1, -1
          count = chain(b) % rows % num_unknowns()
          call frozen_at(chain(b), design, unknowns, inputs)
          call solved_linear(chain(b) % rows, unknowns, inputs, rhs(1:count, b), .true., &
               & systems(b) % mark, one)
          lambda(1:count, b, i) = one
          do k = 1, chain(b) % given * chain(b) % width
             instant = chain(b) % first + ((k - 1) / chain(b) % width) * chain(b) % stride
             d       = mod(k - 1, chain(b) % width)
             call holder_of(chain(1:b - 1), instant, held_by, at)
             if (held_by > 0) rhs(at + d + 1, held_by) = rhs(at + d + 1, held_by) + one(k)
          end do
       end do
    end do

  end subroutine costates

  !===================================================================!
  ! A functional's second partials over one block's owned points, as
  ! a right side over the block's unknowns: f_qq w_j + f_q_pj, the
  ! latter the integrand's mixed partial in the physics' design or,
  ! for a grid design, the gradient weighted by the steps' partial.
  !===================================================================!

  subroutine functional_second(rule, points, state, knobs, weight, at, degrees, count, w, j, &
       & step_partials, b, from, to, node_measure, g2)

    type(expression)           , intent(in) :: rule
    type(stored_directed_graph), intent(in) :: points
    type(stored_field)         , intent(in) :: state, knobs
    real(dp)                   , intent(in) :: weight(:), w(:)
    integer                    , intent(in) :: at(:), degrees, count, j, from, to
    real(dp), intent(in), optional          :: step_partials(:,:), node_measure(:)
    type(chain_block)          , intent(in) :: b
    real(dp), allocatable      , intent(out) :: g2(:)

    real(dp), allocatable :: owned_g(:), mixed(:), dweight(:), along(:)
    integer , allocatable :: varied_at(:)
    integer :: p, d

    along = [((w(at(p) + d + 1), d = 0, degrees - 1), p = 1, size(at))]
    call functional_gradient(rule, points, [state, knobs], weight, size(at), degrees, &
         & points % vertex_set(), owned_g, along_state=along)
    if (j == 1) then
       call functional_gradient(rule, points, [state, knobs], weight, size(at), degrees, &
            & points % vertex_set(), mixed, along_design=spread(1.0_dp, 1, size(at)))
    else
       call owned_points(b, from, to, node_measure, varied_at, dweight, &
            & along=along_of(b, step_partials(:, j - 1)))
       call functional_gradient(rule, points, [state, knobs], dweight, size(at), degrees, &
            & points % vertex_set(), mixed)
    end if
    allocate(g2(count), source=0.0_dp)
    do p = 1, size(at)
       do d = 0, degrees - 1
          g2(at(p) + d + 1) = owned_g((p - 1) * degrees + d + 1) + mixed((p - 1) * degrees + d + 1)
       end do
    end do

  end subroutine functional_second

  !===================================================================!
  ! The transposed contraction (R_qq w_j + R_q_pj)^T lambda over one
  ! block. The physics is nodal, so at each point the row's second
  ! partial in its own degrees along w_j, and its mixed partial in
  ! the design, are contracted with that row's costate; the varied
  ! rows of a grid design are applied transposed.
  !===================================================================!

  subroutine rows_second_transposed(b, degrees, unknowns, inputs, w, lambda, j, step_partials, mu)

    type(chain_block)          , intent(in) :: b
    integer                    , intent(in) :: degrees, j
    type(stored_directed_graph), intent(in) :: unknowns
    type(stored_field)         , intent(in) :: inputs(:)
    real(dp)                   , intent(in) :: w(:), lambda(:)
    real(dp), intent(in), optional          :: step_partials(:,:)
    real(dp), allocatable      , intent(out) :: mu(:)

    integer , allocatable :: at(:)
    real(dp), allocatable :: e(:), second(:)
    integer :: count, d, p, primary

    count   = size(w)
    at      = b % rows % points_at()
    primary = b % primary
    allocate(mu(count), source=0.0_dp)
    allocate(e(count))
    do d = 0, degrees - 1
       ! the unit direction at this degree of every point
       e = 0.0_dp
       do p = 1, size(at)
          e(at(p) + d + 1) = 1.0_dp
       end do
       call varied_by(b % rows, unknowns, inputs, 1, unknowns % vertex_set(), e, &
            & 1, unknowns % vertex_set(), w, second)
       do p = 1, size(at)
          mu(at(p) + d + 1) = mu(at(p) + d + 1) + lambda(at(p) + primary + 1) * second(at(p) + primary + 1)
       end do
       if (j == 1) then
          call varied_by(b % rows, unknowns, inputs, 1, unknowns % vertex_set(), e, &
               & 2, unknowns % vertex_set(), spread(1.0_dp, 1, b % rows % num_points()), second)
          do p = 1, size(at)
             mu(at(p) + d + 1) = mu(at(p) + d + 1) + lambda(at(p) + primary + 1) * second(at(p) + primary + 1)
          end do
       end if
    end do
    if (j > 1) then
       mu = mu + varied_transposed(b, degrees, along_of(b, step_partials(:, j - 1)), lambda)
    end if

  end subroutine rows_second_transposed

  !===================================================================!
  ! A functional's mixed second partial f_pk_q w_j + f_pk_pj over one
  ! block's owned points: in the physics' design, the integrand's
  ! own partials; in a grid design, the gradient or the value
  ! weighted by the steps' partial, and by their mixed second partial
  ! for two grid designs.
  !===================================================================!

  subroutine functional_mixed(rule, points, state, knobs, weight, at, degrees, count, w, j, k, &
       & step_partials, tower, b, from, to, node_measure, values, explicit)

    type(expression)           , intent(in) :: rule
    type(stored_directed_graph), intent(in) :: points
    type(stored_field)         , intent(in) :: state, knobs
    real(dp)                   , intent(in) :: weight(:), w(:), values(:)
    integer                    , intent(in) :: at(:), degrees, count, j, k, from, to
    real(dp)   , intent(in), optional       :: step_partials(:,:), node_measure(:)
    type(expansion)            , intent(in) :: tower
    type(chain_block)          , intent(in) :: b
    real(dp), allocatable      , intent(out) :: explicit(:,:)

    real(dp), allocatable :: along(:), g(:), dweight(:), u(:), uweight(:)
    integer , allocatable :: varied_at(:)
    real(dp) :: part
    integer :: p, d

    allocate(explicit(1, 1), source=0.0_dp)
    along = [((w(at(p) + d + 1), d = 0, degrees - 1), p = 1, size(at))]
    if (k == 1) then
       ! f_nu_q w_j, then f_nu_pj
       call functional_design_partial(rule, points, [state, knobs], weight, size(at), &
            & points % vertex_set(), part, along_state=along)
       explicit(1, 1) = part
       if (j == 1) then
          call functional_design_partial(rule, points, [state, knobs], weight, size(at), &
               & points % vertex_set(), part, along_design=spread(1.0_dp, 1, size(at)))
          explicit(1, 1) = explicit(1, 1) + part
       else
          call owned_points(b, from, to, node_measure, varied_at, dweight, &
               & along=along_of(b, step_partials(:, j - 1)))
          call functional_design_partial(rule, points, [state, knobs], dweight, size(at), &
               & points % vertex_set(), part)
          explicit(1, 1) = explicit(1, 1) + part
       end if
    else
       ! the measure's partial in design k against the gradient along
       ! w_j, then its mixed partial in j and k
       call owned_points(b, from, to, node_measure, varied_at, dweight, &
            & along=along_of(b, step_partials(:, k - 1)))
       call functional_gradient(rule, points, [state, knobs], dweight, size(at), degrees, &
            & points % vertex_set(), g)
       explicit(1, 1) = dot_product(g, along)
       if (j == 1) then
          call functional_design_partial(rule, points, [state, knobs], dweight, size(at), &
               & points % vertex_set(), part)
          explicit(1, 1) = explicit(1, 1) + part
       else
          call tower % step_second_partials(j - 1, k - 1, u)
          call owned_points(b, from, to, node_measure, varied_at, uweight, along=along_of(b, u))
          explicit(1, 1) = explicit(1, 1) + sum(uweight * values)
       end if
    end if
    associate (u1 => count); end associate

  end subroutine functional_mixed

  !===================================================================!
  ! The rows' mixed second partial R_pk_q w_j + R_pk_pj over one
  ! block: in the physics' design the physics' own, in a grid design
  ! the varied rows applied to the tangent, and for two grid designs
  ! the rows along both step directions and along the steps' mixed
  ! second partial, applied to the state.
  !===================================================================!

  subroutine rows_mixed(b, degrees, unknowns, inputs, w, j, k, step_partials, tower, r2)

    type(chain_block)          , intent(in) :: b
    integer                    , intent(in) :: degrees, j, k
    type(stored_directed_graph), intent(in) :: unknowns
    type(stored_field)         , intent(in) :: inputs(:)
    real(dp)                   , intent(in) :: w(:)
    real(dp)   , intent(in), optional       :: step_partials(:,:)
    type(expansion)            , intent(in) :: tower
    real(dp), allocatable      , intent(out) :: r2(:)

    real(dp), allocatable :: second(:), u(:)
    integer :: count, npts

    count = size(w)
    npts  = b % rows % num_points()
    if (k == 1) then
       call varied_by(b % rows, unknowns, inputs, 2, unknowns % vertex_set(), &
            & spread(1.0_dp, 1, npts), 1, unknowns % vertex_set(), w, r2)
       if (j == 1) then
          call varied_by(b % rows, unknowns, inputs, 2, unknowns % vertex_set(), &
               & spread(1.0_dp, 1, npts), 2, unknowns % vertex_set(), spread(1.0_dp, 1, npts), second)
          r2 = r2 + second
       end if
    else
       r2 = varied_applied(b, degrees, along_of(b, step_partials(:, k - 1)), w)
       if (j > 1) then
          r2 = r2 + varied_applied(b, degrees, along_of(b, step_partials(:, j - 1)), &
               & b % state, along2=along_of(b, step_partials(:, k - 1)))
          call tower % step_second_partials(j - 1, k - 1, u)
          r2 = r2 + varied_applied(b, degrees, along_of(b, u), b % state)
       end if
    end if

  end subroutine rows_mixed

  !===================================================================!
  ! The rows along a direction in the steps - and a second, given -
  ! applied to a vector, and applied transposed. A carried row holds
  ! a given number whatever the steps, so the block's residual has no
  ! such row there and neither has its derivative: a time discretization stencil row the
  ! family would have placed at a carried instant is dropped, as the
  ! block drops it.
  !===================================================================!

  function varied_applied(b, degrees, along, x, along2) result(r)

    type(chain_block), intent(in) :: b
    integer          , intent(in) :: degrees
    real(dp)         , intent(in) :: along(:), x(:)
    real(dp)         , intent(in), optional :: along2(:)
    real(dp), allocatable :: r(:)

    type(stencil) :: rows
    real(dp), allocatable :: w(:)
    integer :: e

    rows = b % rows % rows_varied(b % scheme, b % dt, along, along2)
    call rows % weights % real_vector(w)
    allocate(r(size(x)), source=0.0_dp)
    associate (u1 => degrees); end associate
    do e = 1, rows % pattern % num_edges()
       r(rows % pattern % edge_head(e)) = r(rows % pattern % edge_head(e)) &
            & + w(e) * x(rows % pattern % edge_tail(e))
    end do
    r(b % rows % carried_unknowns()) = 0.0_dp

  end function varied_applied

  function varied_transposed(b, degrees, along, y) result(r)

    type(chain_block), intent(in) :: b
    integer          , intent(in) :: degrees
    real(dp)         , intent(in) :: along(:), y(:)
    real(dp), allocatable :: r(:)

    type(stencil) :: rows
    real(dp), allocatable :: w(:), kept(:)
    integer :: e

    rows = b % rows % rows_varied(b % scheme, b % dt, along)
    call rows % weights % real_vector(w)
    kept = y
    associate (u1 => degrees); end associate
    kept(b % rows % carried_unknowns()) = 0.0_dp
    allocate(r(size(y)), source=0.0_dp)
    do e = 1, rows % pattern % num_edges()
       r(rows % pattern % edge_tail(e)) = r(rows % pattern % edge_tail(e)) &
            & + w(e) * kept(rows % pattern % edge_head(e))
    end do

  end function varied_transposed

  !===================================================================!
  ! What the model says an expansion to the given order costs in
  ! substitutions, for the accounting layer to set beside what it
  ! counted: per order, one per block by the forward route at one
  ! design and one functional.
  !===================================================================!

  pure integer function expansion_substitutions(num_blocks, order) result(count)

    integer, intent(in) :: num_blocks, order

    count = num_blocks * route_substitutions(route_of(1, 1, order), 1, 1, order)

  end function expansion_substitutions

  !===================================================================!
  ! Forward: one solve per design. Each block solves for its own
  ! sensitivity, with the rows it carried set to what an earlier block
  ! already found at those instants, and every functional reads it.
  ! The table has one row per functional and one column per design.
  !===================================================================!

  function chain_by_tangent(chain, systems, degrees, design) result(df)

    type(chain_block) , intent(in) :: chain(:)
    type(chain_system), intent(in) :: systems(:)
    integer           , intent(in) :: degrees
    real(dp)          , intent(in) :: design
    real(dp), allocatable :: df(:,:)

    type(stored_directed_graph) :: unknowns
    type(stored_field), allocatable :: inputs(:)
    real(dp), allocatable :: w(:,:), rhs(:), one(:)
    integer :: b, widest, i, j, k, d, held_by, at, num_designs, num_functionals

    num_functionals = size(systems(1) % g, 2)
    num_designs     = size(systems(1) % rate, 2)
    widest = 0
    do b = 1, size(chain)
       widest = max(widest, chain(b) % rows % num_unknowns())
    end do
    allocate(w(widest, size(chain)), source=0.0_dp)
    allocate(df(num_functionals, num_designs), source=0.0_dp)
    do b = 1, size(chain)
       df = df + systems(b) % explicit
    end do

    ! one sensitivity per design, marched along the chain
    do j = 1, num_designs
       w = 0.0_dp
       do b = 1, size(chain)
          rhs = -systems(b) % rate(:, j)
          do i = 1, chain(b) % given * chain(b) % width
             k = chain(b) % first + ((i - 1) / chain(b) % width) * chain(b) % stride
             d = mod(i - 1, chain(b) % width)
             call holder_of(chain(1:b - 1), k, held_by, at)
             if (held_by > 0) rhs(i) = w(at + d + 1, held_by)
          end do
          call frozen_at(chain(b), design, unknowns, inputs)
          call solved_linear(chain(b) % rows, unknowns, inputs, rhs, .false., systems(b) % mark, one)
          w(1:size(one), b) = one
          do i = 1, num_functionals
             df(i, j) = df(i, j) + dot_product(systems(b) % g(:, i), one)
          end do
       end do
    end do

  end function chain_by_tangent

  !===================================================================!
  ! Backward: one solve per functional. Each block solves against its
  ! own transpose, with its own share of that functional's gradient
  ! plus whatever a later block's costate left on the instants they
  ! share, and every design reads the costate. The same table.
  !===================================================================!

  function chain_by_adjoint(chain, systems, degrees, design) result(df)

    type(chain_block) , intent(in) :: chain(:)
    type(chain_system), intent(in) :: systems(:)
    integer           , intent(in) :: degrees
    real(dp)          , intent(in) :: design
    real(dp), allocatable :: df(:,:)

    type(stored_directed_graph) :: unknowns
    type(stored_field), allocatable :: inputs(:)
    real(dp), allocatable :: rhs(:,:), lambda(:)
    integer :: b, widest, i, j, k, d, held_by, at, instant, num_designs, num_functionals

    num_functionals = size(systems(1) % g, 2)
    num_designs     = size(systems(1) % rate, 2)
    widest = 0
    do b = 1, size(chain)
       widest = max(widest, chain(b) % rows % num_unknowns())
    end do
    allocate(rhs(widest, size(chain)), source=0.0_dp)
    allocate(df(num_functionals, num_designs), source=0.0_dp)
    do b = 1, size(chain)
       df = df + systems(b) % explicit
    end do

    ! one costate per functional, marched back along the chain
    do i = 1, num_functionals
       rhs = 0.0_dp
       do b = 1, size(chain)
          rhs(1:size(systems(b) % g, 1), b) = systems(b) % g(:, i)
       end do
       do b = size(chain), 1, -1
          call frozen_at(chain(b), design, unknowns, inputs)
          call solved_linear(chain(b) % rows, unknowns, inputs, &
               & rhs(1:chain(b) % rows % num_unknowns(), b), .true., systems(b) % mark, lambda)
          do j = 1, num_designs
             df(i, j) = df(i, j) - dot_product(lambda, systems(b) % rate(:, j))
          end do
          do k = 1, chain(b) % given * chain(b) % width
             instant = chain(b) % first + ((k - 1) / chain(b) % width) * chain(b) % stride
             d       = mod(k - 1, chain(b) % width)
             call holder_of(chain(1:b - 1), instant, held_by, at)
             if (held_by > 0) rhs(at + d + 1, held_by) = rhs(at + d + 1, held_by) + lambda(k)
          end do
       end do
    end do

  end function chain_by_adjoint

end module gti_chain
