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
  use physics_integrand, only : nodal_integrand
  use gti_expansion    , only : family_holder, marches_by_stages
  use gti_block        , only : block_residual
  use gti_march        , only : imbalance, swept, solved_linear, fresh_stamp, partitioned, horizon_bounds, block_of, solved, &
       & frozen_inputs
  use gti_stage        , only : stage_block_of, instant_at
  use view_directed_stored, only : stored_directed_graph
  use field_calculus   , only : field
  use field_stored     , only : stored_field
  use operation_action , only : variation
  use gti_sweeps       , only : design_partial, route_of, functional_gradient, &
       & route_substitutions, forward_route, reverse_route
  use util_tally            , only : tally_order, tally_enter, tally_leave, &
       & at_horizon, at_block, at_stage
  use gti_taylor       , only : nodal_coefficient

  implicit none

  private
  public :: chain_block, march_chain, chain_expansion, instant_components
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
     integer               :: first = 0
     integer               :: last  = 0
     integer               :: given = 0
     integer               :: primary = 0

  end type chain_block

  !===================================================================!
  ! What one block contributes to a sensitivity: its jacobian in the
  ! state, its partial in the design, and the part of the
  ! functional's gradient it owns.
  !===================================================================!

  type :: chain_system

     real(dp), allocatable :: rate(:)
     real(dp), allocatable :: g(:)

     ! The stamp of the tangent at the frozen state. Every order of
     ! the expansion, the tangent and the adjoint solve against that
     ! one tangent, and a direct solver factorises it once.
     integer :: mark = 0

  end type chain_system

contains

  !===================================================================!
  ! The components a chain holds at one of the horizon's instants,
  ! found in whichever block computed it. A shared instant is held by
  ! both and reads the same either way, so the earlier is taken.
  !===================================================================!

  pure function instant_components(chain, instant, degrees) result(x)

    type(chain_block), intent(in) :: chain(:)
    integer          , intent(in) :: instant, degrees
    real(dp), allocatable :: x(:)

    integer :: b, local, at

    do b = 1, size(chain)
       if (instant < chain(b) % first .or. instant > chain(b) % last) cycle
       local = instant - chain(b) % first + 1
       at    = chain(b) % instants_at(local)
       x     = chain(b) % state(at + 1:at + degrees)
       return
    end do

    error stop 'gti_chain: that instant lies outside the chain'

  end function instant_components

  !===================================================================!
  ! What a block is given at the instants it shares with the one
  ! before it, laid out the way its own carried rows expect: instant
  ! by instant, degrees within an instant.
  !===================================================================!

  pure function handed_over(earlier, first, given, degrees) result(held)

    type(chain_block), intent(in) :: earlier(:)
    integer          , intent(in) :: first, given, degrees
    real(dp), allocatable :: held(:)

    integer :: i

    allocate(held(given * degrees))

    do i = 1, given
       held((i - 1) * degrees + 1:i * degrees) = &
            & instant_components(earlier, first + i - 1, degrees)
    end do

  end function handed_over

  !===================================================================!
  ! One block's statement, and where its instants sit. A stage family
  ! keeps its instants between its stages; every other keeps one set
  ! per instant.
  !===================================================================!

  subroutine built(scheme, physics, degrees, n, dt, held, rows, instants_at)

    class(family)         , intent(in)  :: scheme
    class(nodal_integrand), intent(in)  :: physics
    integer               , intent(in)  :: degrees, n
    real(dp)              , intent(in)  :: dt(:), held(:)
    type(block_residual)  , intent(out) :: rows
    integer, allocatable  , intent(out) :: instants_at(:)

    integer :: k

    if (marches_by_stages(scheme, degrees)) then
       rows        = stage_block_of(scheme, physics, degrees, n, dt, held)
       instants_at = [(instant_at(k, scheme % num_stages(), degrees), k = 1, n)]
    else
       rows        = block_of(scheme, physics, degrees, n, dt, held)
       instants_at = [((k - 1) * degrees, k = 1, n)]
    end if

  end subroutine built

  !===================================================================!
  ! The whole chain built and marched, block after block, each given
  ! what its predecessor computed at the instants they share.
  !===================================================================!

  subroutine march_chain(schemes, added, physics, degrees, steps, &
       & design, initial, chain, dt, t, achieved, grid_design, left)

    type(family_holder)   , intent(in) :: schemes(:)
    integer               , intent(in) :: added(:), degrees
    class(nodal_integrand), intent(in) :: physics
    real(dp)              , intent(in) :: design, initial(:)
    class(grid)           , intent(in) :: steps
    type(chain_block), allocatable, intent(out) :: chain(:)
    real(dp)         , allocatable, intent(out) :: dt(:), t(:)
    real(dp)              , intent(out) :: achieved
    real(dp), intent(in), optional     :: grid_design(:)
    type(imbalance), intent(out), optional :: left

    type(imbalance) :: one_left
    integer , allocatable :: first(:), last(:)
    real(dp) :: one_achieved
    integer :: b

    if (size(added) < 1) then
       error stop 'gti_chain: a chain holds at least one block'
    end if

    call horizon_bounds(schemes, added, degrees - 1, first, last)
    call partitioned(steps, last(size(added)), dt, t, grid_design)

    allocate(chain(size(added)))
    achieved = 0.0_dp

    call tally_enter(at_horizon)
    do b = 1, size(added)
       call one_block(chain, b, schemes(b) % scheme, physics, degrees, &
            & first(b), last(b), dt, design, initial, one_achieved, one_left)
       achieved = max(achieved, one_achieved)

       ! The report kept is the first block's that did not converge:
       ! every block after it reads a state it never reached.
       if (present(left)) then
          if (b == 1) left = one_left
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

  subroutine one_block(chain, b, scheme, physics, degrees, first, last, dt, &
       & design, initial, achieved, left)

    type(chain_block)     , intent(inout) :: chain(:)
    integer               , intent(in)    :: b, degrees, first, last
    class(family)         , intent(in)    :: scheme
    class(nodal_integrand), intent(in)    :: physics
    real(dp)              , intent(in)    :: dt(:), design, initial(:)
    real(dp)              , intent(out)   :: achieved
    type(imbalance)       , intent(out)   :: left

    real(dp), allocatable :: held(:)

    chain(b) % first   = first
    chain(b) % last    = last
    chain(b) % given   = scheme % history_depth(degrees - 1)
    chain(b) % primary = scheme % primary_degree(degrees - 1)

    if (b == 1) then
       held = initial
    else
       held = handed_over(chain(1:b - 1), first, chain(b) % given, degrees)
    end if

    ! A block whose scheme keeps stages within a step is filed under
    ! the stage level, every other under the block level.
    if (scheme % num_stages() > 1) then
       call tally_enter(at_stage)
    else
       call tally_enter(at_block)
    end if

    call built(scheme, physics, degrees, last - first + 1, dt(first:last), &
         & held, chain(b) % rows, chain(b) % instants_at)

    call swept(chain(b) % rows, design, 1, chain(b) % state, achieved, left)

    call tally_leave()

  end subroutine one_block

  !===================================================================!
  ! The instants one block owns: the ones it computed, and for the
  ! first block the ones it was given as well, since nobody else
  ! holds them.
  !===================================================================!

  pure subroutine owned(chain, b, from, to)

    type(chain_block), intent(in)  :: chain(:)
    integer          , intent(in)  :: b
    integer          , intent(out) :: from, to

    from = chain(b) % first + chain(b) % given
    if (b == 1) from = chain(b) % first
    to = chain(b) % last

  end subroutine owned

  !===================================================================!
  ! The functional and every derivative of it in the design, along a
  ! chain. Each order sweeps the blocks forward, handing its
  ! coefficient over at every junction just as the trajectory does.
  !===================================================================!

  subroutine chain_expansion(chain, physics, integrand, degrees, dt, design, &
       & max_order, f)

    type(chain_block)     , intent(in) :: chain(:)
    class(nodal_integrand), intent(in) :: physics, integrand
    integer               , intent(in) :: degrees, max_order
    real(dp)              , intent(in) :: dt(:), design
    real(dp), allocatable , intent(out) :: f(:)

    ! One design and one functional: what this expansion is built for,
    ! and what the gate is asked about at every order.
    integer, parameter :: num_designs = 1, num_functionals = 1

    type(chain_system), allocatable :: systems(:)
    real(dp), allocatable :: series(:,:,:)
    integer :: b, m, widest, route

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
       call chain_systems(chain, integrand, degrees, dt, design, systems)
    end if

    do m = 1, max_order
       call tally_order(m)
       call tally_enter(at_horizon)

       ! THE GATE. The route is chosen from the counts and not assumed.
       ! Only the forward route is built for an expansion, so a choice
       ! of the other stops the program and says so rather than taking
       ! the dearer one silently.
       route = route_of(num_designs, num_functionals, m)
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

    call chain_functional(chain, integrand, degrees, dt, design, max_order, series, f)

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
    class(nodal_integrand), intent(in)    :: physics
    real(dp)              , intent(in)    :: design
    real(dp)              , intent(inout) :: series(0:, :, :)

    type(stored_directed_graph) :: unknowns
    type(stored_field), allocatable :: inputs(:)
    real(dp), allocatable :: frozen(:,:), coefficient(:), r(:), w(:), held(:)
    integer , allocatable :: at(:)
    integer :: count, carried, p

    count   = chain(b) % rows % num_unknowns()
    carried = chain(b) % rows % num_carried()
    at      = chain(b) % rows % points_at()

    allocate(frozen(0:ubound(series, 1), count))
    frozen = series(:, 1:count, b)
    frozen(m, :) = 0.0_dp

    call nodal_coefficient(physics, degrees, at, frozen, design, m, coefficient)

    allocate(r(count), source=0.0_dp)
    do p = 1, size(at)
       r(at(p) + chain(b) % primary + 1) = coefficient(p)
    end do

    if (b == 1) then
       r(1:carried) = 0.0_dp
    else
       held = coefficients_handed(chain(1:b - 1), chain(b) % first, &
            & chain(b) % given, degrees, series, m)
       r(1:carried) = -held
    end if

    call frozen_at(chain(b), design, unknowns, inputs)
    call solved_linear(chain(b) % rows, unknowns, inputs, -r, .false., systems(b) % mark, 1, w)
    series(m, 1:count, b) = w

  end subroutine one_order

  !===================================================================!
  ! What the earlier blocks found at one order, at the instants a
  ! later one carries. A coefficient travels the junction the way the
  ! trajectory does, and comes from whichever block computed that
  ! instant.
  !===================================================================!

  pure function coefficients_handed(earlier, first, given, degrees, series, order) &
       & result(held)

    type(chain_block), intent(in) :: earlier(:)
    integer          , intent(in) :: first, given, degrees, order
    real(dp)         , intent(in) :: series(0:, :, :)
    real(dp), allocatable :: held(:)

    integer :: i, b, local, at

    allocate(held(given * degrees), source=0.0_dp)

    do i = 1, given
       do b = 1, size(earlier)
          if (first + i - 1 < earlier(b) % first) cycle
          if (first + i - 1 > earlier(b) % last) cycle
          local = first + i - 1 - earlier(b) % first + 1
          at    = earlier(b) % instants_at(local)
          held((i - 1) * degrees + 1:i * degrees) = series(order, at + 1:at + degrees, b)
          exit
       end do
    end do

  end function coefficients_handed

  !===================================================================!
  ! The functional at every order: each block's own instants, each
  ! counted once, weighted by the step that ends at it.
  !===================================================================!

  subroutine chain_functional(chain, integrand, degrees, dt, design, max_order, &
       & series, f)

    type(chain_block)     , intent(in) :: chain(:)
    class(nodal_integrand), intent(in) :: integrand
    integer               , intent(in) :: degrees, max_order
    real(dp)              , intent(in) :: dt(:), design, series(0:, :, :)
    real(dp), allocatable , intent(out) :: f(:)

    real(dp), allocatable :: values(:)
    integer , allocatable :: at(:)
    integer :: b, m, from, to, k, count

    allocate(f(0:max_order), source=0.0_dp)

    do b = 1, size(chain)
       call owned(chain, b, from, to)
       count = chain(b) % rows % num_unknowns()
       at    = [(chain(b) % instants_at(k - chain(b) % first + 1), k = from, to)]

       do m = 0, max_order
          call nodal_coefficient(integrand, degrees, at, series(:, 1:count, b), &
               & design, m, values)
          f(m) = f(m) + sum(dt(from:to) * values)
       end do
    end do

  end subroutine chain_functional

  !===================================================================!
  ! Every block's system at the trajectory already marched. The
  ! gradient is shared out by ownership, each block taking the
  ! instants it computed and none of the instants it was given, so
  ! that a shared instant counts once.
  !===================================================================!

  subroutine chain_systems(chain, integrand, degrees, dt, design, systems)

    type(chain_block)     , intent(in) :: chain(:)
    class(nodal_integrand), intent(in) :: integrand
    integer               , intent(in) :: degrees
    real(dp)              , intent(in) :: dt(:), design
    type(chain_system), allocatable, intent(out) :: systems(:)

    type(stored_directed_graph) :: unknowns
    type(stored_field), allocatable :: inputs(:)
    integer :: b, count

    allocate(systems(size(chain)))

    do b = 1, size(chain)
       count = chain(b) % rows % num_unknowns()
       call frozen_at(chain(b), design, unknowns, inputs)

       systems(b) % mark = fresh_stamp()
       call design_partial(chain(b) % rows, unknowns, inputs, &
            & chain(b) % rows % num_points(), unknowns % vertex_set(), &
            & systems(b) % rate)
       call owned_gradient(chain, b, integrand, degrees, dt, design, inputs, &
            & unknowns, systems(b) % g)
    end do

  end subroutine chain_systems

  !===================================================================!
  ! The functional's gradient over the instants one block owns,
  ! weighted by the step that ends at each. The integrand reads one
  ! instant at a time, so one partial action per degree gives the
  ! whole of it rather than one per unknown.
  !===================================================================!

  subroutine owned_gradient(chain, b, integrand, degrees, dt, design, inputs, &
       & unknowns, g)

    type(chain_block)          , intent(in) :: chain(:)
    integer                    , intent(in) :: b, degrees
    class(nodal_integrand)     , intent(in) :: integrand
    real(dp)                   , intent(in) :: dt(:), design
    type(stored_field)         , intent(in) :: inputs(:)
    type(stored_directed_graph), intent(in) :: unknowns
    real(dp), allocatable      , intent(out) :: g(:)

    type(stored_directed_graph) :: instants
    type(stored_field) :: state, knobs
    real(dp), allocatable :: owned_g(:)
    integer :: count, from, to, k, d, at, held

    count = chain(b) % rows % num_unknowns()
    allocate(g(count), source=0.0_dp)
    call owned(chain, b, from, to)

    held = to - from + 1
    instants = stored_directed_graph(held, tails=[integer ::], heads=[integer ::])

    call at_owned_instants(chain, b, degrees, design, inputs, from, to, instants, &
         & state, knobs)

    ! the gradient over the owned instants, then scattered to where
    ! those instants lie in the block
    call functional_gradient(integrand, instants, [state, knobs], dt(from:to), held, degrees, &
         & instants % vertex_set(), owned_g)

    do d = 0, degrees - 1
       do k = from, to
          at = chain(b) % instants_at(k - chain(b) % first + 1)
          g(at + d + 1) = owned_g((k - from) * degrees + d + 1)
       end do
    end do

    associate (u1 => unknowns); end associate

  end subroutine owned_gradient

  !===================================================================!
  ! The trajectory at the instants one block owns, laid out one
  ! instant at a time so that a nodal rule reads them.
  !===================================================================!

  subroutine at_owned_instants(chain, b, degrees, design, inputs, from, to, &
       & instants, state, knobs)

    type(chain_block)          , intent(in)  :: chain(:)
    integer                    , intent(in)  :: b, degrees, from, to
    real(dp)                   , intent(in)  :: design
    type(stored_field)         , intent(in)  :: inputs(:)
    type(stored_directed_graph), intent(in)  :: instants
    type(stored_field)         , intent(out) :: state, knobs

    real(dp), allocatable :: whole(:), v(:)
    integer :: held, k, at

    held = to - from + 1
    call inputs(1) % real_vector(whole)
    allocate(v(held * degrees))

    do k = from, to
       at = chain(b) % instants_at(k - chain(b) % first + 1)
       v((k - from) * degrees + 1:(k - from + 1) * degrees) = whole(at + 1:at + degrees)
    end do

    state = stored_field('state', instants % vertex_set(), held * degrees)
    knobs = stored_field('design', instants % vertex_set(), held)
    call state % set_real_vector(v)
    call knobs % set_real_vector(spread(design, 1, held))

  end subroutine at_owned_instants

  !===================================================================!
  ! Which block holds one global instant, and where in it.
  !===================================================================!

  pure subroutine holder_of(chain, instant, degrees, held_by, at)

    type(chain_block), intent(in)  :: chain(:)
    integer          , intent(in)  :: instant, degrees
    integer          , intent(out) :: held_by, at

    integer :: b

    held_by = 0
    at      = 0

    do b = 1, size(chain)
       if (instant < chain(b) % first .or. instant > chain(b) % last) cycle
       held_by = b
       at      = chain(b) % instants_at(instant - chain(b) % first + 1)
       return
    end do

    associate (u1 => degrees); end associate

  end subroutine holder_of

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
  ! Forward. Each block solves for its own sensitivity, with the rows
  ! it carried set to what an earlier block already found at those
  ! instants.
  !===================================================================!

  real(dp) function chain_by_tangent(chain, systems, degrees, design) result(df)

    type(chain_block) , intent(in) :: chain(:)
    type(chain_system), intent(in) :: systems(:)
    integer           , intent(in) :: degrees
    real(dp)          , intent(in) :: design

    type(stored_directed_graph) :: unknowns
    type(stored_field), allocatable :: inputs(:)

    real(dp), allocatable :: w(:,:), rhs(:), one(:)
    integer :: b, widest, i, k, d, held_by, at

    widest = 0
    do b = 1, size(chain)
       widest = max(widest, chain(b) % rows % num_unknowns())
    end do

    allocate(w(widest, size(chain)), source=0.0_dp)
    df = 0.0_dp

    do b = 1, size(chain)
       rhs = -systems(b) % rate

       do i = 1, chain(b) % given * degrees
          k = chain(b) % first + (i - 1) / degrees
          d = mod(i - 1, degrees)
          call holder_of(chain(1:b - 1), k, degrees, held_by, at)
          if (held_by > 0) rhs(i) = w(at + d + 1, held_by)
       end do

       call frozen_at(chain(b), design, unknowns, inputs)
       call solved_linear(chain(b) % rows, unknowns, inputs, rhs, .false., systems(b) % mark, 1, one)
       w(1:size(one), b) = one
       df = df + dot_product(systems(b) % g, one)
    end do

  end function chain_by_tangent

  !===================================================================!
  ! Backward. Each block solves against its own transpose, with its
  ! own share of the gradient plus whatever a later block's costate
  ! left on the instants they share.
  !===================================================================!

  real(dp) function chain_by_adjoint(chain, systems, degrees, design) result(df)

    type(chain_block) , intent(in) :: chain(:)
    type(chain_system), intent(in) :: systems(:)
    integer           , intent(in) :: degrees
    real(dp)          , intent(in) :: design

    type(stored_directed_graph) :: unknowns
    type(stored_field), allocatable :: inputs(:)

    real(dp), allocatable :: rhs(:,:), lambda(:)
    integer :: b, c, widest, i, k, d, held_by, at

    widest = 0
    do b = 1, size(chain)
       widest = max(widest, chain(b) % rows % num_unknowns())
    end do

    allocate(rhs(widest, size(chain)), source=0.0_dp)
    do b = 1, size(chain)
       rhs(1:size(systems(b) % g), b) = systems(b) % g
    end do

    df = 0.0_dp

    do b = size(chain), 1, -1
       call frozen_at(chain(b), design, unknowns, inputs)
       call solved_linear(chain(b) % rows, unknowns, inputs, &
            & rhs(1:chain(b) % rows % num_unknowns(), b), .true., systems(b) % mark, 1, lambda)
       df = df - dot_product(lambda, systems(b) % rate)

       do i = 1, chain(b) % given * degrees
          k = chain(b) % first + (i - 1) / degrees
          d = mod(i - 1, degrees)
          call holder_of(chain(1:b - 1), k, degrees, held_by, at)
          if (held_by > 0) rhs(at + d + 1, held_by) = rhs(at + d + 1, held_by) + lambda(i)
       end do

       associate (u1 => c); end associate
    end do

  end function chain_by_adjoint

end module gti_chain
