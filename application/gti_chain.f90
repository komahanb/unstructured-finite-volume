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
! unknowns its k-th instant sits. A block hands its successor components, never
! rows, so that question is the whole of the interface between them.
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

  use iso_fortran_env  , only : dp => REAL64
  use operation_family , only : family
  use operation_grid   , only : grid
  use physics_integrand, only : nodal_integrand
  use gti_expansion    , only : family_holder, marches_by_stages
  use gti_block        , only : block_residual
  use gti_march        , only : partitioned, horizon_bounds, block_of, solved
  use gti_stage        , only : stage_block_of, instant_at
  use view_directed_stored, only : stored_directed_graph
  use field_stored     , only : stored_field
  use gti_sweeps       , only : dense_solve, jacobian_of
  use gti_taylor       , only : nodal_coefficient

  implicit none

  private
  public :: chain_block, march_chain, chain_expansion

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

contains

  !===================================================================!
  ! What a block is given at the instants it shares with the one
  ! before it, laid out the way its own carried rows expect: instant
  ! by instant, degrees within an instant.
  !===================================================================!

  pure function handed_over(previous, first, given, degrees, series) result(held)

    type(chain_block), intent(in) :: previous
    integer          , intent(in) :: first, given, degrees
    real(dp)         , intent(in) :: series(:)
    real(dp), allocatable :: held(:)

    integer :: i, at

    allocate(held(given * degrees))

    do i = 1, given
       at = previous % instants_at(first + i - 1 - previous % first + 1)
       held((i - 1) * degrees + 1:i * degrees) = series(at + 1:at + degrees)
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
       & design, initial, chain, dt, t, achieved, grid_design)

    type(family_holder)   , intent(in) :: schemes(:)
    integer               , intent(in) :: added(:), degrees
    class(nodal_integrand), intent(in) :: physics
    real(dp)              , intent(in) :: design, initial(:)
    class(grid)           , intent(in) :: steps
    type(chain_block), allocatable, intent(out) :: chain(:)
    real(dp)         , allocatable, intent(out) :: dt(:), t(:)
    real(dp)              , intent(out) :: achieved
    real(dp), intent(in), optional     :: grid_design(:)

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

    do b = 1, size(added)
       call one_block(chain, b, schemes(b) % scheme, physics, degrees, &
            & first(b), last(b), dt, design, initial, one_achieved)
       achieved = max(achieved, one_achieved)
    end do

  end subroutine march_chain

  !===================================================================!
  ! One block of a chain: given what its predecessor computed at the
  ! instants they share, or the initial conditions if it is first,
  ! then built and solved.
  !===================================================================!

  subroutine one_block(chain, b, scheme, physics, degrees, first, last, dt, &
       & design, initial, achieved)

    type(chain_block)     , intent(inout) :: chain(:)
    integer               , intent(in)    :: b, degrees, first, last
    class(family)         , intent(in)    :: scheme
    class(nodal_integrand), intent(in)    :: physics
    real(dp)              , intent(in)    :: dt(:), design, initial(:)
    real(dp)              , intent(out)   :: achieved

    real(dp), allocatable :: held(:)

    chain(b) % first   = first
    chain(b) % last    = last
    chain(b) % given   = scheme % history_depth(degrees - 1)
    chain(b) % primary = scheme % primary_degree(degrees - 1)

    if (b == 1) then
       held = initial
    else
       held = handed_over(chain(b - 1), first, chain(b) % given, degrees, &
            & chain(b - 1) % state)
    end if

    call built(scheme, physics, degrees, last - first + 1, dt(first:last), &
         & held, chain(b) % rows, chain(b) % instants_at)

    call solved(chain(b) % rows, design, chain(b) % state, achieved)

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

    real(dp), allocatable :: series(:,:,:), a(:,:), jacobians(:,:,:)
    integer :: b, m, widest

    widest = 0
    do b = 1, size(chain)
       widest = max(widest, chain(b) % rows % num_unknowns())
    end do

    allocate(series(0:max_order, widest, size(chain)), source=0.0_dp)
    allocate(jacobians(widest, widest, size(chain)), source=0.0_dp)

    do b = 1, size(chain)
       call block_jacobian(chain(b), degrees, design, a)
       jacobians(1:size(a, 1), 1:size(a, 2), b) = a
       series(0, 1:size(chain(b) % state), b) = chain(b) % state
    end do

    do m = 1, max_order
       do b = 1, size(chain)
          call one_order(chain, b, physics, degrees, design, m, jacobians, series)
       end do
    end do

    call chain_functional(chain, integrand, degrees, dt, design, max_order, series, f)

  end subroutine chain_expansion

  !===================================================================!
  ! One block's jacobian at the trajectory it marched.
  !===================================================================!

  subroutine block_jacobian(b, degrees, design, a)

    type(chain_block), intent(in) :: b
    integer          , intent(in) :: degrees
    real(dp)         , intent(in) :: design
    real(dp), allocatable, intent(out) :: a(:,:)

    type(stored_directed_graph) :: unknowns
    type(stored_field) :: state, knobs
    integer :: count

    associate (u1 => degrees); end associate

    count    = b % rows % num_unknowns()
    unknowns = stored_directed_graph(count, tails=[integer ::], heads=[integer ::])
    state    = stored_field('state', unknowns % vertex_set(), count)
    knobs    = stored_field('design', unknowns % vertex_set(), b % rows % num_points())
    call state % set_real_vector(b % state)
    call knobs % set_real_vector(spread(design, 1, b % rows % num_points()))

    call jacobian_of(b % rows, unknowns, [state, knobs], count, &
         & unknowns % vertex_set(), a)

  end subroutine block_jacobian

  !===================================================================!
  ! One block at one order: its own physics on the rows it governs,
  ! and its predecessor's coefficient on the rows it was given.
  !===================================================================!

  subroutine one_order(chain, b, physics, degrees, design, m, jacobians, series)

    type(chain_block)     , intent(in)    :: chain(:)
    integer               , intent(in)    :: b, degrees, m
    class(nodal_integrand), intent(in)    :: physics
    real(dp)              , intent(in)    :: design, jacobians(:,:,:)
    real(dp)              , intent(inout) :: series(0:, :, :)

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
       held = handed_over(chain(b - 1), chain(b) % first, chain(b) % given, degrees, &
            & series(m, 1:chain(b - 1) % rows % num_unknowns(), b - 1))
       r(1:carried) = -held
    end if

    call dense_solve(jacobians(1:count, 1:count, b), -r, .false., w)
    series(m, 1:count, b) = w

  end subroutine one_order

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

end module gti_chain
