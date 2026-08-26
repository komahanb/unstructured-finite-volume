! Adaptive time stepping, two-phase: an error-controlled forward march
! discovers the grid, then the existing chain runs on that grid
! unchanged and its sensitivities pass the same cross-checks as on a
! fixed grid.
!
!             THE ERROR ESTIMATE
!
! Step doubling: one step of size h against two of h/2, both from the
! same state with the same scheme. Their difference on the solution
! components estimates the local error of the h step, which for an
! order-p scheme is O(h^(p+1)). A step is accepted when the estimate,
! measured the way a solve's tolerance is - relative to the state or
! absolute - is at or below the tolerance; the next step is
! h (tolerance / estimate)^(1/(p+1)), bounded so one step neither
! grows nor shrinks without limit, and clamped so the last lands on
! the duration exactly. The state advanced is the two-half-step one,
! the more accurate of the pair.
!
! This is scheme-agnostic in principle; it is exercised here on the
! diagonally implicit families, which are self-starting and take one
! step at a time, so a step's coefficients do not depend on the steps
! around it. A variable-step multistep family, whose coefficients do,
! is a separate construction.
!
!             THE VERIFICATION
!
! The grid covers the duration to round-off; a tighter tolerance makes
! more steps and a functional nearer the refined fixed-grid value; and
! the discovered grid, handed to the existing expansion, yields a
! functional and its design derivatives whose forward and reverse
! routes agree - the adaptive grid is an ordinary grid to everything
! above the march.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program adaptive_grid

  use util_precision       , only : dp
  use operation_family     , only : family
  use operation_family_dirk, only : implicit_midpoint, crouzeix_two_stage, crouzeix_three_stage
  use operation_grid       , only : designed_grid, uniform_grid
  use gti_block            , only : block_residual
  use gti_expansion        , only : expansion, family_holder
  use gti_march            , only : consistent_state, block_from, solved
  use gti_chain            , only : chain_block, march_chain, chain_expansion, &
       & chain_derivative, chain_system, chain_systems, functional_holder, one_functional
  use gti_sweeps           , only : forward_route, reverse_route
  use physics_vanderpol    , only : van_der_pol, van_der_pol_energy

  implicit none

  integer , parameter :: degrees  = 3          ! van der Pol is degree two
  real(dp), parameter :: duration = 4.0_dp
  real(dp), parameter :: design   = 1.0_dp
  real(dp), parameter :: q0 = 2.0_dp, qd0 = 0.0_dp

  call report('implicit midpoint (order 2)', 2)
  call report('crouzeix two stage (order 3)', 3)
  call report('crouzeix three stage (order 4)', 4)

contains

  function dirk_of(order) result(scheme)

    integer, intent(in) :: order
    class(family), allocatable :: scheme

    select case (order)
    case (2); allocate(scheme, source=implicit_midpoint())
    case (3); allocate(scheme, source=crouzeix_two_stage())
    case (4); allocate(scheme, source=crouzeix_three_stage())
    case default; error stop 'adaptive_grid: order two, three or four'
    end select

  end function dirk_of

  !-------------------------------------------------------------------!
  ! One step of size h from a state, and the same by two of h/2. The
  ! arriving instant's components come back from each; the pair drives
  ! the estimate. A block over n instants of a uniform grid of extent
  ! h has n - 1 steps of h / (n - 1), so n = 2 is one step and n = 3
  ! two half steps.
  !-------------------------------------------------------------------!

  subroutine stepped(scheme, state, h, n, arrived)

    class(family), intent(in)  :: scheme
    real(dp)     , intent(in)  :: state(:), h
    integer      , intent(in)  :: n
    real(dp), allocatable, intent(out) :: arrived(:)

    type(expansion)      :: tower
    type(family_holder)  :: holder(1)
    type(block_residual) :: rows
    integer, allocatable :: at(:)
    real(dp), allocatable :: q(:)
    real(dp) :: achieved
    integer  :: last

    allocate(holder(1) % scheme, source=scheme)
    call tower % build(van_der_pol(degrees - 1), holder, [n], uniform_grid(h), 0, design)
    call block_from(tower, 1, scheme, van_der_pol(degrees - 1), state, rows, at)
    call solved(rows, design, q, achieved)

    last    = at(size(at))
    arrived = q(last + 1:last + degrees)

  end subroutine stepped

  !-------------------------------------------------------------------!
  ! The estimate the way a solve reads its tolerance: over the
  ! solution components below the highest, relative to the state.
  !-------------------------------------------------------------------!

  pure real(dp) function estimate(coarse, fine) result(e)

    real(dp), intent(in) :: coarse(:), fine(:)

    e = norm2(coarse(1:degrees - 1) - fine(1:degrees - 1)) &
         & / max(norm2(fine(1:degrees - 1)), tiny(1.0_dp))

  end function estimate

  !-------------------------------------------------------------------!
  ! The adaptive march: the accepted steps of a scheme of order p over
  ! [0, duration] to a tolerance. The controller is bounded and the
  ! last step is clamped to land on the duration.
  !-------------------------------------------------------------------!

  subroutine adaptive_partition(scheme, p, tol, dt, rejects)

    class(family), intent(in)  :: scheme
    integer      , intent(in)  :: p
    real(dp)     , intent(in)  :: tol
    real(dp), allocatable, intent(out) :: dt(:)
    integer      , intent(out) :: rejects

    real(dp), parameter :: safety = 0.9_dp, grow = 5.0_dp, shrink = 0.2_dp
    real(dp), allocatable :: state(:), coarse(:), fine(:)
    real(dp) :: t, h, e, factor
    integer  :: attempt

    state = consistent_state(van_der_pol(degrees - 1), degrees, [q0, qd0], design)
    dt    = [real(dp) ::]
    t     = 0.0_dp
    h     = duration / 8.0_dp
    rejects = 0

    do while (t < duration * (1.0_dp - 1.0e-12_dp))
       h = min(h, duration - t)
       attempt = 0
       do
          attempt = attempt + 1
          call stepped(scheme, state, h, 2, coarse)
          call stepped(scheme, state, h, 3, fine)
          e      = estimate(coarse, fine)
          factor = safety * (tol / max(e, tiny(1.0_dp))) ** (1.0_dp / real(p + 1, dp))
          factor = min(grow, max(shrink, factor))
          if (e <= tol .or. h <= duration * 1.0e-10_dp) exit
          rejects = rejects + 1
          h = h * factor
          if (attempt > 50) error stop 'adaptive_grid: a step is refused past fifty attempts'
       end do
       dt    = [dt, h]
       t     = t + h
       state = fine
       h     = h * factor
    end do

  end subroutine adaptive_partition

  !-------------------------------------------------------------------!
  ! The functional and its first design derivative on a given grid, by
  ! the forward expansion and by the reverse route, so the two can be
  ! compared: the adaptive grid is an ordinary designed grid here.
  !-------------------------------------------------------------------!

  subroutine on_grid(scheme, dt, f, forward, reverse)

    class(family), intent(in)  :: scheme
    real(dp)     , intent(in)  :: dt(:)
    real(dp)     , intent(out) :: f, forward, reverse

    type(family_holder)     :: schemes(1)
    type(functional_holder) :: functionals(1)
    type(chain_block), allocatable :: chain(:)
    type(chain_system), allocatable :: systems(:)
    type(expansion), allocatable :: tower
    real(dp), allocatable :: grid_dt(:), t(:), fvals(:,:), df(:,:), other(:,:)
    real(dp) :: achieved
    integer  :: n

    n = size(dt) + 1
    allocate(schemes(1) % scheme, source=scheme)
    functionals(1) = one_functional(van_der_pol_energy(degrees - 1))

    call march_chain(schemes, [n - 1], van_der_pol(degrees - 1), degrees, &
         & designed_grid(duration), design, &
         & consistent_state(van_der_pol(degrees - 1), degrees, [q0, qd0], design), &
         & chain, tower, grid_dt, t, achieved, grid_design=dt)

    call chain_expansion(chain, tower, functionals, degrees, 1, fvals)
    f = fvals(0, 1)

    call chain_systems(chain, tower, functionals, degrees, systems)
    call chain_derivative(chain, tower, systems, functionals, degrees, 1, forward_route, df)
    call chain_derivative(chain, tower, systems, functionals, degrees, 1, reverse_route, other)
    forward = df(1, 1)
    reverse = other(1, 1)

  end subroutine on_grid

  subroutine report(title, order)

    character(len=*), intent(in) :: title
    integer         , intent(in) :: order

    class(family), allocatable :: scheme
    real(dp), allocatable :: dt(:)
    real(dp) :: tol, f, forward, reverse, span
    integer  :: rejects, level

    scheme = dirk_of(order)

    write(*,'(a)') ' '
    write(*,'(a)') ' ' // title
    write(*,'(a)') '   tolerance     steps   rejects        sum dt - T          functional     forward-reverse'
    do level = 1, 4
       tol = 10.0_dp ** (-3 - level)
       call adaptive_partition(scheme, order, tol, dt, rejects)
       call on_grid(scheme, dt, f, forward, reverse)
       span = sum(dt) - duration
       write(*,'(a,es9.1,i9,i9,es18.2,f18.9,es18.2)') '   ', tol, size(dt), rejects, span, f, &
            & abs(forward - reverse) / max(1.0_dp, abs(forward))
    end do

  end subroutine report

end program adaptive_grid
