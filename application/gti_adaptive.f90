!=====================================================================!
! Adaptive time stepping, phase one: an error-controlled forward march
! that discovers a grid. The grid it returns is an ordinary partition
! of the duration; the chain runs on it unchanged, so its sensitivities
! are formed exactly as on a fixed grid. The march itself takes no
! derivatives.
!
!             THE ERROR ESTIMATE
!
! Step doubling: one step of size h against two of h/2, both from the
! same state with the same scheme. Their difference on the solution
! components below the highest estimates the local error of the h step,
! which for an order-p scheme is O(h^(p+1)). A step is accepted when the
! estimate, measured relative to the state or absolute as a solve's
! tolerance is, is at or below the tolerance; the next step is
! h (tolerance / estimate)^(1/(p+1)), bounded so one step neither grows
! nor shrinks without limit, and clamped so the last lands on the
! duration exactly. The state advanced is the two-half-step one.
!
! The families marched are the diagonally implicit ones, self-starting
! and one step at a time, so a step's coefficients do not depend on the
! steps around it. A variable-step multistep family, whose coefficients
! do, is a separate construction and is refused here by its history
! reaching past one instant.
!
!             WHAT IS REFUSED
!
! A scheme that reaches back over more than one instant; a nonpositive
! duration, tolerance or first step; a step that stays above the
! tolerance past fifty attempts, where the estimate is not falling.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_adaptive

  use util_precision   , only : dp
  use operation_family , only : family
  use operation_grid   , only : uniform_grid
  use operation_expression, only : expression
  use gti_block        , only : block_residual
  use gti_expansion    , only : expansion, family_holder
  use gti_march        , only : consistent_state, block_from, solved

  implicit none

  private
  public :: adaptive_partition

contains

  !-------------------------------------------------------------------!
  ! One step of size h from a state by n - 1 uniform steps: n = 2 is
  ! the single step, n = 3 the two half steps. A block over n instants
  ! of a uniform grid of extent h has n - 1 steps of h / (n - 1); the
  ! arriving instant's components come back.
  !-------------------------------------------------------------------!

  subroutine stepped(scheme, physics, degrees, state, h, n, design, arrived)

    class(family)   , intent(in)  :: scheme
    type(expression), intent(in)  :: physics
    integer         , intent(in)  :: degrees, n
    real(dp)        , intent(in)  :: state(:), h, design
    real(dp), allocatable, intent(out) :: arrived(:)

    type(expansion)      :: tower
    type(family_holder)  :: holder(1)
    type(block_residual) :: rows
    integer, allocatable :: at(:)
    real(dp), allocatable :: q(:)
    real(dp) :: achieved
    integer  :: last

    allocate(holder(1) % scheme, source=scheme)
    call tower % build(physics, holder, [n], uniform_grid(h), 0, design)
    call block_from(tower, 1, scheme, physics, state, rows, at)
    call solved(rows, design, q, achieved)

    last    = at(size(at))
    arrived = q(last + 1:last + degrees)

  end subroutine stepped

  !-------------------------------------------------------------------!
  ! The estimate over the solution components below the highest, the
  ! highest being algebraically determined: relative to the state, or
  ! absolute.
  !-------------------------------------------------------------------!

  pure real(dp) function estimate(coarse, fine, degrees, relative) result(e)

    real(dp), intent(in) :: coarse(:), fine(:)
    integer , intent(in) :: degrees
    logical , intent(in) :: relative

    e = norm2(coarse(1:degrees - 1) - fine(1:degrees - 1))
    if (relative) e = e / max(norm2(fine(1:degrees - 1)), tiny(1.0_dp))

  end function estimate

  !-------------------------------------------------------------------!
  ! The accepted steps of a scheme of order p over [0, duration] to a
  ! tolerance, from the state consistent with the lower components. The
  ! first step is a fraction of the duration; the controller is bounded
  ! and the last step clamped to the duration.
  !-------------------------------------------------------------------!

  function adaptive_partition(scheme, p, physics, degrees, duration, lower, design, &
       & tolerance, relative, rejects) result(dt)

    class(family)   , intent(in)  :: scheme
    integer         , intent(in)  :: p, degrees
    type(expression), intent(in)  :: physics
    real(dp)        , intent(in)  :: duration, lower(:), design, tolerance
    logical         , intent(in)  :: relative
    integer         , intent(out), optional :: rejects
    real(dp), allocatable :: dt(:)

    real(dp), parameter :: safety = 0.9_dp, grow = 5.0_dp, shrink = 0.2_dp
    real(dp), allocatable :: state(:), coarse(:), fine(:)
    real(dp) :: t, h, e, factor
    integer  :: attempt, rejected

    if (scheme % history_depth(degrees - 1) > 1) then
       error stop 'gti_adaptive: an adaptive march is a self-starting scheme'
    end if
    if (duration <= 0.0_dp .or. tolerance <= 0.0_dp) then
       error stop 'gti_adaptive: the duration and the tolerance are positive'
    end if

    state    = consistent_state(physics, degrees, lower, design)
    dt       = [real(dp) ::]
    t        = 0.0_dp
    h        = duration / 8.0_dp
    rejected = 0

    do while (t < duration * (1.0_dp - 1.0e-12_dp))
       h = min(h, duration - t)
       attempt = 0
       do
          attempt = attempt + 1
          call stepped(scheme, physics, degrees, state, h, 2, design, coarse)
          call stepped(scheme, physics, degrees, state, h, 3, design, fine)
          e      = estimate(coarse, fine, degrees, relative)
          factor = safety * (tolerance / max(e, tiny(1.0_dp))) ** (1.0_dp / real(p + 1, dp))
          factor = min(grow, max(shrink, factor))
          if (e <= tolerance .or. h <= duration * 1.0e-10_dp) exit
          rejected = rejected + 1
          h = h * factor
          if (attempt > 50) then
             error stop 'gti_adaptive: a step stays above the tolerance past fifty attempts'
          end if
       end do
       dt    = [dt, h]
       t     = t + h
       state = fine
       h     = h * factor
    end do

    if (present(rejects)) rejects = rejected

  end function adaptive_partition

end module gti_adaptive
