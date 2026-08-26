!=====================================================================!
! One block of a stage family, as a statement to be driven to zero.
!
! A stage block does not hold one set of components per instant. Each
! step carries s stages and then the instant it arrives at, and the
! block's first instant carries nothing but itself, because no step
! arrives there. So the unknowns run
!
!      slice 1        the instant alone
!      slice k > 1    stage 1 .. stage s, then the arriving instant
!
! and the block is square in them.
!
!             WHICH ROW IS WHICH
!
! Every degree below the highest is determined by the tableau: a
! stage reads the incoming instant at its own degree and the stages
! at or before it one degree up, and the arriving instant reads the
! incoming instant at its own degree and every stage one degree up.
! The highest degree is determined at a stage by the physics and at
! the arriving instant by the recovery, which reads every stage at
! that same degree.
!
! So the physics is evaluated at the stages and nowhere else, and the
! instants are recovered from them. That is what a Runge-Kutta method
! is, and it is why the block residual takes its evaluation points
! rather than assuming one per instant.
!
!             THE WEIGHTS
!
! From the family, never from this module. The family numbers a step
! as the incoming instant, then its stages, then the arriving
! instant, and that numbering is mapped onto the block's unknowns
! here. The step is read at every vertex of it, the stages included,
! because that is where the scaling reads it.
!
!             WHAT IS REFUSED
!
! A block of fewer than two slices, since it would hold no step, and
! a family that reaches further back than one instant, which no
! stage family does.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_stage

  use util_precision  , only : dp
  use field_calculus          , only : field
  use operation_stencil       , only : stencil
  use operation_family        , only : family
  use operation_coupling      , only : weights_of
  use operation_weight        , only : scheme_weight
  use operation_scheme_stencil, only : derived_constraints
  use physics_integrand       , only : nodal_integrand
  use gti_block               , only : block_residual

  implicit none

  private
  public :: stage_block_of, stage_unknowns, stage_points, instant_at

contains

  pure integer function slice_base(kk, s, nd) result(at)

    integer, intent(in) :: kk, s, nd

    if (kk == 1) then
       at = 0
    else
       at = nd + (kk - 2) * (s + 1) * nd
    end if

  end function slice_base

  pure integer function instant_at(kk, s, nd) result(at)

    integer, intent(in) :: kk, s, nd

    at = slice_base(kk, s, nd)
    if (kk > 1) at = at + s * nd

  end function instant_at

  pure integer function stage_at(kk, i, s, nd) result(at)

    integer, intent(in) :: kk, i, s, nd

    at = slice_base(kk, s, nd) + (i - 1) * nd

  end function stage_at

  pure integer function stage_unknowns(n, s, nd) result(count)

    integer, intent(in) :: n, s, nd

    count = nd + (n - 1) * (s + 1) * nd

  end function stage_unknowns

  !===================================================================!
  ! The evaluation points: every stage of every step, in order. The
  ! instants are not among them; they are recovered.
  !===================================================================!

  pure function stage_points(n, s, nd) result(at)

    integer, intent(in) :: n, s, nd
    integer, allocatable :: at(:)

    integer :: kk, i, e

    allocate(at((n - 1) * s))
    e = 0

    do kk = 2, n
       do i = 1, s
          e = e + 1
          at(e) = stage_at(kk, i, s, nd)
       end do
    end do

  end function stage_points

  !===================================================================!
  ! The rows of one step, in the family's own numbering: vertex one
  ! is the incoming instant, then the stages, then the arriving
  ! instant. Counted on the first pass and filled on the second.
  !===================================================================!

  subroutine step_reach(nd, s, tails, heads, source_degree, determines)

    integer, intent(in) :: nd, s
    integer, allocatable, intent(out) :: tails(:), heads(:)
    integer, allocatable, intent(out) :: source_degree(:), determines(:)

    integer :: d, i, j, e

    e = (nd - 1) * (2 * s + s * (s + 1) / 2 + 1) + s
    allocate(tails(e), heads(e), source_degree(e), determines(e))
    e = 0

    do d = 0, nd - 2
       do i = 1, s
          call put(e, 1, 1 + i, d, d)
          do j = 1, i
             call put(e, 1 + j, 1 + i, d + 1, d)
          end do
       end do
       call put(e, 1, 2 + s, d, d)
       do j = 1, s
          call put(e, 1 + j, 2 + s, d + 1, d)
       end do
    end do

    do j = 1, s
       call put(e, 1 + j, 2 + s, nd - 1, nd - 1)
    end do

  contains

    subroutine put(e, tail, head, sd, det)

      integer, intent(in)    :: tail, head, sd, det
      integer, intent(inout) :: e

      e = e + 1

      tails(e)         = tail
      heads(e)         = head
      source_degree(e) = sd
      determines(e)    = det

    end subroutine put

  end subroutine step_reach

  !===================================================================!
  ! A vertex of one step, as an unknown of the block.
  !===================================================================!

  pure integer function unknown_of(vertex, degree, kk, s, nd) result(at)

    integer, intent(in) :: vertex, degree, kk, s, nd

    if (vertex == 1) then
       at = instant_at(kk - 1, s, nd) + degree + 1
    else if (vertex == 2 + s) then
       at = instant_at(kk, s, nd) + degree + 1
    else
       at = stage_at(kk, vertex - 1, s, nd) + degree + 1
    end if

  end function unknown_of

  !===================================================================!
  ! The weights of one step, from the family.
  !===================================================================!

  subroutine step_weights(scheme, s, tails, heads, source_degree, determines, step, w)

    class(family), intent(in) :: scheme
    integer      , intent(in) :: s, tails(:), heads(:), source_degree(:), determines(:)
    real(dp)     , intent(in) :: step
    real(dp), allocatable, intent(out) :: w(:)

    call weights_of(scheme_weight(scheme), s + 2, tails, heads, spread(step, 1, s + 2), &
         & source_degree, determines, w)

  end subroutine step_weights

  !===================================================================!
  ! The derived rows of the whole block: one step's pattern repeated,
  ! its weights taken at each step's own size, and its vertices
  ! mapped onto the block's unknowns.
  !===================================================================!

  function stage_rows(scheme, nd, n, dt) result(rows)

    class(family), intent(in) :: scheme
    integer      , intent(in) :: nd, n
    real(dp)     , intent(in) :: dt(:)
    type(stencil) :: rows

    integer , allocatable :: tails(:), heads(:), source_degree(:), determines(:)
    integer , allocatable :: into(:), from(:)
    real(dp), allocatable :: w(:), weight(:)
    integer :: s, kk, e, at, per

    s = scheme % num_stages()
    call step_reach(nd, s, tails, heads, source_degree, determines)
    per = size(tails)

    allocate(into((n - 1) * per), from((n - 1) * per), weight((n - 1) * per))
    at = 0

    do kk = 2, n
       call step_weights(scheme, s, tails, heads, source_degree, determines, dt(kk), w)
       do e = 1, per
          at = at + 1
          into(at)   = unknown_of(heads(e), determines(e), kk, s, nd)
          from(at)   = unknown_of(tails(e), source_degree(e), kk, s, nd)
          weight(at) = w(e)
       end do
    end do

    rows = derived_constraints(into, from, weight, stage_unknowns(n, s, nd), &
         & 'stage rows')

  end function stage_rows

  !===================================================================!
  ! The whole statement of one stage block. The first instant is
  ! carried; every stage and every later instant is solved for.
  !===================================================================!

  function stage_block_of(scheme, physics, nd, n, dt, held) result(rows)

    class(family)         , intent(in) :: scheme
    class(nodal_integrand), intent(in) :: physics
    integer               , intent(in) :: nd, n
    real(dp)              , intent(in) :: dt(:), held(:)
    type(block_residual) :: rows

    integer :: s, d

    if (n < 2) then
       error stop 'gti_stage: a stage block holds at least one step'
    end if
    if (scheme % history_depth(nd - 1) /= 1) then
       error stop 'gti_stage: a stage family reaches one instant back'
    end if
    if (size(held) /= nd) then
       error stop 'gti_stage: the first instant is carried, one value per degree'
    end if

    s = scheme % num_stages()

    rows = block_residual(stage_rows(scheme, nd, n, dt), physics, &
         & stage_points(n, s, nd), stage_unknowns(n, s, nd), nd, nd - 1, &
         & [(d, d = 1, nd)], held)

    ! where every unknown lies: the first instant is its own member,
    ! and a step's stages and the instant it arrives at are one member
    call rows % placed_in(step_labels(n, s, nd), spread(1, 1, stage_unknowns(n, s, nd)))

  end function stage_block_of

  pure function step_labels(n, s, nd) result(label)

    integer, intent(in) :: n, s, nd
    integer :: label(stage_unknowns(n, s, nd))

    integer :: kk, i, d

    label(1:nd) = 1
    do kk = 2, n
       do i = 1, s
          do d = 1, nd
             label(stage_at(kk, i, s, nd) + d) = kk
          end do
       end do
       do d = 1, nd
          label(instant_at(kk, s, nd) + d) = kk
       end do
    end do

  end function step_labels

end module gti_stage
