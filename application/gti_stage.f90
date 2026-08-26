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
! and the block is square in them. Over a field every stage and every
! instant holds one set of components per node, node by node, so each
! of those is width = nodes x degrees wide; a single node's block is
! the case nodes = 1.
!
!             WHICH ROW IS WHICH
!
! Every degree below the highest is determined by the tableau: a
! stage reads the incoming instant at its own degree and the stages
! at or before it one degree up, and the arriving instant reads the
! incoming instant at its own degree and every stage one degree up.
! The highest degree is determined at a stage by the physics and at
! the arriving instant by the recovery, which reads every stage at
! that same degree. Each node's history is its own, so the rows are
! the same at every node; the level below, laid on the stages, is
! what couples the nodes.
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
  use operation_coupling      , only : weights_of, weights_varied
  use operation_weight        , only : scheme_weight
  use operation_scheme_stencil, only : derived_constraints
  use physics_integrand       , only : nodal_integrand
  use gti_block               , only : block_residual

  implicit none

  private
  public :: stage_block_of, stage_unknowns, stage_points, instant_at, stage_rows

contains

  pure integer function slice_base(kk, s, width) result(at)

    integer, intent(in) :: kk, s, width

    if (kk == 1) then
       at = 0
    else
       at = width + (kk - 2) * (s + 1) * width
    end if

  end function slice_base

  pure integer function instant_at(kk, s, width) result(at)

    integer, intent(in) :: kk, s, width

    at = slice_base(kk, s, width)
    if (kk > 1) at = at + s * width

  end function instant_at

  pure integer function stage_at(kk, i, s, width) result(at)

    integer, intent(in) :: kk, i, s, width

    at = slice_base(kk, s, width) + (i - 1) * width

  end function stage_at

  pure integer function stage_unknowns(n, s, width) result(count)

    integer, intent(in) :: n, s, width

    count = width + (n - 1) * (s + 1) * width

  end function stage_unknowns

  !===================================================================!
  ! The evaluation points: every node at every stage of every step,
  ! in order. The instants are not among them; they are recovered.
  !===================================================================!

  pure function stage_points(n, s, nd, nodes) result(at)

    integer, intent(in) :: n, s, nd, nodes
    integer, allocatable :: at(:)

    integer :: kk, i, p, e

    allocate(at((n - 1) * s * nodes))
    e = 0

    do kk = 2, n
       do i = 1, s
          do p = 1, nodes
             e = e + 1
             at(e) = stage_at(kk, i, s, nd * nodes) + (p - 1) * nd
          end do
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
  ! A vertex of one step at one node, as an unknown of the block.
  !===================================================================!

  pure integer function unknown_of(vertex, degree, kk, s, nd, node, nodes) result(at)

    integer, intent(in) :: vertex, degree, kk, s, nd, node, nodes

    integer :: width

    width = nd * nodes
    if (vertex == 1) then
       at = instant_at(kk - 1, s, width)
    else if (vertex == 2 + s) then
       at = instant_at(kk, s, width)
    else
       at = stage_at(kk, vertex - 1, s, width)
    end if
    at = at + (node - 1) * nd + degree + 1

  end function unknown_of

  !===================================================================!
  ! The weights of one step, from the family.
  !===================================================================!

  subroutine step_weights(scheme, s, tails, heads, source_degree, determines, step, w, along, along2)

    class(family), intent(in) :: scheme
    integer      , intent(in) :: s, tails(:), heads(:), source_degree(:), determines(:)
    real(dp)     , intent(in) :: step
    real(dp), allocatable, intent(out) :: w(:)
    real(dp)     , intent(in), optional :: along, along2

    if (present(along2)) then
       call weights_varied(scheme_weight(scheme), s + 2, tails, heads, spread(step, 1, s + 2), &
            & spread(along, 1, s + 2), source_degree, determines, w, spread(along2, 1, s + 2))
    else if (present(along)) then
       call weights_varied(scheme_weight(scheme), s + 2, tails, heads, spread(step, 1, s + 2), &
            & spread(along, 1, s + 2), source_degree, determines, w)
    else
       call weights_of(scheme_weight(scheme), s + 2, tails, heads, spread(step, 1, s + 2), &
            & source_degree, determines, w)
    end if

  end subroutine step_weights

  !===================================================================!
  ! The derived rows of the whole block: one step's pattern repeated
  ! at every step and every node, its weights taken at each step's
  ! own size, and its vertices mapped onto the block's unknowns. Along
  ! a direction in the steps, the partial of those rows.
  !===================================================================!

  function stage_rows(scheme, nd, n, dt, nodes, along, along2) result(rows)

    class(family), intent(in) :: scheme
    integer      , intent(in) :: nd, n, nodes
    real(dp)     , intent(in) :: dt(:)
    real(dp)     , intent(in), optional :: along(:), along2(:)
    type(stencil) :: rows

    integer , allocatable :: tails(:), heads(:), source_degree(:), determines(:)
    integer , allocatable :: into(:), from(:)
    real(dp), allocatable :: w(:), weight(:)
    integer :: s, kk, p, e, at, per

    s = scheme % num_stages()
    call step_reach(nd, s, tails, heads, source_degree, determines)
    per = size(tails)

    allocate(into((n - 1) * nodes * per), from((n - 1) * nodes * per), &
         &   weight((n - 1) * nodes * per))
    at = 0

    do kk = 2, n
       if (present(along2)) then
          call step_weights(scheme, s, tails, heads, source_degree, determines, dt(kk), w, &
               & along=along(kk), along2=along2(kk))
       else if (present(along)) then
          call step_weights(scheme, s, tails, heads, source_degree, determines, dt(kk), w, &
               & along=along(kk))
       else
          call step_weights(scheme, s, tails, heads, source_degree, determines, dt(kk), w)
       end if
       do p = 1, nodes
          do e = 1, per
             at = at + 1
             into(at)   = unknown_of(heads(e), determines(e), kk, s, nd, p, nodes)
             from(at)   = unknown_of(tails(e), source_degree(e), kk, s, nd, p, nodes)
             weight(at) = w(e)
          end do
       end do
    end do

    ! along a direction in the steps the rows are the partial of the
    ! weights, which the determined component, entering with one,
    ! takes no part in
    if (present(along)) then
       rows = stencil(into, from, -weight, spread(0.0_dp, 1, stage_unknowns(n, s, nd * nodes)), &
            & 'varied stage rows')
    else
       rows = derived_constraints(into, from, weight, stage_unknowns(n, s, nd * nodes), &
            & 'stage rows')
    end if

  end function stage_rows

  !===================================================================!
  ! The whole statement of one stage block. The first instant is
  ! carried at every node; every stage and every later instant is
  ! solved for. Over a field the level below, a stencil over the
  ! nodes, is laid on every stage.
  !===================================================================!

  function stage_block_of(scheme, physics, nd, n, dt, held, nodes, spatial) result(rows)

    class(family)         , intent(in)           :: scheme
    class(nodal_integrand), intent(in)           :: physics
    integer               , intent(in)           :: nd, n
    real(dp)              , intent(in)           :: dt(:), held(:)
    integer               , intent(in), optional :: nodes
    type(stencil)         , intent(in), optional :: spatial
    type(block_residual) :: rows

    integer, allocatable :: slice(:), node(:), moment(:)
    integer :: s, m, u

    m = 1
    if (present(nodes)) m = nodes

    if (n < 2) then
       error stop 'gti_stage: a stage block holds at least one step'
    end if
    if (scheme % history_depth(nd - 1) /= 1) then
       error stop 'gti_stage: a stage family reaches one instant back'
    end if
    if (size(held) /= nd * m) then
       error stop 'gti_stage: the first instant is carried, one value per degree at every node'
    end if

    s = scheme % num_stages()

    rows = block_residual(stage_rows(scheme, nd, n, dt, m), physics, &
         & stage_points(n, s, nd, m), stage_unknowns(n, s, nd * m), nd, nd - 1, &
         & [(u, u = 1, nd * m)], held)

    ! where every unknown lies: the first instant is its own member,
    ! and a step's stages and the instant it arrives at are one member
    call stage_labels(n, s, nd, m, slice, node, moment)
    call rows % placed_in(slice, node, moment)
    if (present(spatial)) call rows % spatial_laid(spatial)

  end function stage_block_of

  !===================================================================!
  ! The member of the time level, the node, and the moment of every
  ! unknown: the first instant is moment one, and each step's stages
  ! and arriving instant follow in order.
  !===================================================================!

  subroutine stage_labels(n, s, nd, nodes, slice, node, moment)

    integer, intent(in) :: n, s, nd, nodes
    integer, allocatable, intent(out) :: slice(:), node(:), moment(:)

    integer :: kk, i, p, d, at, count, width

    width = nd * nodes
    count = stage_unknowns(n, s, width)
    allocate(slice(count), node(count), moment(count))

    do p = 1, nodes
       do d = 1, nd
          call put((p - 1) * nd + d, 1, p, 1)
       end do
    end do
    do kk = 2, n
       do i = 1, s
          do p = 1, nodes
             at = stage_at(kk, i, s, width) + (p - 1) * nd
             do d = 1, nd
                call put(at + d, kk, p, 1 + (kk - 2) * (s + 1) + i)
             end do
          end do
       end do
       do p = 1, nodes
          at = instant_at(kk, s, width) + (p - 1) * nd
          do d = 1, nd
             call put(at + d, kk, p, 1 + (kk - 1) * (s + 1))
          end do
       end do
    end do

  contains

    subroutine put(u, in_slice, at_node, at_moment)

      integer, intent(in) :: u, in_slice, at_node, at_moment

      slice(u)  = in_slice
      node(u)   = at_node
      moment(u) = at_moment

    end subroutine put

  end subroutine stage_labels

end module gti_stage
