!=====================================================================!
! The sensitivity of a functional over a horizon of several blocks.
!
! A block reaches back over its predecessor's last instants and
! carries there what the predecessor computed. That junction is a
! dependency, so a horizon's trajectory depends on its design through
! every block in turn, and a sensitivity has to travel the same way.
!
!      tangent    forward, block after block: each block's own
!                 partial in the design, plus what its predecessor's
!                 sensitivity hands it at the instants they share
!      adjoint    backward, block before block: each block's own
!                 gradient, plus what its successor's costate hands
!                 it at the instants they share
!
! The junction row is the identity less what it was given, so its
! partial in the predecessor's state is minus one, and both hand-overs
! are that one number: a carried entry moves across unchanged. Which
! is why the two directions are so nearly the same code, and why they
! must agree.
!
!             WHOSE INSTANT IS IT
!
! An instant shared by two blocks is computed by the earlier one and
! carried by the later. The functional counts it once, under the
! block that computed it, so each block's gradient is the horizon's
! gradient with its carried rows struck out. Counting it twice would
! double every shared instant's contribution and pass no check.
!
!             WHAT IS REFUSED
!
! A horizon of no blocks, and a set of systems whose offsets do not
! run in order, since a junction only ever hands forward.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_horizon

  use iso_fortran_env       , only : dp => REAL64
  use view_directed_stored  , only : stored_directed_graph
  use field_stored          , only : stored_field
  use operation_family      , only : family
  use physics_integrand     , only : nodal_integrand
  use gti_expansion         , only : family_holder
  use gti_block             , only : block_residual
  use gti_march             , only : horizon_bounds, partition, block_of, unknowns_graph
  use gti_sweeps            , only : functional_gradient, design_partial, jacobian_of, &
       & dense_solve

  implicit none

  private
  public :: block_system, horizon_systems, horizon_by_tangent, horizon_by_adjoint

  !===================================================================!
  ! What one block contributes to a sensitivity: its jacobian in the
  ! state, its partial in the design, the part of the functional's
  ! gradient it owns, where it sits in the horizon, and how much of
  ! it was carried in.
  !===================================================================!

  type :: block_system

     real(dp), allocatable :: a(:,:)
     real(dp), allocatable :: rate(:)
     real(dp), allocatable :: g(:)
     integer               :: at = 0
     integer               :: carried = 0

  end type block_system

contains

  !===================================================================!
  ! Every block's system, at the trajectory already marched. The
  ! gradient is formed once over the whole horizon and then shared
  ! out, each block taking the rows it computed and none of the rows
  ! it was given.
  !===================================================================!

  subroutine horizon_systems(schemes, added, physics, integrand, degrees, &
       & duration, design_value, q, systems, gradient)

    type(family_holder)   , intent(in) :: schemes(:)
    integer               , intent(in) :: added(:), degrees
    class(nodal_integrand), intent(in) :: physics, integrand
    real(dp)              , intent(in) :: duration, design_value, q(:)
    type(block_system), allocatable, intent(out) :: systems(:)
    real(dp)          , allocatable, intent(out) :: gradient(:)

    integer , allocatable :: first(:), last(:)
    real(dp), allocatable :: dt(:), t(:)
    integer :: b, n

    if (size(added) < 1) then
       error stop 'gti_horizon: a horizon holds at least one block'
    end if

    call horizon_bounds(schemes, added, degrees - 1, first, last)
    call partition(duration, last(size(added)), dt, t)
    call whole_gradient(integrand, degrees, last(size(added)), design_value, q, dt, gradient)

    allocate(systems(size(added)))

    do b = 1, size(added)
       n = last(b) - first(b) + 1
       systems(b) % at      = (first(b) - 1) * degrees
       systems(b) % carried = schemes(b) % scheme % history_depth(degrees - 1) * degrees

       call one_system(schemes(b) % scheme, physics, degrees, n, &
            & dt(first(b):last(b)), design_value, q, gradient, systems(b))
    end do

  end subroutine horizon_systems

  !===================================================================!
  ! The functional's gradient over every instant of the horizon.
  !===================================================================!

  subroutine whole_gradient(integrand, degrees, n, design_value, q, dt, gradient)

    class(nodal_integrand), intent(in) :: integrand
    integer               , intent(in) :: degrees, n
    real(dp)              , intent(in) :: design_value, q(:), dt(:)
    real(dp), allocatable , intent(out) :: gradient(:)

    type(stored_directed_graph) :: unknowns, instants
    type(stored_field) :: state, knobs

    unknowns = unknowns_graph(n, degrees)
    instants = stored_directed_graph(n, tails=[integer ::], heads=[integer ::])
    state    = stored_field('state', unknowns % vertex_set(), size(q))
    knobs    = stored_field('design', unknowns % vertex_set(), n)
    call state % set_real_vector(q)
    call knobs % set_real_vector(spread(design_value, 1, n))

    call functional_gradient(integrand, instants, [state, knobs], dt, n, degrees, &
         & unknowns % vertex_set(), gradient)

  end subroutine whole_gradient

  !===================================================================!
  ! One block's jacobian, design partial, and share of the gradient.
  ! The rows it was given are struck out of its share, because the
  ! block before it computed those instants and already counts them.
  !===================================================================!

  subroutine one_system(scheme, physics, degrees, n, dt, design_value, q, gradient, s)

    class(family)         , intent(in)    :: scheme
    class(nodal_integrand), intent(in)    :: physics
    integer               , intent(in)    :: degrees, n
    real(dp)              , intent(in)    :: dt(:), design_value, q(:), gradient(:)
    type(block_system)    , intent(inout) :: s

    type(block_residual) :: rows
    type(stored_directed_graph) :: unknowns
    type(stored_field) :: state, knobs
    integer :: unknown_count

    unknown_count = n * degrees

    rows = block_of(scheme, physics, degrees, n, dt, &
         & q(s % at + 1:s % at + s % carried))

    unknowns = unknowns_graph(n, degrees)
    state    = stored_field('state', unknowns % vertex_set(), unknown_count)
    knobs    = stored_field('design', unknowns % vertex_set(), n)
    call state % set_real_vector(q(s % at + 1:s % at + unknown_count))
    call knobs % set_real_vector(spread(design_value, 1, n))

    call design_partial(rows, unknowns, [state, knobs], n, unknowns % vertex_set(), s % rate)
    call jacobian_of(rows, unknowns, [state, knobs], unknown_count, &
         & unknowns % vertex_set(), s % a)

    s % g = gradient(s % at + 1:s % at + unknown_count)
    s % g(1:s % carried) = 0.0_dp

  end subroutine one_system

  !===================================================================!
  ! Forward. Each block solves for its own sensitivity, with the
  ! rows it carried set to what its predecessor already found there.
  ! The horizon's sensitivity is the gradient read along the whole.
  !===================================================================!

  real(dp) function horizon_by_tangent(systems, gradient) result(df)

    type(block_system), intent(in) :: systems(:)
    real(dp)          , intent(in) :: gradient(:)

    real(dp), allocatable :: w(:), b(:), whole(:)
    integer :: at, i, previous

    allocate(whole(size(gradient)), source=0.0_dp)

    do at = 1, size(systems)
       b = -systems(at) % rate

       do i = 1, systems(at) % carried
          previous = systems(at) % at + i
          b(i) = whole(previous)
       end do

       call dense_solve(systems(at) % a, b, .false., w)
       whole(systems(at) % at + 1:systems(at) % at + size(w)) = w
    end do

    df = dot_product(gradient, whole)

  end function horizon_by_tangent

  !===================================================================!
  ! Backward. Each block solves against its own transpose, with its
  ! own share of the gradient plus whatever its successor's costate
  ! left on the rows they share.
  !===================================================================!

  real(dp) function horizon_by_adjoint(systems) result(df)

    type(block_system), intent(in) :: systems(:)

    real(dp), allocatable :: lambda(:), b(:), handed(:)
    integer :: at, i, into

    df = 0.0_dp
    allocate(handed(0))

    do at = size(systems), 1, -1
       b = systems(at) % g

       if (at < size(systems)) then
          do i = 1, systems(at + 1) % carried
             into = systems(at + 1) % at + i - systems(at) % at
             b(into) = b(into) + handed(i)
          end do
       end if

       call dense_solve(systems(at) % a, b, .true., lambda)
       df = df - dot_product(lambda, systems(at) % rate)
       handed = lambda
    end do

  end function horizon_by_adjoint

end module gti_horizon
