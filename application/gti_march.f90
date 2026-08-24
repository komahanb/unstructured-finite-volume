!=====================================================================!
! Building one block and solving it.
!
! The pieces are the same ones the assembly uses - the rows a family
! reaches over, the weights on them, the stencil they make, and the
! statement that adds the governing and carried rows to it - gathered
! here so that a caller marching a block and a caller differentiating
! one write them once.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_march

  use iso_fortran_env         , only : dp => REAL64
  use view_directed_stored    , only : stored_directed_graph
  use field_calculus          , only : field
  use field_stored            , only : stored_field
  use operation_stencil       , only : stencil
  use operation_newton        , only : newton
  use operation_dense_direct  , only : dense_direct
  use operation_family        , only : family
  use operation_grid          , only : uniform_grid
  use operation_weight        , only : scheme_weight
  use operation_scheme_stencil, only : derived_constraints
  use physics_integrand       , only : nodal_integrand
  use gti_expansion           , only : block_reach
  use gti_block               , only : block_residual

  implicit none

  private
  public :: partition, scheme_rows, block_of, solved, unknowns_graph

contains

  !===================================================================!
  ! A uniform partition of the duration, and the instants it makes.
  !===================================================================!

  subroutine partition(duration, n, dt, t)

    real(dp), intent(in) :: duration
    integer , intent(in) :: n
    real(dp), allocatable, intent(out) :: dt(:), t(:)

    type(stored_directed_graph) :: instants
    type(stored_field) :: knobs
    type(uniform_grid) :: steps
    class(field), allocatable :: out
    integer :: k

    instants = stored_directed_graph(n, tails=[integer ::], heads=[integer ::])
    knobs    = stored_field('design', instants % vertex_set(), 1)
    call knobs % set_real_vector([0.0_dp])

    steps = uniform_grid(duration)
    call steps % apply(instants, [knobs], out)
    call out % real_vector(dt)

    allocate(t(n))
    t(1) = 0.0_dp
    do k = 2, n
       t(k) = t(k - 1) + dt(k)
    end do

  end subroutine partition

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

    type(stored_directed_graph) :: edges
    type(stored_field) :: steps, source_field, condition_field
    type(scheme_weight) :: weights
    class(field), allocatable :: out
    integer , allocatable :: tails(:), heads(:), source_degree(:), determines(:)
    real(dp), allocatable :: w(:)
    integer :: e

    call block_reach(scheme, degrees, n, tails, heads, source_degree, determines)

    edges           = stored_directed_graph(n, tails=tails, heads=heads)
    steps           = stored_field('dt', edges % vertex_set(), n)
    source_field    = stored_field('source degree', edges % edge_set(), size(tails))
    condition_field = stored_field('determines', edges % edge_set(), size(tails))
    call steps           % set_real_vector(dt)
    call source_field    % set_integer_vector(source_degree)
    call condition_field % set_integer_vector(determines)

    weights = scheme_weight(scheme)
    call weights % apply(edges, [steps, source_field, condition_field], out)
    call out % real_vector(w)

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

    h = scheme % history_depth()
    carried = [((unknown(k, d, degrees), d = 0, degrees - 1), k = 1, h)]

    if (size(held) /= size(carried)) then
       error stop 'gti_march: one value per carried component'
    end if

    rows = block_residual(scheme_rows(scheme, degrees, n, dt), physics, n, degrees, &
         & scheme % primary_degree(degrees - 1), carried, held)

  end function block_of

  !===================================================================!
  ! Newton over the whole block. The design is held while the state
  ! varies, which is what a minimizer supplies as an extra input.
  !===================================================================!

  subroutine solved(rows, n, degrees, design_value, q, achieved)

    type(block_residual), intent(in)  :: rows
    integer             , intent(in)  :: n, degrees
    real(dp)            , intent(in)  :: design_value
    real(dp), allocatable, intent(out) :: q(:)
    real(dp)            , intent(out) :: achieved

    type(newton) :: solver
    type(stored_directed_graph) :: unknowns
    type(stored_field) :: design

    unknowns = unknowns_graph(n, degrees)
    design   = stored_field('nu', unknowns % vertex_set(), n)
    call design % set_real_vector(spread(design_value, 1, n))

    allocate(solver % inner, source=dense_direct())
    solver % tolerance = 1.0e-12_dp
    call solver % attach(rows, unknowns, unknowns % vertex_set(), n * degrees, &
         & held_inputs = [design])

    allocate(q(n * degrees), source=0.0_dp)
    call solver % solve(spread(0.0_dp, 1, n * degrees), q, achieved)

  end subroutine solved

end module gti_march
