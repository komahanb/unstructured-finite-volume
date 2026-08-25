!=====================================================================!
! Building one block and solving it.
!
! The pieces are the same ones the assembly uses - the rows a family
! reaches over, the weights on them, the stencil they make, and the
! statement that adds the governing and carried rows to it - gathered
! here so that a caller marching a block and a caller differentiating
! one write them once.
!
!             THE JUNCTION BETWEEN BLOCKS
!
! A horizon is partitioned among blocks, and a block reaches back
! over instants that begin before it does. Those instants are the
! previous block's last, so every block after the first overlaps its
! predecessor by exactly what its family reaches, and what it carries
! there is what the previous block computed. That hand-over is the
! junction, and it is why a horizon can change scheme part way
! through: a block asks its predecessor only for components, never
! for the rows that made them.
!
! The state is laid out instant by instant, degrees within an
! instant, so a block's own numbering and the horizon's differ by a
! single shift and the junction is a contiguous copy. The first
! block has no predecessor and carries the initial conditions
! instead.
!
! A block must add more instants than its family reaches back over,
! or it would consist of nothing but what it was given.
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
  use operation_minimization  , only : minimizer
  use operation_dense_direct  , only : dense_direct
  use operation_gmres         , only : gmres
  use operation_family        , only : family
  use operation_grid          , only : grid, uniform_grid
  use operation_weight        , only : scheme_weight
  use operation_scheme_stencil, only : derived_constraints
  use physics_integrand       , only : nodal_integrand
  use gti_expansion           , only : block_reach, family_holder
  use gti_block               , only : block_residual
  use gti_sweeps              , only : krylov_above

  implicit none

  private
  public :: partition, partitioned, scheme_rows, block_of, solved, unknowns_graph
  public :: horizon_bounds, marched_horizon

contains

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
  !===================================================================!

  subroutine solved(rows, design_value, q, achieved)

    type(block_residual), intent(in)  :: rows
    real(dp)            , intent(in)  :: design_value
    real(dp), allocatable, intent(out) :: q(:)
    real(dp)            , intent(out) :: achieved

    type(newton) :: solver
    type(stored_directed_graph) :: unknowns
    type(stored_field) :: design
    integer :: count

    count    = rows % num_unknowns()
    unknowns = stored_directed_graph(count, tails=[integer ::], heads=[integer ::])
    design   = stored_field('nu', unknowns % vertex_set(), rows % num_points())
    call design % set_real_vector(spread(design_value, 1, rows % num_points()))

    allocate(solver % inner, source=inner_solver(count))
    solver % tolerance = 1.0e-12_dp
    call solver % attach(rows, unknowns, unknowns % vertex_set(), count, &
         & held_inputs = [design])

    q = at_first_instant(rows, count)
    call solver % solve(spread(0.0_dp, 1, count), q, achieved)

  end subroutine solved

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

  function inner_solver(unknowns) result(inner)

    integer, intent(in) :: unknowns
    class(minimizer), allocatable :: inner

    type(gmres) :: krylov

    if (unknowns <= krylov_above()) then
       allocate(inner, source=dense_direct())
    else
       krylov = gmres()
       krylov % restart    = min(unknowns, 200)
       krylov % tolerance  = 1.0e-13_dp
       krylov % max_iterations = 4 * unknowns
       allocate(inner, source=krylov)
    end if

  end function inner_solver

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
    end do

  end subroutine horizon_bounds

  !===================================================================!
  ! The whole horizon marched, block after block. Each block is given
  ! what its predecessor computed at the instants they share, solved
  ! on its own, and its state written back into the horizon.
  !===================================================================!

  subroutine marched_horizon(schemes, added, physics, degrees, duration, &
       & design_value, initial, q, achieved)

    type(family_holder)   , intent(in) :: schemes(:)
    integer               , intent(in) :: added(:), degrees
    class(nodal_integrand), intent(in) :: physics
    real(dp)              , intent(in) :: duration, design_value, initial(:)
    real(dp), allocatable , intent(out) :: q(:)
    real(dp)              , intent(out) :: achieved

    type(block_residual) :: rows
    integer , allocatable :: first(:), last(:)
    real(dp), allocatable :: dt(:), t(:), block_state(:)
    integer :: b, n, shift, held_size
    real(dp) :: block_achieved

    call horizon_bounds(schemes, added, degrees - 1, first, last)
    call partition(duration, last(size(added)), dt, t)

    allocate(q(last(size(added)) * degrees), source=0.0_dp)
    q(1:size(initial)) = initial

    achieved = 0.0_dp

    do b = 1, size(added)
       n         = last(b) - first(b) + 1
       shift     = (first(b) - 1) * degrees
       held_size = schemes(b) % scheme % history_depth(degrees - 1) * degrees

       rows = block_of(schemes(b) % scheme, physics, degrees, n, &
            & dt(first(b):last(b)), q(shift + 1:shift + held_size))

       call solved(rows, design_value, block_state, block_achieved)

       q(shift + 1:shift + n * degrees) = block_state
       achieved = max(achieved, block_achieved)
    end do

  end subroutine marched_horizon

end module gti_march
