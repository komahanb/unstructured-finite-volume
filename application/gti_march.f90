!=====================================================================!
! Building one block and solving it.
!
! The pieces are the same ones the assembly uses - the rows a family
! reaches over, the weights on them, the stencil they make, and the
! statement that adds the governing and carried rows to it - gathered
! here so that a caller marching a block and a caller differentiating
! one write them once.
!
!             WHERE THE BLOCKS OF A HORIZON SIT
!
! horizon_bounds says which instants each block of a chain spans. A
! block reaches back over instants that begin before it does, so
! every block after the first overlaps what came before it by
! exactly what its family reaches. What is done with that overlap -
! the junction, and the layouts either side of it - belongs to
! gti_chain, which marches them.
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
  use operation_minimization  , only : minimizer, relative, absolute, &
       & by_count, by_rate
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

  !-------------------------------------------------------------------!
  ! HOW A MARCH STOPS. A caller that sets nothing gets a tolerance
  ! measured against the imbalance the march began at, and a budget
  ! taken from the rate the march itself shows. The count is a
  ! backstop and not the operative limit.
  !-------------------------------------------------------------------!

  real(dp), save :: stopping_tolerance  = 1.0e-12_dp
  integer , save :: stopping_criterion  = relative
  integer , save :: stopping_budget     = by_rate
  integer , save :: stopping_iterations = 100

  private
  public :: partition, partitioned, scheme_rows, block_of, solved, unknowns_graph
  public :: set_stopping
  public :: horizon_bounds

contains

  !===================================================================!
  ! How every march that follows stops. A criterion or a budget that
  ! is neither of its two stops the program.
  !===================================================================!

  subroutine set_stopping(tolerance, criterion, budget, iterations)

    real(dp), intent(in) :: tolerance
    integer , intent(in) :: criterion, budget, iterations

    if (tolerance <= 0.0_dp) then
       error stop 'gti_march: a tolerance is positive'
    end if
    if (criterion /= relative .and. criterion /= absolute) then
       error stop 'gti_march: a tolerance is measured relative or absolute'
    end if
    if (budget /= by_count .and. budget /= by_rate) then
       error stop 'gti_march: a budget is counted or taken from the rate'
    end if
    if (iterations < 1) then
       error stop 'gti_march: an iteration budget is positive'
    end if

    stopping_tolerance  = tolerance
    stopping_criterion  = criterion
    stopping_budget     = budget
    stopping_iterations = iterations

  end subroutine set_stopping

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
  !
  !             WHAT COUNTS AS SOLVED
  !
  ! A scheme's rows carry a power of the step, so a difference on the
  ! second derivative weighs its sources by the inverse square of it.
  ! Refining the grid therefore raises the size of a residual for the
  ! same trajectory, and the smallest one reachable in the arithmetic
  ! rises with it: at a hundredth of a unit it is near ten to the
  ! minus thirteen, and finer than that it passes any fixed target.
  !
  ! Asked for a fixed one, newton reaches the trajectory in two steps
  ! and then spends its whole budget failing to better it. Measured on a
  ! degree-two problem over three units: a hundred and twenty instants
  ! took a hundred and sixty seconds to produce what forty iterations
  ! produce in a sixth of one, to the same six digits.
  !
  ! So the target is set against the residual the first guess gives,
  ! which is the only scale in the problem that is known before it is
  ! solved, and the budget is a backstop rather than a cost.
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
    call solver % attach(rows, unknowns, unknowns % vertex_set(), count, &
         & held_inputs = [design])

    q = at_first_instant(rows, count)

    solver % max_iterations = stopping_iterations
    solver % tolerance      = stopping_tolerance
    solver % criterion      = stopping_criterion
    solver % budget         = stopping_budget

    call solver % solve(spread(0.0_dp, 1, count), q, achieved)

  end subroutine solved

  !===================================================================!
  ! How large a residual the first guess gives, which is the scale
  ! the target is set against.
  !===================================================================!


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

    type(gmres)        :: krylov
    type(dense_direct) :: factorisation

    if (unknowns <= krylov_above()) then
       ! The matrix here is a tangent frozen at an intermediate newton
       ! iterate, where a singular pivot is a fact about the iterate
       ! rather than a fault, so it is reported and not stopped on.
       factorisation = dense_direct()
       factorisation % singular_reported = .true.
       allocate(inner, source=factorisation)
    else
       krylov = gmres()
       krylov % restart        = min(unknowns, 60)
       krylov % tolerance      = 1.0e-13_dp
       krylov % max_iterations = 4
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
       if (first(b) < 1) then
          error stop 'gti_march: the horizon holds every instant its blocks reach back over'
       end if
    end do

  end subroutine horizon_bounds

end module gti_march
