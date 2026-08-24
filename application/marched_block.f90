! One block solved.
!
! The rows a scheme makes, the row the physics makes, and the rows of
! the instants carried in are added into one statement, and newton
! drives it to zero. The tangent it linearizes with is exact: the
! stencil is its own jacobian and the physics differentiates its own
! rule, so nothing anywhere in the chain is differenced.
!
! The check is exact rather than comparative. Van der Pol at a design
! of zero is
!
!      q" - 0 (1 - q^2) q' + q  =  q" + q  =  0 ,
!
! the harmonic oscillator, whose solution through q(0) = 1, q'(0) = 0
! is the cosine. So the marched value is compared against cos(t)
! directly, and the steps are then halved: the error must fall by two
! raised to the scheme's order, which is what the printed ratio is.
!
! The last part marches the true oscillator, at a design of one,
! where the statement is nonlinear and newton has to iterate.
program marched_block

  use iso_fortran_env       , only : dp => REAL64
  use view_directed_stored  , only : stored_directed_graph
  use field_calculus        , only : field
  use field_stored          , only : stored_field
  use operation_stencil     , only : stencil
  use operation_newton      , only : newton
  use operation_dense_direct, only : dense_direct
  use operation_family      , only : family
  use operation_family_bdf  , only : bdf_family
  use operation_family_adams, only : adams_family
  use operation_grid        , only : uniform_grid
  use operation_weight      , only : scheme_weight
  use operation_scheme_stencil, only : derived_constraints
  use physics_vanderpol     , only : van_der_pol
  use gti_expansion         , only : block_reach
  use gti_block             , only : block_residual

  implicit none

  integer , parameter :: max_state_degree = 2
  integer , parameter :: degrees = max_state_degree + 1
  real(dp), parameter :: duration = 2.0_dp

  call order_of('bdf 2',           bdf_family(2),   2.0_dp)
  call order_of('bdf 3',           bdf_family(3),   3.0_dp)
  call order_of('adams-moulton 3', adams_family(3), 3.0_dp)
  call nonlinear()

contains

  pure integer function unknown(instant, degree)

    integer, intent(in) :: instant, degree

    unknown = (instant - 1) * degrees + degree + 1

  end function unknown

  !-------------------------------------------------------------------!
  ! The instants of a uniform partition of the duration.
  !-------------------------------------------------------------------!

  subroutine partition(n, dt, t)

    integer, intent(in) :: n
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

  !-------------------------------------------------------------------!
  ! The derived rows of a block, as a stencil.
  !-------------------------------------------------------------------!

  function scheme_rows(scheme, n, dt) result(rows)

    class(family), intent(in) :: scheme
    integer      , intent(in) :: n
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
         & [(unknown(heads(e), determines(e)), e = 1, size(heads))], &
         & [(unknown(tails(e), source_degree(e)), e = 1, size(tails))], &
         & w, n * degrees, 'derived rows')

  end function scheme_rows

  !-------------------------------------------------------------------!
  ! March a block: the steps, the rows, the carried instants, and
  ! newton over all of them.
  !-------------------------------------------------------------------!

  subroutine march(scheme, n, design_value, q, t, achieved)

    class(family), intent(in)  :: scheme
    integer      , intent(in)  :: n
    real(dp)     , intent(in)  :: design_value
    real(dp), allocatable, intent(out) :: q(:), t(:)
    real(dp)     , intent(out) :: achieved

    type(block_residual) :: rows
    type(newton) :: solver
    type(stored_directed_graph) :: unknowns
    type(stored_field) :: design
    real(dp), allocatable :: dt(:), held(:)
    integer , allocatable :: carried(:)
    integer :: h, k, d

    call partition(n, dt, t)
    h = scheme % history_depth()

    carried = [(((k - 1) * degrees + d + 1, d = 0, degrees - 1), k = 1, h)]
    held    = [((exact(d, t(k)), d = 0, degrees - 1), k = 1, h)]

    rows = block_residual(scheme_rows(scheme, n, dt), van_der_pol(max_state_degree), &
         & n, degrees, scheme % primary_degree(max_state_degree), carried, held)

    unknowns = stored_directed_graph(n * degrees, tails=[integer ::], heads=[integer ::])
    design   = stored_field('nu', unknowns % vertex_set(), n)
    call design % set_real_vector(spread(design_value, 1, n))

    allocate(solver % inner, source=dense_direct())
    solver % tolerance = 1.0e-12_dp
    call solver % attach(rows, unknowns, unknowns % vertex_set(), n * degrees, &
         & held_inputs = [design])

    allocate(q(n * degrees), source=0.0_dp)
    call solver % solve(spread(0.0_dp, 1, n * degrees), q, achieved)

  end subroutine march

  !-------------------------------------------------------------------!
  ! The d-th derivative of the cosine.
  !-------------------------------------------------------------------!

  pure real(dp) function exact(d, t) result(q)

    integer , intent(in) :: d
    real(dp), intent(in) :: t

    select case (mod(d, 4))
    case (0)
       q =  cos(t)
    case (1)
       q = -sin(t)
    case (2)
       q = -cos(t)
    case default
       q =  sin(t)
    end select

  end function exact

  pure real(dp) function worst(q, t) result(e)

    real(dp), intent(in) :: q(:), t(:)

    integer :: k

    e = 0.0_dp
    do k = 1, size(t)
       e = max(e, abs(q(unknown(k, 0)) - exact(0, t(k))))
    end do

  end function worst

  !-------------------------------------------------------------------!
  ! The error at two resolutions, and the ratio between them.
  !-------------------------------------------------------------------!

  subroutine order_of(title, scheme, expected)

    character(len=*), intent(in) :: title
    class(family)   , intent(in) :: scheme
    real(dp)        , intent(in) :: expected

    real(dp), allocatable :: q(:), t(:)
    real(dp) :: e(4), achieved
    integer :: level, steps

    do level = 1, 4
       steps = 10 * 2 ** (level - 1)
       call march(scheme, steps + 1, 0.0_dp, q, t, achieved)
       e(level) = worst(q, t)
    end do

    write(*,'(a)')          ' '
    write(*,'(a)')          ' ' // title // ' on the harmonic oscillator'
    write(*,'(a,4i11)')     '   steps                      ', [(10 * 2 ** (level - 1), level = 1, 4)]
    write(*,'(a,4es11.3)')  '   worst error                ', e
    write(*,'(a,33x,3f11.3)') '   ratio                    ', e(1:3) / e(2:4)
    write(*,'(a,f11.3)')    '   two to the scheme order    ', 2.0_dp ** expected
    write(*,'(a,es11.3)')   '   residual newton achieved   ', achieved

  end subroutine order_of

  !-------------------------------------------------------------------!
  ! The true oscillator, where the statement is nonlinear.
  !-------------------------------------------------------------------!

  subroutine nonlinear()

    real(dp), allocatable :: q(:), t(:)
    real(dp) :: achieved

    call march(bdf_family(2), 41, 1.0_dp, q, t, achieved)

    write(*,'(a)')        ' '
    write(*,'(a)')        ' van der pol at a design of one, bdf 2, 40 steps'
    write(*,'(a,es11.3)') '   residual newton achieved   ', achieved
    write(*,'(a,3f11.5)') '   the last instant, q q'' q"  ', &
         & q(unknown(size(t), 0)), q(unknown(size(t), 1)), q(unknown(size(t), 2))

  end subroutine nonlinear

end program marched_block
