!=====================================================================!
! The sensitivity of a functional to a design, by the tangent and
! by the adjoint.
!
! A block's statement R(q, x) = 0 determines q from x, and the
! functional is a sum over the instants,
!
!      f  =  sum over k of  dt_k F(q_k, x) .
!
! Differentiating the statement gives J dq/dx = -dR/dx, so
!
!      tangent    solve J w = -dR/dx  once, then  df/dx = f_x + g.w
!      adjoint    solve J^T l = g     once, then  df/dx = f_x - l.dR/dx
!
! with J the jacobian in the state, g the gradient of f in the state
! and f_x its own partial in the design. The two read the same
! three objects and must agree to round-off; that they do is what
! this module exists to make checkable.
!
! Each of the three comes from a partial action that is exact: the
! scheme's rows are linear and are their own jacobian, and the
! physics and the integrand differentiate their own rules. Nothing
! here is differenced.
!
!             HOW THE JACOBIAN IS FORMED
!
! Column by column, one partial action per unknown, and then solved
! densely. That is the same cost dense_direct already pays inside
! newton and it does not scale: the block residual supplies a matvec
! through its partial action, so a krylov solver needs no matrix at
! all, and the adjoint needs the stencil's compiled transpose rather
! than a dense one. Both are what a long horizon needs, and
! neither changes the numbers below.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_sweeps

  use iso_fortran_env       , only : dp => REAL64
  use operation_action      , only : operation, variation
  use view_directed         , only : directed_graph
  use view_directed_stored  , only : stored_directed_graph
  use graph_fractal         , only : graph
  use field_calculus        , only : field
  use field_stored          , only : stored_field
  use operation_stencil     , only : stencil
  use operation_dense_direct, only : dense_direct

  implicit none

  private
  public :: functional_of, functional_gradient
  public :: design_partial, jacobian_of
  public :: by_tangent, by_adjoint, dense_solve

contains

  subroutine applied(action, on, inputs, y)

    class(operation)     , intent(in) :: action
    class(directed_graph), intent(in) :: on
    type(stored_field)   , intent(in) :: inputs(:)
    real(dp), allocatable, intent(out) :: y(:)

    class(field), allocatable :: out

    call action % apply(on, inputs, out)
    call out % real_vector(y)

  end subroutine applied

  !===================================================================!
  ! The partial action of a statement along one direction in one of
  ! its arguments.
  !===================================================================!

  subroutine varied(action, on, inputs, which, domain, v, y)

    class(operation)     , intent(in) :: action
    class(directed_graph), intent(in) :: on
    type(stored_field)   , intent(in) :: inputs(:)
    integer              , intent(in) :: which
    type(graph)          , intent(in) :: domain
    real(dp)             , intent(in) :: v(:)
    real(dp), allocatable, intent(out) :: y(:)

    type(stored_field) :: direction
    class(field), allocatable :: out

    direction = stored_field('direction', domain, size(v))
    call direction % set_real_vector(v)

    call action % partial_action(on, inputs, &
         & [variation(action % argument(which), direction)], out)
    call out % real_vector(y)

  end subroutine varied

  !===================================================================!
  ! The functional: the integrand at every instant, weighted by the
  ! step that ends there. The first instant carries no step and so
  ! contributes nothing.
  !===================================================================!

  real(dp) function functional_of(integrand, instants, inputs, dt) result(f)

    class(operation)     , intent(in) :: integrand
    class(directed_graph), intent(in) :: instants
    type(stored_field)   , intent(in) :: inputs(:)
    real(dp)             , intent(in) :: dt(:)

    real(dp), allocatable :: values(:)

    call applied(integrand, instants, inputs, values)
    f = sum(dt * values)

  end function functional_of

  !===================================================================!
  ! Its gradient in the state. The integrand at one instant reads
  ! only that instant's components, so one partial action per degree
  ! gives the whole gradient rather than one per unknown.
  !===================================================================!

  subroutine functional_gradient(integrand, instants, inputs, dt, n, degrees, &
       & state_domain, g)

    class(operation)     , intent(in) :: integrand
    class(directed_graph), intent(in) :: instants
    type(stored_field)   , intent(in) :: inputs(:)
    real(dp)             , intent(in) :: dt(:)
    integer              , intent(in) :: n, degrees
    type(graph)          , intent(in) :: state_domain
    real(dp), allocatable, intent(out) :: g(:)

    real(dp), allocatable :: v(:), rate(:)
    integer :: d, k

    allocate(g(n * degrees), source=0.0_dp)
    allocate(v(n * degrees))

    do d = 0, degrees - 1
       v = 0.0_dp
       do k = 1, n
          v((k - 1) * degrees + d + 1) = 1.0_dp
       end do

       call varied(integrand, instants, inputs, 1, state_domain, v, rate)

       do k = 1, n
          g((k - 1) * degrees + d + 1) = dt(k) * rate(k)
       end do
    end do

  end subroutine functional_gradient

  !===================================================================!
  ! The statement's partial in the design, along the direction that
  ! varies every instant's design value together, which is a single
  ! design number for the whole block.
  !===================================================================!

  subroutine design_partial(rows, unknowns, inputs, n, design_domain, d)

    class(operation)     , intent(in) :: rows
    class(directed_graph), intent(in) :: unknowns
    type(stored_field)   , intent(in) :: inputs(:)
    integer              , intent(in) :: n
    type(graph)          , intent(in) :: design_domain
    real(dp), allocatable, intent(out) :: d(:)

    call varied(rows, unknowns, inputs, 2, design_domain, spread(1.0_dp, 1, n), d)

  end subroutine design_partial

  !===================================================================!
  ! The jacobian in the state, one column per unknown.
  !===================================================================!

  subroutine jacobian_of(rows, unknowns, inputs, num_unknowns, state_domain, a)

    class(operation)     , intent(in) :: rows
    class(directed_graph), intent(in) :: unknowns
    type(stored_field)   , intent(in) :: inputs(:)
    integer              , intent(in) :: num_unknowns
    type(graph)          , intent(in) :: state_domain
    real(dp), allocatable, intent(out) :: a(:,:)

    real(dp), allocatable :: v(:), column(:)
    integer :: j

    allocate(a(num_unknowns, num_unknowns))
    allocate(v(num_unknowns), source=0.0_dp)

    do j = 1, num_unknowns
       v    = 0.0_dp
       v(j) = 1.0_dp
       call varied(rows, unknowns, inputs, 1, state_domain, v, column)
       a(:, j) = column
    end do

  end subroutine jacobian_of

  !===================================================================!
  ! One dense solve, through the tower's own machinery: the matrix as
  ! a stencil, transposed on request, driven to the given right side.
  !===================================================================!

  subroutine dense_solve(a, b, transposed, x)

    real(dp), intent(in)  :: a(:,:), b(:)
    logical , intent(in)  :: transposed
    real(dp), allocatable, intent(out) :: x(:)

    type(stencil) :: matrix
    type(dense_direct) :: solver
    type(stored_directed_graph) :: on
    real(dp) :: achieved

    matrix = stencil(a, 'jacobian')
    on     = stored_directed_graph(size(b), tails=[integer ::], heads=[integer ::])

    if (transposed) then
       call solver % attach(matrix % transpose(), on, on % vertex_set(), size(b))
    else
       call solver % attach(matrix, on, on % vertex_set(), size(b))
    end if

    allocate(x(size(b)), source=0.0_dp)
    call solver % solve(b, x, achieved)

  end subroutine dense_solve

  !===================================================================!
  ! The tangent: one solve in the state, then the gradient read
  ! along what it gives.
  !===================================================================!

  real(dp) function by_tangent(a, g, design_rate, explicit) result(df)

    real(dp), intent(in) :: a(:,:), g(:), design_rate(:), explicit

    real(dp), allocatable :: w(:)

    call dense_solve(a, -design_rate, .false., w)
    df = explicit + dot_product(g, w)

  end function by_tangent

  !===================================================================!
  ! The adjoint: one solve against the transpose, then the
  ! statement's design partial read along what it gives.
  !===================================================================!

  real(dp) function by_adjoint(a, g, design_rate, explicit) result(df)

    real(dp), intent(in) :: a(:,:), g(:), design_rate(:), explicit

    real(dp), allocatable :: lambda(:)

    call dense_solve(a, g, .true., lambda)
    df = explicit - dot_product(lambda, design_rate)

  end function by_adjoint

end module gti_sweeps
