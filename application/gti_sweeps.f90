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
! For the adjoint, column by column: one partial action per unknown,
! and then solved densely. That does not scale, and what would is
! the stencil's compiled transpose rather than a dense one.
!
! For the tangent no matrix is formed at all past a large block. The
! statement's partial action is already a matvec, so freezing it at
! the trajectory gives a linear operation a krylov solver can be
! driven with directly. Below the same threshold gti_march uses, a
! dense factorisation is the faster of the two and is taken.
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
  use util_tally, only : tally_record, tangent_loops, adjoint_loops
  use util_factorisation, only : dense_factorisation
  use operation_gmres       , only : gmres
  use operation_minimization, only : minimizer
  use operation_linearization, only : linearization, tangent_of

  implicit none

  private
  public :: functional_of, functional_gradient
  public :: design_partial, jacobian_of
  public :: by_tangent, by_adjoint, tangent_solve
  public :: forward_route, reverse_route, route_of, route_substitutions

  integer, parameter :: forward_route = 1
  integer, parameter :: reverse_route = 2
  public :: krylov_above, set_krylov_above

  !===================================================================!
  ! Where the linear solve stops factorising and starts iterating.
  !
  ! A dense factorisation forms the jacobian one column at a time and
  ! then costs the cube of the count. A krylov solver forms no matrix
  ! and costs a matvec per step, so how it fares depends on how many
  ! steps it needs, which is a question about the conditioning of the
  ! block and not about its size.
  !
  ! The two families sit on opposite sides of that. A stage block
  ! weighs its sources by the step, and iterating beats factorising
  ! on one: at six hundred unknowns five seconds against eleven, and
  ! at a thousand ten seconds against sixty. A difference block on a
  ! second derivative weighs its sources by the inverse square of the
  ! step, and iterating does not converge on one at all - at three
  ! hundred and sixty unknowns it does not finish in the time
  ! factorising takes a second to do.
  !
  ! Nothing here preconditions, and without that a krylov solver
  ! cannot be the default. So the default factorises, which always
  ! finishes, and iterating is asked for: set krylov_above to a count
  ! above which to iterate, or to zero to iterate throughout.
  !
  ! It is worth asking for on a large stage block, and worth asking
  ! for where a statement is singular for reasons of its own - a
  ! coarse grid on a diverging problem will do it - because a
  ! factorisation meeting a singular pivot stops the program where an
  ! iteration reports that a row did not converge.
  !
  ! For the second of those to be worth anything the iteration has to
  ! give up rather than grind, so it is held to a few restarts of a
  ! few dozen steps. On a block it suits that is more than it needs;
  ! on one it does not, newton is handed a poor step, fails to
  ! converge, and the row says so - which is the point.
  !===================================================================!

  integer, private, save :: crossing = huge(1)

contains

  pure integer function krylov_above()

    krylov_above = crossing

  end function krylov_above

  subroutine set_krylov_above(count)

    integer, intent(in) :: count

    if (count < 0) then
       error stop 'gti_sweeps: the crossing is not negative'
    end if

    crossing = count

  end subroutine set_krylov_above

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
  ! One solve against the statement's tangent in the state, frozen at
  ! the inputs given. No matrix is formed past a small block: the
  ! statement's own partial action is the matvec.
  !===================================================================!

  subroutine tangent_solve(rows, on, inputs, b, x)

    class(operation)     , intent(in) :: rows
    class(directed_graph), intent(in) :: on
    type(stored_field)   , intent(in) :: inputs(:)
    real(dp)             , intent(in) :: b(:)
    real(dp), allocatable, intent(out) :: x(:)

    type(linearization) :: jacobian
    class(minimizer), allocatable :: solver
    type(gmres) :: krylov
    real(dp) :: achieved

    jacobian = tangent_of(rows, rows % argument(1))
    call jacobian % freeze(inputs)

    if (size(b) <= krylov_above()) then
       allocate(solver, source=dense_direct())
    else
       krylov = gmres()
       krylov % restart        = min(size(b), 60)
       krylov % tolerance      = 1.0e-13_dp
       krylov % max_iterations = 4
       allocate(solver, source=krylov)
    end if

    call solver % attach(jacobian, on, on % vertex_set(), size(b))

    allocate(x(size(b)), source=0.0_dp)
    call solver % solve(b, x, achieved)

  end subroutine tangent_solve

  !===================================================================!
  ! The tangent: one solve in the state, then the gradient read
  ! along what it gives.
  !===================================================================!

  real(dp) function by_tangent(factor, g, design_rate, explicit) result(df)

    type(dense_factorisation), intent(in) :: factor
    real(dp)                 , intent(in) :: g(:), design_rate(:), explicit

    real(dp), allocatable :: w(:)

    call tally_record(tangent_loops)
    call factor % substitute(-design_rate, w, transposed=.false.)
    df = explicit + dot_product(g, w)

  end function by_tangent

  !===================================================================!
  ! The adjoint: one substitution against the transpose, then the
  ! statement's design partial read along what it gives.
  !===================================================================!

  real(dp) function by_adjoint(factor, g, design_rate, explicit) result(df)

    type(dense_factorisation), intent(in) :: factor
    real(dp)                 , intent(in) :: g(:), design_rate(:), explicit

    real(dp), allocatable :: lambda(:)

    call tally_record(adjoint_loops)
    call factor % substitute(g, lambda, transposed=.true.)
    df = explicit - dot_product(lambda, design_rate)

  end function by_adjoint

  !===================================================================!
  ! THE GATE. Which route computes the m-th derivative of n_f
  ! functionals in n_d designs, by the count of substitutions each
  ! costs against one kept factorisation:
  !
  !      forward, m times            C(n_d + m - 1, m)
  !      forward m-1 times over      (1 + n_f) C(n_d + m - 2, m - 1)
  !      one reverse
  !
  ! The ratio of the first to the second is (n_d + m - 1) / (m (1 + n_f)),
  ! so the forward route is the cheaper exactly where n_d <= m n_f. At
  ! m = 1 that is the rule for a gradient, n_d <= n_f. At one design
  ! and one functional the forward route costs m substitutions and the
  ! other 2 m, at every order.
  !
  ! A count below one, or an order below one, stops the program.
  !===================================================================!

  pure integer function route_of(num_designs, num_functionals, order) result(route)

    integer, intent(in) :: num_designs, num_functionals, order

    if (num_designs < 1 .or. num_functionals < 1 .or. order < 1) then
       error stop 'gti_sweeps: a route is chosen for at least one design, one functional and order one'
    end if

    if (num_designs <= order * num_functionals) then
       route = forward_route
    else
       route = reverse_route
    end if

  end function route_of

  !===================================================================!
  ! The substitutions a route costs at one order, per block. What the
  ! accounting layer counts is compared against this.
  !===================================================================!

  pure integer function route_substitutions(route, num_designs, num_functionals, order) &
       & result(count)

    integer, intent(in) :: route, num_designs, num_functionals, order

    select case (route)
    case (forward_route)
       count = choose(num_designs + order - 1, order)
    case (reverse_route)
       count = (1 + num_functionals) * choose(num_designs + order - 2, order - 1)
    case default
       error stop 'gti_sweeps: a route is forward or reverse'
    end select

  end function route_substitutions

  pure integer function choose(n, k) result(c)

    integer, intent(in) :: n, k

    integer :: i

    c = 1
    do i = 1, k
       c = c * (n - k + i) / i
    end do

  end function choose

end module gti_sweeps
