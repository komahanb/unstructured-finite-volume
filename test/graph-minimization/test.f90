!=====================================================================!
! The minimization suite: rung 4's acceptance.
!
! The whole tower is exercised by one problem: heat conduction on a
! chain of three cells, boundaries fixed at 0 and 10,
!
!       0 |--- (1) --- (2) --- (3) ---| 10        k = 1, A = 1
!         d 0.5     1       1     0.5
!
! whose exact solution is q = [5/3, 5, 25/3]. The mesh stores the
! measurements, the constitution supplies the coefficients, the
! calculus builds the rows, and the minimization drives the residual
! down - every level doing only its own job, and the solver speaking
! only its own words: matvec, inner product, norm, diagonal, sweep
! order, each a delegation to an engine entry.
!=====================================================================!

module cubic_statement_fixture

  use iso_fortran_env, only : dp => REAL64
  use operation_action, only : operation, binding, bound_real_vector
  use view_directed, only : directed_graph
  use field_calculus, only : field
  ! An action names a domain and counts it: identity and count is
  ! the whole of it.
  use graph_fractal      , only : graph
  use field_stored  , only : stored_field
  use operation_differential, only : differential_operator

  implicit none

  private
  public :: cubic_statement

  !===================================================================!
  ! A nonlinear statement for newton: the chain's rows with
  ! a small cube AGAINST them - a stable reaction, so the jacobian
  ! keeps the rows' own sign and the root is one.
  !===================================================================!

  type, extends(operation) :: cubic_statement

     type(differential_operator) :: linear_part
     real(dp) :: strength = 0.0_dp

   contains

     procedure :: name   => cubic_name
     procedure :: domain => cubic_domain
     procedure :: apply  => cubic_apply

  end type cubic_statement

  interface cubic_statement
     module procedure create_cubic
  end interface cubic_statement

contains

  ! The constructor declares the one argument, the state; the linear
  ! part and the strength are assigned afterwards.
  function create_cubic() result(this)
    type(cubic_statement) :: this
    call this % declare_arguments(1)
  end function create_cubic

  pure function cubic_name(this) result(name)
    class(cubic_statement), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'cubic statement'
  end function cubic_name

  subroutine cubic_domain(this, input_graph, domain, num_entries)
    class(cubic_statement), intent(in)     :: this
    class(directed_graph), intent(in)               :: input_graph
    type(graph), intent(out) :: domain
    integer        , intent(out) :: num_entries
    domain   = input_graph % vertex_set()
    num_entries = input_graph % num_vertices()
  end subroutine cubic_domain

  subroutine cubic_apply(this, input_graph, inputs, output)

    class(cubic_statement), intent(in)             :: this
    class(directed_graph), intent(in)                       :: input_graph
    type(binding), intent(in), optional       :: inputs(:)
    class(field), allocatable, intent(inout) :: output

    type(graph)   :: cells
    type(stored_field)   :: out
    real(dp), allocatable :: q(:), y(:)
    integer :: nv, v

    nv = input_graph % num_vertices()
    if (present(inputs)) then
       call this % linear_part % apply(input_graph, &
            & this % linear_part % bind(inputs), output)
       call output % real_vector(y)
       call bound_real_vector(inputs, this % argument(1), q)
       do v = 1, min(nv, size(q))
          y(v) = y(v) - this % strength * q(v)**3
       end do
    else
       call this % linear_part % apply(input_graph, output=output)
       call output % real_vector(y)
    end if

    cells = input_graph % vertex_set()
    out = stored_field('cubic', cells, nv)
    call out % set_real_vector(y)
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine cubic_apply

end module cubic_statement_fixture

!=====================================================================!
! Fixtures for solver restriction: a test-only solver composition and
! a two-component implicit march whose residual the temporal
! partition constrains to members.
!=====================================================================!

module restriction_fixture

  use iso_fortran_env, only : dp => REAL64
  use graph_fractal, only : graph
  use view_directed, only : directed_graph
  use field_stored, only : stored_field
  use operation_action, only : operation
  use operation_minimization, only : minimizer, state, restrict, solve_result, SOLVE_SINGULAR
  use operation_stencil, only : stencil
  use operation_residual, only : residual_operator
  use operation_expression, only : unknown, derivative, constant, stated, design
  use operation_expression, only : operator(+), operator(*), operator(**)

  implicit none

  private
  public :: delegating_solver, march_residual, designed_march_residual, num_restrictions

  ! how many restrictions the delegating solvers have received
  integer :: num_restrictions = 0

  !===================================================================!
  ! A solver of no metadata of its own, stating and solving through
  ! the inner minimizer it owns. It restricts by delegation alone and
  ! counts the restrictions it receives. With report_failure set its
  ! solve reports a singular result, so an outer solver's propagation
  ! of an unsuccessful inner solve is observable.
  !===================================================================!

  type, extends(minimizer) :: delegating_solver
     class(minimizer), allocatable :: inner
     logical :: report_failure = .false.
   contains
     procedure :: name     => delegating_name
     procedure :: state    => delegating_state
     procedure :: restrict => delegating_restrict
     procedure :: solve    => delegating_solve
  end type delegating_solver

contains

  pure function delegating_name(this) result(name)
    class(delegating_solver), intent(in) :: this
    character(len=:), allocatable :: name
    associate (u1 => this); end associate
    name = 'delegating solver'
  end function delegating_name

  subroutine delegating_state(this, action, context, unknown_domain, num_unknowns, &
       & num_components, coupling, stored_inputs)
    class(delegating_solver), intent(inout) :: this
    class(operation)        , intent(in)    :: action
    class(directed_graph)   , intent(in)    :: context
    type(graph)             , intent(in)    :: unknown_domain
    integer                 , intent(in)    :: num_unknowns
    integer                 , intent(in), optional :: num_components
    class(directed_graph)   , intent(in), optional :: coupling
    type(stored_field)      , intent(in), optional :: stored_inputs(:)
    call state(this, action, context, unknown_domain, num_unknowns, num_components, coupling, stored_inputs)
    call this % inner % state(action, context, unknown_domain, num_unknowns, num_components, coupling, stored_inputs)
  end subroutine delegating_state

  subroutine delegating_restrict(this, selected)
    class(delegating_solver), intent(inout) :: this
    integer                 , intent(in)    :: selected(:)
    call restrict(this, selected)
    num_restrictions = num_restrictions + 1
    if (allocated(this % inner)) call this % inner % restrict(selected)
  end subroutine delegating_restrict

  subroutine delegating_solve(this, rhs, x, achieved)
    class(delegating_solver), intent(inout) :: this
    real(dp), intent(in)    :: rhs(:)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out)   :: achieved
    real(dp), allocatable :: r(:)
    type(solve_result) :: outcome
    call this % initialize_residual_history()
    if (this % report_failure) then
       call this % imbalance(rhs, x, r)
       achieved = this % norm(r)
       call this % record_residual_norm(achieved)
       call this % record_result(achieved, 0, SOLVE_SINGULAR)
       return
    end if
    call this % inner % solve(rhs, x, achieved)
    outcome = this % inner % result()
    call this % record_residual_norm(this % inner % initial_residual_norm())
    call this % record_result(outcome % residual, outcome % iterations, outcome % reason)
  end subroutine delegating_solve

  !===================================================================!
  ! The implicit march q' = -c q^3 over n instants of step h from
  ! q(0) = q0. Instant p stores the tuple (q_p, q'_p) at unknowns
  ! 2p - 1 and 2p. The physics governs the q row, q'_p + c q_p^3, and
  ! the stencil ties the q' row, q'_p - (q_p - q_{p-1}) / h. The first
  ! instant is fixed at (q0, -c q0^3).
  !===================================================================!

  function march_residual(n, h, c, q0) result(residual)
    integer , intent(in) :: n
    real(dp), intent(in) :: h, c, q0
    type(residual_operator) :: residual
    type(stencil) :: tying
    integer , allocatable :: rows(:), columns(:)
    real(dp), allocatable :: weights(:)
    integer :: p, e
    allocate(rows(3 * (n - 1)), columns(3 * (n - 1)), weights(3 * (n - 1)))
    e = 0
    do p = 2, n
       rows(e + 1:e + 3)    = 2 * p
       columns(e + 1:e + 3) = [2 * p, 2 * p - 1, 2 * p - 3]
       weights(e + 1:e + 3) = [1.0_dp, -1.0_dp / h, 1.0_dp / h]
       e = e + 3
    end do
    tying = stencil(rows, columns, weights, spread(0.0_dp, 1, 2 * n), 'tying rows')
    residual = residual_operator(tying, &
         & stated(derivative(unknown(), 1) + constant(c) * derivative(unknown(), 0) ** 3, 1, 'cubic decay'), &
         & [(2 * (p - 1), p = 1, n)], 2 * n, 2, [0], [1, 2], [q0, -c * q0 ** 3])
  end function march_residual

  !===================================================================!
  ! The same march with a design at every point: q' = -c q^3 - nu q,
  ! the physics row q'_p + c q_p^3 + nu_p q_p, so the design tangent,
  ! the mixed state-design partial and the third state partial are
  ! nonzero at every governed point. The first instant is fixed at
  ! (q0, -c q0^3 - nu_1 q0).
  !===================================================================!

  function designed_march_residual(n, h, c, q0, nu_1) result(residual)
    integer , intent(in) :: n
    real(dp), intent(in) :: h, c, q0, nu_1
    type(residual_operator) :: residual
    type(stencil) :: tying
    integer , allocatable :: rows(:), columns(:)
    real(dp), allocatable :: weights(:)
    integer :: p, e
    allocate(rows(3 * (n - 1)), columns(3 * (n - 1)), weights(3 * (n - 1)))
    e = 0
    do p = 2, n
       rows(e + 1:e + 3)    = 2 * p
       columns(e + 1:e + 3) = [2 * p, 2 * p - 1, 2 * p - 3]
       weights(e + 1:e + 3) = [1.0_dp, -1.0_dp / h, 1.0_dp / h]
       e = e + 3
    end do
    tying = stencil(rows, columns, weights, spread(0.0_dp, 1, 2 * n), 'tying rows')
    residual = residual_operator(tying, &
         & stated(derivative(unknown(), 1) + constant(c) * derivative(unknown(), 0) ** 3 &
         &        + design() * derivative(unknown(), 0), 1, 'designed cubic decay'), &
         & [(2 * (p - 1), p = 1, n)], 2 * n, 2, [0], [1, 2], [q0, -c * q0 ** 3 - nu_1 * q0])
  end function designed_march_residual

end module restriction_fixture

program test_graph_minimization

  use iso_fortran_env, only : dp => REAL64
  use graph_fractal  , only : graph
  use map_set_store, only : set_store
  use view_directed, only : SIDE_EDGE, SIDE_VERTEX
  use view_directed_stored, only : stored_directed_graph
  use field_calculus, only : field
  use field_stored, only : stored_field
  use view_mesh   , only : mesh, values_of
  use view_mesh_builder , only : mesh_from_gmsh
  use operation_robin_condition, only : robin_condition, dirichlet, COEFFICIENT_OPERATOR
  use operation_conduction     , only : conduction
  use operation_differential, only : differential_operator
  use operation_differential, only : laplacian, gradient
  use operation_jacobi   , only : jacobi
  use operation_conjugate_gradient, only : conjugate_gradient
  use operation_gauss_seidel, only : gauss_seidel
  use operation_gmres    , only : gmres
  use operation_newton   , only : newton
  use operation_minimization, only : minimizer, solve_result, absolute, SOLVE_EXHAUSTED, SOLVE_BREAKDOWN, &
       & SOLVE_INNER_FAILED, SOLVE_NOT_STARTED, SOLVE_STAGNATED, state_tuple
  use operation_linearization, only : linearization, tangent_of
  use operation_action, only : varied, variation, dense_of_triples
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf, ieee_is_finite, ieee_is_nan
  use operation_dense_direct, only : dense_direct
  use operation_stencil, only : stencil
  use operation_balance  , only : balance
  use operation_elimination, only : elimination
  use operation_multigrid, only : multigrid
  use operation_temporal_minimization, only : temporal_minimizer
  use operation_residual, only : residual_operator
  use operation_domain, only : continuous_domain, discrete_domain
  use operation_expression, only : expression, unknown, derivative, constant, stated, &
       & operator(+), operator(*), operator(**)
  use field_stored, only : typed_field_domain
  use cubic_statement_fixture, only : cubic_statement
  use restriction_fixture, only : delegating_solver, march_residual, designed_march_residual, num_restrictions

  implicit none

  integer :: num_failures

  num_failures = 0

  call check_minimizer_operations(num_failures)
  call check_residual_arithmetic(num_failures)
  call check_three_cell_solution(num_failures)
  call check_real_mesh(num_failures)
  call check_sweeping_family(num_failures)
  call check_gmres_family(num_failures)
  call check_newton(num_failures)
  call check_iteration_limits(num_failures)
  call check_scaled_gmres(num_failures)
  call check_gmres_exit(num_failures)
  call check_direction_scale(num_failures)
  call check_restriction_maps(num_failures)
  call check_restricted_operator(num_failures)
  call check_solver_restriction(num_failures)
  call check_residual_boundary(num_failures)
  call check_one_law_two_placements(num_failures)

  write(*, '(a)') ' ============================================='
  if (num_failures == 0) then
     write(*, '(a)') ' all minimization checks passed'
  else
     write(*, '(a, i0, a)') ' ', num_failures, ' minimization checks FAILED'
     error stop 1
  end if

contains

  subroutine check_residual_arithmetic(num_failures)
    integer, intent(inout) :: num_failures
    type(jacobi) :: solver
    real(dp) :: least, overflowing
    call report(abs(solver % norm([3.0e200_dp, 4.0e200_dp]) / 5.0e200_dp - 1.0_dp) < 1.0e-14_dp, &
         & 'the residual norm avoids overflow in its sum of squares', num_failures)
    call report(abs(solver % norm([3.0e-200_dp, 4.0e-200_dp]) / 5.0e-200_dp - 1.0_dp) < 1.0e-14_dp, &
         & 'the residual norm avoids underflow in its sum of squares', num_failures)
    call report(abs(solver % norm([3.0e200_dp, 4.0e200_dp, tiny(1.0_dp), -1.0e-200_dp]) / &
         & 5.0e200_dp - 1.0_dp) < 1.0e-14_dp, &
         & 'mixed huge and tiny components do not underflow while scaling the residual norm', num_failures)
    call report(abs(solver % norm([huge(1.0_dp) / 2.0_dp, huge(1.0_dp) / 2.0_dp, tiny(1.0_dp)]) / &
         & (huge(1.0_dp) / 2.0_dp) - sqrt(2.0_dp)) < 1.0e-14_dp, &
         & 'a representable residual norm near the largest number stays finite', num_failures)
    overflowing = solver % norm([huge(1.0_dp), huge(1.0_dp)])
    call report(.not. ieee_is_finite(overflowing) .and. .not. ieee_is_nan(overflowing), &
         & 'a genuinely overflowing residual norm is reported as positive infinity', num_failures)
    least = scale(1.0_dp, minexponent(1.0_dp) - digits(1.0_dp))
    call report(solver % norm([least, least]) == least, &
         & 'a subnormal residual norm rounds without trapping underflow', num_failures)
    call solver % initialize_residual_history()
    call solver % record_residual_norm(1.0_dp)
    solver % tolerance = 0.25_dp
    call report(solver % converged(0.25_dp) .and. .not. solver % converged(nearest(0.25_dp, 1.0_dp)), &
         & 'relative convergence includes its exact boundary and excludes the next number', num_failures)
    call report(.not. solver % converged(ieee_value(0.0_dp, ieee_positive_inf)), &
         & 'an infinite residual is never converged', num_failures)
    call solver % initialize_residual_history()
    call solver % record_residual_norm(scale(1.0_dp, -1000))
    solver % tolerance = scale(1.0_dp, -1000)
    call report(solver % converged(0.0_dp) .and. .not. solver % converged(tiny(1.0_dp)), &
         & 'a relative threshold below the arithmetic range is compared without multiplying it', num_failures)
  end subroutine check_residual_arithmetic

  subroutine check_iteration_limits(num_failures)
    integer, intent(inout) :: num_failures
    class(minimizer), allocatable :: solver
    type(newton) :: nonlinear
    type(stencil) :: identity
    type(solve_result) :: outcome
    real(dp) :: x(2), rhs(2), achieved
    integer :: member

    identity = stencil([1, 2], [1, 2], [1.0_dp, 1.0_dp], [0.0_dp, 0.0_dp], 'identity')
    rhs = [1.0_dp, 2.0_dp]
    allocate(nonlinear % inner, source=dense_direct())
    do member = 1, 5
       select case (member)
       case (1)
          allocate(solver, source=jacobi())
       case (2)
          allocate(solver, source=gauss_seidel())
       case (3)
          allocate(solver, source=conjugate_gradient())
       case (4)
          allocate(solver, source=gmres())
       case (5)
          allocate(solver, source=nonlinear)
       end select
       solver % max_iterations = 1
       solver % tolerance = 1.0e-12_dp
       call solver % state(identity, identity % pattern, identity % pattern % vertex_set(), 2, &
            & coupling=identity % pattern)
       x = 0.0_dp
       call solver % solve(rhs, x, achieved)
       outcome = solver % result()
       call report(maxval(abs(x - rhs)) < 1.0e-12_dp .and. outcome % converged() .and. &
            & outcome % iterations == 1, solver % name() // ': one permitted step is completed', num_failures)
       x = rhs
       call solver % solve(rhs, x, achieved)
       outcome = solver % result()
       call report(achieved == 0.0_dp .and. outcome % converged() .and. outcome % iterations == 0, &
            & solver % name() // ': an initially solved system takes no step', num_failures)
       call solver % state(identity, identity % pattern, identity % pattern % vertex_set(), 2, &
            & coupling=identity % pattern)
       outcome = solver % result()
       call report(outcome % reason == SOLVE_NOT_STARTED, &
            & solver % name() // ': stating another system clears the previous outcome', num_failures)
       solver % max_iterations = 0
       x = 0.0_dp
       call solver % solve(rhs, x, achieved)
       outcome = solver % result()
       call report(all(x == 0.0_dp) .and. outcome % reason == SOLVE_EXHAUSTED .and. &
            & outcome % iterations == 0, solver % name() // ': a zero iteration limit performs no update', num_failures)
       deallocate(solver)
    end do
  end subroutine check_iteration_limits

  subroutine check_scaled_gmres(num_failures)
    integer, intent(inout) :: num_failures
    type(stencil) :: a
    type(gmres) :: solver
    type(jacobi) :: preconditioner
    type(solve_result) :: outcome
    real(dp) :: x(2), rhs(2), achieved, scale
    integer :: exponent

    preconditioner % max_iterations = 1
    allocate(solver % preconditioner, source=preconditioner)
    solver % restart = 2
    solver % max_iterations = 20
    solver % tolerance = 1.0e-8_dp
    rhs = [1.0_dp, 0.0_dp]
    do exponent = -12, 12, 12
       scale = 10.0_dp ** exponent
       a = stencil([1, 1, 2, 2], [1, 2, 1, 2], scale * [2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], &
            & [0.0_dp, 0.0_dp], 'scaled system')
       call solver % state(a, a % pattern, a % pattern % vertex_set(), 2, coupling=a % pattern)
       x = 0.0_dp
       call solver % solve(rhs, x, achieved)
       outcome = solver % result()
       call report(achieved <= solver % tolerance .and. outcome % converged(), &
            & 'GMRES verifies the true residual under scalar preconditioning', num_failures)
    end do
  end subroutine check_scaled_gmres

  subroutine check_gmres_exit(num_failures)
    integer, intent(inout) :: num_failures
    type(stencil) :: a
    type(gmres) :: solver
    type(dense_direct) :: singular
    type(jacobi) :: stopped
    type(solve_result) :: outcome
    real(dp) :: x(2), rhs(2), achieved, true_residual_norm

    a = stencil([1, 1, 2, 2], [1, 2, 1, 2], [2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], &
         & [0.0_dp, 0.0_dp], 'one restart')
    solver % restart = 1
    solver % max_iterations = 1
    solver % tolerance = 1.0e-12_dp
    call solver % state(a, a % pattern, a % pattern % vertex_set(), 2, coupling=a % pattern)
    rhs = [1.0_dp, 0.0_dp]
    x = 0.0_dp
    call solver % solve(rhs, x, achieved)
    outcome = solver % result()
    true_residual_norm = norm2(rhs - [2.0_dp * x(1) + x(2), x(1) + 2.0_dp * x(2)])
    call report(any(x /= 0.0_dp) .and. abs(achieved - true_residual_norm) < 1.0e-14_dp .and. &
         & outcome % reason == SOLVE_EXHAUSTED .and. outcome % iterations == 1, &
         & 'exhausted GMRES reports the residual after its permitted restart', num_failures)

    a = stencil([1, 2], [1, 2], [0.0_dp, 0.0_dp], [0.0_dp, 0.0_dp], 'zero operator')
    call solver % state(a, a % pattern, a % pattern % vertex_set(), 2, coupling=a % pattern)
    x = 0.0_dp
    call solver % solve(rhs, x, achieved)
    outcome = solver % result()
    call report(all(x == 0.0_dp) .and. achieved == 1.0_dp .and. &
         & outcome % reason == SOLVE_BREAKDOWN, 'GMRES reports a zero Arnoldi column without dividing by zero', num_failures)

    stopped % max_iterations = 0
    allocate(solver % preconditioner, source=stopped)
    call solver % state(a, a % pattern, a % pattern % vertex_set(), 2, coupling=a % pattern)
    call solver % solve(rhs, x, achieved)
    outcome = solver % result()
    call report(all(x == 0.0_dp) .and. achieved == 1.0_dp .and. &
         & outcome % reason == SOLVE_BREAKDOWN, 'GMRES reports an annihilated residual without dividing by zero', num_failures)
    deallocate(solver % preconditioner)

    singular % singular_reported = .true.
    allocate(solver % preconditioner, source=singular)
    call solver % state(a, a % pattern, a % pattern % vertex_set(), 2, coupling=a % pattern)
    call solver % solve(rhs, x, achieved)
    outcome = solver % result()
    call report(all(x == 0.0_dp) .and. achieved == 1.0_dp .and. &
         & outcome % reason == SOLVE_INNER_FAILED, 'GMRES propagates a singular preconditioner', num_failures)
  end subroutine check_gmres_exit

  subroutine check_direction_scale(num_failures)
    integer, intent(inout) :: num_failures
    type(stored_directed_graph) :: single_vertex
    type(cubic_statement) :: action
    type(linearization) :: tangent
    type(stored_field) :: state, direction
    class(field), allocatable :: image
    type(gmres) :: solver
    type(solve_result) :: outcome
    real(dp), allocatable :: value(:)
    real(dp) :: unit, factor, x(1), achieved
    integer :: exponent

    single_vertex = stored_directed_graph(1, tails=[integer ::], heads=[integer ::])
    action = cubic_statement()
    action % linear_part = differential_operator(SIDE_VERTEX, 0, coefficient=2.0_dp)
    state = stored_field('state', single_vertex % vertex_set(), 1)
    call state % set_real_vector([1.0_dp])
    tangent = tangent_of(action, action % argument(1))
    call tangent % freeze([state])
    direction = stored_field('direction', single_vertex % vertex_set(), 1)
    call direction % set_real_vector([1.0_dp])
    call tangent % apply(single_vertex, tangent % bind([direction]), image)
    call image % real_vector(value)
    unit = value(1)
    do exponent = -200, 200, 100
       factor = 10.0_dp ** exponent
       call direction % set_real_vector([factor])
       call tangent % apply(single_vertex, tangent % bind([direction]), image)
       call image % real_vector(value)
       call report(abs(value(1) / factor - unit) < 1.0e-14_dp, &
            & 'a scalar directional difference is homogeneous across direction scales', num_failures)
    end do
    call direction % set_real_vector([0.0_dp])
    call tangent % apply(single_vertex, tangent % bind([direction]), image)
    call image % real_vector(value)
    call report(value(1) == 0.0_dp, 'a zero direction has an exact zero numerical derivative', num_failures)

    ! Opposite directions can disagree at a numerical difference's
    ! truncation scale. Exhausting the one-dimensional Krylov space
    ! reports that mismatch without accepting its estimated residual.
    action % strength = 1.0_dp
    tangent = tangent_of(action, action % argument(1))
    call tangent % freeze([state])
    call solver % state(tangent, single_vertex, single_vertex % vertex_set(), 1)
    solver % tolerance = 1.0e-14_dp
    x = 0.0_dp
    call solver % solve([1.0_dp], x, achieved)
    outcome = solver % result()
    call report(outcome % reason == SOLVE_STAGNATED .and. .not. outcome % converged() .and. &
         & achieved > solver % tolerance .and. outcome % iterations == 1, &
         & 'an invariant Krylov space reports a numerical derivative mismatch as stagnation', num_failures)
  end subroutine check_direction_scale

  subroutine report(passed, message, num_failures)

    logical         , intent(in)    :: passed
    character(len=*), intent(in)    :: message
    integer         , intent(inout) :: num_failures

    if (passed) then
       write(*, '(a)') ' PASS : ' // message
    else
       write(*, '(a)') ' FAIL : ' // message
       num_failures = num_failures + 1
    end if

  end subroutine report

  !===================================================================!
  ! The three-cell conduction problem, assembled from the levels.
  !===================================================================!

  type(mesh) function chain_mesh() result(m)

    m = mesh(3, tails=[1, 2, 1, 3], heads=[2, 3, 0, 0], &
         & volumes      = [1.0_dp, 1.0_dp, 1.0_dp], &
         & cell_centres = [0.5_dp, 0.0_dp, 0.0_dp, &
         &                 1.5_dp, 0.0_dp, 0.0_dp, &
         &                 2.5_dp, 0.0_dp, 0.0_dp], &
         & areas        = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
         & deltas       = [1.0_dp, 1.0_dp, 0.5_dp, 0.5_dp], &
         & normals      = [ 1.0_dp, 0.0_dp, 0.0_dp, &
         &                  1.0_dp, 0.0_dp, 0.0_dp, &
         &                 -1.0_dp, 0.0_dp, 0.0_dp, &
         &                  1.0_dp, 0.0_dp, 0.0_dp], &
         & face_centres = [1.0_dp, 0.0_dp, 0.0_dp, &
         &                 2.0_dp, 0.0_dp, 0.0_dp, &
         &                 0.0_dp, 0.0_dp, 0.0_dp, &
         &                 3.0_dp, 0.0_dp, 0.0_dp], &
         & weights      = [0.5_dp, 0.5_dp, 1.0_dp, 1.0_dp], &
         & etags        = [character(len=4) :: '', '', 'west', 'east'])

  end function chain_mesh

  !===================================================================!
  ! The assembly: the constitution's coefficients placed into the
  ! operator's argument list - conduction inside, one condition per
  ! boundary - and nothing else.
  !===================================================================!

  type(differential_operator) function chain_operator(m) result(op)

    type(mesh), intent(in) :: m

    real(dp), allocatable :: c(:), b(:)

    call assemble(m, conduction(1.0_dp), &
         & [dirichlet('west', 0.0_dp), dirichlet('east', 10.0_dp)], &
         & 1.0_dp, c, b)

    op = laplacian(coefficients=c, &
         & spacings=[1.0_dp, 1.0_dp, 0.5_dp, 0.5_dp], boundary_values=b)

  end function chain_operator

  subroutine assemble(m, law, conditions, kappa, c, b)

    type(mesh), intent(in)            :: m
    type(conduction), intent(in)      :: law
    type(robin_condition), intent(in) :: conditions(:)
    real(dp), intent(in)              :: kappa
    real(dp), allocatable, intent(out) :: c(:), b(:)

    type(graph)       :: members
    type(set_store)       :: sets
    integer , allocatable :: face(:)
    real(dp), allocatable :: cw(:), bw(:)
    integer :: k, f, e

    ! the headless faces enter through the conditions, so the
    ! conduction coefficient is zero there
    call law % edge_coefficients(m, .true., c)
    allocate(b(size(c)))
    b = 0.0_dp

    do k = 1, size(conditions)
       call conditions(k) % faces(m, sets, members)
       call conditions(k) % coefficient_values(m, kappa, COEFFICIENT_OPERATOR, cw)
       call conditions(k) % boundary_values(m, bw)
       call sets % members_of(members, face)
       do f = 1, size(face)
          e = face(f)
          c(e) = cw(f)
          b(e) = bw(f)
       end do
    end do

  end subroutine assemble

  !===================================================================!
  ! The solver's words, checked against hand values on the chain:
  ! the diagonal by coloured indicators is the stencil's own, and the colouring
  ! never gives neighbours one colour.
  !===================================================================!

  subroutine check_minimizer_operations(num_failures)

    integer, intent(inout) :: num_failures

    type(mesh) :: m
    type(jacobi) :: js
    real(dp), allocatable :: d(:,:,:)
    integer , allocatable :: colours(:), nbrs(:)
    logical :: proper
    integer :: v, i

    m = chain_mesh()
    call js % state(chain_operator(m), m, m % vertex_set(), m % num_vertices(), coupling = m)

    call js % block_diagonal(d)
    call report(all(abs(d(1, 1, :) - [-3.0_dp, -2.0_dp, -3.0_dp]) < 1.0d-12), &
         & 'the diagonal by coloured indicators is the stencil diagonal, by hand', num_failures)

    call js % sweep_order(colours)
    proper = .true.
    do v = 1, m % num_vertices()
       call m % adjacent_vertices(v, nbrs)
       do i = 1, size(nbrs)
          if (colours(nbrs(i)) == colours(v)) proper = .false.
       end do
    end do
    call report(proper, 'the sweep order never gives neighbours one colour', num_failures)

    call report(abs(js % inner_product([1.0_dp, 2.0_dp, 3.0_dp], &
         &                             [2.0_dp, 1.0_dp, 0.0_dp]) - 4.0_dp) &
         & < 1.0d-14, 'the inner product is the measured sum', num_failures)

    call report(abs(js % norm([3.0_dp, 4.0_dp, 0.0_dp]) - 5.0_dp) < 1.0d-14, &
         & 'the Euclidean norm of (3,4,0) equals five', num_failures)

  end subroutine check_minimizer_operations

  !===================================================================!
  ! The three-cell diffusion equation with an exact discrete solution.
  !===================================================================!

  subroutine check_three_cell_solution(num_failures)

    integer, intent(inout) :: num_failures

    type(mesh) :: m
    type(jacobi) :: js
    real(dp), allocatable :: g(:), rhs(:), x(:), y(:)
    real(dp) :: achieved
    real(dp), parameter :: exact(3) = [5.0_dp/3.0_dp, 5.0_dp, 25.0_dp/3.0_dp]

    m = chain_mesh()
    call js % state(chain_operator(m), m, m % vertex_set(), m % num_vertices(), coupling = m)
    js % max_iterations = 5000
    js % tolerance      = 1.0d-12

    ! The equation action(q) = 0 reads matvec(q) = -affine.
    g   = js % affine
    rhs = -g

    allocate(x(3))
    x = 0.0_dp
    call js % solve(rhs, x, achieved)

    call report(achieved < 1.0d-10, 'jacobi drives the residual down', num_failures)
    call report(all(abs(x - exact) < 1.0d-8), &
         & 'the discrete solution equals 5/3, 5, 25/3', num_failures)

    call js % matvec(x, y)
    call report(js % norm(rhs - y) < 1.0d-10, &
         & 'the solution satisfies the assembled equation', num_failures)

  end subroutine check_three_cell_solution

  !===================================================================!
  ! The real mesh. A symmetric positive operator - the conduction
  ! coefficients negated - and a consistent right hand side
  ! manufactured from a known state; conjugate gradient must return
  ! a state the operator cannot tell apart from it.
  !===================================================================!

  subroutine check_real_mesh(num_failures)

    integer, intent(inout) :: num_failures

    type(mesh) :: m
    type(conjugate_gradient) :: cg
    type(conduction) :: law
    real(dp), allocatable :: c(:), deltas(:), xref(:), b(:), x(:), y(:)
    real(dp) :: achieved
    integer :: v

    m = mesh_from_gmsh('../square-10.msh')

    law = conduction(1.0_dp)
    call law % edge_coefficients(m, .true., c)
    c = -c

    call values_of(m % face_delta(), deltas)

    call cg % state(laplacian(coefficients=c, spacings=deltas), m, m % vertex_set(), &
         & m % num_vertices())
    cg % max_iterations = 2000
    cg % tolerance      = 1.0d-10

    allocate(xref(m % num_vertices()))
    do v = 1, size(xref)
       xref(v) = real(mod(3 * v, 17), dp)
    end do

    call cg % matvec(xref, b)

    allocate(x(size(xref)))
    x = 0.0_dp
    call cg % solve(b, x, achieved)

    call cg % matvec(x, y)
    call report(cg % norm(b - y) < 1.0d-7 * (1.0_dp + cg % norm(b)), &
         & 'conjugate gradient closes the manufactured equation on a real mesh', &
         & num_failures)

  end subroutine check_real_mesh

  !===================================================================!
  ! Gauss-seidel and its omega: the same chain, the same exact
  ! solution, the sweep ordered by colour. SOR is not another type -
  ! it is this one at omega away from one.
  !===================================================================!

  subroutine check_sweeping_family(num_failures)

    integer, intent(inout) :: num_failures

    type(mesh) :: m
    type(gauss_seidel) :: gs
    real(dp), allocatable :: g(:), rhs(:), x(:)
    real(dp) :: achieved
    real(dp), parameter :: exact(3) = [5.0_dp/3.0_dp, 5.0_dp, 25.0_dp/3.0_dp]

    m = chain_mesh()
    call gs % state(chain_operator(m), m, m % vertex_set(), m % num_vertices(), coupling = m)
    gs % max_iterations = 2000
    gs % tolerance      = 1.0d-12

    g   = gs % affine
    rhs = -g

    allocate(x(3))
    x = 0.0_dp
    call gs % solve(rhs, x, achieved)
    call report(all(abs(x - exact) < 1.0d-8), &
         & 'gauss-seidel sweeps by colour to the exact solution', num_failures)

    gs % omega = 1.2_dp
    x = 0.0_dp
    call gs % solve(rhs, x, achieved)
    call report(all(abs(x - exact) < 1.0d-8), &
         & 'and over-relaxed it is sor, a parameter, not a type', num_failures)

  end subroutine check_sweeping_family

  !===================================================================!
  ! GMRES: the chain again, then an unsymmetric statement -
  ! advection added to diffusion through the balance - where two
  ! different solvers must meet on one solution.
  !===================================================================!

  subroutine check_gmres_family(num_failures)

    integer, intent(inout) :: num_failures

    type(mesh) :: m
    type(gmres)  :: gm
    type(jacobi) :: js
    type(balance) :: statement
    real(dp), allocatable :: g(:), rhs(:), x(:), xj(:)
    real(dp) :: achieved
    real(dp), parameter :: exact(3) = [5.0_dp/3.0_dp, 5.0_dp, 25.0_dp/3.0_dp]

    m = chain_mesh()

    call gm % state(chain_operator(m), m, m % vertex_set(), m % num_vertices())
    gm % tolerance = 1.0d-12
    g   = gm % affine
    rhs = -g
    allocate(x(3))
    x = 0.0_dp
    call gm % solve(rhs, x, achieved)
    call report(all(abs(x - exact) < 1.0d-8), &
         & 'gmres lands on the exact solution', num_failures)

    ! Diffusion and upwind advection in one balance: unsymmetric.
    statement = balance(edge_terms=[ &
         & gradient(coefficients=[1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp], &
         &      spacings=[1.0_dp, 1.0_dp, 0.5_dp, 0.5_dp], &
         &      boundary_values=[0.0_dp, 0.0_dp, 0.0_dp, 10.0_dp]), &
         & differential_operator(SIDE_EDGE, 0, &
         &      coefficients=[0.4_dp, 0.4_dp, 0.0_dp, 0.0_dp], &
         &      one_sided=.true.)])

    call gm % state(statement, m, m % vertex_set(), m % num_vertices())
    g   = gm % affine
    rhs = -g
    x = 0.0_dp
    call gm % solve(rhs, x, achieved)
    call report(achieved < 1.0d-10, &
         & 'gmres closes the unsymmetric statement', num_failures)

    call js % state(statement, m, m % vertex_set(), m % num_vertices(), coupling = m)
    js % max_iterations = 5000
    js % tolerance      = 1.0d-12
    allocate(xj(3))
    xj = 0.0_dp
    call js % solve(rhs, xj, achieved)
    call report(all(abs(x - xj) < 1.0d-7), &
         & 'and two different solvers meet on one solution', num_failures)

  end subroutine check_gmres_family

  !===================================================================!
  ! Newton over the linear family. A cubic term is added to the chain
  ! statement; newton linearizes by directional differences, passes
  ! each linear system to the gmres it governs, and drives the
  ! nonlinear residual to zero.
  !===================================================================!

  subroutine check_newton(num_failures)

    integer, intent(inout) :: num_failures

    type(mesh) :: m
    type(newton) :: ns
    type(cubic_statement) :: action
    real(dp), allocatable :: q(:)
    real(dp) :: achieved

    m = chain_mesh()

    action = cubic_statement()
    action % linear_part = chain_operator(m)
    action % strength    = 0.05_dp

    allocate(ns % inner, source=gmres())
    ns % inner % tolerance = 1.0d-12
    ns % tolerance = 1.0d-7
    ns % criterion = absolute

    call ns % state(action, m, m % vertex_set(), m % num_vertices())

    allocate(q(3))
    q = 0.0_dp
    call ns % solve([0.0_dp, 0.0_dp, 0.0_dp], q, achieved)

    call report(achieved < 1.0d-7, &
         & 'newton drives the nonlinear residual to the difference floor', num_failures)
    call report(q(1) > 0.0_dp .and. q(1) < q(2) .and. q(2) < q(3), &
         & 'and the heated chain still rises monotonically', num_failures)

  end subroutine check_newton

  !===================================================================!
  ! SOLVER RESTRICTION. A nested composition over the march of n
  ! instants, each instant a tuple (q, q'):
  !
  !   1  newton > elimination of q' > gmres > multigrid preconditioner
  !      over the retained q, aggregated by instant pairs
  !   2  newton > gmres > elimination preconditioner > dense direct
  !   3  newton > delegating solver > gmres > multigrid preconditioner
  !      over every unknown in blocks of two, aggregated by instant pairs
  !   4  newton > elimination of q' > delegating solver reporting
  !      failure > dense direct
  !   5  newton > elimination of q' within one entry of storage >
  !      dense direct
  !
  ! The temporal engine dispatches on none of these.
  !===================================================================!

  subroutine composed(kind, n, template)

    integer, intent(in) :: kind, n
    class(minimizer), allocatable, intent(out) :: template

    type(newton)       :: ns
    type(elimination)  :: schur
    type(gmres)        :: krylov
    type(multigrid)    :: levels
    type(jacobi)       :: smoother
    type(gauss_seidel) :: sweeps
    type(dense_direct) :: factorisation
    type(delegating_solver) :: delegate
    integer :: p

    ns % tolerance      = 1.0e-11_dp
    ns % criterion      = absolute
    ns % max_iterations = 30
    krylov % tolerance      = 1.0e-13_dp
    krylov % max_iterations = 50
    levels % max_iterations = 1
    factorisation % singular_reported = .true.

    select case (kind)
    case (1)
       smoother % max_iterations = 2
       allocate(levels % smoother, source=smoother)
       allocate(levels % coarse, source=factorisation)
       levels % aggregates = [((p + 1) / 2, p = 1, n)]
       allocate(krylov % preconditioner, source=levels)
       schur % eliminated = [(mod(p, 2) == 0, p = 1, 2 * n)]
       allocate(schur % inner, source=krylov)
       allocate(ns % inner, source=schur)
    case (2)
       schur % eliminated = [(mod(p, 2) == 0, p = 1, 2 * n)]
       allocate(schur % inner, source=factorisation)
       allocate(krylov % preconditioner, source=schur)
       allocate(ns % inner, source=krylov)
    case (3)
       sweeps % max_iterations = 2
       sweeps % block_width    = 2
       allocate(levels % smoother, source=sweeps)
       allocate(levels % coarse, source=factorisation)
       levels % block_width = 2
       levels % aggregates  = [(((p - 1) / 2) / 2 + 1, p = 1, 2 * n)]
       allocate(krylov % preconditioner, source=levels)
       allocate(delegate % inner, source=krylov)
       allocate(ns % inner, source=delegate)
    case (4)
       delegate % report_failure = .true.
       allocate(delegate % inner, source=factorisation)
       schur % eliminated = [(mod(p, 2) == 0, p = 1, 2 * n)]
       allocate(schur % inner, source=delegate)
       allocate(ns % inner, source=schur)
    case (5)
       schur % eliminated = [(mod(p, 2) == 0, p = 1, 2 * n)]
       schur % max_entries = 1
       allocate(schur % inner, source=factorisation)
       allocate(ns % inner, source=schur)
    case default
       error stop 'test: a composition is one of the five'
    end select
    allocate(template, source=ns)

  end subroutine composed

  ! whether a composition's template still holds its whole-domain
  ! metadata at every level
  logical function template_intact(kind, n, template) result(intact)

    integer, intent(in) :: kind, n
    class(minimizer), intent(in) :: template

    integer :: p

    intact = .false.
    select type (template)
    type is (newton)
       select type (inner => template % inner)
       type is (elimination)
          intact = flags_whole(inner, n)
          if (kind == 1) then
             select type (krylov => inner % inner)
             type is (gmres)
                select type (levels => krylov % preconditioner)
                type is (multigrid)
                   intact = intact .and. size(levels % aggregates) == n
                   if (intact) intact = all(levels % aggregates == [((p + 1) / 2, p = 1, n)])
                end select
             end select
          end if
       type is (gmres)
          select type (schur => inner % preconditioner)
          type is (elimination)
             intact = flags_whole(schur, n)
          end select
       type is (delegating_solver)
          select type (krylov => inner % inner)
          type is (gmres)
             select type (levels => krylov % preconditioner)
             type is (multigrid)
                intact = size(levels % aggregates) == 2 * n .and. levels % block_width == 2
                if (intact) intact = all(levels % aggregates == [(((p - 1) / 2) / 2 + 1, p = 1, 2 * n)])
             end select
          end select
       end select
    end select

  end function template_intact

  ! whether an elimination still flags the q' rows of n instants
  logical function flags_whole(schur, n)
    type(elimination), intent(in) :: schur
    integer          , intent(in) :: n
    integer :: p
    flags_whole = size(schur % eliminated) == 2 * n
    if (flags_whole) flags_whole = all(schur % eliminated .eqv. [(mod(p, 2) == 0, p = 1, 2 * n)])
  end function flags_whole

  ! the state the march solves begin from: q0 at every instant, q'
  ! zero, the fixed first instant at its values
  subroutine seeded(residual, n, q0, x)
    type(residual_operator), intent(in) :: residual
    integer , intent(in) :: n
    real(dp), intent(in) :: q0
    real(dp), allocatable, intent(out) :: x(:)
    allocate(x(2 * n), source=0.0_dp)
    x(1::2) = q0
    x(residual % fixed_unknowns()) = residual % fixed_values()
  end subroutine seeded

  ! the partition of 2n unknowns: instants in order, instants
  ! reversed, or the parity members reversed
  subroutine labels(layout, n, member_of, member_order)
    integer, intent(in) :: layout, n
    integer, allocatable, intent(out) :: member_of(:), member_order(:)
    integer :: p
    select case (layout)
    case (1)
       member_of    = [((p + 1) / 2, p = 1, 2 * n)]
       member_order = [(p, p = 1, n)]
    case (2)
       member_of    = [((p + 1) / 2, p = 1, 2 * n)]
       member_order = [(n + 1 - p, p = 1, n)]
    case default
       member_of    = [(2 - mod((p + 1) / 2, 2), p = 1, 2 * n)]
       member_order = [2, 1]
    end select
  end subroutine labels

  !===================================================================!
  ! Each solver maps its own metadata through a reordered,
  ! noncontiguous selection and induces its children's selections:
  ! instants 4, 2 and 6 of six, in that order.
  !===================================================================!

  subroutine check_restriction_maps(num_failures)

    integer, intent(inout) :: num_failures

    integer, parameter :: n = 6
    integer, parameter :: selected(6) = [7, 8, 3, 4, 11, 12]
    class(minimizer), allocatable :: template, copy
    type(temporal_minimizer) :: solver
    logical :: mapped
    integer :: p

    call composed(1, n, template)
    select type (template)
    type is (newton)
       select type (schur => template % inner)
       type is (elimination)
          schur % max_entries = 12345
       end select
    end select
    allocate(copy, source=template)
    call copy % restrict(selected)
    mapped = .false.
    select type (copy)
    type is (newton)
       select type (schur => copy % inner)
       type is (elimination)
          mapped = size(schur % eliminated) == 6 .and. schur % max_entries == 12345
          if (mapped) mapped = all(schur % eliminated .eqv. [.false., .true., .false., .true., .false., .true.])
          select type (krylov => schur % inner)
          type is (gmres)
             select type (levels => krylov % preconditioner)
             type is (multigrid)
                ! retained positions 4, 2, 6 of the whole; their
                ! aggregates 2, 1, 3 relabelled in order of appearance
                mapped = mapped .and. size(levels % aggregates) == 3
                if (mapped) mapped = all(levels % aggregates == [1, 2, 3])
             end select
          end select
       end select
    end select
    call report(mapped, 'elimination flags and storage limit follow the selection and the multigrid below it &
         &coarsens by the retained members', num_failures)
    call report(copy % num_unknowns == 0 .and. .not. allocated(copy % action), &
         & 'a restricted solver discards the whole statement until stated on the selection', num_failures)
    call report(template_intact(1, n, template) .and. template % num_unknowns == 0, &
         & 'the whole template is unchanged by restricting a copy', num_failures)
    deallocate(copy, template)

    num_restrictions = 0
    call composed(3, n, template)
    allocate(copy, source=template)
    call copy % restrict(selected)
    mapped = .false.
    select type (copy)
    type is (newton)
       select type (delegate => copy % inner)
       type is (delegating_solver)
          select type (krylov => delegate % inner)
          type is (gmres)
             select type (levels => krylov % preconditioner)
             type is (multigrid)
                mapped = size(levels % aggregates) == 6 .and. levels % block_width == 2
                if (mapped) mapped = all(levels % aggregates == [1, 1, 2, 2, 3, 3])
             end select
          end select
       end select
    end select
    call report(mapped .and. num_restrictions == 1, &
         & 'a test-only solver delegates restriction and the multigrid below it keeps its block layout', &
         & num_failures)
    call report(template_intact(3, n, template), 'the test-only template is unchanged', num_failures)
    deallocate(copy, template)

    ! the temporal minimizer's own partition: instants in reverse
    ! order, restricted to instants 4, 2, 6
    call solver % partition([((p + 1) / 2, p = 1, 2 * n)], [(n + 1 - p, p = 1, n)])
    call solver % restrict(selected)
    mapped = size(solver % member_of) == 6 .and. size(solver % member_order) == 3
    if (mapped) mapped = all(solver % member_of == [1, 1, 2, 2, 3, 3]) .and. all(solver % member_order == [3, 1, 2])
    call report(mapped, 'a partition restricts to the members the selection meets, in the stated order', &
         & num_failures)

  end subroutine check_restriction_maps

  !===================================================================!
  ! The operator on a selection: the residual constrained to instants
  ! 4, 2 and 6 with the exterior fixed at the state agrees with the
  ! whole residual on those rows, its explicit tangent is the whole
  ! tangent's submatrix, the linearization action reproduces that
  ! assembled tangent, and the transposed stencil's action is its
  ! adjoint under the Euclidean pairing.
  !===================================================================!

  subroutine check_restricted_operator(num_failures)

    integer, intent(inout) :: num_failures

    integer , parameter :: n = 6, k = 6
    integer , parameter :: selected(k) = [7, 8, 3, 4, 11, 12]
    real(dp), parameter :: h = 0.1_dp, c = 0.5_dp, q0 = 1.0_dp
    type(residual_operator) :: residual, sub
    type(stored_directed_graph) :: unknowns, members
    type(stored_field), allocatable :: inputs(:), inputs_sub(:)
    type(linearization) :: tangent
    type(stencil) :: assembled
    type(dense_direct) :: action, adjoint
    class(field), allocatable :: image
    integer , allocatable :: rows(:), columns(:)
    real(dp), allocatable :: weights(:), x(:), r_whole(:), r_sub(:), j_whole(:,:), j_sub(:,:), y(:)
    real(dp), allocatable :: basis(:), u(:), v(:), ju(:), jtv(:)
    real(dp) :: scale, difference
    logical :: tangent_defined
    integer :: p, i, j

    residual = march_residual(n, h, c, q0)
    unknowns = residual % unknown_graph()
    allocate(x(2 * n))
    do p = 1, n
       x(2 * p - 1) = q0 * 0.9_dp ** p
       x(2 * p)     = -0.3_dp / real(p, dp)
    end do
    inputs = residual % frozen_tuple(x, spread(0.0_dp, 1, n))

    call residual % apply(unknowns, residual % bind(inputs), image)
    call image % real_vector(r_whole)
    call residual % explicit_tangent(unknowns, residual % bind(inputs), 1, rows, columns, weights, tangent_defined)
    call report(tangent_defined, 'the march residual states an explicit tangent', num_failures)
    allocate(j_whole(2 * n, 2 * n), source=0.0_dp)
    do i = 1, size(rows)
       j_whole(rows(i), columns(i)) = j_whole(rows(i), columns(i)) + weights(i)
    end do

    sub     = residual % constrain(selected, x)
    members = sub % unknown_graph()
    inputs_sub = sub % frozen_tuple(x(selected), spread(0.0_dp, 1, sub % num_points()))
    call sub % apply(members, sub % bind(inputs_sub), image)
    call image % real_vector(r_sub)
    scale = max(1.0_dp, maxval(abs(r_whole)))
    call report(sub % num_points() == 3 .and. size(r_sub) == k, 'the constrained residual is over the &
         &selected points', num_failures)
    call report(maxval(abs(r_sub - r_whole(selected))) <= 1.0e-13_dp * scale, &
         & 'the constrained residual with the exterior fixed equals the whole residual on the selected rows', &
         & num_failures)

    call sub % explicit_tangent(members, sub % bind(inputs_sub), 1, rows, columns, weights, tangent_defined)
    allocate(j_sub(k, k), source=0.0_dp)
    do i = 1, size(rows)
       j_sub(rows(i), columns(i)) = j_sub(rows(i), columns(i)) + weights(i)
    end do
    difference = 0.0_dp
    do j = 1, k
       do i = 1, k
          difference = max(difference, abs(j_sub(i, j) - j_whole(selected(i), selected(j))))
       end do
    end do
    scale = max(1.0_dp, maxval(abs(j_whole)))
    call report(tangent_defined .and. difference <= 1.0e-13_dp * scale, &
         & 'the assembled tangent of the constrained residual is the whole tangent''s submatrix', num_failures)

    ! the linearization action, column by column, against the
    ! assembled tangent
    tangent = tangent_of(sub, sub % argument(1))
    call tangent % freeze(inputs_sub)
    call action % state(tangent, members, sub % unknown_domain(), k)
    allocate(basis(k))
    difference = 0.0_dp
    do j = 1, k
       basis    = 0.0_dp
       basis(j) = 1.0_dp
       call action % matvec(basis, y)
       difference = max(difference, maxval(abs(y - j_sub(:, j))))
    end do
    call report(tangent % exact() .and. difference <= 1.0e-13_dp * scale, &
         & 'the exact tangent action on the selection reproduces the assembled tangent', num_failures)

    ! the adjoint: <J u, v> = <u, J^T v> through the transposed stencil
    assembled = stencil(rows, columns, weights, spread(0.0_dp, 1, k), 'selected tangent')
    assembled = assembled % transpose()
    call adjoint % state(assembled, assembled % pattern, assembled % pattern % vertex_set(), k)
    u = [(0.3_dp + 0.1_dp * real(i, dp), i = 1, k)]
    v = [(1.0_dp - 0.2_dp * real(i, dp), i = 1, k)]
    call action  % matvec(u, ju)
    call adjoint % matvec(v, jtv)
    difference = abs(dot_product(ju, v) - dot_product(u, jtv))
    call report(difference <= 1.0e-13_dp * max(abs(dot_product(ju, v)), 1.0_dp), &
         & 'the transposed action on the selection is the adjoint of the tangent action', num_failures)

  end subroutine check_restricted_operator

  !===================================================================!
  ! The temporal partition restricts a copy of its inner template to
  ! every member through the solvers' own restrictions. Three
  ! compositions under three partitions - instants in order, instants
  ! reversed, and the parity members reversed - reach the coupled
  ! solution; the template keeps its whole-domain metadata; and a
  ! failing solver nested below an elimination is reported at the top.
  !===================================================================!

  subroutine check_solver_restriction(num_failures)

    integer, intent(inout) :: num_failures

    integer , parameter :: n = 6
    real(dp), parameter :: h = 0.1_dp, c = 0.5_dp, q0 = 1.0_dp
    character(len=44), parameter :: composition(3) = [character(len=44) :: &
         & 'newton > elimination > gmres > multigrid', &
         & 'newton > gmres > elimination > dense direct', &
         & 'newton > delegating > gmres > multigrid']
    character(len=24), parameter :: layout_named(3) = [character(len=24) :: &
         & 'instants in order', 'instants reversed', 'parity members reversed']
    type(residual_operator) :: residual
    type(stored_directed_graph) :: unknowns
    type(stored_field) :: frozen(2)
    type(newton) :: coupled
    type(temporal_minimizer), allocatable :: solver
    class(minimizer), allocatable :: template
    real(dp), allocatable :: reference(:), q(:), zeros(:), scalar_march(:)
    integer , allocatable :: member_of(:), member_order(:)
    logical , allocatable :: fixed(:)
    real(dp) :: achieved, root
    integer :: kind, layout, p, active_members, m
    type(solve_result) :: outcome

    residual = march_residual(n, h, c, q0)
    unknowns = residual % unknown_graph()
    frozen   = residual % frozen_tuple(spread(0.0_dp, 1, 2 * n), spread(0.0_dp, 1, n))
    allocate(zeros(2 * n), source=0.0_dp)
    fixed = residual % fixed_indicator()

    ! the coupled reference, and the scalar recurrence
    ! q_p + h c q_p^3 = q_{p-1} it must reproduce
    allocate(coupled % inner, source=dense_direct())
    coupled % tolerance = 1.0e-12_dp
    coupled % criterion = absolute
    call coupled % state(residual, unknowns, residual % unknown_domain(), 2 * n, stored_inputs=[frozen(2)])
    call seeded(residual, n, q0, reference)
    call coupled % solve(zeros, reference, achieved)
    outcome = coupled % result()
    allocate(scalar_march(n))
    scalar_march(1) = q0
    do p = 2, n
       root = scalar_march(p - 1)
       do m = 1, 50
          root = root - (root + h * c * root ** 3 - scalar_march(p - 1)) / (1.0_dp + 3.0_dp * h * c * root ** 2)
       end do
       scalar_march(p) = root
    end do
    call report(outcome % converged() .and. maxval(abs(reference(1::2) - scalar_march)) <= 1.0e-10_dp, &
         & 'the coupled newton solve of the march reproduces the scalar implicit recurrence', num_failures)

    do kind = 1, 3
       do layout = 1, 3
          call labels(layout, n, member_of, member_order)
          call composed(kind, n, template)
          if (allocated(solver)) deallocate(solver)
          allocate(solver)
          call move_alloc(template, solver % inner)
          call solver % state(residual, unknowns, residual % unknown_domain(), 2 * n, stored_inputs=[frozen(2)])
          solver % tolerance      = 1.0e-10_dp
          solver % criterion      = absolute
          solver % max_iterations = 40
          call solver % partition(member_of, member_order, seed_from_previous=.true.)
          call seeded(residual, n, q0, q)
          num_restrictions = 0
          call solver % solve(zeros, q, achieved)
          outcome = solver % result()
          call report(outcome % converged() .and. maxval(abs(q - reference)) <= 1.0e-9_dp, &
               & trim(composition(kind)) // ' under ' // trim(layout_named(layout)) // &
               & ' reaches the coupled solution', num_failures)
          call report(template_intact(kind, n, solver % inner), &
               & trim(composition(kind)) // ': the template keeps its whole-domain metadata', num_failures)
          if (kind == 3) then
             active_members = 0
             do m = 1, maxval(member_of)
                if (any(.not. fixed .and. member_of == m)) active_members = active_members + 1
             end do
             call report(num_restrictions == active_members * outcome % iterations, &
                  & 'the test-only solver is restricted once per member per pass', num_failures)
          end if
       end do
    end do

    call labels(1, n, member_of, member_order)
    call composed(4, n, template)
    if (allocated(solver)) deallocate(solver)
    allocate(solver)
    call move_alloc(template, solver % inner)
    call solver % state(residual, unknowns, residual % unknown_domain(), 2 * n, stored_inputs=[frozen(2)])
    solver % tolerance = 1.0e-10_dp
    solver % criterion = absolute
    call solver % partition(member_of, member_order, seed_from_previous=.true.)
    call seeded(residual, n, q0, q)
    call solver % solve(zeros, q, achieved)
    outcome = solver % result()
    call report(outcome % failed() .and. outcome % reason == SOLVE_INNER_FAILED .and. outcome % iterations == 0, &
         & 'a failing solver below an elimination below newton is reported by the temporal minimizer', &
         & num_failures)

    ! the storage limit is metadata of every restricted member: an
    ! elimination within one entry is refused in the first member and
    ! reported at the top without an allocation or a stop
    call composed(5, n, template)
    if (allocated(solver)) deallocate(solver)
    allocate(solver)
    call move_alloc(template, solver % inner)
    call solver % state(residual, unknowns, residual % unknown_domain(), 2 * n, stored_inputs=[frozen(2)])
    solver % tolerance = 1.0e-10_dp
    solver % criterion = absolute
    call solver % partition(member_of, member_order, seed_from_previous=.true.)
    call seeded(residual, n, q0, q)
    call solver % solve(zeros, q, achieved)
    outcome = solver % result()
    call report(outcome % failed() .and. outcome % reason == SOLVE_INNER_FAILED .and. outcome % iterations == 0, &
         & 'an elimination refused within its storage limit below newton is reported by the temporal minimizer', &
         & num_failures)

  end subroutine check_solver_restriction

  !===================================================================!
  ! THE RESIDUAL BOUNDARY. One frozen tuple (Q, nu) on U x P is read
  ! by every consumer: the value, the explicit tangent, the tangent
  ! action, the frozen linearization and its transpose, the design
  ! tangent, the mixed and third partials, and Newton by the explicit
  ! and by the tangent action. On the fixed rows F: R_i = Q_i - h_i,
  ! (D_Q R[v])_i = v_i, J has unit rows, and the design and higher
  ! partials vanish. The pairing: under <u, v>_M = u^T M v the
  ! adjoint of J is M^-1 J^T M; the coordinate transpose violates the
  ! adjoint identity by exactly u^T (J^T M - M J^T) v; the sensitivity
  ! of F = <g, Q>_M in the design by the tangent solve (J w = -D_nu R e)
  ! and by the transposed solve (J^T lambda = M g) agree with a central
  ! difference. The temporal partition restricts a point-varying
  ! design to each member's points.
  !===================================================================!

  subroutine check_residual_boundary(num_failures)

    integer, intent(inout) :: num_failures

    integer , parameter :: n = 5, m = 2 * n
    real(dp), parameter :: h = 0.1_dp, c = 0.5_dp, q0 = 1.0_dp
    real(dp), parameter :: eps = 1.0e-3_dp, delta = 1.0e-5_dp
    type(residual_operator) :: residual, lin, lin_t, moved
    type(stored_directed_graph) :: unknowns
    type(graph) :: u_domain, p_domain
    type(typed_field_domain) :: states, designs
    type(stored_field) :: inputs(2), shifted(2), lin_inputs(2)
    type(newton) :: by_entries, by_action
    type(dense_direct) :: factorisation
    type(temporal_minimizer) :: partitioned
    type(solve_result) :: outcome, outcome_action
    class(field), allocatable :: image
    integer , allocatable :: rows(:), columns(:), member_of(:), member_order(:)
    real(dp), allocatable :: weights(:), x(:), nu(:), r(:), j(:,:), jm(:,:), y(:), basis(:), w(:)
    real(dp), allocatable :: v1(:), v2(:), v3(:), wnu(:), y12(:), y21(:), plus(:), minus(:), y1n(:)
    real(dp), allocatable :: y123(:), rhs(:), mm(:), u(:), vv(:), g(:), lambda(:), b(:), q(:), q2(:)
    real(dp), allocatable :: expected(:), zeros(:)
    real(dp) :: scale, difference, achieved, left, right, violation, exact_violation, df_t, df_a, df_d
    logical :: tangent_defined, on_domain
    integer :: p, i, k

    allocate(x(m), nu(n), zeros(m))
    zeros = 0.0_dp
    do p = 1, n
       x(2 * p - 1) = q0 * 0.9_dp ** p
       x(2 * p)     = -0.3_dp / real(p, dp)
       nu(p)        = 0.2_dp + 0.1_dp * real(p, dp)
    end do
    residual = designed_march_residual(n, h, c, q0, nu(1))
    unknowns = residual % unknown_graph()
    u_domain = residual % unknown_domain()
    p_domain = residual % design_domain()
    inputs   = residual % frozen_tuple(x, nu)
    states   = residual % state_fields()
    designs  = residual % design_fields()
    call report(inputs(1) % defined_on(u_domain) .and. inputs(2) % defined_on(p_domain) .and. &
         & inputs(2) % num_entries() == n, 'the frozen tuple is a state on U and a design on P, one value per point', &
         & num_failures)

    ! the value: physics rows, tying rows, fixed rows
    call residual % apply(unknowns, residual % bind(inputs), image)
    on_domain = image % defined_on(u_domain)
    call image % real_vector(r)
    allocate(expected(m))
    expected(1) = x(1) - q0
    expected(2) = x(2) - (-c * q0 ** 3 - nu(1) * q0)
    do p = 2, n
       expected(2 * p - 1) = x(2 * p) + c * x(2 * p - 1) ** 3 + nu(p) * x(2 * p - 1)
       expected(2 * p)     = x(2 * p) - (x(2 * p - 1) - x(2 * p - 3)) / h
    end do
    scale = max(1.0_dp, maxval(abs(expected)))
    call report(on_domain .and. maxval(abs(r - expected)) <= 1.0e-13_dp * scale, &
         & 'the residual on U: the physics with the point design, the tying rows, and Q_i - h_i on F', num_failures)

    ! the explicit tangent: entries, unit rows on F
    call residual % explicit_tangent(unknowns, residual % bind(inputs), 1, rows, columns, weights, tangent_defined)
    call dense_of_triples(m, rows, columns, weights, j)
    difference = 0.0_dp
    do p = 2, n
       difference = max(difference, abs(j(2 * p - 1, 2 * p - 1) - (3.0_dp * c * x(2 * p - 1) ** 2 + nu(p))), &
            & abs(j(2 * p - 1, 2 * p) - 1.0_dp), abs(j(2 * p, 2 * p) - 1.0_dp), &
            & abs(j(2 * p, 2 * p - 1) + 1.0_dp / h), abs(j(2 * p, 2 * p - 3) - 1.0_dp / h))
    end do
    scale = max(1.0_dp, maxval(abs(j)))
    call report(tangent_defined .and. difference <= 1.0e-13_dp * scale .and. &
         & maxval(abs(j(1, :) - [(merge(1.0_dp, 0.0_dp, k == 1), k = 1, m)])) == 0.0_dp .and. &
         & maxval(abs(j(2, :) - [(merge(1.0_dp, 0.0_dp, k == 2), k = 1, m)])) == 0.0_dp, &
         & 'the explicit tangent: the physics and design partials, the tying weights, unit rows on F', num_failures)

    ! the tangent action column by column against the explicit entries
    allocate(basis(m))
    difference = 0.0_dp
    do k = 1, m
       basis    = 0.0_dp
       basis(k) = 1.0_dp
       call varied(residual, unknowns, inputs, 1, u_domain, basis, y)
       difference = max(difference, maxval(abs(y - j(:, k))))
    end do
    call report(difference <= 1.0e-13_dp * scale, &
         & 'the tangent action D_Q R[e_k] is the k-th column of the explicit tangent, F rows included', num_failures)

    ! the frozen linearization and its transpose on the same U and P
    rhs = [(0.1_dp * real(k, dp) - 0.4_dp, k = 1, m)]
    w   = [(0.5_dp - 0.07_dp * real(k, dp), k = 1, m)]
    lin   = residual % linearize(unknowns, residual % bind(inputs), rhs, .false., 7)
    lin_t = residual % linearize(unknowns, residual % bind(inputs), rhs, .true., 8)
    lin_inputs = lin % frozen_tuple(w, nu)
    call lin % apply(lin % unknown_graph(), lin % bind(lin_inputs), image)
    call image % real_vector(y)
    difference = maxval(abs(y - (matmul(j, w) - rhs)))
    call lin_t % apply(lin_t % unknown_graph(), lin_t % bind(lin_inputs), image)
    call image % real_vector(y)
    difference = max(difference, maxval(abs(y - (matmul(transpose(j), w) - rhs))))
    call report(difference <= 1.0e-13_dp * scale .and. lin % version() == 7 .and. lin_t % version() == 8 .and. &
         & lin_t % transpose_version() .and. .not. lin % transpose_version(), &
         & 'the frozen linearization computes J w - rhs and its transpose J^T w - rhs, versioned', num_failures)
    call report(u_domain % same_as(lin % unknown_domain()) .and. p_domain % same_as(lin % design_domain()) .and. &
         & u_domain % same_as(lin_t % unknown_domain()), &
         & 'A = D_Q R maps U to U: the frozen linearization keeps the unknown and point domains', num_failures)

    ! the design tangent: the physics partial at each point, zero on F
    allocate(wnu(n))
    wnu = [(1.0_dp + 0.3_dp * real(p, dp), p = 1, n)]
    call varied(residual, unknowns, inputs, 2, p_domain, wnu, y)
    expected = 0.0_dp
    do p = 2, n
       expected(2 * p - 1) = wnu(p) * x(2 * p - 1)
    end do
    call report(maxval(abs(y - expected)) <= 1.0e-13_dp * scale, &
         & 'the design tangent D_nu R[w] is w_p q_p on the governed rows and zero on F and the tying rows', &
         & num_failures)

    ! second partials: symmetric, exact, and against a central
    ! difference of the first partial along the second direction
    v1 = [(0.3_dp + 0.1_dp * real(k, dp), k = 1, m)]
    v2 = [(1.0_dp - 0.15_dp * real(k, dp), k = 1, m)]
    v3 = [(0.2_dp * real(k, dp) - 0.9_dp, k = 1, m)]
    call varied(residual, unknowns, inputs, 1, u_domain, v1, y12, 1, u_domain, v2)
    call varied(residual, unknowns, inputs, 1, u_domain, v2, y21, 1, u_domain, v1)
    expected = 0.0_dp
    do p = 2, n
       expected(2 * p - 1) = 6.0_dp * c * x(2 * p - 1) * v1(2 * p - 1) * v2(2 * p - 1)
    end do
    shifted = residual % frozen_tuple(x + eps * v2, nu)
    call varied(residual, unknowns, shifted, 1, u_domain, v1, plus)
    shifted = residual % frozen_tuple(x - eps * v2, nu)
    call varied(residual, unknowns, shifted, 1, u_domain, v1, minus)
    call report(maxval(abs(y12 - y21)) <= 1.0e-13_dp * scale .and. maxval(abs(y12 - expected)) <= 1.0e-13_dp * scale &
         & .and. maxval(abs(y12 - (plus - minus) / (2.0_dp * eps))) <= 1.0e-9_dp * scale, &
         & 'D^2 R[v1, v2] is symmetric, equals 6 c q v1 v2 on the governed rows, zero on F, and is the &
         &central difference of D_Q R[v1] along v2', num_failures)

    ! the mixed state-design partial against a difference of the
    ! design tangent along the state direction
    call varied(residual, unknowns, inputs, 1, u_domain, v1, y1n, 2, p_domain, wnu)
    shifted = residual % frozen_tuple(x + eps * v1, nu)
    call varied(residual, unknowns, shifted, 2, p_domain, wnu, plus)
    shifted = residual % frozen_tuple(x - eps * v1, nu)
    call varied(residual, unknowns, shifted, 2, p_domain, wnu, minus)
    expected = 0.0_dp
    do p = 2, n
       expected(2 * p - 1) = wnu(p) * v1(2 * p - 1)
    end do
    call report(maxval(abs(y1n - expected)) <= 1.0e-13_dp * scale .and. &
         & maxval(abs(y1n - (plus - minus) / (2.0_dp * eps))) <= 1.0e-9_dp * scale, &
         & 'the mixed partial D_Q D_nu R[v, w] is w_p v_p on the governed rows and the difference of D_nu R[w] &
         &along v', num_failures)

    ! the third partial in the state, and the third mixed partial
    call residual % partial_action(unknowns, residual % bind(inputs), &
         & [variation(residual % argument(1), states % direction(v1)), &
         &  variation(residual % argument(1), states % direction(v2)), &
         &  variation(residual % argument(1), states % direction(v3))], image)
    call image % real_vector(y123)
    expected = 0.0_dp
    do p = 2, n
       expected(2 * p - 1) = 6.0_dp * c * v1(2 * p - 1) * v2(2 * p - 1) * v3(2 * p - 1)
    end do
    shifted = residual % frozen_tuple(x + eps * v3, nu)
    call varied(residual, unknowns, shifted, 1, u_domain, v1, plus, 1, u_domain, v2)
    shifted = residual % frozen_tuple(x - eps * v3, nu)
    call varied(residual, unknowns, shifted, 1, u_domain, v1, minus, 1, u_domain, v2)
    difference = maxval(abs(y123 - (plus - minus) / (2.0_dp * eps)))
    call residual % partial_action(unknowns, residual % bind(inputs), &
         & [variation(residual % argument(1), states % direction(v1)), &
         &  variation(residual % argument(1), states % direction(v2)), &
         &  variation(residual % argument(2), designs % direction(wnu))], image)
    call image % real_vector(y)
    call report(maxval(abs(y123 - expected)) <= 1.0e-13_dp * scale .and. difference <= 1.0e-9_dp * scale .and. &
         & maxval(abs(y)) == 0.0_dp, &
         & 'D^3 R[v1, v2, v3] is 6 c v1 v2 v3 on the governed rows, the difference of D^2 R[v1, v2] along v3, and &
         &D^2_Q D_nu R vanishes', num_failures)

    ! Newton by the explicit entries and by the tangent action from
    ! one seed reach one iterate
    allocate(by_entries % inner, source=dense_direct())
    allocate(by_action  % inner, source=dense_direct())
    by_entries % explicit = .true.
    by_action  % explicit = .false.
    by_entries % tolerance = 1.0e-12_dp
    by_action  % tolerance = 1.0e-12_dp
    by_entries % criterion = absolute
    by_action  % criterion = absolute
    call by_entries % state(residual, unknowns, u_domain, m, stored_inputs=[inputs(2)])
    call by_action  % state(residual, unknowns, u_domain, m, stored_inputs=[inputs(2)])
    call seeded(residual, n, q0, q)
    call by_entries % solve(zeros, q, achieved)
    outcome = by_entries % result()
    call seeded(residual, n, q0, q2)
    call by_action % solve(zeros, q2, achieved)
    outcome_action = by_action % result()
    call report(outcome % converged() .and. outcome_action % converged() .and. &
         & maxval(abs(q - q2)) <= 1.0e-10_dp, &
         & 'newton by the explicit entries and by the tangent action reach one iterate from one seed', num_failures)

    ! the pairing <u, v>_M = u^T M v: the M-adjoint of J is M^-1 J^T M;
    ! the coordinate transpose violates the identity by u^T (J^T M - M J^T) v
    mm = [(1.0_dp + 0.4_dp * real(k, dp), k = 1, m)]
    u  = [(0.3_dp + 0.1_dp * real(k, dp), k = 1, m)]
    vv = [(1.0_dp - 0.2_dp * real(k, dp), k = 1, m)]
    allocate(jm(m, m))
    do i = 1, m
       do k = 1, m
          jm(i, k) = j(k, i) * mm(k) / mm(i)
       end do
    end do
    left  = dot_product(matmul(j, u), mm * vv)
    right = dot_product(u, mm * matmul(jm, vv))
    violation       = left - dot_product(u, mm * matmul(transpose(j), vv))
    exact_violation = 0.0_dp
    do i = 1, m
       do k = 1, m
          exact_violation = exact_violation + u(i) * (j(k, i) * mm(k) - mm(i) * j(k, i)) * vv(k)
       end do
    end do
    write(*, '(a, es12.4, a, es12.4)') '        <J u, v>_M = ', left, &
         & '   coordinate-transpose violation of the M-adjoint identity = ', violation
    call report(abs(left - right) <= 1.0e-13_dp * max(abs(left), 1.0_dp) .and. &
         & abs(violation - exact_violation) <= 1.0e-13_dp * max(abs(left), 1.0_dp) .and. &
         & abs(violation) > 1.0e-13_dp * max(abs(left), 1.0_dp), &
         & 'under <u, v>_M the adjoint of J is M^-1 J^T M, and the coordinate transpose violates the identity by &
         &u^T (J^T M - M J^T) v', num_failures)

    ! the pairing law on the typed supports: a direction u on U, its
    ! image J u on Y and a costate v on Y pair as <J u, v>_Y = <u, J^T v>_U
    ! under the Euclidean measure, evaluated through the fields'
    ! inner product; a costate pairs with the residual on Y
    block
      type(typed_field_domain) :: residuals
      type(stored_field) :: u_field, v_field, ju_field, jtv_field, lambda_field
      real(dp) :: paired_left, paired_right, paired_lagrangian
      residuals  = residual % residual_fields()
      u_field    = states % direction(u)
      call residual % partial_action(unknowns, residual % bind(inputs), &
           & [variation(residual % argument(1), u_field)], image)
      call image % real_vector(y)
      ju_field   = residuals % residual(y, 'J u')
      v_field    = residuals % costate(vv)
      jtv_field  = states % direction(matmul(transpose(j), vv))
      lambda_field = residuals % costate(vv)
      paired_left  = ju_field % inner_product(v_field)
      paired_right = u_field % inner_product(jtv_field)
      paired_lagrangian = lambda_field % inner_product(image)
      call report(abs(paired_left - dot_product(matmul(j, u), vv)) <= 1.0e-13_dp * max(abs(paired_left), 1.0_dp) &
           & .and. abs(paired_left - paired_right) <= 1.0e-13_dp * max(abs(paired_left), 1.0_dp) &
           & .and. abs(paired_lagrangian - dot_product(vv, y)) <= 1.0e-13_dp * max(abs(paired_lagrangian), 1.0_dp) &
           & .and. ju_field % defined_on(u_domain) .and. u_field % defined_on(u_domain), &
           & 'through the typed supports <J u, v>_Y = <u, J^T v>_U under the Euclidean measure, and the costate &
           &pairs with the residual on Y = U', num_failures)
    end block

    ! the sensitivity of F = <g, Q(nu)>_M along the design direction
    ! wnu: tangent J w = -D_nu R[wnu], adjoint J^T lambda = M g, and a
    ! central difference of the solved F
    g = [(0.1_dp * real(k, dp), k = 1, m)]
    call seeded(residual, n, q0, q)
    call by_entries % solve(zeros, q, achieved)
    inputs = residual % frozen_tuple(q, nu)
    call varied(residual, unknowns, inputs, 2, p_domain, wnu, b)
    lin   = residual % linearize(unknowns, residual % bind(inputs), zeros, .false., 9)
    lin_t = residual % linearize(unknowns, residual % bind(inputs), zeros, .true., 10)
    call factorisation % state(lin, lin % unknown_graph(), lin % unknown_domain(), m, stored_inputs=[inputs(2)])
    allocate(lambda(m))
    w = 0.0_dp
    call factorisation % solve(-b, w, achieved)
    df_t = dot_product(g, mm * w)
    call factorisation % state(lin_t, lin_t % unknown_graph(), lin_t % unknown_domain(), m, stored_inputs=[inputs(2)])
    lambda = 0.0_dp
    call factorisation % solve(mm * g, lambda, achieved)
    df_a = -dot_product(lambda, b)
    ! the moved design is a residual of its own, on its own U and P;
    ! the fixed values h are data of the residual and do not move
    moved   = designed_march_residual(n, h, c, q0, nu(1))
    shifted = moved % frozen_tuple(zeros, nu + delta * wnu)
    call by_entries % state(moved, moved % unknown_graph(), moved % unknown_domain(), m, stored_inputs=[shifted(2)])
    call seeded(moved, n, q0, plus)
    call by_entries % solve(zeros, plus, achieved)
    moved   = designed_march_residual(n, h, c, q0, nu(1))
    shifted = moved % frozen_tuple(zeros, nu - delta * wnu)
    call by_entries % state(moved, moved % unknown_graph(), moved % unknown_domain(), m, stored_inputs=[shifted(2)])
    call seeded(moved, n, q0, minus)
    call by_entries % solve(zeros, minus, achieved)
    df_d = (dot_product(g, mm * plus) - dot_product(g, mm * minus)) / (2.0_dp * delta)
    write(*, '(a, 3es20.12)') '        dF/dnu[w]: tangent, adjoint, difference ', df_t, df_a, df_d
    call report(abs(df_t - df_a) <= 1.0e-11_dp * max(abs(df_t), 1.0_dp) .and. &
         & abs(df_t - df_d) <= 1.0e-7_dp * max(abs(df_t), 1.0_dp), &
         & 'the M-weighted sensitivity by the tangent and by the transposed solve agree, and with a central &
         &difference', num_failures)

    ! the temporal partition with the point-varying design: each
    ! member reads the design at its own points
    call labels(1, n, member_of, member_order)
    allocate(partitioned % inner, source=by_entries)
    call partitioned % state(residual, unknowns, u_domain, m, stored_inputs=[inputs(2)])
    partitioned % tolerance      = 1.0e-10_dp
    partitioned % criterion      = absolute
    partitioned % max_iterations = 40
    call partitioned % partition(member_of, member_order, seed_from_previous=.true.)
    call seeded(residual, n, q0, q2)
    call partitioned % solve(zeros, q2, achieved)
    outcome = partitioned % result()
    call report(outcome % converged() .and. maxval(abs(q2 - q)) <= 1.0e-9_dp, &
         & 'the temporal partition restricts the point-varying design to each member and reaches the coupled &
         &solution', num_failures)

  end subroutine check_residual_boundary
  !===================================================================!
  ! ONE LAW, TWO DISCRETIZATIONS. A continuous law is placed on an
  ! instant graph of five vertices and on a cell graph of nine: both
  ! placements evaluate the same law to the same value at equal
  ! tuples, their state supports are distinct identities of extents
  ! 5 and 9 with the law's component count, and a residual placed on
  ! its own point graph accepts the design typed by the placement on
  ! that graph (the other placement is refused, restriction_refusal
  ! case other_placement).
  !===================================================================!

  subroutine check_one_law_two_placements(num_failures)

    integer, intent(inout) :: num_failures

    integer , parameter :: n = 5, m = 2 * n
    real(dp), parameter :: h = 0.1_dp, c = 0.5_dp, q0 = 1.0_dp
    type(expression) :: law_expression, law_a, law_b
    type(continuous_domain) :: law
    type(discrete_domain) :: on_instants, on_cells, on_points
    type(stored_directed_graph) :: instants, cells, unknowns
    type(graph) :: instant_identity, cell_identity
    type(continuous_domain) :: residual_law
    type(typed_field_domain) :: instant_states, cell_states, instant_designs, cell_designs, point_designs
    type(stored_field) :: state_a, design_a, state_b, design_b, inputs(2), nu_field
    type(residual_operator) :: residual
    class(field), allocatable :: out, image
    real(dp), allocatable :: ra(:), rb(:), r(:), r_placed(:), x(:), nu(:)
    real(dp) :: tuple(2), value
    integer :: k, p

    law_expression = stated(derivative(unknown(), 1) + constant(c) * derivative(unknown(), 0) ** 3, 1, 'cubic decay')
    law      = continuous_domain(law_expression)
    instants = stored_directed_graph(5, tails=[integer ::], heads=[integer ::])
    cells    = stored_directed_graph(9, tails=[integer ::], heads=[integer ::])
    on_instants = law % discrete(instants)
    on_cells    = law % discrete(cells)
    instant_states  = on_instants % state_fields()
    cell_states     = on_cells % state_fields()
    instant_designs = on_instants % design_fields()
    cell_designs    = on_cells % design_fields()
    tuple = [0.8_dp, -0.3_dp]
    value = tuple(2) + c * tuple(1) ** 3
    state_a  = instant_states % state([(tuple, k = 1, 5)])
    design_a = instant_designs % design([(0.2_dp, k = 1, 5)])
    state_b  = cell_states % state([(tuple, k = 1, 9)])
    design_b = cell_designs % design([(0.2_dp, k = 1, 9)])
    law_a = on_instants % law()
    law_b = on_cells % law()
    call law_a % apply(instants, law_a % bind([state_a, design_a]), out)
    call out % real_vector(ra)
    call law_b % apply(cells, law_b % bind([state_b, design_b]), out)
    call out % real_vector(rb)
    call report(on_instants % num_points() == 5 .and. on_cells % num_points() == 9 .and. &
         & on_instants % equation_degree() == on_cells % equation_degree() .and. &
         & on_instants % num_components() == on_cells % num_components() .and. &
         & size(ra) == 5 .and. size(rb) == 9 .and. &
         & maxval(abs(ra - value)) <= 1.0e-14_dp .and. maxval(abs(rb - value)) <= 1.0e-14_dp .and. &
         & out % defined_on(cells % vertex_set()), &
         & 'one continuous law placed on five instants and on nine cells evaluates to dq/dt + c q^3 at every &
         &point of either placement', num_failures)
    instant_identity = instant_states % domain()
    cell_identity    = cell_states % domain()
    call report(.not. instant_identity % same_as(cell_identity) .and. &
         & instant_identity % same_as(instants % vertex_set()) .and. &
         & cell_identity % same_as(cells % vertex_set()) .and. &
         & instant_states % num_entries() == 5 .and. cell_states % num_entries() == 9 .and. &
         & instant_states % num_components() == on_instants % num_components() .and. &
         & cell_states % num_components() == on_cells % num_components() .and. &
         & instant_designs % num_components() == 1 .and. cell_designs % num_entries() == 9, &
         & 'the two placements keep their own discrete-domain identities: extents 5 and 9, the law''s component &
         &count on the state, one value per point on the design', num_failures)

    ! the residual placed on its own points accepts the design typed
    ! by the placement of its law on that same point graph
    allocate(x(m), nu(n))
    do p = 1, n
       x(2 * p - 1) = q0 * 0.9_dp ** p
       x(2 * p)     = -0.3_dp / real(p, dp)
       nu(p)        = 0.2_dp + 0.1_dp * real(p, dp)
    end do
    residual  = designed_march_residual(n, h, c, q0, nu(1))
    unknowns  = residual % unknown_graph()
    inputs    = residual % frozen_tuple(x, nu)
    residual_law = continuous_domain(residual % rule())
    on_points    = residual_law % discrete(residual % point_domain())
    point_designs = on_points % design_fields()
    nu_field  = point_designs % design(nu)
    call residual % apply(unknowns, residual % bind(inputs), image)
    call image % real_vector(r)
    call residual % apply(unknowns, residual % bind([inputs(1), nu_field]), image)
    call image % real_vector(r_placed)
    call report(nu_field % defined_on(residual % design_domain()) .and. on_points % num_points() == n .and. &
         & maxval(abs(r - r_placed)) == 0.0_dp, &
         & 'the law placed on the residual''s own point graph types a design the residual accepts, with the &
         &same value', num_failures)

  end subroutine check_one_law_two_placements


end program test_graph_minimization
