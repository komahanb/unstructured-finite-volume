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
       & SOLVE_INNER_FAILED, SOLVE_NOT_STARTED, SOLVE_STAGNATED
  use operation_linearization, only : linearization, tangent_of
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf, ieee_is_finite, ieee_is_nan
  use operation_dense_direct, only : dense_direct
  use operation_stencil, only : stencil
  use operation_balance  , only : balance
  use cubic_statement_fixture, only : cubic_statement

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

end program test_graph_minimization
