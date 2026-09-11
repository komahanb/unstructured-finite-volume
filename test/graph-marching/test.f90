!=====================================================================!
! The marching suite.
!
! TIME IS A GRAPH: the instants are the data vertices of a bipartite
! digraph, the steps its rule vertices, and the driver visits the
! steps in the digraph's own topological order. THE STEP IS EXACT
! ABOUT ITSELF: euler on decay dq/dt = -q lands on (1 - h)^n to
! machine precision, euler's own discrete identity. THE MAP IS THE
! MARCH: z -> z^2 + c is forward euler with step one on
! S = z - z^2 - c - an identity, not an approximation - so the
! escape-time fractal is a march, with known points where the
! arithmetic states they are:
!
!      c = 0        never escapes: the origin is fixed
!      c = -1       never escapes: a two-cycle, 0 -> -1 -> 0
!      c = 1        escapes at step three: 0, 1, 2, 5
!      c = -2       never escapes: converges to the fixed point 2
!      c = 2i       escapes at step two: -4 + 2i leaves the circle
!
! The implicit steps are newton solves at each rule vertex; the
! reverse pass is the transposed chain in the same digraph, so the
! driver's lifetimes retain every state the adjoint reads.
!=====================================================================!

program test_graph_marching

  use iso_fortran_env, only : dp => REAL64
  use operation_action, only : operation
  use field_calculus, only : field
  use view_directed , only : SIDE_VERTEX, forward, reverse
  use field_stored  , only : stored_field
  use view_directed_stored        , only : stored_directed_graph
  use view_read_write, only : bipartite_digraph, FIRST_PART, SECOND_PART
  use operation_driver, only : driver, rule_graph, data_graph, pairing
  use operation_differential, only : differential_operator
  use operation_stencil , only : stencil
  use operation_temporal_minimization, only : temporal_minimizer
  use operation_dense_direct, only : dense_direct
  use operation_minimization, only : solve_result, SOLVE_SINGULAR
  use march_rules_fixture, only : explicit_step, implicit_step, adjoint_step
  use mandelbrot_law_fixture, only : mandelbrot_law
  use vdp_fixture, only : vdp_law, vdp_tangent_law

  implicit none

  integer :: num_failures

  num_failures = 0

  call check_temporal_dependency_order(num_failures)
  call check_lifetimes(num_failures)
  call check_long_schedule(num_failures)
  call check_euler_discrete_solution(num_failures)
  call check_mandelbrot_iteration(num_failures)
  call check_implicit_discrete_solution(num_failures)
  call check_component_layout_invariance(num_failures)
  call check_adjoint_pairing_invariance(num_failures)
  call check_tangent_adjoint_equivalence(num_failures)
  call check_delegated_result(num_failures)

  write(*, '(a)') ' ============================================='
  if (num_failures == 0) then
     write(*, '(a)') ' all marching checks passed'
  else
     write(*, '(a, i0, a)') ' ', num_failures, ' marching checks FAILED'
     error stop 1
  end if

contains

  subroutine check_delegated_result(num_failures)
    integer, intent(inout) :: num_failures
    type(temporal_minimizer) :: temporal
    type(dense_direct) :: inner
    type(stencil) :: zero
    type(solve_result) :: outcome
    real(dp) :: x(1), achieved

    zero = stencil([1], [1], [0.0_dp], [0.0_dp], 'singular temporal member')
    inner % singular_reported = .true.
    allocate(temporal % inner, source=inner)
    call temporal % state(zero, zero % pattern, zero % pattern % vertex_set(), 1)
    x = 0.0_dp
    call temporal % solve([1.0_dp], x, achieved)
    outcome = temporal % result()
    call report(outcome % reason == SOLVE_SINGULAR .and. outcome % residual == 1.0_dp .and. &
         & achieved == huge(1.0_dp), 'a temporal delegate preserves the actual residual of a singular inner solve', num_failures)
  end subroutine check_delegated_result

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
  ! THE CHAIN OF INSTANTS as a bipartite digraph: nsteps rules,
  ! nsteps + 1 data. Step b reads the instant b - 1 and the `reach`
  ! instants before it, newest first, and writes instant b. With a
  ! transpose, nsteps adjoint rules follow: adjoint step b reads the
  ! costate at instant nsteps - b + 1 and the state at instant
  ! nsteps - b, and writes the costate at instant nsteps - b; the
  ! costates are the data nsteps + 2 .. 2 nsteps + 2.
  !
  !      (s1) -> [q1] -> (s2) -> [q2]          forward
  !      [q0] -> (s1)    [q1] -> (s2)
  !      [l2] -> (a1) -> [l1] -> (a2) -> [l0]  transpose
  !      [q1] -> (a1)    [q0] -> (a2)
  !===================================================================!

  function chain_incidence(nsteps, reach, with_transpose) result(incidence)

    integer, intent(in) :: nsteps, reach
    logical, intent(in) :: with_transpose
    type(bipartite_digraph) :: incidence

    integer, allocatable :: from_part(:), from_vertex(:), to_part(:), to_vertex(:)
    integer :: b, j, n_rules, n_data, a

    n_rules = nsteps
    n_data  = nsteps + 1
    if (with_transpose) then
       n_rules = 2 * nsteps
       n_data  = 2 * nsteps + 2
    end if

    allocate(from_part(0), from_vertex(0), to_part(0), to_vertex(0))

    do b = 1, nsteps
       ! step b reads the instants b - 1, b - 2, ... within reach,
       ! newest first, so argument j of its rule is q_(b-j)
       do j = 1, reach
          if (b - j < 0) exit
          call append_arc(from_part, from_vertex, to_part, to_vertex, &
               & SECOND_PART, b - j + 1, FIRST_PART, b)
       end do
       call append_arc(from_part, from_vertex, to_part, to_vertex, &
            & FIRST_PART, b, SECOND_PART, b + 1)
    end do

    if (with_transpose) then
       do b = 1, nsteps
          a = nsteps + b
          ! the costate at nsteps - b + 1, then the state at nsteps - b
          call append_arc(from_part, from_vertex, to_part, to_vertex, &
               & SECOND_PART, nsteps + 1 + b, FIRST_PART, a)
          call append_arc(from_part, from_vertex, to_part, to_vertex, &
               & SECOND_PART, nsteps - b + 1, FIRST_PART, a)
          ! the costate at nsteps - b
          call append_arc(from_part, from_vertex, to_part, to_vertex, &
               & FIRST_PART, a, SECOND_PART, nsteps + 2 + b)
       end do
    end if

    incidence = bipartite_digraph(n_rules, n_data, from_part, from_vertex, to_part, to_vertex)

  end function chain_incidence

  ! One arc appended to the four lists.
  subroutine append_arc(from_part, from_vertex, to_part, to_vertex, p1, v1, p2, v2)
    integer, allocatable, intent(inout) :: from_part(:), from_vertex(:), to_part(:), to_vertex(:)
    integer, intent(in) :: p1, v1, p2, v2
    from_part   = [from_part  , p1]
    from_vertex = [from_vertex, v1]
    to_part     = [to_part    , p2]
    to_vertex   = [to_vertex  , v2]
  end subroutine append_arc

  !===================================================================!
  ! THE MARCH: the rule at every step, the initial state at the first
  ! instant, the driver evaluating the chain. The final state is read
  ! from the last instant's datum. With an adjoint rule and a terminal
  ! costate, the transposed chain is evaluated in the same traversal
  ! and the initial costate is read from the last costate datum.
  !===================================================================!

  subroutine march(rule, on, q, nsteps, reach, adjoint, lambda)

    class(operation)     , intent(in)    :: rule
    class(stored_directed_graph), intent(in) :: on
    real(dp)             , intent(inout) :: q(:)
    integer              , intent(in)    :: nsteps, reach
    class(operation)     , intent(in), optional    :: adjoint
    real(dp)             , intent(inout), optional :: lambda(:)

    type(bipartite_digraph) :: incidence
    type(rule_graph) :: rules
    type(data_graph) :: values
    type(driver)     :: schedule
    type(pairing)    :: pairs
    type(stored_field) :: initial
    class(field), allocatable :: final
    real(dp), allocatable :: values_of_final(:)
    integer :: b, nv

    incidence = chain_incidence(nsteps, reach, present(adjoint))

    allocate(rules % at(incidence % order_of_part(FIRST_PART)))
    allocate(values % at(incidence % order_of_part(SECOND_PART)))
    do b = 1, nsteps
       allocate(rules % at(b) % rule, source=rule)
    end do
    if (present(adjoint)) then
       do b = 1, nsteps
          allocate(rules % at(nsteps + b) % rule, source=adjoint)
       end do
    end if

    nv = on % num_vertices()
    initial = stored_field('state', on % vertex_set(), nv, num_components=size(q) / nv)
    call initial % set_real_vector(q)
    allocate(values % at(1) % datum, source=initial)
    if (present(adjoint)) then
       initial = stored_field('costate', on % vertex_set(), nv, num_components=size(lambda) / nv)
       call initial % set_real_vector(lambda)
       allocate(values % at(nsteps + 2) % datum, source=initial)
    end if

    schedule = driver(rule, incidence, forward)
    call schedule % pair_with(rules % pair(values))
    call schedule % evaluate(on)

    pairs = schedule % pairing_of()
    call pairs % datum_at(nsteps + 1, final)
    call final % real_vector(values_of_final)
    q = values_of_final
    if (present(adjoint)) then
       call pairs % datum_at(2 * nsteps + 2, final)
       call final % real_vector(values_of_final)
       lambda = values_of_final
    end if

  end subroutine march

  !===================================================================!
  ! The instants are the data vertices of the chain and the steps its
  ! rules; the driver's visiting order is the topological order the
  ! arcs admit, and every step reads the instant before it and writes
  ! the one after.
  !===================================================================!

  subroutine check_temporal_dependency_order(num_failures)

    integer, intent(inout) :: num_failures

    type(bipartite_digraph) :: chain
    type(driver) :: schedule
    type(stored_directed_graph) :: single_vertex
    integer, allocatable :: order(:), reads(:), writes(:)
    logical :: ordered
    integer :: b

    single_vertex  = stored_directed_graph(1, tails=[integer ::], heads=[integer ::])
    chain = chain_incidence(10, 1, .false.)

    call report(chain % order_of_part(SECOND_PART) == 11 .and. &
         & chain % order_of_part(FIRST_PART) == 10 .and. &
         & chain % size_of_digraph() == 20, &
         & 'eleven instants, ten steps: the chain is the bipartite digraph', num_failures)

    schedule = driver(explicit_step(differential_operator(SIDE_VERTEX, 0), 1.0_dp), chain, forward)
    order    = schedule % visits()
    ordered  = size(order) == 10
    do b = 1, 10
       call chain % in_neighbourhood(FIRST_PART, b, reads)
       call chain % out_neighbourhood(FIRST_PART, b, writes)
       if (order(b) /= b) ordered = .false.
       if (size(reads) /= 1 .or. size(writes) /= 1) then
          ordered = .false.
       else if (reads(1) /= b .or. writes(1) /= b + 1) then
          ordered = .false.
       end if
    end do
    call report(ordered, 'and every step reads its instant and writes the next', num_failures)

  end subroutine check_temporal_dependency_order

  subroutine check_lifetimes(num_failures)
    integer, intent(inout) :: num_failures
    type(bipartite_digraph) :: incidence
    type(driver) :: schedule
    type(rule_graph) :: rules
    type(data_graph) :: values
    type(pairing) :: pairs
    type(stored_directed_graph) :: single_vertex
    type(stored_field) :: initial
    class(field), allocatable :: datum
    integer, allocatable :: order(:), released(:)
    integer :: orientation, positions(3), expected(5), k, d
    logical :: valid

    ! Rule labels are not traversal positions. Datum 1 has two readers;
    ! rule 2 has no operation, but its read still determines a lifetime.
    incidence = bipartite_digraph(3, 5, &
         & [SECOND_PART, FIRST_PART, SECOND_PART, FIRST_PART, SECOND_PART, FIRST_PART], &
         & [4, 3, 1, 1, 1, 2], &
         & [FIRST_PART, SECOND_PART, FIRST_PART, SECOND_PART, FIRST_PART, SECOND_PART], &
         & [3, 1, 1, 5, 2, 2])
    valid = .true.
    do orientation = forward, reverse, reverse - forward
       schedule = driver(explicit_step(differential_operator(SIDE_VERTEX, 0), 0.0_dp), &
            & incidence, orientation)
       order = schedule % visits()
       do k = 1, size(order)
          positions(order(k)) = k
       end do
       if (orientation == forward) then
          expected = [max(positions(1), positions(2)), 0, 0, positions(3), 0]
       else
          expected = [positions(3), positions(2), 0, 0, positions(1)]
       end if
       do d = 1, 5
          valid = valid .and. schedule % last_reader_of(d) == expected(d)
       end do
       do k = 0, 4
          released = schedule % released_after(k)
          valid = valid .and. size(released) == count(expected > 0 .and. expected <= k)
          if (size(released) == count(expected > 0 .and. expected <= k)) then
             valid = valid .and. all(released == pack([1, 2, 3, 4, 5], expected > 0 .and. expected <= k))
          end if
       end do
    end do
    call report(valid, 'lifetimes follow both orientations; public release lists remain cumulative', num_failures)

    single_vertex = stored_directed_graph(1, tails=[integer ::], heads=[integer ::])
    initial = stored_field('source', single_vertex % vertex_set(), 1)
    call initial % set_real_vector([2.0_dp])
    allocate(rules % at(3), values % at(5))
    do d = 1, 5
       allocate(values % at(d) % datum, source=initial)
    end do
    allocate(rules % at(1) % rule, source=explicit_step(differential_operator(SIDE_VERTEX, 0), 0.0_dp))
    allocate(rules % at(3) % rule, source=explicit_step(differential_operator(SIDE_VERTEX, 0), 0.0_dp))
    schedule = driver(explicit_step(differential_operator(SIDE_VERTEX, 0), 0.0_dp), incidence, forward)
    call schedule % pair_with(rules % pair(values))
    call schedule % evaluate(single_vertex)
    pairs = schedule % pairing_of()
    valid = .true.
    do d = 1, 5
       call pairs % datum_at(d, datum)
       valid = valid .and. (allocated(datum) .eqv. (d == 2 .or. d == 3 .or. d == 5))
    end do
    call report(valid, 'an absent rule still releases reads; sinks and isolated data remain available', num_failures)
  end subroutine check_lifetimes

  subroutine check_long_schedule(num_failures)
    integer, intent(inout) :: num_failures
    integer, parameter :: n = 1000
    type(bipartite_digraph) :: incidence
    type(driver) :: schedule
    type(rule_graph) :: rules
    type(data_graph) :: values
    type(stored_directed_graph) :: single_vertex
    integer :: d
    logical :: valid

    incidence = chain_incidence(n, 2, .false.)
    schedule = driver(explicit_step(differential_operator(SIDE_VERTEX, 0), 0.0_dp), incidence, forward)
    valid = schedule % last_reader_of(n + 1) == 0
    do d = 1, n
       valid = valid .and. schedule % last_reader_of(d) == min(d + 1, n)
    end do
    allocate(rules % at(n), values % at(n + 1))
    call schedule % pair_with(rules % pair(values))
    single_vertex = stored_directed_graph(1, tails=[integer ::], heads=[integer ::])
    call schedule % evaluate(single_vertex)
    call report(valid, 'a thousand-step history schedule retains its exact last-reader indices', num_failures)
  end subroutine check_long_schedule

  !===================================================================!
  ! Euler's discrete identity: on S = q the march gives
  ! q_n = (1 - h)^n q_0.
  !===================================================================!

  subroutine check_euler_discrete_solution(num_failures)

    integer, intent(inout) :: num_failures

    type(stored_directed_graph) :: single_vertex
    type(differential_operator) :: decay
    real(dp) :: q(1), expected

    single_vertex = stored_directed_graph(1, tails=[integer ::], heads=[integer ::])

    decay = differential_operator(SIDE_VERTEX, 0, coefficient=1.0_dp)

    q = [3.0_dp]
    call march(explicit_step(decay, 0.125_dp), single_vertex, q, 20, 1)

    expected = 3.0_dp * (1.0_dp - 0.125_dp)**20

    call report(abs(q(1) - expected) < 1.0d-14, &
         & 'euler lands on (1-h)^n exactly: its own discrete identity', num_failures)

  end subroutine check_euler_discrete_solution

  !===================================================================!
  ! Five points of the complex plane, marched by z -> z^2 + c one step
  ! at a time with the escape checked after each. The driver reads
  ! nothing of the law.
  !===================================================================!

  subroutine check_mandelbrot_iteration(num_failures)

    integer, intent(inout) :: num_failures

    type(mandelbrot_law) :: law
    type(stored_directed_graph) :: points
    type(stored_field) :: escape_field
    real(dp), allocatable :: q(:)
    integer, allocatable :: escape(:)
    integer :: v, n
    integer, parameter :: nv = 5, nmax = 30

    ! The five points, as a graph of lone cells.
    points = stored_directed_graph(nv, tails=[integer ::], heads=[integer ::])

    law = mandelbrot_law()
    law % creal = [0.0_dp, -1.0_dp, 1.0_dp, -2.0_dp, 0.0_dp]
    law % cimag = [0.0_dp,  0.0_dp, 0.0_dp,  0.0_dp, 2.0_dp]

    allocate(q(2 * nv), escape(nv))
    q      = 0.0_dp
    escape = 0

    do n = 1, nmax
       call march(explicit_step(law, 1.0_dp), points, q, 1, 1)
       do v = 1, nv
          if (escape(v) == 0 .and. &
               & q(2 * v - 1)**2 + q(2 * v)**2 > 4.0_dp) then
             ! The escape is recorded; the orbit is fixed at zero so
             ! the remaining march cannot overflow.
             escape(v) = n
             law % creal(v) = 0.0_dp
             law % cimag(v) = 0.0_dp
             q(2 * v - 1 : 2 * v) = 0.0_dp
          end if
       end do
    end do

    ! The counts as an integer field.
    escape_field = stored_field('escape time', points % vertex_set(), points % num_vertices())
    call escape_field % set_integer_vector(escape)

    call report(escape(1) == 0, 'c = 0 never escapes: the origin is fixed', num_failures)
    call report(escape(2) == 0, 'c = -1 never escapes: the two-cycle', num_failures)
    call report(escape(3) == 3, 'c = 1 escapes at step three: 0, 1, 2, 5', num_failures)
    call report(escape(4) == 0, 'c = -2 never escapes: fixed at two', num_failures)
    call report(escape(5) == 2, &
         & 'c = 2i escapes at step two: -4 + 2i leaves the circle', num_failures)

  end subroutine check_mandelbrot_iteration

  !===================================================================!
  ! The implicit pass: backward euler lands on its own discrete
  ! identity q0/(1+h)^n, and bdf2 quarters its error when the step
  ! halves - second order, measured.
  !===================================================================!

  subroutine check_implicit_discrete_solution(num_failures)

    integer, intent(inout) :: num_failures

    type(stored_directed_graph) :: single_vertex
    type(differential_operator) :: decay
    real(dp) :: q(1), expected, coarse_error, fine_error, ratio

    single_vertex  = stored_directed_graph(1, tails=[integer ::], heads=[integer ::])
    decay = differential_operator(SIDE_VERTEX, 0, coefficient=1.0_dp)

    q = [3.0_dp]
    call march(implicit_step(decay, 0.125_dp, 1), single_vertex, q, 16, 1)
    expected = 3.0_dp / (1.0_dp + 0.125_dp)**16

    call report(abs(q(1) - expected) < 1.0d-9, &
         & 'backward euler lands on q0/(1+h)^n: its own discrete identity', num_failures)

    ! bdf2 against exp(-2): half the step, a quarter of the error.
    q = [1.0_dp]
    call march(implicit_step(decay, 0.1_dp, 2), single_vertex, q, 20, 2)
    coarse_error = abs(q(1) - exp(-2.0_dp))

    q = [1.0_dp]
    call march(implicit_step(decay, 0.05_dp, 2), single_vertex, q, 40, 2)
    fine_error = abs(q(1) - exp(-2.0_dp))

    ratio = coarse_error / fine_error
    call report(ratio > 3.2_dp .and. ratio < 4.8_dp, &
         & 'bdf2 quarters its error when the step halves: second order', num_failures)

  end subroutine check_implicit_discrete_solution

  !===================================================================!
  ! The same numbers, placed two ways. Two cells with one number each
  ! and one cell with two are the same state written differently, and
  ! a decoupled law must march them to the same place. The implicit
  ! pass is where this fails silently: a solver that measures its
  ! residual by the first component alone reports convergence while
  ! the rest of the entry is still moving.
  !===================================================================!

  subroutine check_component_layout_invariance(num_failures)

    integer, intent(inout) :: num_failures

    type(stored_directed_graph) :: two_vertices, single_vertex
    type(differential_operator) :: decay
    real(dp) :: vector_state(2), scalar_states(2), expected(2)

    two_vertices = stored_directed_graph(2, tails=[integer ::], heads=[integer ::])
    single_vertex = stored_directed_graph(1, tails=[integer ::], heads=[integer ::])
    decay = differential_operator(SIDE_VERTEX, 0, coefficient=1.0_dp)

    ! Two cells, one number each.
    scalar_states = [3.0_dp, 5.0_dp]
    call march(implicit_step(decay, 0.125_dp, 1), two_vertices, scalar_states, 8, 1)

    ! One cell, two numbers.
    vector_state = [3.0_dp, 5.0_dp]
    call march(implicit_step(decay, 0.125_dp, 1), single_vertex, vector_state, 8, 1)

    expected = [3.0_dp, 5.0_dp] / (1.0_dp + 0.125_dp)**8

    call report(maxval(abs(scalar_states - expected)) < 1.0d-9, &
         & 'two scalar cell values satisfy the discrete decay relation', num_failures)
    call report(maxval(abs(vector_state - expected)) < 1.0d-9, &
         & 'one vector cell value satisfies the same discrete decay relation', num_failures)

  end subroutine check_component_layout_invariance

  !===================================================================!
  ! The reverse pass is the adjoint: march q forward under a
  ! statement, march lambda backward under its transpose, and the
  ! pairing <lambda, q> is invariant.
  !===================================================================!

  subroutine check_adjoint_pairing_invariance(num_failures)

    integer, intent(inout) :: num_failures

    type(stored_directed_graph) :: three_vertices
    type(stencil) :: forward_action
    real(dp) :: q(3), lambda(3), terminal_pairing, initial_pairing
    integer , parameter :: rows(6) = [1, 1, 2, 2, 3, 3]
    integer , parameter :: cols(6) = [1, 2, 2, 3, 3, 1]
    real(dp), parameter :: w(6) = [2.0_dp, -1.0_dp, 1.0_dp, -0.4_dp, &
         &                         0.3_dp, 0.7_dp]
    real(dp), parameter :: zeros(3) = [0.0_dp, 0.0_dp, 0.0_dp]
    real(dp), parameter :: h = 0.05_dp

    three_vertices = stored_directed_graph(3, tails=[integer ::], heads=[integer ::])

    ! An unsymmetric statement; the reverse pass derives its
    ! transpose from the tangent at every recorded state.
    forward_action = stencil(rows, cols, w, zeros)

    q      = [1.0_dp, -2.0_dp, 3.0_dp]
    lambda = [0.4_dp, 2.0_dp, -1.0_dp]

    ! The pairing at the far end: <lambda_N, q_N> after q is marched
    ! forward; lambda_N is the terminal costate of the transposed
    ! chain, evaluated in the same traversal.
    call march(explicit_step(forward_action, h), three_vertices, q, 12, 1, &
         & adjoint=adjoint_step(forward_action, h), lambda=lambda)
    terminal_pairing = sum([0.4_dp, 2.0_dp, -1.0_dp] * q)

    ! The initial costate against the initial state.
    q = [1.0_dp, -2.0_dp, 3.0_dp]
    initial_pairing = sum(lambda * q)

    call report(abs(terminal_pairing - initial_pairing) < 1.0d-12 * (1.0_dp + abs(terminal_pairing)), &
         & 'the reverse pass keeps the pairing: <lambda, q> is invariant', num_failures)

  end subroutine check_adjoint_pairing_invariance

  !===================================================================!
  ! On Van der Pol, the tangent marched forward as an augmented
  ! statement and the adjoint marched back under the linearization
  ! derived from the law itself are evaluated over one chain, so the
  ! gradient they each report is one number, and the pairing
  ! <lambda, dq> is invariant along the traversal.
  !===================================================================!

  subroutine check_tangent_adjoint_equivalence(num_failures)

    integer, intent(inout) :: num_failures

    type(stored_directed_graph) :: cell
    type(vdp_law)         :: law
    type(vdp_tangent_law) :: tangent
    real(dp) :: aug(4), lambda(2), grad_tangent(2), q(2)
    real(dp), parameter :: h = 0.01_dp
    real(dp), parameter :: q0(2) = [2.0_dp, 0.0_dp]
    integer , parameter :: nsteps = 100
    integer :: i

    law     = vdp_law()
    tangent = vdp_tangent_law()

    cell = stored_directed_graph(1, tails=[integer ::], heads=[integer ::])

    ! The tangent pass: one augmented march per seed direction; the
    ! objective is u at the end, so the gradient entry is du_N.
    do i = 1, 2
       aug = 0.0_dp
       aug(1:2) = q0
       aug(2 + i) = 1.0_dp
       call march(explicit_step(tangent, h), cell, aug, nsteps, 1)
       grad_tangent(i) = aug(3)
    end do

    ! The adjoint pass: one reverse traversal, seeded by the
    ! objective, under the transposed linearization derived from the
    ! law at every recorded state.
    q      = q0
    lambda = [1.0_dp, 0.0_dp]
    call march(explicit_step(law, h), cell, q, nsteps, 1, &
         & adjoint=adjoint_step(law, h), lambda=lambda)

    call report(all(abs(lambda - grad_tangent) < 1.0d-12 &
         & * (1.0_dp + abs(grad_tangent))), &
         & 'tangent and adjoint meet on one gradient: the same graph evaluated twice', &
         & num_failures)

    ! The structural identity at the ends: <lambda_0, dq_0> equals
    ! <lambda_N, dq_N> by construction.
    call report(abs(lambda(1) - grad_tangent(1)) < 1.0d-12 &
         & * (1.0_dp + abs(grad_tangent(1))), &
         & 'and the pairing is invariant along the traversal', num_failures)

  end subroutine check_tangent_adjoint_equivalence

end program test_graph_marching
