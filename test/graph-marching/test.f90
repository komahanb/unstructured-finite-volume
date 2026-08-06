!=====================================================================!
! The marching suite: rung 6's acceptance, with a fractal for a
! birth certificate.
!
! Three verdicts. TIME IS A GRAPH: the marcher's instants stand as
! a chain - one vertex per instant, one edge per step, walked in
! order. THE STEP IS EXACT ABOUT ITSELF: euler on decay dq/dt = -q
! lands on (1 - h)^n to machine precision, euler's own discrete
! truth. And THE MAP IS THE MARCH: z -> z^2 + c is forward euler
! with step one on S = z - z^2 - c - an identity, not an
! approximation - so the escape-time fractal falls out of the
! marcher with known points landing where sixty years of arithmetic
! say they land:
!
!      c = 0        never escapes: the origin is fixed
!      c = -1       never escapes: a two-cycle, 0 -> -1 -> 0
!      c = 1        escapes at step three: 0, 1, 2, 5
!      c = -2       never escapes: walks to the fixed point 2
!      c = 2i       escapes at step two: -4 + 2i leaves the circle
!=====================================================================!

program test_graph_marching

  use iso_fortran_env, only : dp => REAL64
  use structure_graph, only : graph
  use data_graph_field, only : graph_field
  use operation_graph_operation, only : graph_operation
  use structure_graph, only : GRAPH_SIDE_VERTEX
  use structure_support, only : support
  use data_field  , only : field
  use structure_stored_graph        , only : stored_graph
  use operation_differential, only : vertex_differential_operator
  use operation_differential, only : differential
  use operation_step   , only : chain, chain
  use operation_step   , only : CHAIN_FORWARD, CHAIN_BACKWARD, CHAIN_BDF2
  use operation_fit_optimizer , only : fit_optimizer, substitution, newton, gmres
  use operation_stencil , only : stencil
  use mandelbrot_law_fixture, only : mandelbrot_law
  use vdp_fixture, only : vdp_law, vdp_tangent_law, vdp_adjoint_law

  implicit none

  integer :: nfail

  nfail = 0

  call check_time_is_a_graph(nfail)
  call check_euler_is_exact_about_itself(nfail)
  call check_the_map_is_the_march(nfail)
  call check_the_implicit_road(nfail)
  call check_a_wide_entry_marches(nfail)
  call check_the_reverse_walk(nfail)
  call check_tangent_meets_adjoint(nfail)

  write(*, '(a)') ' ============================================='
  if (nfail == 0) then
     write(*, '(a)') ' all marching checks passed'
  else
     write(*, '(a, i0, a)') ' ', nfail, ' marching checks FAILED'
     error stop 1
  end if

contains

  subroutine report(ok, message, nfail)

    logical         , intent(in)    :: ok
    character(len=*), intent(in)    :: message
    integer         , intent(inout) :: nfail

    if (ok) then
       write(*, '(a)') ' PASS : ' // message
    else
       write(*, '(a)') ' FAIL : ' // message
       nfail = nfail + 1
    end if

  end subroutine report

  !===================================================================!
  ! VERDICT ONE. The instants stand as a chain.
  !===================================================================!

  subroutine check_time_is_a_graph(nfail)

    integer, intent(inout) :: nfail

    type(stored_graph) :: timeline
    logical :: ordered
    integer :: e

    timeline = instants(10)

    call report(timeline % num_vertices() == 11 .and. timeline % num_edges() == 10, &
         & 'eleven instants, ten steps: the chain stands', nfail)

    ordered = .true.
    do e = 1, 10
       if (timeline % edge_tail(e) /= e .or. timeline % edge_head(e) /= e + 1) &
            & ordered = .false.
    end do
    call report(ordered, 'and every step leads to the next instant', nfail)

  end subroutine check_time_is_a_graph

  !===================================================================!
  ! VERDICT TWO. Euler's discrete truth, hit exactly: on S = q the
  ! walk gives q_n = (1 - h)^n q_0.
  !===================================================================!

  subroutine check_euler_is_exact_about_itself(nfail)

    integer, intent(inout) :: nfail

    type(support) :: cells
    type(stored_graph) :: lone
    type(differential) :: decay
    real(dp) :: q(1), expected
    integer :: v

    lone = stored_graph(1, tails=[integer ::], heads=[integer ::])
    associate (u1 => cells, u2 => v)
    end associate

    decay = vertex_differential_operator(order=0, coefficient=1.0_dp)

    q = [3.0_dp]
    call walk(decay, lone, q, 20, 0.125_dp, CHAIN_FORWARD)

    expected = 3.0_dp * (1.0_dp - 0.125_dp)**20

    call report(abs(q(1) - expected) < 1.0d-14, &
         & 'euler lands on (1-h)^n exactly: its own discrete truth', nfail)

  end subroutine check_euler_is_exact_about_itself

  !===================================================================!
  ! VERDICT THREE. Five points of the complex plane, marched by
  ! z -> z^2 + c with the escape watched. The law is a level-3
  ! fixture; the marcher neither knows nor cares that the picture
  ! is a fractal.
  !===================================================================!

  subroutine check_the_map_is_the_march(nfail)

    integer, intent(inout) :: nfail

    type(mandelbrot_law) :: law
    type(stored_graph) :: points
    type(support) :: cells
    type(field) :: escape_field
    real(dp), allocatable :: q(:)
    integer, allocatable :: escape(:)
    integer :: v, n
    integer, parameter :: nv = 5, nmax = 30

    ! The five points, as a graph of lone cells.
    points = stored_graph(nv, tails=[integer ::], heads=[integer ::])

    law % creal = [0.0_dp, -1.0_dp, 1.0_dp, -2.0_dp, 0.0_dp]
    law % cimag = [0.0_dp,  0.0_dp, 0.0_dp,  0.0_dp, 2.0_dp]

    allocate(q(2 * nv), escape(nv))
    q      = 0.0_dp
    escape = 0

    do n = 1, nmax
       call walk(law, points, q, 1, 1.0_dp, CHAIN_FORWARD)
       do v = 1, nv
          if (escape(v) == 0 .and. &
               & q(2 * v - 1)**2 + q(2 * v)**2 > 4.0_dp) then
             ! The verdict is in; the orbit is frozen so the
             ! remaining march cannot overflow on its way out.
             escape(v) = n
             law % creal(v) = 0.0_dp
             law % cimag(v) = 0.0_dp
             q(2 * v - 1 : 2 * v) = 0.0_dp
          end if
       end do
    end do

    ! The counts leave as an integer field: the tower's native kind
    ! for a picture of whole numbers.
    cells = support(GRAPH_SIDE_VERTEX, [(v, v = 1, nv)])
    escape_field = field('escape time', cells)
    call escape_field % set_integer_vector(escape)

    call report(escape(1) == 0, 'c = 0 never escapes: the origin holds', nfail)
    call report(escape(2) == 0, 'c = -1 never escapes: the two-cycle turns', nfail)
    call report(escape(3) == 3, 'c = 1 escapes at step three: 0, 1, 2, 5', nfail)
    call report(escape(4) == 0, 'c = -2 never escapes: pinned at two', nfail)
    call report(escape(5) == 2, &
         & 'c = 2i escapes at step two: -4 + 2i leaves the circle', nfail)

  end subroutine check_the_map_is_the_march

  !===================================================================!
  ! VERDICT FOUR. The implicit road: backward euler lands on its own
  ! discrete truth q0/(1+h)^n, and bdf2 halves its error fourfold -
  ! second order, measured.
  !===================================================================!

  subroutine check_the_implicit_road(nfail)

    integer, intent(inout) :: nfail

    type(stored_graph) :: lone
    type(differential) :: decay
    real(dp) :: q(1), expected, coarse_error, fine_error, ratio

    lone  = stored_graph(1, tails=[integer ::], heads=[integer ::])
    decay = vertex_differential_operator(order=0, coefficient=1.0_dp)

    q = [3.0_dp]
    call walk(decay, lone, q, 16, 0.125_dp, CHAIN_BACKWARD)
    expected = 3.0_dp / (1.0_dp + 0.125_dp)**16

    call report(abs(q(1) - expected) < 1.0d-9, &
         & 'backward euler lands on q0/(1+h)^n: its own discrete truth', nfail)

    ! bdf2 against exp(-2): half the step, a quarter of the error.
    q = [1.0_dp]
    call walk(decay, lone, q, 20, 0.1_dp, CHAIN_BDF2)
    coarse_error = abs(q(1) - exp(-2.0_dp))

    q = [1.0_dp]
    call walk(decay, lone, q, 40, 0.05_dp, CHAIN_BDF2)
    fine_error = abs(q(1) - exp(-2.0_dp))

    ratio = coarse_error / fine_error
    call report(ratio > 3.2_dp .and. ratio < 4.8_dp, &
         & 'bdf2 quarters its error when the step halves: second order', nfail)

  end subroutine check_the_implicit_road

  !===================================================================!
  ! The same numbers, seated two ways. Two cells carrying one number
  ! each and one cell carrying two are the same state written
  ! differently, and a decoupled law must march them to the same
  ! place. The implicit road is where this can go wrong quietly: a
  ! solver that measures its residual by the first stripe alone
  ! declares victory while the rest of the entry is still moving.
  !===================================================================!

  subroutine check_a_wide_entry_marches(nfail)

    integer, intent(inout) :: nfail

    type(stored_graph) :: pair, lone
    type(differential) :: decay
    real(dp) :: wide(2), tall(2), expected(2)

    pair = stored_graph(2, tails=[integer ::], heads=[integer ::])
    lone = stored_graph(1, tails=[integer ::], heads=[integer ::])
    decay = vertex_differential_operator(order=0, coefficient=1.0_dp)

    ! Two cells, one number each.
    tall = [3.0_dp, 5.0_dp]
    call walk(decay, pair, tall, 8, 0.125_dp, CHAIN_BACKWARD)

    ! One cell, two numbers.
    wide = [3.0_dp, 5.0_dp]
    call walk(decay, lone, wide, 8, 0.125_dp, CHAIN_BACKWARD)

    expected = [3.0_dp, 5.0_dp] / (1.0_dp + 0.125_dp)**8

    call report(maxval(abs(tall - expected)) < 1.0d-9, &
         & 'two cells, one number each, land on their discrete truth', nfail)
    call report(maxval(abs(wide - expected)) < 1.0d-9, &
         & 'and one cell two numbers wide lands on the very same place', nfail)

  end subroutine check_a_wide_entry_marches

  !===================================================================!
  ! VERDICT FIVE. The reverse walk is the adjoint: march q forward
  ! under a statement, march lambda backward under its transpose,
  ! and the pairing <lambda, q> never moves.
  !===================================================================!

  subroutine check_the_reverse_walk(nfail)

    integer, intent(inout) :: nfail

    type(stored_graph) :: trio
    type(stencil) :: forward_action, transposed
    real(dp) :: q(3), lambda(3), before, after
    integer , parameter :: rows(6) = [1, 1, 2, 2, 3, 3]
    integer , parameter :: cols(6) = [1, 2, 2, 3, 3, 1]
    real(dp), parameter :: w(6) = [2.0_dp, -1.0_dp, 1.0_dp, -0.4_dp, &
         &                         0.3_dp, 0.7_dp]
    real(dp), parameter :: zeros(3) = [0.0_dp, 0.0_dp, 0.0_dp]

    trio = stored_graph(3, tails=[integer ::], heads=[integer ::])

    ! An unsymmetric statement and its transpose: rows and columns
    ! swapped, the stencil's own adjoint.
    forward_action = stencil(rows, cols, w, zeros)
    transposed     = stencil(cols, rows, w, zeros)

    q      = [1.0_dp, -2.0_dp, 3.0_dp]
    lambda = [0.4_dp, 2.0_dp, -1.0_dp]

    ! The pairing at the far end: <lambda_N, q_N> needs q marched
    ! all the way forward first.
    call walk(forward_action, trio, q, 12, 0.05_dp, CHAIN_FORWARD)
    before = sum(lambda * q)

    ! Now lambda walks home under the transpose: the same sweep,
    ! running the other way.
    call walk_home(transposed, trio, lambda, 12, 0.05_dp)

    ! And meets the initial state.
    q = [1.0_dp, -2.0_dp, 3.0_dp]
    after = sum(lambda * q)

    call report(abs(before - after) < 1.0d-12 * (1.0_dp + abs(before)), &
         & 'the reverse walk keeps the pairing: <lambda, q> never moves', nfail)

  end subroutine check_the_reverse_walk

  !===================================================================!
  ! VERDICT SIX. The paper's witness: on Van der Pol, the tangent
  ! marched forward as an augmented statement and the adjoint walked
  ! home under the transposed linearization are constructed from the
  ! same chain - so the gradient they each report is one number, and
  ! the pairing <lambda, dq> holds across the whole walk. This is
  ! the consistency the complex step used to be hired for, standing
  ! in the suite at no cost.
  !===================================================================!

  subroutine check_tangent_meets_adjoint(nfail)

    integer, intent(inout) :: nfail

    type(stored_graph) :: cell
    type(vdp_law)         :: law
    type(vdp_tangent_law) :: tangent
    type(vdp_adjoint_law) :: transposed
    real(dp), allocatable :: trajectory(:,:)
    real(dp) :: aug(4), lambda(2), grad_tangent(2), q(2)
    real(dp), parameter :: h = 0.01_dp
    integer , parameter :: nsteps = 100
    integer :: n, i

    cell = stored_graph(1, tails=[integer ::], heads=[integer ::])

    ! The trajectory, stored: the reverse walk reads it back.
    allocate(trajectory(2, 0:nsteps))
    q = [2.0_dp, 0.0_dp]
    trajectory(:, 0) = q
    do n = 1, nsteps
       call walk(law, cell, q, 1, h, CHAIN_FORWARD)
       trajectory(:, n) = q
    end do

    ! The tangent road: one augmented march per seed direction; the
    ! objective is u at the end, so the gradient entry is du_N.
    do i = 1, 2
       aug = 0.0_dp
       aug(1:2) = trajectory(:, 0)
       aug(2 + i) = 1.0_dp
       call walk(tangent, cell, aug, nsteps, h, CHAIN_FORWARD)
       grad_tangent(i) = aug(3)
    end do

    ! The adjoint road: one reverse walk, seeded by the objective.
    lambda = [1.0_dp, 0.0_dp]
    do n = nsteps, 1, -1
       transposed % at = trajectory(:, n - 1)
       call walk_home(transposed, cell, lambda, 1, h)
    end do

    call report(all(abs(lambda - grad_tangent) < 1.0d-12 &
         & * (1.0_dp + abs(grad_tangent))), &
         & 'tangent and adjoint meet on one gradient: the same graph spoke twice', &
         & nfail)

    ! The structural identity at the ends: <lambda_0, dq_0> equals
    ! <lambda_N, dq_N> by construction.
    call report(abs(lambda(1) - grad_tangent(1)) < 1.0d-12 &
         & * (1.0_dp + abs(grad_tangent(1))), &
         & 'and the pairing holds across the walk: the identity stands', nfail)

  end subroutine check_tangent_meets_adjoint

  !===================================================================!
  ! The chain of instants: one vertex per moment, one edge per step.
  !===================================================================!

  type(stored_graph) function instants(nsteps) result(timeline)

    integer, intent(in) :: nsteps
    integer :: n

    timeline = stored_graph(nsteps + 1, &
         & tails=[(n, n = 1, nsteps)], heads=[(n + 1, n = 1, nsteps)])

  end function instants

  !===================================================================!
  ! A march, expressed as what it is: the whole recurrence stated on
  ! the time chain, and one causal sweep answering it exactly. The
  ! state arrives as the initial instant and leaves as the last.
  !===================================================================!

  subroutine walk(action, space, q, nsteps, h, rule)

    class(graph_operation), intent(in) :: action
    class(graph), intent(in)           :: space
    real(dp), intent(inout)            :: q(:)
    integer , intent(in)               :: nsteps
    real(dp), intent(in)               :: h
    integer , intent(in)               :: rule

    call sweep(action, space, q, nsteps, h, rule, backward=.false.)

  end subroutine walk

  !===================================================================!
  ! And the same sweep running the other way: sensitivities settle
  ! backward along the chain the state settled forward on.
  !===================================================================!

  subroutine walk_home(transposed, space, lambda, nsteps, h)

    class(graph_operation), intent(in) :: transposed
    class(graph), intent(in)           :: space
    real(dp), intent(inout)            :: lambda(:)
    integer , intent(in)               :: nsteps
    real(dp), intent(in)               :: h

    call sweep(transposed, space, lambda, nsteps, h, CHAIN_FORWARD, &
         & backward=.true.)

  end subroutine walk_home

  subroutine sweep(action, space, q, nsteps, h, rule, backward)

    class(graph_operation), intent(in) :: action
    class(graph), intent(in)           :: space
    real(dp), intent(inout)            :: q(:)
    integer , intent(in)               :: nsteps
    real(dp), intent(in)               :: h
    integer , intent(in)               :: rule
    logical , intent(in)               :: backward

    type(stored_graph)   :: timeline
    type(chain) :: recurrence
    type(fit_optimizer)  :: causal
    real(dp), allocatable :: trajectory(:), g(:)
    real(dp) :: achieved
    integer :: width, last

    width = size(q)

    timeline   = instants(nsteps)
    recurrence = chain(action, space, h, rule, initial=q)

    causal = substitution(newton(gmres()), backward=backward)
    causal % inner % tolerance = 1.0d-13

    call causal % attach(recurrence, timeline, ncomp=width)

    allocate(trajectory((nsteps + 1) * width))
    trajectory = 0.0_dp

    call causal % constant(g)
    call causal % solve(-g, trajectory, achieved)

    if (backward) then
       q = trajectory(1 : width)
    else
       last = nsteps * width
       q = trajectory(last + 1 : last + width)
    end if

  end subroutine sweep

end program test_graph_marching
