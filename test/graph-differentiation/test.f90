!=====================================================================!
! Tests for the differentiation stack:
!
!      set_bdf            the nonuniform and uniform order-2 rows
!                         and the order-1 row
!      tangent_of         exact/difference dispatch, and the value
!                         of the exact tangent
!      the scheme         its max_degree and partial actions, the
!                         exact tangent of the step equation, the
!                         chain rule over it
!      chain_rule         total derivatives to degree 8
!      the marcher        march with trajectory recording,
!                         march_directional to degree 8,
!                         march_adjoint with terminal and
!                         per-instant seeds, march_adaptive, a
!                         nonuniform implicit march, and the
!                         adjoint on the difference road
!
! Every expected value is an exact rational or a closed form
! derived independently of the code under test; the derivation is
! stated in each section comment.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_graph_differentiation

  use iso_fortran_env     , only : dp => REAL64
  use field_calculus, only : field
  use graph_fractal       , only : graph
  use view_directed_stored         , only : stored_directed_graph
  use field_stored   , only : stored_field
  use operation_action  , only : variation
  use operation_step    , only : scheme, backward_euler, bdf_variable
  use operation_chain_rule, only : chain_rule, argument_path
  use operation_linearization, only : linearization, tangent_of, dual_by_basis
  use operation_stencil , only : stencil
  use operation_marching , only : marcher, MARCH_BACKWARD
  use operation_step_policy, only : halving_policy
  use operation_newton  , only : newton
  use operation_gmres   , only : gmres
  use toy_differentiable_forms, only : quartic_form, power8_form, &
       & equilibrium_law, linear_law, scalar_pair, fill_path

  implicit none

  type(quartic_form)    :: quartic
  type(power8_form)     :: p8
  type(equilibrium_law) :: equil
  type(linear_law)      :: lin

  type(stored_directed_graph) :: lone
  type(graph)             :: cells
  type(marcher)               :: clock

  integer :: nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "graph differentiation tower suite"
  write(*,'(1x,a)') "============================================="

  lone  = stored_directed_graph(1, tails=[integer ::], heads=[integer ::])
  cells = lone % vertex_set()

  ! every fixture is built by its constructor, which declares its
  ! arguments: q first, xi second
  quartic = quartic_form()
  p8      = power8_form()
  equil   = equilibrium_law()
  lin     = linear_law()

  clock % rule = MARCH_BACKWARD
  clock % step = 1.0_dp
  allocate(clock % inner, source=newton())
  select type (nw => clock % inner)
  type is (newton)
     allocate(nw % inner, source=gmres())
     nw % inner % tolerance = 1.0d-14
     nw % tolerance = 1.0d-12
  end select

  call check_the_table_seat(nfail)
  call check_the_tangent_chooser(nfail)
  call check_the_step_tangent(nfail)
  call check_the_chain_rule(nfail)
  call check_the_taylor_convolution(nfail)
  call check_the_derivative_walks(nfail)
  call check_the_adaptive_walk(nfail)
  call check_the_argument_calculus(nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all differentiation tower checks passed"
  else
     error stop
  end if

contains

  !===================================================================!
  ! set_bdf order 2 at steps (2, 3): with h0 = 2 and h1 = 3 the
  ! coefficients are a0 = (2 h0 + h1)/(h0 + h1) = 7/5,
  ! a1 = -(h0 + h1)/h1 = -5/3, a2 = h0^2/(h1 (h0 + h1)) = 4/15.
  ! Equal steps must give the constants (1.5, -2, 0.5) exactly,
  ! because set_bdf assigns them without arithmetic. Order 1 is
  ! (1, -1, 0) for any step.
  !===================================================================!

  subroutine check_the_table_seat(nfail)

    integer, intent(inout) :: nfail

    type(scheme) :: statement

    statement = bdf_variable(2, lin, [2.0_dp, 3.0_dp])
    call report( &
         & near(statement % a0,  7.0_dp /  5.0_dp, 1.0e-15_dp) .and. &
         & near(statement % a1, -5.0_dp /  3.0_dp, 1.0e-15_dp) .and. &
         & near(statement % a2,  4.0_dp / 15.0_dp, 1.0e-15_dp) .and. &
         & near(statement % hs,  2.0_dp, 0.0_dp), &
         & "set_bdf at steps (2, 3) gives [7/5, -5/3, 4/15], h = 2", nfail)

    call statement % set_bdf(2, [0.7_dp, 0.7_dp])
    call report( &
         & statement % a0 ==  1.5_dp .and. &
         & statement % a1 == -2.0_dp .and. &
         & statement % a2 ==  0.5_dp, &
         & "equal steps give the uniform coefficients exactly", nfail)

    statement = backward_euler(lin, 0.25_dp)
    call report( &
         & statement % a0 ==  1.0_dp .and. &
         & statement % a1 == -1.0_dp .and. &
         & statement % a2 ==  0.0_dp .and. &
         & statement % reach == 1, &
         & "backward euler gives [1, -1, 0], reach 1", nfail)

  end subroutine check_the_table_seat

  !===================================================================!
  ! tangent_of must take the exact road when the operation's
  ! max_degree is at least one and the difference road otherwise;
  ! the two are distinguished here by their name() prefixes. The
  ! exact tangent of the quartic frozen at
  ! q = 1 (xi defaulting to 2) is Phi_q(1, 2) = 26, so applying it
  ! to the direction v = 3 must return 78.
  !===================================================================!

  subroutine check_the_tangent_chooser(nfail)

    integer, intent(inout) :: nfail

    type(linearization) :: tangent, slow
    class(field), allocatable :: output
    type(stored_field) :: direction
    real(dp), allocatable :: rv(:)

    tangent = tangent_of(quartic)
    slow    = tangent_of(lin)

    call report(lin % max_degree() == 0 .and. quartic % max_degree() == 4, &
         & "max_degree is 0 unless the operation differentiates itself", nfail)
    call report(index(tangent % name(), 'exact derivative of') == 1, &
         & "tangent_of picks the exact linearization when differentiable", &
         & nfail)
    call report(index(slow % name(), 'derivative of') == 1 .and. &
         & index(slow % name(), 'exact') == 0, &
         & "tangent_of picks the difference linearization otherwise", nfail)

    call tangent % freeze([1.0_dp])

    direction = stored_field('v', cells, 1, num_components=1)
    call direction % set_real_vector([3.0_dp])
    call tangent % apply(lone, [direction], output)
    call output % real_vector(rv)

    call report(size(rv) == 1 .and. near(rv(1), 78.0_dp, 1.0e-12_dp), &
         & "the exact tangent of the quartic at 1: J v = 26 v", nfail)

  end subroutine check_the_tangent_chooser

  !===================================================================!
  ! The scheme's own calculus. Its max_degree is the action's. Its
  ! partial actions follow from a0 q + a1 qold + a2 qolder + hs S(q):
  ! hs times the action's partial, plus a0 v for the first
  ! derivative in slot 1 only. With the quartic at (q, xi) = (1, 2)
  ! - Phi_q = 26, Phi_xi = 49, Phi_qq = 32 - and backward euler at
  ! h = 1/2 (a0 = 1, hs = 1/2):
  !
  !      D_q [3]      = 1*3 + 26*3/2   = 42
  !      D_xi [1]     = 49/2           = 24.5
  !      D_qq [1, 1]  = 32/2           = 16
  !
  ! and after set_bdf(2, [2, 3]) (a0 = 7/5, hs = 2) the tangent at
  ! [3] is 21/5 + 156 = 160.2. The chain rule over the scheme with
  ! the quartic's paths q^(k) = (1, 2, 3, 4), xi^(k) = (5, 7, 11,
  ! 13) gives a0 q^(n) + hs d_n = n + d_n/2 for n = 1..4 and, at
  ! degree 0 with qold = 0, 1 + 31/2 = 16.5. A scheme over
  ! linear_law keeps the difference tangent, which reads qold: at
  ! h = 1 it is a0 + h = 2.
  !===================================================================!

  subroutine check_the_step_tangent(nfail)

    integer, intent(inout) :: nfail

    type(scheme) :: statement, shallow, slowstep
    type(linearization) :: tangent, slow
    type(chain_rule)    :: composer
    type(argument_path) :: full(2)
    type(stored_field)  :: inputs(2), v, w, one
    class(field), allocatable :: output
    real(dp), allocatable :: rv(:)
    real(dp) :: expected(0:4)
    logical  :: ok
    integer  :: n

    statement = backward_euler(equil, 1.0_dp)
    shallow   = bdf_variable(2, lin, [1.0_dp, 1.0_dp])
    call report(statement % max_degree() == 8 .and. &
         & shallow % max_degree() == 0, &
         & "a scheme's max_degree is its action's", nfail)

    statement = backward_euler(quartic, 0.5_dp)
    statement % qold = [0.0_dp]
    call scalar_pair(1.0_dp, 2.0_dp, cells, inputs)

    v = stored_field('v', cells, 1, num_components=1)
    call v % set_real_vector([3.0_dp])
    w = stored_field('w', cells, 1, num_components=1)
    call w % set_real_vector([1.0_dp])
    one = stored_field('one', cells, 1, num_components=1)
    call one % set_real_vector([1.0_dp])

    call statement % partial_action(lone, inputs, &
         & [variation(statement % state(), v)], output)
    call output % real_vector(rv)
    ok = size(rv) == 1 .and. near(rv(1), 42.0_dp, 1.0e-12_dp)
    call statement % partial_action(lone, inputs, &
         & [variation(statement % auxiliary(1), w)], output)
    call output % real_vector(rv)
    ok = ok .and. near(rv(1), 24.5_dp, 1.0e-12_dp)
    call statement % partial_action(lone, inputs, &
         & [variation(statement % state(), one), variation(statement % state(), one)], output)
    call output % real_vector(rv)
    ok = ok .and. near(rv(1), 16.0_dp, 1.0e-12_dp)
    call report(ok, &
         & "the scheme's partials a0 v + hs D_q S, hs D_xi S, hs D_qq S: &
         &42, 24.5, 16", nfail)

    tangent = tangent_of(statement)
    call tangent % freeze([1.0_dp])
    call tangent % apply(lone, [v], output)
    call output % real_vector(rv)
    ok = index(tangent % name(), 'exact derivative of') == 1 .and. &
         & near(rv(1), 42.0_dp, 1.0e-12_dp)

    call statement % set_bdf(2, [2.0_dp, 3.0_dp])
    tangent = tangent_of(statement)
    call tangent % freeze([1.0_dp])
    call tangent % apply(lone, [v], output)
    call output % real_vector(rv)
    ok = ok .and. near(rv(1), 160.2_dp, 1.0e-12_dp)
    call report(ok, &
         & "tangent_of(scheme) is exact: 42 at backward euler, 160.2 after &
         &set_bdf(2, [2, 3])", nfail)

    call statement % set_bdf(1, [0.5_dp])
    call fill_path(full(1), statement % state(), [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], cells)
    call fill_path(full(2), statement % auxiliary(1), [5.0_dp, 7.0_dp, 11.0_dp, 13.0_dp], cells)
    expected = [16.5_dp, 136.5_dp, 1105.5_dp, 8347.0_dp, 59129.5_dp]
    ok = .true.
    do n = 0, 4
       call composer % assemble(statement, lone, inputs, n, full, output)
       call output % real_vector(rv)
       ok = ok .and. size(rv) == 1 .and. near(rv(1), expected(n), 1.0e-10_dp)
    end do
    call report(ok, &
         & "the chain rule over the scheme, a0 q^(n) + hs d_n: 16.5, 136.5, &
         &1105.5, 8347, 59129.5", nfail)

    slowstep = backward_euler(lin, 1.0_dp)
    slowstep % qold = [0.0_dp]
    slow = tangent_of(slowstep)
    call slow % freeze([1.0_dp])
    call slow % apply(lone, [w], output)
    call output % real_vector(rv)
    call report(index(slow % name(), 'exact') == 0 .and. &
         & near(rv(1), 2.0_dp, 1.0e-6_dp), &
         & "a scheme over a max_degree-0 action keeps the difference tangent: &
         &a0 + h = 2", nfail)

  end subroutine check_the_step_tangent

  !===================================================================!
  ! chain_rule % assemble on the quartic at (q, xi) = (1, 2) with
  ! path derivatives q^(k) = (1, 2, 3, 4) and xi^(k) = (5, 7, 11,
  ! 13). Expected totals were computed by rational Taylor
  ! composition, independently of the pattern generator:
  !
  !      d0 = 31,  d1 = 271,  d2 = 2207,  d3 = 16688,  d4 = 118251
  !
  ! With xi's derivatives past the first unoccupied, those terms
  ! are read as zero: d3 = 9156, d4 = 27908.
  !===================================================================!

  subroutine check_the_chain_rule(nfail)

    integer, intent(inout) :: nfail

    type(chain_rule)    :: composer
    type(argument_path) :: full(2), sparse(2)
    type(stored_field)         :: inputs(2)
    class(field), allocatable :: output
    real(dp), allocatable :: rv(:)
    real(dp) :: expected(0:4)
    logical  :: degrees_ok
    integer  :: n

    call scalar_pair(1.0_dp, 2.0_dp, cells, inputs)

    call fill_path(full(1), quartic % argument(1), [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], cells)
    call fill_path(full(2), quartic % argument(2), [5.0_dp, 7.0_dp, 11.0_dp, 13.0_dp], cells)

    sparse(1) = full(1)
    call fill_path(sparse(2), quartic % argument(2), [5.0_dp], cells)

    expected = [31.0_dp, 271.0_dp, 2207.0_dp, 16688.0_dp, 118251.0_dp]

    degrees_ok = .true.
    do n = 0, 4
       call composer % assemble(quartic, lone, inputs, n, full, output)
       call output % real_vector(rv)
       degrees_ok = degrees_ok .and. size(rv) == 1 .and. &
            & near(rv(1), expected(n), 1.0e-10_dp)
    end do
    call report(degrees_ok, &
         & "quartic degrees 0..4: 31, 271, 2207, 16688, 118251", nfail)

    call composer % assemble(quartic, lone, inputs, 3, sparse, output)
    call output % real_vector(rv)
    call report(near(rv(1), 9156.0_dp, 1.0e-10_dp), &
         & "an unoccupied derivative reads as zero: sparse degree 3 is 9156", &
         & nfail)

    call composer % assemble(quartic, lone, inputs, 4, sparse, output)
    call output % real_vector(rv)
    call report(near(rv(1), 27908.0_dp, 1.0e-10_dp), &
         & "an unoccupied derivative reads as zero: sparse degree 4 is 27908", &
         & nfail)

  end subroutine check_the_chain_rule

  !===================================================================!
  ! Degrees 1..8 of Phi = (q + xi)^8 at (1, 2) with path
  ! derivatives q^(k) = k! and xi^(k) = k k!. Then u(e) = q(e) +
  ! xi(e) has the normalized Taylor coefficients a_0 = 3 and
  ! a_k = k + 1, all exact in floating point. The expected values
  ! are computed here by raising u(e) to the eighth power with
  ! truncated polynomial convolution, independently of the pattern
  ! generator; the degree-n total is n! times the coefficient of
  ! e^n.
  !===================================================================!

  subroutine check_the_taylor_convolution(nfail)

    integer, intent(inout) :: nfail

    type(chain_rule)    :: composer
    type(argument_path) :: paths(2)
    type(stored_field)         :: inputs(2)
    class(field), allocatable :: output
    real(dp), allocatable :: rv(:)
    real(dp) :: qseats(8), xiseats(8)   ! the path derivatives of q and xi
    real(dp) :: acoef(0:8), upow(0:8), convolved(0:8)
    real(dp) :: fact, expected
    logical  :: degrees_ok
    integer  :: i, j, k, n

    call scalar_pair(1.0_dp, 2.0_dp, cells, inputs)

    fact = 1.0_dp
    do k = 1, 8
       fact       = fact * real(k, dp)
       qseats(k)  = fact
       xiseats(k) = real(k, dp) * fact
    end do

    call fill_path(paths(1), p8 % argument(1), qseats, cells)
    call fill_path(paths(2), p8 % argument(2), xiseats, cells)

    acoef(0) = 3.0_dp
    do k = 1, 8
       acoef(k) = real(k + 1, dp)
    end do
    upow = acoef
    do k = 1, 7
       convolved = 0.0_dp
       do i = 0, 8
          do j = 0, 8 - i
             convolved(i + j) = convolved(i + j) + upow(i) * acoef(j)
          end do
       end do
       upow = convolved
    end do

    degrees_ok = .true.
    fact = 1.0_dp
    do n = 1, 8
       fact     = fact * real(n, dp)
       expected = fact * upow(n)
       call composer % assemble(p8, lone, inputs, n, paths, output)
       call output % real_vector(rv)
       degrees_ok = degrees_ok .and. size(rv) == 1 .and. &
            & abs(rv(1) - expected) <= 1.0e-12_dp * abs(expected)
    end do

    call report(degrees_ok, &
         & "every degree 1..8 of (q+xi)^8 matches the convolution oracle", &
         & nfail)

  end subroutine check_the_taylor_convolution

  !===================================================================!
  ! Marcher derivative tests on S = q^2 - xi with xi = 1, backward
  ! euler at h = 1, started at the fixed point q = 1: the
  ! trajectory is constant and every step Jacobian equals
  ! a0 + 2 h q = 3. With the parameter path xi^(1) = 1 the
  ! second-instant sensitivities satisfy
  !
  !      u_1 = 1/3,
  !      u_s = -( sum_{i=1}^{s-1} C(s,i) u_i u_(s-i) ) / 3
  !
  ! and the third-instant ones the same recursion plus the
  ! -a1 u_s term carried from the previous instant. The fractions
  ! below were computed by hand from these recursions.
  ! march_adjoint is backward substitution with the same Jacobian:
  ! a terminal seed of 1 becomes (1/3)^2 = 1/9 after two edges,
  ! and per-instant seeds of 1 give 1 + (1 + 1/3)/3 = 13/9. The
  ! last check marches S = q with steps (1, 1/2): backward euler
  ! gives q1 = 1/2, then q2 = q1 / (3/2) = 1/3.
  !===================================================================!

  subroutine check_the_derivative_walks(nfail)

    integer, intent(inout) :: nfail

    type(argument_path) :: xipath(1)
    type(stored_field) :: xifield(1)
    real(dp), allocatable :: trajectory(:,:), sensitivities(:,:,:)
    real(dp) :: q(1), lambda(1), seeds(1,3)
    real(dp) :: instant2(8), instant3(8)
    logical  :: orders_ok
    integer  :: s

    !----------------------------------------------------------------!
    ! march with trajectory recording: 3 instants, all equal to 1.
    !----------------------------------------------------------------!

    q = [1.0_dp]
    call clock % march(equil, lone, q, 2, trajectory=trajectory)

    call report(size(trajectory, 2) == 3 .and. &
         & maxval(abs(trajectory - 1.0_dp)) < 1.0e-12_dp, &
         & "march records the trajectory: constant at the fixed point", nfail)

    !----------------------------------------------------------------!
    ! march_directional to order 8, the parameter path supplied on
    ! input slot 2.
    !----------------------------------------------------------------!

    xifield(1) = stored_field('xi', cells, 1, num_components=1)
    call xifield(1) % set_real_vector([1.0_dp])
    call fill_path(xipath(1), equil % argument(2), [1.0_dp], cells)

    call clock % march_directional(equil, lone, 2, trajectory, 8, &
         & sensitivities, parameters=xifield, paths=xipath)

    instant2 = [ 1.0_dp / 3.0_dp,            -2.0_dp / 27.0_dp,       &
         &       4.0_dp / 81.0_dp,          -40.0_dp / 729.0_dp,      &
         &     560.0_dp / 6561.0_dp,      -1120.0_dp / 6561.0_dp,     &
         &   24640.0_dp / 59049.0_dp,   -640640.0_dp / 531441.0_dp ]

    instant3 = [ 4.0_dp / 9.0_dp,           -38.0_dp / 243.0_dp,      &
         &     340.0_dp / 2187.0_dp,     -14848.0_dp / 59049.0_dp,    &
         &  897680.0_dp / 1594323.0_dp, -95200.0_dp / 59049.0_dp,     &
         &  726772480.0_dp / 129140163.0_dp,                          &
         &  -80862642560.0_dp / 3486784401.0_dp ]

    orders_ok = .true.
    do s = 1, 8
       orders_ok = orders_ok .and. &
            & near(sensitivities(1, s, 1), 0.0_dp,       1.0e-14_dp) .and. &
            & near(sensitivities(1, s, 2), instant2(s),  1.0e-10_dp) .and. &
            & near(sensitivities(1, s, 3), instant3(s),  1.0e-10_dp)
    end do

    call report(orders_ok, &
         & "march_directional matches the exact rationals to degree 8", nfail)

    !----------------------------------------------------------------!
    ! march_adjoint over the recorded trajectory.
    !----------------------------------------------------------------!

    lambda = [1.0_dp]
    call clock % march_adjoint(equil, lone, lambda, 2, trajectory)
    call report(near(lambda(1), 1.0_dp / 9.0_dp, 1.0e-12_dp), &
         & "a terminal seed crosses two edges as (1/3)^2 = 1/9", nfail)

    lambda = [1.0_dp]
    seeds  = 1.0_dp
    call clock % march_adjoint(equil, lone, lambda, 2, trajectory, &
         & seeds=seeds)
    call report(near(lambda(1), 13.0_dp / 9.0_dp, 1.0e-12_dp), &
         & "per-instant seeds accumulate to 13/9", nfail)

    !----------------------------------------------------------------!
    ! Implicit march with per-edge steps (1, 1/2).
    !----------------------------------------------------------------!

    q = [1.0_dp]
    call clock % march(lin, lone, q, 2, steps=[1.0_dp, 0.5_dp])
    call report(near(q(1), 1.0_dp / 3.0_dp, 1.0e-9_dp), &
         & "S = q on steps (1, 1/2) gives exactly 1/3", nfail)

    !----------------------------------------------------------------!
    ! The difference road: S = q keeps the difference tangent, so
    ! the step Jacobian a0 + h = 2 at h = 1 is differenced from the
    ! whole residual; a terminal seed crosses two edges as (1/2)^2.
    !----------------------------------------------------------------!

    q = [1.0_dp]
    call clock % march(lin, lone, q, 2, trajectory=trajectory)
    lambda = [1.0_dp]
    call clock % march_adjoint(lin, lone, lambda, 2, trajectory)
    call report(near(lambda(1), 0.25_dp, 1.0e-7_dp), &
         & "on the difference road a terminal seed crosses two edges as (1/2)^2", &
         & nfail)

  end subroutine check_the_derivative_walks

  !===================================================================!
  ! march_adaptive with the halving policy, q0 = 2, duration 1/2.
  ! Each backward-euler trial solves h q^2 + q - (qold + h) = 0,
  ! whose positive root is
  !
  !      q = ( -1 + sqrt(1 + 4 h (qold + h)) ) / (2 h).
  !
  ! Tolerance 1: the first proposal h = 1/2 is accepted (error
  ! estimate |trial - q0| is about 0.5505), giving q = -1 + sqrt 6.
  ! Tolerance 0.45: h = 1/2 is rejected, the halved retry h = 1/4
  ! is accepted (estimate about 0.3944), and a second h = 1/4 edge
  ! completes the duration. Tolerance 1e-30 with max_attempts = 3:
  ! no trial is accepted, so completed must be false, steps_taken
  ! empty, and q unchanged, because a rejected trial never writes
  ! to q.
  !===================================================================!

  subroutine check_the_adaptive_walk(nfail)

    integer, intent(inout) :: nfail

    type(halving_policy) :: generous, strict, impossible
    real(dp), allocatable :: taken(:)
    real(dp) :: q(1), first, second
    logical  :: completed

    generous % first_step = 0.5_dp
    generous % tolerance  = 1.0_dp

    q = [2.0_dp]
    call clock % march_adaptive(equil, lone, q, 0.5_dp, generous, 5, &
         & taken, completed)
    call report(completed .and. size(taken) == 1 .and. &
         & taken(1) == 0.5_dp .and. &
         & near(q(1), -1.0_dp + sqrt(6.0_dp), 1.0e-9_dp), &
         & "a generous tolerance accepts the first trial: -1 + sqrt 6", nfail)

    strict % first_step = 0.5_dp
    strict % tolerance  = 0.45_dp

    first  = -2.0_dp + sqrt(13.0_dp)
    second = 2.0_dp * (sqrt(first + 1.25_dp) - 1.0_dp)

    q = [2.0_dp]
    call clock % march_adaptive(equil, lone, q, 0.5_dp, strict, 5, &
         & taken, completed)
    call report(completed .and. size(taken) == 2 .and. &
         & all(taken == 0.25_dp) .and. near(q(1), second, 1.0e-9_dp), &
         & "a rejected proposal halves: the march takes steps (1/4, 1/4)", &
         & nfail)

    impossible % first_step = 0.5_dp
    impossible % tolerance  = 1.0e-30_dp

    q = [2.0_dp]
    call clock % march_adaptive(equil, lone, q, 0.5_dp, impossible, 3, &
         & taken, completed)
    call report((.not. completed) .and. size(taken) == 0 .and. &
         & q(1) == 2.0_dp, &
         & "a spent attempt budget: completed false, no steps, q unchanged", &
         & nfail)

  end subroutine check_the_adaptive_walk

  pure function near(value, expected, tolerance) result(ok)

    real(dp), intent(in) :: value, expected, tolerance
    logical :: ok

    ok = abs(value - expected) <= tolerance

  end function near

  subroutine report(ok, label, nfail)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: label
    integer, intent(inout) :: nfail
    if (ok) then
       write(*,'(1x,a,a)') "PASS : ", label
    else
       write(*,'(1x,a,a)') "FAIL : ", label
       nfail = nfail + 1
    end if
  end subroutine report


  !===================================================================!
  ! The argument calculus. (1) The scheme's partials at theta = 1/2
  ! over the quartic at (q, xi) = (1, 2), qold = 0, h = 1/2. With
  ! Phi_q(1,2) = 26, Phi_xi(1,2) = 49, Phi_qq(1,2) = 32 and, at the
  ! previous state, Phi_q(0,2) = 8, Phi_xi(0,2) = 32, Phi_qq(0,2) = 8:
  !
  !      state, v = 3      a0 v + h theta Phi_q(1,2) v        = 22.5
  !      history(1), v = 3 a1 v + h (1-theta) Phi_q(0,2) v    = 3
  !      auxiliary, w = 1  h [theta Phi_xi(1,2) + (1-theta) Phi_xi(0,2)] = 20.25
  !      state x2          h theta Phi_qq(1,2)                = 8
  !      history(1) x2     h (1-theta) Phi_qq(0,2)            = 2
  !      state x history   0
  !      the residual      q - qold + h [theta Phi(1,2) + (1-theta) Phi(0,2)]
  !                        = 1 + (31 + 16)/4                  = 12.75
  !
  ! (2) A history state supplied as input 3 wins over the stored
  ! qold: with qold_in = 1 and theta = 1 the residual is
  ! 0 + h Phi(1,2) = 15.5, not 16.5. (3) The dual by basis under the
  ! Euclidean pairing, on S = q^2 - xi over three vertices at
  ! q = (1,2,3): in the state the block is diag(2q), so the dual of
  ! lambda is 2 q lambda and equals the compiled transpose; in xi the
  ! block is the column (-1,-1,-1), so the dual of lambda is the
  ! scalar -sum(lambda); both satisfy <J v, lambda> = <v, J^T lambda>.
  !===================================================================!

  subroutine check_the_argument_calculus(nfail)

    integer, intent(inout) :: nfail

    type(scheme) :: statement
    type(stored_directed_graph) :: three
    type(graph) :: points
    type(linearization) :: tangent_q, tangent_xi
    type(stencil) :: compiled, adjoint
    type(stored_field) :: inputs(2), with_history(3), v, w, one
    type(stored_field) :: qf, xif, vf, wf, lf
    class(field), allocatable :: output
    real(dp), allocatable :: rv(:), g(:), jv(:), gt(:)
    real(dp) :: q3(3), lambda(3), lhs, rhs
    logical  :: ok

    statement = backward_euler(quartic, 0.5_dp)
    statement % qold  = [0.0_dp]
    statement % theta = 0.5_dp
    call scalar_pair(1.0_dp, 2.0_dp, cells, inputs)

    v = stored_field('v', cells, 1, num_components=1)
    call v % set_real_vector([3.0_dp])
    w = stored_field('w', cells, 1, num_components=1)
    call w % set_real_vector([1.0_dp])
    one = stored_field('one', cells, 1, num_components=1)
    call one % set_real_vector([1.0_dp])

    call statement % partial_action(lone, inputs, [variation(statement % state(), v)], output)
    call output % real_vector(rv)
    ok = near(rv(1), 22.5_dp, 1.0e-12_dp)
    call statement % partial_action(lone, inputs, [variation(statement % history(1), v)], output)
    call output % real_vector(rv)
    ok = ok .and. near(rv(1), 3.0_dp, 1.0e-12_dp)
    call statement % partial_action(lone, inputs, [variation(statement % auxiliary(1), w)], output)
    call output % real_vector(rv)
    ok = ok .and. near(rv(1), 20.25_dp, 1.0e-12_dp)
    call statement % partial_action(lone, inputs, &
         & [variation(statement % state(), one), variation(statement % state(), one)], output)
    call output % real_vector(rv)
    ok = ok .and. near(rv(1), 8.0_dp, 1.0e-12_dp)
    call statement % partial_action(lone, inputs, &
         & [variation(statement % history(1), one), variation(statement % history(1), one)], output)
    call output % real_vector(rv)
    ok = ok .and. near(rv(1), 2.0_dp, 1.0e-12_dp)
    call statement % partial_action(lone, inputs, &
         & [variation(statement % state(), one), variation(statement % history(1), one)], output)
    call output % real_vector(rv)
    ok = ok .and. near(rv(1), 0.0_dp, 1.0e-12_dp)
    call statement % apply(lone, inputs, output)
    call output % real_vector(rv)
    ok = ok .and. near(rv(1), 12.75_dp, 1.0e-12_dp)
    call report(ok, &
         & "theta = 1/2: the scheme's partials in state, history and auxiliary &
         &are 22.5, 3, 20.25, 8, 2, 0 and the residual 12.75", nfail)

    statement % theta = 1.0_dp
    with_history(1:2) = inputs
    with_history(3) = stored_field('qold', cells, 1, num_components=1)
    call with_history(3) % set_real_vector([1.0_dp])
    call statement % apply(lone, with_history, output)
    call output % real_vector(rv)
    call report(near(rv(1), 15.5_dp, 1.0e-12_dp), &
         & "a history state supplied as input 3 is read before the stored qold: 15.5", &
         & nfail)

    three  = stored_directed_graph(3, tails=[integer ::], heads=[integer ::])
    points = three % vertex_set()
    q3     = [1.0_dp, 2.0_dp, 3.0_dp]
    lambda = [0.5_dp, -1.0_dp, 2.0_dp]

    qf = stored_field('q', points, 3, num_components=1)
    call qf % set_real_vector(q3)
    xif = stored_field('xi', cells, 1, num_components=1)
    call xif % set_real_vector([1.0_dp])
    vf = stored_field('v', points, 3, num_components=1)
    call vf % set_real_vector([1.0_dp, -2.0_dp, 0.5_dp])
    wf = stored_field('w', cells, 1, num_components=1)
    call wf % set_real_vector([0.75_dp])
    lf = stored_field('lambda', points, 3, num_components=1)
    call lf % set_real_vector(lambda)

    ! the state block, square: dual by basis against the compiled transpose
    tangent_q = tangent_of(equil, equil % argument(1), at_inputs=[qf, xif])
    call tangent_q % apply(three, [vf], output)
    call output % real_vector(jv)
    lhs = dot_product(jv, lambda)
    call dual_by_basis(tangent_q, three, lambda, g)
    rhs = dot_product([1.0_dp, -2.0_dp, 0.5_dp], g)
    compiled = stencil(tangent_q, three, 3)
    adjoint  = compiled % transpose()
    call adjoint % apply(three, [lf], output)
    call output % real_vector(gt)
    call report(near(lhs, rhs, 1.0e-12_dp) .and. maxval(abs(g - 2.0_dp * q3 * lambda)) < 1.0e-12_dp &
         & .and. maxval(abs(g - gt)) < 1.0e-12_dp, &
         & "state block: <J v, lambda> = <v, J^T lambda>, J^T lambda = 2 q lambda, &
         &and the dual by basis equals the compiled transpose", nfail)

    ! the auxiliary block, rectangular: three residuals, one parameter
    tangent_xi = tangent_of(equil, equil % argument(2), at_inputs=[qf, xif])
    call tangent_xi % apply(three, [wf], output)
    call output % real_vector(jv)
    lhs = dot_product(jv, lambda)
    call dual_by_basis(tangent_xi, three, lambda, g)
    rhs = 0.75_dp * g(1)
    call report(size(g) == 1 .and. near(lhs, rhs, 1.0e-12_dp) .and. &
         & near(g(1), -sum(lambda), 1.0e-12_dp), &
         & "auxiliary block: the dual of lambda in xi is -sum(lambda), and the &
         &pairing law holds on a rectangular block", nfail)

  end subroutine check_the_argument_calculus

end program test_graph_differentiation
