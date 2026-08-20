!=====================================================================!
! Tests for the differentiation stack:
!
!      set_bdf            the nonuniform and uniform order-2 rows
!                         and the order-1 row
!      tangent_of         exact/difference dispatch, and the value
!                         of the exact tangent
!      chain_rule         total derivatives to degree 8
!      the marcher        march with trajectory recording,
!                         march_directional to degree 8,
!                         march_adjoint with terminal and
!                         per-instant seeds, march_adaptive, and a
!                         nonuniform implicit march
!
! Every expected value is an exact rational or a closed form
! derived independently of the code under test; the derivation is
! stated in each section comment.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_graph_differentiation

  use iso_fortran_env     , only : dp => REAL64
  use graph_field_calculus, only : graph_field
  use graph_discretization      , only : linearization_operator
  use fractal_graph       , only : set_graph => graph
  use class_graph         , only : directed_stored_graph
  use class_graph_field   , only : field
  use class_graph_step    , only : step_operator, backward_euler, bdf_variable
  use class_graph_chain_rule, only : chain_rule, argument_path
  use class_graph_linearization, only : tangent_of
  use class_graph_marcher , only : marcher, MARCH_BACKWARD
  use class_graph_step_policy, only : halving_policy
  use class_graph_newton  , only : newton
  use class_graph_gmres   , only : gmres
  use toy_differentiable_forms, only : quartic_form, power8_form, &
       & equilibrium_law, linear_law, scalar_pair, fill_path

  implicit none

  type(quartic_form)    :: quartic
  type(power8_form)     :: p8
  type(equilibrium_law) :: equil
  type(linear_law)      :: lin

  type(directed_stored_graph) :: lone
  type(set_graph)             :: cells
  type(marcher)               :: clock

  integer :: nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "graph differentiation tower suite"
  write(*,'(1x,a)') "============================================="

  lone  = directed_stored_graph(1, tails=[integer ::], heads=[integer ::])
  cells = lone % vertex_set()

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
  call check_the_chain_rule(nfail)
  call check_the_taylor_convolution(nfail)
  call check_the_derivative_walks(nfail)
  call check_the_adaptive_walk(nfail)

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

    type(step_operator) :: statement

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
  ! tangent_of must return exact_linearization for a
  ! differentiable_operation and difference_linearization for any
  ! other operation; the two are distinguished here by their
  ! name() prefixes. The exact tangent of the quartic frozen at
  ! q = 1 (xi defaulting to 2) is Phi_q(1, 2) = 26, so applying it
  ! to the direction v = 3 must return 78.
  !===================================================================!

  subroutine check_the_tangent_chooser(nfail)

    integer, intent(inout) :: nfail

    class(linearization_operator), allocatable :: tangent, slow
    class(graph_field), allocatable :: output
    type(field) :: direction
    real(dp), allocatable :: rv(:)

    tangent = tangent_of(quartic)
    slow    = tangent_of(lin)

    call report(index(tangent % name(), 'exact derivative of') == 1, &
         & "tangent_of picks the exact linearization when differentiable", &
         & nfail)
    call report(index(slow % name(), 'derivative of') == 1 .and. &
         & index(slow % name(), 'exact') == 0, &
         & "tangent_of picks the difference linearization otherwise", nfail)

    call tangent % freeze([1.0_dp])

    direction = field('v', cells, 1, ncomp=1)
    call direction % set_real_vector([3.0_dp])
    call tangent % apply(lone, [direction], output)
    call output % get_real_vector(rv)

    call report(size(rv) == 1 .and. near(rv(1), 78.0_dp, 1.0e-12_dp), &
         & "the exact tangent of the quartic at 1: J v = 26 v", nfail)

  end subroutine check_the_tangent_chooser

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
    type(field)         :: inputs(2)
    class(graph_field), allocatable :: output
    real(dp), allocatable :: rv(:)
    real(dp) :: expected(0:4)
    logical  :: degrees_ok
    integer  :: n

    call scalar_pair(1.0_dp, 2.0_dp, cells, inputs)

    call fill_path(full(1), 1, [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], cells)
    call fill_path(full(2), 2, [5.0_dp, 7.0_dp, 11.0_dp, 13.0_dp], cells)

    sparse(1) = full(1)
    call fill_path(sparse(2), 2, [5.0_dp], cells)

    expected = [31.0_dp, 271.0_dp, 2207.0_dp, 16688.0_dp, 118251.0_dp]

    degrees_ok = .true.
    do n = 0, 4
       call composer % assemble(quartic, lone, inputs, n, full, output)
       call output % get_real_vector(rv)
       degrees_ok = degrees_ok .and. size(rv) == 1 .and. &
            & near(rv(1), expected(n), 1.0e-10_dp)
    end do
    call report(degrees_ok, &
         & "quartic degrees 0..4: 31, 271, 2207, 16688, 118251", nfail)

    call composer % assemble(quartic, lone, inputs, 3, sparse, output)
    call output % get_real_vector(rv)
    call report(near(rv(1), 9156.0_dp, 1.0e-10_dp), &
         & "an unoccupied derivative reads as zero: sparse degree 3 is 9156", &
         & nfail)

    call composer % assemble(quartic, lone, inputs, 4, sparse, output)
    call output % get_real_vector(rv)
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
    type(field)         :: inputs(2)
    class(graph_field), allocatable :: output
    real(dp), allocatable :: rv(:)
    real(dp) :: qseats(8), xiseats(8)
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

    call fill_path(paths(1), 1, qseats, cells)
    call fill_path(paths(2), 2, xiseats, cells)

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
       call output % get_real_vector(rv)
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
    type(field) :: xifield(1)
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

    xifield(1) = field('xi', cells, 1, ncomp=1)
    call xifield(1) % set_real_vector([1.0_dp])
    call fill_path(xipath(1), 2, [1.0_dp], cells)

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
    call clock % march_adjoint(equil, lone, lambda, 2, action=equil, &
         & trajectory=trajectory)
    call report(near(lambda(1), 1.0_dp / 9.0_dp, 1.0e-12_dp), &
         & "a terminal seed crosses two edges as (1/3)^2 = 1/9", nfail)

    lambda = [1.0_dp]
    seeds  = 1.0_dp
    call clock % march_adjoint(equil, lone, lambda, 2, action=equil, &
         & trajectory=trajectory, seeds=seeds)
    call report(near(lambda(1), 13.0_dp / 9.0_dp, 1.0e-12_dp), &
         & "per-instant seeds accumulate to 13/9", nfail)

    !----------------------------------------------------------------!
    ! Implicit march with per-edge steps (1, 1/2).
    !----------------------------------------------------------------!

    q = [1.0_dp]
    call clock % march(lin, lone, q, 2, steps=[1.0_dp, 0.5_dp])
    call report(near(q(1), 1.0_dp / 3.0_dp, 1.0e-9_dp), &
         & "S = q on steps (1, 1/2) gives exactly 1/3", nfail)

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

end program test_graph_differentiation
