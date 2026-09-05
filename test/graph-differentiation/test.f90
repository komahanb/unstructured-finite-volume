!=====================================================================!
! Tests for the differentiation stack:
!
!      slope_at_zero      the nonuniform and uniform order-2 BDF
!                         coefficients and the order-1 coefficients,
!                         read from the family's Lagrange slopes
!      tangent_of         exact/difference dispatch, and the value
!                         of the exact tangent
!      total_derivative   total derivatives to degree 8
!      the pairing law    <J v, lambda> = <v, J^T lambda> with J^T
!                         the compiled transpose of the tangent
!      stored inputs      a minimizer's fixed inputs enter the
!                         residual and the frozen tangent alike
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
  use operation_chain_rule, only : total_derivative, derivative_of, argument_path
  use operation_linearization, only : linearization, tangent_of
  use operation_stencil , only : stencil
  use operation_family  , only : slope_at_zero
  use util_derivative_terms, only : derivative_terms, value
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

  call check_the_bdf_coefficients(nfail)
  call check_the_tangent_chooser(nfail)
  call check_the_chain_rule(nfail)
  call check_the_taylor_convolution(nfail)
  call check_the_pairing_law(nfail)
  call check_stored_inputs(nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all differentiation tower checks passed"
  else
     error stop
  end if

contains

  !===================================================================!
  ! The BDF coefficients are the slopes at zero of the Lagrange basis
  ! through the nodes u_j = -(t_k - t_(k-j)) / dt_k. Order 2 at the
  ! steps (2, 3), the newest step first: the nodes are (0, -1, -5/2)
  ! and the slopes are a0 = (2 h0 + h1)/(h0 + h1) = 7/5,
  ! a1 = -(h0 + h1)/h1 = -5/3, a2 = h0^2/(h1 (h0 + h1)) = 4/15.
  ! Equal steps, nodes (0, -1, -2), give (3/2, -2, 1/2) exactly,
  ! because every quotient formed is a dyadic rational. Order 1,
  ! nodes (0, -1), gives (1, -1) for any step.
  !===================================================================!

  subroutine check_the_bdf_coefficients(nfail)

    integer, intent(inout) :: nfail

    type(derivative_terms) :: nodes(0:2)
    real(dp) :: a(0:2)
    integer  :: j

    nodes(0) = derivative_terms(0.0_dp, 0)
    nodes(1) = derivative_terms(-1.0_dp, 0)
    nodes(2) = derivative_terms(-2.5_dp, 0)
    do j = 0, 2
       a(j) = value(slope_at_zero(nodes, j))
    end do
    call report( &
         & near(a(0),  7.0_dp /  5.0_dp, 1.0e-15_dp) .and. &
         & near(a(1), -5.0_dp /  3.0_dp, 1.0e-15_dp) .and. &
         & near(a(2),  4.0_dp / 15.0_dp, 1.0e-15_dp), &
         & "bdf 2 at steps (2, 3) gives [7/5, -5/3, 4/15]", nfail)

    nodes(2) = derivative_terms(-2.0_dp, 0)
    do j = 0, 2
       a(j) = value(slope_at_zero(nodes, j))
    end do
    call report( &
         & a(0) ==  1.5_dp .and. &
         & a(1) == -2.0_dp .and. &
         & a(2) ==  0.5_dp, &
         & "equal steps give the uniform coefficients exactly", nfail)

    do j = 0, 1
       a(j) = value(slope_at_zero(nodes(0:1), j))
    end do
    call report( &
         & a(0) ==  1.0_dp .and. &
         & a(1) == -1.0_dp, &
         & "bdf 1 gives [1, -1]", nfail)

  end subroutine check_the_bdf_coefficients

  !===================================================================!
  ! tangent_of must take the exact mode when the operation's
  ! max_degree is at least one and the difference mode otherwise;
  ! the two are distinguished here by their name() prefixes. The
  ! exact tangent of the quartic frozen at
  ! q = 1 (xi defaulting to 2) is Phi_q(1, 2) = 26, so applying it
  ! to the direction v = 3 must return 78.
  !===================================================================!

  subroutine check_the_tangent_chooser(nfail)

    integer, intent(inout) :: nfail

    type(linearization) :: tangent, slow
    class(field), allocatable :: output
    type(stored_field) :: direction, state
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

    state = stored_field('q', cells, 1, num_components=1)
    call state % set_real_vector([1.0_dp])
    call tangent % freeze([state])

    direction = stored_field('v', cells, 1, num_components=1)
    call direction % set_real_vector([3.0_dp])
    call tangent % apply(lone, tangent % bind([direction]), output)
    call output % real_vector(rv)

    call report(size(rv) == 1 .and. near(rv(1), 78.0_dp, 1.0e-12_dp), &
         & "the exact tangent of the quartic at 1: J v = 26 v", nfail)

  end subroutine check_the_tangent_chooser

  !===================================================================!
  ! total_derivative % assemble on the quartic at (q, xi) = (1, 2)
  ! with path derivatives q^(k) = (1, 2, 3, 4) and xi^(k) = (5, 7,
  ! 11, 13). Expected totals were computed by rational Taylor
  ! composition, independently of the pattern generator:
  !
  !      d0 = 31,  d1 = 271,  d2 = 2207,  d3 = 16688,  d4 = 118251
  !
  ! With xi's derivatives past the first unoccupied, those terms
  ! are read as zero: d3 = 9156, d4 = 27908. The sparse degrees are
  ! taken through derivative_of, the composed rule applied to the
  ! statement's own inputs.
  !===================================================================!

  subroutine check_the_chain_rule(nfail)

    integer, intent(inout) :: nfail

    type(total_derivative) :: composer, composed
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

    composer = total_derivative(4)
    degrees_ok = .true.
    do n = 0, 4
       call composer % assemble(quartic, lone, quartic % bind(inputs), n, full, output)
       call output % real_vector(rv)
       degrees_ok = degrees_ok .and. size(rv) == 1 .and. &
            & near(rv(1), expected(n), 1.0e-10_dp)
    end do
    call report(degrees_ok, &
         & "quartic degrees 0..4: 31, 271, 2207, 16688, 118251", nfail)

    composed = derivative_of(quartic, 3, sparse)
    call composed % apply(lone, composed % bind(inputs), output)
    call output % real_vector(rv)
    call report(near(rv(1), 9156.0_dp, 1.0e-10_dp), &
         & "an unoccupied derivative reads as zero: sparse degree 3 is 9156", &
         & nfail)

    composed = derivative_of(quartic, 4, sparse)
    call composed % apply(lone, composed % bind(inputs), output)
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

    type(total_derivative) :: composer
    type(argument_path) :: paths(2)
    type(stored_field)         :: inputs(2)
    class(field), allocatable :: output
    real(dp), allocatable :: rv(:)
    real(dp) :: qpath(8), xipath(8)   ! the path derivatives of q and xi
    real(dp) :: acoef(0:8), upow(0:8), convolved(0:8)
    real(dp) :: fact, expected
    logical  :: degrees_ok
    integer  :: i, j, k, n

    call scalar_pair(1.0_dp, 2.0_dp, cells, inputs)

    fact = 1.0_dp
    do k = 1, 8
       fact      = fact * real(k, dp)
       qpath(k)  = fact
       xipath(k) = real(k, dp) * fact
    end do

    call fill_path(paths(1), p8 % argument(1), qpath, cells)
    call fill_path(paths(2), p8 % argument(2), xipath, cells)

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

    composer = total_derivative(8)
    degrees_ok = .true.
    fact = 1.0_dp
    do n = 1, 8
       fact     = fact * real(n, dp)
       expected = fact * upow(n)
       call composer % assemble(p8, lone, p8 % bind(inputs), n, paths, output)
       call output % real_vector(rv)
       degrees_ok = degrees_ok .and. size(rv) == 1 .and. &
            & abs(rv(1) - expected) <= 1.0e-12_dp * abs(expected)
    end do

    call report(degrees_ok, &
         & "every degree 1..8 of (q+xi)^8 matches the convolution oracle", &
         & nfail)

  end subroutine check_the_taylor_convolution

  pure function near(value, expected, tolerance) result(passed)

    real(dp), intent(in) :: value, expected, tolerance
    logical :: passed

    passed = abs(value - expected) <= tolerance

  end function near

  subroutine report(passed, label, nfail)
    logical, intent(in) :: passed
    character(len=*), intent(in) :: label
    integer, intent(inout) :: nfail
    if (passed) then
       write(*,'(1x,a,a)') "PASS : ", label
    else
       write(*,'(1x,a,a)') "FAIL : ", label
       nfail = nfail + 1
    end if
  end subroutine report

  !===================================================================!
  ! The pairing law under the Euclidean pairing, on S = q^2 - xi
  ! over three vertices at q = (1,2,3). In the state the block is
  ! diag(2q), so J^T lambda is 2 q lambda and the compiled transpose
  ! of the exact tangent returns it, with <J v, lambda> = <v, J^T
  ! lambda>. In xi the block is the column (-1,-1,-1), so the tangent
  ! along w = 3/4 is -3/4 at every vertex and <J w, lambda> is
  ! -(3/4) sum(lambda): the dual of lambda in xi is -sum(lambda).
  !===================================================================!

  subroutine check_the_pairing_law(nfail)

    integer, intent(inout) :: nfail

    type(stored_directed_graph) :: three
    type(graph) :: points
    type(linearization) :: tangent_q, tangent_xi
    type(stencil) :: compiled, adjoint
    type(stored_field) :: qf, xif, vf, wf, lf
    class(field), allocatable :: output
    real(dp), allocatable :: jv(:), gt(:)
    real(dp) :: q3(3), lambda(3), v3(3), lhs, rhs

    three  = stored_directed_graph(3, tails=[integer ::], heads=[integer ::])
    points = three % vertex_set()
    q3     = [1.0_dp, 2.0_dp, 3.0_dp]
    lambda = [0.5_dp, -1.0_dp, 2.0_dp]
    v3     = [1.0_dp, -2.0_dp, 0.5_dp]

    qf = stored_field('q', points, 3, num_components=1)
    call qf % set_real_vector(q3)
    xif = stored_field('xi', cells, 1, num_components=1)
    call xif % set_real_vector([1.0_dp])
    vf = stored_field('v', points, 3, num_components=1)
    call vf % set_real_vector(v3)
    wf = stored_field('w', cells, 1, num_components=1)
    call wf % set_real_vector([0.75_dp])
    lf = stored_field('lambda', points, 3, num_components=1)
    call lf % set_real_vector(lambda)

    ! the state block, square: the compiled transpose of the tangent
    tangent_q = tangent_of(equil, equil % argument(1))
    call tangent_q % freeze([qf, xif])
    call tangent_q % apply(three, tangent_q % bind([vf]), output)
    call output % real_vector(jv)
    lhs = dot_product(jv, lambda)
    compiled = stencil(tangent_q, three, 3)
    adjoint  = compiled % transpose()
    call adjoint % apply(three, adjoint % bind([lf]), output)
    call output % real_vector(gt)
    rhs = dot_product(v3, gt)
    call report(near(lhs, rhs, 1.0e-12_dp) .and. &
         & maxval(abs(gt - 2.0_dp * q3 * lambda)) < 1.0e-12_dp, &
         & "state block: <J v, lambda> = <v, J^T lambda>, and the compiled &
         &transpose returns J^T lambda = 2 q lambda", nfail)

    ! the auxiliary block, rectangular: three residuals, one parameter
    tangent_xi = tangent_of(equil, equil % argument(2))
    call tangent_xi % freeze([qf, xif])
    call tangent_xi % apply(three, tangent_xi % bind([wf]), output)
    call output % real_vector(jv)
    lhs = dot_product(jv, lambda)
    rhs = 0.75_dp * (-sum(lambda))
    call report(size(jv) == 3 .and. maxval(abs(jv + 0.75_dp)) < 1.0e-12_dp .and. &
         & near(lhs, rhs, 1.0e-12_dp), &
         & "auxiliary block: the dual of lambda in xi is -sum(lambda), and the &
         &pairing law is satisfied on a rectangular block", nfail)

  end subroutine check_the_pairing_law

  !===================================================================!
  ! Stored inputs. Newton stated on the quartic with xi = 3 stored:
  ! the residual is evaluated on [q, xi] and the Jacobian is frozen
  ! at the same tuple, so the tangent at q = 1 is Phi_q(1, 3) = 58
  ! (it would read 26 at the default xi = 2), it agrees with the
  ! differenced matvec, and the solve of Phi(q, 3) = Phi(2, 3) = 211
  ! from q = 1 lands on 2.
  !===================================================================!

  subroutine check_stored_inputs(nfail)

    integer, intent(inout) :: nfail

    type(newton) :: solver
    type(linearization) :: jacobian
    type(stored_field) :: stored(1), v
    type(stored_field), allocatable :: inputs(:)
    class(field), allocatable :: output
    real(dp), allocatable :: rv(:), y0(:), y1(:)
    real(dp) :: x(1), achieved, differenced
    real(dp), parameter :: eps = 1.0e-7_dp

    stored(1) = stored_field('xi', cells, 1, num_components=1)
    call stored(1) % set_real_vector([3.0_dp])
    v = stored_field('v', cells, 1, num_components=1)
    call v % set_real_vector([1.0_dp])

    allocate(solver % inner, source=gmres())
    solver % inner % tolerance = 1.0d-14
    solver % tolerance = 1.0d-12
    call solver % state(quartic, lone, cells, 1, stored_inputs=stored)

    call solver % evaluate([1.0_dp], y0, inputs)
    jacobian = tangent_of(quartic, quartic % argument(1))
    call jacobian % freeze(inputs)
    call jacobian % apply(lone, jacobian % bind([v]), output)
    call output % real_vector(rv)

    call solver % matvec([1.0_dp], y0)
    call solver % matvec([1.0_dp + eps], y1)
    differenced = (y1(1) - y0(1)) / eps

    call report(size(inputs) == 2 .and. near(rv(1), 58.0_dp, 1.0e-12_dp) .and. &
         & near(differenced, 58.0_dp, 1.0e-4_dp), &
         & "the Jacobian frozen at the evaluation tuple reads the stored xi: &
         &Phi_q(1, 3) = 58, and agrees with the differenced residual", nfail)

    x = [1.0_dp]
    call solver % solve([211.0_dp], x, achieved)
    call report(near(x(1), 2.0_dp, 1.0e-9_dp), &
         & "newton solves Phi(q, 3) = 211 with xi stored: q = 2", nfail)

  end subroutine check_stored_inputs

end program test_graph_differentiation
