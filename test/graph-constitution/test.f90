!=====================================================================!
! Boundary-condition and constitutive-law consistency tests.
!
! A mesh of two cells and two named boundary faces:
!
!         in                e1                 out
!       o------ (1) -----------------> (2) ------o
!       a 0.7       a 2.0, d 0.8           a 1.1
!       d 0.3                              d 0.6
!
! Two properties are tested. FORMULA FIDELITY: the Robin
! condition computes, per tagged face, exactly the four numbers the
! old class_boundary_condition computes - dirichlet, neumann, and
! mixed, to machine precision. OPERATOR CONSISTENCY: for a > 0,
! the coefficients and stored value the condition supplies make the
! calculus reproduce the eliminated boundary flux on the row,
!
!      row with the boundary - row without it = lhs*q_p - rhs
!
! with lhs and rhs evaluated from the reference formulas, so both
! representations produce the same scalar.
!=====================================================================!

program test_graph_constitution

  use iso_fortran_env, only : dp => REAL64
  use view_directed, only : directed_graph
  use field_calculus, only : field
  use view_directed , only : SIDE_VERTEX, SIDE_EDGE
  use graph_fractal  , only : graph
  use map_set_store, only : set_store
  use field_stored  , only : stored_field
  use view_mesh   , only : mesh
  use operation_robin_condition, only : robin_condition, robin, dirichlet, neumann
  use operation_robin_condition, only : COEFFICIENT_LHS, COEFFICIENT_RHS
  use operation_robin_condition, only : COEFFICIENT_ADVECTION_LHS, COEFFICIENT_ADVECTION_RHS
  use operation_robin_condition, only : COEFFICIENT_OPERATOR
  use operation_differential, only : differential_operator
  use operation_balance  , only : balance
  use operation_conduction     , only : conduction, advection
  use operation_diffusion, only : diffusion_stencil
  use operation_stencil  , only : stencil
  use operation_gmres    , only : gmres

  implicit none

  integer :: num_failures

  num_failures = 0

  call check_tag_resolves_once(num_failures)
  call check_formula_fidelity(num_failures)
  call check_boundary_operator(num_failures)
  call check_conduction_law(num_failures)
  call check_advection_law(num_failures)
  call check_diffusion_solution(num_failures)
  call check_robin_affine_relation(num_failures)

  write(*, '(a)') ' ============================================='
  if (num_failures == 0) then
     write(*, '(a)') ' all constitution checks passed'
  else
     write(*, '(a, i0, a)') ' ', num_failures, ' constitution checks FAILED'
     error stop 1
  end if

contains

  subroutine report(satisfied, message, num_failures)

    logical         , intent(in)    :: satisfied
    character(len=*), intent(in)    :: message
    integer         , intent(inout) :: num_failures

    if (satisfied) then
       write(*, '(a)') ' PASS : ' // message
    else
       write(*, '(a)') ' FAIL : ' // message
       num_failures = num_failures + 1
    end if

  end subroutine report

  !===================================================================!
  ! The hand mesh, one place.
  !===================================================================!

  type(mesh) function two_cell_mesh() result(m)

    m = mesh(2, tails=[1, 1, 2], heads=[2, 0, 0], &
         & volumes      = [1.5_dp, 2.5_dp], &
         & cell_centres = [0.5_dp, 0.0_dp, 0.0_dp, &
         &                 1.5_dp, 0.0_dp, 0.0_dp], &
         & areas        = [2.0_dp, 0.7_dp, 1.1_dp], &
         & deltas       = [0.8_dp, 0.3_dp, 0.6_dp], &
         & normals      = [ 1.0_dp, 0.0_dp, 0.0_dp, &
         &                 -1.0_dp, 0.0_dp, 0.0_dp, &
         &                  1.0_dp, 0.0_dp, 0.0_dp], &
         & face_centres = [1.0_dp, 0.0_dp, 0.0_dp, &
         &                 0.0_dp, 0.0_dp, 0.0_dp, &
         &                 2.0_dp, 0.0_dp, 0.0_dp], &
         & weights      = [0.5_dp, 1.0_dp, 1.0_dp], &
         & etags        = [character(len=4) :: '', 'in', 'out'])

  end function two_cell_mesh

  !===================================================================!
  ! The tag becomes a member set, once.
  !===================================================================!

  subroutine check_tag_resolves_once(num_failures)

    integer, intent(inout) :: num_failures

    type(mesh) :: m
    type(robin_condition) :: bc

    !----------------------------------------------------------------!
    ! The tagged-face subset is constructed and evaluated inside this
    ! check, so their interpretation is local to it.
    !----------------------------------------------------------------!

    type(graph)     :: members
    type(set_store) :: sets

    m  = two_cell_mesh()
    bc = dirichlet('in', 5.0_dp)

    call bc % faces(m, sets, members)
    call report(sets % num_members_of(members) == 1, &
         & 'the tag names one face', num_failures)
    call report(sets % member_of(members, 1) == 2, &
         & 'and it is the second face', num_failures)

  end subroutine check_tag_resolves_once

  !===================================================================!
  ! FORMULA FIDELITY. Three conditions, four numbers each, against
  ! the old world's own functions on the same face.
  !===================================================================!

  subroutine check_formula_fidelity(num_failures)

    integer, intent(inout) :: num_failures

    type(mesh) :: m
    real(dp), parameter :: kappa = 1.7_dp, vn = 0.9_dp

    m = two_cell_mesh()

    call fidelity_case(m, dirichlet('in', 5.0_dp), &
         & 1.0_dp, 0.0_dp, 5.0_dp, 0.7_dp, 0.3_dp, &
         & kappa, vn, 'dirichlet', num_failures)

    call fidelity_case(m, neumann('out', 2.0_dp), &
         & 0.0_dp, 1.0_dp, 2.0_dp, 1.1_dp, 0.6_dp, &
         & kappa, vn, 'neumann', num_failures)

    call fidelity_case(m, robin('out', 2.0_dp, 3.0_dp, 4.0_dp), &
         & 2.0_dp, 3.0_dp, 4.0_dp, 1.1_dp, 0.6_dp, &
         & kappa, vn, 'mixed robin', num_failures)

  end subroutine check_formula_fidelity

  subroutine fidelity_case(m, bc, a, b, c, area, delta, kappa, vn, name, num_failures)

    type(mesh)            , intent(in)    :: m
    type(robin_condition) , intent(in)    :: bc
    real(dp)              , intent(in)    :: a, b, c, area, delta, kappa, vn
    character(len=*)      , intent(in)    :: name
    integer               , intent(inout) :: num_failures

    real(dp), allocatable :: computed_values(:)
    real(dp) :: denom
    real(dp), parameter :: tol = 1.0d-15

    ! The reference formulas for comparison
    ! now that the old module is gone: denom = a + b/delta.
    denom = a + b / delta

    call bc % coefficient_values(m, kappa, COEFFICIENT_LHS, computed_values)
    call report(size(computed_values) == 1 .and. &
         & abs(computed_values(1) - (-kappa * area * a / (delta * denom))) < tol, &
         & name // ': the diffusive diagonal matches the recorded formula', num_failures)

    call bc % coefficient_values(m, kappa, COEFFICIENT_RHS, computed_values)
    call report(abs(computed_values(1) - (-kappa * area * c / (delta * denom))) < tol, &
         & name // ': the diffusive constant matches', num_failures)

    call bc % coefficient_values(m, vn, COEFFICIENT_ADVECTION_LHS, computed_values)
    call report(abs(computed_values(1) - (-vn * area * (b / delta) / denom)) < tol, &
         & name // ': the advective diagonal matches', num_failures)

    call bc % coefficient_values(m, vn, COEFFICIENT_ADVECTION_RHS, computed_values)
    call report(abs(computed_values(1) - (vn * area * c / denom)) < tol, &
         & name // ': the advective constant matches', num_failures)

  end subroutine fidelity_case

  !===================================================================!
  ! OPERATOR CONSISTENCY. The condition's coefficients and stored value
  ! determine the boundary flux: the difference between
  ! the row with the condition and the row without it equals
  ! lhs*q_p - rhs - the old row's own boundary term, sign and all -
  ! both sides computed independently.
  !===================================================================!

  subroutine check_boundary_operator(num_failures)

    integer, intent(inout) :: num_failures

    type(mesh) :: m
    type(robin_condition) :: bc
    type(differential_operator) :: with_boundary, without
    type(stored_field) :: state
    type(graph) :: cells
    class(field), allocatable :: y
    real(dp), allocatable :: cin(:), bin(:), rows_with(:), rows_without(:)
    real(dp) :: c(3), b(3), expected
    real(dp), parameter :: kappa = 1.7_dp, kint = 1.3_dp
    real(dp), parameter :: q1 = 2.0_dp, q2 = 4.0_dp

    m   = two_cell_mesh()
    bc  = dirichlet('in', 5.0_dp)

    call bc % coefficient_values(m, kappa, COEFFICIENT_OPERATOR, cin)
    call bc % boundary_values(m, bin)

    ! Interior conductivity on e1; the condition on e2; nothing on e3.
    c = [kint * 2.0_dp, cin(1), 0.0_dp]
    b = [0.0_dp, bin(1), 0.0_dp]

    cells = m % vertex_set()
    state = stored_field('q', cells, m % num_vertices())
    call state % set_real_vector([q1, q2])

    with_boundary = differential_operator(SIDE_VERTEX, 2, coefficients=c, &
         & spacings=[0.8_dp, 0.3_dp, 0.6_dp], boundary_values=b)
    call with_boundary % apply(m, with_boundary % bind([state]), y)
    call y % real_vector(rows_with)

    without = differential_operator(SIDE_VERTEX, 2, &
         & coefficients=[c(1), 0.0_dp, 0.0_dp], &
         & spacings=[0.8_dp, 0.3_dp, 0.6_dp])
    call without % apply(m, without % bind([state]), y)
    call y % real_vector(rows_without)

    ! lhs*q1 - rhs, from the recorded formulas at a=1, b=0, c=5.
    expected = (-kappa * 0.7_dp / 0.3_dp) * q1 &
         &     - (-kappa * 0.7_dp * 5.0_dp / 0.3_dp)

    call report(abs((rows_with(1) - rows_without(1)) - expected) < 1.0d-11, &
         & 'the row difference is the eliminated boundary flux', num_failures)

    call report(abs(rows_with(2) - rows_without(2)) < 1.0d-14, &
         & 'and the far cell feels nothing from it', num_failures)

  end subroutine check_boundary_operator

  !===================================================================!
  ! THE CONDUCTION LAW. keff = n^T K n from the mesh's own normals:
  ! isotropic k evaluates to k on every unit normal; a diagonal tensor
  ! evaluates to the component its normal selects. The dictionary
  ! coefficient keff*area is defined on interior faces only.
  !===================================================================!

  subroutine check_conduction_law(num_failures)

    integer, intent(inout) :: num_failures

    type(mesh) :: m
    type(conduction) :: law
    real(dp), allocatable :: keff(:), c(:)
    real(dp) :: k(3, 3)

    m = two_cell_mesh()

    law = conduction(2.5_dp)
    call law % normal_value(m, keff)
    call report(all(abs(keff - 2.5_dp) < 1.0d-14), &
         & 'an isotropic material evaluates to k on every unit normal', num_failures)

    k = 0.0_dp
    k(1, 1) = 2.0_dp
    k(2, 2) = 3.0_dp
    k(3, 3) = 4.0_dp
    law = conduction(k)
    call law % normal_value(m, keff)
    call report(all(abs(keff - 2.0_dp) < 1.0d-14), &
         & 'a diagonal tensor evaluates to the component selected by its normal', num_failures)

    call law % edge_coefficients(m, .true., c)
    call report(abs(c(1) - 2.0_dp * 2.0_dp) < 1.0d-14, &
         & 'the interior coefficient is keff times the area', num_failures)
    call report(c(2) == 0.0_dp .and. c(3) == 0.0_dp, &
         & 'headless face coefficients are determined by boundary conditions', num_failures)

  end subroutine check_conduction_law

  !===================================================================!
  ! THE ADVECTION LAW. vn = v.n and the coefficient vn*area; the
  ! scheme is the calculus's one_sided flag, and both settings must
  ! reproduce the old assembler's weights,
  !
  !      upwind    wp = max(vn,0), wn = min(vn,0)
  !      central   wp = wn = vn/2
  !
  ! through the balance rows -A*(wp*q_p + wn*q_n), computed here by
  ! hand, both flow directions.
  !===================================================================!

  subroutine check_advection_law(num_failures)

    integer, intent(inout) :: num_failures

    type(mesh) :: m
    type(conduction) :: flow
    type(balance) :: sums
    type(stored_field) :: state
    type(graph) :: cells
    class(field), allocatable :: y
    real(dp), allocatable :: vn(:), c(:), computed_values(:)
    real(dp), parameter :: q1 = 2.0_dp, q2 = 4.0_dp
    real(dp) :: area, wp, wn

    m = two_cell_mesh()
    area = 2.0_dp

    cells = m % vertex_set()
    state = stored_field('q', cells, m % num_vertices())
    call state % set_real_vector([q1, q2])

    flow = advection([1.5_dp, 0.0_dp, 0.0_dp])
    call flow % normal_value(m, vn)
    call report(abs(vn(1) - 1.5_dp) < 1.0d-14, &
         & 'the normal speed is v dot n', num_failures)

    call flow % edge_coefficients(m, .true., c)
    call report(abs(c(1) - 1.5_dp * area) < 1.0d-14 .and. c(2) == 0.0_dp, &
         & 'the coefficient is vn times area, interior only', num_failures)

    ! Upwind, flow along the edge: the tail value is selected.
    sums = balance(edge_terms=[differential_operator(SIDE_EDGE, 0, &
         & coefficients=c, one_sided=.true.)])
    call sums % apply(m, sums % bind([state]), y)
    call y % real_vector(computed_values)

    wp = max(vn(1), 0.0_dp)
    wn = min(vn(1), 0.0_dp)
    call report(abs(computed_values(1) - (-area * (wp * q1 + wn * q2))) < 1.0d-12 .and. &
         &      abs(computed_values(2) - (+area * (wp * q1 + wn * q2))) < 1.0d-12, &
         & 'upwind with the flow: the old rows, both cells', num_failures)

    ! Upwind, flow against the edge: the head value is selected.
    flow = advection([-1.5_dp, 0.0_dp, 0.0_dp])
    call flow % edge_coefficients(m, .true., c)
    sums = balance(edge_terms=[differential_operator(SIDE_EDGE, 0, &
         & coefficients=c, one_sided=.true.)])
    call sums % apply(m, sums % bind([state]), y)
    call y % real_vector(computed_values)

    wp = max(-1.5_dp, 0.0_dp)
    wn = min(-1.5_dp, 0.0_dp)
    call report(abs(computed_values(1) - (-area * (wp * q1 + wn * q2))) < 1.0d-12, &
         & 'upwind against the flow: the head is upstream', num_failures)

    ! Central: both ends, evenly - the old half weights.
    flow = advection([1.5_dp, 0.0_dp, 0.0_dp])
    call flow % edge_coefficients(m, .true., c)
    sums = balance(edge_terms=[differential_operator(SIDE_EDGE, 0, &
         & coefficients=c, one_sided=.false.)])
    call sums % apply(m, sums % bind([state]), y)
    call y % real_vector(computed_values)

    wp = 0.5_dp * 1.5_dp
    wn = 0.5_dp * 1.5_dp
    call report(abs(computed_values(1) - (-area * (wp * q1 + wn * q2))) < 1.0d-12, &
         & 'central: the old half weights, exactly', num_failures)

  end subroutine check_advection_law

  !===================================================================!
  ! THE STATEMENT. The constitution's whole sentence in one call:
  ! a material, two boundaries, a mesh - and the compiled operator converges
  ! under GMRES, with a monotone solution between its prescribed values.
  !===================================================================!

  subroutine check_diffusion_solution(num_failures)

    integer, intent(inout) :: num_failures

    type(mesh) :: m
    type(stencil) :: op
    type(gmres) :: gm
    real(dp), allocatable :: g(:), rhs(:), x(:)
    real(dp) :: achieved

    m = two_cell_mesh()

    op = diffusion_stencil(m, conduction(1.7_dp), &
         & [dirichlet('in', 0.0_dp), dirichlet('out', 10.0_dp)])

    call gm % state(op, m, m % vertex_set(), m % num_vertices())
    gm % tolerance = 1.0d-12

    g = gm % affine
    rhs = -g

    allocate(x(2))
    x = 0.0_dp
    call gm % solve(rhs, x, achieved)

    call report(achieved < 1.0d-10, &
         & 'GMRES converges for the diffusion statement', num_failures)
    call report(x(1) > 0.0_dp .and. x(1) < x(2) .and. x(2) < 10.0_dp, &
         & 'the solution is monotone between the prescribed boundary values', num_failures)

  end subroutine check_diffusion_solution

  !===================================================================!
  ! Prescribed boundary values and prescribed boundary gradients define
  ! distinct conditions. Eliminating the face gives an affine relation,
  !
  !      phi_b = (1 - w)*phi_p + v
  !
  ! A statement that represents only v omits the coefficient of phi_p:
  ! it would reduce an insulated boundary to a zero Dirichlet boundary
  ! and a Robin boundary to a Dirichlet boundary c/a. Two exact solutions test
  ! the difference.
  !===================================================================!

  subroutine check_robin_affine_relation(num_failures)

    integer, intent(inout) :: num_failures

    type(mesh) :: m
    type(stencil) :: op
    real(dp), allocatable :: x(:), mixed(:), dirichlet_solution(:)
    real(dp) :: achieved

    m = two_cell_mesh()

    ! A prescribed value of five and an insulated opposite boundary give
    ! zero flux and a constant solution of five.
    op = diffusion_stencil(m, conduction(1.7_dp), &
         & [dirichlet('in', 5.0_dp), neumann('out', 0.0_dp)])

    call solve_statement(op, m, x, achieved)

    call report(achieved < 1.0d-10, &
         & 'the insulated statement closes', num_failures)
    call report(all(abs(x - 5.0_dp) < 1.0d-9), &
         & 'an insulated boundary gives the prescribed constant solution', num_failures)

    ! The Robin condition robin(2, 3, 5) has the same c/a
    ! with dirichlet(2.5) and must not compile to the same operator.
    op = diffusion_stencil(m, conduction(1.7_dp), &
         & [dirichlet('in', 0.0_dp), robin('out', 2.0_dp, 3.0_dp, 5.0_dp)])
    call solve_statement(op, m, mixed, achieved)

    op = diffusion_stencil(m, conduction(1.7_dp), &
         & [dirichlet('in', 0.0_dp), dirichlet('out', 2.5_dp)])
    call solve_statement(op, m, dirichlet_solution, achieved)

    call report(maxval(abs(mixed - dirichlet_solution)) > 1.0d-6, &
         & 'a Robin boundary differs from a Dirichlet boundary with value c/a', num_failures)

    ! The Robin flux is smaller than the Dirichlet flux at c/a,
    ! so the Robin solution has a smaller boundary-cell value.
    call report(mixed(2) < dirichlet_solution(2) .and. mixed(2) > 0.0_dp, &
         & 'the Robin solution lies between zero and the Dirichlet solution', num_failures)

  end subroutine check_robin_affine_relation

  !===================================================================!
  ! Solve one residual equation with GMRES.
  !===================================================================!

  subroutine solve_statement(op, m, x, achieved)

    type(stencil), intent(in) :: op
    type(mesh), intent(in)             :: m
    real(dp), allocatable, intent(out) :: x(:)
    real(dp), intent(out)              :: achieved

    type(gmres) :: gm
    real(dp), allocatable :: g(:), rhs(:)

    call gm % state(op, m, m % vertex_set(), m % num_vertices())
    gm % tolerance      = 1.0d-12
    gm % max_iterations = 200

    g = gm % affine
    rhs = -g

    allocate(x(m % num_vertices()))
    x = 0.0_dp
    call gm % solve(rhs, x, achieved)

  end subroutine solve_statement

end program test_graph_constitution
