!=====================================================================!
! The constitution suite: rung 3's acceptance.
!
! A hand mesh of two cells and two named boundary faces:
!
!         in                e1                 out
!       o------ (1) -----------------> (2) ------o
!       a 0.7       a 2.0, d 0.8           a 1.1
!       d 0.3                              d 0.6
!
! Two things are proven here. FORMULA FIDELITY: the new robin
! condition computes, per tagged face, exactly the four numbers the
! old class_boundary_condition computes - dirichlet, neumann, and
! mixed, to machine precision. And THE OPERATOR ROAD: for a > 0,
! the coefficients and stored value the condition supplies make the
! calculus reproduce the eliminated boundary flux on the row,
!
!      row with the wall  -  row without it  =  lhs*q_p - rhs
!
! with lhs and rhs taken from the OLD functions, so the two worlds
! meet on one number.
!=====================================================================!

program test_graph_constitution

  use iso_fortran_env, only : dp => REAL64
  use graph_grammar  , only : graph, graph_field
  use graph_calculus , only : GRAPH_SIDE_VERTEX
  use class_graph_support, only : support
  use class_graph_field  , only : field
  use class_graph_mesh   , only : mesh
  use class_robin_condition, only : robin_condition, robin, dirichlet, neumann
  use class_boundary_condition, only : old_condition => boundary_condition
  use class_boundary_condition, only : old_robin => robin
  use class_graph_differential_operator, only : differential_operator
  use class_graph_differential_operator, only : vertex_differential_operator
  use class_graph_differential_operator, only : edge_differential_operator
  use class_graph_balance  , only : balance
  use class_conduction     , only : conduction
  use class_advection      , only : advection
  use class_diffusion_statement, only : diffusion_statement
  use class_graph_stencil  , only : stencil_operator
  use class_graph_gmres    , only : gmres

  implicit none

  integer :: nfail

  nfail = 0

  call check_tag_resolves_once(nfail)
  call check_formula_fidelity(nfail)
  call check_operator_road(nfail)
  call check_conduction_law(nfail)
  call check_advection_law(nfail)
  call check_the_statement_speaks(nfail)

  write(*, '(a)') ' ============================================='
  if (nfail == 0) then
     write(*, '(a)') ' all constitution checks passed'
  else
     write(*, '(a, i0, a)') ' ', nfail, ' constitution checks FAILED'
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
  ! The hand mesh, one place.
  !===================================================================!

  type(mesh) function hand_mesh() result(m)

    m = mesh(2, tails=[1, 1, 2], heads=[2, 0, 0], &
         & volumes      = [1.5_dp, 2.5_dp], &
         & cell_centers = [0.5_dp, 0.0_dp, 0.0_dp, &
         &                 1.5_dp, 0.0_dp, 0.0_dp], &
         & areas        = [2.0_dp, 0.7_dp, 1.1_dp], &
         & deltas       = [0.8_dp, 0.3_dp, 0.6_dp], &
         & normals      = [ 1.0_dp, 0.0_dp, 0.0_dp, &
         &                 -1.0_dp, 0.0_dp, 0.0_dp, &
         &                  1.0_dp, 0.0_dp, 0.0_dp], &
         & face_centers = [1.0_dp, 0.0_dp, 0.0_dp, &
         &                 0.0_dp, 0.0_dp, 0.0_dp, &
         &                 2.0_dp, 0.0_dp, 0.0_dp], &
         & weights      = [0.5_dp, 1.0_dp, 1.0_dp], &
         & etags        = [character(len=4) :: '', 'in', 'out'])

  end function hand_mesh

  !===================================================================!
  ! The tag becomes a member set, once.
  !===================================================================!

  subroutine check_tag_resolves_once(nfail)

    integer, intent(inout) :: nfail

    type(mesh) :: m
    type(robin_condition) :: bc
    class(graph), allocatable :: members

    m  = hand_mesh()
    bc = dirichlet('in', 5.0_dp)

    call bc % faces(m, members)
    call report(members % num_vertices() == 1, &
         & 'the tag names one face', nfail)
    call report(members % global_vertex_index(1) == 2, &
         & 'and it is the second face', nfail)

  end subroutine check_tag_resolves_once

  !===================================================================!
  ! FORMULA FIDELITY. Three conditions, four numbers each, against
  ! the old world's own functions on the same face.
  !===================================================================!

  subroutine check_formula_fidelity(nfail)

    integer, intent(inout) :: nfail

    type(mesh) :: m
    real(dp), parameter :: kappa = 1.7_dp, vn = 0.9_dp

    m = hand_mesh()

    call fidelity_case(m, dirichlet('in', 5.0_dp), &
         & old_robin(1.0_dp, 0.0_dp, 5.0_dp), 0.7_dp, 0.3_dp, &
         & kappa, vn, 'dirichlet', nfail)

    call fidelity_case(m, neumann('out', 2.0_dp), &
         & old_robin(0.0_dp, 1.0_dp, 2.0_dp), 1.1_dp, 0.6_dp, &
         & kappa, vn, 'neumann', nfail)

    call fidelity_case(m, robin('out', 2.0_dp, 3.0_dp, 4.0_dp), &
         & old_robin(2.0_dp, 3.0_dp, 4.0_dp), 1.1_dp, 0.6_dp, &
         & kappa, vn, 'mixed robin', nfail)

  end subroutine check_formula_fidelity

  subroutine fidelity_case(m, bc, old, area, delta, kappa, vn, name, nfail)

    type(mesh)            , intent(in)    :: m
    type(robin_condition) , intent(in)    :: bc
    type(old_condition)   , intent(in)    :: old
    real(dp)              , intent(in)    :: area, delta, kappa, vn
    character(len=*)      , intent(in)    :: name
    integer               , intent(inout) :: nfail

    real(dp), allocatable :: got(:)
    real(dp), parameter :: tol = 1.0d-15

    call bc % lhs_coefficients(m, kappa, got)
    call report(size(got) == 1 .and. &
         & abs(got(1) - old % lhs_coeff(area, delta, kappa)) < tol, &
         & name // ': the diffusive diagonal matches the old world', nfail)

    call bc % rhs_coefficients(m, kappa, got)
    call report(abs(got(1) - old % rhs_coeff(area, delta, kappa)) < tol, &
         & name // ': the diffusive constant matches', nfail)

    call bc % adv_lhs_coefficients(m, vn, got)
    call report(abs(got(1) - old % adv_lhs_coeff(area, delta, vn)) < tol, &
         & name // ': the advective diagonal matches', nfail)

    call bc % adv_rhs_coefficients(m, vn, got)
    call report(abs(got(1) - old % adv_rhs_coeff(area, delta, vn)) < tol, &
         & name // ': the advective constant matches', nfail)

  end subroutine fidelity_case

  !===================================================================!
  ! THE OPERATOR ROAD. The condition's coefficients and stored value
  ! make the calculus carry the boundary flux: the difference between
  ! the row with the condition and the row without it equals
  ! lhs*q_p - rhs - the old row's own boundary term, sign and all -
  ! both sides computed independently.
  !===================================================================!

  subroutine check_operator_road(nfail)

    integer, intent(inout) :: nfail

    type(mesh) :: m
    type(robin_condition) :: bc
    type(old_condition)   :: old
    type(differential_operator) :: with_wall, without
    type(field) :: state
    type(support) :: cells
    class(graph_field), allocatable :: y
    real(dp), allocatable :: cin(:), bin(:), rows_with(:), rows_without(:)
    real(dp) :: c(3), b(3), expected
    real(dp), parameter :: kappa = 1.7_dp, kint = 1.3_dp
    real(dp), parameter :: q1 = 2.0_dp, q2 = 4.0_dp
    integer :: v

    m   = hand_mesh()
    bc  = dirichlet('in', 5.0_dp)
    old = old_robin(1.0_dp, 0.0_dp, 5.0_dp)

    call bc % operator_coefficients(m, kappa, cin)
    call bc % boundary_values(m, bin)

    ! Interior conductivity on e1; the condition on e2; nothing on e3.
    c = [kint * 2.0_dp, cin(1), 0.0_dp]
    b = [0.0_dp, bin(1), 0.0_dp]

    cells = support(GRAPH_SIDE_VERTEX, [(v, v = 1, 2)])
    state = field('q', cells)
    call state % set_real_vector([q1, q2])

    with_wall = vertex_differential_operator(order=2, coefficients=c, &
         & spacings=[0.8_dp, 0.3_dp, 0.6_dp], boundary_values=b)
    call with_wall % apply(m, [state], y)
    call y % get_real_vector(rows_with)

    without = vertex_differential_operator(order=2, &
         & coefficients=[c(1), 0.0_dp, 0.0_dp], &
         & spacings=[0.8_dp, 0.3_dp, 0.6_dp])
    call without % apply(m, [state], y)
    call y % get_real_vector(rows_without)

    expected = old % lhs_coeff(0.7_dp, 0.3_dp, kappa) * q1 &
         &     - old % rhs_coeff(0.7_dp, 0.3_dp, kappa)

    call report(abs((rows_with(1) - rows_without(1)) - expected) < 1.0d-11, &
         & 'the row difference is the eliminated boundary flux', nfail)

    call report(abs(rows_with(2) - rows_without(2)) < 1.0d-14, &
         & 'and the far cell feels nothing from it', nfail)

  end subroutine check_operator_road

  !===================================================================!
  ! THE CONDUCTION LAW. keff = n^T K n from the mesh's own normals:
  ! isotropic k answers k on every unit normal; a diagonal tensor
  ! answers the component its normal selects. The dictionary
  ! coefficient keff*area lives on interior faces only.
  !===================================================================!

  subroutine check_conduction_law(nfail)

    integer, intent(inout) :: nfail

    type(mesh) :: m
    type(conduction) :: law
    real(dp), allocatable :: keff(:), c(:)
    real(dp) :: k(3, 3)

    m = hand_mesh()

    law = conduction(2.5_dp)
    call law % normal_conductivity(m, keff)
    call report(all(abs(keff - 2.5_dp) < 1.0d-14), &
         & 'an isotropic material answers k on every unit normal', nfail)

    k = 0.0_dp
    k(1, 1) = 2.0_dp
    k(2, 2) = 3.0_dp
    k(3, 3) = 4.0_dp
    law = conduction(k)
    call law % normal_conductivity(m, keff)
    call report(all(abs(keff - 2.0_dp) < 1.0d-14), &
         & 'a diagonal tensor answers the part its normal selects', nfail)

    call law % edge_coefficients(m, c)
    call report(abs(c(1) - 2.0_dp * 2.0_dp) < 1.0d-14, &
         & 'the interior coefficient is keff times the area', nfail)
    call report(c(2) == 0.0_dp .and. c(3) == 0.0_dp, &
         & 'and the headless faces carry none: theirs come from a condition', nfail)

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

  subroutine check_advection_law(nfail)

    integer, intent(inout) :: nfail

    type(mesh) :: m
    type(advection) :: flow
    type(balance) :: sums
    type(field) :: state
    type(support) :: cells
    class(graph_field), allocatable :: y
    real(dp), allocatable :: vn(:), c(:), got(:)
    real(dp), parameter :: q1 = 2.0_dp, q2 = 4.0_dp
    real(dp) :: area, wp, wn
    integer :: v

    m = hand_mesh()
    area = 2.0_dp

    cells = support(GRAPH_SIDE_VERTEX, [(v, v = 1, 2)])
    state = field('q', cells)
    call state % set_real_vector([q1, q2])

    flow = advection([1.5_dp, 0.0_dp, 0.0_dp])
    call flow % normal_speed(m, vn)
    call report(abs(vn(1) - 1.5_dp) < 1.0d-14, &
         & 'the normal speed is v dot n', nfail)

    call flow % edge_coefficients(m, c)
    call report(abs(c(1) - 1.5_dp * area) < 1.0d-14 .and. c(2) == 0.0_dp, &
         & 'the coefficient is vn times area, interior only', nfail)

    ! Upwind, flow along the edge: the tail carries it.
    sums = balance(edge_terms=[edge_differential_operator(order=0, &
         & coefficients=c, one_sided=.true.)])
    call sums % apply(m, [state], y)
    call y % get_real_vector(got)

    wp = max(vn(1), 0.0_dp)
    wn = min(vn(1), 0.0_dp)
    call report(abs(got(1) - (-area * (wp * q1 + wn * q2))) < 1.0d-12 .and. &
         &      abs(got(2) - (+area * (wp * q1 + wn * q2))) < 1.0d-12, &
         & 'upwind with the flow: the old rows, both cells', nfail)

    ! Upwind, flow against the edge: the head carries it.
    flow = advection([-1.5_dp, 0.0_dp, 0.0_dp])
    call flow % edge_coefficients(m, c)
    sums = balance(edge_terms=[edge_differential_operator(order=0, &
         & coefficients=c, one_sided=.true.)])
    call sums % apply(m, [state], y)
    call y % get_real_vector(got)

    wp = max(-1.5_dp, 0.0_dp)
    wn = min(-1.5_dp, 0.0_dp)
    call report(abs(got(1) - (-area * (wp * q1 + wn * q2))) < 1.0d-12, &
         & 'upwind against the flow: the head is upstream', nfail)

    ! Central: both ends, evenly - the old half weights.
    flow = advection([1.5_dp, 0.0_dp, 0.0_dp])
    call flow % edge_coefficients(m, c)
    sums = balance(edge_terms=[edge_differential_operator(order=0, &
         & coefficients=c, one_sided=.false.)])
    call sums % apply(m, [state], y)
    call y % get_real_vector(got)

    wp = 0.5_dp * 1.5_dp
    wn = 0.5_dp * 1.5_dp
    call report(abs(got(1) - (-area * (wp * q1 + wn * q2))) < 1.0d-12, &
         & 'central: the old half weights, exactly', nfail)

  end subroutine check_advection_law

  !===================================================================!
  ! THE STATEMENT. The constitution's whole sentence in one call:
  ! a material, two walls, a mesh - and the compiled operator closes
  ! under the house solver, monotone between its held values.
  !===================================================================!

  subroutine check_the_statement_speaks(nfail)

    integer, intent(inout) :: nfail

    type(mesh) :: m
    type(stencil_operator) :: op
    type(gmres) :: gm
    real(dp), allocatable :: g(:), rhs(:), x(:)
    real(dp) :: achieved

    m = hand_mesh()

    op = diffusion_statement(m, conduction(1.7_dp), &
         & [dirichlet('in', 0.0_dp), dirichlet('out', 10.0_dp)])

    call gm % attach(op, m)
    gm % tolerance = 1.0d-12

    call gm % constant(g)
    rhs = -g

    allocate(x(2))
    x = 0.0_dp
    call gm % solve(rhs, x, achieved)

    call report(achieved < 1.0d-10, &
         & 'the statement closes under the house solver', nfail)
    call report(x(1) > 0.0_dp .and. x(1) < x(2) .and. x(2) < 10.0_dp, &
         & 'and the answer sits monotone between the held walls', nfail)

  end subroutine check_the_statement_speaks

end program test_graph_constitution
