!=====================================================================!
! The minimization suite: rung 4's acceptance.
!
! The whole tower stands in one problem: heat conduction on a chain
! of three cells, walls held at 0 and 10,
!
!       0 |--- (1) --- (2) --- (3) ---| 10        k = 1, A = 1
!         d 0.5     1       1     0.5
!
! whose exact answer is q = [5/3, 5, 25/3]. The mesh seats the
! measurements, the constitution supplies the coefficients, the
! calculus builds the rows, and the minimization drives the residual
! down - every level doing only its own job, and the solver speaking
! only its own words: matvec, inner product, norm, diagonal, sweep
! order, each a delegation to an engine seat.
!=====================================================================!

module cubic_statement_fixture

  use iso_fortran_env, only : dp => REAL64
  use graph_grammar  , only : graph, graph_operation
  use graph_field_calculus, only : graph_field
  use graph_calculus , only : GRAPH_SIDE_VERTEX
  ! An action names a domain and counts it. It asks no membership,
  ! so it holds no map: identity and count is the whole of it.
  use fractal_graph      , only : set_graph => graph
  use class_graph_field  , only : field
  use class_graph_differential_operator, only : differential_operator

  implicit none

  private
  public :: cubic_statement

  !===================================================================!
  ! A nonlinear statement for newton to chew: the chain's rows with
  ! a small cube AGAINST them - a stable reaction, so the jacobian
  ! keeps the rows' own sign and the root is one.
  !===================================================================!

  type, extends(graph_operation) :: cubic_statement

     type(differential_operator) :: linear_part
     real(dp) :: strength = 0.0_dp

   contains

     procedure :: name   => cubic_name
     procedure :: domain => cubic_domain
     procedure :: apply  => cubic_apply

  end type cubic_statement

contains

  pure function cubic_name(this) result(name)
    class(cubic_statement), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'cubic statement'
  end function cubic_name

  subroutine cubic_domain(this, input_graph, domain, nentries)
    class(cubic_statement), intent(in)     :: this
    class(graph), intent(in)               :: input_graph
    type(set_graph), intent(out) :: domain
    integer        , intent(out) :: nentries
    domain   = input_graph % all_vertices()
    nentries = input_graph % num_vertices()
  end subroutine cubic_domain

  subroutine cubic_apply(this, input_graph, input_data, output)

    class(cubic_statement), intent(in)             :: this
    class(graph), intent(in)                       :: input_graph
    class(graph_field), intent(in), optional       :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output

    type(set_graph)   :: cells
    type(field)   :: out
    real(dp), allocatable :: q(:), y(:)
    integer :: nv, v

    call this % linear_part % apply(input_graph, input_data, output)
    call output % get_real_vector(y)

    nv = input_graph % num_vertices()
    if (present(input_data)) then
       call input_data(1) % get_real_vector(q)
       do v = 1, min(nv, size(q))
          y(v) = y(v) - this % strength * q(v)**3
       end do
    end if

    cells = input_graph % vertex_set()
    out = field('cubic', cells, nv)
    call out % set_real_vector(y)
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine cubic_apply

end module cubic_statement_fixture

program test_graph_minimization

  use iso_fortran_env, only : dp => REAL64
  use fractal_graph  , only : set_graph => graph
  use graph_set_map  , only : set_map
  use graph_label_map, only : label_map
  use graph_inclusion_map, only : inclusion_map
  use graph_grammar  , only : graph
  use class_graph_mesh   , only : mesh
  use class_mesh_builder , only : mesh_from_gmsh
  use class_robin_condition, only : robin_condition, dirichlet
  use class_conduction     , only : conduction
  use class_graph_differential_operator, only : differential_operator
  use class_graph_differential_operator, only : vertex_differential_operator
  use class_graph_jacobi   , only : jacobi
  use class_graph_conjugate_gradient, only : conjugate_gradient
  use class_graph_gauss_seidel, only : gauss_seidel
  use class_graph_gmres    , only : gmres
  use class_graph_newton   , only : newton
  use class_graph_differential_operator, only : edge_differential_operator
  use class_graph_balance  , only : balance
  use cubic_statement_fixture, only : cubic_statement

  implicit none

  integer :: nfail

  nfail = 0

  call check_solver_words(nfail)
  call check_hand_problem(nfail)
  call check_real_mesh(nfail)
  call check_sweeping_family(nfail)
  call check_gmres_family(nfail)
  call check_newton(nfail)

  write(*, '(a)') ' ============================================='
  if (nfail == 0) then
     write(*, '(a)') ' all minimization checks passed'
  else
     write(*, '(a, i0, a)') ' ', nfail, ' minimization checks FAILED'
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
  ! The three-cell conduction problem, assembled from the levels.
  !===================================================================!

  type(mesh) function chain_mesh() result(m)

    m = mesh(3, tails=[1, 2, 1, 3], heads=[2, 3, 0, 0], &
         & volumes      = [1.0_dp, 1.0_dp, 1.0_dp], &
         & cell_centers = [0.5_dp, 0.0_dp, 0.0_dp, &
         &                 1.5_dp, 0.0_dp, 0.0_dp, &
         &                 2.5_dp, 0.0_dp, 0.0_dp], &
         & areas        = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
         & deltas       = [1.0_dp, 1.0_dp, 0.5_dp, 0.5_dp], &
         & normals      = [ 1.0_dp, 0.0_dp, 0.0_dp, &
         &                  1.0_dp, 0.0_dp, 0.0_dp, &
         &                 -1.0_dp, 0.0_dp, 0.0_dp, &
         &                  1.0_dp, 0.0_dp, 0.0_dp], &
         & face_centers = [1.0_dp, 0.0_dp, 0.0_dp, &
         &                 2.0_dp, 0.0_dp, 0.0_dp, &
         &                 0.0_dp, 0.0_dp, 0.0_dp, &
         &                 3.0_dp, 0.0_dp, 0.0_dp], &
         & weights      = [0.5_dp, 0.5_dp, 1.0_dp, 1.0_dp], &
         & etags        = [character(len=4) :: '', '', 'west', 'east'])

  end function chain_mesh

  !===================================================================!
  ! The assembly: the constitution's coefficients placed into the
  ! operator's argument list - conduction inside, one condition per
  ! wall - and nothing else.
  !===================================================================!

  type(differential_operator) function chain_operator(m) result(op)

    type(mesh), intent(in) :: m

    type(robin_condition) :: west, east
    real(dp), allocatable :: c(:), b(:)

    call assemble(m, conduction(1.0_dp), &
         & [dirichlet('west', 0.0_dp), dirichlet('east', 10.0_dp)], &
         & 1.0_dp, c, b)

    op = vertex_differential_operator(order=2, coefficients=c, &
         & spacings=[1.0_dp, 1.0_dp, 0.5_dp, 0.5_dp], boundary_values=b)

  end function chain_operator

  subroutine assemble(m, law, conditions, kappa, c, b)

    type(mesh), intent(in)            :: m
    type(conduction), intent(in)      :: law
    type(robin_condition), intent(in) :: conditions(:)
    real(dp), intent(in)              :: kappa
    real(dp), allocatable, intent(out) :: c(:), b(:)

    type(set_graph)       :: members
    type(set_map)         :: sets
    type(label_map)       :: labels
    type(inclusion_map)   :: inclusions
    integer , allocatable :: face(:)
    real(dp), allocatable :: cw(:), bw(:)
    integer :: k, f, e

    call law % edge_coefficients(m, c)
    allocate(b(size(c)))
    b = 0.0_dp

    do k = 1, size(conditions)
       call conditions(k) % faces(m, sets, labels, inclusions, members)
       call conditions(k) % operator_coefficients(m, kappa, cw)
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
  ! the probed diagonal is the stencil's own, and the colouring
  ! never gives neighbours one colour.
  !===================================================================!

  subroutine check_solver_words(nfail)

    integer, intent(inout) :: nfail

    type(mesh) :: m
    type(jacobi) :: js
    real(dp), allocatable :: d(:)
    integer , allocatable :: colours(:), nbrs(:)
    logical :: proper
    integer :: v, i

    m = chain_mesh()
    call js % attach(chain_operator(m), m, m % vertex_set(), m % num_vertices(), coupling = m)

    call js % diagonal(d)
    call report(all(abs(d - [-3.0_dp, -2.0_dp, -3.0_dp]) < 1.0d-12), &
         & 'the probed diagonal is the stencil diagonal, by hand', nfail)

    call js % sweep_order(colours)
    proper = .true.
    do v = 1, m % num_vertices()
       call m % adjacent_vertices(v, nbrs)
       do i = 1, size(nbrs)
          if (colours(nbrs(i)) == colours(v)) proper = .false.
       end do
    end do
    call report(proper, 'the sweep order never gives neighbours one colour', nfail)

    call report(abs(js % inner_product([1.0_dp, 2.0_dp, 3.0_dp], &
         &                             [2.0_dp, 1.0_dp, 0.0_dp]) - 4.0_dp) &
         & < 1.0d-14, 'the inner product is the measured sum', nfail)

    call report(abs(js % norm([3.0_dp, 4.0_dp, 0.0_dp]) - 5.0_dp) < 1.0d-14, &
         & 'the norm is the norm', nfail)

  end subroutine check_solver_words

  !===================================================================!
  ! The hand problem, driven to its exact answer.
  !===================================================================!

  subroutine check_hand_problem(nfail)

    integer, intent(inout) :: nfail

    type(mesh) :: m
    type(jacobi) :: js
    real(dp), allocatable :: g(:), rhs(:), x(:), y(:)
    real(dp) :: achieved
    real(dp), parameter :: exact(3) = [5.0_dp/3.0_dp, 5.0_dp, 25.0_dp/3.0_dp]

    m = chain_mesh()
    call js % attach(chain_operator(m), m, m % vertex_set(), m % num_vertices(), coupling = m)
    js % max_iterations = 5000
    js % tolerance      = 1.0d-12

    ! The equation action(q) = 0 reads matvec(q) = -constant.
    call js % constant(g)
    rhs = -g

    allocate(x(3))
    x = 0.0_dp
    call js % solve(rhs, x, achieved)

    call report(achieved < 1.0d-10, 'jacobi drives the residual down', nfail)
    call report(all(abs(x - exact) < 1.0d-8), &
         & 'and lands on the exact answer: 5/3, 5, 25/3', nfail)

    call js % matvec(x, y)
    call report(js % norm(rhs - y) < 1.0d-10, &
         & 'the answer satisfies the assembled equation', nfail)

  end subroutine check_hand_problem

  !===================================================================!
  ! The real mesh. A symmetric positive operator - the conduction
  ! coefficients negated - and a consistent right hand side
  ! manufactured from a known state; conjugate gradient must return
  ! a state the operator cannot tell apart from it.
  !===================================================================!

  subroutine check_real_mesh(nfail)

    integer, intent(inout) :: nfail

    type(mesh) :: m
    type(conjugate_gradient) :: cg
    type(conduction) :: law
    real(dp), allocatable :: c(:), deltas(:), xref(:), b(:), x(:), y(:)
    real(dp) :: achieved
    integer :: v

    m = mesh_from_gmsh('../square-10.msh')

    law = conduction(1.0_dp)
    call law % edge_coefficients(m, c)
    c = -c

    block
      use class_graph_field, only : field
      type(field) :: fd
      fd = m % face_delta()
      call fd % get_real_vector(deltas)
    end block

    call cg % attach(vertex_differential_operator(order=2, &
         & coefficients=c, spacings=deltas), m, m % vertex_set(), &
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
         & nfail)

  end subroutine check_real_mesh

  !===================================================================!
  ! Gauss-seidel and its omega: the same chain, the same exact
  ! answer, the sweep ordered by colour. SOR is not another type -
  ! it is this one at omega away from one.
  !===================================================================!

  subroutine check_sweeping_family(nfail)

    integer, intent(inout) :: nfail

    type(mesh) :: m
    type(gauss_seidel) :: gs
    real(dp), allocatable :: g(:), rhs(:), x(:)
    real(dp) :: achieved
    real(dp), parameter :: exact(3) = [5.0_dp/3.0_dp, 5.0_dp, 25.0_dp/3.0_dp]

    m = chain_mesh()
    call gs % attach(chain_operator(m), m, m % vertex_set(), m % num_vertices(), coupling = m)
    gs % max_iterations = 2000
    gs % tolerance      = 1.0d-12

    call gs % constant(g)
    rhs = -g

    allocate(x(3))
    x = 0.0_dp
    call gs % solve(rhs, x, achieved)
    call report(all(abs(x - exact) < 1.0d-8), &
         & 'gauss-seidel sweeps by colour to the exact answer', nfail)

    gs % omega = 1.2_dp
    x = 0.0_dp
    call gs % solve(rhs, x, achieved)
    call report(all(abs(x - exact) < 1.0d-8), &
         & 'and over-relaxed it is sor, a parameter, not a type', nfail)

  end subroutine check_sweeping_family

  !===================================================================!
  ! GMRES: the chain again, then an unsymmetric statement -
  ! advection riding diffusion through the balance - where two
  ! different solvers must meet on one answer.
  !===================================================================!

  subroutine check_gmres_family(nfail)

    integer, intent(inout) :: nfail

    type(mesh) :: m
    type(gmres)  :: gm
    type(jacobi) :: js
    type(balance) :: statement
    real(dp), allocatable :: g(:), rhs(:), x(:), xj(:)
    real(dp) :: achieved
    real(dp), parameter :: exact(3) = [5.0_dp/3.0_dp, 5.0_dp, 25.0_dp/3.0_dp]

    m = chain_mesh()

    call gm % attach(chain_operator(m), m, m % vertex_set(), m % num_vertices())
    gm % tolerance = 1.0d-12
    call gm % constant(g)
    rhs = -g
    allocate(x(3))
    x = 0.0_dp
    call gm % solve(rhs, x, achieved)
    call report(all(abs(x - exact) < 1.0d-8), &
         & 'gmres lands on the exact answer', nfail)

    ! Diffusion and upwind advection in one balance: unsymmetric.
    statement = balance(edge_terms=[ &
         & edge_differential_operator(order=1, &
         &      coefficients=[1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp], &
         &      spacings=[1.0_dp, 1.0_dp, 0.5_dp, 0.5_dp], &
         &      boundary_values=[0.0_dp, 0.0_dp, 0.0_dp, 10.0_dp]), &
         & edge_differential_operator(order=0, &
         &      coefficients=[0.4_dp, 0.4_dp, 0.0_dp, 0.0_dp], &
         &      one_sided=.true.)])

    call gm % attach(statement, m, m % vertex_set(), m % num_vertices())
    call gm % constant(g)
    rhs = -g
    x = 0.0_dp
    call gm % solve(rhs, x, achieved)
    call report(achieved < 1.0d-10, &
         & 'gmres closes the unsymmetric statement', nfail)

    call js % attach(statement, m, m % vertex_set(), m % num_vertices(), coupling = m)
    js % max_iterations = 5000
    js % tolerance      = 1.0d-12
    allocate(xj(3))
    xj = 0.0_dp
    call js % solve(rhs, xj, achieved)
    call report(all(abs(x - xj) < 1.0d-7), &
         & 'and two different solvers meet on one answer', nfail)

  end subroutine check_gmres_family

  !===================================================================!
  ! Newton over the linear family. A cubic term rides the chain
  ! statement; newton linearizes by directional differences, hands
  ! each linear question to the gmres it governs, and drives the
  ! nonlinear residual to zero.
  !===================================================================!

  subroutine check_newton(nfail)

    integer, intent(inout) :: nfail

    type(mesh) :: m
    type(newton) :: ns
    type(cubic_statement) :: action
    real(dp), allocatable :: q(:)
    real(dp) :: achieved

    m = chain_mesh()

    action % linear_part = chain_operator(m)
    action % strength    = 0.05_dp

    allocate(ns % inner, source=gmres())
    ns % inner % tolerance = 1.0d-12
    ns % tolerance = 1.0d-7

    call ns % attach(action, m, m % vertex_set(), m % num_vertices())

    allocate(q(3))
    q = 0.0_dp
    call ns % solve([0.0_dp, 0.0_dp, 0.0_dp], q, achieved)

    call report(achieved < 1.0d-7, &
         & 'newton drives the nonlinear residual to the difference floor', nfail)
    call report(q(1) > 0.0_dp .and. q(1) < q(2) .and. q(2) < q(3), &
         & 'and the heated chain still rises monotonically', nfail)

  end subroutine check_newton

end program test_graph_minimization
