!=====================================================================!
! Multigrid consistency and convergence tests.
!
! An eight-cell conduction chain, with endpoint values 0 and 10, compiled
! by the fitted balance into one stencil operator; blocks of two
! from pairwise coarsening. Four assertions:
!
!   1. the aggregates cover every cell exactly once
!   2. THE COMMUTATION SQUARE: a coarse vector is
!      reproduced identically by the coarsened statement and by the
!      fine one - solve_coarse(R(A(P e))) returns e
!   3. the two-grid cycle and a fine-grid GMRES solve
!      converge to the same field
!   4. the two-grid residual after one cycle is
!      far below one smoothing pass alone
!=====================================================================!

module coarse_solve_fixture
  use iso_fortran_env, only : dp => REAL64
  use operation_minimization, only : minimizer, state
  use operation_action      , only : operation
  use view_directed         , only : directed_graph
  use graph_fractal         , only : graph
  use field_stored          , only : stored_field
  implicit none
  integer :: num_coarse_solves = 0
  ! the versions of the coarse statements received, in order of statement
  integer :: num_statements_received = 0
  integer :: statement_versions(16) = 0
  type, extends(minimizer) :: counted_coarse_solver
  contains
    procedure :: state => record_statement
    procedure :: solve => solve_coarse
  end type counted_coarse_solver
contains
  subroutine record_statement(this, action, context, unknown_domain, num_unknowns, &
       & num_components, coupling, stored_inputs)
    class(counted_coarse_solver), intent(inout)        :: this
    class(operation)            , intent(in)           :: action
    class(directed_graph)       , intent(in)           :: context
    type(graph)                 , intent(in)           :: unknown_domain
    integer                     , intent(in)           :: num_unknowns
    integer                     , intent(in), optional :: num_components
    class(directed_graph)       , intent(in), optional :: coupling
    type(stored_field)          , intent(in), optional :: stored_inputs(:)
    num_statements_received = num_statements_received + 1
    statement_versions(num_statements_received) = action % version()
    call state(this, action, context, unknown_domain, num_unknowns, num_components, coupling, stored_inputs)
  end subroutine record_statement
  subroutine solve_coarse(this, rhs, x, achieved)
    class(counted_coarse_solver), intent(inout) :: this
    real(dp), intent(in) :: rhs(:)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out) :: achieved
    real(dp), allocatable :: residual(:)
    num_coarse_solves = num_coarse_solves + 1
    call this % initialize_residual_history()
    call this % record_residual_norm(this % norm(rhs))
    x = rhs / 6.0_dp
    call this % imbalance(rhs, x, residual)
    achieved = this % norm(residual)
    call this % record_result(achieved, 1)
  end subroutine solve_coarse
end module coarse_solve_fixture

program test_graph_multigrid

  use iso_fortran_env, only : dp => REAL64
  use view_directed, only : directed_graph
  use view_mesh   , only : mesh
  use operation_stencil, only : stencil
  use operation_fitted_balance, only : fitted_balance_stencil
  use field_forms       , only : polynomial_form
  use transform_coarsener, only : coarsener, COARSEN_PAIRWISE
  use operation_multigrid, only : multigrid
  use operation_jacobi   , only : jacobi
  use operation_gmres    , only : gmres
  use operation_minimization, only : solve_result, SOLVE_EXHAUSTED
  use coarse_solve_fixture, only : counted_coarse_solver, num_coarse_solves, &
       & num_statements_received, statement_versions

  implicit none

  integer :: num_failures

  num_failures = 0

  call check_blocks_cover_once(num_failures)
  call check_commutation_square(num_failures)
  call check_multigrid_gmres_equivalence(num_failures)
  call check_one_cycle(num_failures)
  call check_statement_versions(num_failures)

  write(*, '(a)') ' ============================================='
  if (num_failures == 0) then
     write(*, '(a)') ' all multigrid checks passed'
  else
     write(*, '(a, i0, a)') ' ', num_failures, ' multigrid checks FAILED'
     error stop 1
  end if

contains

  subroutine check_one_cycle(num_failures)
    integer, intent(inout) :: num_failures
    type(stencil) :: a
    type(multigrid) :: solver
    type(jacobi) :: smoother
    type(solve_result) :: outcome
    real(dp) :: x(2), rhs(2), achieved

    a = stencil([1, 1, 2, 2], [1, 2, 1, 2], [2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], &
         & [0.0_dp, 0.0_dp], 'one cycle')
    smoother % max_iterations = 1
    allocate(solver % smoother, source=smoother)
    allocate(solver % coarse, source=counted_coarse_solver())
    call solver % state(a, a % pattern, a % pattern % vertex_set(), 2, coupling=a % pattern)
    call solver % setup([1, 1])
    solver % max_iterations = 1
    solver % tolerance = 1.0e-12_dp
    rhs = [1.0_dp, 0.0_dp]
    x = 0.0_dp
    num_coarse_solves = 0
    call solver % solve(rhs, x, achieved)
    outcome = solver % result()
    call report(num_coarse_solves == 1 .and. outcome % iterations == 1 .and. achieved < 0.2_dp, &
         & 'one multigrid cycle reaches the coarse solve and completes both smoothing passes', num_failures)

    x = 0.0_dp
    num_coarse_solves = 0
    call solver % solve([0.0_dp, 0.0_dp], x, achieved)
    outcome = solver % result()
    call report(num_coarse_solves == 0 .and. outcome % iterations == 0 .and. outcome % converged(), &
         & 'an initially solved multigrid system takes no cycle', num_failures)

    solver % max_iterations = 0
    call solver % solve(rhs, x, achieved)
    outcome = solver % result()
    call report(num_coarse_solves == 0 .and. all(x == 0.0_dp) .and. outcome % reason == SOLVE_EXHAUSTED, &
         & 'a zero multigrid iteration limit performs no cycle', num_failures)
  end subroutine check_one_cycle

  !===================================================================!
  ! The version of the k-th coarse statement of one multigrid object is
  ! k: the sequence belongs to the object, so two objects stated
  ! alternately receive 1, 1, 2, 2, 3 and a copy continues its own
  ! sequence from the count it was copied with.
  !===================================================================!

  subroutine check_statement_versions(num_failures)
    integer, intent(inout) :: num_failures
    type(stencil) :: a
    type(multigrid) :: first, second, copy
    type(jacobi) :: smoother

    a = stencil([1, 1, 2, 2], [1, 2, 1, 2], [2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], &
         & [0.0_dp, 0.0_dp], 'statement versions')
    smoother % max_iterations = 1
    allocate(first % smoother, source=smoother)
    allocate(first % coarse, source=counted_coarse_solver())
    allocate(second % smoother, source=smoother)
    allocate(second % coarse, source=counted_coarse_solver())
    num_statements_received = 0
    statement_versions = 0

    call first % state(a, a % pattern, a % pattern % vertex_set(), 2, coupling=a % pattern)
    call first % setup([1, 1])
    call second % state(a, a % pattern, a % pattern % vertex_set(), 2, coupling=a % pattern)
    call second % setup([1, 1])
    ! aggregates are stored, so a further state re-forms the coarse statement
    call first % state(a, a % pattern, a % pattern % vertex_set(), 2, coupling=a % pattern)
    call second % state(a, a % pattern, a % pattern % vertex_set(), 2, coupling=a % pattern)
    call first % state(a, a % pattern, a % pattern % vertex_set(), 2, coupling=a % pattern)

    call report(num_statements_received == 5 .and. all(statement_versions(1:5) == [1, 1, 2, 2, 3]), &
         & 'two multigrid objects stated alternately number their coarse statements independently', num_failures)
    call report(first % num_statements == 3 .and. second % num_statements == 2, &
         & 'each multigrid object counts its own coarse statements', num_failures)

    copy = first
    call copy % state(a, a % pattern, a % pattern % vertex_set(), 2, coupling=a % pattern)
    call first % state(a, a % pattern, a % pattern % vertex_set(), 2, coupling=a % pattern)
    call report(num_statements_received == 7 .and. all(statement_versions(6:7) == [4, 4]) .and. &
         & copy % num_statements == 4 .and. first % num_statements == 4, &
         & 'a copied multigrid continues its own statement sequence from the copied count', num_failures)
  end subroutine check_statement_versions

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
  ! The chain: eight cells, unit faces, endpoint values 0 and 10.
  !===================================================================!

  type(mesh) function chain_mesh() result(m)

    integer :: kf

    m = mesh(8, &
         & tails=[1, 2, 3, 4, 5, 6, 7,  1, 8], &
         & heads=[2, 3, 4, 5, 6, 7, 8,  0, 0], &
         & volumes      = [(1.0_dp, kf = 1, 8)], &
         & cell_centres = [(real(kf - 1, dp) + 0.5_dp, 0.0_dp, 0.0_dp, &
         &                  kf = 1, 8)], &
         & areas        = [(1.0_dp, kf = 1, 9)], &
         & deltas       = [(1.0_dp, kf = 1, 7), 0.5_dp, 0.5_dp], &
         & normals      = [(1.0_dp, 0.0_dp, 0.0_dp, kf = 1, 7), &
         &                 -1.0_dp, 0.0_dp, 0.0_dp, &
         &                  1.0_dp, 0.0_dp, 0.0_dp], &
         & face_centres = [(real(kf, dp), 0.0_dp, 0.0_dp, kf = 1, 7), &
         &                  0.0_dp, 0.0_dp, 0.0_dp, &
         &                  8.0_dp, 0.0_dp, 0.0_dp], &
         & weights      = [(0.5_dp, kf = 1, 9)])

  end function chain_mesh

  type(stencil) function chain_statement(m) result(op)

    type(mesh), intent(in) :: m

    real(dp) :: vb(9), farea(9)
    integer :: e

    do e = 1, 9
       farea(e) = 1.0_dp
       vb(e)    = 0.0_dp
    end do
    vb(9) = 10.0_dp

    op = fitted_balance_stencil(m, polynomial_form(), farea, &
         & boundary_values=vb)

  end function chain_statement

  !===================================================================!
  ! The aggregates partition the cells.
  !===================================================================!

  subroutine check_blocks_cover_once(num_failures)

    integer, intent(inout) :: num_failures

    type(mesh) :: m
    type(coarsener) :: c
    integer, allocatable :: assignment(:)
    integer :: nb

    m = chain_mesh()
    c = coarsener(COARSEN_PAIRWISE)
    call c % blocks(m, assignment, nb)

    call report(size(assignment) == 8 .and. nb >= 2, &
         & 'every cell belongs to one block, with at least two blocks', num_failures)
    call report(minval(assignment) >= 1 .and. maxval(assignment) == nb, &
         & 'the blocks cover once: no gaps, no strays', num_failures)

  end subroutine check_blocks_cover_once

  !===================================================================!
  ! The commutation square, exact: for any block
  ! field e, push it through the fine statement, gather, and the
  ! coarsened statement must reproduce e.
  !===================================================================!

  subroutine check_commutation_square(num_failures)

    integer, intent(inout) :: num_failures

    type(mesh) :: m
    type(multigrid) :: mg
    type(coarsener) :: c
    integer, allocatable :: assignment(:)
    real(dp), allocatable :: e_blocks(:), e_fine(:), r_fine(:), rc(:), ec(:)
    real(dp) :: coarse_residual_norm
    integer :: nb, v, b

    m = chain_mesh()
    c = coarsener(COARSEN_PAIRWISE)
    call c % blocks(m, assignment, nb)

    allocate(mg % smoother, source=jacobi())
    allocate(mg % coarse  , source=gmres())
    call mg % state(chain_statement(m), m, m % vertex_set(), &
         & m % num_vertices())
    call mg % setup(assignment)
    mg % coarse % tolerance = 1.0d-13

    allocate(e_blocks(nb), e_fine(8), rc(nb), ec(nb))
    do b = 1, nb
       e_blocks(b) = real(mod(3 * b, 5) + 1, dp)
    end do

    ! P e: every cell takes its block's value. A(P e): the fine
    ! statement. R: gathered onto blocks.
    do v = 1, 8
       e_fine(v) = e_blocks(assignment(v))
    end do
    call mg % matvec(e_fine, r_fine)
    rc = 0.0_dp
    do v = 1, 8
       rc(assignment(v)) = rc(assignment(v)) + r_fine(v)
    end do

    ec = 0.0_dp
    call mg % coarse % solve(rc, ec, coarse_residual_norm)

    call report(all(abs(ec - e_blocks) < 1.0d-9), &
         & 'the commutation square is exact: coarsen-then-solve returns e', num_failures)

  end subroutine check_commutation_square

  !===================================================================!
  ! Multigrid and fine-grid GMRES solution equivalence, and the
  ! detour pays.
  !===================================================================!

  subroutine check_multigrid_gmres_equivalence(num_failures)

    integer, intent(inout) :: num_failures

    type(mesh) :: m
    type(multigrid) :: mg
    type(gmres)     :: fine_solver
    type(coarsener) :: c
    integer, allocatable :: assignment(:)
    real(dp), allocatable :: g(:), rhs(:), x(:), xd(:)
    real(dp) :: achieved
    integer :: nb

    m = chain_mesh()
    c = coarsener(COARSEN_PAIRWISE)
    call c % blocks(m, assignment, nb)

    allocate(mg % smoother, source=jacobi())
    allocate(mg % coarse  , source=gmres())
    call mg % state(chain_statement(m), m, m % vertex_set(), &
         & m % num_vertices())
    call mg % setup(assignment)

    mg % smoother % max_iterations = 3
    mg % smoother % tolerance      = 0.0_dp
    select type (s => mg % smoother)
    type is (jacobi)
       s % omega = 0.7_dp
    end select
    mg % coarse   % tolerance      = 1.0d-13
    mg % tolerance                 = 1.0d-11
    mg % max_iterations            = 60

    g = mg % affine
    rhs = -g

    allocate(x(8))
    x = 0.0_dp
    call mg % solve(rhs, x, achieved)

    call report(achieved < 1.0d-9, &
         & 'the two-grid cycle closes the chain statement', num_failures)

    call fine_solver % state(chain_statement(m), m, m % vertex_set(), &
         & m % num_vertices())
    fine_solver % tolerance = 1.0d-12
    allocate(xd(8))
    xd = 0.0_dp
    call fine_solver % solve(rhs, xd, achieved)

    call report(all(abs(x - xd) < 1.0d-7), &
         & 'the multigrid and fine-grid GMRES solutions agree', num_failures)

    !----------------------------------------------------------------!
    ! State replacement and default component regressions.
    !
    ! STATE MAY BE REPEATED. Newton calls it once per iteration, and
    ! it declares the solver's number domain each time - which a graph
    ! refuses, because graph identity is immutable. The old counted_set
    ! constructor declared a distinct number domain per state, so
    ! that semantics is retained: reset to an undeclared graph, then declare.
    ! Stating twice must therefore succeed and return the same.
    !----------------------------------------------------------------!

    call fine_solver % state(chain_statement(m), m, m % vertex_set(), &
         & m % num_vertices())
    xd = 0.0_dp
    call fine_solver % solve(rhs, xd, achieved)

    call report(all(abs(x - xd) < 1.0d-7), &
         & 'state is re-enterable: the second state solves the same', &
         & num_failures)

    !----------------------------------------------------------------!
    ! A COUNT COMPONENT NEEDS A DEFAULT. jacobi() and gmres() name no
    ! component, so every component of minimizer must have a value
    ! without being named - a count with no value is not a different
    ! state from a count of zero.
    !----------------------------------------------------------------!

    default_components: block
      type(jacobi) :: j
      type(gmres)  :: q
      j = jacobi()
      q = gmres()
      call report(j % num_unknowns .eq. 0 .and. &
           &      q % num_unknowns .eq. 0 .and. &
           &      j % num_residuals .eq. 0, &
           & 'an unstated solver counts no unknowns', &
           & num_failures)
    end block default_components

  end subroutine check_multigrid_gmres_equivalence

end program test_graph_multigrid
