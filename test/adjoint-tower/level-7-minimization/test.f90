!=====================================================================!
! ADJOINT TOWER . LEVEL 7 . THE SOLVER IS NEUTRAL
!
! The level answers one question: CAN ONE ORDINARY MINIMIZER GOVERN
! BOTH ORIENTATIONS. Two equations are supplied opaquely - the
! primal R(q,2) = 0 running Q -> Y, and the transposed
! A^T lambda = c running Y -> Q - and the SAME GMRES type solves
! both:
!
!      unknown Q, residual Y   ->  q      = [2, 4]
!      unknown Y, residual Q   ->  lambda = [-0.4, 0.6]
!
! The solver is never told which is which. It does not know the
! word adjoint; it knows an unknown domain, a residual domain, and
! an operation that answers on the latter. That is the whole claim
! of the rung, and there is no adjoint_solver, transpose_gmres or
! reverse_gmres anywhere.
!
! Gate A's theorem becomes numerical here. |Q| = |Y| = 2, so a
! right-hand side on the WRONG domain has exactly the right size,
! and only identity can reject it - the production minimizer does,
! by message, and that refusal is proved in check_refusals.sh.
! Likewise the non-symmetric A convicts a reversed orientation: had
! the adjoint been solved in the primal orientation it would answer
! [0.4, 0.2], and this test pins that it does not.
!
! The compatibility host is a five-vertex graph that is NEITHER Q
! NOR Y and contributes nothing to either answer - the legacy
! graph_operation face requires it, so it is supplied honestly and
! left alone. What Level 7 does NOT yet prove: that these two
! equations are one truth. They are independently supplied here, on
! purpose. Level 8 removes the duplication.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program adjoint_level_7

  use iso_fortran_env  , only : dp => REAL64
  use adjoint_assert   , only : report, verdict
  use adjoint_assert   , only : VAR_P, VAR_U, VAR_V
  use adjoint_assert   , only : TGT_R1, TGT_R2, TGT_F
  use graph_carrier    , only : counted_set, subset_set, member_set
  use graph_grammar    , only : graph_field
  use class_graph      , only : stored_graph
  use class_graph_field, only : field
  use class_graph_gmres, only : gmres
  use opaque_equation_fixture, only : opaque_primal, opaque_adjoint

  implicit none

  type(counted_set)               :: v, t, hv
  type(subset_set)                :: q_dom, y_dom
  type(stored_graph)              :: host
  type(opaque_primal)             :: primal_eq
  type(opaque_adjoint)            :: adjoint_eq
  type(gmres)                     :: primal_solver, adjoint_solver
  type(field)                     :: rhs_y, rhs_q
  class(member_set), allocatable  :: dom
  class(graph_field), allocatable :: q_sol, lam_sol
  real(dp), allocatable           :: gv(:), qv(:), lv(:)
  integer                         :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "adjoint tower . level 7 . solver neutrality"
  write(*,'(1x,a)') "============================================="

  v = counted_set('variables', 3)
  t = counted_set('targets'  , 3)
  q_dom = subset_set('state'   , v, [VAR_U, VAR_V])
  y_dom = subset_set('residual', t, [TGT_R1, TGT_R2])

  ! The compatibility host: five vertices, nobody's domain.
  host = stored_graph(5, tails=[1,2,3,4], heads=[2,3,4,5])

  primal_eq  = opaque_primal(q_dom, y_dom)
  adjoint_eq = opaque_adjoint(y_dom, q_dom)

  call check_host_is_nobodys_domain(nfail)
  call check_primal_solve(nfail)
  call check_adjoint_solve(nfail)
  call check_orientation_convicted(nfail)
  call check_one_solver_family(nfail)

  call verdict(nfail, "level 7")

contains

  !===================================================================!
  ! Three identities, all different: the solver host, the unknown
  ! domain and the residual domain. The legacy face demands a graph;
  ! neither Q nor Y pretends to be its vertex set to satisfy it.
  !===================================================================!

  subroutine check_host_is_nobodys_domain(nfail)

    integer, intent(inout) :: nfail

    hv = host % vertex_set()

    call report(.not. hv % same_as(q_dom) .and. &
         &      .not. hv % same_as(y_dom), &
         & "the host's vertex set is neither Q nor Y", nfail)
    call report(host % num_vertices() .ne. q_dom % size() .and. &
         &      host % num_vertices() .ne. y_dom % size(), &
         & "and it is not even their size: solver host, unknown " // &
         & "domain and residual domain are three different things", &
         & nfail)

  end subroutine check_host_is_nobodys_domain

  !===================================================================!
  ! The primal: unknown Q, residual Y, solved through the solver's
  ! own operation face - a right-hand side on Y in, a state on Q out.
  !===================================================================!

  subroutine check_primal_solve(nfail)

    integer, intent(inout) :: nfail

    call primal_solver % attach(primal_eq, host, q_dom)
    primal_solver % tolerance      = 1.0d-12
    primal_solver % max_iterations = 50

    call primal_solver % domain(host, dom)
    call report(dom % same_as(q_dom), &
         & "the primal solver answers on Q, by identity", nfail)
    call primal_eq % domain(host, dom)
    call report(dom % same_as(y_dom) .and. .not. dom % same_as(q_dom), &
         & "and the equation answers on Y - which is not Q", nfail)

    call primal_solver % constant(gv)
    call report(size(gv) .eq. 2 .and. &
         &      abs(gv(y_dom % local_index(TGT_R1)) + 8.0_dp) < 1.0d-12 &
         & .and. abs(gv(y_dom % local_index(TGT_R2)) + 22.0_dp) &
         &      < 1.0d-12, &
         & "R(0) = [-8, -22] on the residual rows", nfail)

    rhs_y = field('rhs', y_dom)
    call rhs_y % set_real_vector(-gv)
    call primal_solver % apply(host, [rhs_y], q_sol)

    call q_sol % domain(dom)
    call report(dom % same_as(q_dom), &
         & "the solution field lives on Q", nfail)

    call q_sol % get_real_vector(qv)
    call report(abs(qv(q_dom % local_index(VAR_U)) - 2.0_dp) < 1.0d-9, &
         & "u = 2, read through Q's enumeration", nfail)
    call report(abs(qv(q_dom % local_index(VAR_V)) - 4.0_dp) < 1.0d-9, &
         & "v = 4", nfail)

  end subroutine check_primal_solve

  !===================================================================!
  ! The adjoint: the orientation exchanged - unknown Y, residual Q -
  ! and the SAME solver type governing it.
  !===================================================================!

  subroutine check_adjoint_solve(nfail)

    integer, intent(inout) :: nfail

    call adjoint_solver % attach(adjoint_eq, host, y_dom)
    adjoint_solver % tolerance      = 1.0d-12
    adjoint_solver % max_iterations = 50

    call adjoint_solver % domain(host, dom)
    call report(dom % same_as(y_dom), &
         & "the adjoint solver answers on Y, by identity", nfail)
    call adjoint_eq % domain(host, dom)
    call report(dom % same_as(q_dom), &
         & "and its equation answers on Q: the orientation is " // &
         & "exchanged, end for end", nfail)

    call adjoint_solver % constant(gv)
    call report(abs(gv(q_dom % local_index(VAR_U)) + 1.0_dp) < 1.0d-12 &
         & .and. abs(gv(q_dom % local_index(VAR_V)) + 2.0_dp) &
         &      < 1.0d-12, &
         & "its constant is -c = [-1, -2], indexed by STATE slots", &
         & nfail)

    rhs_q = field('rhs', q_dom)
    call rhs_q % set_real_vector(-gv)
    call adjoint_solver % apply(host, [rhs_q], lam_sol)

    call lam_sol % domain(dom)
    call report(dom % same_as(y_dom), &
         & "the adjoint field lives on Y", nfail)

    call lam_sol % get_real_vector(lv)
    call report(abs(lv(y_dom % local_index(TGT_R1)) + 0.4_dp) < 1.0d-9, &
         & "lambda(r1) = -0.4", nfail)
    call report(abs(lv(y_dom % local_index(TGT_R2)) - 0.6_dp) < 1.0d-9, &
         & "lambda(r2) = 0.6", nfail)

  end subroutine check_adjoint_solve

  !===================================================================!
  ! The non-symmetric A doing its job: the answer is the TRANSPOSE
  ! solution, and demonstrably not the primal-orientation one.
  !===================================================================!

  subroutine check_orientation_convicted(nfail)

    integer, intent(inout) :: nfail

    call report(abs(lv(y_dom % local_index(TGT_R1)) - 0.4_dp) > 0.1_dp &
         & .and. abs(lv(y_dom % local_index(TGT_R2)) - 0.2_dp) &
         &      > 0.1_dp, &
         & "and it is NOT [0.4, 0.2] - what solving A lambda = c " // &
         & "would have given: a reversed transpose cannot hide " // &
         & "behind equal sizes", nfail)

  end subroutine check_orientation_convicted

  !===================================================================!
  ! One family, two orientations. The solvers differ only in the
  ! domains they were handed; neither knows the word adjoint.
  !===================================================================!

  subroutine check_one_solver_family(nfail)

    integer, intent(inout) :: nfail

    class(member_set), allocatable :: a, b

    call primal_solver % domain(host, a)
    call adjoint_solver % domain(host, b)

    call report(.not. a % same_as(b) .and. a % size() .eq. b % size(), &
         & "the two solvers answer on different domains of equal " // &
         & "size - the same machinery, told apart only by identity", &
         & nfail)

  end subroutine check_one_solver_family

end program adjoint_level_7
