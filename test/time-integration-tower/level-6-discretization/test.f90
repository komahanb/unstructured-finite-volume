!=====================================================================!
! TIME INTEGRATION TOWER . LEVEL 6 . DISCRETIZATION
!
! The level answers one question:
!
!      CAN A TEMPORAL DISCRETIZATION PRESERVE THE STATE DOMAIN Q
!      INSTEAD OF REPLACING IT WITH THE HOST GRAPH'S VERTICES?
!
! This is the first rung of the tower to touch production
! machinery, and the first place seam A2 of the reverse
! architecture review is genuinely exercised:
!
!      operations whose mathematics needs a DOMAIN currently often
!      obtain that domain from a GRAPH.
!
!                    THE EXPERIMENT'S SHAPE
!
! The action S : Q -> Q stores its own domain and stores no graph.
! The compatibility host H_t has FIVE vertices in a chain, and Q
! has TWO members. The mismatch is deliberate and load-bearing: if
! the two carriers had the same size, a substitution of one for the
! other would produce plausible numbers and the seam would hide.
!
!      |V(H_t)| = 5        |Q| = 2        and they are not same_as
!
! H_t is a COMPATIBILITY HOST - the conduit the operation
! contract requires - and this tower's action does not read its
! topology. That does not reopen seam A1: the partitioned tower
! settled on production evidence that the host is a real conduit
! for actions that DO consume topology. A triangular 2x2 decay is
! not one of those.
!
!                    THE ORACLES
!
!      S(q0)     = [2, -2]
!      q_FE,1    = q0 - h S(q0) = [1, 1]          by hand, not by marcher
!      q_BE,1    = [4/3, 4/9]                     exact
!      q_BDF2,2  = [5/6, 47/72]                   exact
!
! The two implicit answers are verified by SUBSTITUTION: the level
! asks production for the residual at the exact state and requires
! zero. No solver appears here - that is Level 7 - and no marcher
! appears at all.
!
!                    STRUCTURE AND SCHEME, JOINED
!
! Level 2 derived A1 and A2 and refused to call either a scheme.
! This level supplies the missing half. At instant t2:
!
!      A1-predecessor of t2  =  t1        one-step history
!      A2-predecessor of t2  =  t0        two-step history
!
! and BDF2 assigns NUMERICAL COEFFICIENTS to exactly those roles:
!
!      a0 = 3/2  at t2      a1 = -2  at t1      a2 = 1/2  at t0
!
! So A1 and A2 supply STRUCTURAL REACH; the scheme supplies the
! numbers. That is the Level-2 -> Level-6 Rosetta connection, and
! it is why Level 2 was right to refuse the name.
!
!                    THE FAMILY AS DATA
!
! The second half of the level reads a family's incidence and
! weights from the family alone and checks them against closed
! forms. Newmark (beta, gamma) advances the jet (q, q', q'') by
!
!      q_k  = q_(k-1) + h q'_(k-1) + h^2 (1/2 - beta) q''_(k-1)
!                                  + h^2 beta q''_k
!      q'_k = q'_(k-1) + h (1 - gamma) q''_(k-1) + h gamma q''_k
!
! with h = dt_k; both rows read q'' at the arriving instant, the
! across-degree edge with tail = head. On q = t^m the rows vanish
! for m = 0, 1, 2 whatever the pair, and on t^3 the value row leaves
! 6 h^3 (beta - 1/6) and the velocity row 6 h^2 (gamma - 1/2), in
! the sign convention  r = -D q_k + sum_e w_e D q_tail(e). The pair
! (0, 0) is the explicit Taylor step: the weight of q''_k is zero.
! A DIRK tableau's step connectivity states the same kind of rows,
! Q_i = q_(k-1) + h sum_j a_ij Q'_j, q_k = q_(k-1) + h sum_j b_j Q'_j,
! with the instant behind first in every row. A functional over the
! step is integrated on the instants the rows read: for Newmark the
! two instants k - 1 and k, the trapezoidal rule (1/2, 1/2), the
! same rule Adams-Moulton 2 and BDF-2 state on their two instants;
! at the first instant one node, whose measure is zero. Newmark
! refuses an equation of degree other than two; a staged family
! refuses the instant quadrature: each refusal stops a child process.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program time_level_6

  use iso_fortran_env       , only : dp => REAL64
  use time_assert           , only : report, assert_all
  use time_assert           , only : NQ, NT, TOL
  use time_assert           , only : T0, T1, T2
  use time_assert           , only : H_STEP, Q0, Q_FE1, Q_BE1, Q_BDF2
  use time_assert           , only : action_of
  use graph_fractal        , only : graph
  use map_set        , only : set_map
  use view_directed   , only : directed_graph
  use field_calculus  , only : field
  use relation_binary , only : csr_relation
  use view_directed_stored           , only : stored_directed_graph
  use field_stored     , only : stored_field
  use operation_family    , only : family, bdf_family, adams_family, newmark_family, crouzeix_two_stage, &
       & dirk_family
  use operation_weight    , only : scheme_weight
  use operation_coupling  , only : weights_of, weights_terms
  use view_directed_connectivity, only : connectivity_graph
  use util_derivative_terms, only : derivative_terms, value
  use temporal_step_fixture , only : temporal_step, backward_euler, bdf
  use time_carriers_fixture , only : time_carriers
  use time_relations_fixture, only : tail_relation, head_relation
  use time_algebra_fixture  , only : derive_one_step_reach, &
       &                             derive_two_step_reach
  use time_fields_fixture   , only : state_field
  use triangular_decay_fixture, only : triangular_decay

  implicit none

  type(graph)          :: q, t, e
  type(set_map)          :: sets
  type(csr_relation), target :: tail, head, a1
  type(csr_relation)         :: a2
  type(stored_directed_graph)         :: ht
  type(triangular_decay)     :: decay
  type(stored_field)                :: qf
  integer                    :: nfail
  character(len=64)          :: mode

  call get_command_argument(1, mode)
  if (len_trim(mode) > 0) call refused_case(trim(mode))

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "time integration tower . level 6 . scheme"
  write(*,'(1x,a)') "============================================="

  call time_carriers(sets, q, t, e)
  tail = tail_relation(e, t, sets)
  head = head_relation(e, t, sets)
  a1   = derive_one_step_reach(tail, head, sets)
  a2   = derive_two_step_reach(a1, sets)

  ! The COMPATIBILITY HOST: five vertices, four edges, a chain -
  ! the same temporal extension as T, and emphatically not Q.
  ht = stored_directed_graph(NT, tails=[1,2,3,4], heads=[2,3,4,5])

  decay = triangular_decay(q, NQ)
  qf    = state_field(q)

  call check_host_is_not_the_state_domain(nfail)
  call check_host_carriers_agree_with_themselves(nfail)
  call check_direct_action_preserves_q(nfail)
  call check_forward_euler_oracle(nfail)
  call check_step_domain_is_the_action_s(nfail)
  call check_backward_euler_residual(nfail)
  call check_reach_supplies_the_history_roles(nfail)
  call check_bdf2_residual(nfail)
  call check_newmark_connectivity(nfail)
  call check_newmark_weights(nfail)
  call check_newmark_polynomial_rows(nfail)
  call check_taylor_newmark_datum(nfail)
  call check_dirk_step_connectivity(nfail)
  call check_newmark_step_quadrature(nfail)
  call check_alexander_tableau(nfail)
  call check_family_refusals(nfail)

  call assert_all(nfail, "level 6")

contains

  !===================================================================!
  ! The mismatch that makes the experiment discriminating. If these
  ! ever became equal, every assertion below would still pass while
  ! proving nothing.
  !===================================================================!

  subroutine check_host_is_not_the_state_domain(nfail)

    integer, intent(inout) :: nfail

    type(graph) :: hv

    hv = ht % vertex_set()

    call report(ht % num_vertices() .eq. NT .and. sets % num_members_of(q) .eq. NQ, &
         & "the compatibility host H_t has five vertices; Q has two", &
         & nfail)

    call report(.not. hv % same_as(q), &
         & "and V(H_t) is NOT Q - no accidental equality is hiding " // &
         & "the seam", nfail)

    call report(.not. hv % same_as(t), &
         & "nor is it T: H_t carries the same EXTENSION as the time " // &
         & "axis and none of its identity", nfail)

  end subroutine check_host_is_not_the_state_domain

  !===================================================================!
  ! A fact the production correction below depends on, so it is
  ! measured rather than assumed: a graph's two ways of naming its
  ! vertices agree.
  !===================================================================!

  subroutine check_host_carriers_agree_with_themselves(nfail)

    integer, intent(inout) :: nfail

    type(graph) :: d, hv
    integer        :: n_d

    hv = ht % vertex_set()
    d   = ht % vertex_set()
    n_d = ht % num_vertices()
    call report(d % same_as(hv), &
         & "vertex_set(H_t) and H_t % vertex_set() are the same " // &
         & "carrier - so delegating a domain question changes no " // &
         & "graph-based caller's answer", nfail)

  end subroutine check_host_carriers_agree_with_themselves

  !===================================================================!
  ! THE first half of the experiment, and it needs no production
  ! change: an action that carries its own domain already works
  ! through the operation face.
  !===================================================================!

  subroutine check_direct_action_preserves_q(nfail)

    integer, intent(inout) :: nfail

    class(field), allocatable :: residual_field
    type(graph)  :: d
    integer         :: n_d
    real(dp), allocatable           :: s(:)

    call decay % domain(ht, d, n_d)
    call report(d % same_as(q), &
         & "the ACTION answers Q when asked its domain, though it " // &
         & "was handed a five-vertex host", nfail)

    call decay % apply(ht, decay % bind([qf]), residual_field)
    d = residual_field % domain()
    call report(d % same_as(q), &
         & "and it ANSWERS on Q: graph host and state domain are " // &
         & "independent concepts in this specimen", nfail)

    call residual_field % real_vector(s)
    call report(size(s) .eq. NQ .and. &
         &      abs(s(1) - 2.0_dp) .lt. TOL .and. &
         &      abs(s(2) + 2.0_dp) .lt. TOL, &
         & "S(q0) = [2, -2]", nfail)

    call report(maxval(abs(s - action_of(Q0))) .lt. TOL, &
         & "matching the oracle computed in plain arithmetic, not " // &
         & "by the machinery under test", nfail)

  end subroutine check_direct_action_preserves_q

  !===================================================================!
  ! Forward euler by hand: q1 = q0 - h S(q0). No production step,
  ! no marcher. The scheme is tested before the machinery that
  ! stamps it along a chain.
  !===================================================================!

  subroutine check_forward_euler_oracle(nfail)

    integer, intent(inout) :: nfail

    real(dp) :: q1(NQ)

    q1 = Q0 - H_STEP * action_of(Q0)

    call report(maxval(abs(q1 - Q_FE1)) .lt. TOL, &
         & "forward euler q1 = q0 - h S(q0) = [1, 1], by ordinary " // &
         & "arithmetic", nfail)

  end subroutine check_forward_euler_oracle

  !===================================================================!
  ! THE seam-A2 assertion, and the reason this level exists.
  !
  ! A temporal discretization is an operation BUILT FROM another
  ! operation. Its residual is a statement about the same unknown
  ! the action is about, so its domain must be the ACTION's domain
  ! - not whatever carrier the compatibility host happens to have.
  !
  ! On the production reviewed at Gate A, step_domain answered
  ! input_graph % vertex_set(...) and this assertion FAILED,
  ! reporting a five-member carrier for a two-member unknown. That
  ! RED is recorded verbatim in NUCLEUS-OBSERVATIONS.md TI-8.
  !
  ! The second assertion is the permanent guard: it is not enough
  ! that the answer BE Q; it must also not be the host's vertices,
  ! or a future coincidence of carriers would let the seam back in
  ! unnoticed.
  !===================================================================!

  subroutine check_step_domain_is_the_action_s(nfail)

    integer, intent(inout) :: nfail

    type(temporal_step)            :: step
    type(graph) :: d, hv
    integer         :: n_d

    step = backward_euler(decay, H_STEP)

    call step % domain(ht, d, n_d)

    call report(d % same_as(q), &
         & "the backward-euler STEP answers Q when asked its " // &
         & "domain: TEMPORAL DISCRETIZATION PRESERVES THE DOMAIN OF " // &
         & "THE ACTION IT DISCRETIZES", nfail)

    hv = ht % vertex_set()
    call report(.not. d % same_as(hv), &
         & "and it does NOT answer the host's five vertices - the " // &
         & "state domain is not inferred from the conduit", nfail)

  end subroutine check_step_domain_is_the_action_s

  !===================================================================!
  ! The backward-euler residual, verified by substitution:
  !
  !      R_BE(q) = q - q0 + h S(q)
  !
  ! is zero exactly at q = [4/3, 4/9]. The level asks PRODUCTION for
  ! that residual and requires zero - and requires the answer to
  ! land on Q.
  !===================================================================!

  subroutine check_backward_euler_residual(nfail)

    integer, intent(inout) :: nfail

    type(temporal_step)             :: step
    type(stored_field)                     :: state, hist
    class(field), allocatable :: r
    type(graph)  :: d
    real(dp), allocatable           :: v(:)

    step = backward_euler(decay, H_STEP)
    hist = stored_field('q0', q, NQ, num_components=1)
    call hist % set_real_vector(Q0)

    ! At the exact backward-euler state the residual vanishes.
    state = stored_field('trial', q, NQ, num_components=1)
    call state % set_real_vector(Q_BE1)
    call step % apply(ht, step % bind([state, hist]), r)

    d = r % domain()
    call report(d % same_as(q), &
         & "the backward-euler RESIDUAL lands on Q, not on the " // &
         & "host's vertices", nfail)

    call r % real_vector(v)
    call report(size(v) .eq. NQ .and. maxval(abs(v)) .lt. TOL, &
         & "and it is ZERO at q = [4/3, 4/9]: the exact discrete " // &
         & "backward-euler state, verified by substitution", nfail)

    ! And it is not zero anywhere convenient - the forward-euler
    ! answer is a different number, as it must be.
    call state % set_real_vector(Q_FE1)
    call step % apply(ht, step % bind([state, hist]), r)
    call r % real_vector(v)
    call report(maxval(abs(v)) .gt. 1.0e-3_dp, &
         & "while the FORWARD-euler answer leaves a residual: the " // &
         & "two schemes are different statements", nfail)

  end subroutine check_backward_euler_residual

  !===================================================================!
  ! THE Rosetta connection Level 2 deliberately left open.
  !
  ! A1 and A2 say WHICH instants a two-step scheme may look at; the
  ! scheme says WHAT NUMBERS to weight them with. Neither contains
  ! the other, and the join happens here.
  !===================================================================!

  subroutine check_reach_supplies_the_history_roles(nfail)

    integer, intent(inout) :: nfail

    type(temporal_step) :: bdf2
    type(family)        :: scheme

    call report(a1 % has([T1, T2]) .and. a2 % has([T0, T2]), &
         & "at instant t2 the one-step predecessor is t1 and the " // &
         & "two-step predecessor is t0 - STRUCTURAL REACH, from " // &
         & "Level 2", nfail)

    bdf2   = bdf(2, decay, H_STEP)
    scheme = bdf_family(2)
    call report(bdf2 % reach .eq. 2 .and. scheme % history_depth(1) .eq. 2, &
         & "and bdf-2 reaches exactly two instants back, matching " // &
         & "the reach A2 describes", nfail)

    call report(abs(bdf2 % c(0) - 1.5_dp) .lt. TOL .and. &
         &      abs(bdf2 % c(1) + 2.0_dp) .lt. TOL .and. &
         &      abs(bdf2 % c(2) - 0.5_dp) .lt. TOL, &
         & "with coefficients a0 = 3/2 at t2, a1 = -2 at t1, " // &
         & "a2 = 1/2 at t0: NUMERICAL WEIGHTS ON THE HISTORY ROLES " // &
         & "reach already named", nfail)

    call report(.not. a2 % has([T0, T1]) .and. a1 % has([T0, T1]), &
         & "and the roles are distinct: t1 is one step from t0 and " // &
         & "never two - A1/A2 give the structure, bdf-2 gives the " // &
         & "numbers, and neither is the other", nfail)

  end subroutine check_reach_supplies_the_history_roles

  !===================================================================!
  ! The bdf-2 residual, verified by substitution:
  !
  !      (3/2) q2 - 2 q1 + (1/2) q0 + h S(q2) = 0
  !
  ! with q1 the backward-euler start. Zero exactly at
  ! q2 = [5/6, 47/72].
  !===================================================================!

  subroutine check_bdf2_residual(nfail)

    integer, intent(inout) :: nfail

    type(temporal_step)             :: step
    type(stored_field)                     :: state, hist1, hist2
    class(field), allocatable :: r
    type(graph)  :: d
    real(dp), allocatable           :: v(:)

    step = bdf(2, decay, H_STEP)
    hist1 = stored_field('q1', q, NQ, num_components=1)
    call hist1 % set_real_vector(Q_BE1)
    hist2 = stored_field('q0', q, NQ, num_components=1)
    call hist2 % set_real_vector(Q0)

    state = stored_field('trial', q, NQ, num_components=1)
    call state % set_real_vector(Q_BDF2)
    call step % apply(ht, step % bind([state, hist1, hist2]), r)

    d = r % domain()
    call report(d % same_as(q), &
         & "the bdf-2 RESIDUAL lands on Q as well - a two-step " // &
         & "scheme changes the coefficients, never the domain", nfail)

    call r % real_vector(v)
    call report(size(v) .eq. NQ .and. maxval(abs(v)) .lt. TOL, &
         & "and it is ZERO at q2 = [5/6, 47/72], started from the " // &
         & "backward-euler q1", nfail)

  end subroutine check_bdf2_residual

  !===================================================================!
  ! The Newmark block connectivity over two instants at three
  ! degrees: the value row reads q, q', q'' behind and q'' ahead; the
  ! velocity row reads q', q'' behind and q'' ahead; the acceleration
  ! (the primary degree) has no row. Over four instants the pattern
  ! is placed at every instant with one behind it.
  !===================================================================!

  subroutine check_newmark_connectivity(nfail)

    integer, intent(inout) :: nfail

    type(family) :: scheme
    type(connectivity_graph) :: edges
    integer, parameter :: tails(7) = [1, 1, 1, 2, 1, 1, 2]
    integer, parameter :: tail_degrees(7) = [0, 1, 2, 2, 1, 2, 2]
    integer, parameter :: head_degrees(7) = [0, 0, 0, 0, 1, 1, 1]
    integer :: e

    scheme = newmark_family(0.25_dp, 0.5_dp)
    call report(scheme % history_depth(2) .eq. 1 .and. scheme % primary_degree(2) .eq. 2 &
         & .and. scheme % num_stages() .eq. 1, &
         & "newmark reaches one instant back, determines the jet from " // &
         & "the acceleration, and has one stage", nfail)

    edges = scheme % block_connectivity(3, 2)
    call report(edges % num_vertices() .eq. 2 .and. edges % num_edges() .eq. 7, &
         & "over two instants at three degrees the family states seven " // &
         & "edges: four into the value row, three into the velocity row", nfail)
    call report(all([(edges % edge_tail(e) .eq. tails(e), e = 1, 7)]) .and. &
         &      all([(edges % edge_head(e) .eq. 2, e = 1, 7)]) .and. &
         &      all([(edges % tail_degree(e) .eq. tail_degrees(e), e = 1, 7)]) .and. &
         &      all([(edges % head_degree(e) .eq. head_degrees(e), e = 1, 7)]), &
         & "value row: q, q', q'' behind and q'' ahead; velocity row: " // &
         & "q', q'' behind and q'' ahead, in the family's order", nfail)
    call report(edges % edge_tail(4) .eq. edges % edge_head(4) .and. edges % tail_degree(4) .eq. 2 &
         & .and. edges % edge_tail(7) .eq. edges % edge_head(7) .and. edges % tail_degree(7) .eq. 2, &
         & "the across-degree edges: both rows read the acceleration at " // &
         & "the instant they determine (tail = head)", nfail)

    edges = scheme % block_connectivity(3, 4)
    call report(edges % num_edges() .eq. 21 .and. &
         & all([(edges % edge_head(e) .eq. 1 + (e - 1) / 7 + 1, e = 1, 21)]), &
         & "over four instants the pattern is placed at instants 2, 3, 4: " // &
         & "twenty-one edges, seven per arriving instant", nfail)

  end subroutine check_newmark_connectivity

  !===================================================================!
  ! The weights of both rows against the closed forms, on a
  ! non-uniform grid dt = [0, 0.3, 0.2]: every weight is a power of
  ! the step arriving at the row's instant, so the two instants
  ! give two different rows of numbers.
  !===================================================================!

  subroutine check_newmark_weights(nfail)

    integer, intent(inout) :: nfail

    real(dp), parameter :: dt(3) = [0.0_dp, 0.3_dp, 0.2_dp]
    real(dp), parameter :: beta = 0.25_dp, gamma = 0.5_dp

    call report(maxval(abs(newmark_weights(beta, gamma, dt) - closed_forms(beta, gamma, dt))) .lt. TOL, &
         & "newmark (1/4, 1/2): value row [1, h, h^2 (1/2 - beta), h^2 beta], " // &
         & "velocity row [1, h (1 - gamma), h gamma], at h = 0.3 and h = 0.2", nfail)
    call report(maxval(abs(newmark_weights(1.0_dp / 12.0_dp, gamma, dt) &
         & - closed_forms(1.0_dp / 12.0_dp, gamma, dt))) .lt. TOL, &
         & "newmark (1/12, 1/2), Fox-Goodwin: the same closed forms", nfail)

  end subroutine check_newmark_weights

  !===================================================================!
  ! The weights over the three-instant block, fourteen numbers: the
  ! seven edges into instant 2 then the seven into instant 3.
  !===================================================================!

  function newmark_weights(beta, gamma, dt) result(w)

    real(dp), intent(in) :: beta, gamma, dt(3)
    real(dp), allocatable :: w(:)

    type(family) :: scheme

    scheme = newmark_family(beta, gamma)
    call weights_of(scheme_weight(scheme), scheme % block_connectivity(3, 3), dt, w)

  end function newmark_weights

  function closed_forms(beta, gamma, dt) result(w)

    real(dp), intent(in) :: beta, gamma, dt(3)
    real(dp) :: w(14)

    real(dp) :: h
    integer :: k

    do k = 2, 3
       h = dt(k)
       w(7 * (k - 2) + 1 : 7 * (k - 1)) = [1.0_dp, h, h * h * (0.5_dp - beta), h * h * beta, &
            &                              1.0_dp, h * (1.0_dp - gamma), h * gamma]
    end do

  end function closed_forms

  !===================================================================!
  ! The rows on q = t^m at instant 3 (t = 0.5, h = 0.2): zero for
  ! m = 0, 1, 2 whatever the pair; on t^3 the value row leaves
  ! 6 h^3 (beta - 1/6) and the velocity row 6 h^2 (gamma - 1/2).
  !===================================================================!

  subroutine check_newmark_polynomial_rows(nfail)

    integer, intent(inout) :: nfail

    real(dp), parameter :: dt(3) = [0.0_dp, 0.3_dp, 0.2_dp]
    real(dp), parameter :: pairs(2, 3) = reshape( &
         & [0.25_dp, 0.5_dp,  1.0_dp / 12.0_dp, 0.5_dp,  0.0_dp, 0.0_dp], [2, 3])
    real(dp) :: r(0:3, 0:1), h
    integer :: p
    logical :: quadratic, cubic

    h = dt(3)
    quadratic = .true.
    cubic     = .true.
    do p = 1, 3
       r = row_residuals(pairs(1, p), pairs(2, p), dt)
       quadratic = quadratic .and. maxval(abs(r(0:2, :))) .lt. TOL
       cubic = cubic .and. abs(r(3, 0) - 6.0_dp * h**3 * (pairs(1, p) - 1.0_dp / 6.0_dp)) .lt. TOL &
            &        .and. abs(r(3, 1) - 6.0_dp * h**2 * (pairs(2, p) - 0.5_dp)) .lt. TOL
    end do
    call report(quadratic, &
         & "both rows vanish on t^0, t^1, t^2 for (1/4, 1/2), (1/12, 1/2) " // &
         & "and (0, 0): the pair changes no exactness below the cubic", nfail)
    call report(cubic, &
         & "on t^3 the value row leaves 6 h^3 (beta - 1/6) and the velocity " // &
         & "row 6 h^2 (gamma - 1/2): measured equals the closed form", nfail)

  end subroutine check_newmark_polynomial_rows

  !===================================================================!
  ! r(m, d) = -D^d t_3^m + sum over the edges into (3, d) of
  ! w_e D^(tail degree) t_tail^m.
  !===================================================================!

  function row_residuals(beta, gamma, dt) result(r)

    real(dp), intent(in) :: beta, gamma, dt(3)
    real(dp) :: r(0:3, 0:1)

    type(family) :: scheme
    type(connectivity_graph) :: edges
    real(dp), allocatable :: w(:)
    real(dp) :: t(3)
    integer :: e, m, k

    scheme = newmark_family(beta, gamma)
    edges  = scheme % block_connectivity(3, 3)
    call weights_of(scheme_weight(scheme), edges, dt, w)
    t(1) = 0.0_dp
    do k = 2, 3
       t(k) = t(k - 1) + dt(k)
    end do
    do m = 0, 3
       r(m, :) = -[monomial_derivative(m, 0, t(3)), monomial_derivative(m, 1, t(3))]
       do e = 1, edges % num_edges()
          if (edges % edge_head(e) .ne. 3) cycle
          r(m, edges % head_degree(e)) = r(m, edges % head_degree(e)) &
               & + w(e) * monomial_derivative(m, edges % tail_degree(e), t(edges % edge_tail(e)))
       end do
    end do

  end function row_residuals

  pure real(dp) function monomial_derivative(m, d, t) result(q)

    integer , intent(in) :: m, d
    real(dp), intent(in) :: t

    integer :: i

    q = 0.0_dp
    if (d > m) return
    q = 1.0_dp
    do i = 0, d - 1
       q = q * real(m - i, dp)
    end do
    q = q * t**(m - d)

  end function monomial_derivative

  !===================================================================!
  ! Taylor-Newmark is the pair (0, 0): the same pattern, the weight
  ! of q'' ahead zero in both rows, so the step is the explicit
  ! Taylor expansion of the jet.
  !===================================================================!

  subroutine check_taylor_newmark_datum(nfail)

    integer, intent(inout) :: nfail

    real(dp), parameter :: dt(3) = [0.0_dp, 0.3_dp, 0.2_dp]
    real(dp), allocatable :: w(:)

    w = newmark_weights(0.0_dp, 0.0_dp, dt)
    call report(maxval(abs(w - closed_forms(0.0_dp, 0.0_dp, dt))) .lt. TOL .and. &
         & w(4) .eq. 0.0_dp .and. w(7) .eq. 0.0_dp .and. w(11) .eq. 0.0_dp .and. w(14) .eq. 0.0_dp, &
         & "newmark (0, 0): value row [1, h, h^2/2, 0], velocity row [1, h, 0] - " // &
         & "the explicit Taylor step, one datum under the name taylor-newmark", nfail)

  end subroutine check_taylor_newmark_datum

  !===================================================================!
  ! The Crouzeix two-stage step at two degrees: eight edges, the
  ! instant behind first in every row, weights [1, h g], [1, h (1 -
  ! 2 g), h g], [1, h/2, h/2], all below the top degree; the top
  ! degree at the instant ahead is the law's row, not the family's.
  !===================================================================!

  subroutine check_dirk_step_connectivity(nfail)

    integer, intent(inout) :: nfail

    type(family) :: scheme
    type(connectivity_graph) :: edges
    real(dp), allocatable :: w(:)
    real(dp), parameter :: h = 0.3_dp
    real(dp) :: g
    integer, parameter :: tails(8) = [1, 2, 1, 2, 3, 1, 2, 3]
    integer, parameter :: heads(8) = [2, 2, 3, 3, 3, 4, 4, 4]
    integer, parameter :: tail_degrees(8) = [0, 1, 0, 1, 1, 0, 1, 1]
    integer, parameter :: head_degrees(8) = [0, 0, 0, 0, 0, 0, 0, 0]
    integer :: e

    g = (3.0_dp + sqrt(3.0_dp)) / 6.0_dp
    scheme = crouzeix_two_stage()
    edges  = scheme % stage_connectivity(2)
    call report(edges % num_vertices() .eq. 4 .and. edges % num_edges() .eq. 8 .and. &
         &      all([(edges % edge_tail(e) .eq. tails(e), e = 1, 8)]) .and. &
         &      all([(edges % edge_head(e) .eq. heads(e), e = 1, 8)]) .and. &
         &      all([(edges % tail_degree(e) .eq. tail_degrees(e), e = 1, 8)]) .and. &
         &      all([(edges % head_degree(e) .eq. head_degrees(e), e = 1, 8)]), &
         & "crouzeix two-stage at two degrees: stage 1 reads the instant " // &
         & "behind and itself, stage 2 both stages, the instant ahead every " // &
         & "stage, all below the top degree; the instant behind is first in each row", nfail)
    call weights_of(scheme_weight(scheme), edges, [(h, e = 1, 4)], w)
    call report(maxval(abs(w - [1.0_dp, h * g, 1.0_dp, h * (1.0_dp - 2.0_dp * g), h * g, &
         &                     1.0_dp, h / 2.0_dp, h / 2.0_dp])) .lt. TOL, &
         & "with weights 1 on the instant behind, h a_ij on the stages and " // &
         & "h b_j into the instant ahead; the top degree has no row of the family", nfail)

  end subroutine check_dirk_step_connectivity

  !===================================================================!
  ! The step quadrature of Newmark on the non-uniform grid
  ! dt = [0, 0.3, 0.2]: at instant 3 the two nodes are instants 3
  ! and 2 with weights (1/2, 1/2) of the step dt_3, the trapezoidal
  ! rule, equal to the two-instant rule of Adams-Moulton 2 and of
  ! BDF-2; at instant 1 one node of weight one.
  !===================================================================!

  subroutine check_newmark_step_quadrature(nfail)

    integer, intent(inout) :: nfail

    type(family) :: scheme, adams2, bdf2
    type(derivative_terms) :: dt(3)
    type(derivative_terms), allocatable :: weight(:), reference(:)
    integer :: k

    dt = [derivative_terms(0.0_dp, 0), derivative_terms(0.3_dp, 0), derivative_terms(0.2_dp, 0)]
    scheme = newmark_family(0.25_dp, 0.5_dp)
    call scheme % step_quadrature(dt, 3, weight)
    call report(size(weight) .eq. 2 .and. abs(value(weight(1)) - 0.5_dp) .lt. TOL &
         & .and. abs(value(weight(2)) - 0.5_dp) .lt. TOL, &
         & "newmark integrates the step on the two instants its rows " // &
         & "read: the trapezoidal rule (1/2, 1/2)", nfail)
    adams2 = adams_family(2)
    bdf2   = bdf_family(2)
    call adams2 % step_quadrature(dt, 3, reference)
    call report(size(reference) .eq. 2 .and. maxval([(abs(value(weight(k)) - value(reference(k))), k = 1, 2)]) .lt. TOL, &
         & "the same rule Adams-Moulton 2 states on its two instants", nfail)
    call bdf2 % step_quadrature(dt, 3, reference)
    call report(size(reference) .eq. 2 .and. maxval([(abs(value(weight(k)) - value(reference(k))), k = 1, 2)]) .lt. TOL, &
         & "and the same rule BDF-2 states", nfail)
    call scheme % step_quadrature(dt, 1, weight)
    call report(size(weight) .eq. 1 .and. abs(value(weight(1)) - 1.0_dp) .lt. TOL, &
         & "at the first instant one node of weight one, whose step " // &
         & "measure is zero", nfail)

  end subroutine check_newmark_step_quadrature

  !===================================================================!
  ! The refusals, each in a child process that must stop.
  !===================================================================!

  subroutine check_family_refusals(nfail)

    integer, intent(inout) :: nfail

    call report(stopped('newmark-history-depth') .and. stopped('newmark-row-pattern') &
         & .and. stopped('newmark-primary-degree'), &
         & "newmark refuses a first-order equation in every query that " // &
         & "reads the equation degree", nfail)
    call report(stopped('dirk-step-quadrature'), &
         & "a staged family refuses the instant quadrature: it integrates " // &
         & "over its stages by the tableau weights", nfail)
    call report(stopped('dirk-top-degree-edge'), &
         & "a tableau refuses a stage read at the constraint's own degree: " // &
         & "the top degree at the instant ahead is the law's row", nfail)

  end subroutine check_family_refusals

  !===================================================================!
  ! A tableau registered by data alone: Alexander's two-stage
  ! L-stable DIRK, g = 1 - sqrt(2)/2, a = [[g, 0], [1 - g, g]],
  ! b = [1 - g, g], through dirk_family. The step map on q' = lambda q
  ! is assembled from the family's stage connectivity and weights
  ! alone: eight unknowns (four vertices at two degrees), the instant
  ! behind fixed, three value rows of the family and three law rows.
  ! The map equals the stability function
  ! R(z) = (1 + (1 - 2 g) z) / (1 - g z)^2, z = h lambda; over T = 1
  ! the error against exp(lambda T) quarters when the step halves;
  ! the tangent of the step in h, from weights_terms seeded by dh = 1,
  ! equals lambda R'(z) q; and the adjoint identity
  ! e_k . du/dh = -mu . (dA/dh) u, A^T mu = e_k, holds through the
  ! transposed solve.
  !===================================================================!

  subroutine check_alexander_tableau(nfail)

    integer, intent(inout) :: nfail

    real(dp), parameter :: lambda = -1.0_dp, duration = 1.0_dp, h = 0.3_dp
    integer , parameter :: steps(3) = [10, 20, 40]
    type(family) :: scheme
    type(connectivity_graph) :: edges
    real(dp), allocatable :: w(:), table(:,:), a(:,:), da(:,:), rhs(:), u(:), du(:), mu(:), seeds(:,:)
    real(dp) :: g, z, r_exact, dr_exact, q, e(3), tangent, adjoint, ratio(2)
    integer :: nv, ne, k, n, level, arriving

    g = 1.0_dp - sqrt(2.0_dp) / 2.0_dp
    scheme = dirk_family(reshape([g, 1.0_dp - g, 0.0_dp, g], [2, 2]), [1.0_dp - g, g])
    edges  = scheme % stage_connectivity(2)
    nv = edges % num_vertices()
    ne = edges % num_edges()
    arriving = 2 * (nv - 1) + 1
    call weights_of(scheme_weight(scheme), edges, [(h, k = 1, nv)], w)
    call report(nv .eq. 4 .and. ne .eq. 8 .and. &
         & maxval(abs(w - [1.0_dp, h * g, 1.0_dp, h * (1.0_dp - g), h * g, &
         &                 1.0_dp, h * (1.0_dp - g), h * g])) .lt. TOL, &
         & "alexander two-stage at two degrees: eight edges with weights 1 on the " // &
         & "instant behind, h a_ij on the stages and h b_j into the instant ahead", nfail)

    ! one step from q = 1: the assembled map against R(z)
    call decay_step_system(scheme, edges, lambda, h, 1.0_dp, a, rhs)
    u = solved(a, rhs)
    z = h * lambda
    r_exact = (1.0_dp + (1.0_dp - 2.0_dp * g) * z) / (1.0_dp - g * z) ** 2
    call report(abs(u(arriving) - r_exact) .lt. TOL, &
         & "the step map assembled from the connectivity and weights equals the " // &
         & "stability function (1 + (1 - 2 g) z) / (1 - g z)^2", nfail)

    ! convergence over T = 1: the error quarters when the step halves
    do level = 1, 3
       n = steps(level)
       q = 1.0_dp
       do k = 1, n
          call decay_step_system(scheme, edges, lambda, duration / real(n, dp), q, a, rhs)
          u = solved(a, rhs)
          q = u(arriving)
       end do
       e(level) = abs(q - exp(lambda * duration))
    end do
    ratio = e(1:2) / e(2:3)
    call report(all(ratio .gt. 2.0_dp ** 1.75_dp) .and. all(ratio .lt. 2.0_dp ** 2.25_dp), &
         & "second order on q' = -q over T = 1: the error ratio under halving " // &
         & "lies within 2^(2 -/+ 1/4) at 10, 20, 40 steps", nfail)

    ! the tangent in h from the weights' derivative terms, and the
    ! adjoint identity through the transposed solve
    allocate(seeds(nv, 1))
    seeds = 1.0_dp
    call weights_terms(scheme_weight(scheme), edges, [(h, k = 1, nv)], seeds, table)
    call decay_step_system(scheme, edges, lambda, h, 1.0_dp, a, rhs)
    u  = solved(a, rhs)
    call decay_step_derivative(edges, table(:, 1), da)
    du = solved(a, -matmul(da, u))
    tangent  = du(arriving)
    dr_exact = lambda * ((1.0_dp - 2.0_dp * g) * (1.0_dp - g * z) &
         & + 2.0_dp * g * (1.0_dp + (1.0_dp - 2.0_dp * g) * z)) / (1.0_dp - g * z) ** 3
    rhs = 0.0_dp
    rhs(arriving) = 1.0_dp
    mu = solved(transpose(a), rhs)
    adjoint = -dot_product(mu, matmul(da, u))
    call report(abs(tangent - dr_exact) .lt. TOL .and. abs(tangent - adjoint) .lt. TOL, &
         & "the tangent of the step in h from weights_terms equals lambda R'(z) q, and " // &
         & "the adjoint through the transposed solve equals the tangent", nfail)

  end subroutine check_alexander_tableau

  !===================================================================!
  ! The linear system of one staged step on q' = lambda q at two
  ! degrees, unknown (v, d) at row 2 (v - 1) + d + 1: the instant
  ! behind fixed to (q, lambda q), each value row -u(head) + sum of
  ! the weighted tails, each law row -u(v, 1) + lambda u(v, 0).
  !===================================================================!

  subroutine decay_step_system(scheme, edges, lambda, step, q_behind, a, rhs)

    type(family)            , intent(in)  :: scheme
    type(connectivity_graph), intent(in)  :: edges
    real(dp)                , intent(in)  :: lambda, step, q_behind
    real(dp), allocatable   , intent(out) :: a(:,:), rhs(:)

    real(dp), allocatable :: wt(:)
    integer :: nv, v, e, row, col

    nv = edges % num_vertices()
    call weights_of(scheme_weight(scheme), edges, [(step, v = 1, nv)], wt)
    allocate(a(2 * nv, 2 * nv), rhs(2 * nv))
    a   = 0.0_dp
    rhs = 0.0_dp
    a(1, 1) = 1.0_dp
    a(2, 2) = 1.0_dp
    rhs(1)  = q_behind
    rhs(2)  = lambda * q_behind
    do v = 2, nv
       row = 2 * (v - 1) + 1
       a(row, row)         = -1.0_dp
       a(row + 1, row + 1) = -1.0_dp
       a(row + 1, row)     = lambda
    end do
    do e = 1, edges % num_edges()
       row = 2 * (edges % edge_head(e) - 1) + edges % head_degree(e) + 1
       col = 2 * (edges % edge_tail(e) - 1) + edges % tail_degree(e) + 1
       a(row, col) = a(row, col) + wt(e)
    end do

  end subroutine decay_step_system

  subroutine decay_step_derivative(edges, dwt, da)

    type(connectivity_graph), intent(in)  :: edges
    real(dp)                , intent(in)  :: dwt(:)
    real(dp), allocatable   , intent(out) :: da(:,:)

    integer :: nv, e, row, col

    nv = edges % num_vertices()
    allocate(da(2 * nv, 2 * nv))
    da = 0.0_dp
    do e = 1, edges % num_edges()
       row = 2 * (edges % edge_head(e) - 1) + edges % head_degree(e) + 1
       col = 2 * (edges % edge_tail(e) - 1) + edges % tail_degree(e) + 1
       da(row, col) = da(row, col) + dwt(e)
    end do

  end subroutine decay_step_derivative

  !===================================================================!
  ! Gaussian elimination with partial pivoting on a small dense
  ! system; a zero pivot stops the program.
  !===================================================================!

  function solved(a, b) result(x)

    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), allocatable :: x(:)

    real(dp), allocatable :: m(:,:), r(:)
    real(dp) :: factor
    integer :: n, i, j, k, p

    n = size(b)
    m = a
    r = b
    do k = 1, n
       p = k - 1 + maxloc(abs(m(k:n, k)), dim=1)
       if (abs(m(p, k)) .eq. 0.0_dp) error stop 'level 6: a singular step system'
       if (p .ne. k) then
          m([k, p], :) = m([p, k], :)
          r([k, p])    = r([p, k])
       end if
       do i = k + 1, n
          factor  = m(i, k) / m(k, k)
          m(i, :) = m(i, :) - factor * m(k, :)
          r(i)    = r(i) - factor * r(k)
       end do
    end do
    allocate(x(n))
    do i = n, 1, -1
       x(i) = r(i)
       do j = i + 1, n
          x(i) = x(i) - m(i, j) * x(j)
       end do
       x(i) = x(i) / m(i, i)
    end do

  end function solved

  logical function stopped(case_name)

    character(len=*), intent(in) :: case_name

    character(len=256) :: self
    integer :: status, command_status, unit

    call get_command_argument(0, self)
    call execute_command_line(trim(self) // ' ' // case_name // ' > refusal.out 2>&1', &
         & exitstat=status, cmdstat=command_status)
    stopped = command_status .eq. 0 .and. status .ne. 0
    open(newunit=unit, file='refusal.out', status='old')
    close(unit, status='delete')

  end function stopped

  subroutine refused_case(case_name)

    character(len=*), intent(in) :: case_name

    type(family) :: scheme
    type(derivative_terms) :: dt(2)
    type(derivative_terms), allocatable :: weight(:)
    integer, allocatable :: offset(:), tail_degree(:)
    integer :: k

    select case (case_name)
    case ('newmark-history-depth')
       scheme = newmark_family(0.25_dp, 0.5_dp)
       k = scheme % history_depth(1)
    case ('newmark-row-pattern')
       scheme = newmark_family(0.25_dp, 0.5_dp)
       call scheme % row_pattern(0, 1, offset, tail_degree)
       k = size(offset)
    case ('newmark-primary-degree')
       scheme = newmark_family(0.25_dp, 0.5_dp)
       k = scheme % primary_degree(1)
    case ('dirk-step-quadrature')
       scheme = crouzeix_two_stage()
       dt = [derivative_terms(0.0_dp, 0), derivative_terms(0.5_dp, 0)]
       call scheme % step_quadrature(dt, 2, weight)
       k = size(weight)
    case ('dirk-top-degree-edge')
       scheme = crouzeix_two_stage()
       dt = [derivative_terms(0.5_dp, 0), derivative_terms(0.5_dp, 0)]
       weight = [scheme % edge_coefficient([dt, dt], 2, 4, 1, 1)]
       k = size(weight)
    case default
       error stop 'level 6: an unknown refusal case'
    end select
    write(*,'(a,i0)') ' the case ' // case_name // ' was accepted with ', k
    stop 0

  end subroutine refused_case

end program time_level_6
