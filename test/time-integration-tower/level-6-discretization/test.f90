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
! The action S : Q -> Q carries its own domain and stores no graph.
! The compatibility host H_t has FIVE vertices in a chain, and Q
! has TWO members. The mismatch is deliberate and load-bearing: if
! the two carriers had the same size, a substitution of one for the
! other would produce plausible numbers and the seam would hide.
!
!      |V(H_t)| = 5        |Q| = 2        and they are not same_as
!
! H_t is a COMPATIBILITY HOST - the conduit the graph_operation
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
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program time_level_6

  use iso_fortran_env       , only : dp => REAL64
  use time_assert           , only : report, verdict
  use time_assert           , only : NQ, NT, TOL
  use time_assert           , only : T0, T1, T2
  use time_assert           , only : H_STEP, Q0, Q_FE1, Q_BE1, Q_BDF2
  use time_assert           , only : action_of
  use fractal_graph        , only : set_graph => graph
  use graph_set_map        , only : set_map
  use graph_ordinary_view   , only : ordinary_graph
  use graph_field_calculus  , only : graph_field
  use graph_binary_relation , only : csr_relation
  use class_graph           , only : ordinary_stored_graph
  use class_graph_field     , only : field
  use class_graph_step      , only : step_operator, backward_euler, bdf
  use time_carriers_fixture , only : time_carriers
  use time_relations_fixture, only : tail_relation, head_relation
  use time_algebra_fixture  , only : derive_one_step_reach, &
       &                             derive_two_step_reach
  use time_fields_fixture   , only : state_field
  use triangular_decay_fixture, only : triangular_decay

  implicit none

  type(set_graph)          :: q, t, e
  type(set_map)          :: sets
  type(csr_relation), target :: tail, head, a1
  type(csr_relation)         :: a2
  type(ordinary_stored_graph)         :: ht
  type(triangular_decay)     :: decay
  type(field)                :: qf
  integer                    :: nfail

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
  ht = ordinary_stored_graph(NT, tails=[1,2,3,4], heads=[2,3,4,5])

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

  call verdict(nfail, "level 6")

contains

  !===================================================================!
  ! The mismatch that makes the experiment discriminating. If these
  ! ever became equal, every assertion below would still pass while
  ! proving nothing.
  !===================================================================!

  subroutine check_host_is_not_the_state_domain(nfail)

    integer, intent(inout) :: nfail

    type(set_graph) :: hv

    hv = ht % vertex_set()

    call report(ht % num_vertices() .eq. NT .and. sets % size_of(q) .eq. NQ, &
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

    type(set_graph) :: d, hv
    integer        :: n_d

    hv = ht % vertex_set()
    d   = ht % all_vertices()
    n_d = ht % num_vertices()
    call report(d % same_as(hv), &
         & "all_vertices(H_t) and H_t % vertex_set() are the same " // &
         & "carrier - so delegating a domain question changes no " // &
         & "graph-based caller's answer", nfail)

  end subroutine check_host_carriers_agree_with_themselves

  !===================================================================!
  ! THE first half of the experiment, and it needs no production
  ! change: an action that carries its own domain already works
  ! through the graph_operation face.
  !===================================================================!

  subroutine check_direct_action_preserves_q(nfail)

    integer, intent(inout) :: nfail

    class(graph_field), allocatable :: answer
    type(set_graph)  :: d
    integer         :: n_d
    real(dp), allocatable           :: s(:)

    call decay % domain(ht, d, n_d)
    call report(d % same_as(q), &
         & "the ACTION answers Q when asked its domain, though it " // &
         & "was handed a five-vertex host", nfail)

    call decay % apply(ht, [qf], answer)
    d = answer % domain()
    call report(d % same_as(q), &
         & "and it ANSWERS on Q: graph host and state domain are " // &
         & "independent concepts in this specimen", nfail)

    call answer % get_real_vector(s)
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
  ! input_graph % all_vertices(...) and this assertion FAILED,
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

    type(step_operator)            :: step
    type(set_graph) :: d, hv
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

    type(step_operator)             :: step
    type(field)                     :: state
    class(graph_field), allocatable :: r
    type(set_graph)  :: d
    real(dp), allocatable           :: v(:)

    step = backward_euler(decay, H_STEP)
    step % qold = Q0

    ! At the exact backward-euler state the residual vanishes.
    state = field('trial', q, NQ, ncomp=1)
    call state % set_real_vector(Q_BE1)
    call step % apply(ht, [state], r)

    d = r % domain()
    call report(d % same_as(q), &
         & "the backward-euler RESIDUAL lands on Q, not on the " // &
         & "host's vertices", nfail)

    call r % get_real_vector(v)
    call report(size(v) .eq. NQ .and. maxval(abs(v)) .lt. TOL, &
         & "and it is ZERO at q = [4/3, 4/9]: the exact discrete " // &
         & "backward-euler state, verified by substitution", nfail)

    ! And it is not zero anywhere convenient - the forward-euler
    ! answer is a different number, as it must be.
    call state % set_real_vector(Q_FE1)
    call step % apply(ht, [state], r)
    call r % get_real_vector(v)
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

    type(step_operator) :: scheme

    call report(a1 % has([T1, T2]) .and. a2 % has([T0, T2]), &
         & "at instant t2 the one-step predecessor is t1 and the " // &
         & "two-step predecessor is t0 - STRUCTURAL REACH, from " // &
         & "Level 2", nfail)

    scheme = bdf(2, decay, H_STEP)
    call report(scheme % reach .eq. 2, &
         & "and bdf-2 reaches exactly two instants back, matching " // &
         & "the reach A2 describes", nfail)

    call report(abs(scheme % a0 - 1.5_dp) .lt. TOL .and. &
         &      abs(scheme % a1 + 2.0_dp) .lt. TOL .and. &
         &      abs(scheme % a2 - 0.5_dp) .lt. TOL, &
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

    type(step_operator)             :: step
    type(field)                     :: state
    class(graph_field), allocatable :: r
    type(set_graph)  :: d
    real(dp), allocatable           :: v(:)

    step = bdf(2, decay, H_STEP)
    step % qold   = Q_BE1
    step % qolder = Q0

    state = field('trial', q, NQ, ncomp=1)
    call state % set_real_vector(Q_BDF2)
    call step % apply(ht, [state], r)

    d = r % domain()
    call report(d % same_as(q), &
         & "the bdf-2 RESIDUAL lands on Q as well - a two-step " // &
         & "scheme changes the coefficients, never the domain", nfail)

    call r % get_real_vector(v)
    call report(size(v) .eq. NQ .and. maxval(abs(v)) .lt. TOL, &
         & "and it is ZERO at q2 = [5/6, 47/72], started from the " // &
         & "backward-euler q1", nfail)

  end subroutine check_bdf2_residual

end program time_level_6
