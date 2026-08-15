!=====================================================================!
! ADJOINT TOWER . LEVEL 9 . THE COMPLETE STATEMENT
!
! Gate C. The last question: WHAT SENSITIVITY PROBLEM IS BEING
! ASKED. The statement is: at p = 2, solve R(q,p) = 0 and compute
! df/dp by the adjoint method. It SELECTS - it invents nothing, and
! it adds no coefficient of its own:
!
!      structure       R_dep and the four blocks, DERIVED as at
!                      Gate A and then OWNED by the model graph
!      constitution    the Level-8 law table, reused, never redone
!      parameter       p = 2, a field on P
!      requested       df/dp
!
! The primary road runs adjoint-only:
!
!      p -> primal solve -> q -> f_q^T -> adjoint solve -> lambda
!        -> f_p - lambda^T R_p -> 7
!
! Every derivative input on that road is GENERATED, never written:
! f_q^T is the response block read in reverse on a unit seed, R_p
! and f_p are the parameter blocks acting on a unit direction. The
! literals [1,2], [-4,-11] and 2 appear in assertions only. The
! tangent road runs afterwards, independently, and certifies the
! same 7 without ever having contributed to it.
!
! Two graphs stand here and they are NOT interchangeable. The MODEL
! is a related_graph owning the mathematical structure: the
! external selectors are located inside it by identity and then
! DESTROYED, and every number below comes through model-owned
! relations. The HOST is a seven-vertex stored_graph that exists
! only because the legacy graph_operation face demands a
! class(ordinary_graph); it is provably neither Q nor Y nor their size, and
! it contributes no domain, no coefficient and no topology.
!
! And the roles are enumerated HOSTILELY on purpose - Q = {v,u},
! Y = {r2,r1} - the configuration Level 8 proved a positional
! implementation cannot survive. The sealed tower inherits the hard
! case, not the easy one.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program adjoint_level_9

  use iso_fortran_env  , only : dp => REAL64
  use adjoint_assert   , only : report, verdict
  use adjoint_assert   , only : VAR_P, VAR_U, VAR_V
  use adjoint_assert   , only : TGT_R1, TGT_R2, TGT_F
  use graph_set    , only : index_set, subset, set
  use graph_relation   , only : stored_relation, relation
  use graph_relation_algebra, only : compose_binary
  use graph_binary_relation , only : csr_relation, transposed_view, &
       &                             transpose_of, inclusion_of
  use graph_structure  , only : related_graph, declared_set, declared_relation
  use graph_grammar    , only : graph_field
  use class_graph      , only : stored_graph
  use class_graph_field, only : field
  use class_graph_gmres, only : gmres
  use adjoint_constitution_fixture, only : constituted_primal, &
       & constituted_adjoint, constituted_tangent, &
       & response_of, rq_forward, rp_forward, fq_forward, fq_reverse, &
       & fp_forward

  implicit none

  type(index_set)                  :: v, t, hv
  type(subset)                   :: p_dom, q_dom, y_dom, z_dom
  type(stored_relation), allocatable :: dep
  type(csr_relation)   , allocatable, target :: inc_y, inc_z, inc_q, inc_p
  type(csr_relation)   , allocatable :: jq, jp, fq, fp
  type(related_graph), target     :: model
  type(stored_graph)                 :: host
  class(relation), pointer           :: rp   => null()
  class(relation), pointer           :: gdep => null()
  class(relation), pointer           :: gjq => null(), gjp => null()
  class(relation), pointer           :: gfq => null(), gfp => null()
  type(constituted_primal)           :: primal_eq
  type(constituted_adjoint)          :: adjoint_eq
  type(constituted_tangent)          :: tangent_eq
  type(gmres)                        :: primal_solver, adjoint_solver
  type(gmres)                        :: tangent_solver
  type(field)                        :: p_field, rhs_y, rhs_q
  class(set) , allocatable    :: dom, dom2
  class(graph_field), allocatable    :: sol
  real(dp), allocatable              :: pv(:), gv(:), qv(:), lv(:), qp(:)
  real(dp)                           :: fqt(2), rp_dir(2), resp(1)
  real(dp)                           :: direct(1), contrib(1), zbar(1)
  real(dp)                           :: lambda_rp, sensitivity
  real(dp)                           :: tangent_sensitivity
  integer                            :: table(2, 9)
  integer                            :: nfail, k

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "adjoint tower . level 9 . the statement"
  write(*,'(1x,a)') "============================================="

  !-- the structure the statement selects, HOSTILELY enumerated -----
  v = index_set('variables', 3)
  t = index_set('targets'  , 3)

  p_dom = subset('parameter', v, [VAR_P])
  q_dom = subset('state'    , v, [VAR_V, VAR_U])   ! reversed
  y_dom = subset('residual' , t, [TGT_R2, TGT_R1]) ! reversed
  z_dom = subset('response' , t, [TGT_F])

  table(:, 1) = [TGT_R1, VAR_P]
  table(:, 2) = [TGT_R1, VAR_U]
  table(:, 3) = [TGT_R1, VAR_V]
  table(:, 4) = [TGT_R2, VAR_P]
  table(:, 5) = [TGT_R2, VAR_U]
  table(:, 6) = [TGT_R2, VAR_V]
  table(:, 7) = [TGT_F , VAR_P]
  table(:, 8) = [TGT_F , VAR_U]
  table(:, 9) = [TGT_F , VAR_V]
  allocate(dep)
  dep = stored_relation('dependency', [t, v], table)

  !-- the four blocks, derived exactly as Gate A derives them -------
  allocate(inc_y, inc_z, inc_q, inc_p)
  inc_y = inclusion_of(y_dom)
  inc_z = inclusion_of(z_dom)
  inc_q = inclusion_of(q_dom)
  inc_p = inclusion_of(p_dom)

  allocate(jq, jp, fq, fp)
  jq = compose_binary(compose_binary(inc_y, dep), transpose_of(inc_q))
  jp = compose_binary(compose_binary(inc_y, dep), transpose_of(inc_p))
  fq = compose_binary(compose_binary(inc_z, dep), transpose_of(inc_q))
  fp = compose_binary(compose_binary(inc_z, dep), transpose_of(inc_p))

  !-- the MODEL graph owns the mathematics --------------------------
  model = related_graph('adjoint model', &
       & [declared_set(v), declared_set(t), declared_set(p_dom), declared_set(q_dom), &
       &  declared_set(y_dom), declared_set(z_dom)], &
       & [declared_relation(dep), declared_relation(jq), declared_relation(jp), &
       &  declared_relation(fq), declared_relation(fp)])

  !-- locate the model's own citizens by identity, then destroy the
  !   construction selectors: from here on the statement runs on
  !   graph-owned structure alone
  do k = 1, model % num_relations()
     rp => model % relation_at(k)
     if (rp % equals(dep)) gdep => rp
     if (rp % equals(jq))  gjq  => rp
     if (rp % equals(jp))  gjp  => rp
     if (rp % equals(fq))  gfq  => rp
     if (rp % equals(fp))  gfp  => rp
  end do

  deallocate(dep, jq, jp, fq, fp)
  deallocate(inc_y, inc_z, inc_q, inc_p)

  !-- the compatibility host: nobody's domain -----------------------
  host = stored_graph(7, tails=[1,2,3,4,5,6], heads=[2,3,4,5,6,7])

  call check_hostile_enumeration(nfail)
  call check_model_ownership(nfail)
  call check_host_is_not_the_model(nfail)

  call state_the_parameter(nfail)
  call solve_the_primal(nfail)
  call evaluate_the_response(nfail)
  call generate_the_response_covector(nfail)
  call solve_the_adjoint(nfail)
  call generate_the_parameter_actions(nfail)
  call assemble_the_sensitivity(nfail)
  call certify_with_the_tangent(nfail)

  ! The tower's one result: a real scalar, serialized as it stands.
  write(*,'(1x,a)', advance='no') "ADJOINT_RESULT ="
  write(*,'(es24.16)') sensitivity

  call verdict(nfail, "level 9")

contains

  !===================================================================!
  ! The sealed tower inherits Level 8's hard configuration, not the
  ! easy alignment: both two-member roles run backwards.
  !===================================================================!

  subroutine check_hostile_enumeration(nfail)

    integer, intent(inout) :: nfail

    call report(q_dom % local_index(VAR_V) .eq. 1 .and. &
         &      q_dom % local_index(VAR_U) .eq. 2, &
         & "Q is enumerated {v, u}: the state runs backwards", nfail)
    call report(y_dom % local_index(TGT_R2) .eq. 1 .and. &
         &      y_dom % local_index(TGT_R1) .eq. 2, &
         & "Y is enumerated {r2, r1}: and so do the residual rows", &
         & nfail)
    call report(.not. q_dom % equals(y_dom) .and. &
         &      q_dom % size() .eq. y_dom % size(), &
         & "and they are still two different domains of equal size", &
         & nfail)

  end subroutine check_hostile_enumeration

  !===================================================================!
  ! The lifetime proof: the construction selectors are gone, and the
  ! model's own relations answer with the signatures the statement
  ! needs. This is a LIFETIME truth, not another extension test.
  !===================================================================!

  subroutine check_model_ownership(nfail)

    integer, intent(inout) :: nfail

    call report(.not. allocated(dep) .and. .not. allocated(jq) .and. &
         &      .not. allocated(jp) .and. .not. allocated(fq) .and. &
         &      .not. allocated(fp), &
         & "every construction selector is destroyed", nfail)

    call report(associated(gdep) .and. associated(gjq) .and. &
         &      associated(gjp) .and. associated(gfq) .and. &
         &      associated(gfp), &
         & "and the model's own five relations were located by " // &
         & "identity before they died", nfail)

    call report(model % num_relations() .eq. 5 .and. &
         &      model % num_sets() .eq. 6, &
         & "the model still owns six sets and five relations", &
         & nfail)

    dom  = gjq % domain(1)
    dom2 = gjq % domain(2)
    call report(dom % equals(y_dom) .and. dom2 % equals(q_dom), &
         & "model-owned J_Q still runs Y x Q", nfail)

    dom  = gjp % domain(2)
    call report(dom % equals(p_dom), &
         & "model-owned J_P still answers for the parameter", nfail)

    dom  = gfq % domain(1)
    dom2 = gfq % domain(2)
    call report(dom % equals(z_dom) .and. dom2 % equals(q_dom), &
         & "model-owned F_Q still runs Z x Q", nfail)

    dom  = gfp % domain(2)
    call report(dom % equals(p_dom), &
         & "and model-owned F_P for the parameter", nfail)

  end subroutine check_model_ownership

  !===================================================================!
  ! Two graphs, two roles. The model owns the mathematics; the host
  ! satisfies a legacy signature and owns nothing that matters.
  !===================================================================!

  subroutine check_host_is_not_the_model(nfail)

    integer, intent(inout) :: nfail

    hv = host % vertex_set()

    call report(.not. hv % equals(q_dom) .and. &
         &      .not. hv % equals(y_dom), &
         & "the host's vertex set is neither Q nor Y", nfail)
    call report(host % num_vertices() .ne. q_dom % size() .and. &
         &      host % num_vertices() .ne. y_dom % size(), &
         & "nor of their size: the solver host is not the model", &
         & nfail)
    call report(.not. model % holds_set(hv), &
         & "and the model does not own it at all", nfail)

  end subroutine check_host_is_not_the_model

  !===================================================================!
  ! The statement's one given: p = 2, as a field on P, read through
  ! P's own map rather than passed as a naked literal.
  !===================================================================!

  subroutine state_the_parameter(nfail)

    integer, intent(inout) :: nfail

    p_field = field('parameter', p_dom)
    call p_field % set_real_vector([2.0_dp])
    call p_field % get_real_vector(pv)

    call p_field % domain(dom)
    call report(dom % equals(p_dom) .and. &
         &      abs(pv(p_dom % local_index(VAR_P)) - 2.0_dp) < 1.0d-14, &
         & "the statement states p = 2 on P", nfail)

  end subroutine state_the_parameter

  !===================================================================!
  ! R(q, p) = 0, through model-owned supports and the Level-8 law.
  !===================================================================!

  subroutine solve_the_primal(nfail)

    integer, intent(inout) :: nfail

    primal_eq = constituted_primal(gjq, gjp, q_dom, y_dom, p_dom, pv)

    call primal_solver % attach(primal_eq, host, q_dom)
    primal_solver % tolerance      = 1.0d-12
    primal_solver % max_iterations = 50

    call primal_solver % domain(host, dom)
    call report(dom % equals(q_dom), &
         & "the primal solver answers on Q", nfail)

    call primal_solver % constant(gv)
    rhs_y = field('rhs', y_dom)
    call rhs_y % set_real_vector(-gv)
    call primal_solver % apply(host, [rhs_y], sol)

    call sol % domain(dom)
    call sol % get_real_vector(qv)
    call report(dom % equals(q_dom), &
         & "and the state comes back as a field on Q", nfail)
    call report(abs(qv(q_dom % local_index(VAR_U)) - 2.0_dp) < 1.0d-9, &
         & "u = 2, by member", nfail)
    call report(abs(qv(q_dom % local_index(VAR_V)) - 4.0_dp) < 1.0d-9, &
         & "v = 4, by member", nfail)

  end subroutine solve_the_primal

  !===================================================================!
  ! f(q, p) at the solved state, from the model's response blocks.
  !===================================================================!

  subroutine evaluate_the_response(nfail)

    integer, intent(inout) :: nfail

    call response_of(gfq, gfp, z_dom, q_dom, p_dom, qv, pv, resp)

    call report(abs(resp(z_dom % local_index(TGT_F)) - 14.0_dp) &
         &      < 1.0d-9, &
         & "f = 14 at the solved state", nfail)

  end subroutine evaluate_the_response

  !===================================================================!
  ! f_q^T, GENERATED: the response block read in reverse on a unit
  ! seed. The vector [1, 2] is an oracle here, never an input - and
  ! the right-hand side the adjoint solver actually consumes is
  ! proved to be this generated field.
  !===================================================================!

  subroutine generate_the_response_covector(nfail)

    integer, intent(inout) :: nfail

    zbar = 1.0_dp
    call fq_reverse(gfq, z_dom, q_dom, zbar, fqt)

    call report(abs(fqt(q_dom % local_index(VAR_U)) - 1.0_dp) < 1.0d-12 &
         & .and. abs(fqt(q_dom % local_index(VAR_V)) - 2.0_dp) &
         &      < 1.0d-12, &
         & "f_q^T = [1, 2] on Q - generated by reverse action, not " // &
         & "written down", nfail)

  end subroutine generate_the_response_covector

  !===================================================================!
  ! R_q^T lambda = f_q^T, through the SAME solver family with the
  ! orientation exchanged. The equation generates its own constant
  ! from the same reverse reading, so the right-hand side the solver
  ! consumes IS the covector generated above - proved, not assumed.
  !===================================================================!

  subroutine solve_the_adjoint(nfail)

    integer, intent(inout) :: nfail

    adjoint_eq = constituted_adjoint(gjq, gfq, y_dom, q_dom, z_dom)

    call adjoint_solver % attach(adjoint_eq, host, y_dom)
    adjoint_solver % tolerance      = 1.0d-12
    adjoint_solver % max_iterations = 50

    call adjoint_solver % domain(host, dom)
    call report(dom % equals(y_dom), &
         & "the adjoint solver answers on Y - the orientation is " // &
         & "exchanged, the machinery is not", nfail)

    call adjoint_solver % constant(gv)
    call report(maxval(abs(-gv - fqt)) < 1.0d-12, &
         & "the right-hand side it consumes IS the generated " // &
         & "f_q^T", nfail)

    rhs_q = field('rhs', q_dom)
    call rhs_q % set_real_vector(-gv)
    call adjoint_solver % apply(host, [rhs_q], sol)

    call sol % domain(dom)
    call sol % get_real_vector(lv)
    call report(dom % equals(y_dom), &
         & "lambda comes back as a field on Y", nfail)
    call report(abs(lv(y_dom % local_index(TGT_R1)) + 0.4_dp) < 1.0d-9, &
         & "lambda(r1) = -0.4, by member", nfail)
    call report(abs(lv(y_dom % local_index(TGT_R2)) - 0.6_dp) < 1.0d-9, &
         & "lambda(r2) = 0.6, by member", nfail)
    call report(abs(lv(y_dom % local_index(TGT_R1)) - 0.4_dp) > 0.1_dp, &
         & "and not the primal-orientation [0.4, 0.2]", nfail)

  end subroutine solve_the_adjoint

  !===================================================================!
  ! R_p and f_p, GENERATED by acting on the unit parameter
  ! direction. The literals [-4, -11] and 2 are oracles only.
  !===================================================================!

  subroutine generate_the_parameter_actions(nfail)

    integer, intent(inout) :: nfail

    call rp_forward(gjp, y_dom, p_dom, [1.0_dp], rp_dir)
    call report(abs(rp_dir(y_dom % local_index(TGT_R1)) + 4.0_dp) &
         &      < 1.0d-12 .and. &
         &      abs(rp_dir(y_dom % local_index(TGT_R2)) + 11.0_dp) &
         &      < 1.0d-12, &
         & "R_p = [-4, -11] on Y, generated by the parameter action", &
         & nfail)

    call fp_forward(gfp, z_dom, p_dom, [1.0_dp], direct)
    call report(abs(direct(z_dom % local_index(TGT_F)) - 2.0_dp) &
         &      < 1.0d-12, &
         & "f_p = 2 on Z, generated the same way", nfail)

  end subroutine generate_the_parameter_actions

  !===================================================================!
  ! THE STATEMENT'S ANSWER: df/dp = f_p - lambda^T R_p, adjoint-only.
  ! The pairing is taken over Y, and both operands are proved to
  ! live there before a single product is formed.
  !===================================================================!

  subroutine assemble_the_sensitivity(nfail)

    integer, intent(inout) :: nfail

    call sol % domain(dom)
    call report(dom % equals(y_dom) .and. size(lv) .eq. size(rp_dir), &
         & "lambda and R_p share the residual domain: the pairing " // &
         & "is legitimate", nfail)

    lambda_rp = dot_product(lv, rp_dir)
    call report(abs(lambda_rp + 5.0_dp) < 1.0d-9, &
         & "lambda^T R_p = -5", nfail)

    sensitivity = direct(z_dom % local_index(TGT_F)) - lambda_rp
    call report(abs(sensitivity - 7.0_dp) < 1.0d-9, &
         & "df/dp = f_p - lambda^T R_p = 7 - the adjoint road, " // &
         & "alone", nfail)

  end subroutine assemble_the_sensitivity

  !===================================================================!
  ! The tangent road, run AFTERWARDS and never consulted above:
  ! R_q q_p = -R_p, then f_q q_p + f_p. It certifies the answer; it
  ! does not produce it.
  !===================================================================!

  subroutine certify_with_the_tangent(nfail)

    integer, intent(inout) :: nfail

    tangent_eq = constituted_tangent(gjq, gjp, q_dom, y_dom, p_dom, &
         &                           [1.0_dp])

    call tangent_solver % attach(tangent_eq, host, q_dom)
    tangent_solver % tolerance      = 1.0d-12
    tangent_solver % max_iterations = 50

    call tangent_solver % constant(gv)
    rhs_y = field('rhs', y_dom)
    call rhs_y % set_real_vector(-gv)
    call tangent_solver % apply(host, [rhs_y], sol)
    call sol % get_real_vector(qp)

    call report(abs(qp(q_dom % local_index(VAR_U)) - 1.0_dp) < 1.0d-9 &
         & .and. abs(qp(q_dom % local_index(VAR_V)) - 2.0_dp) &
         &      < 1.0d-9, &
         & "q_p = [1, 2] on Q", nfail)

    call fq_forward(gfq, z_dom, q_dom, qp, contrib)
    tangent_sensitivity = contrib(z_dom % local_index(TGT_F)) + &
         &                direct(z_dom % local_index(TGT_F))

    call report(abs(tangent_sensitivity - 7.0_dp) < 1.0d-9, &
         & "the tangent road independently reaches 7", nfail)
    call report(abs(tangent_sensitivity - sensitivity) < 1.0d-9, &
         & "and it agrees with the adjoint answer it never touched", &
         & nfail)

  end subroutine certify_with_the_tangent

end program adjoint_level_9
