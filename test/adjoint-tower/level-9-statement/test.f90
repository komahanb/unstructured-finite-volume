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
! is a graph read as (S, P), with a binding to the objects its
! elements denote, and it carries the mathematical structure: the
! external selectors are located inside it by identity and then
! DESTROYED, and every number below comes through bound relations. The HOST is a seven-vertex stored_graph that exists
! only because the legacy graph_operation face demands a
! class(graph); it is provably neither Q nor Y nor their size, and
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
  use fractal_graph        , only : set_graph => graph
  use graph_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use graph_set_map        , only : set_map
  use graph_label_map      , only : label_map
  use graph_inclusion_map  , only : inclusion_map, declared_subobject
  use graph_relation   , only : stored_relation, relation
  use graph_relation_algebra, only : compose_binary
  use graph_binary_relation , only : csr_relation, transposed_view, &
       &                             transpose_of, inclusion_of
  use fractal_graph        , only : graph, known_branch, null_branch
  use graph_relational_view, only : relational_binding, &
       & num_member_sets, member_set_at, num_relations, relation_at, &
       & holds_set
  use graph_field_calculus, only : graph_field
  use class_graph      , only : stored_graph
  use class_graph_field, only : field
  use class_graph_gmres, only : gmres
  use adjoint_constitution_fixture, only : constituted_primal, &
       & constituted_adjoint, constituted_tangent, &
       & response_of, rq_forward, rp_forward, fq_forward, fq_reverse, &
       & fp_forward

  implicit none

  type(set_graph)                  :: v, t, hv
  type(set_graph)                   :: p_dom, q_dom, y_dom, z_dom
  type(stored_relation), allocatable :: dep
  type(csr_relation)   , allocatable, target :: inc_y, inc_z, inc_q, inc_p
  type(csr_relation)   , allocatable :: jq, jp, fq, fp
  type(graph)             , target :: model
  type(graph)             , target :: scell(6), selem(6)
  type(graph)             , target :: rcell(5), relem(5)
  type(relational_binding)         :: bnd
  integer                          :: kcell
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
  type(set_graph) , allocatable    :: dom, dom2
  class(graph_field), allocatable    :: sol
  real(dp), allocatable              :: pv(:), gv(:), qv(:), lv(:), qp(:)
  real(dp)                           :: fqt(2), rp_dir(2), resp(1)
  real(dp)                           :: direct(1), contrib(1), zbar(1)
  real(dp)                           :: lambda_rp, sensitivity
  real(dp)                           :: tangent_sensitivity
  integer                            :: table(2, 9)
  integer                            :: nfail, k
  type(set_map)     :: sets
  type(label_map)     :: labels
  type(inclusion_map)     :: inclusions

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "adjoint tower . level 9 . the statement"
  write(*,'(1x,a)') "============================================="

  !-- the structure the statement selects, HOSTILELY enumerated -----
  call v % declare()
  call sets % bind(v, counted_set_representation(3))
  call t % declare()
  call sets % bind(t, counted_set_representation(3))

  call p_dom % declare()
  call sets       % bind(p_dom, listed_set_representation([VAR_P]))
  call inclusions % include_in(p_dom, v)
  call q_dom % declare()
  call sets       % bind(q_dom, listed_set_representation([VAR_V, VAR_U]))
  call inclusions % include_in(q_dom, v)
  call y_dom % declare()
  call sets       % bind(y_dom, listed_set_representation([TGT_R2, TGT_R1]))
  call inclusions % include_in(y_dom, t)
  call z_dom % declare()
  call sets       % bind(z_dom, listed_set_representation([TGT_F]))
  call inclusions % include_in(z_dom, t)

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
  dep = stored_relation('dependency', [t, v], table, sets)

  !-- the four blocks, derived exactly as Gate A derives them -------
  allocate(inc_y, inc_z, inc_q, inc_p)
  inc_y = inclusion_of(y_dom, t, sets, labels)
  inc_z = inclusion_of(z_dom, t, sets, labels)
  inc_q = inclusion_of(q_dom, v, sets, labels)
  inc_p = inclusion_of(p_dom, v, sets, labels)

  allocate(jq, jp, fq, fp)
  jq = compose_binary(compose_binary(inc_y, dep, sets), transpose_of(inc_q), sets)
  jp = compose_binary(compose_binary(inc_y, dep, sets), transpose_of(inc_p), sets)
  fq = compose_binary(compose_binary(inc_z, dep, sets), transpose_of(inc_q), sets)
  fp = compose_binary(compose_binary(inc_z, dep, sets), transpose_of(inc_p), sets)

  !-- the MODEL graph owns the mathematics --------------------------
  ! 'adjoint model': (S, P) as one sequence on each branch.
  call model % declare()
  do kcell = 1, 6
     call scell(kcell) % declare()
     call selem(kcell) % declare()
  end do
  do kcell = 1, 5
     call rcell(kcell) % declare()
     call relem(kcell) % declare()
  end do

  call bnd % bind_set(selem(1), v)
  call bnd % bind_set(selem(2), t)
  call bnd % bind_set(selem(3), p_dom)
  call bnd % bind_set(selem(4), q_dom)
  call bnd % bind_set(selem(5), y_dom)
  call bnd % bind_set(selem(6), z_dom)
  call bnd % bind_relation(relem(1), dep)
  call bnd % bind_relation(relem(2), jq)
  call bnd % bind_relation(relem(3), jp)
  call bnd % bind_relation(relem(4), fq)
  call bnd % bind_relation(relem(5), fp)

  do kcell = 1, 6
     scell(kcell) % branch(1) = known_branch(selem(kcell))
     if (kcell .lt. 6) scell(kcell) % branch(2) = &
          & known_branch(scell(kcell + 1))
  end do
  do kcell = 1, 5
     rcell(kcell) % branch(1) = known_branch(relem(kcell))
     if (kcell .lt. 5) rcell(kcell) % branch(2) = &
          & known_branch(rcell(kcell + 1))
  end do

  model % branch(1) = known_branch(scell(1))
  model % branch(2) = known_branch(rcell(1))

  !-- locate the model's own citizens by identity, then destroy the
  !   construction selectors: from here on the statement runs on
  !   graph-owned structure alone
  do k = 1, num_relations(model)
     rp => relation_at(model, bnd, k)
     if (rp % same_as(dep)) gdep => rp
     if (rp % same_as(jq))  gjq  => rp
     if (rp % same_as(jp))  gjp  => rp
     if (rp % same_as(fq))  gfq  => rp
     if (rp % same_as(fp))  gfp  => rp
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

    call report(sets % index_in(q_dom, VAR_V) .eq. 1 .and. &
         &      sets % index_in(q_dom, VAR_U) .eq. 2, &
         & "Q is enumerated {v, u}: the state runs backwards", nfail)
    call report(sets % index_in(y_dom, TGT_R2) .eq. 1 .and. &
         &      sets % index_in(y_dom, TGT_R1) .eq. 2, &
         & "Y is enumerated {r2, r1}: and so do the residual rows", &
         & nfail)
    call report(.not. q_dom % same_as(y_dom) .and. &
         &      sets % size_of(q_dom) .eq. sets % size_of(y_dom), &
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

    call report(num_relations(model) .eq. 5 .and. &
         &      num_member_sets(model) .eq. 6, &
         & "the model still owns six carriers and five relations", &
         & nfail)

    dom  = gjq % domain(1)
    dom2 = gjq % domain(2)
    call report(dom % same_as(y_dom) .and. dom2 % same_as(q_dom), &
         & "model-owned J_Q still runs Y x Q", nfail)

    dom  = gjp % domain(2)
    call report(dom % same_as(p_dom), &
         & "model-owned J_P still answers for the parameter", nfail)

    dom  = gfq % domain(1)
    dom2 = gfq % domain(2)
    call report(dom % same_as(z_dom) .and. dom2 % same_as(q_dom), &
         & "model-owned F_Q still runs Z x Q", nfail)

    dom  = gfp % domain(2)
    call report(dom % same_as(p_dom), &
         & "and model-owned F_P for the parameter", nfail)

  end subroutine check_model_ownership

  !===================================================================!
  ! Two graphs, two roles. The model owns the mathematics; the host
  ! satisfies a legacy signature and owns nothing that matters.
  !===================================================================!

  subroutine check_host_is_not_the_model(nfail)

    integer, intent(inout) :: nfail

    hv = host % vertex_set()

    call report(.not. hv % same_as(q_dom) .and. &
         &      .not. hv % same_as(y_dom), &
         & "the host's vertex set is neither Q nor Y", nfail)
    call report(host % num_vertices() .ne. sets % size_of(q_dom) .and. &
         &      host % num_vertices() .ne. sets % size_of(y_dom), &
         & "nor of their size: the solver host is not the model", &
         & nfail)
    call report(.not. holds_set(model, bnd, hv), &
         & "and the model does not own it at all", nfail)

  end subroutine check_host_is_not_the_model

  !===================================================================!
  ! The statement's one given: p = 2, as a field on P, read through
  ! P's own map rather than passed as a naked literal.
  !===================================================================!

  subroutine state_the_parameter(nfail)

    integer, intent(inout) :: nfail

    p_field = field('parameter', p_dom, sets % size_of(p_dom))
    call p_field % set_real_vector([2.0_dp])
    call p_field % get_real_vector(pv)

    dom = p_field % domain()
    call report(dom % same_as(p_dom) .and. &
         &      abs(pv(sets % index_in(p_dom, VAR_P)) - 2.0_dp) < 1.0d-14, &
         & "the statement states p = 2 on P", nfail)

  end subroutine state_the_parameter

  !===================================================================!
  ! R(q, p) = 0, through model-owned supports and the Level-8 law.
  !===================================================================!

  subroutine solve_the_primal(nfail)

    integer, intent(inout) :: nfail
    integer         :: n_dom

    primal_eq = constituted_primal(gjq, gjp, q_dom, y_dom, p_dom, pv, sets)

    call primal_solver % attach(primal_eq, host, q_dom, sets % size_of(q_dom))
    primal_solver % tolerance      = 1.0d-12
    primal_solver % max_iterations = 50

    call primal_solver % domain(host, dom, n_dom)
    call report(dom % same_as(q_dom), &
         & "the primal solver answers on Q", nfail)

    call primal_solver % constant(gv)
    rhs_y = field('rhs', y_dom, sets % size_of(y_dom))
    call rhs_y % set_real_vector(-gv)
    call primal_solver % apply(host, [rhs_y], sol)

    dom = sol % domain()
    call sol % get_real_vector(qv)
    call report(dom % same_as(q_dom), &
         & "and the state comes back as a field on Q", nfail)
    call report(abs(qv(sets % index_in(q_dom, VAR_U)) - 2.0_dp) < 1.0d-9, &
         & "u = 2, by member", nfail)
    call report(abs(qv(sets % index_in(q_dom, VAR_V)) - 4.0_dp) < 1.0d-9, &
         & "v = 4, by member", nfail)

  end subroutine solve_the_primal

  !===================================================================!
  ! f(q, p) at the solved state, from the model's response blocks.
  !===================================================================!

  subroutine evaluate_the_response(nfail)

    integer, intent(inout) :: nfail

    call response_of(gfq, gfp, z_dom, q_dom, p_dom, sets, qv, pv, resp)

    call report(abs(resp(sets % index_in(z_dom, TGT_F)) - 14.0_dp) &
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
    call fq_reverse(gfq, z_dom, q_dom, sets, zbar, fqt)

    call report(abs(fqt(sets % index_in(q_dom, VAR_U)) - 1.0_dp) < 1.0d-12 &
         & .and. abs(fqt(sets % index_in(q_dom, VAR_V)) - 2.0_dp) &
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
    integer         :: n_dom

    adjoint_eq = constituted_adjoint(gjq, gfq, y_dom, q_dom, z_dom, sets)

    call adjoint_solver % attach(adjoint_eq, host, y_dom, sets % size_of(y_dom))
    adjoint_solver % tolerance      = 1.0d-12
    adjoint_solver % max_iterations = 50

    call adjoint_solver % domain(host, dom, n_dom)
    call report(dom % same_as(y_dom), &
         & "the adjoint solver answers on Y - the orientation is " // &
         & "exchanged, the machinery is not", nfail)

    call adjoint_solver % constant(gv)
    call report(maxval(abs(-gv - fqt)) < 1.0d-12, &
         & "the right-hand side it consumes IS the generated " // &
         & "f_q^T", nfail)

    rhs_q = field('rhs', q_dom, sets % size_of(q_dom))
    call rhs_q % set_real_vector(-gv)
    call adjoint_solver % apply(host, [rhs_q], sol)

    dom = sol % domain()
    call sol % get_real_vector(lv)
    call report(dom % same_as(y_dom), &
         & "lambda comes back as a field on Y", nfail)
    call report(abs(lv(sets % index_in(y_dom, TGT_R1)) + 0.4_dp) < 1.0d-9, &
         & "lambda(r1) = -0.4, by member", nfail)
    call report(abs(lv(sets % index_in(y_dom, TGT_R2)) - 0.6_dp) < 1.0d-9, &
         & "lambda(r2) = 0.6, by member", nfail)
    call report(abs(lv(sets % index_in(y_dom, TGT_R1)) - 0.4_dp) > 0.1_dp, &
         & "and not the primal-orientation [0.4, 0.2]", nfail)

  end subroutine solve_the_adjoint

  !===================================================================!
  ! R_p and f_p, GENERATED by acting on the unit parameter
  ! direction. The literals [-4, -11] and 2 are oracles only.
  !===================================================================!

  subroutine generate_the_parameter_actions(nfail)

    integer, intent(inout) :: nfail

    call rp_forward(gjp, y_dom, p_dom, sets, [1.0_dp], rp_dir)
    call report(abs(rp_dir(sets % index_in(y_dom, TGT_R1)) + 4.0_dp) &
         &      < 1.0d-12 .and. &
         &      abs(rp_dir(sets % index_in(y_dom, TGT_R2)) + 11.0_dp) &
         &      < 1.0d-12, &
         & "R_p = [-4, -11] on Y, generated by the parameter action", &
         & nfail)

    call fp_forward(gfp, z_dom, p_dom, sets, [1.0_dp], direct)
    call report(abs(direct(sets % index_in(z_dom, TGT_F)) - 2.0_dp) &
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

    dom = sol % domain()
    call report(dom % same_as(y_dom) .and. size(lv) .eq. size(rp_dir), &
         & "lambda and R_p share the residual domain: the pairing " // &
         & "is legitimate", nfail)

    lambda_rp = dot_product(lv, rp_dir)
    call report(abs(lambda_rp + 5.0_dp) < 1.0d-9, &
         & "lambda^T R_p = -5", nfail)

    sensitivity = direct(sets % index_in(z_dom, TGT_F)) - lambda_rp
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
         &                           [1.0_dp], sets)

    call tangent_solver % attach(tangent_eq, host, q_dom, sets % size_of(q_dom))
    tangent_solver % tolerance      = 1.0d-12
    tangent_solver % max_iterations = 50

    call tangent_solver % constant(gv)
    rhs_y = field('rhs', y_dom, sets % size_of(y_dom))
    call rhs_y % set_real_vector(-gv)
    call tangent_solver % apply(host, [rhs_y], sol)
    call sol % get_real_vector(qp)

    call report(abs(qp(sets % index_in(q_dom, VAR_U)) - 1.0_dp) < 1.0d-9 &
         & .and. abs(qp(sets % index_in(q_dom, VAR_V)) - 2.0_dp) &
         &      < 1.0d-9, &
         & "q_p = [1, 2] on Q", nfail)

    call fq_forward(gfq, z_dom, q_dom, sets, qp, contrib)
    tangent_sensitivity = contrib(sets % index_in(z_dom, TGT_F)) + &
         &                direct(sets % index_in(z_dom, TGT_F))

    call report(abs(tangent_sensitivity - 7.0_dp) < 1.0d-9, &
         & "the tangent road independently reaches 7", nfail)
    call report(abs(tangent_sensitivity - sensitivity) < 1.0d-9, &
         & "and it agrees with the adjoint answer it never touched", &
         & nfail)

  end subroutine certify_with_the_tangent

end program adjoint_level_9
