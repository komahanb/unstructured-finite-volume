!=====================================================================!
! ADJOINT TOWER . LEVEL 8 . ONE CONSTITUTION, BOTH DIRECTIONS
!
! The level answers one question: CAN ONE LAW SOURCE GENERATE
! EVERYTHING. Level 7's two independently supplied equations are
! retired here. A single coefficient table - A, b, c, d, each
! entry written once and keyed by the members it relates - now
! produces:
!
!      R(q,p)      the primal residual        Q -> Y
!      f(q,p)      the response               Q,P -> Z
!      Rq v        the forward state action   Q -> Y
!      Rq^T lam    the reverse state action   Y -> Q
!      Rp dp       the parameter action       P -> Y
!      fq v, fq^T  the response's two readings
!      fp dp       the direct response action P -> Z
!
! and no A^T is written anywhere: the reverse action walks the same
! J_Q with the same coeff_state, accumulating into the state slots.
! Every action is driven by the Gate-A structural supports rather
! than by a 2x2 loop, and every vector is indexed through its
! domain's own local_index.
!
! Three roads then meet on one number, and none may call another:
!
!      primal solve      ->  q = [2,4],  f = 14
!      adjoint solve     ->  lambda = [-0.4, 0.6]
!      tangent solve     ->  q_p = [1,2]
!
!      tangent  df/dp = fq q_p + fp        = 7
!      adjoint  df/dp = fp - lambda^T Rp   = 7
!
! The whole battery runs TWICE: once with Q = [u,v], Y = [r1,r2],
! and once with the roles INDEPENDENTLY permuted to Q' = [v,u],
! Y' = [r2,r1]. Every answer is checked BY MEMBER, so the stored
! vectors come back in different orders and the truths do not move.
! An implementation that quietly assumed the i-th state matches the
! i-th residual row would pass the first run and fail the second.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program adjoint_level_8

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
  use graph_grammar    , only : graph_field
  use class_graph      , only : stored_graph
  use class_graph_field, only : field
  use class_graph_gmres, only : gmres
  use adjoint_constitution_fixture, only : constituted_primal, &
       & constituted_adjoint, constituted_tangent, &
       & response_of, rq_forward, rq_reverse, rp_forward, &
       & fq_forward, fp_forward

  implicit none

  type(set_graph)     :: v, t
  type(set_graph)      :: p_dom, z_dom
  type(set_graph)      :: q_can, y_can, q_perm, y_perm
  type(stored_graph)    :: host
  type(stored_relation) :: dep
  integer               :: table(2, 9)
  integer               :: nfail
  type(set_map)     :: sets
  type(label_map)     :: labels
  type(inclusion_map)     :: inclusions

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "adjoint tower . level 8 . one constitution"
  write(*,'(1x,a)') "============================================="

  call v % declare()
  call sets % bind(v, counted_set_representation(3))
  call t % declare()
  call sets % bind(t, counted_set_representation(3))

  ! The one-member roles cannot be permuted; the two-member ones can.
  call p_dom % declare()
  call sets       % bind(p_dom, listed_set_representation([VAR_P]))
  call inclusions % include_in(p_dom, v)
  call z_dom % declare()
  call sets       % bind(z_dom, listed_set_representation([TGT_F]))
  call inclusions % include_in(z_dom, t)

  call q_can % declare()
  call sets       % bind(q_can, listed_set_representation([VAR_U, VAR_V]))
  call inclusions % include_in(q_can, v)
  call y_can % declare()
  call sets       % bind(y_can, listed_set_representation([TGT_R1, TGT_R2]))
  call inclusions % include_in(y_can, t)
  call q_perm % declare()
  call sets       % bind(q_perm, listed_set_representation([VAR_V, VAR_U]))
  call inclusions % include_in(q_perm, v)
  call y_perm % declare()
  call sets       % bind(y_perm, listed_set_representation([TGT_R2, TGT_R1]))
  call inclusions % include_in(y_perm, t)

  table(:, 1) = [TGT_R1, VAR_P]
  table(:, 2) = [TGT_R1, VAR_U]
  table(:, 3) = [TGT_R1, VAR_V]
  table(:, 4) = [TGT_R2, VAR_P]
  table(:, 5) = [TGT_R2, VAR_U]
  table(:, 6) = [TGT_R2, VAR_V]
  table(:, 7) = [TGT_F , VAR_P]
  table(:, 8) = [TGT_F , VAR_U]
  table(:, 9) = [TGT_F , VAR_V]
  dep = stored_relation('dependency', [t, v], table, sets)

  host = stored_graph(5, tails=[1,2,3,4], heads=[2,3,4,5])

  ! The same battery, twice: once aligned, once with BOTH two-member
  ! roles independently reversed. The expected storage positions are
  ! passed in so the run must prove it really is enumerated the way
  ! its name claims - otherwise "permuted" could quietly be a second
  ! canonical run.
  call run_battery(q_can , y_can , sets, "canonical", 1, 1, nfail)
  call run_battery(q_perm, y_perm, sets, "permuted" , 2, 2, nfail)

  call verdict(nfail, "level 8")

contains

  !===================================================================!
  ! The whole of Level 8, parameterized by how the two-member roles
  ! happen to be enumerated. Nothing below reads a raw position.
  !===================================================================!

  subroutine run_battery(q_dom, y_dom, sets, tag, u_at, r1_at, nfail)

    type(set_graph) , intent(in)    :: q_dom, y_dom
    type(set_map)   , intent(in)    :: sets
    character(len=*) , intent(in)    :: tag
    integer          , intent(in)    :: u_at, r1_at
    integer          , intent(inout) :: nfail

    type(csr_relation), target      :: inc_y, inc_z, inc_q, inc_p
    type(transposed_view)           :: inc_q_t, inc_p_t
    type(csr_relation)              :: jq, jp, fq, fp
    type(constituted_primal)        :: primal_eq
    type(constituted_adjoint)       :: adjoint_eq
    type(constituted_tangent)       :: tangent_eq
    type(gmres)                     :: primal_solver, adjoint_solver
    type(gmres)                     :: tangent_solver
    type(field)                     :: rhs_y, rhs_q, state
    type(set_graph)  :: dom
    class(graph_field), allocatable :: sol
    real(dp), allocatable           :: gv(:), qv(:), lv(:), qp(:)
    real(dp), allocatable           :: rp_dir(:)
    real(dp), allocatable                        :: resp(:)
    real(dp), allocatable                        :: contrib(:)
    real(dp), allocatable                        :: direct(:)
    real(dp)                        :: tangent_sens, adjoint_sens
    real(dp)                        :: lam_rp

    allocate(rp_dir(sets % size_of(y_dom)))
    allocate(direct(sets % size_of(z_dom)))
    allocate(contrib(sets % size_of(z_dom)))
    allocate(resp(sets % size_of(z_dom)))
    write(*,'(1x,a,a,a)') "--- ", tag, " enumeration ---"

    ! Prove the enumeration is what the tag claims, so a "permuted"
    ! run cannot silently be a second canonical one.
    call report(sets % index_in(q_dom, VAR_U) .eq. u_at .and. &
         &      sets % index_in(y_dom, TGT_R1) .eq. r1_at, &
         & tag // ": u is stored at position " // &
         & achar(48 + u_at) // " and r1 at position " // &
         & achar(48 + r1_at), nfail)

    ! The Gate-A supports, re-derived for THESE role identities.
    inc_y   = inclusion_of(y_dom, t, sets, labels)
    inc_z   = inclusion_of(z_dom, t, sets, labels)
    inc_q   = inclusion_of(q_dom, v, sets, labels)
    inc_p   = inclusion_of(p_dom, v, sets, labels)
    inc_q_t = transpose_of(inc_q)
    inc_p_t = transpose_of(inc_p)

    jq = compose_binary(compose_binary(inc_y, dep, sets), inc_q_t, sets)
    jp = compose_binary(compose_binary(inc_y, dep, sets), inc_p_t, sets)
    fq = compose_binary(compose_binary(inc_z, dep, sets), inc_q_t, sets)
    fp = compose_binary(compose_binary(inc_z, dep, sets), inc_p_t, sets)

    primal_eq  = constituted_primal(jq, jp, q_dom, y_dom, p_dom, &
         &                          [2.0_dp], sets)
    adjoint_eq = constituted_adjoint(jq, fq, y_dom, q_dom, z_dom, sets)
    tangent_eq = constituted_tangent(jq, jp, q_dom, y_dom, p_dom, &
         &                           [1.0_dp], sets)

    !-- the primal, through the constituted residual ----------------
    call primal_solver % attach(primal_eq, host, q_dom, sets % size_of(q_dom))
    primal_solver % tolerance      = 1.0d-12
    primal_solver % max_iterations = 50

    call primal_solver % constant(gv)
    rhs_y = field('rhs', y_dom, sets % size_of(y_dom))
    call rhs_y % set_real_vector(-gv)
    call primal_solver % apply(host, [rhs_y], sol)

    dom = sol % domain()
    call sol % get_real_vector(qv)
    call report(dom % same_as(q_dom) .and. &
         &      abs(qv(sets % index_in(q_dom, VAR_U)) - 2.0_dp) < 1.0d-9 &
         & .and. abs(qv(sets % index_in(q_dom, VAR_V)) - 4.0_dp) &
         &      < 1.0d-9, &
         & tag // ": the constituted primal gives u = 2, v = 4 on Q", &
         & nfail)

    !-- the response at the solved state ----------------------------
    call response_of(fq, fp, z_dom, q_dom, p_dom, sets, qv, [2.0_dp], resp)
    call report(abs(resp(sets % index_in(z_dom, TGT_F)) - 14.0_dp) &
         &      < 1.0d-9, &
         & tag // ": f = 14 at the solved state", nfail)

    !-- the adjoint, through the reverse reading of the same law ----
    call adjoint_solver % attach(adjoint_eq, host, y_dom, sets % size_of(y_dom))
    adjoint_solver % tolerance      = 1.0d-12
    adjoint_solver % max_iterations = 50

    call adjoint_solver % constant(gv)
    rhs_q = field('rhs', q_dom, sets % size_of(q_dom))
    call rhs_q % set_real_vector(-gv)
    call adjoint_solver % apply(host, [rhs_q], sol)

    dom = sol % domain()
    call sol % get_real_vector(lv)
    call report(dom % same_as(y_dom) .and. &
         &      abs(lv(sets % index_in(y_dom, TGT_R1)) + 0.4_dp) < 1.0d-9 &
         & .and. abs(lv(sets % index_in(y_dom, TGT_R2)) - 0.6_dp) &
         &      < 1.0d-9, &
         & tag // ": the constituted adjoint gives lambda = " // &
         & "[-0.4, 0.6] on Y", nfail)
    call report(abs(lv(sets % index_in(y_dom, TGT_R1)) - 0.4_dp) &
         &      > 0.1_dp, &
         & tag // ": and not the primal-orientation [0.4, 0.2]", &
         & nfail)

    !-- the tangent, an INDEPENDENT road to the same sensitivity ----
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
         & tag // ": q_p = [1, 2] on Q", nfail)

    call fq_forward(fq, z_dom, q_dom, sets, qp, contrib)
    call fp_forward(fp, z_dom, p_dom, sets, [1.0_dp], direct)
    tangent_sens = contrib(sets % index_in(z_dom, TGT_F)) + &
         &         direct(sets % index_in(z_dom, TGT_F))
    call report(abs(tangent_sens - 7.0_dp) < 1.0d-9, &
         & tag // ": tangent df/dp = fq q_p + fp = 7", nfail)

    !-- the adjoint road, meeting it without ever calling it -------
    call rp_forward(jp, y_dom, p_dom, sets, [1.0_dp], rp_dir)
    lam_rp = dot_product(lv, rp_dir)
    call report(abs(lam_rp + 5.0_dp) < 1.0d-9, &
         & tag // ": lambda^T Rp = -5", nfail)
    call report(abs(direct(sets % index_in(z_dom, TGT_F)) - 2.0_dp) &
         &      < 1.0d-12, &
         & tag // ": fp = 2, the response's DIRECT parameter term", &
         & nfail)

    adjoint_sens = direct(sets % index_in(z_dom, TGT_F)) - lam_rp
    call report(abs(adjoint_sens - 7.0_dp) < 1.0d-9, &
         & tag // ": adjoint df/dp = fp - lambda^T Rp = 7", nfail)
    call report(abs(tangent_sens - adjoint_sens) < 1.0d-9, &
         & tag // ": the two roads meet on one number, neither " // &
         & "having called the other", nfail)

    !-- duality, through the constituted actions only ---------------
    call check_duality(jq, q_dom, y_dom, sets, tag, nfail)

  end subroutine run_battery

  !===================================================================!
  ! < mu, A v >_Y = < A^T mu, v >_Q, both sides computed from the
  ! one law table - the forward reading and the reverse reading of
  ! the same coefficients. No dense transpose stands in the path.
  !===================================================================!

  subroutine check_duality(jq, q_dom, y_dom, sets, tag, nfail)

    class(relation)  , intent(in)    :: jq
    type(set_graph) , intent(in)    :: q_dom, y_dom
    type(set_map)   , intent(in)    :: sets
    character(len=*) , intent(in)    :: tag
    integer          , intent(inout) :: nfail

    real(dp), allocatable :: vq(:), mu(:), av(:), atmu(:)
    real(dp) :: lhs, rhs

    allocate(vq(sets % size_of(q_dom)), mu(sets % size_of(y_dom)))
    allocate(av(sets % size_of(y_dom)), atmu(sets % size_of(q_dom)))
    ! Seeded BY MEMBER: v = [-1 at u, 4 at v], mu = [2 at r1, -3 at r2]
    vq(sets % index_in(q_dom, VAR_U)) = -1.0_dp
    vq(sets % index_in(q_dom, VAR_V)) =  4.0_dp
    mu(sets % index_in(y_dom, TGT_R1)) =  2.0_dp
    mu(sets % index_in(y_dom, TGT_R2)) = -3.0_dp

    call rq_forward(jq, y_dom, q_dom, sets, vq, av)
    call rq_reverse(jq, y_dom, q_dom, sets, mu, atmu)

    call report(abs(av(sets % index_in(y_dom, TGT_R1)) - 2.0_dp) < 1.0d-12 &
         & .and. abs(av(sets % index_in(y_dom, TGT_R2)) - 13.0_dp) &
         &      < 1.0d-12, &
         & tag // ": A v = [2, 13] on Y", nfail)
    call report(abs(atmu(sets % index_in(q_dom, VAR_U)) + 5.0_dp) &
         &      < 1.0d-12 .and. &
         &      abs(atmu(sets % index_in(q_dom, VAR_V)) + 10.0_dp) &
         &      < 1.0d-12, &
         & tag // ": A^T mu = [-5, -10] on Q", nfail)

    lhs = dot_product(mu, av)
    rhs = dot_product(atmu, vq)
    call report(abs(lhs - rhs) < 1.0d-12 .and. &
         &      abs(lhs + 35.0_dp) < 1.0d-12, &
         & tag // ": < mu, Av >_Y = < A^T mu, v >_Q = -35", nfail)

  end subroutine check_duality

end program adjoint_level_8
