!=====================================================================!
! LEARNING TOWER . LEVEL 8 . CONSTITUTION
!
! The level answers one question: WHAT THE SYMBOLS MEAN - and
! whether meaning plus the existing structure GENERATES the very
! residual Level 7 was handed as an oracle. Only here bind:
!
!      predict(w, x) = w * x        error(yhat, y) = yhat - y
!
! and the semantics are the learning tower's own: Theta is
! trainable state, K = {y, x} is observed state, U = {e, yhat} is
! COMPUTED - so evaluation runs the laws INTO the computed slots
! along the DERIVED order [predict, error], and the residual is the
! VALUE at the home L locates: r = e = w*x - y. Never q(e) - law:
! there is no independent q(e). The structure supplies execution
! order (restrict, project, compose, admit, interpret, sort); the
! constitution supplies laws; L supplies the home; the evaluator
! hard-codes no slot, no order, no location. The generated map
! reproduces Level 7's oracle - w = 0 gives -6, w = 3 gives 0 - and
! a second data instance (8, 4) proves data is not constitution.
! Meaning does not bend topology: the trainable slots the
! constituted evaluation actually READS equal Level 6's structural
! support {w}, re-derived here without one numerical perturbation.
! No solver runs; no parameter changes; no derivative exists.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program learning_level_8

  use iso_fortran_env, only : dp => REAL64
  use learning_assert, only : report, verdict
  use learning_assert, only : SLOT_W, SLOT_X, SLOT_YHAT, SLOT_Y, SLOT_E
  use learning_assert, only : OP_PREDICT, OP_ERROR
  use learning_assert, only : PORT_IN1, PORT_IN2, PORT_OUT
  use graph_fractal        , only : graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set        , only : set_map
  use map_inclusion  , only : inclusion_map, declared_subobject
  use relation_finitary , only : stored_relation, relation
  use relation_algebra, only : restrict_slot, project_slots, &
       &                             compose_binary
  use relation_algorithms, only : reachable, topological_order
  use field_stored, only : stored_field
  use learning_constitution_fixture, only : apply_law, slot_for_port, &
       &                                    located_slot, generated_residual
  use graph_fractal        , only : graph, known_branch, null_branch
  use view_relational, only : relational_binding, &
       & num_member_sets, member_set_at, num_relations, relation_at, &
       & has_set

  implicit none


  ! Residual rows exist only from Level 6 upward.
  integer, parameter :: ROW_R = 1

  type(graph)              :: v, o, p, y
  type(graph)               :: k, theta, u, p_in, p_out
  type(stored_relation)          :: flow, backwards, located
  type(stored_relation)          :: consumes, produces
  class(relation), allocatable   :: d, a, d2
  type(graph)             , target :: g
  type(graph)             , target :: scell(3), selem(3)
  type(graph)             , target :: rcell(2), relem(2)
  type(relational_binding)         :: bnd
  integer                          :: kcell
  type(graph)             , target :: g_a
  type(graph)             , target :: scell2(1), selem2(1)
  type(graph)             , target :: rcell2(1), relem2(1)
  type(relational_binding)         :: bnd2
  integer                          :: kcell2
  type(graph)             , target :: g2
  type(graph)             , target :: scell3(3), selem3(3)
  type(graph)             , target :: rcell3(2), relem3(2)
  type(relational_binding)         :: bnd3
  integer                          :: kcell3
  type(stored_field)                    :: q_k, q_k2
  integer, allocatable           :: order(:), order2(:)
  real(dp), allocatable          :: obs(:)
  integer                        :: table(3, 6)
  integer                        :: nfail
  type(set_map)     :: sets
  type(inclusion_map)     :: inclusions

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "learning tower . level 8 . constitution"
  write(*,'(1x,a)') "============================================="

  call v % declare()
  call sets % bind(v, counted_set_representation(5))
  call o % declare()
  call sets % bind(o, counted_set_representation(2))
  call p % declare()
  call sets % bind(p, counted_set_representation(3))
  call y % declare()
  call sets % bind(y, counted_set_representation(1))
  call k % declare()
  call sets       % bind(k, listed_set_representation([SLOT_Y, SLOT_X]))
  call inclusions % include_in(k, v)
  call theta % declare()
  call sets       % bind(theta, listed_set_representation([SLOT_W]))
  call inclusions % include_in(theta, v)
  call u % declare()
  call sets       % bind(u, listed_set_representation([SLOT_E, SLOT_YHAT]))
  call inclusions % include_in(u, v)

  table(:, 1) = [OP_PREDICT, SLOT_W   , PORT_IN1]
  table(:, 2) = [OP_PREDICT, SLOT_X   , PORT_IN2]
  table(:, 3) = [OP_PREDICT, SLOT_YHAT, PORT_OUT]
  table(:, 4) = [OP_ERROR  , SLOT_YHAT, PORT_IN1]
  table(:, 5) = [OP_ERROR  , SLOT_Y   , PORT_IN2]
  table(:, 6) = [OP_ERROR  , SLOT_E   , PORT_OUT]
  flow = stored_relation('flow', [o, v, p], table, sets)

  located = stored_relation('located', [y, v], &
       & reshape([ROW_R, SLOT_E], [2, 1]), sets)

  call p_in % declare()
  call sets       % bind(p_in, listed_set_representation([PORT_IN1, PORT_IN2]))
  call inclusions % include_in(p_in, p)
  call p_out % declare()
  call sets       % bind(p_out, listed_set_representation([PORT_OUT]))
  call inclusions % include_in(p_out, p)

  ! The Level-5 observed field, carried forward: [6, 2] on K = {y, x}
  ! in declaration order. U gets NO field - before evaluation the
  ! computed domain owns no numbers, and that stays true.
  q_k = stored_field('observations', k, sets % num_members_of(k))
  call q_k % set_real_vector([6.0_dp, 2.0_dp])
  call q_k % real_vector(obs)

  ! The lower road, walked again: D derived, admitted, interpreted,
  ! sorted - the structure supplies the execution order.
  consumes = project_slots(restrict_slot(flow, 3, p_in , sets, inclusions), [2, 1], sets)
  produces = project_slots(restrict_slot(flow, 3, p_out, sets, inclusions), [1, 2], sets)
  d = compose_binary(produces, consumes, sets)
  ! 'learning': (S, P) as one sequence on each branch.
  call g % declare()
  do kcell = 1, 3
     call scell(kcell) % declare()
     call selem(kcell) % declare()
  end do
  do kcell = 1, 2
     call rcell(kcell) % declare()
     call relem(kcell) % declare()
  end do

  call bnd % bind_set(selem(1), v)
  call bnd % bind_set(selem(2), o)
  call bnd % bind_set(selem(3), p)
  call bnd % bind_relation(relem(1), flow)
  call bnd % bind_relation(relem(2), d)

  do kcell = 1, 3
     scell(kcell) % branch(1) = known_branch(selem(kcell))
     if (kcell .lt. 3) scell(kcell) % branch(2) = &
          & known_branch(scell(kcell + 1))
  end do
  do kcell = 1, 2
     rcell(kcell) % branch(1) = known_branch(relem(kcell))
     if (kcell .lt. 2) rcell(kcell) % branch(2) = &
          & known_branch(rcell(kcell + 1))
  end do

  g % branch(1) = known_branch(scell(1))
  g % branch(2) = known_branch(rcell(1))
  call topological_order(d, sets, order)

  call check_derived_order(nfail)
  call check_laws(nfail)
  call check_observed_readback(nfail)
  call check_generated_map(nfail)
  call check_intermediates(nfail)
  call check_topology_preserved(nfail)
  call check_data_witness(nfail)
  call check_order_invariance(nfail)

  call verdict(nfail, "level 8")

contains

  !===================================================================!
  ! The execution order is DERIVED, never handed over: the evaluator
  ! receives what the topological sort answered.
  !===================================================================!

  subroutine check_derived_order(nfail)

    integer, intent(inout) :: nfail

    call report(size(order) .eq. 2 .and. &
         &      order(1) .eq. OP_PREDICT .and. order(2) .eq. OP_ERROR, &
         & "the derived order is [predict, error], exactly", nfail)

  end subroutine check_derived_order

  !===================================================================!
  ! The law table, certified independently: this is the exact place
  ! where numerical meaning enters. The unbound-symbol refusal dies
  ! in the refusal suite.
  !===================================================================!

  subroutine check_laws(nfail)

    integer, intent(inout) :: nfail

    call report(abs(apply_law(OP_PREDICT, 3.0_dp, 2.0_dp) - 6.0_dp) &
         &      < 1.0d-14, &
         & "predict is multiplication: predict(3, 2) = 6", nfail)
    call report(abs(apply_law(OP_ERROR, 6.0_dp, 6.0_dp)) < 1.0d-14, &
         & "error is subtraction: error(6, 6) = 0", nfail)
    call report(abs(apply_law(OP_ERROR, 4.0_dp, 6.0_dp) + 2.0_dp) &
         &      < 1.0d-14, &
         & "and error(4, 6) = -2: signed, not a distance", nfail)

  end subroutine check_laws

  !===================================================================!
  ! The observed values, read back through K's own enumeration -
  ! Level 5's truth, one line each, not re-proved wholesale.
  !===================================================================!

  subroutine check_observed_readback(nfail)

    integer, intent(inout) :: nfail

    call report(abs(obs(sets % index_in(k, SLOT_Y)) - 6.0_dp) < 1.0d-14 &
         & .and. abs(obs(sets % index_in(k, SLOT_X)) - 2.0_dp) < 1.0d-14, &
         & "K still answers y = 6 and x = 2, by enumeration", nfail)

  end subroutine check_observed_readback

  !===================================================================!
  ! The generated map, probed OFF-SOLUTION: structure + laws + data
  ! + L reproduce r(w) = 2w - 6 - the very map Level 7 received
  ! opaquely. The formula lives in this comment and these expected
  ! numbers; it does not live in the evaluator.
  !===================================================================!

  subroutine check_generated_map(nfail)

    integer, intent(inout) :: nfail

    real(dp) :: r(1)

    call generated_residual(flow, located, v, sets, y, k, obs, &
         & theta, [0.0_dp], u, order, r)
    call report(abs(r(1) + 6.0_dp) < 1.0d-12, &
         & "w = 0 generates r = -6: Level 7's affine constant", nfail)

    call generated_residual(flow, located, v, sets, y, k, obs, &
         & theta, [1.0_dp], u, order, r)
    call report(abs(r(1) + 4.0_dp) < 1.0d-12, &
         & "w = 1 generates r = -4", nfail)

    call generated_residual(flow, located, v, sets, y, k, obs, &
         & theta, [-1.0_dp], u, order, r)
    call report(abs(r(1) + 8.0_dp) < 1.0d-12, &
         & "w = -1 generates r = -8", nfail)

    call generated_residual(flow, located, v, sets, y, k, obs, &
         & theta, [3.0_dp], u, order, r)
    call report(abs(r(1)) < 1.0d-12, &
         & "w = 3 generates r = 0: the solution, evaluated - " // &
         & "never solved for", nfail)

  end subroutine check_generated_map

  !===================================================================!
  ! The computed slots acquire their values DURING evaluation: at
  ! w = 3, the laws produce yhat = 6 and e = 0 on U - slots that
  ! owned no number before the run.
  !===================================================================!

  subroutine check_intermediates(nfail)

    integer, intent(inout) :: nfail

    real(dp)              :: r(1)
    real(dp), allocatable :: trace(:)

    call generated_residual(flow, located, v, sets, y, k, obs, &
         & theta, [3.0_dp], u, order, r, trace=trace)

    call report(abs(trace(sets % index_in(v, SLOT_YHAT)) - 6.0_dp) &
         &      < 1.0d-12, &
         & "the laws computed yhat = 6 on the computed domain", nfail)
    call report(abs(trace(sets % index_in(v, SLOT_E))) < 1.0d-12, &
         & "and e = 0 - produced, never preseeded", nfail)

  end subroutine check_intermediates

  !===================================================================!
  ! The Level-8 preservation law: meaning added /= topology changed.
  ! The structural support is re-derived exactly as Level 6 derived
  ! it - A from the flow, reachability to the home read from L -
  ! and held against the trainable slots the constituted evaluation
  ! actually READ. No numerical perturbation anywhere.
  !===================================================================!

  subroutine check_topology_preserved(nfail)

    integer, intent(inout) :: nfail

    real(dp)             :: r(1)
    integer, allocatable :: touched(:)
    integer, allocatable              :: sup(:)
    integer              :: ti, nsup, home
    logical              :: same

    allocate(sup(sets % num_members_of(theta)))
    a = compose_binary(consumes, produces, sets)
    ! 'value dependency': (S, P) as one sequence on each branch.
    call g_a % declare()
    do kcell2 = 1, 1
       call scell2(kcell2) % declare()
       call selem2(kcell2) % declare()
    end do
    do kcell2 = 1, 1
       call rcell2(kcell2) % declare()
       call relem2(kcell2) % declare()
    end do

    call bnd2 % bind_set(selem2(1), v)
    call bnd2 % bind_relation(relem2(1), a)

    do kcell2 = 1, 1
       scell2(kcell2) % branch(1) = known_branch(selem2(kcell2))
       if (kcell2 .lt. 1) scell2(kcell2) % branch(2) = &
            & known_branch(scell2(kcell2 + 1))
    end do
    do kcell2 = 1, 1
       rcell2(kcell2) % branch(1) = known_branch(relem2(kcell2))
       if (kcell2 .lt. 1) rcell2(kcell2) % branch(2) = &
            & known_branch(rcell2(kcell2 + 1))
    end do

    g_a % branch(1) = known_branch(scell2(1))
    g_a % branch(2) = known_branch(rcell2(1))

    home = located_slot(located, v, sets, ROW_R)
    nsup = 0
    do ti = 1, sets % num_members_of(theta)
       if (reachable(a, sets, sets % member_of(theta, ti), home)) then
          nsup      = nsup + 1
          sup(nsup) = sets % member_of(theta, ti)
       end if
    end do

    call report(nsup .eq. 1 .and. sup(1) .eq. SLOT_W, &
         & "the structural support of r is { w }, re-derived", nfail)

    call generated_residual(flow, located, v, sets, y, k, obs, &
         & theta, [1.0_dp], u, order, r, touched=touched)

    call report(size(touched) .eq. 1 .and. touched(1) .eq. SLOT_W, &
         & "the constituted evaluation read exactly { w }", nfail)

    same = size(touched) .eq. nsup
    do ti = 1, nsup
       same = same .and. any(touched .eq. sup(ti))
    end do
    do ti = 1, size(touched)
       same = same .and. any(sup(1:nsup) .eq. touched(ti))
    end do
    call report(same, &
         & "meaning added, topology unchanged: the two supports " // &
         & "agree extensionally", nfail)

  end subroutine check_topology_preserved

  !===================================================================!
  ! Data is not constitution: the SAME flow, laws and evaluator,
  ! handed observations (8, 4), generate the other model's map.
  !===================================================================!

  subroutine check_data_witness(nfail)

    integer, intent(inout) :: nfail

    real(dp)              :: r(1)
    real(dp), allocatable :: obs2(:)

    q_k2 = stored_field('observations again', k, sets % num_members_of(k))
    call q_k2 % set_real_vector([8.0_dp, 4.0_dp])
    call q_k2 % real_vector(obs2)

    call generated_residual(flow, located, v, sets, y, k, obs2, &
         & theta, [2.0_dp], u, order, r)
    call report(abs(r(1)) < 1.0d-12, &
         & "with (y, x) = (8, 4): w = 2 generates r = 0", nfail)

    call generated_residual(flow, located, v, sets, y, k, obs2, &
         & theta, [1.0_dp], u, order, r)
    call report(abs(r(1) + 4.0_dp) < 1.0d-12, &
         & "and w = 1 generates r = -4: data changed, " // &
         & "constitution did not", nfail)

  end subroutine check_data_witness

  !===================================================================!
  ! The six facts handed backwards: their OWN derivation - D2,
  ! order2 - drives the same constitution to the same residual and
  ! the same trainable support at a nontrivial probe.
  !===================================================================!

  subroutine check_order_invariance(nfail)

    integer, intent(inout) :: nfail

    integer              :: rev(3, 6), j
    integer, allocatable :: touched(:), touched2(:)
    real(dp)             :: r_fwd(1), r_back(1)

    do j = 1, 6
       rev(:, j) = table(:, 7 - j)
    end do
    backwards = stored_relation('flow backwards', [o, v, p], rev, sets)

    d2 = compose_binary( &
         & project_slots(restrict_slot(backwards, 3, p_out, sets, inclusions), [1, 2], sets), &
         & project_slots(restrict_slot(backwards, 3, p_in , sets, inclusions), [2, 1], sets), sets)
    ! 'learning backwards': (S, P) as one sequence on each branch.
    call g2 % declare()
    do kcell3 = 1, 3
       call scell3(kcell3) % declare()
       call selem3(kcell3) % declare()
    end do
    do kcell3 = 1, 2
       call rcell3(kcell3) % declare()
       call relem3(kcell3) % declare()
    end do

    call bnd3 % bind_set(selem3(1), v)
    call bnd3 % bind_set(selem3(2), o)
    call bnd3 % bind_set(selem3(3), p)
    call bnd3 % bind_relation(relem3(1), backwards)
    call bnd3 % bind_relation(relem3(2), d2)

    do kcell3 = 1, 3
       scell3(kcell3) % branch(1) = known_branch(selem3(kcell3))
       if (kcell3 .lt. 3) scell3(kcell3) % branch(2) = &
            & known_branch(scell3(kcell3 + 1))
    end do
    do kcell3 = 1, 2
       rcell3(kcell3) % branch(1) = known_branch(relem3(kcell3))
       if (kcell3 .lt. 2) rcell3(kcell3) % branch(2) = &
            & known_branch(rcell3(kcell3 + 1))
    end do

    g2 % branch(1) = known_branch(scell3(1))
    g2 % branch(2) = known_branch(rcell3(1))
    call topological_order(d2, sets, order2)

    call generated_residual(flow, located, v, sets, y, k, obs, &
         & theta, [-2.0_dp], u, order, r_fwd, touched=touched)
    call generated_residual(backwards, located, v, sets, y, k, obs, &
         & theta, [-2.0_dp], u, order2, r_back, touched=touched2)

    call report(abs(r_fwd(1) - r_back(1)) < 1.0d-14 .and. &
         &      abs(r_fwd(1) + 10.0_dp) < 1.0d-12, &
         & "w = -2: both enumerations generate r = -10, " // &
         & "identically", nfail)
    call report(size(touched2) .eq. size(touched) .and. &
         &      touched2(1) .eq. touched(1), &
         & "and the same trainable support, either way written", &
         & nfail)

  end subroutine check_order_invariance

end program learning_level_8
