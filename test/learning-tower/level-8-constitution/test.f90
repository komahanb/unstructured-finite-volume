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
  use graph_carrier  , only : counted_set, subset_set, member_set
  use graph_relation , only : stored_relation, relation
  use graph_relation_algebra, only : restrict_slot, project_slots, &
       &                             compose_binary
  use graph_structure, only : relational_graph, held_set, held_relation
  use graph_profile  , only : directed_adjacency_view
  use graph_algorithms, only : reachable, topological_order
  use class_graph_field, only : field
  use learning_constitution_fixture, only : apply_law, slot_for_port, &
       &                                    located_slot, generated_residual
  use fractal_graph        , only : graph
  use graph_relational_view, only : relational_binding
  use relational_fixture   , only : fractal_fixture

  implicit none

  type(fractal_fixture)             :: fx_
  type(graph)             , pointer :: fg_
  type(relational_binding), pointer :: fb_

  ! Residual rows exist only from Level 6 upward.
  integer, parameter :: ROW_R = 1

  type(counted_set)              :: v, o, p, y
  type(subset_set)               :: k, theta, u, p_in, p_out
  type(stored_relation)          :: flow, backwards, located
  type(stored_relation)          :: consumes, produces
  class(relation), allocatable   :: d, a, d2
  type(relational_graph), target :: g, g_a, g2
  type(directed_adjacency_view)  :: view, dep_view, view2
  type(field)                    :: q_k, q_k2
  integer, allocatable           :: order(:), order2(:)
  real(dp), allocatable          :: obs(:)
  integer                        :: table(3, 6)
  integer                        :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "learning tower . level 8 . constitution"
  write(*,'(1x,a)') "============================================="

  v     = counted_set('value-slots'  , 5)
  o     = counted_set('operations'   , 2)
  p     = counted_set('ports'        , 3)
  y     = counted_set('residual-rows', 1)
  k     = subset_set('observed' , v, [SLOT_Y, SLOT_X])
  theta = subset_set('trainable', v, [SLOT_W])
  u     = subset_set('computed' , v, [SLOT_E, SLOT_YHAT])

  table(:, 1) = [OP_PREDICT, SLOT_W   , PORT_IN1]
  table(:, 2) = [OP_PREDICT, SLOT_X   , PORT_IN2]
  table(:, 3) = [OP_PREDICT, SLOT_YHAT, PORT_OUT]
  table(:, 4) = [OP_ERROR  , SLOT_YHAT, PORT_IN1]
  table(:, 5) = [OP_ERROR  , SLOT_Y   , PORT_IN2]
  table(:, 6) = [OP_ERROR  , SLOT_E   , PORT_OUT]
  flow = stored_relation('flow', [o, v, p], table)

  located = stored_relation('located', [y, v], &
       & reshape([ROW_R, SLOT_E], [2, 1]))

  p_in  = subset_set('input-ports', p, [PORT_IN1, PORT_IN2])
  p_out = subset_set('output-port', p, [PORT_OUT])

  ! The Level-5 observed field, carried forward: [6, 2] on K = {y, x}
  ! in declaration order. U gets NO field - before evaluation the
  ! computed domain owns no numbers, and that stays true.
  q_k = field('observations', k)
  call q_k % set_real_vector([6.0_dp, 2.0_dp])
  call q_k % get_real_vector(obs)

  ! The lower road, walked again: D derived, admitted, interpreted,
  ! sorted - the structure supplies the execution order.
  consumes = project_slots(restrict_slot(flow, 3, p_in ), [2, 1])
  produces = project_slots(restrict_slot(flow, 3, p_out), [1, 2])
  d = compose_binary(produces, consumes)
  g = relational_graph('learning', &
       & [held_set(v), held_set(o), held_set(p)], &
       & [held_relation(flow), held_relation(d)])
  call fx_ % to_fractal(g, fg_, fb_)
  view = directed_adjacency_view(fg_, fb_, d)
  call topological_order(view, order)

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

    call report(abs(obs(k % local_index(SLOT_Y)) - 6.0_dp) < 1.0d-14 &
         & .and. abs(obs(k % local_index(SLOT_X)) - 2.0_dp) < 1.0d-14, &
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

    call generated_residual(flow, located, v, y, k, obs, &
         & theta, [0.0_dp], u, order, r)
    call report(abs(r(1) + 6.0_dp) < 1.0d-12, &
         & "w = 0 generates r = -6: Level 7's affine constant", nfail)

    call generated_residual(flow, located, v, y, k, obs, &
         & theta, [1.0_dp], u, order, r)
    call report(abs(r(1) + 4.0_dp) < 1.0d-12, &
         & "w = 1 generates r = -4", nfail)

    call generated_residual(flow, located, v, y, k, obs, &
         & theta, [-1.0_dp], u, order, r)
    call report(abs(r(1) + 8.0_dp) < 1.0d-12, &
         & "w = -1 generates r = -8", nfail)

    call generated_residual(flow, located, v, y, k, obs, &
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

    call generated_residual(flow, located, v, y, k, obs, &
         & theta, [3.0_dp], u, order, r, trace=trace)

    call report(abs(trace(v % local_index(SLOT_YHAT)) - 6.0_dp) &
         &      < 1.0d-12, &
         & "the laws computed yhat = 6 on the computed domain", nfail)
    call report(abs(trace(v % local_index(SLOT_E))) < 1.0d-12, &
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
    integer              :: sup(theta % size())
    integer              :: ti, nsup, home
    logical              :: same

    a = compose_binary(consumes, produces)
    g_a = relational_graph('value dependency', [held_set(v)], &
         & [held_relation(a)])
    call fx_ % to_fractal(g_a, fg_, fb_)
    dep_view = directed_adjacency_view(fg_, fb_, a)

    home = located_slot(located, v, ROW_R)
    nsup = 0
    do ti = 1, theta % size()
       if (reachable(dep_view, theta % member(ti), home)) then
          nsup      = nsup + 1
          sup(nsup) = theta % member(ti)
       end if
    end do

    call report(nsup .eq. 1 .and. sup(1) .eq. SLOT_W, &
         & "the structural support of r is { w }, re-derived", nfail)

    call generated_residual(flow, located, v, y, k, obs, &
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

    q_k2 = field('observations again', k)
    call q_k2 % set_real_vector([8.0_dp, 4.0_dp])
    call q_k2 % get_real_vector(obs2)

    call generated_residual(flow, located, v, y, k, obs2, &
         & theta, [2.0_dp], u, order, r)
    call report(abs(r(1)) < 1.0d-12, &
         & "with (y, x) = (8, 4): w = 2 generates r = 0", nfail)

    call generated_residual(flow, located, v, y, k, obs2, &
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
    backwards = stored_relation('flow backwards', [o, v, p], rev)

    d2 = compose_binary( &
         & project_slots(restrict_slot(backwards, 3, p_out), [1, 2]), &
         & project_slots(restrict_slot(backwards, 3, p_in ), [2, 1]))
    g2 = relational_graph('learning backwards', &
         & [held_set(v), held_set(o), held_set(p)], &
         & [held_relation(backwards), held_relation(d2)])
    call fx_ % to_fractal(g2, fg_, fb_)
    view2 = directed_adjacency_view(fg_, fb_, d2)
    call topological_order(view2, order2)

    call generated_residual(flow, located, v, y, k, obs, &
         & theta, [-2.0_dp], u, order, r_fwd, touched=touched)
    call generated_residual(backwards, located, v, y, k, obs, &
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
