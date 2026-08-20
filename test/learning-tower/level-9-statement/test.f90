!=====================================================================!
! LEARNING TOWER . LEVEL 9 . THE STATEMENT
!
! The level asks the only question left: WHAT LEARNING PROBLEM IS
! BEING ASKED. The statement is: given the observation (x, y) =
! (2, 6) and the model structure predict -> yhat, error -> e with
! the Level-8 constitution, learn w from w0 = 0 such that r = 0.
! It SELECTS - it invents nothing:
!
!      structure       T_flow; D DERIVED; order DERIVED
!      observation     x = 2, y = 6 on K = { y, x }
!      trainable       Theta = { w }
!      computed        U = { e, yhat } - valueless until laws run
!      constitution    the Level-8 law table, reused, never redone
!      location        L = { (r, e) }
!      requested       w, on the trainable domain
!
! The residual operation is a test-local adapter holding the
! GRAPH-OWNED flow (the external selector dies before the solve),
! delegating every number to the Level-8 evaluator; ordinary GMRES
! solves through its own operation face - rhs on Y in, solution on
! Theta out - discovering the linear action from the constituted
! model alone. Level 7's affine oracle appears NOWHERE in this
! path. The learned parameter is read only through Theta's own
! enumeration; the literal 3 stands only in final assertions; and
! training changes a VALUE, never the topology. Afterward, the
! SAME predict law - its operands chosen by structure - serves
! inference: x* = 4 gives yhat* = 12, no retraining, no target.
!
! The tower's PRIMARY result is the learned parameter, 3 - not the
! inference 12. This is structure-driven supervised parameter
! identification through residual minimization: not gradient
! descent, not SGD, not backpropagation - and the Level-6 reverse
! structure J_Theta^T is deliberately not consumed by this solve.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program learning_level_9

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
  use relation_algorithms, only : topological_order
  use field_calculus, only : field
  use view_directed_stored    , only : stored_directed_graph
  use field_stored, only : stored_field
  use operation_gmres, only : gmres
  use learning_constitution_fixture, only : apply_law, slot_for_port
  use constituted_residual_fixture , only : constituted_learning_residual
  use graph_fractal        , only : graph, known_branch, null_branch
  use view_relational, only : relational_binding, &
       & num_member_sets, member_set_at, num_relations, relation_at, &
       & has_set

  implicit none


  integer, parameter :: ROW_R = 1

  type(graph)                  :: v, o, p, y
  type(graph)                   :: k, theta, u, p_in, p_out
  type(stored_relation), allocatable :: flow
  type(stored_relation)              :: located, consumes, produces
  class(relation), allocatable       :: d
  class(relation), pointer           :: gflow
  type(graph)             , target :: g
  type(graph)             , target :: scell(3), selem(3)
  type(graph)             , target :: rcell(2), relem(2)
  type(relational_binding)         :: bnd
  integer                          :: kcell
  type(stored_directed_graph)                 :: host
  type(constituted_learning_residual) :: residual_op
  type(gmres)                        :: solver
  type(stored_field)                        :: q_k, rhsf
  type(graph)     :: dom
  integer     :: n_dom
  class(field), allocatable    :: sol, rf
  integer, allocatable               :: order(:)
  real(dp), allocatable              :: obs(:), gv(:), solval(:), rv(:)
  real(dp)                           :: w_learned
  integer                            :: table(3, 6)
  integer                            :: nfail
  type(inclusion_map)     :: inclusions
  type(set_map)     :: sets

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "learning tower . level 9 . the statement"
  write(*,'(1x,a)') "============================================="

  ! -- the structure, and its one relational model graph
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
  allocate(flow)
  flow = stored_relation('flow', [o, v, p], table, sets)

  call p_in % declare()
  call sets       % bind(p_in, listed_set_representation([PORT_IN1, PORT_IN2]))
  call inclusions % include_in(p_in, p)
  call p_out % declare()
  call sets       % bind(p_out, listed_set_representation([PORT_OUT]))
  call inclusions % include_in(p_out, p)
  consumes = project_slots(restrict_slot(flow, 3, p_in , sets, inclusions), [2, 1], sets)
  produces = project_slots(restrict_slot(flow, 3, p_out, sets, inclusions), [1, 2], sets)
  d        = compose_binary(produces, consumes, sets)

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

  ! -- the discretization, distinct from the graph on purpose
  located = stored_relation('located', [y, v], &
       & reshape([ROW_R, SLOT_E], [2, 1]), sets)

  ! -- the execution order, from structure - never from the adapter
  call topological_order(d, sets, order)
  call report(size(order) .eq. 2 .and. &
       &      order(1) .eq. OP_PREDICT .and. order(2) .eq. OP_ERROR, &
       & "the derived order is [predict, error], exactly", nfail)

  ! -- the observation: Level 5's field on K, and NO field on U
  q_k = stored_field('observations', k, sets % num_members_of(k))
  call q_k % set_real_vector([6.0_dp, 2.0_dp])
  call q_k % real_vector(obs)

  ! -- the adapter keeps the GRAPH-OWNED flow; the selector dies.
  residual_op = constituted_learning_residual(g, bnd, flow, located, &
       & v, y, sets, k, obs, theta, u, order)
  deallocate(flow)

  ! -- the compatibility scenery: seven vertices, nobody's trainables
  host = stored_directed_graph(7, tails=[1,2,3,4,5,6], heads=[2,3,4,5,6,7])

  call check_structure(nfail, "before training")

  ! -- the ordinary solver, through its own operation face
  call solver % attach(residual_op, host, theta, sets % num_members_of(theta))
  solver % tolerance      = 1.0d-12
  solver % max_iterations = 50

  call solver % domain(host, dom, n_dom)
  call report(dom % same_as(theta), &
       & "the solver answers on Theta, by identity", nfail)
  call residual_op % domain(host, dom, n_dom)
  call report(dom % same_as(y) .and. .not. dom % same_as(theta), &
       & "the residual answers on Y - and Y is not Theta", nfail)

  ! -- the affine constant, through the ENTIRE constituted model
  call solver % constant(gv)
  call report(size(gv) .eq. 1 .and. abs(gv(1) + 6.0_dp) < 1.0d-12, &
       & "R(0) = -6: the full model evaluated at the untrained w", &
       & nfail)

  rhsf = stored_field('rhs', y, sets % num_members_of(y))
  call rhsf % set_real_vector(-gv)

  ! -- the solve: rhs on Y in, solution on Theta out; the internal
  !    state starts from zero, so w0 = 0 is the statement's own
  !    initial parameter - never injected.
  call solver % apply(host, [rhsf], sol)

  dom = sol % domain()
  call report(dom % same_as(theta), &
       & "the learned state is a field on Theta, by identity", nfail)

  call sol % real_vector(solval)
  w_learned = solval(sets % index_in(theta, SLOT_W))
  call report(abs(w_learned - 3.0_dp) < 1.0d-9, &
       & "learned w = 3, read through Theta's enumeration", nfail)
  call report(abs(w_learned) > 1.0_dp, &
       & "and the parameter MOVED: 3 is not the initial 0", nfail)

  ! -- the full symbolic model judges its own learned solution
  select type (sol)
  type is (stored_field)
     call residual_op % apply(host, [sol], rf)
  end select
  call rf % real_vector(rv)
  call report(sqrt(sum(rv * rv)) < 1.0d-9, &
       & "the constituted residual vanishes at the learned w", nfail)

  call check_structure(nfail, "after training")
  call check_inference(nfail)

  ! -- the tower's one primary result, computed - never a literal
  write(*,'(1x,a,i0)') "LEARNING_RESULT = ", nint(w_learned)

  call verdict(nfail, "level 9")

contains

  !===================================================================!
  ! Training changes a VALUE, never the topology: the graph-owned
  ! flow keeps its six facts, the derived dependency keeps its one,
  ! and the role subdomains keep their extensions - held both
  ! before the solve and after it.
  !===================================================================!

  subroutine check_structure(nfail, when)

    integer         , intent(inout) :: nfail
    character(len=*), intent(in)    :: when

    class(relation), pointer :: rp
    integer :: kk
    logical :: flow_ok, d_ok, roles_ok

    flow_ok = .false.
    d_ok    = .false.
    do kk = 1, num_relations(g)
       rp => relation_at(g, bnd, kk)
       if (rp % arity() .eq. 3) then
          gflow => rp
          flow_ok = rp % num_tuples() .eq. 6 .and. &
               & rp % has([OP_PREDICT, SLOT_W   , PORT_IN1]) .and. &
               & rp % has([OP_PREDICT, SLOT_X   , PORT_IN2]) .and. &
               & rp % has([OP_PREDICT, SLOT_YHAT, PORT_OUT]) .and. &
               & rp % has([OP_ERROR  , SLOT_YHAT, PORT_IN1]) .and. &
               & rp % has([OP_ERROR  , SLOT_Y   , PORT_IN2]) .and. &
               & rp % has([OP_ERROR  , SLOT_E   , PORT_OUT])
       else if (rp % arity() .eq. 2) then
          d_ok = rp % num_tuples() .eq. 1 .and. &
               & rp % has([OP_PREDICT, OP_ERROR])
       end if
    end do

    roles_ok = sets % num_members_of(theta) .eq. 1 .and. sets % has(theta, SLOT_W)  .and. &
         &     sets % num_members_of(k) .eq. 2 .and. sets % has(k, SLOT_Y)          .and. &
         &     sets % has(k, SLOT_X)                                  .and. &
         &     sets % num_members_of(u) .eq. 2 .and. sets % has(u, SLOT_E)          .and. &
         &     sets % has(u, SLOT_YHAT)

    call report(flow_ok .and. d_ok .and. roles_ok, &
         & "six facts, one dependency, three roles - intact " // when, &
         & nfail)

  end subroutine check_structure

  !===================================================================!
  ! Inference, separate from training: learned w plus a new x* = 4,
  ! and NO target - prediction needs no y. Structure still chooses
  ! the operands: the ports of predict are discovered from the
  ! graph-owned flow, the roles decide which value stands where,
  ! and the SAME Level-8 law computes. No retraining, no 3*4
  ! formula, no second constitution.
  !===================================================================!

  subroutine check_inference(nfail)

    integer, intent(inout) :: nfail

    integer  :: in1, in2, out1
    real(dp) :: a1, a2, prediction
    real(dp), parameter :: xstar = 4.0_dp

    in1  = slot_for_port(gflow, v, sets, OP_PREDICT, PORT_IN1)
    in2  = slot_for_port(gflow, v, sets, OP_PREDICT, PORT_IN2)
    out1 = slot_for_port(gflow, v, sets, OP_PREDICT, PORT_OUT)

    call report(in1 .eq. SLOT_W .and. in2 .eq. SLOT_X .and. &
         &      out1 .eq. SLOT_YHAT, &
         & "structure chooses the prediction operands: " // &
         & "w in1, x in2, yhat out", nfail)

    ! Positions dictated by the discovered ports; roles decide the
    ! values: the trainable slot carries the learned state, the
    ! other carries the fresh input.
    if (sets % has(theta, in1)) then
       a1 = w_learned
    else
       a1 = xstar
    end if
    if (sets % has(theta, in2)) then
       a2 = w_learned
    else
       a2 = xstar
    end if

    prediction = apply_law(OP_PREDICT, a1, a2)

    call report(abs(prediction - 12.0_dp) < 1.0d-9, &
         & "inference: x* = 4 gives yhat* = 12 - the same law, " // &
         & "no retraining, no target", nfail)

  end subroutine check_inference

end program learning_level_9
