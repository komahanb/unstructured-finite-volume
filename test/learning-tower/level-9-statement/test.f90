!=====================================================================!
! LEARNING TOWER . LEVEL 9 . THE STATEMENT
!
! The level asks the only question left: WHAT LEARNING PROBLEM IS
! BEING ASKED. The statement is: given the observation (x, y) =
! (2, 6) and the model structure predict -> yhat, error -> e with
! the Level-8 constitution, learn w from w0 = 0 such that r = 0.
! It SELECTS - it invents nothing:
!
!      structure       R_flow; D DERIVED; order DERIVED
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
  use graph_set  , only : index_set, subset, set
  use graph_relation , only : stored_relation, relation
  use graph_relation_algebra, only : restrict_slot, project_slots, &
       &                             compose_binary
  use graph_structure, only : related_graph, declared_set, declared_relation
  use graph_interpretation  , only : directed_adjacency_view
  use graph_algorithms, only : topological_order
  use graph_grammar  , only : graph_field
  use class_graph    , only : stored_graph
  use class_graph_field, only : field
  use class_graph_gmres, only : gmres
  use learning_constitution_fixture, only : apply_law, slot_for_port
  use constituted_residual_fixture , only : constituted_learning_residual

  implicit none

  integer, parameter :: ROW_R = 1

  type(index_set)                  :: v, o, p, y
  type(subset)                   :: k, theta, u, p_in, p_out
  type(stored_relation), allocatable :: flow
  type(stored_relation)              :: located, consumes, produces
  class(relation), allocatable       :: d
  class(relation), pointer           :: gflow
  type(related_graph), target     :: g
  type(directed_adjacency_view)      :: view
  type(stored_graph)                 :: host
  type(constituted_learning_residual) :: residual_op
  type(gmres)                        :: solver
  type(field)                        :: q_k, rhsf
  class(set), allocatable     :: dom
  class(graph_field), allocatable    :: sol, rf
  integer, allocatable               :: order(:)
  real(dp), allocatable              :: obs(:), gv(:), solval(:), rv(:)
  real(dp)                           :: w_learned
  integer                            :: table(3, 6)
  integer                            :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "learning tower . level 9 . the statement"
  write(*,'(1x,a)') "============================================="

  ! -- the structure, and its one relational model graph
  v     = index_set('value-slots'  , 5)
  o     = index_set('operations'   , 2)
  p     = index_set('ports'        , 3)
  y     = index_set('residual-rows', 1)
  k     = subset('observed' , v, [SLOT_Y, SLOT_X])
  theta = subset('trainable', v, [SLOT_W])
  u     = subset('computed' , v, [SLOT_E, SLOT_YHAT])

  table(:, 1) = [OP_PREDICT, SLOT_W   , PORT_IN1]
  table(:, 2) = [OP_PREDICT, SLOT_X   , PORT_IN2]
  table(:, 3) = [OP_PREDICT, SLOT_YHAT, PORT_OUT]
  table(:, 4) = [OP_ERROR  , SLOT_YHAT, PORT_IN1]
  table(:, 5) = [OP_ERROR  , SLOT_Y   , PORT_IN2]
  table(:, 6) = [OP_ERROR  , SLOT_E   , PORT_OUT]
  allocate(flow)
  flow = stored_relation('flow', [o, v, p], table)

  p_in     = subset('input-ports', p, [PORT_IN1, PORT_IN2])
  p_out    = subset('output-port', p, [PORT_OUT])
  consumes = project_slots(restrict_slot(flow, 3, p_in ), [2, 1])
  produces = project_slots(restrict_slot(flow, 3, p_out), [1, 2])
  d        = compose_binary(produces, consumes)

  g = related_graph('learning', &
       & [declared_set(v), declared_set(o), declared_set(p)], &
       & [declared_relation(flow), declared_relation(d)])

  ! -- the discretization, distinct from the graph on purpose
  located = stored_relation('located', [y, v], &
       & reshape([ROW_R, SLOT_E], [2, 1]))

  ! -- the execution order, from structure - never from the adapter
  view = directed_adjacency_view(g, d)
  call topological_order(view, order)
  call report(size(order) .eq. 2 .and. &
       &      order(1) .eq. OP_PREDICT .and. order(2) .eq. OP_ERROR, &
       & "the derived order is [predict, error], exactly", nfail)

  ! -- the observation: Level 5's field on K, and NO field on U
  q_k = field('observations', k)
  call q_k % set_real_vector([6.0_dp, 2.0_dp])
  call q_k % get_real_vector(obs)

  ! -- the adapter keeps the GRAPH-OWNED flow; the selector dies.
  residual_op = constituted_learning_residual(g, flow, located, &
       & v, y, k, obs, theta, u, order)
  deallocate(flow)

  ! -- the compatibility scenery: seven vertices, nobody's trainables
  host = stored_graph(7, tails=[1,2,3,4,5,6], heads=[2,3,4,5,6,7])

  call check_structure(nfail, "before training")

  ! -- the ordinary solver, through its own operation face
  call solver % attach(residual_op, host, theta)
  solver % tolerance      = 1.0d-12
  solver % max_iterations = 50

  call solver % domain(host, dom)
  call report(dom % equals(theta), &
       & "the solver answers on Theta, by identity", nfail)
  call residual_op % domain(host, dom)
  call report(dom % equals(y) .and. .not. dom % equals(theta), &
       & "the residual answers on Y - and Y is not Theta", nfail)

  ! -- the affine constant, through the ENTIRE constituted model
  call solver % constant(gv)
  call report(size(gv) .eq. 1 .and. abs(gv(1) + 6.0_dp) < 1.0d-12, &
       & "R(0) = -6: the full model evaluated at the untrained w", &
       & nfail)

  rhsf = field('rhs', y)
  call rhsf % set_real_vector(-gv)

  ! -- the solve: rhs on Y in, solution on Theta out; the internal
  !    state starts from zero, so w0 = 0 is the statement's own
  !    initial parameter - never injected.
  call solver % apply(host, [rhsf], sol)

  call sol % domain(dom)
  call report(dom % equals(theta), &
       & "the learned state is a field on Theta, by identity", nfail)

  call sol % get_real_vector(solval)
  w_learned = solval(theta % local_index(SLOT_W))
  call report(abs(w_learned - 3.0_dp) < 1.0d-9, &
       & "learned w = 3, read through Theta's enumeration", nfail)
  call report(abs(w_learned) > 1.0_dp, &
       & "and the parameter MOVED: 3 is not the initial 0", nfail)

  ! -- the full symbolic model judges its own learned solution
  select type (sol)
  type is (field)
     call residual_op % apply(host, [sol], rf)
  end select
  call rf % get_real_vector(rv)
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
    do kk = 1, g % num_relations()
       rp => g % relation_at(kk)
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

    roles_ok = theta % size() .eq. 1 .and. theta % has(SLOT_W)  .and. &
         &     k % size() .eq. 2 .and. k % has(SLOT_Y)          .and. &
         &     k % has(SLOT_X)                                  .and. &
         &     u % size() .eq. 2 .and. u % has(SLOT_E)          .and. &
         &     u % has(SLOT_YHAT)

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

    in1  = slot_for_port(gflow, v, OP_PREDICT, PORT_IN1)
    in2  = slot_for_port(gflow, v, OP_PREDICT, PORT_IN2)
    out1 = slot_for_port(gflow, v, OP_PREDICT, PORT_OUT)

    call report(in1 .eq. SLOT_W .and. in2 .eq. SLOT_X .and. &
         &      out1 .eq. SLOT_YHAT, &
         & "structure chooses the prediction operands: " // &
         & "w in1, x in2, yhat out", nfail)

    ! Positions dictated by the discovered ports; roles decide the
    ! values: the trainable slot carries the learned state, the
    ! other carries the fresh input.
    if (theta % has(in1)) then
       a1 = w_learned
    else
       a1 = xstar
    end if
    if (theta % has(in2)) then
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
