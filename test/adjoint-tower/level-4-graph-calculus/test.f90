!=====================================================================!
! ADJOINT TOWER . LEVEL 4 . THE IMPLICIT SYSTEM IS NOT A PROGRAM
!
! The level answers one question: WHAT DOES GRAPH CALCULUS HONESTLY
! SAY ABOUT AN IMPLICIT SYSTEM. Its siblings all found an execution
! order here - the calculator walked [plus, times], the learner
! walked [predict, error], the derivative specimen walked
! [product, sum]. This tower must NOT manufacture one.
!
! The state is implicitly coupled, and the coupling is derived from
! the same J_Q, never authored:
!
!      C_Q = J_Q o J_Q^T  <=  Q x Q            the path Q -> Y -> Q
!          = { (u,u), (u,v), (v,u), (v,v) }
!
! written in code as compose_binary(J_Q^T, J_Q); as a Boolean
! matrix pattern the same object reads J_Q^T J_Q. The two state
! slots share both residual rows, so each depends on the OTHER: the
! off-diagonal pair (u,v) and (v,u) is what makes the coupling
! mutual. Self-couplings are present too, but they are not the
! reason.
!
! Interpreted as a directed graph this is a perfectly VALID
! structure - the view builds, the domain is Q, reachability
! answers in both directions - and it is emphatically not a DAG:
!
!      a valid directed graph  /=  a DAG
!
! A topological order therefore does not exist, and the nucleus
! says so by REFUSING rather than inventing one; that refusal is
! this rung's central truth and it is checked, by message, in
! check_refusals.sh. What solves an implicit system is minimization
! (Level 7), never a walk. No acyclic fiction was invented here to
! make graph calculus produce an order it has no right to produce.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program adjoint_level_4

  use adjoint_assert, only : report, verdict
  use adjoint_assert, only : VAR_P, VAR_U, VAR_V
  use adjoint_assert, only : TGT_R1, TGT_R2, TGT_F
  use graph_set , only : index_set, subset, set
  use graph_relation, only : stored_relation, relation
  use graph_relation_algebra, only : compose_binary
  use graph_binary_relation , only : csr_relation, transposed_view, &
       &                             transpose_of, inclusion_of
  use graph_structure, only : related_graph, declared_set, declared_relation
  use graph_interpretation  , only : directed_adjacency_view
  use graph_algorithms, only : reachable, sources, sinks

  implicit none

  type(index_set)              :: v, t
  type(subset)               :: p_dom, q_dom, y_dom, z_dom
  type(stored_relation)          :: dep
  type(csr_relation), target     :: inc_y, inc_q, jq
  type(transposed_view)          :: inc_q_t, jq_t
  type(csr_relation)             :: coupling
  type(related_graph), target :: g
  type(directed_adjacency_view)  :: view
  integer                        :: table(2, 9)
  integer                        :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "adjoint tower . level 4 . implicit coupling"
  write(*,'(1x,a)') "============================================="

  v = index_set('variables', 3)
  t = index_set('targets'  , 3)

  p_dom = subset('parameter', v, [VAR_P])
  q_dom = subset('state'    , v, [VAR_U, VAR_V])
  y_dom = subset('residual' , t, [TGT_R1, TGT_R2])
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
  dep = stored_relation('dependency', [t, v], table)

  inc_y   = inclusion_of(y_dom)
  inc_q   = inclusion_of(q_dom)
  inc_q_t = transpose_of(inc_q)
  jq      = compose_binary(compose_binary(inc_y, dep), inc_q_t)

  ! The coupling, derived from the SAME J_Q - state to residual row
  ! and back again. No second structure is authored for it.
  jq_t     = transpose_of(jq)
  coupling = compose_binary(jq_t, jq)

  g = related_graph('state coupling', [declared_set(q_dom)], &
       & [declared_relation(coupling)])
  view = directed_adjacency_view(g, coupling)

  call check_coupling_extension(nfail)
  call check_view_is_valid(nfail)
  call check_the_cycle(nfail)
  call check_no_source_or_sink(nfail)

  ! topological_order is NOT called here: it would - correctly - stop
  ! the program. Its refusal is the rung's truth, and it is proved
  ! in refusal.f90 where dying loudly is the expected outcome.

  call verdict(nfail, "level 4")

contains

  !===================================================================!
  ! The coupling, derived and exact: both state slots share both
  ! residual rows, so all four pairs stand - including the two
  ! self-couplings.
  !===================================================================!

  subroutine check_coupling_extension(nfail)

    integer, intent(inout) :: nfail

    class(set), allocatable :: dom

    dom = coupling % domain(1)
    call report(dom % equals(q_dom), &
         & "the coupling runs from the state slots", nfail)
    dom = coupling % domain(2)
    call report(dom % equals(q_dom), &
         & "back into the state slots", nfail)

    call report(coupling % num_tuples() .eq. 4 .and. &
         &      coupling % has([VAR_U, VAR_U]) .and. &
         &      coupling % has([VAR_U, VAR_V]) .and. &
         &      coupling % has([VAR_V, VAR_U]) .and. &
         &      coupling % has([VAR_V, VAR_V]), &
         & "C_Q = { (u,u), (u,v), (v,u), (v,v) } - derived from " // &
         & "J_Q alone", nfail)

  end subroutine check_coupling_extension

  !===================================================================!
  ! The interpretation is perfectly legitimate: a view builds over
  ! the graph-owned coupling and walks the state domain.
  !===================================================================!

  subroutine check_view_is_valid(nfail)

    integer, intent(inout) :: nfail

    class(set), allocatable :: dom

    dom = view % domain()
    call report(dom % equals(q_dom) .and. dom % size() .eq. 2, &
         & "the view is valid and walks the state slots", nfail)

  end subroutine check_view_is_valid

  !===================================================================!
  ! And it holds a genuine cycle. The MUTUAL coupling is the
  ! off-diagonal pair: u depends on v and v depends on u, because
  ! they share residual rows. The self-couplings exist as stored
  ! facts too - worth pinning, since a walk would trip on them
  ! either way - but they are not why the state is coupled.
  !===================================================================!

  subroutine check_the_cycle(nfail)

    integer, intent(inout) :: nfail

    call report(coupling % has([VAR_U, VAR_V]) .and. &
         &      coupling % has([VAR_V, VAR_U]), &
         & "the off-diagonal pair stands: u depends on v and v on " // &
         & "u - THAT is the mutual coupling", nfail)
    call report(reachable(view, VAR_U, VAR_V) .and. &
         &      reachable(view, VAR_V, VAR_U), &
         & "and each reaches the other, so the walk finds a cycle", &
         & nfail)
    call report(coupling % has([VAR_U, VAR_U]) .and. &
         &      coupling % has([VAR_V, VAR_V]), &
         & "self-couplings are present too, by stored fact - not " // &
         & "the reason, but not fiction either", nfail)

  end subroutine check_the_cycle

  !===================================================================!
  ! The structural signature of an implicit system: the derived
  ! state coupling has no source and no sink STATE SLOT - C_Q
  ! relates slots, not equations - so there is nowhere for a walk
  ! to begin or end, which is exactly why a walk cannot solve it.
  !===================================================================!

  subroutine check_no_source_or_sink(nfail)

    integer, intent(inout) :: nfail

    type(subset) :: src, snk

    src = sources(view)
    snk = sinks(view)

    call report(src % size() .eq. 0 .and. snk % size() .eq. 0, &
         & "the state coupling has no source or sink state slot: " // &
         & "nowhere for a walk to begin", nfail)

  end subroutine check_no_source_or_sink

end program adjoint_level_4
