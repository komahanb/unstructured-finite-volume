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
! written in code as compose_binary(J_Q^T, J_Q, sets); as a Boolean
! matrix pattern the same object reads J_Q^T J_Q. The two state
! slots share both residual rows, so each depends on the OTHER: the
! off-diagonal pair (u,v) and (v,u) is what makes the coupling
! mutual. Self-couplings are present too, but they are not the
! reason.
!
! Interpreted as a directed graph this is a perfectly VALID
! structure - the reading holds, the domain is Q, reachability
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
  use graph_fractal        , only : set_graph => graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set        , only : set_map
  use map_label      , only : label_map
  use map_inclusion  , only : inclusion_map, declared_subobject
  use relation_finitary, only : stored_relation, relation
  use relation_algebra, only : compose_binary
  use relation_binary , only : csr_relation, transposed_view, &
       &                             transpose_of, inclusion_of
  use graph_algorithms, only : reachable, sources, sinks
  use graph_fractal        , only : graph, known_branch, null_branch
  use graph_relational_view, only : relational_binding, &
       & num_member_sets, member_set_at, num_relations, relation_at, &
       & holds_set

  implicit none


  type(set_graph)              :: v, t
  type(set_graph)               :: p_dom, q_dom, y_dom, z_dom
  type(stored_relation)          :: dep
  type(csr_relation), target     :: inc_y, inc_q, jq
  type(transposed_view)          :: inc_q_t, jq_t
  type(csr_relation)             :: coupling
  type(graph)             , target :: g
  type(graph)             , target :: scell(1), selem(1)
  type(graph)             , target :: rcell(1), relem(1)
  type(relational_binding)         :: bnd
  integer                          :: kcell
  integer                        :: table(2, 9)
  integer                        :: nfail
  type(set_map)     :: sets
  type(inclusion_map)     :: inclusions
  type(label_map)     :: labels

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "adjoint tower . level 4 . implicit coupling"
  write(*,'(1x,a)') "============================================="

  call v % declare()
  call sets % bind(v, counted_set_representation(3))
  call t % declare()
  call sets % bind(t, counted_set_representation(3))

  call p_dom % declare()
  call sets       % bind(p_dom, listed_set_representation([VAR_P]))
  call inclusions % include_in(p_dom, v)
  call q_dom % declare()
  call sets       % bind(q_dom, listed_set_representation([VAR_U, VAR_V]))
  call inclusions % include_in(q_dom, v)
  call y_dom % declare()
  call sets       % bind(y_dom, listed_set_representation([TGT_R1, TGT_R2]))
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
  dep = stored_relation('dependency', [t, v], table, sets)

  inc_y   = inclusion_of(y_dom, t, sets, labels)
  inc_q   = inclusion_of(q_dom, v, sets, labels)
  inc_q_t = transpose_of(inc_q)
  jq      = compose_binary(compose_binary(inc_y, dep, sets), inc_q_t, sets)

  ! The coupling, derived from the SAME J_Q - state to residual row
  ! and back again. No second structure is authored for it.
  jq_t     = transpose_of(jq)
  coupling = compose_binary(jq_t, jq, sets)

  ! 'state coupling': (S, P) as one sequence on each branch.
  call g % declare()
  do kcell = 1, 1
     call scell(kcell) % declare()
     call selem(kcell) % declare()
  end do
  do kcell = 1, 1
     call rcell(kcell) % declare()
     call relem(kcell) % declare()
  end do

  call bnd % bind_set(selem(1), q_dom)
  call bnd % bind_relation(relem(1), coupling)

  do kcell = 1, 1
     scell(kcell) % branch(1) = known_branch(selem(kcell))
     if (kcell .lt. 1) scell(kcell) % branch(2) = &
          & known_branch(scell(kcell + 1))
  end do
  do kcell = 1, 1
     rcell(kcell) % branch(1) = known_branch(relem(kcell))
     if (kcell .lt. 1) rcell(kcell) % branch(2) = &
          & known_branch(rcell(kcell + 1))
  end do

  g % branch(1) = known_branch(scell(1))
  g % branch(2) = known_branch(rcell(1))

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

    type(set_graph) :: dom

    dom = coupling % domain(1)
    call report(dom % same_as(q_dom), &
         & "the coupling runs from the state slots", nfail)
    dom = coupling % domain(2)
    call report(dom % same_as(q_dom), &
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
  ! The interpretation is perfectly legitimate: the reading runs over
  ! the graph-owned coupling and walks the state domain.
  !===================================================================!

  subroutine check_view_is_valid(nfail)

    integer, intent(inout) :: nfail

    type(set_graph) :: dom

    dom = coupling % source()
    call report(dom % same_as(q_dom) .and. sets % size_of(dom) .eq. 2, &
         & "the coupling is valid and walks the state slots", nfail)

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
    call report(reachable(coupling, sets, VAR_U, VAR_V) .and. &
         &      reachable(coupling, sets, VAR_V, VAR_U), &
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

    type(set_graph) :: src, snk

    call sources(coupling, sets, labels, inclusions, src)
    call sinks(coupling, sets, labels, inclusions, snk)

    call report(sets % size_of(src) .eq. 0 .and. sets % size_of(snk) .eq. 0, &
         & "the state coupling has no source or sink state slot: " // &
         & "nowhere for a walk to begin", nfail)

  end subroutine check_no_source_or_sink

end program adjoint_level_4
