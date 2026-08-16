!=====================================================================!
! TIME INTEGRATION TOWER . LEVEL 4 . GRAPH CALCULUS
!
! The level answers one question: WHAT DOES CAUSAL TRAVERSAL OF THE
! TIME DEPENDENCY MEAN.
!
!      sources(A1)           = { t0 }
!      sinks(A1)             = { t4 }
!      reachable(t0, t4)     = true
!      reachable(t4, t0)     = false
!      topological_order(A1) = [ t0 t1 t2 t3 t4 ]
!
! That last list is the FORWARD CAUSAL ORDER, and it is worth being
! exact about where it comes from: the walk takes the first READY
! member in the carrier's own declaration order, through
! local_index, and never compares member values arithmetically.
! Time's order here is a consequence of structure, not of integers
! happening to ascend.
!
!                    CAUSALITY IS AN INTERPRETATION
!
! A1 is a relation; it does not know it is time. Nothing in Level 2
! called t0 "first". What makes this a CAUSAL order is the reading
! imposed here, by a profile and four algorithms that would say the
! same things about a dependency between calculator operations. The
! interpretation is the level's content, and it is separable from
! the structure it interprets.
!
!                    REVERSE CAUSAL ORDER
!
! [ t4 t3 t2 t1 t0 ] - stated as the OPPOSITE TRAVERSAL of the
! established forward order, and nothing more. No new algorithm is
! written for it. No adjoint is implemented. No derivative or
! adjoint machinery is imported, at this level or any level below.
!
! The observation is only that reverse causal order EXISTS
! STRUCTURALLY BEFORE AN ADJOINT EXISTS: whoever later builds one
! will find the ordering already here, and will not have to invent
! it alongside the mathematics.
!
!                    TWO VIEWS, ONE CARRIER
!
! A2 is interpreted the same way, over the SAME T. That is the
! Rosetta truth a later scheme can consume - one instant carrier,
! several structural readings of it - and it is checked here
! together with the refusal that keeps A2 honest: t0 does NOT reach
! t3 under two-step traversal.
!
! The views BORROW graph-owned storage. The selectors are destroyed
! immediately after the views are built, and every assertion below
! runs afterwards, so the borrow is demonstrated rather than
! described.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program time_level_4

  use time_assert           , only : report, verdict
  use time_assert           , only : NT
  use time_assert           , only : T0, T1, T2, T3, T4
  use graph_carrier         , only : counted_set, subset_set, member_set
  use graph_relation        , only : relation
  use graph_binary_relation , only : csr_relation
  use graph_profile         , only : directed_adjacency_view
  use graph_algorithms      , only : sources, sinks, reachable, &
       &                             topological_order
  use time_carriers_fixture , only : time_carriers
  use time_relations_fixture, only : tail_relation, head_relation
  use time_algebra_fixture  , only : derive_one_step_reach, &
       &                             derive_two_step_reach
  use fractal_graph        , only : graph, known_branch, null_branch
  use graph_relational_view, only : relational_binding, &
       & num_member_sets, member_set_at, num_relations, relation_at, &
       & holds_set

  implicit none


  type(counted_set)              :: q, t, e
  type(csr_relation), target     :: tail, head, a1
  type(csr_relation)             :: a2
  class(relation), allocatable   :: sel_a1, sel_a2
  type(graph)             , target :: g
  type(graph)             , target :: scell(3), selem(3)
  type(graph)             , target :: rcell(4), relem(4)
  type(relational_binding)         :: bnd
  integer                          :: kcell
  type(directed_adjacency_view)  :: view_a1, view_a2
  integer                        :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "time integration tower . level 4 . calculus"
  write(*,'(1x,a)') "============================================="

  call time_carriers(q, t, e)
  tail = tail_relation(e, t)
  head = head_relation(e, t)
  a1   = derive_one_step_reach(tail, head)
  a2   = derive_two_step_reach(a1)

  ! 'time': (S, P) as one sequence on each branch.
  call g % declare()
  do kcell = 1, 3
     call scell(kcell) % declare()
     call selem(kcell) % declare()
  end do
  do kcell = 1, 4
     call rcell(kcell) % declare()
     call relem(kcell) % declare()
  end do

  call bnd % bind_set(selem(1), q)
  call bnd % bind_set(selem(2), t)
  call bnd % bind_set(selem(3), e)
  call bnd % bind_relation(relem(1), tail)
  call bnd % bind_relation(relem(2), head)
  call bnd % bind_relation(relem(3), a1)
  call bnd % bind_relation(relem(4), a2)

  do kcell = 1, 3
     scell(kcell) % branch(1) = known_branch(selem(kcell))
     if (kcell .lt. 3) scell(kcell) % branch(2) = &
          & known_branch(scell(kcell + 1))
  end do
  do kcell = 1, 4
     rcell(kcell) % branch(1) = known_branch(relem(kcell))
     if (kcell .lt. 4) rcell(kcell) % branch(2) = &
          & known_branch(rcell(kcell + 1))
  end do

  g % branch(1) = known_branch(scell(1))
  g % branch(2) = known_branch(rcell(1))

  ! The interpretations read GRAPH-OWNED relations. A SELECTOR IS AN
  ! IDENTITY, NOT A THING: these two are separate objects that merely
  ! carry A1's and A2's identities, they are used only to find the
  ! graph's own storage, and they are destroyed before a single
  ! question is asked of the views.
  allocate(sel_a1, source=a1)
  allocate(sel_a2, source=a2)
  view_a1 = directed_adjacency_view(g, bnd, sel_a1)
  view_a2 = directed_adjacency_view(g, bnd, sel_a2)
  deallocate(sel_a1)
  deallocate(sel_a2)

  call check_causal_ends(nfail)
  call check_causal_reachability(nfail)
  call check_forward_causal_order(nfail)
  call check_two_step_view(nfail)
  call check_two_views_one_carrier(nfail)

  call verdict(nfail, "level 4")

contains

  !===================================================================!
  ! The ends of time, as subobjects of the instants - proved after
  ! the selectors' death, so the borrow demonstrably lives in the
  ! graph.
  !===================================================================!

  subroutine check_causal_ends(nfail)

    integer, intent(inout) :: nfail

    type(subset_set) :: src, snk

    src = sources(view_a1)
    snk = sinks(view_a1)

    call report(src % size() .eq. 1 .and. src % has(T0), &
         & "sources(A1) = { t0 }: one instant nothing leads to, the " // &
         & "selector long dead", nfail)
    call report(snk % size() .eq. 1 .and. snk % has(T4), &
         & "sinks(A1) = { t4 }: one instant leading nowhere", nfail)

    call report(src % is_subobject_of(t) .and. snk % is_subobject_of(t), &
         & "and both stand embedded in T - the ends of time are " // &
         & "instants, not a new kind of thing", nfail)

  end subroutine check_causal_ends

  !===================================================================!
  ! Time's asymmetry, measured: the first instant reaches the last,
  ! and the last reaches the first not at all.
  !===================================================================!

  subroutine check_causal_reachability(nfail)

    integer, intent(inout) :: nfail

    call report(reachable(view_a1, T0, T4), &
         & "t0 reaches t4 along the causal chain", nfail)
    call report(.not. reachable(view_a1, T4, T0), &
         & "and t4 reaches t0 not at all: THE ASYMMETRY IS " // &
         & "STRUCTURAL, not a convention about which way to loop", &
         & nfail)

    call report(reachable(view_a1, T2, T2), &
         & "every instant reaches itself by the zero-length path", &
         & nfail)
    call report(reachable(view_a1, T1, T3) .and. &
         &      .not. reachable(view_a1, T3, T1), &
         & "and the asymmetry holds in the interior too, not only " // &
         & "at the ends", nfail)

  end subroutine check_causal_reachability

  !===================================================================!
  ! THE rung's order. Checked against the CARRIER's declaration
  ! order through member(), never against a written [1,2,3,4,5]: the
  ! two coincide here, and a literal would hide the dependency that
  ! makes the coincidence meaningful.
  !===================================================================!

  subroutine check_forward_causal_order(nfail)

    integer, intent(inout) :: nfail

    integer, allocatable :: order(:)
    integer              :: i
    logical              :: ok

    call topological_order(view_a1, order)

    ok = size(order) .eq. NT
    do i = 1, min(size(order), t % size())
       ok = ok .and. (order(i) .eq. t % member(i))
    end do
    call report(ok, &
         & "topological_order(A1) = [t0 t1 t2 t3 t4], read from the " // &
         & "carrier's own declaration order - THE FORWARD CAUSAL " // &
         & "ORDER", nfail)

    ! The order is causal, not merely a permutation: every edge of
    ! A1 runs forwards in it.
    ok = .true.
    do i = 1, size(order) - 1
       ok = ok .and. .not. reachable(view_a1, order(i + 1), order(i))
    end do
    call report(ok, &
         & "and no later instant reaches an earlier one: the order " // &
         & "respects every dependency, which is what makes it causal", &
         & nfail)

    ! REVERSE CAUSAL ORDER, stated as the opposite traversal of the
    ! order just established. No new algorithm; no adjoint.
    ok = .true.
    do i = 1, size(order)
       ok = ok .and. (order(size(order) - i + 1) .eq. &
            &         t % member(t % size() - i + 1))
    end do
    call report(ok, &
         & "and read backwards it is [t4 t3 t2 t1 t0]: REVERSE " // &
         & "CAUSAL ORDER EXISTS STRUCTURALLY BEFORE ANY ADJOINT " // &
         & "DOES - no algorithm was written for it here", nfail)

  end subroutine check_forward_causal_order

  !===================================================================!
  ! A2 read the same way, and the refusal that keeps it honest.
  !===================================================================!

  subroutine check_two_step_view(nfail)

    integer, intent(inout) :: nfail

    integer, pointer :: fibre(:)
    logical          :: ok

    fibre => view_a2 % successors_view(T0)
    ok = size(fibre) .eq. 1
    if (ok) ok = fibre(1) .eq. T2
    fibre => view_a2 % successors_view(T1)
    ok = ok .and. size(fibre) .eq. 1
    if (ok) ok = ok .and. fibre(1) .eq. T3
    fibre => view_a2 % successors_view(T2)
    ok = ok .and. size(fibre) .eq. 1
    if (ok) ok = ok .and. fibre(1) .eq. T4
    call report(ok, &
         & "successors under A2: t0->{t2}, t1->{t3}, t2->{t4}", nfail)

    call report(reachable(view_a2, T0, T4), &
         & "t0 reaches t4 under repeated two-step traversal, by way " // &
         & "of t2", nfail)

    call report(.not. reachable(view_a2, T0, T3), &
         & "but t0 does NOT reach t3: two-step reach lands only on " // &
         & "even offsets - REACH IS NOT A SCHEME, and a BDF2 " // &
         & "dependency would have had to see t3", nfail)

  end subroutine check_two_step_view

  !===================================================================!
  ! THE Rosetta truth a later scheme can consume: not two carriers
  ! and not two graphs, but two readings of ONE instant carrier.
  !===================================================================!

  subroutine check_two_views_one_carrier(nfail)

    integer, intent(inout) :: nfail

    class(member_set), allocatable :: d1, d2

    d1 = view_a1 % domain()
    d2 = view_a2 % domain()

    call report(d1 % same_as(t) .and. d2 % same_as(t), &
         & "both views run over T itself, by identity", nfail)

    call report(d1 % same_as(d2), &
         & "so A1 and A2 are DIFFERENT STRUCTURAL VIEWS OVER THE " // &
         & "SAME CARRIER - one time axis, several readings of it", &
         & nfail)

    call report(.not. d1 % same_as(q) .and. .not. d1 % same_as(e), &
         & "and neither reading has quietly become the state axis " // &
         & "or the steps", nfail)

  end subroutine check_two_views_one_carrier

end program time_level_4
