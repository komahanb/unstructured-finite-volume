!=====================================================================!
! The graph algorithms suite: the laws of the four earned walks
! (level 4) over the directed adjacency interpretation, pinned
! generically, apart from the calculator that earned them.
!
! The domain here is built to EXPOSE numeric-order cheating: its
! declaration order is [3, 1, 5, 4, 2], so every canonical answer
! differs from ascending integers. The DAG:
!
!      3 --> 4        1 --> 4        4 --> 2        5 alone
!
! sources [3, 1, 5], sinks [5, 2] - in declaration order - and one
! topological order: [3, 1, 5, 4, 2]. An implementation that sorts
! by member value fails here.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_graph_algorithms

  use graph_fractal           , only : graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set           , only : set_map
  use map_label         , only : label_map
  use map_inclusion     , only : inclusion_map, declared_subobject
  use relation_binary, only : csr_relation
  use graph_algorithms     , only : sources, sinks, reachable, &
       &                            topological_order
  use graph_fractal        , only : graph, known_branch, null_branch
  use graph_relational_view, only : relational_binding, &
       & num_member_sets, member_set_at, num_relations, relation_at, &
       & holds_set

  implicit none


  type(graph)              :: ground
  type(graph)               :: v
  type(set_map)                    :: sets
  type(label_map)                  :: labels
  type(inclusion_map)              :: inclusions
  type(csr_relation)             :: after
  type(graph)             , target :: g
  type(graph)             , target :: scell(1), selem(1)
  type(graph)             , target :: rcell(1), relem(1)
  type(relational_binding)         :: bnd
  integer                          :: kcell
  integer                        :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "graph algorithms suite (level 4)"
  write(*,'(1x,a)') "============================================="

  call ground % declare()
  call sets % bind(ground, counted_set_representation(5))
  call labels % bind(ground, 'ground')
  call v % declare()
  call sets % bind(v, listed_set_representation([3, 1, 5, 4, 2]))
  call labels % bind(v, 'ordered-domain')
  call inclusions % include_in(v, ground)

  after = csr_relation('after', v, v, &
       & reshape([3,4,  1,4,  4,2], [2, 3]), sets)

  ! 'dag': (S, P) as one sequence on each branch.
  call g % declare()
  do kcell = 1, 1
     call scell(kcell) % declare()
     call selem(kcell) % declare()
  end do
  do kcell = 1, 1
     call rcell(kcell) % declare()
     call relem(kcell) % declare()
  end do

  call bnd % bind_set(selem(1), v)
  call bnd % bind_relation(relem(1), after)

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

  call check_sources_and_sinks(sets, labels, inclusions, nfail)
  call check_reachability(sets, nfail)
  call check_topological_order(sets, nfail)
  call check_tuple_order_invariance(v, sets, labels, inclusions, nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all graph algorithm checks passed"
  else
     write(*,'(1x,a,i0,a)') "FAILED: ", nfail, " algorithm check(s)"
     error stop
  end if

contains

  subroutine report(ok, label, nfail)

    logical         , intent(in)    :: ok
    character(len=*), intent(in)    :: label
    integer         , intent(inout) :: nfail

    if (ok) then
       write(*,'(1x,a,a)') "PASS : ", label
    else
       write(*,'(1x,a,a)') "FAIL : ", label
       nfail = nfail + 1
    end if

  end subroutine report

  !===================================================================!
  ! Sources and sinks come back as subobjects of the domain, in the
  ! DOMAIN'S declaration order - and the isolated member 5 sits
  ! lawfully in both.
  !===================================================================!

  subroutine check_sources_and_sinks(sets, labels, inclusions, nfail)

    type(set_map)                , intent(inout) :: sets
    type(label_map)              , intent(inout) :: labels
    type(inclusion_map)          , intent(inout) :: inclusions
    integer                      , intent(inout) :: nfail

    type(graph)          :: src, snk
    integer, allocatable :: idx(:)

    ! sources and sinks CARVE - a fresh declared set each - so the
    ! caller's maps are where they say what they are.
    call sources(after, sets, labels, inclusions, src)
    call sinks(after, sets, labels, inclusions, snk)

    call sets % members_of(src, idx)
    call report(size(idx) .eq. 3 .and. all(idx .eq. [3, 1, 5]), &
         & "sources are [3, 1, 5], by declaration, not by number", nfail)

    call sets % members_of(snk, idx)
    call report(size(idx) .eq. 2 .and. all(idx .eq. [5, 2]), &
         & "sinks are [5, 2], likewise", nfail)

    call report(sets % has_in(src, 5) .and. sets % has_in(snk, 5), &
         & "the isolated member is both source and sink", nfail)

    call report(declared_subobject(src, v, inclusions) .and. declared_subobject(snk, v, inclusions), &
         & "and both answers are subobjects of the domain", nfail)

  end subroutine check_sources_and_sinks

  !===================================================================!
  ! Reachability is path existence, with the zero-length path
  ! included and outsiders excluded.
  !===================================================================!

  subroutine check_reachability(sets, nfail)
    type(set_map)                , intent(in)    :: sets
    integer                      , intent(inout) :: nfail

    call report(reachable(after, sets, 3, 2), &
         & "3 reaches 2 through 4", nfail)
    call report(.not. reachable(after, sets, 2, 3), &
         & "and never backwards", nfail)
    call report(reachable(after, sets, 5, 5), &
         & "every member reaches itself by the zero-length path", nfail)
    call report(.not. reachable(after, sets, 5, 2), &
         & "the isolated member reaches nothing else", nfail)
    call report(.not. reachable(after, sets, 7, 2) .and. &
         &      .not. reachable(after, sets, 3, 0), &
         & "an outsider endpoint answers false, never a crash", nfail)

  end subroutine check_reachability

  !===================================================================!
  ! One graph, one order: the deterministic Kahn walk answers
  ! [3, 1, 5, 4, 2] - declaration-order ties, member values out.
  !===================================================================!

  subroutine check_topological_order(sets, nfail)
    type(set_map)                , intent(in)    :: sets
    integer                      , intent(inout) :: nfail

    integer, allocatable :: order(:)

    call topological_order(after, sets, order)

    call report(size(order) .eq. 5 .and. &
         &      all(order .eq. [3, 1, 5, 4, 2]), &
         & "the canonical order is [3, 1, 5, 4, 2]", nfail)

  end subroutine check_topological_order

  !===================================================================!
  ! The adjacency does not remember how its tuples were written:
  ! the same extension handed backwards answers the same sources,
  ! sinks, reachability and order.
  !===================================================================!

  subroutine check_tuple_order_invariance(v, sets, labels, inclusions, nfail)

    type(graph)        , intent(in)    :: v
    type(set_map)      , intent(inout) :: sets
    type(label_map)    , intent(inout) :: labels
    type(inclusion_map), intent(inout) :: inclusions
    integer            , intent(inout) :: nfail

    type(csr_relation)             :: backwards
    type(graph)             , target :: g2
    type(graph)             , target :: scell2(1), selem2(1)
    type(graph)             , target :: rcell2(1), relem2(1)
    type(relational_binding)         :: bnd2
    integer                          :: kcell2
    type(graph)               :: src, snk
    integer, allocatable           :: idx(:), order(:)
    logical                        :: ok

    backwards = csr_relation('after backwards', v, v, &
         & reshape([4,2,  1,4,  3,4], [2, 3]), sets)

    ! 'dag again': (S, P) as one sequence on each branch.
    call g2 % declare()
    do kcell2 = 1, 1
       call scell2(kcell2) % declare()
       call selem2(kcell2) % declare()
    end do
    do kcell2 = 1, 1
       call rcell2(kcell2) % declare()
       call relem2(kcell2) % declare()
    end do

    call bnd2 % bind_set(selem2(1), v)
    call bnd2 % bind_relation(relem2(1), backwards)

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

    g2 % branch(1) = known_branch(scell2(1))
    g2 % branch(2) = known_branch(rcell2(1))

    call sources(backwards, sets, labels, inclusions, src)
    call sinks(backwards, sets, labels, inclusions, snk)
    call sets % members_of(src, idx)
    ok = all(idx .eq. [3, 1, 5])
    call sets % members_of(snk, idx)
    ok = ok .and. all(idx .eq. [5, 2])
    call report(ok, &
         & "sources and sinks stand, tuples shuffled", nfail)

    call report(reachable(backwards, sets, 3, 2) .and. &
         &      .not. reachable(backwards, sets, 2, 3), &
         & "and so does reachability", nfail)

    call topological_order(backwards, sets, order)
    call report(all(order .eq. [3, 1, 5, 4, 2]), &
         & "and so does the one canonical order", nfail)

  end subroutine check_tuple_order_invariance

end program test_graph_algorithms
