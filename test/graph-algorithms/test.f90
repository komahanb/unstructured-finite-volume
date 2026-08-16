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

  use graph_carrier        , only : counted_set, subset_set, member_set
  use graph_binary_relation, only : csr_relation
  use graph_structure      , only : relational_graph, held_set, &
       &                            held_relation
  use graph_profile        , only : directed_adjacency_view
  use graph_algorithms     , only : sources, sinks, reachable, &
       &                            topological_order
  use fractal_graph        , only : graph
  use graph_relational_view, only : relational_binding
  use relational_fixture   , only : fractal_fixture

  implicit none

  type(fractal_fixture)             :: fx_
  type(graph)             , pointer :: fg_
  type(relational_binding), pointer :: fb_

  type(counted_set)              :: ground
  type(subset_set)               :: v
  type(csr_relation)             :: after
  type(relational_graph), target :: g
  type(directed_adjacency_view)  :: view
  integer                        :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "graph algorithms suite (level 4)"
  write(*,'(1x,a)') "============================================="

  ground = counted_set('ground', 5)
  v      = subset_set('ordered-domain', ground, [3, 1, 5, 4, 2])

  after = csr_relation('after', v, v, &
       & reshape([3,4,  1,4,  4,2], [2, 3]))

  g = relational_graph('dag', [held_set(v)], [held_relation(after)])

  call fx_ % to_fractal(g, fg_, fb_)
  view = directed_adjacency_view(fg_, fb_, after)

  call check_sources_and_sinks(view, nfail)
  call check_reachability(view, nfail)
  call check_topological_order(view, nfail)
  call check_tuple_order_invariance(nfail)

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

  subroutine check_sources_and_sinks(view, nfail)

    type(directed_adjacency_view), intent(in)    :: view
    integer                      , intent(inout) :: nfail

    type(subset_set)     :: src, snk
    integer, allocatable :: idx(:)

    src = sources(view)
    snk = sinks(view)

    call src % members(idx)
    call report(size(idx) .eq. 3 .and. all(idx .eq. [3, 1, 5]), &
         & "sources are [3, 1, 5], by declaration, not by number", nfail)

    call snk % members(idx)
    call report(size(idx) .eq. 2 .and. all(idx .eq. [5, 2]), &
         & "sinks are [5, 2], likewise", nfail)

    call report(src % has(5) .and. snk % has(5), &
         & "the isolated member is both source and sink", nfail)

    call report(src % is_subobject_of(v) .and. snk % is_subobject_of(v), &
         & "and both answers are subobjects of the domain", nfail)

  end subroutine check_sources_and_sinks

  !===================================================================!
  ! Reachability is path existence, with the zero-length path
  ! included and outsiders excluded.
  !===================================================================!

  subroutine check_reachability(view, nfail)

    type(directed_adjacency_view), intent(in)    :: view
    integer                      , intent(inout) :: nfail

    call report(reachable(view, 3, 2), &
         & "3 reaches 2 through 4", nfail)
    call report(.not. reachable(view, 2, 3), &
         & "and never backwards", nfail)
    call report(reachable(view, 5, 5), &
         & "every member reaches itself by the zero-length path", nfail)
    call report(.not. reachable(view, 5, 2), &
         & "the isolated member reaches nothing else", nfail)
    call report(.not. reachable(view, 7, 2) .and. &
         &      .not. reachable(view, 3, 0), &
         & "an outsider endpoint answers false, never a crash", nfail)

  end subroutine check_reachability

  !===================================================================!
  ! One graph, one order: the deterministic Kahn walk answers
  ! [3, 1, 5, 4, 2] - declaration-order ties, member values out.
  !===================================================================!

  subroutine check_topological_order(view, nfail)

    type(directed_adjacency_view), intent(in)    :: view
    integer                      , intent(inout) :: nfail

    integer, allocatable :: order(:)

    call topological_order(view, order)

    call report(size(order) .eq. 5 .and. &
         &      all(order .eq. [3, 1, 5, 4, 2]), &
         & "the canonical order is [3, 1, 5, 4, 2]", nfail)

  end subroutine check_topological_order

  !===================================================================!
  ! The adjacency does not remember how its tuples were written:
  ! the same extension handed backwards answers the same sources,
  ! sinks, reachability and order.
  !===================================================================!

  subroutine check_tuple_order_invariance(nfail)

    integer, intent(inout) :: nfail

    type(csr_relation)             :: backwards
    type(relational_graph), target :: g2
    type(directed_adjacency_view)  :: view2
    type(subset_set)               :: src, snk
    integer, allocatable           :: idx(:), order(:)
    logical                        :: ok

    backwards = csr_relation('after backwards', v, v, &
         & reshape([4,2,  1,4,  3,4], [2, 3]))

    g2 = relational_graph('dag again', [held_set(v)], &
         & [held_relation(backwards)])
    call fx_ % to_fractal(g2, fg_, fb_)
    view2 = directed_adjacency_view(fg_, fb_, backwards)

    src = sources(view2)
    snk = sinks(view2)
    call src % members(idx)
    ok = all(idx .eq. [3, 1, 5])
    call snk % members(idx)
    ok = ok .and. all(idx .eq. [5, 2])
    call report(ok, &
         & "sources and sinks stand, tuples shuffled", nfail)

    call report(reachable(view2, 3, 2) .and. &
         &      .not. reachable(view2, 2, 3), &
         & "and so does reachability", nfail)

    call topological_order(view2, order)
    call report(all(order .eq. [3, 1, 5, 4, 2]), &
         & "and so does the one canonical order", nfail)

  end subroutine check_tuple_order_invariance

end program test_graph_algorithms
