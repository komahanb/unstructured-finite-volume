!=====================================================================!
! The ordinary profile suite (AGENTS.md, level 4, phase 4B).
!
! The ordinary directed graph is one schema over the relational
! container: T <= E x V total, H <= E x V partial, the boundary an
! absence in H. This suite's centrepiece is the compatibility law
! (AGENTS.md 62): on every topology below, the profile - deriving
! every answer from T and H alone - is held against the old
! stored_graph, query for query, vertex for vertex, edge for edge.
! One source of truth on the new road, and the old road as its
! oracle until the day it retires.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_graph_ordinary

  use class_graph          , only : stored_graph
  use graph_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use graph_set_map           , only : set_map
  use graph_binary_relation, only : csr_relation
  use graph_profile        , only : ordinary_graph_view
  use fractal_graph        , only : graph, known_branch
  use graph_relational_view, only : relational_binding

  implicit none

  integer :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "graph ordinary profile suite (AGENTS phase 4B)"
  write(*,'(1x,a)') "============================================="

  ! A diamond with a wall and an island: 1->2->4, 1->3->4, 4->wall,
  ! vertex 5 alone.
  call compare_topology(5, [1, 1, 2, 3, 4], [2, 3, 4, 4, 0], &
       & "diamond with a wall and an island", nfail)

  ! Two edges joining one pair.
  call compare_topology(2, [1, 1], [2, 2], &
       & "a parallel pair", nfail)

  ! The chain with a wall at its end.
  call compare_topology(3, [1, 2, 3], [2, 3, 0], &
       & "a chain into a wall", nfail)

  ! A self-loop beside ordinary edges: 1->2, 2->2, 2->3, 3->wall.
  ! Old law: incident counts the loop twice, adjacency never names
  ! the vertex to itself.
  call compare_topology(3, [1, 2, 2, 3], [2, 2, 3, 0], &
       & "a self-loop beside a wall", nfail)

  ! A relation is a set: hand the same topologies in scrambled tuple
  ! order and every answer must stand.
  call compare_topology(5, [1, 1, 2, 3, 4], [2, 3, 4, 4, 0], &
       & "the diamond, tuples shuffled", nfail, scrambled=.true.)
  call compare_topology(3, [1, 2, 2, 3], [2, 2, 3, 0], &
       & "the self-loop, tuples shuffled", nfail, scrambled=.true.)

  call check_wall_is_absence(nfail)
  call check_declared_order(nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all ordinary profile checks passed"
  else
     write(*,'(1x,a,i0,a)') "FAILED: ", nfail, " ordinary profile check(s)"
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
  ! The oracle bout: build the old stored_graph and the relational
  ! road on one topology, then ask both everything and demand one
  ! answer. T gets every edge; H gets the headed ones; both tables
  ! in ascending edge order, as every mesh builder hands them.
  !===================================================================!

  subroutine compare_topology(nv, tails, heads, what, nfail, scrambled)

    integer         , intent(in)           :: nv
    integer         , intent(in)           :: tails(:), heads(:)
    character(len=*), intent(in)           :: what
    integer         , intent(inout)        :: nfail
    logical         , intent(in), optional :: scrambled

    type(stored_graph)              :: old
    type(graph)               :: verts, edges
    type(set_map)               :: sets
    type(csr_relation)              :: t, h
    type(graph)            , target :: g
    type(graph)            , target :: scell(2), selem(2), rcell(2), relem(2)
    type(relational_binding)        :: bnd
    type(ordinary_graph_view)       :: view
    integer, allocatable            :: ttab(:,:), htab(:,:)
    integer, allocatable            :: a(:), b(:)
    integer                         :: ne, nh, e, v, q, k
    logical                         :: ok

    ne = size(tails)

    old = stored_graph(nv, tails=tails, heads=heads)

    call verts % declare()
    call sets % bind(verts, counted_set_representation(nv))
    call edges % declare()
    call sets % bind(edges, counted_set_representation(ne))

    allocate(ttab(2, ne))
    nh = 0
    do e = 1, ne
       ttab(:, e) = [e, tails(e)]
       if (heads(e) >= 1) nh = nh + 1
    end do
    allocate(htab(2, nh))
    nh = 0
    do e = 1, ne
       if (heads(e) >= 1) then
          nh = nh + 1
          htab(:, nh) = [e, heads(e)]
       end if
    end do

    ! A set does not remember the order it was written down in:
    ! reverse the columns and nothing downstream may notice.
    if (present(scrambled)) then
       if (scrambled) then
          ttab = ttab(:, size(ttab, 2) : 1 : -1)
          htab = htab(:, size(htab, 2) : 1 : -1)
       end if
    end if

    t = csr_relation('tail', edges, verts, ttab, sets)
    h = csr_relation('head', edges, verts, htab, sets)

    ! (S, P) = ({V, E}, {T, H}), one sequence on each branch.
    call g % declare()
    do k = 1, 2
       call scell(k) % declare(); call selem(k) % declare()
       call rcell(k) % declare(); call relem(k) % declare()
    end do

    call bnd % bind_set(selem(1), verts)
    call bnd % bind_set(selem(2), edges)
    call bnd % bind_relation(relem(1), t)
    call bnd % bind_relation(relem(2), h)

    do k = 1, 2
       scell(k) % branch(1) = known_branch(selem(k))
       rcell(k) % branch(1) = known_branch(relem(k))
       if (k .lt. 2) then
          scell(k) % branch(2) = known_branch(scell(k + 1))
          rcell(k) % branch(2) = known_branch(rcell(k + 1))
       end if
    end do

    g % branch(1) = known_branch(scell(1))
    g % branch(2) = known_branch(rcell(1))

    view = ordinary_graph_view(g, bnd, sets, tail_at=1, head_at=2)

    ok = (view % num_vertices() .eq. old % num_vertices()) .and. &
         & (view % num_edges() .eq. old % num_edges())

    do e = 1, ne
       ok = ok .and. (view % edge_tail(e) .eq. old % edge_tail(e))
       ok = ok .and. (view % edge_head(e) .eq. old % edge_head(e))
       ok = ok .and. (view % edge_has_head(e) .eqv. old % edge_has_head(e))
    end do

    do v = 1, nv
       do q = 1, 6
          select case (q)
          case (1)
             call view % outgoing_edges(v, a)
             call old % outgoing_edges(v, b)
          case (2)
             call view % incoming_edges(v, a)
             call old % incoming_edges(v, b)
          case (3)
             call view % incident_edges(v, a)
             call old % incident_edges(v, b)
          case (4)
             call view % adjacent_vertices(v, a)
             call old % adjacent_vertices(v, b)
          case (5)
             call view % outgoing_vertices(v, a)
             call old % outgoing_vertices(v, b)
          case (6)
             call view % incoming_vertices(v, a)
             call old % incoming_vertices(v, b)
          end select
          ok = ok .and. (size(a) .eq. size(b))
          if (size(a) .eq. size(b)) ok = ok .and. all(a .eq. b)
       end do
    end do

    call report(ok, "the profile matches the old graph on " // what, nfail)

  end subroutine compare_topology

  !===================================================================!
  ! The boundary law, said directly: the wall edge stands in T and
  ! is an absence in H - no imaginary member, no ghost tuple.
  !===================================================================!

  subroutine check_wall_is_absence(nfail)

    integer, intent(inout) :: nfail

    type(graph)               :: verts, edges
    type(set_map)                   :: sets
    type(csr_relation)              :: t, h
    type(graph)            , target :: g
    type(graph)            , target :: scell(2), selem(2), rcell(2), relem(2)
    type(relational_binding)        :: bnd
    type(ordinary_graph_view)       :: view
    integer                         :: v, k
    logical                         :: ok

    call verts % declare()
    call sets % bind(verts, counted_set_representation(3))
    call edges % declare()
    call sets % bind(edges, counted_set_representation(3))

    t = csr_relation('tail', edges, verts, &
         & reshape([1,1,  2,2,  3,3], [2, 3]), sets)
    h = csr_relation('head', edges, verts, &
         & reshape([1,2,  2,3], [2, 2]), sets)

    call g % declare()
    do k = 1, 2
       call scell(k) % declare(); call selem(k) % declare()
       call rcell(k) % declare(); call relem(k) % declare()
    end do

    call bnd % bind_set(selem(1), verts)
    call bnd % bind_set(selem(2), edges)
    call bnd % bind_relation(relem(1), t)
    call bnd % bind_relation(relem(2), h)

    do k = 1, 2
       scell(k) % branch(1) = known_branch(selem(k))
       rcell(k) % branch(1) = known_branch(relem(k))
       if (k .lt. 2) then
          scell(k) % branch(2) = known_branch(scell(k + 1))
          rcell(k) % branch(2) = known_branch(rcell(k + 1))
       end if
    end do

    g % branch(1) = known_branch(scell(1))
    g % branch(2) = known_branch(rcell(1))

    view = ordinary_graph_view(g, bnd, sets, tail_at=1, head_at=2)

    call report(h % num_tuples() .eq. 2 .and. t % num_tuples() .eq. 3, &
         & "the wall edge stands in T and is an absence in H", nfail)

    ok = .true.
    do v = 1, 3
       ok = ok .and. .not. h % has([3, v])
    end do
    call report(ok, &
         & "no ghost tuple stands across the wall", nfail)

    call report(.not. view % edge_has_head(3) .and. &
         &      view % edge_head(3) .eq. 0, &
         & "and the profile reads the absence as the old zero", nfail)

  end subroutine check_wall_is_absence

  !===================================================================!
  ! Canonical order is the carrier's DECLARED enumeration, not the
  ! members' numeric order. Edges { 30 10 20 } stand in that
  ! declared order in every derived list, however the tuples were
  ! shuffled - a numeric sort would answer 10 20 30 and be wrong.
  !===================================================================!

  subroutine check_declared_order(nfail)

    integer, intent(inout) :: nfail

    type(graph)                     :: edges
    type(set_map)                   :: sets
    type(graph)               :: verts
    type(csr_relation)              :: t, h
    type(graph)            , target :: g
    type(graph)            , target :: scell(2), selem(2), rcell(2), relem(2)
    type(relational_binding)        :: bnd
    type(ordinary_graph_view)       :: view
    integer, allocatable            :: idx(:)
    integer                         :: k

    ! A listed domain: sparse, unordered-world indices. It was a
    ! second member_set concretion once, to prove the relation generic
    ! over carriers; the generality now lives in the representation,
    ! so the fixture is gone and this is one bind.
    call edges % declare()
    call sets  % bind(edges, listed_set_representation([30, 10, 20]))
    call verts % declare()
    call sets % bind(verts, counted_set_representation(2))

    ! Every edge runs 1 -> 2; tuples handed in scrambled order.
    t = csr_relation('tail', edges, verts, &
         & reshape([10,1,  30,1,  20,1], [2, 3]), sets)
    h = csr_relation('head', edges, verts, &
         & reshape([20,2,  10,2,  30,2], [2, 3]), sets)

    call g % declare()
    do k = 1, 2
       call scell(k) % declare(); call selem(k) % declare()
       call rcell(k) % declare(); call relem(k) % declare()
    end do

    call bnd % bind_set(selem(1), verts)
    call bnd % bind_set(selem(2), edges)
    call bnd % bind_relation(relem(1), t)
    call bnd % bind_relation(relem(2), h)

    do k = 1, 2
       scell(k) % branch(1) = known_branch(selem(k))
       rcell(k) % branch(1) = known_branch(relem(k))
       if (k .lt. 2) then
          scell(k) % branch(2) = known_branch(scell(k + 1))
          rcell(k) % branch(2) = known_branch(rcell(k + 1))
       end if
    end do

    g % branch(1) = known_branch(scell(1))
    g % branch(2) = known_branch(rcell(1))

    view = ordinary_graph_view(g, bnd, sets, tail_at=1, head_at=2)

    call view % outgoing_edges(1, idx)
    call report(all(idx .eq. [30, 10, 20]), &
         & "outgoing edges stand in declared order, not numeric", nfail)

    call view % incoming_edges(2, idx)
    call report(all(idx .eq. [30, 10, 20]), &
         & "and so do incoming edges", nfail)

    call view % incident_edges(1, idx)
    call report(all(idx .eq. [30, 10, 20]), &
         & "and the merged incidence follows the same declaration", nfail)

    call report(view % edge_tail(20) .eq. 1 .and. &
         &      view % edge_head(30) .eq. 2, &
         & "the ends answer by member value, wherever it stands", nfail)

  end subroutine check_declared_order

end program test_graph_ordinary
