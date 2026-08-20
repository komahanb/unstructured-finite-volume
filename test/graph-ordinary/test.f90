!=====================================================================!
! The directed-reading equivalence suite. The directed graph is one
! schema over relations: T <= E x V total, H <= E x V partial, the
! boundary an absence in H. On every topology below the stored
! graph is held against externally built T and H relations, edge
! for edge and vertex for vertex, with per-vertex lists compared as
! sets (a relation does not remember tuple order). The graph's own
! derived relations (tail_relation, head_relation) are checked
! against its stored answers on a 400-vertex pseudo-random
! topology.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_graph_ordinary

  use view_directed_stored          , only : directed_stored_graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set           , only : set_map
  use relation_binary, only : csr_relation
  use graph_fractal        , only : graph

  implicit none

  integer :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "directed reading equivalence suite"
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

  call compare_large_topology(nfail)
  call check_relations_of_the_graph(nfail)

  call check_wall_is_absence(nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all directed reading checks passed"
  else
     write(*,'(1x,a,i0,a)') "FAILED: ", nfail, " directed reading check(s)"
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
  ! Sort a small integer list ascending, in place: the set-equality
  ! comparisons above need one canonical order on both sides.
  !===================================================================!

  pure subroutine ascending(a)

    integer, intent(inout) :: a(:)

    integer :: i, j, key

    do i = 2, size(a)
       key = a(i)
       j = i - 1
       do while (j >= 1)
          if (a(j) <= key) exit
          a(j + 1) = a(j)
          j = j - 1
       end do
       a(j + 1) = key
    end do

  end subroutine ascending

  !===================================================================!
  ! The oracle bout: build the old directed_stored_graph and the relational
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

    type(directed_stored_graph)     :: old
    type(graph)               :: verts, edges
    type(set_map)               :: sets
    type(csr_relation)              :: t, h
    integer, allocatable            :: ttab(:,:), htab(:,:)
    integer, allocatable            :: a(:), b(:)
    integer                         :: ne, nh, e, v
    logical                         :: ok

    ne = size(tails)

    old = directed_stored_graph(nv, tails=tails, heads=heads)

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

    ! The stored graph against the relations' fibres. A relation is
    ! a set, so the per-vertex lists are compared as sets: sorted
    ! copies, equal membership. The per-edge fibres have at most one
    ! member and need no order.
    ok = .true.

    do e = 1, ne
       call t % image(e, a)
       ok = ok .and. size(a) == 1
       if (ok) ok = a(1) == old % edge_tail(e)
       call h % image(e, a)
       if (old % edge_has_head(e)) then
          ok = ok .and. size(a) == 1
          if (ok) ok = a(1) == old % edge_head(e)
       else
          ok = ok .and. size(a) == 0
       end if
    end do

    do v = 1, nv
       call t % preimage(v, a)
       call old % outgoing_edges(v, b)
       call ascending(a); call ascending(b)
       ok = ok .and. size(a) == size(b)
       if (ok .and. size(a) > 0) ok = all(a == b)
       call h % preimage(v, a)
       call old % incoming_edges(v, b)
       call ascending(a); call ascending(b)
       ok = ok .and. size(a) == size(b)
       if (ok .and. size(a) > 0) ok = all(a == b)
    end do

    call report(ok, "the relations match the stored graph on " // what, nfail)

  end subroutine compare_topology

  !===================================================================!
  ! The boundary law, said directly: the wall edge stands in T and
  ! is an absence in H - no imaginary member, no ghost tuple.
  !===================================================================!

  !===================================================================!
  ! The equivalence law at scale (consolidation plan, phase 1): a
  ! 400-vertex pseudo-random topology - three edges per vertex to
  ! deterministic LCG targets, roughly one in seven headless - is
  ! compared question by question between the two constructions,
  ! tuples shuffled. The remaining contract questions (the set,
  ! tag, and ownership battery) exist only on the stored graph and
  ! are characterized by the wider battery.
  !===================================================================!

  subroutine compare_large_topology(nfail)

    integer, intent(inout) :: nfail

    integer :: tails(1200), heads(1200)

    call pseudo_random_topology(400, tails, heads)

    call compare_topology(400, tails, heads, &
         & "a 400-vertex pseudo-random topology, tuples shuffled", &
         & nfail, scrambled=.true.)

  end subroutine compare_large_topology

  !===================================================================!
  ! Three edges per vertex to deterministic LCG targets, roughly
  ! one in seven headless. Seedless, so every run builds the same
  ! topology.
  !===================================================================!

  subroutine pseudo_random_topology(nv, tails, heads)

    integer, intent(in)  :: nv
    integer, intent(out) :: tails(:), heads(:)

    integer :: v, k, e, state

    state = 12345
    e = 0
    do v = 1, nv
       do k = 1, 3
          e = e + 1
          tails(e) = v
          state    = mod(1103515245 * state + 12345, 2147483647)
          if (mod(state, 7) == 0) then
             heads(e) = 0                       ! a wall
          else
             heads(e) = 1 + mod(state, nv)      ! any vertex, self allowed
          end if
       end do
    end do

  end subroutine pseudo_random_topology

  !===================================================================!
  ! The graph's own relation reading (consolidation plan, phase 3):
  ! tail_relation and head_relation, derived on request, must agree
  ! with the graph's stored answers - T's fibre at an edge is its
  ! tail, H's fibre its head or empty at a wall, and the preimages
  ! at a vertex are exactly the outgoing and incoming edge lists.
  !===================================================================!

  subroutine check_relations_of_the_graph(nfail)

    integer, intent(inout) :: nfail

    integer :: tails(1200), heads(1200)
    type(directed_stored_graph) :: g
    type(csr_relation) :: t, h
    integer, allocatable :: idx(:), ref(:)
    integer :: e, v
    logical :: ok

    call pseudo_random_topology(400, tails, heads)
    g = directed_stored_graph(400, tails=tails, heads=heads)

    t = g % tail_relation()
    h = g % head_relation()

    ok = .true.

    do e = 1, g % num_edges()
       call t % image(e, idx)
       ok = ok .and. size(idx) == 1
       if (ok) ok = idx(1) == g % edge_tail(e)
       call h % image(e, idx)
       if (g % edge_has_head(e)) then
          ok = ok .and. size(idx) == 1
          if (ok) ok = idx(1) == g % edge_head(e)
       else
          ok = ok .and. size(idx) == 0
       end if
    end do

    call report(ok, &
         & "the derived T and H fibre every edge to its own ends", nfail)

    ok = .true.
    do v = 1, g % num_vertices()
       call t % preimage(v, idx)
       call g % outgoing_edges(v, ref)
       ok = ok .and. size(idx) == size(ref)
       if (ok .and. size(idx) > 0) ok = all(idx == ref)
       call h % preimage(v, idx)
       call g % incoming_edges(v, ref)
       ok = ok .and. size(idx) == size(ref)
       if (ok .and. size(idx) > 0) ok = all(idx == ref)
    end do

    call report(ok, &
         & "the derived preimages are the outgoing and incoming lists", nfail)

  end subroutine check_relations_of_the_graph

  subroutine check_wall_is_absence(nfail)

    integer, intent(inout) :: nfail

    type(directed_stored_graph)     :: old
    type(graph)               :: verts, edges
    type(set_map)                   :: sets
    type(csr_relation)              :: t, h
    integer                         :: v
    logical                         :: ok

    call verts % declare()
    call sets % bind(verts, counted_set_representation(3))
    call edges % declare()
    call sets % bind(edges, counted_set_representation(3))

    t = csr_relation('tail', edges, verts, &
         & reshape([1,1,  2,2,  3,3], [2, 3]), sets)
    h = csr_relation('head', edges, verts, &
         & reshape([1,2,  2,3], [2, 2]), sets)

    call report(h % num_tuples() .eq. 2 .and. t % num_tuples() .eq. 3, &
         & "the wall edge stands in T and is an absence in H", nfail)

    ok = .true.
    do v = 1, 3
       ok = ok .and. .not. h % has([3, v])
    end do
    call report(ok, &
         & "no ghost tuple stands across the wall", nfail)

    old = directed_stored_graph(3, tails=[1, 2, 3], heads=[2, 3, 0])
    call report(.not. old % edge_has_head(3) .and. &
         &      old % edge_head(3) .eq. 0, &
         & "and the stored graph reads the absence as zero", nfail)

  end subroutine check_wall_is_absence



end program test_graph_ordinary
