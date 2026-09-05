!=====================================================================!
! The directed-reading equivalence suite. The directed graph is one
! schema over relations: T <= E x V total, H <= E x V partial, the
! boundary an absence in H. On every topology below the stored
! graph is compared with externally built T and H relations, edge
! for edge and vertex for vertex, with per-vertex lists compared as
! sets (a relation does not remember tuple order), up to a
! 400-vertex pseudo-random topology.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_graph_ordinary

  use view_directed_stored          , only : stored_directed_graph
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

  ! A diamond with a boundary edge and an isolated vertex: 1->2->4,
  ! 1->3->4, 4->boundary, vertex 5 alone.
  call compare_topology(5, [1, 1, 2, 3, 4], [2, 3, 4, 4, 0], &
       & "diamond with a boundary edge and an isolated vertex", nfail)

  ! Two edges joining one pair.
  call compare_topology(2, [1, 1], [2, 2], &
       & "a parallel pair", nfail)

  ! The chain with a boundary edge at its end.
  call compare_topology(3, [1, 2, 3], [2, 3, 0], &
       & "a chain into a boundary edge", nfail)

  ! A self-loop beside ordinary edges: 1->2, 2->2, 2->3, 3->boundary.
  ! incident counts the loop twice, adjacency never names the vertex
  ! to itself.
  call compare_topology(3, [1, 2, 2, 3], [2, 2, 3, 0], &
       & "a self-loop beside a boundary edge", nfail)

  ! A relation is a set: the same topologies in reversed tuple order
  ! give the same fibres.
  call compare_topology(5, [1, 1, 2, 3, 4], [2, 3, 4, 4, 0], &
       & "the diamond, tuples reversed", nfail, reversed=.true.)
  call compare_topology(3, [1, 2, 2, 3], [2, 2, 3, 0], &
       & "the self-loop, tuples reversed", nfail, reversed=.true.)

  call compare_large_topology(nfail)

  call check_boundary_is_absence(nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all directed reading checks passed"
  else
     write(*,'(1x,a,i0,a)') "FAILED: ", nfail, " directed reading check(s)"
     error stop
  end if

contains

  subroutine report(passes, label, nfail)

    logical         , intent(in)    :: passes
    character(len=*), intent(in)    :: label
    integer         , intent(inout) :: nfail

    if (passes) then
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
  ! Build the stored_directed_graph and the relations T and H on one
  ! topology, then compare every fibre. T contains every edge; H the
  ! headed ones; both tables in ascending edge order, as every mesh
  ! builder passes them.
  !===================================================================!

  subroutine compare_topology(nv, tails, heads, what, nfail, reversed)

    integer         , intent(in)           :: nv
    integer         , intent(in)           :: tails(:), heads(:)
    character(len=*), intent(in)           :: what
    integer         , intent(inout)        :: nfail
    logical         , intent(in), optional :: reversed

    type(stored_directed_graph)     :: stored
    type(graph)               :: verts, edges
    type(set_map)               :: sets
    type(csr_relation)              :: t, h
    integer, allocatable            :: ttab(:,:), htab(:,:)
    integer, allocatable            :: a(:), b(:)
    integer                         :: ne, nh, e, v
    logical                         :: passes

    ne = size(tails)

    stored = stored_directed_graph(nv, tails=tails, heads=heads)

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

    ! A set does not remember the order its tuples were listed in:
    ! reverse the columns and every fibre is unchanged.
    if (present(reversed)) then
       if (reversed) then
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
    passes = .true.

    do e = 1, ne
       call t % image(e, a)
       passes = passes .and. size(a) == 1
       if (passes) passes = a(1) == stored % edge_tail(e)
       call h % image(e, a)
       if (stored % edge_has_head(e)) then
          passes = passes .and. size(a) == 1
          if (passes) passes = a(1) == stored % edge_head(e)
       else
          passes = passes .and. size(a) == 0
       end if
    end do

    do v = 1, nv
       call t % preimage(v, a)
       call stored % outgoing_edges(v, b)
       call ascending(a); call ascending(b)
       passes = passes .and. size(a) == size(b)
       if (passes .and. size(a) > 0) passes = all(a == b)
       call h % preimage(v, a)
       call stored % incoming_edges(v, b)
       call ascending(a); call ascending(b)
       passes = passes .and. size(a) == size(b)
       if (passes .and. size(a) > 0) passes = all(a == b)
    end do

    call report(passes, "the relations match the stored graph on " // what, nfail)

  end subroutine compare_topology

  !===================================================================!
  ! The equivalence law at scale: a 400-vertex pseudo-random
  ! topology - three edges per vertex to deterministic LCG targets,
  ! roughly one in seven headless - is compared fibre by fibre
  ! between the two constructions, tuples reversed. The remaining
  ! contract (the set, tag, and ownership queries) exists only on
  ! the stored graph and is checked by the wider suites.
  !===================================================================!

  subroutine compare_large_topology(nfail)

    integer, intent(inout) :: nfail

    integer :: tails(1200), heads(1200)

    call pseudo_random_topology(400, tails, heads)

    call compare_topology(400, tails, heads, &
         & "a 400-vertex pseudo-random topology, tuples reversed", &
         & nfail, reversed=.true.)

  end subroutine compare_large_topology

  !===================================================================!
  ! Three edges per vertex to deterministic LCG targets, roughly
  ! one in seven headless. The seed is fixed, so every run builds
  ! the same topology.
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
             heads(e) = 0                       ! a boundary edge
          else
             heads(e) = 1 + mod(state, nv)      ! any vertex, self allowed
          end if
       end do
    end do

  end subroutine pseudo_random_topology

  subroutine check_boundary_is_absence(nfail)

    integer, intent(inout) :: nfail

    type(stored_directed_graph)     :: stored
    type(graph)               :: verts, edges
    type(set_map)                   :: sets
    type(csr_relation)              :: t, h
    integer                         :: v
    logical                         :: passes

    call verts % declare()
    call sets % bind(verts, counted_set_representation(3))
    call edges % declare()
    call sets % bind(edges, counted_set_representation(3))

    t = csr_relation('tail', edges, verts, &
         & reshape([1,1,  2,2,  3,3], [2, 3]), sets)
    h = csr_relation('head', edges, verts, &
         & reshape([1,2,  2,3], [2, 2]), sets)

    call report(h % num_tuples() .eq. 2 .and. t % num_tuples() .eq. 3, &
         & "the boundary edge is in T and is an absence in H", nfail)

    passes = .true.
    do v = 1, 3
       passes = passes .and. .not. h % has([3, v])
    end do
    call report(passes, &
         & "no invented tuple crosses the boundary", nfail)

    stored = stored_directed_graph(3, tails=[1, 2, 3], heads=[2, 3, 0])
    call report(.not. stored % edge_has_head(3) .and. &
         &      stored % edge_head(3) .eq. 0, &
         & "and the stored graph reads the absence as zero", nfail)

  end subroutine check_boundary_is_absence

end program test_graph_ordinary
