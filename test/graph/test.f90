!=====================================================================!
! Standalone graph suite: exercises the abstract graph by itself - no
! mesh, no solver. Checks, per the plan:
!   1. adjacency queries against a small graph with known answers
!   2. forward and reverse traversal orders (dag/chain + general graph)
!   3. the chain tenant's rule-generated adjacency (nothing materialized)
!   4. partition invariants: every vertex owned exactly once, ghosts
!      consistent, edge cut as reported
!   5. the retained adjacency matches the construction stamp_bfs used to
!      scratch-build (hand-computed xadj/adj on the known graph)
!=====================================================================!

module class_test_graph

  ! a minimal stored tenant: vertices and an edge list set directly,
  ! exactly how a tenant consumes the ancestor

  use interface_graph, only : graph, vertex, edge

  implicit none

  private
  public :: test_graph

  type, extends(graph) :: test_graph
  end type test_graph

  interface test_graph
     module procedure create
  end interface test_graph

contains

  pure type(test_graph) function create(nv, tails, heads) result(this)

    integer, intent(in) :: nv
    integer, intent(in) :: tails(:), heads(:)

    integer :: i

    this % num_vertices = nv
    this % num_edges    = size(tails)

    allocate(this % vertices(nv))
    do i = 1, nv
       this % vertices(i) % number = i
       this % vertices(i) % part   = 1
    end do

    allocate(this % edges(this % num_edges))
    do i = 1, this % num_edges
       this % edges(i) % tail = tails(i)
       this % edges(i) % head = heads(i)
    end do

    call this % build_adjacency()

  end function create

end module class_test_graph

program test_graph_suite

  use iso_fortran_env  , only : dp => REAL64
  use interface_graph  , only : graph
  use class_test_graph , only : test_graph
  use class_chain      , only : chain
  use module_solve_mode, only : FORWARD, REVERSE

  implicit none

  integer :: nfail

  nfail = 0

  call check_known_adjacency(nfail)
  call check_traversal_orders(nfail)
  call check_chain_rule(nfail)
  call check_partition_invariants(nfail)
  call check_dof_map(nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all graph checks passed"
  else
     write(*,'(1x,a,i0,a)') "FAILED: ", nfail, " graph check(s)"
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
  ! The known graph (the chapter's formal-definition example, 1-based):
  ! V = {1..5}, E = {(1,2),(1,3),(1,4),(2,4),(3,4),(4,5)}
  !===================================================================!

  type(test_graph) function known_graph() result(g)
    g = test_graph(5, tails=[1,1,1,2,3,4], heads=[2,3,4,4,4,5])
  end function known_graph

  !===================================================================!
  ! 1 + 5: neighbours/degree with known answers, and the retained
  ! adjacency (xadj/adj) exactly as the scratch construction built it.
  !===================================================================!

  subroutine check_known_adjacency(nfail)

    integer, intent(inout) :: nfail
    type(test_graph) :: g

    g = known_graph()

    call report(all(g % neighbours(1) .eq. [2,3,4]),   "neighbours(1) = [2,3,4]",   nfail)
    call report(all(g % neighbours(4) .eq. [1,2,3,5]), "neighbours(4) = [1,2,3,5]", nfail)
    call report(all(g % neighbours(5) .eq. [4]),       "neighbours(5) = [4]",       nfail)
    call report(g % degree(1) .eq. 3 .and. g % degree(4) .eq. 4 .and. &
         &      g % degree(5) .eq. 1,                  "degrees 3/4/1",             nfail)

    ! retained adjacency, hand-computed: degrees [3,2,2,4,1] ->
    ! xadj = [1,4,6,8,12,13]; adj in edge insertion order
    call report(all(g % xadj .eq. [1,4,6,8,12,13]),    "retained xadj",             nfail)
    call report(all(g % adj  .eq. [2,3,4, 1,4, 1,4, 1,2,3,5, 4]), "retained adj",   nfail)

  end subroutine check_known_adjacency

  !===================================================================!
  ! 2: traversal orders. bfs on the known graph from vertex 1; the
  ! chain's forward order is 1..n and reverse is n..1.
  !===================================================================!

  subroutine check_traversal_orders(nfail)

    integer, intent(inout) :: nfail
    type(test_graph) :: g
    type(chain)      :: c
    integer :: i

    g = known_graph()
    call report(all(g % traversal_order(FORWARD) .eq. [1,2,3,4,5]), &
         & "bfs order on known graph", nfail)
    call report(all(g % traversal_order(REVERSE) .eq. [5,4,3,2,1]), &
         & "reverse order flips it", nfail)

    c = chain(6)
    call report(all(c % traversal_order(FORWARD) .eq. [(i, i=1,6)]), &
         & "chain forward order 1..n", nfail)
    call report(all(c % traversal_order(REVERSE) .eq. [(7-i, i=1,6)]), &
         & "chain reverse order n..1", nfail)

  end subroutine check_traversal_orders

  !===================================================================!
  ! 3: the chain answers by rule and materializes nothing.
  !===================================================================!

  subroutine check_chain_rule(nfail)

    integer, intent(inout) :: nfail
    type(chain) :: c

    c = chain(6)

    call report(all(c % neighbours(1) .eq. [2]),   "chain neighbours(1) = [2]",   nfail)
    call report(all(c % neighbours(3) .eq. [2,4]), "chain neighbours(3) = [2,4]", nfail)
    call report(all(c % neighbours(6) .eq. [5]),   "chain neighbours(6) = [5]",   nfail)
    call report(c % degree(1) .eq. 1 .and. c % degree(3) .eq. 2, &
         & "chain degrees by rule", nfail)
    call report(.not. allocated(c % xadj) .and. .not. allocated(c % edges), &
         & "chain materializes no adjacency and no edge list", nfail)
    call report(c % num_edges .eq. 5, "chain edge count n-1", nfail)

  end subroutine check_chain_rule

  !===================================================================!
  ! 4: partition invariants, on the stored tenant and on the rule
  ! tenant (the inherited partitioner runs on the rule directly).
  !===================================================================!

  subroutine check_partition_invariants(nfail)

    integer, intent(inout) :: nfail
    type(test_graph) :: g
    type(chain)      :: c

    g = known_graph()
    call g % partition(2)
    call assert_partition(g, 2, "known graph", nfail)

    c = chain(10)
    call c % partition(2)
    call assert_partition(c, 2, "chain", nfail)

    ! the chain's 2-part cut is exactly one edge, ghosts one each side
    call report(c % edge_cut() .eq. 1, "chain 2-part cut = 1", nfail)
    call report(c % n_ghosts(1) .eq. 1 .and. c % n_ghosts(2) .eq. 1, &
         & "chain ghosts one per side", nfail)

  end subroutine check_partition_invariants

  subroutine assert_partition(g, nparts, label, nfail)

    class(graph)    , intent(in)    :: g
    integer         , intent(in)    :: nparts
    character(len=*), intent(in)    :: label
    integer         , intent(inout) :: nfail

    integer, allocatable :: own(:), gh(:), nbrs(:), cover(:)
    integer :: k, i, j, v, w, cut
    logical :: ok

    ! every vertex owned exactly once
    allocate(cover(g % num_vertices)); cover = 0
    do k = 1, nparts
       own = g % owned(k)
       do i = 1, size(own)
          cover(own(i)) = cover(own(i)) + 1
       end do
    end do
    call report(all(cover .eq. 1), label//": every vertex owned exactly once", nfail)

    ! each ghost of k is not owned by k and touches an owned vertex of k
    ok = .true.
    do k = 1, nparts
       gh = g % ghosts(k)
       do i = 1, size(gh)
          v = gh(i)
          if (g % part_of(v) .eq. k) ok = .false.
          nbrs = g % neighbours(v)
          if (.not. any([(g % part_of(nbrs(j)) .eq. k, j=1,size(nbrs))])) ok = .false.
       end do
    end do
    call report(ok, label//": ghosts consistent", nfail)

    ! edge cut as reported (independent recount over neighbours)
    cut = 0
    do v = 1, g % num_vertices
       nbrs = g % neighbours(v)
       do j = 1, size(nbrs)
          w = nbrs(j)
          if (w .gt. v .and. g % part_of(v) .ne. g % part_of(w)) cut = cut + 1
       end do
    end do
    call report(cut .eq. g % ncut, label//": edge cut as reported", nfail)

  end subroutine assert_partition

  !===================================================================!
  ! dof interleaving on the ancestor (variable-fastest)
  !===================================================================!

  subroutine check_dof_map(nfail)

    integer, intent(inout) :: nfail
    type(chain) :: c

    c = chain(4, num_variables=2)

    call report(c % dof(1,1) .eq. 1 .and. c % dof(3,2) .eq. 6 .and. &
         &      c % num_dofs() .eq. 8, "dof map variable-fastest", nfail)

  end subroutine check_dof_map

end program test_graph_suite
