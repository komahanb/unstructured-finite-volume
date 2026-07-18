!=====================================================================!
! Standalone graph suite: exercises the abstract graph by itself - no
! mesh, no solver. Checks, per the plan:
!   1. adjacency queries against a small graph with known answers
!   2. forward and reverse traversal orders (dag/chain + general graph)
!   3. the chain subclass's rule-generated adjacency (nothing materialized)
!   4. partition invariants: every vertex owned exactly once, ghosts
!      consistent, edge cut as reported
!   5. the retained adjacency matches the construction stamp_bfs used to
!      scratch-build (hand-computed xadj/adj on the known graph)
!   6. directed structure: dependency order on the diamond dag, cycle
!      detection on a 3-cycle, and the adjoint walk certified through
!      the type-bound witness (chain analytic + diamond nudge)
!   7. orbits under a successor rule: escape time, cycle closure, the
!      tail-into-cycle shape, and the step limit
!=====================================================================!

module class_test_graph

  ! a minimal stored subclass: vertices and an edge list set directly,
  ! exactly how a subclass consumes the ancestor

  use interface_graph, only : graph, digraph, vertex, edge

  implicit none

  private
  public :: test_graph
  public :: test_digraph

  type, extends(graph) :: test_graph
   contains
     ! the deferred contract, delegated to the stored mechanism
     procedure :: neighbours => test_neighbours
     procedure :: degree     => test_degree
  end type test_graph

  ! a minimal stored directed fixture for the structure checks
  type, extends(digraph) :: test_digraph
   contains
     procedure :: out_neighbours => test_out_neighbours
     procedure :: in_neighbours  => test_in_neighbours
  end type test_digraph

  interface test_graph
     module procedure create
  end interface test_graph

  interface test_digraph
     module procedure create_directed
  end interface test_digraph

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

  pure function test_neighbours(this, v) result(nbrs)
    class(test_graph), intent(in) :: this
    integer          , intent(in) :: v
    integer, allocatable :: nbrs(:)
    nbrs = this % stored_neighbours(v)
  end function test_neighbours

  pure integer function test_degree(this, v)
    class(test_graph), intent(in) :: this
    integer          , intent(in) :: v
    test_degree = this % stored_degree(v)
  end function test_degree

  pure type(test_digraph) function create_directed(nv, tails, heads) result(this)

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

    call this % build_directed_adjacency()

  end function create_directed

  pure function test_out_neighbours(this, v) result(nbrs)
    class(test_digraph), intent(in) :: this
    integer            , intent(in) :: v
    integer, allocatable :: nbrs(:)
    nbrs = this % stored_out_neighbours(v)
  end function test_out_neighbours

  pure function test_in_neighbours(this, v) result(nbrs)
    class(test_digraph), intent(in) :: this
    integer            , intent(in) :: v
    integer, allocatable :: nbrs(:)
    nbrs = this % stored_in_neighbours(v)
  end function test_in_neighbours

end module class_test_graph

program test_graph_suite

  use iso_fortran_env  , only : dp => REAL64
  use interface_graph  , only : graph
  use class_test_graph , only : test_graph, test_digraph
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
  call check_directed_structure(nfail)
  call check_orbit(nfail)

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
  ! 4: partition invariants, on the stored subclass and on the rule
  ! subclass (the inherited partitioner runs on the rule directly).
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
  ! Directed structure: the dependency order on the diamond dag, cycle
  ! detection on a 3-cycle (without dying), and the adjoint walk
  ! certified through the type-bound witness - the chain fixture judged
  ! against the analytic derivative and the diamond fixture against a
  ! central-difference nudge, each at its own tolerance inside the
  ! witness. A certification value below one passes both.
  !===================================================================!

  subroutine check_directed_structure(nfail)

    integer, intent(inout) :: nfail

    ! a witness value below this certifies both fixtures at their own
    ! tolerances (the witness normalizes each defect by its tolerance)
    real(dp), parameter :: certification_threshold = 1.0_dp

    type(test_digraph) :: diamond, cycle3
    type(chain)        :: c
    integer, allocatable :: order(:)
    real(dp) :: certification

    ! the diamond dag: 1->2, 1->3, 2->4, 3->4
    diamond = test_digraph(4, tails=[1,1,2,3], heads=[2,3,4,4])

    call report(diamond % is_acyclic(), "diamond dag is acyclic", nfail)

    order = diamond % dependency_order()
    call report(all(order .eq. [1,2,3,4]), &
         & "dependency order on the diamond", nfail)

    call report(all(diamond % out_neighbours(1) .eq. [2,3]) .and. &
         &      all(diamond % in_neighbours(4)  .eq. [2,3]), &
         & "directed queries on the diamond", nfail)

    ! a 3-cycle must be detected without dying
    cycle3 = test_digraph(3, tails=[1,2,3], heads=[2,3,1])
    call report(.not. cycle3 % is_acyclic(), &
         & "cycle detection refuses the 3-cycle", nfail)

    ! the chain's dependency order is 1..n by construction
    c = chain(5)
    order = c % dependency_order()
    call report(all(order .eq. [1,2,3,4,5]), &
         & "chain dependency order is 1..n", nfail)

    ! the walk's witness: chain analytic + diamond nudge, type-bound
    certification = c % verify_adjoint_accumulation()
    call report(certification .lt. certification_threshold, &
         & "adjoint accumulation certified by the witness", nfail)

  end subroutine check_directed_structure

  !===================================================================!
  ! 7: orbits under a successor rule - one vertex, one arrow out,
  ! repeated. Escape (the rule leaves the vertex set; the length is
  ! the escape time), cycle closure (the repeated vertex is kept as
  ! the final entry), the tail-into-cycle shape, and the step limit.
  !===================================================================!

  subroutine check_orbit(nfail)

    integer, intent(inout) :: nfail
    type(test_graph) :: g
    integer, allocatable :: visited(:)

    g = known_graph()

    ! the rule v -> v+1 escapes past vertex 5: escape time 5
    visited = g % orbit(1, next_vertex)
    call report(all(visited .eq. [1,2,3,4,5]), &
         & "orbit escapes; length is the escape time", nfail)

    ! from an interior start the same rule escapes in two visits
    visited = g % orbit(4, next_vertex)
    call report(all(visited .eq. [4,5]), "orbit from an interior start", nfail)

    ! a vertex that is its own successor closes the shortest cycle
    visited = g % orbit(3, same_vertex)
    call report(all(visited .eq. [3,3]), "fixed vertex closes a 1-cycle", nfail)

    ! tail 1 -> 2 into the cycle 3 -> 4 -> 5 -> 3: the repeated vertex
    ! ends the sequence, so the cycle is the slice between its two
    ! appearances; the length is num_vertices + 1, the default limit
    ! the pigeonhole argument promises is enough
    visited = g % orbit(1, tail_then_cycle)
    call report(all(visited .eq. [1,2,3,4,5,3]), &
         & "tail into a cycle; repeated vertex ends the sequence", nfail)

    ! the step limit caps the orbit before the cycle closes
    visited = g % orbit(1, tail_then_cycle, limit=3)
    call report(all(visited .eq. [1,2,3]), "the step limit caps the orbit", nfail)

  end subroutine check_orbit

  !===================================================================!
  ! successor rules for the orbit checks: one arrow out of every vertex
  !===================================================================!

  pure integer function next_vertex(v)
    integer, intent(in) :: v
    next_vertex = v + 1
  end function next_vertex

  pure integer function same_vertex(v)
    integer, intent(in) :: v
    same_vertex = v
  end function same_vertex

  ! the tail 1 -> 2 feeding the cycle 3 -> 4 -> 5 -> 3
  pure integer function tail_then_cycle(v)
    integer, intent(in) :: v
    if (v .lt. 5) then
       tail_then_cycle = v + 1
    else
       tail_then_cycle = 3
    end if
  end function tail_then_cycle

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
