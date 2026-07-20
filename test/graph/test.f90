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
!   8. the squint and the zoom: the partition read back as the quotient
!      graph, and the aggregation partitioner discovering its own
!      parts (the chain of 10 huddles into 4, deterministically)
!   9. escape times resolved in one pass agree with orbit-by-orbit
!      painting, capped and uncapped
!  10. the local frame: a part's dofs in local order (owned first,
!      then the halo), the frame read backwards, and a matrix's rows
!      re-expressed in it - the stage-1 gate of dereplication
!=====================================================================!

program test_graph_suite

  use iso_fortran_env  , only : dp => REAL64
  use interface_graph  , only : graph
  use class_stored_graph, only : stored_graph, stored_digraph
  use class_chain      , only : chain
  use class_csr        , only : csr_matrix
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
  call check_quotient(nfail)
  call check_escape_times(nfail)
  call check_local_frame(nfail)

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

  type(stored_graph) function known_graph() result(g)
    g = stored_graph(5, tails=[1,1,1,2,3,4], heads=[2,3,4,4,4,5])
  end function known_graph

  !===================================================================!
  ! 1 + 5: neighbours/degree with known answers, and the retained
  ! adjacency (xadj/adj) exactly as the scratch construction built it.
  !===================================================================!

  subroutine check_known_adjacency(nfail)

    integer, intent(inout) :: nfail
    type(stored_graph) :: g
    integer, allocatable :: xadj(:), adj(:)

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

    ! the same adjacency at dof granularity (one variable, so dof = v):
    ! every row leads with its own dof, then its neighbours in order
    call g % dof_adjacency(xadj, adj)
    call report(all(xadj .eq. [1,5,8,11,16,18]), "dof xadj: rows are 1+degree wide", nfail)
    call report(all(adj  .eq. [1,2,3,4, 2,1,4, 3,1,4, 4,1,2,3,5, 5,4]), &
         & "dof adj: self-loop first, neighbours after", nfail)

    ! greedy coloring, hand-computed: no edge inside a color
    call report(all(g % coloring() .eq. [1,2,2,3,1]), &
         & "greedy coloring leaves no edge inside a color", nfail)

  end subroutine check_known_adjacency

  !===================================================================!
  ! 2: traversal orders. bfs on the known graph from vertex 1; the
  ! chain's forward order is 1..n and reverse is n..1.
  !===================================================================!

  subroutine check_traversal_orders(nfail)

    integer, intent(inout) :: nfail
    type(stored_graph) :: g
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

    ! the chain raised to a power: edge m -> k whenever k - m <= 2,
    ! still all by rule (this is the shape of a time stencil's dag)
    c = chain(6, power = 2)
    call report(all(c % out_neighbours(3) .eq. [4,5]) .and. &
         &      all(c % out_neighbours(5) .eq. [6]), &
         & "powered chain reaches forward to power", nfail)
    call report(all(c % in_neighbours(3) .eq. [1,2]) .and. &
         &      all(c % in_neighbours(2) .eq. [1]), &
         & "powered chain reaches back to power", nfail)
    call report(all(c % neighbours(3) .eq. [1,2,4,5]) .and. &
         &      c % degree(1) .eq. 2 .and. c % degree(3) .eq. 4, &
         & "powered chain unions both reaches", nfail)
    call report(c % num_edges .eq. 9, "powered chain edge count", nfail)
    call report(all(c % dependency_order() .eq. [1,2,3,4,5,6]), &
         & "powered chain dependency order stays 1..n", nfail)
    call report(all(c % coloring() .eq. [1,2,3,1,2,3]), &
         & "a chain at power p colors with p+1 colors", nfail)

  end subroutine check_chain_rule

  !===================================================================!
  ! 4: partition invariants, on the stored subclass and on the rule
  ! subclass (the inherited partitioner runs on the rule directly).
  !===================================================================!

  subroutine check_partition_invariants(nfail)

    integer, intent(inout) :: nfail
    type(stored_graph) :: g
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

    type(stored_digraph) :: diamond, cycle3, schedule
    type(chain)        :: c
    integer, allocatable :: order(:)
    real(dp) :: certification

    ! the diamond dag: 1->2, 1->3, 2->4, 3->4
    diamond = stored_digraph(4, tails=[1,1,2,3], heads=[2,3,4,4])

    call report(diamond % is_acyclic(), "diamond dag is acyclic", nfail)

    order = diamond % dependency_order()
    call report(all(order .eq. [1,2,3,4]), &
         & "dependency order on the diamond", nfail)

    call report(all(diamond % out_neighbours(1) .eq. [2,3]) .and. &
         &      all(diamond % in_neighbours(4)  .eq. [2,3]), &
         & "directed queries on the diamond", nfail)

    ! a 3-cycle must be detected without dying
    cycle3 = stored_digraph(3, tails=[1,2,3], heads=[2,3,1])
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

    ! the directed citizen carries caller-given numbers, so a small
    ! graph can name things beyond itself. The fixture avoids every
    ! coincidence: no number equals its own vertex index and the list
    ! reads differently reversed, so ignoring or reversing the given
    ! numbers both fail. The neighbour comparisons check their sizes
    ! first - all() on an accidentally empty return says yes to
    ! anything.
    schedule = stored_digraph(5, tails=[1,2,3,4], heads=[2,3,4,5], &
         &                    numbers=[9,8,7,6,4])
    call report(schedule % vertices(2) % number .eq. 8 .and. &
         &      schedule % vertices(5) % number .eq. 4 .and. &
         &      size(schedule % out_neighbours(2)) .eq. 1 .and. &
         &      all(schedule % out_neighbours(2) .eq. [3]) .and. &
         &      size(schedule % out_neighbours(5)) .eq. 0, &
         & "stored digraph carries caller-given numbers", nfail)

    ! and the graph reads its own route: follow the arrows from the
    ! source, hand back the numbers in visit order (size checked
    ! first - all() on an accidentally empty return says yes)
    call report(size(schedule % source_path()) .eq. 5 .and. &
         &      all(schedule % source_path() .eq. [9,8,7,6,4]), &
         & "source path returns the numbers in trip order", nfail)

  end subroutine check_directed_structure

  !===================================================================!
  ! 7: orbits under a successor rule - one vertex, one arrow out,
  ! repeated. Escape (the rule leaves the vertex set; the length is
  ! the escape time), cycle closure (the repeated vertex is kept as
  ! the final entry), the tail-into-cycle shape, and the step limit.
  !===================================================================!

  subroutine check_orbit(nfail)

    integer, intent(inout) :: nfail
    type(stored_graph) :: g
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
  ! 8: the squint. Partition the known graph in two and read it back
  ! as the quotient - two coarse vertices, one coarse edge, however
  ! many fine edges cross. Then the aggregation partitioner on the
  ! chain of 10: no part count given, and Vanek's passes discover
  ! exactly four huddles, deterministically.
  !===================================================================!

  subroutine check_quotient(nfail)

    integer, intent(inout) :: nfail
    type(stored_graph) :: g, coarse, refined
    type(chain)        :: c
    integer  :: v
    real(dp) :: x(12)

    ! bfs partition puts {1,2,3} and {4,5} apart; three fine edges
    ! cross, and the squint sees one coarse edge
    g = known_graph()
    call g % partition(2)
    coarse = stored_graph(g)

    call report(coarse % num_vertices .eq. 2 .and. coarse % num_edges .eq. 1, &
         & "quotient: two parts squint to two vertices, one edge", nfail)
    call report(all(coarse % neighbours(1) .eq. [2]) .and. &
         &      all(coarse % neighbours(2) .eq. [1]), &
         & "quotient adjacency reads back", nfail)

    ! the chain of 10 huddles into 4 parts: {1,2} {3,4,5} {6,7,8} {9,10}
    c = chain(10)
    call c % partition_aggregate()
    call report(c % nparts .eq. 4, "aggregation discovers 4 parts on chain(10)", nfail)
    call report(all([(c % part_of(v), v = 1, 10)] .eq. [1,1,2,2,2,3,3,3,4,4]), &
         & "aggregation huddles the chain as expected", nfail)

    ! the quotient of the aggregated chain is again a chain
    coarse = stored_graph(c)
    call report(coarse % num_vertices .eq. 4 .and. coarse % num_edges .eq. 3 .and. &
         &      all(coarse % neighbours(2) .eq. [1,3]), &
         & "the squinted chain is again a chain", nfail)

    ! the zoom: every vertex of the known graph splits into three
    ! children - chains inside each split vertex, one bridge per
    ! original edge - and the parent map arrives as the partition
    g       = known_graph()
    refined = stored_graph(g, 3)
    call report(refined % num_vertices .eq. 15 .and. &
         &      refined % num_edges .eq. 16, &
         & "refinement: 5 vertices split into 15, 10 chain + 6 bridge edges", nfail)
    call report(refined % part_of(4) .eq. 2 .and. refined % part_of(15) .eq. 5, &
         & "refinement carries the parent map as its partition", nfail)

    ! the squint undoes the zoom: the quotient of the refinement is
    ! the original graph
    coarse = stored_graph(refined)
    call report(coarse % num_vertices .eq. 5 .and. coarse % num_edges .eq. 6 .and. &
         &      size(coarse % neighbours(4)) .eq. 4 .and. &
         &      all(coarse % neighbours(4) .eq. [1,2,3,5]), &
         & "the squint undoes the zoom", nfail)

    ! a rule graph refines as readily as a stored one, and an adopted
    ! partition gathers like any other
    c = chain(6, 2)
    call c % set_partition([1,1,2,2,3,3])
    call report(c % nparts .eq. 3 .and. size(c % owned(2)) .eq. 2 .and. &
         &      all(c % owned(2) .eq. [3,4]), &
         & "an adopted partition gathers like any other", nfail)

    ! values ride the partition: part 2 owns vertices 3 and 4, so with
    ! two variables its dofs are 5..8; gather pulls those entries out,
    ! scatter pushes them back and touches nothing else
    call report(all(c % dofs_of(c % owned(2)) .eq. [5,6,7,8]), &
         & "a vertex list expands to its dofs, variable-fastest", nfail)
    x = [(real(v, dp), v = 1, 12)]
    call report(all(c % gather(2, x) .eq. [5,6,7,8]), &
         & "gather pulls the values at a part's owned dofs", nfail)
    call c % scatter(2, [50.0_dp, 60.0_dp, 70.0_dp, 80.0_dp], x)
    call report(all(x .eq. [1,2,3,4,50,60,70,80,9,10,11,12]), &
         & "scatter pushes them back and leaves the rest alone", nfail)

    ! dot measures where gather reaches: part 2's dofs are 5..8, so
    ! its dot of x=(1..12) with itself is 5^2+..+8^2 = 174, and the
    ! parts' dots sum to the whole graph's dot (each dof owned once)
    x = [(real(v, dp), v = 1, 12)]
    call report(abs(c % dot(2, x, x) - 174.0_dp) .lt. 1.0e-12_dp, &
         & "dot measures a part's owned dofs", nfail)
    call report(abs(c % dot(1, x, x) + c % dot(2, x, x) + c % dot(3, x, x) &
         &          - dot_product(x, x)) .lt. 1.0e-12_dp, &
         & "the parts' dots sum to the whole graph's dot", nfail)
    refined = stored_graph(c, 2)
    call report(refined % num_vertices .eq. 12 .and. refined % num_edges .eq. 11, &
         & "the rule-generated chain refines by its neighbour queries", nfail)

  end subroutine check_quotient

  !===================================================================!
  ! 9: escape times in one pass. The same rules the orbit checks use,
  ! resolved for every vertex at once, must agree with what orbit-by-
  ! orbit painting gives: [5,4,3,2,1] for the escaping rule, all home
  ! at the limit for the tail-into-cycle rule, and the cap biting.
  !===================================================================!

  subroutine check_escape_times(nfail)

    integer, intent(inout) :: nfail
    type(stored_graph) :: g
    integer, allocatable :: times(:)

    g = known_graph()

    times = g % escape_times(next_vertex, 40)
    call report(all(times .eq. [5,4,3,2,1]), &
         & "escape times count down the escaping rule", nfail)

    times = g % escape_times(next_vertex, 3)
    call report(all(times .eq. [3,3,3,2,1]), &
         & "the cap bites the long escapes", nfail)

    times = g % escape_times(tail_then_cycle, 40)
    call report(all(times .eq. 40), &
         & "cycle and its tail are all home at the limit", nfail)

    times = g % escape_times(same_vertex, 7)
    call report(all(times .eq. 7), &
         & "fixed vertices are home", nfail)

    ! the cap may be huge(0) - "no cap" - without the unwind wrapping
    ! negative: escapers keep their exact times, home is painted the cap
    times = g % escape_times(next_vertex, huge(0))
    call report(all(times .eq. [5,4,3,2,1]), &
         & "an uncapped limit leaves escape times exact", nfail)

    times = g % escape_times(tail_then_cycle, huge(0))
    call report(all(times .eq. huge(0)), &
         & "an uncapped limit paints home without overflowing", nfail)

  end subroutine check_escape_times

  !===================================================================!
  ! 10: the local frame. On the chain with two variables and parts
  ! [1,1,2,2,3,3], part 2 owns vertices 3,4 (dofs 5..8) and its halo
  ! is vertices 2,5 (dofs 3,4,9,10):
  !
  !    frame(2)  =  [ 5 6 7 8 | 3 4 9 10 ]
  !                   owned     ghosts
  !
  ! and a hand matrix re-expressed in a hand frame computes exactly
  ! its owned rows through the local block.
  !===================================================================!

  subroutine check_local_frame(nfail)

    integer, intent(inout) :: nfail
    type(chain)      :: c
    type(csr_matrix) :: A, B
    real(dp)         :: y(2)

    c = chain(6, 2)
    call c % set_partition([1,1,2,2,3,3])

    call report(all(c % frame(2) .eq. [5,6,7,8, 3,4,9,10]), &
         & "the frame lists owned dofs first, then the halo", nfail)
    call report(all(c % frame_inverse(2) .eq. [0,0,5,6,1,2,3,4,7,8,0,0]), &
         & "the frame reads backwards, 0 where the part is blind", nfail)

    ! a 4x4 hand matrix; a part owning rows 2,3 reaches columns
    ! 1..4, so its hand frame is [2,3 | 1,4] and loc = [3,1,2,4]
    A = csr_matrix(4, 4, [1,3,6,9,11], [1,2, 1,2,3, 2,3,4, 3,4], &
         & [10.0_dp,1.0_dp, 2.0_dp,20.0_dp,3.0_dp, 4.0_dp,30.0_dp,5.0_dp, 6.0_dp,40.0_dp])
    B = A % local_block([2,3], [3,1,2,4], 4)

    call report(B % num_vertices .eq. 2 .and. B % ncols .eq. 4 .and. &
         &      all(B % out_xadj .eq. [1,4,7]) .and. &
         &      all(B % out_adj  .eq. [3,1,2, 1,2,4]) .and. &
         &      all(B % vals     .eq. [2.0_dp,20.0_dp,3.0_dp, 4.0_dp,30.0_dp,5.0_dp]), &
         & "the local block keeps rows whole, far ends in the frame", nfail)

    ! the block's matvec on a frame-ordered vector IS the owned rows
    ! of the global product: x = [1,2,3,4] framed as [2,3,1,4]
    call B % matvec([2.0_dp,3.0_dp,1.0_dp,4.0_dp], y)
    call report(all(abs(y - [51.0_dp,118.0_dp]) .lt. 1.0e-14_dp), &
         & "the local matvec computes exactly the owned rows", nfail)

  end subroutine check_local_frame

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
