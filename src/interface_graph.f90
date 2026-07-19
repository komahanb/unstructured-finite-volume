!=====================================================================!
! Abstract graph: vertices, edges, the retained adjacency, traversal
! orders, and the partition contract. This base class holds only
! graph-structured data - never a mesh, a face, or any payload.
! Dependencies point downward only: a subclass knows its own domain
! (class_mesh IS the mesh graph, class_chain knows its recurrence rule);
! the base class knows neither.
!
! The adjacency is a compressed neighbour list (xadj/adj), retained and
! publicly queryable - previously buried as scratch inside the bfs
! partitioner and rebuilt on every use. Every piece of machinery here
! consumes only the neighbour queries, so a subclass may either store the
! adjacency (build_adjacency from an edge list) or generate it by rule
! (the chain overrides neighbours/degree and stores nothing).
!
! Traversal: traversal_order(mode) is a VISIT order (breadth-first,
! FORWARD; the same order reversed, REVERSE) - locality for the
! partitioner, with no dependency claim. Dependency semantics are true
! only on the directed graph: digraph carries out/in neighbour queries,
! a topological dependency_order with a cycle refusal, and the discrete
! adjoint's reverse-mode accumulation, whose direction is read from
! structure rather than recovered from a visit order.
!
! Orbit: orbit(start, successor, limit) follows a caller-supplied rule
! giving every vertex one successor - one vertex, one arrow out,
! repeated. The sequence ends when the rule leaves the vertex set (the
! orbit escapes; its length is the escape time), when a vertex repeats
! (a cycle has closed), or at the step limit.
!
! Partition: stamped in place (partition / partition_rcb /
! partition_aggregate), then queried (part_of, owned, ghosts, balance,
! edge_cut). A partition of a graph is still graph-shaped data.
!
! Squint: quotient_edges reads the stamped partition back as a graph -
! every part one coarse vertex, every pair of touching parts one
! coarse edge. The quotient of a graph is the same kind of animal, so
! a hierarchy is just repeated squinting (class_stored_graph gives the
! quotient a constructor). partition_aggregate is the squint algebraic
! multigrid coarsens by: no part count asked for - the graph discovers
! its own aggregates by Vanek's three passes over the neighbour
! queries, strength having already chosen the edges.
!
! Escape times: escape_times resolves every vertex's escape time under
! a successor rule in one pass - orbits that merge share their tails
! instead of re-walking them - one visit per vertex where orbit-by-
! orbit painting spends one whole traversal per vertex.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module interface_graph

  use iso_fortran_env  , only : dp => REAL64
  use module_solve_mode, only : FORWARD, REVERSE

  implicit none

  private
  public :: graph, digraph, vertex, edge
  public :: counting_sort, transpose_adjacency

  type :: vertex
     integer :: number      ! vertex label
     integer :: part = 1    ! owning part (image) after partitioning
  end type vertex

  type :: edge
     integer :: tail, head  ! the two vertices the edge connects
  end type edge

  !===================================================================!
  ! Abstract graph
  !===================================================================!

  type, abstract :: graph

     integer :: num_vertices  = 0
     integer :: num_edges     = 0
     integer :: num_variables = 1   ! fields per vertex (dof interleaving)

     type(vertex), allocatable :: vertices(:)
     type(edge)  , allocatable :: edges(:)

     ! Retained adjacency: compressed neighbour list. Vertex v's
     ! neighbours are adj(xadj(v) : xadj(v+1)-1).
     integer, allocatable :: xadj(:)
     integer, allocatable :: adj(:)

     ! Partition bookkeeping, populated by partition / partition_rcb.
     ! part_of(v) = vertices(v) % part; owned/ghost vertices are stored
     ! csr-style: part k owns own_list(own_ptr(k):own_ptr(k+1)-1).
     integer              :: nparts = 1
     integer              :: ncut   = 0   ! edges crossing parts (cut quality)
     integer, allocatable :: own_ptr(:)   ! (nparts+1)
     integer, allocatable :: own_list(:)  ! (num_vertices)
     integer, allocatable :: gh_ptr(:)    ! (nparts+1)
     integer, allocatable :: gh_list(:)   ! deduped ghosts, by part

   contains

     ! dof indexing (variable-fastest over vertices), and the
     ! adjacency read at dof granularity
     procedure :: dof
     procedure :: num_dofs
     procedure :: dof_adjacency

     ! The contract: every graph answers its neighbour queries, and the
     ! compiler enforces it. A stored subclass implements them as
     ! one-line delegations to the stored_* helpers below; a
     ! rule-generated subclass answers by arithmetic.
     procedure(neighbours_interface), deferred :: neighbours
     procedure(degree_interface)    , deferred :: degree

     ! Shared stored-adjacency mechanism, written once: build the
     ! compressed neighbour list from the edge list, then index it.
     procedure :: build_adjacency
     procedure :: stored_neighbours
     procedure :: stored_degree

     ! visit order (breadth-first) over the whole graph
     procedure :: traversal_order

     ! orbit of a vertex under a caller-supplied successor rule
     procedure :: orbit

     ! every vertex's escape time under the rule, tails shared
     procedure :: escape_times

     ! partitioners - stamp vertex % part and gather the csr, in place.
     ! set_partition adopts parts computed elsewhere (a parent map, a
     ! metis call) through the same gathering.
     procedure :: partition
     procedure :: partition_rcb
     procedure :: partition_aggregate
     procedure :: set_partition

     ! the squint: the stamped partition read back as a coarse edge
     ! list. the zoom runs the other way: every vertex splits into
     ! children, and the refined edge list comes back.
     procedure :: quotient_edges
     procedure :: refine_edges

     ! partition queries
     procedure :: part_of
     procedure :: owned
     procedure :: ghosts
     procedure :: n_owned
     procedure :: n_ghosts
     procedure :: balance
     procedure :: edge_cut

     ! values riding the partition: any vertex list expands to its
     ! dofs, and the gather/scatter pair moves vector values through
     ! a part's owned dofs
     procedure :: dofs_of
     procedure :: gather
     procedure :: scatter

     procedure :: print
     procedure :: print_partition

     ! internal: stamp vertex % part, then gather the owned/ghost csr
     procedure, private :: stamp_bfs
     procedure, private :: stamp_rcb
     procedure, private :: stamp_aggregate
     procedure, private :: gather_partition

  end type graph

  !===================================================================!
  ! Deferred interfaces
  !===================================================================!

  abstract interface

     ! the neighbours of vertex v
     pure function neighbours_interface(this, v) result(nbrs)
       import :: graph
       class(graph), intent(in) :: this
       integer     , intent(in) :: v
       integer, allocatable     :: nbrs(:)
     end function neighbours_interface

     ! the number of neighbours of vertex v
     pure integer function degree_interface(this, v)
       import :: graph
       class(graph), intent(in) :: this
       integer     , intent(in) :: v
     end function degree_interface

  end interface

  !===================================================================!
  ! Directed graph: directedness is a genuine axis - it changes the
  ! protocol (out/in queries, dependency order, adjoint accumulation)
  ! and is orthogonal to storage-versus-rule, so it earns a type, not a
  ! flag and not a runtime refusal. A coupling graph (the mesh) stays on
  ! the base; calling the directed protocol on it does not compile.
  !===================================================================!

  type, abstract, extends(graph) :: digraph

     ! stored directed adjacency: two compressed lists built from the
     ! tail -> head edge list (a rule-generated subclass stores nothing)
     integer, allocatable :: out_xadj(:), out_adj(:)
     integer, allocatable :: in_xadj(:),  in_adj(:)

   contains

     ! the directed contract, compiler-enforced
     procedure(directed_neighbours_interface), deferred :: out_neighbours
     procedure(directed_neighbours_interface), deferred :: in_neighbours

     ! the base contract, provided: the deduplicated union of out and
     ! in, so partitioners and the visit order keep working
     procedure :: neighbours => union_neighbours
     procedure :: degree     => union_degree

     ! shared stored-directed mechanism (mirrors the base stored pattern)
     procedure :: build_directed_adjacency
     procedure :: stored_out_neighbours
     procedure :: stored_in_neighbours

     ! dependency semantics - true only here
     procedure :: dependency_order
     procedure :: is_acyclic

     ! the path from the single source, read as vertex numbers - the
     ! cargo of a route-shaped digraph
     procedure :: source_path

     ! reverse-mode accumulation of the discrete adjoint, direction
     ! read from structure
     procedure :: accumulate_adjoint

     ! the walk's witness, a contract method
     procedure, nopass :: verify_adjoint_accumulation

  end type digraph

  !===================================================================!
  ! Deferred interface of the directed contract
  !===================================================================!

  abstract interface

     ! the out- (or in-) neighbours of vertex v
     pure function directed_neighbours_interface(this, v) result(nbrs)
       import :: digraph
       class(digraph), intent(in) :: this
       integer       , intent(in) :: v
       integer, allocatable       :: nbrs(:)
     end function directed_neighbours_interface

  end interface

  !===================================================================!
  ! Apparatus of the walk's witness (verify_adjoint_accumulation), not
  ! demo code: a minimal stored digraph, private to this module, so the
  ! witness's fixtures cannot be chains-by-import.
  !===================================================================!

  type, extends(digraph) :: witness_digraph
   contains
     procedure :: out_neighbours => witness_out_neighbours
     procedure :: in_neighbours  => witness_in_neighbours
  end type witness_digraph

  ! witness fixture constants: the canonical recurrence chain and the
  ! diamond dag, each judged against its own named tolerance; the
  ! witness returns the worst defect normalized by these, so a value
  ! below one certifies both claims
  integer , parameter :: witness_chain_length      = 8
  real(dp), parameter :: witness_recurrence_factor = 0.7_dp
  real(dp), parameter :: witness_start_value       = 1.3_dp
  real(dp), parameter :: witness_chain_tolerance   = 1.0d-13
  real(dp), parameter :: witness_diamond_tolerance = 1.0d-6
  real(dp), parameter :: witness_factor_12 = 0.3_dp
  real(dp), parameter :: witness_factor_13 = 0.5_dp
  real(dp), parameter :: witness_factor_24 = 0.6_dp
  real(dp), parameter :: witness_factor_34 = 0.9_dp

contains

  !===================================================================!
  ! Global dof for a (vertex, variable) pair. Variable-fastest so a
  ! single field maps dof(v,1) = v.
  !===================================================================!

  pure type(integer) function dof(this, v, ivar)

    class(graph), intent(in) :: this
    integer     , intent(in) :: v, ivar

    dof = (v - 1)*this % num_variables + ivar

  end function dof

  !===================================================================!
  ! Total number of degrees of freedom in the graph
  !===================================================================!

  pure type(integer) function num_dofs(this)

    class(graph), intent(in) :: this

    num_dofs = this % num_vertices * this % num_variables

  end function num_dofs

  !===================================================================!
  ! The adjacency read at dof granularity: one row per dof, and each
  ! row leads with its own dof (the self-loop), then the matching
  ! variable's dof on every neighbour, in neighbour order. Emit the
  ! self pairs, emit the neighbour pairs, and the counting kernel
  ! groups them into rows. Consumes only the neighbour queries, so a
  ! rule-generated graph answers as readily as a stored one.
  !===================================================================!

  pure subroutine dof_adjacency(this, xadj, adj)

    class(graph)        , intent(in)  :: this
    integer, allocatable, intent(out) :: xadj(:), adj(:)

    integer, allocatable :: keys(:), values(:), nbrs(:)
    integer :: v, k, ivar, p, n, ndof

    ndof = this % num_dofs()

    n = ndof
    do v = 1, this % num_vertices
       n = n + this % degree(v) * this % num_variables
    end do
    allocate(keys(n), values(n))

    ! self pairs first - the kernel keeps arrival order within a row
    do p = 1, ndof
       keys(p)   = p
       values(p) = p
    end do

    n = ndof
    do v = 1, this % num_vertices
       nbrs = this % neighbours(v)
       do k = 1, size(nbrs)
          do ivar = 1, this % num_variables
             n = n + 1
             keys(n)   = this % dof(v, ivar)
             values(n) = this % dof(nbrs(k), ivar)
          end do
       end do
    end do

    call counting_sort(ndof, keys, values, xadj, adj)

  end subroutine dof_adjacency

  !===================================================================!
  ! Build the retained adjacency from the undirected edge list: per-
  ! vertex degrees accumulate into the row pointer xadj, then the
  ! neighbours scatter into adj in edge order (both directions).
  !===================================================================!

  pure subroutine build_adjacency(this)

    class(graph), intent(inout) :: this

    integer, allocatable :: keys(:), values(:)
    integer :: e

    ! every edge contributes both ways, in edge order - the same
    ! interleaving the hand-rolled scatter always produced
    allocate(keys(2*this % num_edges), values(2*this % num_edges))
    do e = 1, this % num_edges
       keys(2*e-1)   = this % edges(e) % tail
       values(2*e-1) = this % edges(e) % head
       keys(2*e)     = this % edges(e) % head
       values(2*e)   = this % edges(e) % tail
    end do

    call counting_sort(this % num_vertices, keys, values, this % xadj, this % adj)

  end subroutine build_adjacency

  !===================================================================!
  ! Counting sort - the one counting trick every retained structure
  ! in this module is built with. Given pairs (key, value) in order,
  ! group the values by key: ptr(k) points at key k's first value in
  ! sorted, and within a key the values keep their arrival order.
  ! Count, prefix-sum, scatter - that is the whole algorithm.
  !===================================================================!

  pure subroutine counting_sort(nkeys, keys, values, ptr, sorted)

    integer             , intent(in)  :: nkeys
    integer             , intent(in)  :: keys(:)
    integer             , intent(in)  :: values(:)
    integer, allocatable, intent(out) :: ptr(:)
    integer, allocatable, intent(out) :: sorted(:)

    integer, allocatable :: cursor(:)
    integer :: i, k

    allocate(ptr(nkeys+1))
    ptr = 0
    do i = 1, size(keys)
       ptr(keys(i)+1) = ptr(keys(i)+1) + 1
    end do

    ptr(1) = 1
    do k = 1, nkeys
       ptr(k+1) = ptr(k+1) + ptr(k)
    end do

    allocate(sorted(ptr(nkeys+1)-1))
    allocate(cursor(nkeys))
    cursor = ptr(1:nkeys)
    do i = 1, size(keys)
       sorted(cursor(keys(i))) = values(i)
       cursor(keys(i)) = cursor(keys(i)) + 1
    end do

  end subroutine counting_sort

  !===================================================================!
  ! Transpose an adjacency. One side of a bipartite graph lists its
  ! neighbours: for each key, the values it touches, as a padded
  ! array with per-key counts. This returns the other side's lists -
  ! for each value, the keys that touch it. The transpose is built on
  ! counting_sort, like every retained structure in this module.
  ! Values are numbered 1..n_values.
  !===================================================================!

  pure subroutine transpose_adjacency(forward, num_forward, n_values, reverse, num_reverse)

    integer             , intent(in)  :: forward(:,:)
    integer             , intent(in)  :: num_forward(:)
    integer             , intent(in)  :: n_values
    integer, allocatable, intent(out) :: reverse(:,:)
    integer, allocatable, intent(out) :: num_reverse(:)

    integer, allocatable :: keys(:), values(:), ptr(:), sorted(:)
    integer :: key, k, n, v

    ! flatten the padded lists into (value, key) pairs, in key order
    allocate(keys(sum(num_forward)), values(sum(num_forward)))
    n = 0
    do key = 1, size(num_forward)
       do k = 1, num_forward(key)
          n = n + 1
          keys(n)   = forward(k, key)
          values(n) = key
       end do
    end do

    call counting_sort(n_values, keys, values, ptr, sorted)

    ! pad the grouped rows back into the rectangular shape
    allocate(num_reverse(n_values))
    do v = 1, n_values
       num_reverse(v) = ptr(v+1) - ptr(v)
    end do
    allocate(reverse(maxval(num_reverse), n_values))
    reverse = 0
    do v = 1, n_values
       reverse(1:num_reverse(v), v) = sorted(ptr(v) : ptr(v+1)-1)
    end do

  end subroutine transpose_adjacency

  !===================================================================!
  ! Neighbours of vertex v from the retained adjacency: the shared
  ! mechanism a stored subclass delegates its deferred neighbours to.
  ! The stop below is a broken-invariant report, not an extension hook:
  ! a stored graph's constructor must call build_adjacency.
  !===================================================================!

  pure function stored_neighbours(this, v) result(nbrs)

    class(graph), intent(in) :: this
    integer     , intent(in) :: v

    integer, allocatable :: nbrs(:)

    if (.not. allocated(this % xadj)) then
       error stop "graph: stored graph constructed without adjacency - " // &
            & "the constructor must call build_adjacency"
    end if

    nbrs = this % adj(this % xadj(v) : this % xadj(v+1)-1)

  end function stored_neighbours

  !===================================================================!
  ! Degree of vertex v from the retained adjacency (see
  ! stored_neighbours for the invariant).
  !===================================================================!

  pure integer function stored_degree(this, v)

    class(graph), intent(in) :: this
    integer     , intent(in) :: v

    if (.not. allocated(this % xadj)) then
       error stop "graph: stored graph constructed without adjacency - " // &
            & "the constructor must call build_adjacency"
    end if

    stored_degree = this % xadj(v+1) - this % xadj(v)

  end function stored_degree

  !===================================================================!
  ! Visit order over all vertices: breadth-first, restarting for
  ! disconnected components (FORWARD); the same order reversed
  ! (REVERSE). This is locality for the partitioner, not a dependency
  ! claim - dependency semantics live on digraph. Consumes only the
  ! neighbour queries.
  !===================================================================!

  pure function traversal_order(this, direction) result(order)

    class(graph), intent(in)           :: this
    integer     , intent(in), optional :: direction   ! FORWARD (default) / REVERSE

    integer, allocatable :: order(:)

    integer, allocatable :: queue(:), nbrs(:)
    logical, allocatable :: seen(:)
    integer :: nv, i, v, w, pos, qh, qt, start, dir

    dir = FORWARD
    if (present(direction)) dir = direction

    nv = this % num_vertices
    allocate(order(nv), queue(nv), seen(nv))
    seen = .false.

    pos = 0
    qh  = 1
    qt  = 0

    do start = 1, nv

       if (seen(start)) cycle

       qt = qt + 1
       queue(qt)   = start
       seen(start) = .true.

       do while (qh .le. qt)

          v  = queue(qh)
          qh = qh + 1

          pos = pos + 1
          order(pos) = v

          nbrs = this % neighbours(v)
          do i = 1, size(nbrs)
             w = nbrs(i)
             if (.not. seen(w)) then
                seen(w) = .true.
                qt = qt + 1
                queue(qt) = w
             end if
          end do

       end do

    end do

    if (dir .eq. REVERSE) order = order(nv:1:-1)

  end function traversal_order

  !===================================================================!
  ! Orbit of a vertex under a successor rule: one vertex, one arrow
  ! out, repeated. The rule gives every vertex its single successor; a
  ! value outside 1..num_vertices means the arrow leaves the graph.
  !
  ! The returned sequence starts at start and ends when the rule
  ! leaves the vertex set (the orbit escapes; the length is the escape
  ! time), when a vertex reappears (a cycle has closed; the repeated
  ! vertex is kept as the final entry, so the cycle is the slice
  ! between its two appearances), or when limit entries have been
  ! recorded (a capped sequence makes no claim about the fate beyond
  ! it). The default limit num_vertices + 1 guarantees termination on
  ! its own: a rule that never escapes a finite vertex set must repeat
  ! a vertex by then. Consumes only the rule - the graph contributes
  ! the vertex set the orbit may escape.
  !===================================================================!

  pure function orbit(this, start, successor, limit) result(visited)

    class(graph), intent(in) :: this
    integer     , intent(in) :: start
    interface
       ! the single successor of vertex v under the rule
       pure integer function successor(v)
         integer, intent(in) :: v
       end function successor
    end interface
    integer, intent(in), optional :: limit

    integer, allocatable :: visited(:)

    integer, allocatable :: scratch(:)
    logical, allocatable :: seen(:)
    integer :: nv, cap, v, n

    nv = this % num_vertices

    if (start .lt. 1 .or. start .gt. nv) then
       error stop "graph: orbit must start at a vertex of the graph"
    end if

    cap = nv + 1
    if (present(limit)) cap = limit

    allocate(scratch(cap), seen(nv))
    seen = .false.

    v = start
    n = 0
    do while (n .lt. cap)
       n = n + 1
       scratch(n) = v
       if (seen(v)) exit                    ! a cycle has closed
       seen(v) = .true.
       v = successor(v)
       if (v .lt. 1 .or. v .gt. nv) exit    ! the orbit escapes the graph
    end do

    visited = scratch(1:n)

  end function orbit

  !===================================================================!
  ! Every vertex's escape time under a successor rule, in one pass.
  ! The trick: orbits overlap. Two orbits that merge share their whole
  ! tail, so nothing is walked twice - follow the rule from an
  ! unresolved vertex until it steps off the graph (that last vertex
  ! escapes in 1, and the times count back up the path), reaches a
  ! vertex already resolved (count from its time), or bites its own
  ! path (a cycle: everyone on it is home). Times are capped at limit
  ! and home is painted limit, exactly the picture the orbit method
  ! paints one vertex at a time - but at one visit per vertex.
  !===================================================================!

  pure function escape_times(this, successor, limit) result(times)

    class(graph), intent(in) :: this
    interface
       ! the single successor of vertex v under the rule
       pure integer function successor(v)
         integer, intent(in) :: v
       end function successor
    end interface
    integer, intent(in) :: limit

    integer, allocatable :: times(:)

    integer, allocatable :: path(:), pos(:)
    integer :: nv, start, v, depth, base, i

    if (limit .lt. 1) then
       error stop "graph: escape_times requires a positive limit"
    end if

    nv = this % num_vertices
    allocate(times(nv), path(nv), pos(nv))
    times = 0            ! 0 marks unresolved; every resolved time is >= 1
    pos   = 0

    do start = 1, nv

       if (times(start) .ne. 0) cycle

       depth = 0
       v     = start

       follow: do
          if (v .lt. 1 .or. v .gt. nv) then
             base = 0                    ! stepped off: the path's last vertex escapes in 1
             exit follow
          end if
          if (times(v) .ne. 0) then
             base = times(v)             ! already resolved: count on from its time
             exit follow
          end if
          if (pos(v) .ne. 0) then
             ! the path bit itself: a cycle - everyone on it is home
             do i = pos(v), depth
                times(path(i)) = limit
             end do
             base  = limit
             depth = pos(v) - 1          ! the tail feeding the cycle remains
             exit follow
          end if
          depth       = depth + 1
          path(depth) = v
          pos(v)      = depth
          v = successor(v)
       end do follow

       ! stale pos marks need no clearing: every marked vertex leaves
       ! this pass resolved, and the resolved check runs first. the cap
       ! is applied to the step count before adding, so the sum never
       ! exceeds limit - base + limit would overflow near huge(0)
       do i = depth, 1, -1
          times(path(i)) = base + min(depth - i + 1, limit - base)
       end do

    end do

  end function escape_times

  !===================================================================!
  ! Partition into nparts parts by breadth-first ordering, in place
  ! (vertex % part stamped, owned/ghost csr gathered).
  !
  ! This is the serial placeholder: bfs keeps neighbours close so the
  ! chunks stay mostly connected and the cut is reasonable. (KB: swap the
  ! stamp_bfs body for a parmetis call - same contract.) Sibling to
  ! partition_rcb (geometric).
  !===================================================================!

  pure subroutine partition(this, nparts)

    class(graph), intent(inout) :: this
    integer     , intent(in)    :: nparts

    call this % stamp_bfs(nparts)
    call this % gather_partition(nparts)

  end subroutine partition

  !===================================================================!
  ! Geometric partitioner: recursive coordinate bisection (RCB), in
  ! place. Uses the vertex coordinates coords(:,v) (e.g. cell centroids)
  ! to recursively split by the median coordinate along the longest-
  ! spread axis - balanced, COMPACT subdomains (small edge cut / halo)
  ! where the BFS chunks would go stringy. (KB: METIS could still slot
  ! in behind this same contract.)
  !===================================================================!

  impure subroutine partition_rcb(this, coords, nparts)

    class(graph), intent(inout) :: this
    real(dp)    , intent(in)    :: coords(:,:)   ! (ndim, num_vertices)
    integer     , intent(in)    :: nparts

    call this % stamp_rcb(coords, nparts)
    call this % gather_partition(nparts)

  end subroutine partition_rcb

  !===================================================================!
  ! Partition by greedy aggregation (Vanek's three passes), in place.
  ! No part count is asked for: the graph discovers its own. A vertex
  ! whose whole neighbourhood is still free seeds a part made of
  ! itself and its neighbours; leftovers sweep into a neighbouring
  ! part; whatever remains seeds its own. This is the squint algebraic
  ! multigrid coarsens by - strength has already chosen the edges by
  ! the time this graph exists, so the passes consume only the
  ! neighbour queries. Sibling of partition (breadth-first) and
  ! partition_rcb (geometric).
  !===================================================================!

  pure subroutine partition_aggregate(this)

    class(graph), intent(inout) :: this

    integer :: nparts

    call this % stamp_aggregate(nparts)
    call this % gather_partition(nparts)

  end subroutine partition_aggregate

  !===================================================================!
  ! Adopt a partition computed elsewhere: parts(v) names each vertex's
  ! part, numbered from 1. The same stamping and gathering the built-in
  ! partitioners use - this is the doorway for a parent map from
  ! refinement, or for an external partitioner like metis.
  !===================================================================!

  pure subroutine set_partition(this, parts)

    class(graph), intent(inout) :: this
    integer     , intent(in)    :: parts(:)

    integer :: v

    if (size(parts) .ne. this % num_vertices) then
       error stop "graph: set_partition needs one part per vertex"
    end if
    if (minval(parts) .lt. 1) then
       error stop "graph: parts are numbered from 1"
    end if

    do v = 1, this % num_vertices
       this % vertices(v) % part = parts(v)
    end do

    call this % gather_partition(maxval(parts))

  end subroutine set_partition

  !===================================================================!
  ! Zoom in: every vertex splits into the given number of children,
  ! and the refined edge list comes back. Children of vertex v are
  ! numbered (v-1)*children+1 .. v*children - the parent of any child
  ! is arithmetic. Within a split vertex the children form a chain;
  ! each original edge becomes one bridge between the first children
  ! of its two ends. Consumes only the neighbour queries, so a
  ! rule-generated graph refines as readily as a stored one. The
  ! quotient of the refinement by the parent map is the original
  ! graph - the zoom and the squint undo each other.
  !===================================================================!

  pure subroutine refine_edges(this, children, tails, heads)

    class(graph)        , intent(in)  :: this
    integer             , intent(in)  :: children
    integer, allocatable, intent(out) :: tails(:), heads(:)

    integer, allocatable :: nbrs(:)
    integer              :: v, w, i, c, ne, n_refined

    if (children .lt. 1) then
       error stop "graph: refine_edges needs at least one child per vertex"
    end if

    ! chains within each split vertex, one bridge per original edge
    n_refined = this % num_vertices*(children - 1) + this % num_edges
    allocate(tails(n_refined), heads(n_refined))

    ne = 0
    do v = 1, this % num_vertices
       do c = 1, children - 1
          ne        = ne + 1
          tails(ne) = (v-1)*children + c
          heads(ne) = (v-1)*children + c + 1
       end do
       nbrs = this % neighbours(v)
       do i = 1, size(nbrs)
          w = nbrs(i)
          if (w .gt. v) then
             ne        = ne + 1
             tails(ne) = (v-1)*children + 1
             heads(ne) = (w-1)*children + 1
          end if
       end do
    end do

    tails = tails(1:ne)
    heads = heads(1:ne)

  end subroutine refine_edges

  !===================================================================!
  ! Stand back and squint: read the stamped partition back as a graph.
  ! Every part becomes one coarse vertex; two parts joined by at least
  ! one fine edge become coarse neighbours - one coarse edge, however
  ! many fine edges cross between them. Returned as an edge list
  ! because a partition of a graph is still graph-shaped data; the
  ! stored_graph constructor turns it back into the same animal, ready
  ! to be squinted again.
  !===================================================================!

  pure subroutine quotient_edges(this, tails, heads)

    class(graph)        , intent(in)  :: this
    integer, allocatable, intent(out) :: tails(:), heads(:)

    integer, allocatable :: mark(:), nbrs(:)
    integer :: k, i, v, w, kw, n_coarse, pass

    if (.not. allocated(this % own_ptr)) then
       error stop "graph: quotient requires a stamped partition"
    end if

    ! two passes over each part's owned vertices and their neighbours:
    ! count the coarse edges, then fill them, deduped by stamping
    ! mark(kw) = k (k is monotone, so an old stamp reads as unmarked).
    ! each coarse edge is recorded once, from its smaller part.
    allocate(mark(this % nparts))

    do pass = 1, 2
       mark     = 0
       n_coarse = 0
       do k = 1, this % nparts
          do i = this % own_ptr(k), this % own_ptr(k+1)-1
             v    = this % own_list(i)
             nbrs = this % neighbours(v)
             do w = 1, size(nbrs)
                kw = this % vertices(nbrs(w)) % part
                if (kw .gt. k .and. mark(kw) .ne. k) then
                   mark(kw) = k
                   n_coarse = n_coarse + 1
                   if (pass .eq. 2) then
                      tails(n_coarse) = k
                      heads(n_coarse) = kw
                   end if
                end if
             end do
          end do
       end do
       if (pass .eq. 1) allocate(tails(n_coarse), heads(n_coarse))
    end do

  end subroutine quotient_edges

  !===================================================================!
  ! Stamp vertex % part by cutting the forward traversal order into
  ! nparts balanced contiguous chunks.
  !===================================================================!

  pure subroutine stamp_bfs(this, nparts)

    class(graph) , intent(inout) :: this
    integer      , intent(in)    :: nparts

    integer, allocatable :: order(:)
    integer :: nv, pos

    nv = this % num_vertices
    allocate(order(nv))
    order = this % traversal_order(FORWARD)

    do pos = 1, nv
       this % vertices(order(pos)) % part = (pos-1)*nparts/nv + 1
    end do

  end subroutine stamp_bfs

  !===================================================================!
  ! Stamp vertex % part by Vanek's three greedy passes, discovering
  ! the part count. Pass 1: a vertex whose whole neighbourhood is
  ! still free seeds a part of itself and its neighbours. Pass 2: each
  ! leftover sweeps into the first neighbouring part it sees. Pass 3:
  ! anything still alone seeds its own singleton part.
  !===================================================================!

  pure subroutine stamp_aggregate(this, nparts)

    class(graph), intent(inout) :: this
    integer     , intent(out)   :: nparts

    integer, allocatable :: part(:), nbrs(:)
    integer :: nv, v, i
    logical :: free

    nv = this % num_vertices
    allocate(part(nv))
    part   = 0
    nparts = 0

    ! pass 1: seed where the whole neighbourhood is free
    do v = 1, nv
       if (part(v) .ne. 0) cycle
       nbrs = this % neighbours(v)
       free = .true.
       do i = 1, size(nbrs)
          if (part(nbrs(i)) .ne. 0) then
             free = .false.
             exit
          end if
       end do
       if (.not. free) cycle
       nparts  = nparts + 1
       part(v) = nparts
       do i = 1, size(nbrs)
          part(nbrs(i)) = nparts
       end do
    end do

    ! pass 2: sweep each leftover into a neighbouring part
    do v = 1, nv
       if (part(v) .ne. 0) cycle
       nbrs = this % neighbours(v)
       do i = 1, size(nbrs)
          if (part(nbrs(i)) .ne. 0) then
             part(v) = part(nbrs(i))
             exit
          end if
       end do
    end do

    ! pass 3: whoever is still alone seeds a singleton
    do v = 1, nv
       if (part(v) .eq. 0) then
          nparts  = nparts + 1
          part(v) = nparts
       end if
    end do

    do v = 1, nv
       this % vertices(v) % part = part(v)
    end do

  end subroutine stamp_aggregate

  !===================================================================!
  ! Stamp vertex % part by recursive coordinate bisection of the
  ! coordinates.
  !===================================================================!

  impure subroutine stamp_rcb(this, coords, nparts)

    class(graph), intent(inout) :: this
    real(dp)    , intent(in)    :: coords(:,:)
    integer     , intent(in)    :: nparts

    integer, allocatable :: idx(:), part(:)
    integer :: nv, c

    nv = this % num_vertices
    if (nparts .lt. 1)  error stop "partition_rcb: nparts < 1"
    if (nparts .gt. nv) error stop "partition_rcb: more parts than vertices"

    allocate(idx(nv), part(nv))
    do c = 1, nv
       idx(c) = c
    end do
    part = 1

    call rcb(idx, 1, nv, nparts, 1, coords, part)

    do c = 1, nv
       this % vertices(c) % part = part(c)
    end do

  end subroutine stamp_rcb

  !===================================================================!
  ! Gather the owned/ghost csr from the already-stamped vertex % part.
  ! ghost(k) = vertices not in k but adjacent to a k-owned vertex - the
  ! halo part k must pull from the vertices' owners. Consumes only the
  ! neighbour queries.
  !===================================================================!

  pure subroutine gather_partition(this, nparts)

    class(graph), intent(inout) :: this
    integer     , intent(in)    :: nparts

    integer, allocatable :: cnt(:), mark(:), nbrs(:)
    integer :: nv, v, w, i, k, kv, pos

    nv = this % num_vertices
    this % nparts = nparts

    ! ---- owned vertices per part: counting sort by part ----
    call counting_sort(nparts, [(this % vertices(v) % part, v = 1, nv)], &
         &             [(v, v = 1, nv)], this % own_ptr, this % own_list)

    ! ---- edge cut ----
    this % ncut = this % edge_cut()

    ! ---- ghost vertices per part ----
    ! Two passes over each part's owned vertices and their neighbours:
    ! count (size the csr) then fill, deduped by stamping mark(w)=k (k is
    ! monotone, so an old stamp from part k-1 reads as unmarked for k).
    allocate(cnt(nparts)); cnt = 0
    allocate(mark(nv));    mark = 0
    do k = 1, nparts
       do i = this % own_ptr(k), this % own_ptr(k+1)-1
          v    = this % own_list(i)
          nbrs = this % neighbours(v)
          do kv = 1, size(nbrs)
             w = nbrs(kv)
             if (this % vertices(w) % part .ne. k .and. mark(w) .ne. k) then
                mark(w) = k
                cnt(k)  = cnt(k) + 1
             end if
          end do
       end do
    end do

    if (allocated(this % gh_ptr))  deallocate(this % gh_ptr)
    if (allocated(this % gh_list)) deallocate(this % gh_list)
    allocate(this % gh_ptr(nparts+1))
    this % gh_ptr(1) = 1
    do k = 1, nparts
       this % gh_ptr(k+1) = this % gh_ptr(k) + cnt(k)
    end do
    allocate(this % gh_list(this % gh_ptr(nparts+1)-1))

    mark = 0
    do k = 1, nparts
       pos = this % gh_ptr(k)
       do i = this % own_ptr(k), this % own_ptr(k+1)-1
          v    = this % own_list(i)
          nbrs = this % neighbours(v)
          do kv = 1, size(nbrs)
             w = nbrs(kv)
             if (this % vertices(w) % part .ne. k .and. mark(w) .ne. k) then
                mark(w) = k
                this % gh_list(pos) = w
                pos = pos + 1
             end if
          end do
       end do
    end do

  end subroutine gather_partition

  !===================================================================!
  ! Owning part of a vertex
  !===================================================================!

  pure integer function part_of(this, v)

    class(graph), intent(in) :: this
    integer     , intent(in) :: v

    part_of = this % vertices(v) % part

  end function part_of

  !===================================================================!
  ! Vertices owned by part k
  !===================================================================!

  pure function owned(this, k) result(list)

    class(graph), intent(in) :: this
    integer     , intent(in) :: k

    integer, allocatable :: list(:)

    list = this % own_list(this % own_ptr(k) : this % own_ptr(k+1)-1)

  end function owned

  !===================================================================!
  ! The dofs a list of vertices carries, variable-fastest - the same
  ! interleaving dof(v, ivar) uses. Works for any list the graph
  ! hands out: a part's owned vertices, its ghosts, an aggregate.
  !===================================================================!

  pure function dofs_of(this, verts) result(dofs)

    class(graph), intent(in) :: this
    integer     , intent(in) :: verts(:)

    integer, allocatable :: dofs(:)
    integer :: i, ivar, pos

    allocate(dofs(size(verts) * this % num_variables))

    pos = 0
    do i = 1, size(verts)
       do ivar = 1, this % num_variables
          pos = pos + 1
          dofs(pos) = this % dof(verts(i), ivar)
       end do
    end do

  end function dofs_of

  !===================================================================!
  ! Gather: pull the values x carries at part k's owned dofs into a
  ! dense local vector.
  !===================================================================!

  pure function gather(this, k, x) result(xk)

    class(graph), intent(in) :: this
    integer     , intent(in) :: k
    real(dp)    , intent(in) :: x(:)

    real(dp), allocatable :: xk(:)

    xk = x(this % dofs_of(this % owned(k)))

  end function gather

  !===================================================================!
  ! Scatter: the same sentence read backwards - push a local vector
  ! into x at part k's owned dofs, leaving every other entry alone.
  !===================================================================!

  pure subroutine scatter(this, k, xk, x)

    class(graph), intent(in)    :: this
    integer     , intent(in)    :: k
    real(dp)    , intent(in)    :: xk(:)
    real(dp)    , intent(inout) :: x(:)

    x(this % dofs_of(this % owned(k))) = xk

  end subroutine scatter

  !===================================================================!
  ! Halo (ghost) vertices part k needs from other parts
  !===================================================================!

  pure function ghosts(this, k) result(list)

    class(graph), intent(in) :: this
    integer     , intent(in) :: k

    integer, allocatable :: list(:)

    list = this % gh_list(this % gh_ptr(k) : this % gh_ptr(k+1)-1)

  end function ghosts

  pure integer function n_owned(this, k)

    class(graph), intent(in) :: this
    integer     , intent(in) :: k

    n_owned = this % own_ptr(k+1) - this % own_ptr(k)

  end function n_owned

  pure integer function n_ghosts(this, k)

    class(graph), intent(in) :: this
    integer     , intent(in) :: k

    n_ghosts = this % gh_ptr(k+1) - this % gh_ptr(k)

  end function n_ghosts

  !===================================================================!
  ! Load-balance ratio max(owned)/min(owned) over the parts (1.0 = perfect)
  !===================================================================!

  pure real(dp) function balance(this)

    class(graph), intent(in) :: this

    integer :: k, lo, hi, m

    lo = huge(1); hi = 0
    do k = 1, this % nparts
       m = this % n_owned(k)
       lo = min(lo, m); hi = max(hi, m)
    end do

    if (lo .le. 0) then
       balance = huge(1.0_dp)
    else
       balance = real(hi, dp)/real(lo, dp)
    end if

  end function balance

  !===================================================================!
  ! Number of edges whose endpoints lie in different parts. Counted over
  ! the neighbour queries (each undirected edge visited from both ends;
  ! the w > v guard counts it once).
  !===================================================================!

  pure type(integer) function edge_cut(this)

    class(graph), intent(in) :: this

    integer, allocatable :: nbrs(:)
    integer :: v, i, w

    edge_cut = 0

    do v = 1, this % num_vertices
       nbrs = this % neighbours(v)
       do i = 1, size(nbrs)
          w = nbrs(i)
          if (w .gt. v .and. &
               & this % vertices(v) % part .ne. this % vertices(w) % part) then
             edge_cut = edge_cut + 1
          end if
       end do
    end do

  end function edge_cut

  !===================================================================!
  ! Print a one line summary of the graph
  !===================================================================!

  impure subroutine print(this)

    class(graph), intent(in) :: this

    write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') &
         & "graph: ", this % num_vertices, " vertices, ", this % num_edges, &
         & " edges, ", this % num_variables, " variables/vertex, ", this % num_dofs(), " dofs"

  end subroutine print

  !===================================================================!
  ! Print a summary of the partition (parts, vertices, edge cut, balance).
  !===================================================================!

  impure subroutine print_partition(this)

    class(graph), intent(in) :: this

    integer :: k

    write(*,'(1x,a,i0,a,i0,a,i0,a,f6.3)') "partition: ", this % nparts, &
         & " parts, ", this % num_vertices, " cells, edge cut ", this % ncut, &
         & ", balance ", this % balance()
    do k = 1, this % nparts
       write(*,'(3x,a,i0,a,i0,a,i0,a)') "part ", k, ": ", this % n_owned(k), &
            & " owned, ", this % n_ghosts(k), " ghost"
    end do

  end subroutine print_partition

  !===================================================================!
  ! The base contract on a digraph: neighbours is the deduplicated
  ! union of out- and in-neighbours, so the partitioners and the visit
  ! order consume a directed graph unchanged.
  !===================================================================!

  pure function union_neighbours(this, v) result(nbrs)

    class(digraph), intent(in) :: this
    integer       , intent(in) :: v

    integer, allocatable :: nbrs(:), outn(:), inn(:)
    integer :: n, i

    outn = this % out_neighbours(v)
    inn  = this % in_neighbours(v)

    ! deduplicate within each list as well as across them (a multigraph
    ! may carry repeated edges)
    allocate(nbrs(size(outn) + size(inn)))
    n = 0
    do i = 1, size(outn)
       if (.not. any(nbrs(1:n) .eq. outn(i))) then
          n = n + 1
          nbrs(n) = outn(i)
       end if
    end do
    do i = 1, size(inn)
       if (.not. any(nbrs(1:n) .eq. inn(i))) then
          n = n + 1
          nbrs(n) = inn(i)
       end if
    end do

    nbrs = nbrs(1:n)

  end function union_neighbours

  pure integer function union_degree(this, v)

    class(digraph), intent(in) :: this
    integer       , intent(in) :: v

    integer, allocatable :: nbrs(:)

    nbrs = union_neighbours(this, v)
    union_degree = size(nbrs)

  end function union_degree

  !===================================================================!
  ! Build the stored directed adjacency: two compressed lists from the
  ! tail -> head edge list, out-rows by tail and in-rows by head, both
  ! in edge order.
  !===================================================================!

  pure subroutine build_directed_adjacency(this)

    class(digraph), intent(inout) :: this

    integer, allocatable :: tails(:), heads(:)

    ! the out-lists sort the heads by tail; the in-lists sort the
    ! tails by head - counting sort both ways
    tails = this % edges % tail
    heads = this % edges % head

    call counting_sort(this % num_vertices, tails, heads, this % out_xadj, this % out_adj)
    call counting_sort(this % num_vertices, heads, tails, this % in_xadj, this % in_adj)

  end subroutine build_directed_adjacency

  !===================================================================!
  ! Stored directed queries: the shared mechanism a stored directed
  ! subclass delegates its deferred procedures to. The stop is a
  ! broken-invariant report - a stored digraph's constructor must call
  ! build_directed_adjacency - not an extension hook.
  !===================================================================!

  pure function stored_out_neighbours(this, v) result(nbrs)

    class(digraph), intent(in) :: this
    integer       , intent(in) :: v

    integer, allocatable :: nbrs(:)

    if (.not. allocated(this % out_xadj)) then
       error stop "digraph: stored digraph constructed without directed adjacency - " // &
            & "the constructor must call build_directed_adjacency"
    end if

    nbrs = this % out_adj(this % out_xadj(v) : this % out_xadj(v+1)-1)

  end function stored_out_neighbours

  pure function stored_in_neighbours(this, v) result(nbrs)

    class(digraph), intent(in) :: this
    integer       , intent(in) :: v

    integer, allocatable :: nbrs(:)

    if (.not. allocated(this % in_xadj)) then
       error stop "digraph: stored digraph constructed without directed adjacency - " // &
            & "the constructor must call build_directed_adjacency"
    end if

    nbrs = this % in_adj(this % in_xadj(v) : this % in_xadj(v+1)-1)

  end function stored_in_neighbours

  !===================================================================!
  ! Topological order over the out-edges (Kahn's construction): vertices
  ! with no incoming edges enter first, in ascending vertex order for
  ! determinism. n_ordered < num_vertices signals a directed cycle.
  ! Shared by dependency_order (which refuses on a cycle) and is_acyclic
  ! (which reports without dying).
  !===================================================================!

  pure subroutine topological_order(this, order, n_ordered)

    class(digraph), intent(in)  :: this
    integer, allocatable, intent(out) :: order(:)
    integer              , intent(out) :: n_ordered

    integer, allocatable :: indeg(:), queue(:), nbrs(:)
    integer :: nv, v, w, i, qh, qt

    nv = this % num_vertices
    allocate(order(nv), indeg(nv), queue(nv))

    do v = 1, nv
       nbrs = this % in_neighbours(v)
       indeg(v) = size(nbrs)
    end do

    qh = 1; qt = 0
    do v = 1, nv
       if (indeg(v) .eq. 0) then
          qt = qt + 1
          queue(qt) = v
       end if
    end do

    n_ordered = 0
    do while (qh .le. qt)
       v  = queue(qh)
       qh = qh + 1
       n_ordered = n_ordered + 1
       order(n_ordered) = v
       nbrs = this % out_neighbours(v)
       do i = 1, size(nbrs)
          w = nbrs(i)
          indeg(w) = indeg(w) - 1
          if (indeg(w) .eq. 0) then
             qt = qt + 1
             queue(qt) = w
          end if
       end do
    end do

  end subroutine topological_order

  !===================================================================!
  ! Dependency order over the out-edges: FORWARD as computed, REVERSE
  ! flipped. On a cycle it refuses - standing wiring must be directed
  ! and loop-free before a number moves; this stop is the structural
  ! audit, reachable only by an ill-formed graph, not a stub. Use
  ! is_acyclic to test without dying.
  !===================================================================!

  pure function dependency_order(this, direction) result(order)

    class(digraph), intent(in)           :: this
    integer       , intent(in), optional :: direction   ! FORWARD (default) / REVERSE

    integer, allocatable :: order(:)
    integer :: n_ordered, dir

    dir = FORWARD
    if (present(direction)) dir = direction

    call topological_order(this, order, n_ordered)

    if (n_ordered .lt. this % num_vertices) then
       error stop "digraph: dependency order requires acyclic directed wiring"
    end if

    if (dir .eq. REVERSE) order = order(this % num_vertices : 1 : -1)

  end function dependency_order

  !===================================================================!
  ! The path from the source, read as numbers. Some digraphs are
  ! shaped like a route: one vertex with no arrow in, at most one
  ! arrow out of every vertex. Follow the arrows from the source and
  ! return each vertex's number in the order visited. The numbers are
  ! the cargo - a multigrid schedule names levels with them - but the
  ! path neither knows nor cares what they mean. A digraph that is
  ! not a route stops with a report: two sources or none, a fork, or
  ! a vertex seen twice (with one arrow out of each, a second visit
  ! means the path never ends).
  !===================================================================!

  pure function source_path(this) result(numbers)

    class(digraph), intent(in) :: this

    integer, allocatable :: numbers(:)

    integer, allocatable :: visits(:), nbrs(:)
    integer              :: v, s, n_source, n

    ! the single source: the one vertex with no arrow in
    n_source = 0
    s        = 0
    do v = 1, this % num_vertices
       if (size(this % in_neighbours(v)) .eq. 0) then
          n_source = n_source + 1
          s        = v
       end if
    end do
    if (n_source .ne. 1) then
       error stop "digraph: source_path needs exactly one source"
    end if

    allocate(visits(this % num_vertices))
    n = 0
    do
       n = n + 1
       if (n .gt. this % num_vertices) then
          error stop "digraph: source_path revisits a vertex - the path never ends"
       end if
       visits(n) = this % vertices(s) % number
       nbrs = this % out_neighbours(s)
       if (size(nbrs) .gt. 1) then
          error stop "digraph: source_path hit a fork - one arrow out, not several"
       end if
       if (size(nbrs) .eq. 0) exit
       s = nbrs(1)
    end do

    ! the path must cover the graph: vertices it never reached sit in
    ! a detached component with no source of its own - a cycle the
    ! single-source audit cannot see
    if (n .lt. this % num_vertices) then
       error stop "digraph: source_path left vertices unvisited - the route does not cover the graph"
    end if

    numbers = visits(1:n)

  end function source_path

  !===================================================================!
  ! True when the directed wiring has no cycle (pure report, no stop) -
  ! for constructors and tests to assert without dying.
  !===================================================================!

  pure logical function is_acyclic(this)

    class(digraph), intent(in) :: this

    integer, allocatable :: order(:)
    integer :: n_ordered

    call topological_order(this, order, n_ordered)
    is_acyclic = (n_ordered .eq. this % num_vertices)

  end function is_acyclic

  !===================================================================!
  ! Reverse-mode accumulation - the structure of the discrete adjoint,
  ! implemented once, on the directed graph. Direction is read from the
  ! structure: vertices are taken in dependency order (which proves
  ! acyclicity before a number moves), the adjoint is seeded at the
  ! objective vertex `at` (no hidden assumption about where the
  ! objective sits), and each earlier vertex accumulates its
  ! out-neighbours' adjoints through the caller-supplied edge_apply.
  ! Vertices after `at` in the dependency order correctly remain zero -
  ! they do not influence the objective. The graph never evaluates a
  ! derivative itself; edge_apply is typically an internal procedure of
  ! the caller with access to the required system state.
  !
  ! adjoint is (block, num_vertices): one block-sized vector per vertex
  ! (block = 1 for scalar recurrences, the state size for iterate
  ! chains).
  !===================================================================!

  subroutine accumulate_adjoint(this, seed, edge_apply, at, adjoint)

    class(digraph), intent(in)  :: this
    real(dp)      , intent(in)  :: seed(:)
    interface
       subroutine edge_apply(tail, head, adjoint_head, contribution)
         import :: dp
         integer , intent(in)  :: tail, head
         real(dp), intent(in)  :: adjoint_head(:)
         real(dp), intent(out) :: contribution(:)
       end subroutine edge_apply
    end interface
    integer       , intent(in)  :: at          ! the objective vertex
    real(dp)      , intent(out) :: adjoint(:,:)

    integer , allocatable :: order(:), pos(:), nbrs(:)
    real(dp), allocatable :: contribution(:)
    integer :: nv, i, j, v, w

    nv = this % num_vertices
    allocate(order(nv))
    order = this % dependency_order(FORWARD)

    allocate(pos(nv))
    do i = 1, nv
       pos(order(i)) = i
    end do

    allocate(contribution(size(seed)))

    adjoint = 0.0_dp
    adjoint(:, at) = seed

    ! against the dependency order, from just before the objective
    do i = pos(at) - 1, 1, -1
       v    = order(i)
       nbrs = this % out_neighbours(v)
       do j = 1, size(nbrs)
          w = nbrs(j)
          call edge_apply(v, w, adjoint(:, w), contribution)
          adjoint(:, v) = adjoint(:, v) + contribution
       end do
    end do

  end subroutine accumulate_adjoint

  !===================================================================!
  ! The walk's witness, a contract method (nopass): two self-contained
  ! fixtures, each judged against its own named tolerance. The
  ! canonical recurrence chain (x_{k+1} = a x_k, objective x_n^2 / 2)
  ! is checked against the analytic derivative at machine tolerance;
  ! the diamond dag (1->2, 1->3, 2->4, 3->4) - what proves the walk
  ! beyond chains - is checked against a central-difference nudge.
  ! Returns the worst defect normalized by its fixture's tolerance:
  ! a value below one certifies both claims.
  !===================================================================!

  impure real(dp) function verify_adjoint_accumulation() result(worst)

    type(witness_digraph) :: g
    real(dp), allocatable :: adjoint(:,:)
    real(dp) :: x(witness_chain_length), analytic, nudged, jp, jm, h
    real(dp) :: defect_chain, defect_diamond
    integer  :: k, n

    ! ---- fixture 1: the recurrence chain, analytic at machine tolerance
    n = witness_chain_length
    g = make_witness_digraph(n, [(k, k = 1, n-1)], [(k+1, k = 1, n-1)])

    x(1) = witness_start_value
    do k = 1, n-1
       x(k+1) = witness_recurrence_factor*x(k)
    end do

    allocate(adjoint(1, n))
    call g % accumulate_adjoint([x(n)], chain_edge_apply, n, adjoint)

    analytic    = witness_recurrence_factor**(2*(n-1)) * witness_start_value
    defect_chain = abs(adjoint(1,1) - analytic)/max(abs(analytic), 1.0_dp)

    deallocate(adjoint)

    ! ---- fixture 2: the diamond dag, central-difference nudge
    g = make_witness_digraph(4, [1, 1, 2, 3], [2, 3, 4, 4])

    allocate(adjoint(1, 4))
    call g % accumulate_adjoint([diamond_terminal(witness_start_value)], &
         & diamond_edge_apply, 4, adjoint)

    h  = 1.0d-6
    jp = 0.5_dp*diamond_terminal(witness_start_value + h)**2
    jm = 0.5_dp*diamond_terminal(witness_start_value - h)**2
    nudged = (jp - jm)/(2.0_dp*h)

    defect_diamond = abs(adjoint(1,1) - nudged)/max(abs(nudged), 1.0_dp)

    worst = max(defect_chain/witness_chain_tolerance, &
         &      defect_diamond/witness_diamond_tolerance)

  contains

    ! every chain edge carries the recurrence factor
    subroutine chain_edge_apply(tail, head, adjoint_head, contribution)
      integer , intent(in)  :: tail, head
      real(dp), intent(in)  :: adjoint_head(:)
      real(dp), intent(out) :: contribution(:)
      contribution = witness_recurrence_factor*adjoint_head
    end subroutine chain_edge_apply

    ! the diamond's per-edge factors
    subroutine diamond_edge_apply(tail, head, adjoint_head, contribution)
      integer , intent(in)  :: tail, head
      real(dp), intent(in)  :: adjoint_head(:)
      real(dp), intent(out) :: contribution(:)
      contribution = diamond_factor(tail, head)*adjoint_head
    end subroutine diamond_edge_apply

    pure real(dp) function diamond_factor(tail, head)
      integer, intent(in) :: tail, head
      if (tail .eq. 1 .and. head .eq. 2) then
         diamond_factor = witness_factor_12
      else if (tail .eq. 1 .and. head .eq. 3) then
         diamond_factor = witness_factor_13
      else if (tail .eq. 2 .and. head .eq. 4) then
         diamond_factor = witness_factor_24
      else
         diamond_factor = witness_factor_34
      end if
    end function diamond_factor

    ! the diamond's terminal value x4 as a function of the start value
    pure real(dp) function diamond_terminal(start)
      real(dp), intent(in) :: start
      diamond_terminal = witness_factor_24*(witness_factor_12*start) &
           &           + witness_factor_34*(witness_factor_13*start)
    end function diamond_terminal

  end function verify_adjoint_accumulation

  !===================================================================!
  ! Witness apparatus: stored directed queries and the fixture builder.
  !===================================================================!

  pure function witness_out_neighbours(this, v) result(nbrs)
    class(witness_digraph), intent(in) :: this
    integer               , intent(in) :: v
    integer, allocatable :: nbrs(:)
    nbrs = this % stored_out_neighbours(v)
  end function witness_out_neighbours

  pure function witness_in_neighbours(this, v) result(nbrs)
    class(witness_digraph), intent(in) :: this
    integer               , intent(in) :: v
    integer, allocatable :: nbrs(:)
    nbrs = this % stored_in_neighbours(v)
  end function witness_in_neighbours

  pure function make_witness_digraph(nv, tails, heads) result(g)

    integer, intent(in) :: nv
    integer, intent(in) :: tails(:), heads(:)
    type(witness_digraph) :: g

    integer :: i

    g % num_vertices = nv
    g % num_edges    = size(tails)

    allocate(g % vertices(nv))
    do i = 1, nv
       g % vertices(i) % number = i
       g % vertices(i) % part   = 1
    end do

    allocate(g % edges(g % num_edges))
    do i = 1, g % num_edges
       g % edges(i) % tail = tails(i)
       g % edges(i) % head = heads(i)
    end do

    call g % build_directed_adjacency()

  end function make_witness_digraph

  !===================================================================!
  ! RCB recursion: assign vertices idx(lo:hi) to k parts numbered
  ! base..base+k-1. Halve the part count, send a proportional share of
  ! the vertices (by count) to each side, split at the median along the
  ! axis of largest spread.
  !===================================================================!

  pure recursive subroutine rcb(idx, lo, hi, k, base, coords, part)

    integer , intent(inout) :: idx(:)
    integer , intent(in)    :: lo, hi, k, base
    real(dp), intent(in)    :: coords(:,:)
    integer , intent(inout) :: part(:)

    integer  :: n, d, axis, kL, kR, nL, mid
    real(dp) :: spread, best

    n = hi - lo + 1
    if (k .le. 1 .or. n .le. 0) then
       if (n .gt. 0) part(idx(lo:hi)) = base
       return
    end if

    ! axis of largest coordinate spread over this subset (zero-spread axes,
    ! e.g. z in 2d, are never chosen)
    axis = 1; best = -1.0_dp
    do d = 1, size(coords, 1)
       spread = maxval(coords(d, idx(lo:hi))) - minval(coords(d, idx(lo:hi)))
       if (spread .gt. best) then
          best = spread; axis = d
       end if
    end do

    ! split the part count, and the vertices proportionally (both sides nonempty)
    kL = k/2; kR = k - kL
    nL = nint(real(n,dp)*real(kL,dp)/real(k,dp))
    nL = max(1, min(n-1, nL))

    ! order idx(lo:hi) by the chosen axis so the nL smallest go left
    call qsort_axis(idx, lo, hi, coords, axis)
    mid = lo + nL - 1

    call rcb(idx, lo,    mid, kL, base,      coords, part)
    call rcb(idx, mid+1, hi,  kR, base + kL, coords, part)

  end subroutine rcb

  !===================================================================!
  ! Quicksort idx(lo:hi) ascending by coords(axis, idx(:)). Middle pivot so
  ! structured (already-ordered) meshes don't hit the O(n^2) worst case.
  !===================================================================!

  pure recursive subroutine qsort_axis(idx, lo, hi, coords, axis)

    integer , intent(inout) :: idx(:)
    integer , intent(in)    :: lo, hi, axis
    real(dp), intent(in)    :: coords(:,:)

    integer  :: i, j, tmp
    real(dp) :: pivot

    if (lo .ge. hi) return
    pivot = coords(axis, idx((lo+hi)/2))
    i = lo; j = hi
    do
       do while (coords(axis, idx(i)) .lt. pivot); i = i + 1; end do
       do while (coords(axis, idx(j)) .gt. pivot); j = j - 1; end do
       if (i .le. j) then
          tmp = idx(i); idx(i) = idx(j); idx(j) = tmp
          i = i + 1; j = j - 1
       end if
       if (i .gt. j) exit
    end do
    if (lo .lt. j) call qsort_axis(idx, lo, j, coords, axis)
    if (i .lt. hi) call qsort_axis(idx, i, hi, coords, axis)

  end subroutine qsort_axis

end module interface_graph
