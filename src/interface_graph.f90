!=====================================================================!
! Abstract graph: vertices, edges, the retained adjacency, traversal
! orders, and the partition contract. The computation itself is a graph,
! so this ancestor knows only graph-shaped data - never a mesh, a face,
! or any payload. Dependencies point downward only: a tenant knows its
! own world (class_graph knows the mesh, class_chain knows its rule);
! the ancestor knows nothing about either.
!
! The adjacency is a compressed neighbour list (xadj/adj), retained and
! publicly queryable - previously buried as scratch inside the bfs
! partitioner and rebuilt on every use. Every piece of machinery here
! consumes only the neighbour queries, so a tenant may either store the
! adjacency (build_adjacency from an edge list) or generate it by rule
! (the chain overrides neighbours/degree and stores nothing).
!
! Traversal: traversal_order(mode) visits vertices breadth-first in
! dependency order (FORWARD), or the same order flipped (REVERSE) - the
! direction tags come from module_solve_mode, the one tag vocabulary.
!
! Partition: stamped in place (partition / partition_rcb), then queried
! (part_of, owned, ghosts, balance, edge_cut). A partition of a graph is
! still graph-shaped data.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module interface_graph

  use iso_fortran_env  , only : dp => REAL64
  use module_solve_mode, only : FORWARD, REVERSE

  implicit none

  private
  public :: graph, vertex, edge

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

     ! dof indexing (variable-fastest over vertices)
     procedure :: dof
     procedure :: num_dofs

     ! adjacency: build once from the edge list, query many times.
     ! A rule-generated tenant overrides the queries and builds nothing.
     procedure :: build_adjacency
     procedure :: neighbours
     procedure :: degree

     ! traversal order over the whole graph (FORWARD/REVERSE)
     procedure :: traversal_order

     ! partitioners - stamp vertex % part and gather the csr, in place
     procedure :: partition
     procedure :: partition_rcb

     ! partition queries
     procedure :: part_of
     procedure :: owned
     procedure :: ghosts
     procedure :: n_owned
     procedure :: n_ghosts
     procedure :: balance
     procedure :: edge_cut

     procedure :: print
     procedure :: print_partition

     ! internal: stamp vertex % part, then gather the owned/ghost csr
     procedure, private :: stamp_bfs
     procedure, private :: stamp_rcb
     procedure, private :: gather_partition

  end type graph

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
  ! Build the retained adjacency from the undirected edge list: per-
  ! vertex degrees accumulate into the row pointer xadj, then the
  ! neighbours scatter into adj in edge order (both directions).
  !===================================================================!

  pure subroutine build_adjacency(this)

    class(graph), intent(inout) :: this

    integer, allocatable :: ptr(:)
    integer :: nv, e, i

    nv = this % num_vertices

    if (allocated(this % xadj)) deallocate(this % xadj)
    if (allocated(this % adj))  deallocate(this % adj)

    allocate(this % xadj(nv+1))
    this % xadj = 0

    do e = 1, this % num_edges
       this % xadj(this % edges(e) % tail + 1) = this % xadj(this % edges(e) % tail + 1) + 1
       this % xadj(this % edges(e) % head + 1) = this % xadj(this % edges(e) % head + 1) + 1
    end do

    this % xadj(1) = 1
    do i = 1, nv
       this % xadj(i+1) = this % xadj(i+1) + this % xadj(i)
    end do

    allocate(this % adj(this % xadj(nv+1)-1))
    allocate(ptr(nv))
    ptr = this % xadj(1:nv)

    do e = 1, this % num_edges

       associate(t => this % edges(e) % tail, h => this % edges(e) % head)

       this % adj(ptr(t)) = h
       ptr(t)             = ptr(t) + 1

       this % adj(ptr(h)) = t
       ptr(h)             = ptr(h) + 1

       end associate

    end do

  end subroutine build_adjacency

  !===================================================================!
  ! Neighbours of vertex v, from the retained adjacency. Precondition:
  ! the adjacency exists (a stored tenant builds it at creation; a
  ! rule-generated tenant overrides this query instead).
  !===================================================================!

  pure function neighbours(this, v) result(nbrs)

    class(graph), intent(in) :: this
    integer     , intent(in) :: v

    integer, allocatable :: nbrs(:)

    if (.not. allocated(this % xadj)) then
       error stop "graph: adjacency not built - build_adjacency or override neighbours"
    end if

    nbrs = this % adj(this % xadj(v) : this % xadj(v+1)-1)

  end function neighbours

  !===================================================================!
  ! Degree of vertex v (the number of its neighbours)
  !===================================================================!

  pure integer function degree(this, v)

    class(graph), intent(in) :: this
    integer     , intent(in) :: v

    if (.not. allocated(this % xadj)) then
       error stop "graph: adjacency not built - build_adjacency or override degree"
    end if

    degree = this % xadj(v+1) - this % xadj(v)

  end function degree

  !===================================================================!
  ! Traversal order over all vertices: breadth-first in dependency
  ! order, restarting for disconnected components (FORWARD); the same
  ! order flipped (REVERSE). Consumes only the neighbour queries.
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
  ! Partition into nparts pieces by breadth-first ordering, in place
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

    integer, allocatable :: ptr(:), cnt(:), mark(:), nbrs(:)
    integer :: nv, v, w, i, k, kv, pos

    nv = this % num_vertices
    this % nparts = nparts

    ! ---- owned vertices per part (counting sort by part) ----
    if (allocated(this % own_ptr))  deallocate(this % own_ptr)
    if (allocated(this % own_list)) deallocate(this % own_list)
    allocate(this % own_ptr(nparts+1)); this % own_ptr = 0
    do v = 1, nv
       this % own_ptr(this % vertices(v) % part + 1) = &
            & this % own_ptr(this % vertices(v) % part + 1) + 1
    end do
    this % own_ptr(1) = 1
    do k = 1, nparts
       this % own_ptr(k+1) = this % own_ptr(k+1) + this % own_ptr(k)
    end do
    allocate(this % own_list(nv))
    allocate(ptr(nparts)); ptr = this % own_ptr(1:nparts)
    do v = 1, nv
       k = this % vertices(v) % part
       this % own_list(ptr(k)) = v
       ptr(k) = ptr(k) + 1
    end do
    deallocate(ptr)

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
  ! Vertices owned by part k (its matrix rows)
  !===================================================================!

  pure function owned(this, k) result(cells)

    class(graph), intent(in) :: this
    integer     , intent(in) :: k

    integer, allocatable :: cells(:)

    cells = this % own_list(this % own_ptr(k) : this % own_ptr(k+1)-1)

  end function owned

  !===================================================================!
  ! Halo (ghost) vertices part k needs from other parts
  !===================================================================!

  pure function ghosts(this, k) result(cells)

    class(graph), intent(in) :: this
    integer     , intent(in) :: k

    integer, allocatable :: cells(:)

    cells = this % gh_list(this % gh_ptr(k) : this % gh_ptr(k+1)-1)

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
