!=====================================================================!
! The dof/connectivity graph. Cells are vertices, a shared interior
! face is an edge. It owns the degree-of-freedom map
!
!   dof(cell, variable) -> global dof
!
! and the adjacency (who couples to whom). With num_variables fields per
! cell the dofs are interleaved variable-fastest, so for a single field
! dof(cell,1) = cell and everything reduces to the scalar case.
!
! Lots of other uses live here later: dof reordering, coloring, and -
! the long game - partitioning with parmetis and coarray halos.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_graph

  use iso_fortran_env, only : dp => REAL64
  use class_mesh     , only : mesh

  implicit none

  private
  public :: graph, vertex, edge

  type :: vertex
     integer :: number      ! cell number
     integer :: part = 1    ! owning image (partitioning, later)
  end type vertex

  type :: edge
     integer :: tail, head  ! coupled vertices (cells sharing a face)
     integer :: face        ! the interior face that made this edge
  end type edge

  type :: graph

     integer :: num_vertices
     integer :: num_edges
     integer :: num_variables = 1
     type(vertex), allocatable :: vertices(:)
     type(edge)  , allocatable :: edges(:)

   contains

     procedure :: dof
     procedure :: num_dofs
     procedure :: partition
     procedure :: partition_rcb
     procedure :: edge_cut
     procedure :: print

  end type graph

  interface graph
     module procedure create
  end interface graph

contains

  !===================================================================!
  ! Build the graph from a mesh - cells become vertices, interior faces
  ! (num_face_cells == 2) become edges.
  !===================================================================!

  type(graph) function create(grid, num_variables) result(this)

    type(mesh), intent(in)           :: grid
    integer   , intent(in), optional :: num_variables

    integer :: icell, iface, ctr

    this % num_variables = 1
    if (present(num_variables)) this % num_variables = num_variables

    ! Vertices - one per cell
    this % num_vertices = grid % num_cells
    allocate(this % vertices(this % num_vertices))
    do icell = 1, this % num_vertices
       this % vertices(icell) % number = grid % cell_numbers(icell)
       this % vertices(icell) % part   = 1
    end do

    ! Edges - one per interior face
    this % num_edges = count(grid % num_face_cells .eq. 2)
    allocate(this % edges(this % num_edges))
    ctr = 0
    do iface = 1, grid % num_faces
       if (grid % num_face_cells(iface) .eq. 2) then
          ctr = ctr + 1
          this % edges(ctr) % tail = grid % face_cells(1, iface)
          this % edges(ctr) % head = grid % face_cells(2, iface)
          this % edges(ctr) % face = iface
       end if
    end do

  end function create

  !===================================================================!
  ! Global dof for a (cell, variable) pair. Variable-fastest so a single
  ! field maps dof(cell,1) = cell.
  !===================================================================!

  pure type(integer) function dof(this, cell, ivar)

    class(graph), intent(in) :: this
    integer     , intent(in) :: cell, ivar

    dof = (cell - 1)*this % num_variables + ivar

  end function dof

  !===================================================================!
  ! Total number of degrees of freedom in the graph
  !===================================================================!

  pure type(integer) function num_dofs(this)

    class(graph), intent(in) :: this

    num_dofs = this % num_vertices * this % num_variables

  end function num_dofs

  !===================================================================!
  ! Partition the vertices into nparts pieces and stamp vertex % part.
  !
  ! This is the serial placeholder: a breadth-first ordering of the
  ! graph cut into nparts balanced contiguous chunks. bfs keeps
  ! neighbours close so the chunks stay mostly connected and the cut is
  ! reasonable. (KB: swap this body for a parmetis call - same in/out,
  ! the graph already holds the adjacency parmetis wants.)
  !===================================================================!

  subroutine partition(this, nparts)

    class(graph) , intent(inout) :: this
    integer      , intent(in)    :: nparts

    integer, allocatable :: xadj(:), adj(:), ptr(:), order(:), queue(:)
    logical, allocatable :: seen(:)
    integer :: nv, e, i, v, w, pos, qh, qt, start

    nv = this % num_vertices

    ! Compressed adjacency (csr) from the undirected edge list. First the
    ! per-vertex degree, accumulated into the row pointer xadj ...
    allocate(xadj(nv+1))
    xadj = 0

    do e = 1, this % num_edges
       xadj(this % edges(e) % tail + 1) = xadj(this % edges(e) % tail + 1) + 1
       xadj(this % edges(e) % head + 1) = xadj(this % edges(e) % head + 1) + 1
    end do

    xadj(1) = 1
    do i = 1, nv
       xadj(i+1) = xadj(i+1) + xadj(i)
    end do

    ! ... then scatter the neighbours into adj using a running pointer
    allocate(adj(xadj(nv+1)-1))
    allocate(ptr(nv))
    ptr = xadj(1:nv)

    do e = 1, this % num_edges

       associate(t => this % edges(e) % tail, h => this % edges(e) % head)

       adj(ptr(t)) = h
       ptr(t)      = ptr(t) + 1

       adj(ptr(h)) = t
       ptr(h)      = ptr(h) + 1

       end associate

    end do

    ! Breadth-first ordering, restarting for disconnected components
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

          do i = xadj(v), xadj(v+1)-1
             w = adj(i)
             if (.not. seen(w)) then
                seen(w) = .true.
                qt = qt + 1
                queue(qt) = w
             end if
          end do

       end do

    end do

    ! Cut the bfs order into nparts balanced contiguous chunks
    do pos = 1, nv
       this % vertices(order(pos)) % part = (pos-1)*nparts/nv + 1
    end do

  end subroutine partition

  !===================================================================!
  ! Geometric partitioner: recursive coordinate bisection (RCB).
  !
  ! Sibling to partition() above (which is adjacency-only BFS). Uses the
  ! cell centroids coords(:,cell) (e.g. mesh % cell_centers) to recursively
  ! split the cells by the median coordinate along the longest-spread axis -
  ! balanced, COMPACT subdomains (small edge cut / halo) where the BFS chunks
  ! would go stringy. pure fortran, no external deps; stamps the same
  ! vertex % part contract so everything downstream is unchanged.
  ! (KB: METIS could still slot in behind this same contract.)
  !===================================================================!

  subroutine partition_rcb(this, coords, nparts)

    class(graph), intent(inout) :: this
    real(dp)    , intent(in)    :: coords(:,:)   ! (ndim, num_vertices) centroids
    integer     , intent(in)    :: nparts

    integer, allocatable :: idx(:), part(:)
    integer :: nv, c

    nv = this % num_vertices
    if (nparts .lt. 1)  error stop "partition_rcb: nparts < 1"
    if (nparts .gt. nv) error stop "partition_rcb: more parts than cells"

    allocate(idx(nv), part(nv))
    do c = 1, nv
       idx(c) = c
    end do
    part = 1

    call rcb(idx, 1, nv, nparts, 1, coords, part)

    do c = 1, nv
       this % vertices(c) % part = part(c)
    end do

  end subroutine partition_rcb

  !===================================================================!
  ! RCB recursion: assign cells idx(lo:hi) to k parts numbered
  ! base..base+k-1. Halve the part count, send a proportional share of the
  ! cells (by count - CG work is ~one row per cell) to each side, split at
  ! the median along the axis of largest spread.
  !===================================================================!

  recursive subroutine rcb(idx, lo, hi, k, base, coords, part)

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

    ! split the part count, and the cells proportionally (both sides nonempty)
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

  recursive subroutine qsort_axis(idx, lo, hi, coords, axis)

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

  !===================================================================!
  ! Number of edges whose endpoints lie in different partitions.
  !===================================================================!

  pure type(integer) function edge_cut(this)

    class(graph), intent(in) :: this

    integer :: e

    edge_cut = 0

    do e = 1, this % num_edges
       if (this % vertices(this % edges(e) % tail) % part .ne. &
            & this % vertices(this % edges(e) % head) % part) then
          edge_cut = edge_cut + 1
       end if
    end do

  end function edge_cut

  !===================================================================!
  ! Print a one line summary of the graph
  !===================================================================!

  subroutine print(this)

    class(graph), intent(in) :: this

    write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') &
         & "graph: ", this % num_vertices, " vertices, ", this % num_edges, &
         & " edges, ", this % num_variables, " variables/cell, ", this % num_dofs(), " dofs"

  end subroutine print

end module class_graph
