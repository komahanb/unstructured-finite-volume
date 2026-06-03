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

    ! Compressed adjacency (csr) from the undirected edge list
    allocate(xadj(nv+1)); xadj = 0
    do e = 1, this % num_edges
       xadj(this % edges(e) % tail + 1) = xadj(this % edges(e) % tail + 1) + 1
       xadj(this % edges(e) % head + 1) = xadj(this % edges(e) % head + 1) + 1
    end do
    xadj(1) = 1
    do i = 1, nv
       xadj(i+1) = xadj(i+1) + xadj(i)
    end do
    allocate(adj(xadj(nv+1)-1))
    allocate(ptr(nv)); ptr = xadj(1:nv)
    do e = 1, this % num_edges
       associate(t => this % edges(e) % tail, h => this % edges(e) % head)
       adj(ptr(t)) = h; ptr(t) = ptr(t) + 1
       adj(ptr(h)) = t; ptr(h) = ptr(h) + 1
       end associate
    end do

    ! Breadth-first ordering (restart for disconnected components)
    allocate(order(nv), queue(nv), seen(nv)); seen = .false.
    pos = 0; qh = 1; qt = 0
    do start = 1, nv
       if (seen(start)) cycle
       qt = qt + 1; queue(qt) = start; seen(start) = .true.
       do while (qh .le. qt)
          v = queue(qh); qh = qh + 1
          pos = pos + 1; order(pos) = v
          do i = xadj(v), xadj(v+1)-1
             w = adj(i)
             if (.not. seen(w)) then
                seen(w) = .true.; qt = qt + 1; queue(qt) = w
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
