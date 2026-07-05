!=====================================================================!
! The mesh tenant of the abstract graph: cells are vertices, a shared
! interior face is an edge. This class knows the mesh (its constructor
! takes one, and each edge remembers the interior face that made it);
! the ancestor knows only graph-shaped data - adjacency, traversals,
! and the partition contract all live there.
!
! It owns the degree-of-freedom map dof(cell, variable) -> global dof
! through the ancestor: with num_variables fields per cell the dofs are
! interleaved variable-fastest, so for a single field dof(cell,1) = cell.
!
! A partition of a graph is still graph-shaped data, so partition /
! partition_rcb (inherited) stamp this same object in place and the
! per-part owned and ghost (halo) cell lists are queried from it -
! exactly the bookkeeping a distributed solver needs. Pure serial
! integer bookkeeping, so every image computes the identical partition
! from the replicated graph.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_graph

  use iso_fortran_env, only : dp => REAL64
  use interface_graph, only : graph, vertex, edge
  use class_mesh     , only : mesh

  implicit none

  private
  public :: graph, vertex, edge   ! re-export the ancestor's vocabulary
  public :: mesh_graph

  !===================================================================!
  ! Concrete mesh graph
  !===================================================================!

  type, extends(graph) :: mesh_graph

     ! provenance: edge_face(e) is the interior face that made edge e
     integer, allocatable :: edge_face(:)

  end type mesh_graph

  interface mesh_graph
     module procedure create
  end interface mesh_graph

contains

  !===================================================================!
  ! Build the graph from a mesh - cells become vertices, interior faces
  ! (num_face_cells == 2) become edges. The retained adjacency is built
  ! here, once, and queried ever after.
  !===================================================================!

  pure type(mesh_graph) function create(grid, num_variables) result(this)

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
    allocate(this % edge_face(this % num_edges))
    ctr = 0
    do iface = 1, grid % num_faces
       if (grid % num_face_cells(iface) .eq. 2) then
          ctr = ctr + 1
          this % edges(ctr) % tail = grid % face_cells(1, iface)
          this % edges(ctr) % head = grid % face_cells(2, iface)
          this % edge_face(ctr)    = iface
       end if
    end do

    call this % build_adjacency()

  end function create

end module class_graph
