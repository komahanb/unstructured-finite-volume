!=====================================================================!
! Concrete graph supports.
!
! A support is a chosen set of indices and nothing more. Each of
! these two types holds its set as a plain integer array, which is
! exactly what a set of indices is.
!
!      all_vertices              tagged_edges('wall')
!      +-------------------+     +--------------+
!      | 1  2  3  4  5  6  |     | 11  14  19   |
!      +-------------------+     +--------------+
!         a vertex_support          an edge_support
!
! The two are separate types rather than one type holding a flag.
! This separation distinguishes a vertex field from an edge field
! everywhere else in the library, and the compiler enforces it.
!
! A support is requested once, when an operation begins, and then
! looped over. It is never built inside a loop over cells or faces.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_graph_support

  use abstract_graph_types, only : graph_vertex_support, graph_edge_support
  use abstract_graph_types, only : GRAPH_SUPPORT_VERTEX, GRAPH_SUPPORT_EDGE

  implicit none

  private
  public :: vertex_support, edge_support

  !===================================================================!
  ! A chosen set of vertices.
  !===================================================================!

  type, extends(graph_vertex_support) :: vertex_support

     integer, allocatable :: indices(:)

   contains

     !----------------------------------------------------------------!
     ! The contract, and nothing beyond it.
     !----------------------------------------------------------------!

     procedure :: kind       => vertex_support_kind
     procedure :: size       => vertex_support_size
     procedure :: vertex_indices => vertex_support_indices

  end type vertex_support

  !===================================================================!
  ! A chosen set of edges.
  !===================================================================!

  type, extends(graph_edge_support) :: edge_support

     integer, allocatable :: indices(:)

   contains

     !----------------------------------------------------------------!
     ! The contract, and nothing beyond it.
     !----------------------------------------------------------------!

     procedure :: kind     => edge_support_kind
     procedure :: size     => edge_support_size
     procedure :: edge_indices => edge_support_indices

  end type edge_support

  !===================================================================!
  ! Constructors. Hand in the indices; the support keeps a copy.
  !===================================================================!

  interface vertex_support
     module procedure create_vertex_support
  end interface vertex_support

  interface edge_support
     module procedure create_edge_support
  end interface edge_support

contains

  !===================================================================!
  ! Build a vertex support from a list of vertex numbers.
  !===================================================================!

  pure type(vertex_support) function create_vertex_support(indices) result(this)

    integer, intent(in) :: indices(:)

    allocate(this % indices, source=indices)

  end function create_vertex_support

  !===================================================================!
  ! Build an edge support from a list of edge numbers.
  !===================================================================!

  pure type(edge_support) function create_edge_support(indices) result(this)

    integer, intent(in) :: indices(:)

    allocate(this % indices, source=indices)

  end function create_edge_support

  !===================================================================!
  ! A vertex support holds vertex indices. It says so.
  !===================================================================!

  pure integer function vertex_support_kind(this)

    class(vertex_support), intent(in) :: this

    vertex_support_kind = GRAPH_SUPPORT_VERTEX

  end function vertex_support_kind

  !===================================================================!
  ! How many vertices are in the set.
  !===================================================================!

  pure integer function vertex_support_size(this)

    class(vertex_support), intent(in) :: this

    if (allocated(this % indices)) then
       vertex_support_size = size(this % indices)
    else
       vertex_support_size = 0
    end if

  end function vertex_support_size

  !===================================================================!
  ! Hand over the vertex numbers. An empty set returns a
  ! zero-length array, never with something unallocated, so a caller
  ! can loop without asking first.
  !===================================================================!

  pure subroutine vertex_support_indices(this, indices)

    class(vertex_support), intent(in)  :: this
    integer, allocatable, intent(out)  :: indices(:)

    if (allocated(this % indices)) then
       indices = this % indices
    else
       allocate(indices(0))
    end if

  end subroutine vertex_support_indices

  !===================================================================!
  ! An edge support holds edge indices. It says so.
  !===================================================================!

  pure integer function edge_support_kind(this)

    class(edge_support), intent(in) :: this

    edge_support_kind = GRAPH_SUPPORT_EDGE

  end function edge_support_kind

  !===================================================================!
  ! How many edges are in the set.
  !===================================================================!

  pure integer function edge_support_size(this)

    class(edge_support), intent(in) :: this

    if (allocated(this % indices)) then
       edge_support_size = size(this % indices)
    else
       edge_support_size = 0
    end if

  end function edge_support_size

  !===================================================================!
  ! Hand over the edge numbers, with the same empty-set promise.
  !===================================================================!

  pure subroutine edge_support_indices(this, indices)

    class(edge_support), intent(in)   :: this
    integer, allocatable, intent(out) :: indices(:)

    if (allocated(this % indices)) then
       indices = this % indices
    else
       allocate(indices(0))
    end if

  end subroutine edge_support_indices

end module class_graph_support
