!=====================================================================!
! The plainest citizen of the graph world: a stored graph built from
! an edge list handed straight to the constructor - no mesh, no rule,
! no matrix behind it. Whoever can name vertices and edges gets the
! whole ancestry: adjacency, traversals, orbits, partitions.
!
! Its second constructor is the squint: give it any partitioned graph
! and it becomes the quotient - parts as vertices, touching parts as
! edges (quotient_edges on the ancestor does the reading). The
! quotient of a stored graph is again a stored graph, so the squint
! composes: a hierarchy is just repeated construction.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_stored_graph

  use interface_graph, only : graph, vertex, edge

  implicit none

  private
  public :: stored_graph

  type, extends(graph) :: stored_graph

   contains

     ! the deferred contract, delegated to the shared stored mechanism
     procedure :: neighbours
     procedure :: degree

  end type stored_graph

  interface stored_graph
     module procedure create
     module procedure create_quotient
  end interface stored_graph

contains

  !===================================================================!
  ! Build from an edge list: vertices numbered 1..nv, one edge per
  ! (tail, head) pair, the retained adjacency built once.
  !===================================================================!

  pure type(stored_graph) function create(nv, tails, heads, num_variables) result(this)

    integer, intent(in)           :: nv
    integer, intent(in)           :: tails(:), heads(:)
    integer, intent(in), optional :: num_variables

    integer :: i

    this % num_variables = 1
    if (present(num_variables)) this % num_variables = num_variables

    this % num_vertices = nv
    allocate(this % vertices(nv))
    do i = 1, nv
       this % vertices(i) % number = i
       this % vertices(i) % part   = 1
    end do

    this % num_edges = size(tails)
    allocate(this % edges(this % num_edges))
    do i = 1, this % num_edges
       this % edges(i) % tail = tails(i)
       this % edges(i) % head = heads(i)
    end do

    call this % build_adjacency()

  end function create

  !===================================================================!
  ! The squint: the quotient of a partitioned graph. Parts become
  ! vertices, touching parts become edges, and the result is the same
  ! kind of animal as its input - ready to be squinted again.
  !===================================================================!

  pure type(stored_graph) function create_quotient(fine) result(this)

    class(graph), intent(in) :: fine

    integer, allocatable :: tails(:), heads(:)

    call fine % quotient_edges(tails, heads)
    this = create(fine % nparts, tails, heads, fine % num_variables)

  end function create_quotient

  !===================================================================!
  ! The deferred neighbour queries: one-line delegations to the stored
  ! adjacency the constructor built.
  !===================================================================!

  pure function neighbours(this, v) result(nbrs)

    class(stored_graph), intent(in) :: this
    integer            , intent(in) :: v

    integer, allocatable :: nbrs(:)

    nbrs = this % stored_neighbours(v)

  end function neighbours

  pure integer function degree(this, v)

    class(stored_graph), intent(in) :: this
    integer            , intent(in) :: v

    degree = this % stored_degree(v)

  end function degree

end module class_stored_graph
