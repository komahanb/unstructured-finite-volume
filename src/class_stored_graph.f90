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
! The directed twin, stored_digraph, is the same citizen with arrows:
! an edge list read as tail -> head. Vertices may carry caller-given
! numbers, so a small directed graph can name things beyond itself -
! a multigrid cycle, for instance, is stations and moves: each vertex
! numbered by the level it visits, each arrow the next leg of the
! trip. Build the graph, hand it in, and the machinery follows it.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_stored_graph

  use interface_graph, only : graph, digraph

  implicit none

  private
  public :: stored_graph
  public :: stored_digraph

  type, extends(graph) :: stored_graph

   contains

     ! the deferred contract, delegated to the shared stored mechanism
     procedure :: neighbours
     procedure :: degree

  end type stored_graph

  type, extends(digraph) :: stored_digraph

   contains

     ! the directed contract, delegated to the stored mechanism
     procedure :: out_neighbours
     procedure :: in_neighbours

  end type stored_digraph

  interface stored_graph
     module procedure create
     module procedure create_quotient
     module procedure create_refined
  end interface stored_graph

  interface stored_digraph
     module procedure create_directed
  end interface stored_digraph

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
  ! The zoom: the refinement of any graph. Every vertex splits into
  ! the given number of children (refine_edges on the ancestor draws
  ! the refined edges), and the parent map is adopted as the new
  ! graph's partition - so the quotient of the refinement is the
  ! original graph, and the squint undoes the zoom.
  !===================================================================!

  pure type(stored_graph) function create_refined(coarse, children) result(this)

    class(graph), intent(in) :: coarse
    integer     , intent(in) :: children

    integer, allocatable :: tails(:), heads(:)
    integer              :: i

    call coarse % refine_edges(children, tails, heads)
    this = create(coarse % num_vertices*children, tails, heads, coarse % num_variables)

    ! the parent of child i is arithmetic: parts adopted as a partition
    call this % set_partition([((i-1)/children + 1, i = 1, this % num_vertices)])

  end function create_refined

  !===================================================================!
  ! The directed twin, from a tail -> head edge list. Vertices carry
  ! caller-given numbers when present (a schedule names levels with
  ! them); otherwise they number themselves 1..nv.
  !===================================================================!

  pure type(stored_digraph) function create_directed(nv, tails, heads, numbers) result(this)

    integer, intent(in)           :: nv
    integer, intent(in)           :: tails(:), heads(:)
    integer, intent(in), optional :: numbers(:)

    integer :: i

    this % num_vertices = nv
    allocate(this % vertices(nv))
    do i = 1, nv
       this % vertices(i) % number = i
       if (present(numbers)) this % vertices(i) % number = numbers(i)
       this % vertices(i) % part = 1
    end do

    this % num_edges = size(tails)
    allocate(this % edges(this % num_edges))
    do i = 1, this % num_edges
       this % edges(i) % tail = tails(i)
       this % edges(i) % head = heads(i)
    end do

    call this % build_directed_adjacency()

  end function create_directed

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

  pure function out_neighbours(this, v) result(nbrs)

    class(stored_digraph), intent(in) :: this
    integer              , intent(in) :: v

    integer, allocatable :: nbrs(:)

    nbrs = this % stored_out_neighbours(v)

  end function out_neighbours

  pure function in_neighbours(this, v) result(nbrs)

    class(stored_digraph), intent(in) :: this
    integer              , intent(in) :: v

    integer, allocatable :: nbrs(:)

    nbrs = this % stored_in_neighbours(v)

  end function in_neighbours

end module class_stored_graph
