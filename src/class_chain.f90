!=====================================================================!
! The chain subclass of the abstract graph: vertices 1..n connected by
! the neighbour rule i -> i+1. The iterate sequence of a solver and
! the step sequence of a time integrator are both instances of this
! class.
!
! The adjacency is rule-generated, never materialized: neighbours and
! degree are answered by arithmetic, no edge list and no compressed
! adjacency are stored. Everything inherited from the ancestor
! (traversal orders, partitioning, the queries) consumes only the
! overridden neighbour queries, so it all operates on the rule
! directly - the forward traversal of a chain is 1..n, the reverse is
! n..1.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_chain

  use interface_graph, only : graph, vertex

  implicit none

  private
  public :: chain

  !===================================================================!
  ! Concrete chain graph
  !===================================================================!

  type, extends(graph) :: chain

   contains

     ! rule-generated adjacency: override the queries, store nothing
     procedure :: build_adjacency
     procedure :: neighbours
     procedure :: degree

  end type chain

  interface chain
     module procedure create
  end interface chain

contains

  !===================================================================!
  ! A chain of n vertices: n-1 edges by rule, no edge list stored.
  !===================================================================!

  pure type(chain) function create(n, num_variables) result(this)

    integer, intent(in)           :: n
    integer, intent(in), optional :: num_variables

    integer :: i

    this % num_variables = 1
    if (present(num_variables)) this % num_variables = num_variables

    this % num_vertices = n
    this % num_edges    = n - 1

    ! vertex labels (and part stamps for the inherited partitioners)
    allocate(this % vertices(n))
    do i = 1, n
       this % vertices(i) % number = i
       this % vertices(i) % part   = 1
    end do

  end function create

  !===================================================================!
  ! Nothing to materialize: the adjacency is the rule i -> i+1. The
  ! edge count restates the rule.
  !===================================================================!

  pure subroutine build_adjacency(this)

    class(chain), intent(inout) :: this

    this % num_edges = max(this % num_vertices - 1, 0)

  end subroutine build_adjacency

  !===================================================================!
  ! Neighbours of vertex i by rule: i-1 and i+1, within 1..n.
  !===================================================================!

  pure function neighbours(this, v) result(nbrs)

    class(chain), intent(in) :: this
    integer     , intent(in) :: v

    integer, allocatable :: nbrs(:)

    if (v .gt. 1 .and. v .lt. this % num_vertices) then
       nbrs = [v-1, v+1]
    else if (this % num_vertices .le. 1) then
       allocate(nbrs(0))
    else if (v .eq. 1) then
       nbrs = [2]
    else
       nbrs = [v-1]
    end if

  end function neighbours

  !===================================================================!
  ! Degree by rule: interior vertices 2, end vertices 1.
  !===================================================================!

  pure integer function degree(this, v)

    class(chain), intent(in) :: this
    integer     , intent(in) :: v

    if (this % num_vertices .le. 1) then
       degree = 0
    else if (v .eq. 1 .or. v .eq. this % num_vertices) then
       degree = 1
    else
       degree = 2
    end if

  end function degree

end module class_chain
