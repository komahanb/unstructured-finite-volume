!=====================================================================!
! The chain subclass of the directed graph: vertices 1..n wired by the
! rule i -> i+1. The iterate sequence of a solver and the step sequence
! of a time integrator are both instances of this class; its dependency
! order is 1..n by construction, and the discrete adjoint traverses it
! in reverse through the inherited accumulation.
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

  use interface_graph, only : digraph

  implicit none

  private
  public :: chain

  !===================================================================!
  ! Concrete chain graph
  !===================================================================!

  type, extends(digraph) :: chain

   contains

     ! the directed contract, answered by the rule - nothing stored
     procedure :: out_neighbours
     procedure :: in_neighbours

     ! cheap rule overrides of the provided union queries
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
  ! The directed rule: the one out-edge goes to i+1, the one in-edge
  ! comes from i-1, within 1..n.
  !===================================================================!

  pure function out_neighbours(this, v) result(nbrs)

    class(chain), intent(in) :: this
    integer     , intent(in) :: v

    integer, allocatable :: nbrs(:)

    if (v .lt. this % num_vertices) then
       nbrs = [v+1]
    else
       allocate(nbrs(0))
    end if

  end function out_neighbours

  pure function in_neighbours(this, v) result(nbrs)

    class(chain), intent(in) :: this
    integer     , intent(in) :: v

    integer, allocatable :: nbrs(:)

    if (v .gt. 1) then
       nbrs = [v-1]
    else
       allocate(nbrs(0))
    end if

  end function in_neighbours

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
