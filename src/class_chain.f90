!=====================================================================!
! The chain subclass of the directed graph: vertices 1..n wired by
! the rule i -> i+1, optionally raised to a power - the formal graph
! power, an edge m -> k whenever k - m <= power:
!
!      .-----------.-----------.
!      |           v           v
!      1 --> 2 --> 3 --> 4 --> 5        chain(5, power = 2)
!            |           ^
!            '-----------'
!
! The iterate sequence of a solver and the step sequence of a time
! integrator are both chains; a time stencil that reads deeper than
! one step is the same chain at a higher power (class_bdf carries
! one). Dependency order is 1..n by construction, and the discrete
! adjoint traverses the chain in reverse through the inherited
! accumulation.
!
! The adjacency is rule-generated, never materialized: neighbours and
! degree are answered by arithmetic, no edge list and no compressed
! adjacency are stored. Everything inherited from the ancestor
! (traversal orders, partitioning, the queries) consumes only the
! overridden neighbour queries, so it all operates on the rule
! directly.
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

     integer :: power = 1   ! edge m -> k whenever k - m <= power

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
  ! A chain of n vertices at the given power (default 1, the plain
  ! chain): edges by rule, none stored.
  !===================================================================!

  pure type(chain) function create(n, num_variables, power) result(this)

    integer, intent(in)           :: n
    integer, intent(in), optional :: num_variables
    integer, intent(in), optional :: power

    integer :: i

    this % num_variables = 1
    if (present(num_variables)) this % num_variables = num_variables

    this % power = 1
    if (present(power)) this % power = max(1, power)

    this % num_vertices = n
    this % num_edges    = 0
    do i = 1, n
       this % num_edges = this % num_edges + max(0, min(this % power, n - i))
    end do

    ! vertex labels (and part stamps for the inherited partitioners)
    allocate(this % vertices(n))
    do i = 1, n
       this % vertices(i) % number = i
       this % vertices(i) % part   = 1
    end do

  end function create

  !===================================================================!
  ! The directed rule: out-edges reach forward to the next power
  ! vertices, in-edges reach back the same way, within 1..n.
  !===================================================================!

  pure function out_neighbours(this, v) result(nbrs)

    class(chain), intent(in) :: this
    integer     , intent(in) :: v

    integer, allocatable :: nbrs(:)
    integer              :: k

    nbrs = [(k, k = v + 1, min(this % num_vertices, v + this % power))]

  end function out_neighbours

  pure function in_neighbours(this, v) result(nbrs)

    class(chain), intent(in) :: this
    integer     , intent(in) :: v

    integer, allocatable :: nbrs(:)
    integer              :: k

    nbrs = [(k, k = max(1, v - this % power), v - 1)]

  end function in_neighbours

  !===================================================================!
  ! Neighbours of vertex v by rule: everything within power of v.
  !===================================================================!

  pure function neighbours(this, v) result(nbrs)

    class(chain), intent(in) :: this
    integer     , intent(in) :: v

    integer, allocatable :: nbrs(:)
    integer              :: k

    nbrs = [(k, k = max(1, v - this % power), v - 1), &
         &  (k, k = v + 1, min(this % num_vertices, v + this % power))]

  end function neighbours

  !===================================================================!
  ! Degree by rule: how far the power reaches back plus how far it
  ! reaches forward.
  !===================================================================!

  pure integer function degree(this, v)

    class(chain), intent(in) :: this
    integer     , intent(in) :: v

    degree = min(v - 1, this % power) &
         & + min(this % num_vertices - v, this % power)

  end function degree

end module class_chain
