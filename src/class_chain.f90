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
  ! The concrete chain graph.
  !===================================================================!

  type, extends(digraph) :: chain

     integer :: power = 1   ! An edge m -> k exists whenever k - m <= power.

   contains

     ! The directed contract is answered by the rule; nothing is stored.
     procedure :: out_neighbours
     procedure :: in_neighbours

     ! These are cheap rule overrides of the provided union queries.
     procedure :: neighbours
     procedure :: degree

     ! Form sets THIS object up as a chain of n links, in place - the
     ! seat a type that IS a chain (a marcher, an integrator) uses to
     ! form its own structure once its length is known.
     procedure :: form

  end type chain

  interface chain
     module procedure create
  end interface chain

contains

  !===================================================================!
  ! Create a chain of n vertices at the given power (default 1, the
  ! plain chain): the edges follow the rule, and none are stored.
  !===================================================================!

  pure type(chain) function create(n, num_variables, power) result(this)

    integer, intent(in)           :: n
    integer, intent(in), optional :: num_variables
    integer, intent(in), optional :: power

    call this % form(n, num_variables, power)

  end function create

  !===================================================================!
  ! Form: become a chain of n links, in place.
  !
  !    before:  ( me )              an object with no structure yet
  !    after :  ( me ) = 1 --> 2 --> ... --> n     at the given power
  !
  ! The constructor above builds a fresh chain through this; a type
  ! that extends chain calls it on itself the moment it learns how
  ! long its own sequence is.
  !===================================================================!

  pure subroutine form(this, n, num_variables, power)

    class(chain), intent(inout)        :: this
    integer     , intent(in)           :: n
    integer     , intent(in), optional :: num_variables
    integer     , intent(in), optional :: power

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

    !----------------------------------------------------------------!
    ! The vertex labels carry part stamps for the inherited
    ! partitioners. Form is re-entrant: forming again replaces the
    ! old labels.
    !----------------------------------------------------------------!

    if (allocated(this % vertices)) deallocate(this % vertices)
    allocate(this % vertices(n))
    do i = 1, n
       this % vertices(i) % number = i
       this % vertices(i) % part   = 1
    end do

  end subroutine form

  !===================================================================!
  ! The directed rule, forward: out-edges reach forward to the next
  ! power vertices, within 1..n.
  !===================================================================!

  pure function out_neighbours(this, v) result(nbrs)

    class(chain), intent(in) :: this
    integer     , intent(in) :: v

    integer, allocatable :: nbrs(:)
    integer              :: k

    nbrs = [(k, k = v + 1, min(this % num_vertices, v + this % power))]

  end function out_neighbours

  !===================================================================!
  ! The directed rule, backward: in-edges reach back to the previous
  ! power vertices, within 1..n.
  !===================================================================!

  pure function in_neighbours(this, v) result(nbrs)

    class(chain), intent(in) :: this
    integer     , intent(in) :: v

    integer, allocatable :: nbrs(:)
    integer              :: k

    nbrs = [(k, k = max(1, v - this % power), v - 1)]

  end function in_neighbours

  !===================================================================!
  ! Return the neighbours of vertex v by rule: everything within
  ! power of v.
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
  ! Return the degree by rule: how far the power reaches back plus
  ! how far it reaches forward.
  !===================================================================!

  pure integer function degree(this, v)

    class(chain), intent(in) :: this
    integer     , intent(in) :: v

    degree = min(v - 1, this % power) &
         & + min(this % num_vertices - v, this % power)

  end function degree

end module class_chain
