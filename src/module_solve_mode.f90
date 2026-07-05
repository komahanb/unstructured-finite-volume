!=====================================================================!
! The tag vocabulary of the computation graph, in one home.
!
! mode - the direction an edge is traversed: FORWARD drives the primal
! system, REVERSE the adjoint / transpose. A solver exposes ONE solve
! (an integrator ONE integrate) that takes this mode - not a second
! procedure.
!
! part - the subgraph a product acts on: the WHOLE operator, the
! self-loops (DIAGONAL), or one triangle of the neighbour edges
! (LOWER_TRIANGLE / UPPER_TRIANGLE). Stationary sweeps are traversals
! restricted to these subgraphs.
!
! Declared enumerations, compiler-visible - never bare integers at a
! call site.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module module_solve_mode

  implicit none

  private
  public :: FORWARD, REVERSE
  public :: WHOLE, DIAGONAL, LOWER_TRIANGLE, UPPER_TRIANGLE

  ! direction of edge traversal
  integer, parameter :: FORWARD = 1   ! primal:  A x = b
  integer, parameter :: REVERSE = 2   ! adjoint: A^T x = b

  ! subgraph selector for the one product
  integer, parameter :: WHOLE          =  2   ! the full operator
  integer, parameter :: DIAGONAL       =  0   ! self-loops
  integer, parameter :: LOWER_TRIANGLE = -1   ! one triangle of the neighbour edges
  integer, parameter :: UPPER_TRIANGLE =  1   ! the other triangle

end module module_solve_mode
