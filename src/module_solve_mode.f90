!=====================================================================!
! Named constants shared by the solvers, the integrator and the graph.
!
! mode - the direction of a solve or a product: FORWARD drives the
! primal system, REVERSE the adjoint / transpose. A solver exposes one
! solve (an integrator one integrate) that takes this mode - not a
! second procedure.
!
! part - the portion of the operator a product acts on: the WHOLE
! operator, the DIAGONAL, or one triangle (LOWER_TRIANGLE /
! UPPER_TRIANGLE). Stationary methods work on these parts.
!
! Use the named constants at call sites, never their integer values.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module module_solve_mode

  implicit none

  private
  public :: FORWARD, REVERSE
  public :: WHOLE, DIAGONAL, LOWER_TRIANGLE, UPPER_TRIANGLE

  ! solve / product direction
  integer, parameter :: FORWARD = 1   ! primal:  A x = b
  integer, parameter :: REVERSE = 2   ! adjoint: A^T x = b

  ! operator-part selector for the jacobian-vector product
  integer, parameter :: WHOLE          =  2   ! the full operator
  integer, parameter :: DIAGONAL       =  0   ! the diagonal
  integer, parameter :: LOWER_TRIANGLE = -1   ! the strictly lower triangle
  integer, parameter :: UPPER_TRIANGLE =  1   ! the strictly upper triangle

end module module_solve_mode
