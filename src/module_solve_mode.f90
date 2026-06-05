!=====================================================================!
! Solve direction shared by the linear solvers and the time integrator:
! FORWARD drives the primal system, REVERSE the adjoint / transpose. A
! solver exposes ONE solve (an integrator ONE integrate) that takes this
! mode - not a second procedure.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module module_solve_mode

  implicit none

  private
  public :: FORWARD, REVERSE

  integer, parameter :: FORWARD = 1   ! primal:  A x = b
  integer, parameter :: REVERSE = 2   ! adjoint: A^T x = b

end module module_solve_mode
