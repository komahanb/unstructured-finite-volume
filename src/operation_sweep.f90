!=====================================================================!
! This module holds named constants shared by the solvers, the
! integrator and the graph.
!
! The mode is the direction of a solve or a product: FORWARD drives
! the primal system, and REVERSE drives the adjoint / transpose. A
! solver exposes one solve (an integrator, one integrate) that takes
! this mode rather than a second procedure.
!
! The part is the portion of the operator a product acts on: the
! WHOLE operator, the DIAGONAL, or one triangle (LOWER_TRIANGLE /
! UPPER_TRIANGLE). Stationary methods work on these parts.
!
! Use the named constants at call sites, never their integer values.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_sweep

  implicit none

  private
  public :: FORWARD, REVERSE
  public :: WHOLE, DIAGONAL, LOWER_TRIANGLE, UPPER_TRIANGLE
  public :: is_valid_mode, is_valid_part

  ! These constants name the solve or product direction.
  integer, parameter :: FORWARD = 1   ! The primal solve:  A x = b.
  integer, parameter :: REVERSE = 2   ! The adjoint solve: A^T x = b.

  !-------------------------------------------------------------------!
  ! These constants select the part of the operator a jacobian-vector
  ! product acts on. The range is disjoint from the mode values, so a
  ! part passed where a mode is expected (or the reverse) can never be
  ! silently reinterpreted - the entry validators below refuse it by
  ! name. The part values are only ever compared for equality; the
  ! triangle semantics come from row and column comparisons in the
  ! discretization, not from these values.
  !-------------------------------------------------------------------!

  integer, parameter :: WHOLE          = 11   ! The full operator.
  integer, parameter :: DIAGONAL       = 12   ! The diagonal.
  integer, parameter :: LOWER_TRIANGLE = 13   ! The strictly lower triangle.
  integer, parameter :: UPPER_TRIANGLE = 14   ! The strictly upper triangle.

contains

  !===================================================================!
  ! This entry validator accepts a mode tag: it is either FORWARD or
  ! REVERSE, or it is refused loudly at the door - an illegal value
  ! is never reinterpreted.
  !===================================================================!

  pure logical function is_valid_mode(m)

    integer, intent(in) :: m

    is_valid_mode = (m .eq. FORWARD) .or. (m .eq. REVERSE)

  end function is_valid_mode

  !===================================================================!
  ! This entry validator accepts a part tag: it is one of WHOLE,
  ! DIAGONAL, LOWER_TRIANGLE or UPPER_TRIANGLE, or it is refused
  ! loudly at the door - an illegal value is never reinterpreted.
  !===================================================================!

  pure logical function is_valid_part(p)

    integer, intent(in) :: p

    is_valid_part = (p .eq. WHOLE) .or. (p .eq. DIAGONAL) .or. &
         &          (p .eq. LOWER_TRIANGLE) .or. (p .eq. UPPER_TRIANGLE)

  end function is_valid_part

end module operation_sweep
