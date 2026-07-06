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
  public :: is_valid_mode, is_valid_part

  ! solve / product direction
  integer, parameter :: FORWARD = 1   ! primal:  A x = b
  integer, parameter :: REVERSE = 2   ! adjoint: A^T x = b

  ! operator-part selector for the jacobian-vector product. the range
  ! is disjoint from the mode values, so a part passed where a mode is
  ! expected (or the reverse) can never be silently reinterpreted - the
  ! entry validators below refuse it by name. the part values are only
  ! ever compared for equality (the triangle semantics come from row/
  ! column comparisons in the discretization, not from these values).
  integer, parameter :: WHOLE          = 11   ! the full operator
  integer, parameter :: DIAGONAL       = 12   ! the diagonal
  integer, parameter :: LOWER_TRIANGLE = 13   ! the strictly lower triangle
  integer, parameter :: UPPER_TRIANGLE = 14   ! the strictly upper triangle

contains

  !===================================================================!
  ! Entry validators: a tag is either one of the named constants or it
  ! is refused loudly at the door - an illegal value is never
  ! reinterpreted.
  !===================================================================!

  pure logical function is_valid_mode(m)

    integer, intent(in) :: m

    is_valid_mode = (m .eq. FORWARD) .or. (m .eq. REVERSE)

  end function is_valid_mode

  pure logical function is_valid_part(p)

    integer, intent(in) :: p

    is_valid_part = (p .eq. WHOLE) .or. (p .eq. DIAGONAL) .or. &
         &          (p .eq. LOWER_TRIANGLE) .or. (p .eq. UPPER_TRIANGLE)

  end function is_valid_part

end module module_solve_mode
