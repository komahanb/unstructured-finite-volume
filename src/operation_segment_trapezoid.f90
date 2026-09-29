!=====================================================================!
! THE TRAPEZOID RULE ON THE SEGMENTS OF A LINE OF POINTS: the
! derivatives along a space coordinate represented by relations
! between the jets at the two ends of each segment,
!
!      c(i+1) - c(i) - h/2 (c'(i+1) + c'(i)) = 0,
!
! one for each consecutive pair (c, c') of the chain w, w_x, w_xx, ..
! of the unknown, of order two in the segment length h. The relations
! of a line of P points and a chain of degree D number D (P - 1), so D
! rows stay open: the conditions on the two faces of the line take
! them. Where the values of a moment are determined without the
! physics - the initial instant, an explicit instant of a family - the
! relations would only differentiate known values and leave the chain
! short of data, so there the derivatives are the closed fitted
! differences of the values, of the order given (four when none is).
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_segment_trapezoid

  use view_mesh                         , only : mesh
  use operation_stencil                 , only : stencil
  use operation_derivative_approximation, only : derivative_approximation
  use operation_finite_difference       , only : finite_difference

  implicit none

  private
  public :: segment_trapezoid

  type, extends(derivative_approximation) :: segment_trapezoid
     ! the closed differences at the moments whose values are known
     type(finite_difference) :: closed
   contains
     procedure :: derivative_rows
  end type segment_trapezoid

  interface segment_trapezoid
     module procedure create
  end interface segment_trapezoid

contains

  function create(closed_order) result(this)

    integer, intent(in), optional :: closed_order
    type(segment_trapezoid) :: this

    this % order  = 2
    this % closed = finite_difference(order=4)
    if (present(closed_order)) this % closed = finite_difference(order=closed_order)

  end function create

  !===================================================================!
  ! The closed rows of a derivative at every point, for the moments
  ! whose values are known: those of the fitted differences.
  !===================================================================!

  function derivative_rows(this, cells, axis, derivative) result(rows)

    class(segment_trapezoid), intent(in) :: this
    type(mesh)              , intent(in) :: cells
    integer                 , intent(in) :: axis, derivative
    type(stencil) :: rows

    rows = this % closed % derivative_rows(cells, axis, derivative)

  end function derivative_rows

end module operation_segment_trapezoid
