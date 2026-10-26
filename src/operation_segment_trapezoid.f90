!=====================================================================!
! Module that states the trapezoid rule on the segments of a line of
! points. The derivatives along a space coordinate are represented by
! relations between the jets at the two ends of each segment,
!
!      c(i+1) - c(i) - h/2 (c'(i+1) + c'(i)) = 0,
!
! one for each consecutive pair (c, c') of the chain w, w_x, w_xx, ...
! of the unknown, of order two in the segment length h:
!
!          c'
!          |                o c'(i+1)
!          |           _--- |
!   c'(i)  o------ ---      |      the area under c' over the segment,
!          |                |      c(i+1) - c(i), taken as the area of
!          |                |      the trapezoid, h/2 (c'(i) + c'(i+1))
!          +----------------+----- x
!         x(i)    h      x(i+1)
!
! A line of P points and a chain of degree D have D (P - 1) relations,
! so D rows stay open; the conditions on the two faces of the line
! take them:
!
!   points     1      2      3     ...     P
!              o------o------o---- ... ----o
!   relations     D      D      D   ...       = D (P - 1)
!   faces      ^                           ^  the D open rows
!
! Where the values of a moment are determined without the physics,
! as at the initial instant and at an explicit instant of a family,
! the relations would only differentiate known values and leave the
! chain short of data. There the derivatives are the closed fitted
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

  !-------------------------------------------------------------------!
  ! Type that models the trapezoid relations on a line of points, with
  ! the closed differences for the moments whose values are known
  !-------------------------------------------------------------------!

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

  !===================================================================!
  ! Creates the rule, of order two, with closed differences of the
  ! order given, or of order four
  !===================================================================!

  function create(closed_order) result(this)

    integer, intent(in), optional :: closed_order
    type(segment_trapezoid)       :: this

    this % order  = 2
    this % closed = finite_difference(order=4)
    if (present(closed_order)) this % closed = finite_difference(order=closed_order)

  end function create

  !===================================================================!
  ! Returns the closed rows of a derivative at every point, for the
  ! moments whose values are known: those of the fitted differences
  !===================================================================!

  function derivative_rows(this, cells, axis, derivative) result(rows)

    class(segment_trapezoid), intent(in) :: this
    type(mesh)              , intent(in) :: cells
    integer                 , intent(in) :: axis, derivative
    type(stencil)                        :: rows

    rows = this % closed % derivative_rows(cells, axis, derivative)

  end function derivative_rows

end module operation_segment_trapezoid
