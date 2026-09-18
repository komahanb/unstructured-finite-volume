!=====================================================================!
! THE FINITE-DIFFERENCE METHOD ON A MESH: at every cell, the formula
! over the cell's neighbours that is exact on the polynomials of
! degree p, one formula per derivative the residual names. The
! weights are those of the fit of the polynomial form on the
! neighbourhood (operation_fitting, operation_fitted_balance), grown
! ring by ring until the form's members are independent on it; this
! module names the method and its order p.
!
! THE ORDER OF ACCURACY. A derivative of order m read from the fit of
! degree p has a truncation error O(h^(p + 1 - m)): O(h^p) for a first
! derivative, O(h^(p - 1)) for a second. On a neighbourhood symmetric
! about the cell the odd and even members are orthogonal, so an even
! p gives O(h^p) for both, and an odd p gives the second derivatives
! of p - 1 on the wider neighbourhood p needs. Any p from 2 is
! admitted; the cost and the accuracy follow from this statement.
!
! At order two the form is restricted to the pure powers 1, x, x^2,
! y, y^2: on a cell and its face neighbours their second derivatives
! are the central differences and the laplacian's kernel is the
! constants alone, whereas the form with the mixed member over two
! rings has a grid mode in its kernel.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_finite_difference

  use util_precision          , only : dp
  use field_forms             , only : polynomial_form
  use view_mesh               , only : mesh
  use operation_stencil       , only : stencil
  use operation_fitted_balance, only : fitted_derivative_stencil

  implicit none

  private
  public :: finite_difference

  type :: finite_difference

     integer :: order = 2

   contains

     procedure :: derivative_rows

  end type finite_difference

  interface finite_difference
     module procedure create
  end interface finite_difference

contains

  function create(order) result(this)

    integer, intent(in) :: order
    type(finite_difference) :: this

    character(len=250) :: message

    if (order < 2) then
       write(message,'(a,i0)') 'operation_finite_difference: a second derivative requires order two at &
            &least; order = ', order
       error stop trim(message)
    end if
    this % order = order

  end function create

  !===================================================================!
  ! The rows of the derivative of one order along one axis of the
  ! mesh, at every cell: a stencil with one row per cell over its
  ! neighbours and no constant. An axis outside the mesh's dimension,
  ! or a derivative above the order, stops the program.
  !===================================================================!

  function derivative_rows(this, cells, axis, derivative) result(rows)

    class(finite_difference), intent(in) :: this
    type(mesh)              , intent(in) :: cells
    integer                 , intent(in) :: axis, derivative
    type(stencil) :: rows

    type(polynomial_form) :: shape
    integer, allocatable :: orders(:), pure(:)
    character(len=250) :: message

    if (axis < 1 .or. axis > cells % dimension) then
       write(message,'(a,i0,a,i0)') 'operation_finite_difference: the axis must be one of the mesh; &
            &axis = ', axis, ', dimension = ', cells % dimension
       error stop trim(message)
    end if
    if (derivative < 1 .or. derivative > this % order) then
       write(message,'(a,i0,a,i0)') 'operation_finite_difference: the derivative must lie within the &
            &order; derivative = ', derivative, ', order = ', this % order
       error stop trim(message)
    end if

    shape = polynomial_form(this % order, cells % dimension)
    if (this % order == 2) then
       call shape % pure_members(pure)
       call shape % restrict(pure)
    end if
    allocate(orders(cells % dimension), source=0)
    orders(axis) = derivative
    rows = fitted_derivative_stencil(cells, shape, orders)

  end function derivative_rows

end module operation_finite_difference
