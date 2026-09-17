!=====================================================================!
! THE FINITE-DIFFERENCE METHOD ON A MESH: at every cell, the formula
! over the cell's neighbours that is exact on the polynomials of a
! degree, one formula per derivative the residual names. The weights
! are those of the fit of the polynomial form on the neighbourhood
! (operation_fitting, operation_fitted_balance); this module names
! the method and its degree.
!
! At degree two the form is restricted to the pure powers 1, x, x^2,
! y, y^2: on a cell and its face neighbours their second derivatives
! are the central differences and the laplacian's kernel is the
! constants alone, whereas the form with the mixed member over two
! rings has a grid mode in its kernel. An odd degree fits the second
! derivatives of the even degree below it over a neighbourhood with
! that mode, and is refused.
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

     integer :: degree = 2

   contains

     procedure :: derivative_rows

  end type finite_difference

  interface finite_difference
     module procedure create
  end interface finite_difference

contains

  function create(degree) result(this)

    integer, intent(in) :: degree
    type(finite_difference) :: this

    character(len=250) :: message

    if (degree < 2) then
       write(message,'(a,i0)') 'operation_finite_difference: a second derivative requires degree two at &
            &least; degree = ', degree
       error stop trim(message)
    end if
    if (mod(degree, 2) == 1) then
       write(message,'(a,i0)') 'operation_finite_difference: an odd degree fits the second derivatives of &
            &the even degree below it over a neighbourhood whose laplacian has a grid mode in its &
            &kernel; degree = ', degree
       error stop trim(message)
    end if
    this % degree = degree

  end function create

  !===================================================================!
  ! The rows of the derivative of one order along one axis of the
  ! mesh, at every cell: a stencil with one row per cell over its
  ! neighbours and no constant. An axis outside the mesh's dimension,
  ! or an order above the degree, stops the program.
  !===================================================================!

  function derivative_rows(this, cells, axis, order) result(rows)

    class(finite_difference), intent(in) :: this
    type(mesh)              , intent(in) :: cells
    integer                 , intent(in) :: axis, order
    type(stencil) :: rows

    type(polynomial_form) :: shape
    integer, allocatable :: orders(:), pure(:)
    character(len=250) :: message

    if (axis < 1 .or. axis > cells % dimension) then
       write(message,'(a,i0,a,i0)') 'operation_finite_difference: the axis must be one of the mesh; &
            &axis = ', axis, ', dimension = ', cells % dimension
       error stop trim(message)
    end if
    if (order < 1 .or. order > this % degree) then
       write(message,'(a,i0,a,i0)') 'operation_finite_difference: the order must lie within the &
            &degree; order = ', order, ', degree = ', this % degree
       error stop trim(message)
    end if

    shape = polynomial_form(this % degree, cells % dimension)
    if (this % degree == 2) then
       call shape % pure_members(pure)
       call shape % restrict(pure)
    end if
    allocate(orders(cells % dimension), source=0)
    orders(axis) = order
    rows = fitted_derivative_stencil(cells, shape, orders)

  end function derivative_rows

end module operation_finite_difference
