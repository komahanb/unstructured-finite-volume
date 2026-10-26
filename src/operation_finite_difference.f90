!=====================================================================!
! Module that states the finite-difference method on a mesh. At every
! cell it forms, over the cell's neighbours, the formula that is exact
! on the polynomials of degree p, one formula per derivative the
! residual names:
!
!                 o   o   o
!                   \ | /          a ring of neighbours, grown ring
!             o --- cell --- o     by ring until the form's members
!                   / | \          are independent on it
!                 o   o   o
!
! The weights are those of the fit of the polynomial form on the
! neighbourhood (operation_fitting, operation_fitted_balance). This
! module names the method and its order p.
!
! The order of accuracy p, and the degree of the fit each derivative
! reads. A derivative of order m from a fit of degree d has a
! truncation error O(h^(d + 1 - m)). On a neighbourhood symmetric
! about the cell the odd and even members are orthogonal, so an even d
! gives O(h^d) for the first and the second derivative alike:
!
!   derivative   degree of the fit read          error
!   first        p                               O(h^p)
!   second       p, rounded up to even           O(h^p)
!
! With that, every order from two is O(h^p) in both, and the second
! derivatives never come from an odd degree. The laplacian of an odd
! degree on the quadrilateral box has the odd-even grid mode in its
! kernel, so that the discrete equations are singular:
!
!   +1  -1  +1  -1
!   -1  +1  -1  +1     the odd-even grid mode
!   +1  -1  +1  -1
!
! Measured: at order 3, as a plain cubic fit, every linear solve
! reaches its cycle limit with the residual unreduced.
!
! At degree two the form is restricted to the pure powers 1, x, x^2,
! y, y^2. On a cell and its face neighbours their second derivatives
! are the central differences, and the laplacian's kernel is the
! constants alone; the form with the mixed member over two rings has a
! grid mode in its kernel:
!
!           o
!           |
!      o -- c -- o      1, x, x^2, y, y^2 on the cell c and its four
!           |           face neighbours: the central differences
!           o
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_finite_difference

  use util_precision                    , only : dp
  use field_forms                       , only : polynomial_form
  use view_mesh                         , only : mesh
  use operation_stencil                 , only : stencil
  use operation_fitted_balance          , only : fitted_derivative_stencil
  use operation_derivative_approximation, only : derivative_approximation

  implicit none

  private

  public :: finite_difference

  !-------------------------------------------------------------------!
  ! Type that models the finite-difference method of an order p
  !-------------------------------------------------------------------!

  type, extends(derivative_approximation) :: finite_difference

   contains

     procedure :: derivative_rows

  end type finite_difference

  interface finite_difference
     module procedure create
  end interface finite_difference

contains

  !===================================================================!
  ! Creates the method of the given order, at least two, so that it
  ! can form a second derivative
  !===================================================================!

  function create(order) result(this)

    integer, intent(in)     :: order
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
  ! Returns the rows of the derivative of one order along one axis of
  ! the mesh, at every cell: a stencil with one row per cell over its
  ! neighbours and no constant. An axis outside the mesh's dimension,
  ! or a derivative above the order, stops the program.
  !===================================================================!

  function derivative_rows(this, cells, axis, derivative) result(rows)

    class(finite_difference), intent(in) :: this
    type(mesh)              , intent(in) :: cells
    integer                 , intent(in) :: axis, derivative
    type(stencil)                        :: rows

    type(polynomial_form) :: shape
    integer, allocatable  :: orders(:), pure(:)
    integer               :: degree
    character(len=250)    :: message

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

    ! The degree read: the order for a first derivative, the order
    ! rounded up to even for a second
    degree = this % order
    if (derivative >= 2) degree = degree + mod(degree, 2)

    ! The form, of the pure powers alone at degree two
    shape = polynomial_form(degree, cells % dimension)
    if (degree == 2) then
       call shape % pure_members(pure)
       call shape % restrict(pure)
    end if
    allocate(orders(cells % dimension), source=0)
    orders(axis) = derivative
    rows = fitted_derivative_stencil(cells, shape, orders)

  end function derivative_rows

end module operation_finite_difference
