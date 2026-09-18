!=====================================================================!
! THE APPROXIMATION OF THE DERIVATIVES ALONG THE COORDINATES OF A
! REGION: on a mesh of cells, for each axis and each derivative the
! residual names, a stencil with one row per cell over the cells it
! reads, exact on the polynomials of a degree set by the order of
! accuracy. The methods extend this one type and differ in how the
! rows are formed: finite differences fit a polynomial at the cell;
! finite volumes integrate over the faces of the cell. The residual
! reads every method through derivative_rows alone.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_derivative_approximation

  use view_mesh        , only : mesh
  use operation_stencil, only : stencil

  implicit none

  private
  public :: derivative_approximation

  type, abstract :: derivative_approximation

     ! the order of accuracy p: a derivative from the rows is O(h^p)
     integer :: order = 2

   contains

     procedure(derivative_rows_interface), deferred :: derivative_rows

  end type derivative_approximation

  abstract interface
     !================================================================!
     ! The rows of the derivative of one order along one axis of the
     ! mesh, at every cell: a stencil with one row per cell over the
     ! cells it reads and no constant. An axis outside the mesh's
     ! dimension, or a derivative the method does not form, stops
     ! the program.
     !================================================================!
     function derivative_rows_interface(this, cells, axis, derivative) result(rows)
       import :: derivative_approximation, mesh, stencil
       class(derivative_approximation), intent(in) :: this
       type(mesh)                     , intent(in) :: cells
       integer                        , intent(in) :: axis, derivative
       type(stencil) :: rows
     end function derivative_rows_interface
  end interface

end module operation_derivative_approximation
