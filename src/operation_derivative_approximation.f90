!=====================================================================!
! THE APPROXIMATION OF THE DERIVATIVES ALONG THE COORDINATES OF A
! REGION: on a mesh of cells, how the jet an unknown's law reads is
! formed, and how the law is tested. The methods extend this one type,
! and the residual reads every method through its three procedures:
!
!   derivative_rows   the relations among the components a point
!                     stores: for each axis and each derivative the
!                     law names, a stencil with one row per cell over
!                     the cells it reads, exact on the polynomials of
!                     a degree set by the order of accuracy. Finite
!                     differences fit a polynomial at the cell; finite
!                     volumes integrate over the faces of the cell.
!   projects          whether the law is tested against a basis on
!                     the cells, or at the points.
!   elements          the cells as finite elements, for a method that
!                     projects: their nodes, their quadrature points
!                     and weights, and the jet of the basis at each.
!
! A method that places the law at the points has no elements: each of
! its points is an element of one node, whose basis is the identity
! on the stored jet and whose one quadrature point, of weight one, is
! the node.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_derivative_approximation

  use util_precision   , only : dp
  use view_mesh        , only : mesh
  use operation_stencil, only : stencil

  implicit none

  private
  public :: derivative_approximation, element_quadrature

  !===================================================================!
  ! THE ELEMENTS OF A MESH on which a law is projected. A cell has
  ! nodes, points of the mesh, at each of which an unknown stores its
  ! coefficients: its value and its pure derivatives along each axis
  ! to the element's order. A cell has quadrature points, each with a
  ! weight, a position and the jet of every basis function there:
  ! basis(d, k, i, q, e) is the derivative d of the basis function of
  ! the coefficient k of node i of cell e, at its quadrature point q.
  ! A derivative and a coefficient are indexed alike: zero the value,
  ! (a - 1) * top + o the order o along the axis a, top the highest
  ! order stored (derivatives for the jet, order for the
  ! coefficients). The jet of an unknown at a quadrature point is the
  ! sum over the cell's nodes and coefficients of the basis's jet
  ! times the coefficient.
  !
  ! THE TEST FUNCTION of a coefficient is its basis function divided
  ! by its moment of the coefficient's order about the node, the
  ! integral of the basis function times (x - x_i)^o / o! along the
  ! coefficient's axis: test(d, k, i, q, e) is the weight times the
  ! test function's derivative d. The row of a coefficient is then the
  ! law's component of the coefficient's order at the node, in the weak
  ! sense: the law's average about the node for the value, its slope
  ! for the slope. Tested with a point mass at the node, the row would
  ! be the law's component at the point, the row of a collocation.
  !===================================================================!

  type :: element_quadrature

     integer :: dimension   = 1
     ! the order of the coefficients along each axis, and the highest
     ! order of the basis's jet stored at the quadrature points
     integer :: order       = 0
     integer :: derivatives = 0
     integer , allocatable :: node(:,:)
     real(dp), allocatable :: weight(:,:)
     real(dp), allocatable :: position(:,:,:)
     real(dp), allocatable :: basis(:,:,:,:,:)
     real(dp), allocatable :: test(:,:,:,:,:)

   contains

     procedure :: num_cells             => element_num_cells
     procedure :: num_nodes             => element_num_nodes
     procedure :: num_quadrature_points => element_num_quadrature_points
     procedure :: num_coefficients      => element_num_coefficients
     procedure :: num_derivatives       => element_num_derivatives

  end type element_quadrature

  type, abstract :: derivative_approximation

     ! the order of accuracy p: a derivative from the rows is O(h^p)
     integer :: order = 2

   contains

     procedure(derivative_rows_interface), deferred :: derivative_rows
     procedure :: projects
     procedure :: elements

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

contains

  !===================================================================!
  ! A method places the law at the points unless it states otherwise.
  !===================================================================!

  pure logical function projects(this)
    class(derivative_approximation), intent(in) :: this
    associate (u1 => this); end associate
    projects = .false.
  end function projects

  !===================================================================!
  ! The elements of a method that places the law at the points: none.
  ! A method that projects states its cells, with the coefficients to
  ! the order given and the basis's jet to the derivatives given.
  !===================================================================!

  function elements(this, cells, order, derivatives) result(set)
    class(derivative_approximation), intent(in) :: this
    type(mesh)                     , intent(in) :: cells
    integer                        , intent(in) :: order, derivatives
    type(element_quadrature) :: set
    associate (u1 => this, u2 => cells); end associate
    set % dimension   = 1
    set % order       = order
    set % derivatives = derivatives
    allocate(set % node(0, 0), set % weight(0, 0), set % position(1, 0, 0), set % basis(0:0, 0:0, 0, 0, 0))
    allocate(set % test(0:0, 0:0, 0, 0, 0))
  end function elements

  pure integer function element_num_cells(this)
    class(element_quadrature), intent(in) :: this
    element_num_cells = 0
    if (allocated(this % node)) element_num_cells = size(this % node, 2)
  end function element_num_cells

  pure integer function element_num_nodes(this)
    class(element_quadrature), intent(in) :: this
    element_num_nodes = size(this % node, 1)
  end function element_num_nodes

  pure integer function element_num_quadrature_points(this)
    class(element_quadrature), intent(in) :: this
    element_num_quadrature_points = size(this % weight, 1)
  end function element_num_quadrature_points

  ! the coefficients a node stores of one unknown: the value and the
  ! pure derivatives to the order along each axis
  pure integer function element_num_coefficients(this)
    class(element_quadrature), intent(in) :: this
    element_num_coefficients = 1 + this % dimension * this % order
  end function element_num_coefficients

  ! the derivatives of the basis's jet stored at a quadrature point
  pure integer function element_num_derivatives(this)
    class(element_quadrature), intent(in) :: this
    element_num_derivatives = 1 + this % dimension * this % derivatives
  end function element_num_derivatives

end module operation_derivative_approximation
