!=====================================================================!
! THE FINITE ELEMENT METHOD ON A MESH: the law tested against a basis
! on the cells, in place of being placed at the points. The unknown
! and its multiplier are expanded in the same basis, and each row is
! the law tested with one basis function: an inner product, which is
! an integral, which the quadrature of the cells approximates.
!
! THE ELEMENT of degree p on a cell: the polynomials of degree p, with
! the jet of the unknown at the cell's nodes as its coefficients, to
! the order the law requires. The basis is dual to the coefficients:
! each basis function is the polynomial that its own coefficient reads
! as one and every other as zero, a column of the inverse of the
! square matrix of the coefficients of the polynomials. On a line the
! cells are the segments between consecutive points and the nodes are
! their ends: the value alone at degree one, the linear element; the
! value and the slope at degree three, the Hermite cubic. The basis is
! formed once on the reference segment [0, 1] and mapped to each
! segment of length h: the coefficient of order k along the line is
! the reference's divided by h^k, and a derivative of order o is the
! reference's divided by h^o.
!
! THE QUADRATURE of a segment: Gauss-Legendre of p + 1 points, exact
! on the polynomials of degree 2p + 1, among them the product of two
! basis functions.
!
! THE TEST FUNCTION of a coefficient: its basis function divided by
! its moment of the coefficient's order about the node, summed over
! the node's cells - for the value the integral of the basis function,
! the node's volume; for the slope the integral of the basis function
! times (x - x_i), h^3 / 30 on each side of the node for the Hermite
! cubic. A row is then the law's value, or its slope, about the node.
!
! The coefficients of a node are independent unknowns, each determined
! by the law tested with its own basis function, so the method states
! no relation among the components a point stores.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_finite_element

  use util_precision    , only : dp
  use util_factorisation, only : dense_factorisation
  use field_forms       , only : polynomial_form
  use view_mesh         , only : mesh, values_of
  use operation_stencil , only : stencil
  use operation_grid    , only : gauss_node_weight
  use operation_derivative_approximation, only : derivative_approximation, element_quadrature

  implicit none

  private
  public :: finite_element

  type, extends(derivative_approximation) :: finite_element

     ! the degree of the polynomials of an element
     integer :: degree = 1

   contains

     procedure :: derivative_rows
     procedure :: projects
     procedure :: elements

  end type finite_element

  interface finite_element
     module procedure create
  end interface finite_element

contains

  function create(degree) result(this)

    integer, intent(in) :: degree
    type(finite_element) :: this

    character(len=250) :: message

    if (degree < 1) then
       write(message,'(a,i0)') 'operation_finite_element: an element is of degree one at least; degree = ', degree
       error stop trim(message)
    end if
    this % degree = degree
    this % order  = degree + 1

  end function create

  !===================================================================!
  ! No relation among the components a point stores: a stencil with
  ! no entries over the cells. An axis outside the mesh's dimension
  ! stops the program.
  !===================================================================!

  function derivative_rows(this, cells, axis, derivative) result(rows)

    class(finite_element), intent(in) :: this
    type(mesh)           , intent(in) :: cells
    integer              , intent(in) :: axis, derivative
    type(stencil) :: rows

    character(len=250) :: message

    associate (u1 => this, u2 => derivative); end associate
    if (axis < 1 .or. axis > cells % dimension) then
       write(message,'(a,i0,a,i0)') 'operation_finite_element: the axis must be one of the mesh; &
            &axis = ', axis, ', dimension = ', cells % dimension
       error stop trim(message)
    end if
    rows = stencil([integer ::], [integer ::], [real(dp) ::], spread(0.0_dp, 1, cells % num_vertices()), &
         & 'no relation')

  end function derivative_rows

  pure logical function projects(this)
    class(finite_element), intent(in) :: this
    associate (u1 => this); end associate
    projects = .true.
  end function projects

  !===================================================================!
  ! THE ELEMENTS OF A LINE: each segment between consecutive points a
  ! cell, its two ends its nodes, each storing the coefficients of
  ! orders zero to order, and its p + 1 Gauss-Legendre points, at each
  ! the jet of the basis to the order derivatives. Invalid input: a
  ! mesh of several coordinates, or a degree other than 2 order + 1,
  ! the degree of the polynomials the coefficients of two nodes
  ! determine.
  !===================================================================!

  function elements(this, cells, order, derivatives) result(set)

    class(finite_element), intent(in) :: this
    type(mesh)           , intent(in) :: cells
    integer              , intent(in) :: order, derivatives
    type(element_quadrature) :: set

    real(dp), allocatable :: centre(:), reference(:,:,:), xi(:), w(:), moment(:,:)
    integer , allocatable :: lower(:), upper(:)
    integer :: ncells, nq, e, q, i, k, o, n, p
    real(dp) :: h
    character(len=250) :: message

    if (cells % dimension /= 1) then
       error stop 'operation_finite_element: the elements of a region of several coordinates are not formed; &
            &the elements of a line are'
    end if
    if (this % degree /= 2 * order + 1) then
       write(message,'(a,i0,a,i0,a,i0)') 'operation_finite_element: on a line the coefficients of order ', &
            & order, ' at two nodes determine the polynomials of degree ', 2 * order + 1, '; degree = ', &
            & this % degree
       error stop trim(message)
    end if

    ! the segments: every face with a cell at each end, from the lower
    ! point to the upper
    call values_of(cells % cell_centres, centre)
    allocate(lower(0), upper(0))
    do e = 1, cells % num_edges()
       if (.not. cells % edge_has_head(e)) cycle
       i = cells % edge_tail(e)
       k = cells % edge_head(e)
       if (centre(i) < centre(k)) then
          lower = [lower, i]
          upper = [upper, k]
       else
          lower = [lower, k]
          upper = [upper, i]
       end if
    end do
    ncells = size(lower)

    ! the reference basis, its jet at the reference quadrature points
    nq = this % degree + 1
    allocate(xi(nq), w(nq))
    do q = 1, nq
       call gauss_node_weight(q, nq, 1.0_dp, xi(q), w(q))
    end do
    call reference_jet(order, derivatives, xi, reference)

    set % dimension   = 1
    set % order       = order
    set % derivatives = derivatives
    n = order + 1
    allocate(set % node(2, ncells), set % weight(nq, ncells), set % position(1, nq, ncells))
    allocate(set % basis(0:derivatives, 0:order, 2, nq, ncells))
    do e = 1, ncells
       h = centre(upper(e)) - centre(lower(e))
       set % node(:, e) = [lower(e), upper(e)]
       do q = 1, nq
          set % weight(q, e)      = h * w(q)
          set % position(1, q, e) = centre(lower(e)) + h * xi(q)
          do i = 1, 2
             do k = 0, order
                do o = 0, derivatives
                   set % basis(o, k, i, q, e) = reference((i - 1) * n + k + 1, o, q) * h**(k - o)
                end do
             end do
          end do
       end do
    end do

    ! the moment of each coefficient's basis function about its node,
    ! and the test functions
    allocate(moment(0:order, size(centre)), source=0.0_dp)
    do e = 1, ncells
       do i = 1, 2
          p = set % node(i, e)
          do k = 0, order
             do q = 1, nq
                moment(k, p) = moment(k, p) + set % weight(q, e) * set % basis(0, k, i, q, e) &
                     & * (set % position(1, q, e) - centre(p))**k / real(product([(o, o = 1, k)]), dp)
             end do
          end do
       end do
    end do
    allocate(set % test, mold=set % basis)
    do e = 1, ncells
       do q = 1, nq
          do i = 1, 2
             do k = 0, order
                set % test(:, k, i, q, e) = set % weight(q, e) * set % basis(:, k, i, q, e) / moment(k, set % node(i, e))
             end do
          end do
       end do
    end do

  end function elements

  !===================================================================!
  ! THE BASIS ON THE REFERENCE SEGMENT [0, 1], dual to the coefficients
  ! of orders zero to order at its ends: jet(l, o, q) is the derivative
  ! of order o, at the point xi(q), of the polynomial of degree
  ! 2 order + 1 that the coefficient l reads as one and every other as
  ! zero, the coefficients numbered by node and then by order. The
  ! polynomials are the members of polynomial_form about zero, and the
  ! square matrix of their coefficients is inverted by its
  ! factorisation.
  !===================================================================!

  subroutine reference_jet(order, derivatives, xi, jet)

    integer , intent(in) :: order, derivatives
    real(dp), intent(in) :: xi(:)
    real(dp), allocatable, intent(out) :: jet(:,:,:)

    type(polynomial_form) :: shape
    type(dense_factorisation) :: factors
    real(dp), allocatable :: v(:,:), member(:), unit(:), a(:), dphi(:)
    real(dp) :: ends(2)
    integer :: n, l, i, k, o, q

    n     = 2 * (order + 1)
    shape = polynomial_form(n - 1, 1)
    ends  = [0.0_dp, 1.0_dp]
    allocate(v(n, n), member(n), dphi(n), unit(n))
    do i = 1, 2
       do k = 0, order
          call shape % derivatives([ends(i)], [0.0_dp], [k], member)
          v((i - 1) * (order + 1) + k + 1, :) = member
       end do
    end do
    call factors % factorise(v, 0.0_dp)
    if (factors % singular()) then
       error stop 'operation_finite_element: the coefficients at the two ends do not determine the polynomials'
    end if

    allocate(jet(n, 0:derivatives, size(xi)))
    do l = 1, n
       unit    = 0.0_dp
       unit(l) = 1.0_dp
       call factors % substitute(unit, a, transposed=.false.)
       do q = 1, size(xi)
          do o = 0, derivatives
             call shape % derivatives([xi(q)], [0.0_dp], [o], dphi)
             jet(l, o, q) = dot_product(a, dphi)
          end do
       end do
    end do

  end subroutine reference_jet

end module operation_finite_element
