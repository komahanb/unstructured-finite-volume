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
! IN SEVERAL COORDINATES the cells are those of a mesh read from a file
! and its corners are the nodes, the points of the region; the element
! is of degree one, the value at its corners its coefficients: the
! polynomials of degree one on a triangle and a tetrahedron, the
! bilinear and trilinear ones on a quadrangle and a hexahedron, on the
! reference cell of unit side mapped to each cell by its own basis. The
! quadrature is the tensor product of Gauss-Legendre rules of three
! points, collapsed onto the reference triangle and tetrahedron.
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

    real(dp), allocatable :: centre(:), reference(:,:,:), xi(:), w(:)
    integer , allocatable :: lower(:), upper(:)
    integer :: ncells, nq, e, q, i, k, o, n
    real(dp) :: h
    character(len=250) :: message

    if (cells % dimension > 1) then
       set = corner_elements(this, cells, order, derivatives)
       return
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

    call tested(set, centre)

  end function elements

  !===================================================================!
  ! THE TEST FUNCTIONS of a set of elements: each coefficient's basis
  ! function divided by its moment of the coefficient's order about its
  ! node, summed over the node's cells, times the weight. centre is
  ! the positions of the points, the dimension's coordinates of each.
  !===================================================================!

  subroutine tested(set, centre)

    type(element_quadrature), intent(inout) :: set
    real(dp)                , intent(in)    :: centre(:)

    real(dp), allocatable :: moment(:,:)
    real(dp) :: factorial
    integer :: e, q, i, k, o, p, axis, d, m

    d = set % dimension
    allocate(moment(0:set % num_coefficients() - 1, size(centre) / d), source=0.0_dp)
    do e = 1, set % num_cells()
       do i = 1, set % num_nodes()
          p = set % node(i, e)
          do k = 0, set % num_coefficients() - 1
             axis = 1
             o    = 0
             if (k > 0) then
                axis = (k - 1) / set % order + 1
                o    = k - (axis - 1) * set % order
             end if
             factorial = 1.0_dp
             do m = 2, o
                factorial = factorial * real(m, dp)
             end do
             do q = 1, set % num_quadrature_points()
                moment(k, p) = moment(k, p) + set % weight(q, e) * set % basis(0, k, i, q, e) &
                     & * (set % position(axis, q, e) - centre(d * (p - 1) + axis))**o / factorial
             end do
          end do
       end do
    end do
    allocate(set % test, mold=set % basis)
    do e = 1, set % num_cells()
       do q = 1, set % num_quadrature_points()
          do i = 1, set % num_nodes()
             do k = 0, set % num_coefficients() - 1
                set % test(:, k, i, q, e) = set % weight(q, e) * set % basis(:, k, i, q, e) / moment(k, set % node(i, e))
             end do
          end do
       end do
    end do

  end subroutine tested

  !===================================================================!
  ! THE ELEMENTS OF A MESH OF SEVERAL COORDINATES: the cells the mesh of
  ! corners keeps, all of one kind - triangles, quadrangles, tetrahedra
  ! or hexahedra - each with its corners as its nodes, the value its
  ! one coefficient, and the jet of its basis to the first derivatives
  ! at the quadrature points. The basis on the reference cell is dual to
  ! the values at the reference corners; a cell is the image of the
  ! reference cell under its own basis, x(xi) = sum over the corners of
  ! psi_a(xi) x_a, and the derivatives along x are those along xi
  ! through the inverse transpose of that map's jacobian. Invalid input:
  ! coefficients of a derivative, a degree other than one, second
  ! derivatives, a mesh without its cells, or cells of several kinds or
  ! of another kind.
  !===================================================================!

  function corner_elements(this, cells, order, derivatives) result(set)

    class(finite_element), intent(in) :: this
    type(mesh)           , intent(in) :: cells
    integer              , intent(in) :: order, derivatives
    type(element_quadrature) :: set

    type(polynomial_form) :: shape
    type(dense_factorisation) :: factors
    real(dp), allocatable :: corner(:,:), xi(:,:), w(:), member(:), a(:,:), unit(:), coefficients(:), centre(:)
    real(dp), allocatable :: value(:,:), slope(:,:,:), dphi(:)
    integer , allocatable :: active(:)
    real(dp) :: jacobian(3, 3), inverse(3, 3), determinant, gauss(3), gauss_weight(3), x(3)
    integer :: d, nn, ne, nq, kind, e, q, i, j, r, l, m, u, v, s
    logical :: simplex

    if (order /= 0 .or. this % degree /= 1) then
       error stop 'operation_finite_element: in several coordinates the element is of degree one, the values at &
            &the corners its coefficients: a law of the first derivatives along space'
    end if
    if (derivatives > 1) then
       error stop 'operation_finite_element: in several coordinates the law in flux form reads the first &
            &derivatives along space at most'
    end if
    if (.not. cells % has_corners()) then
       error stop 'operation_finite_element: the elements of a region of several coordinates are the cells of a &
            &mesh read from a file, its corners the points: corners(mesh(file))'
    end if
    d    = cells % dimension
    ne   = size(cells % num_cell_corners)
    kind = cells % cell_kinds(1)
    if (any(cells % cell_kinds /= kind)) then
       error stop 'operation_finite_element: the cells of a mesh are elements of one kind'
    end if

    ! the reference cell's corners, and whether it is a simplex
    select case (kind)
    case (2)
       corner  = reshape([0, 0, 1, 0, 0, 1], [2, 3]) * 1.0_dp
       simplex = .true.
    case (3)
       corner  = reshape([0, 0, 1, 0, 1, 1, 0, 1], [2, 4]) * 1.0_dp
       simplex = .false.
    case (4)
       corner  = reshape([0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 4]) * 1.0_dp
       simplex = .true.
    case (5)
       corner  = reshape([0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1], [3, 8]) * 1.0_dp
       simplex = .false.
    case default
       error stop 'operation_finite_element: the elements are triangles, quadrangles, tetrahedra or hexahedra'
    end select
    if (size(corner, 1) /= d) then
       error stop 'operation_finite_element: the cells of the mesh are of its dimension'
    end if
    nn = size(corner, 2)

    ! the polynomials: of degree one on a simplex; on a box, each power
    ! at most one, the members whose second derivative along every
    ! coordinate vanishes
    if (simplex) then
       shape  = polynomial_form(1, d)
       active = [(j, j = 1, d + 1)]
    else
       shape = polynomial_form(d, d)
       allocate(dphi(shape % num_members()))
       allocate(active(0))
       do j = 1, shape % num_members()
          s = 0
          do r = 1, d
             call shape % derivatives(spread(0.5_dp, 1, d), spread(0.0_dp, 1, d), 2 * unit_vector(r, d), dphi)
             if (dphi(j) /= 0.0_dp) s = s + 1
          end do
          if (s == 0) active = [active, j]
       end do
    end if
    if (size(active) /= nn) then
       error stop 'operation_finite_element: the polynomials of the element number its corners'
    end if

    ! the basis dual to the values at the corners
    allocate(member(shape % num_members()), a(nn, nn), unit(nn))
    do i = 1, nn
       call shape % values(corner(:, i), spread(0.0_dp, 1, d), member)
       a(i, :) = member(active)
    end do
    call factors % factorise(a, 0.0_dp)
    if (factors % singular()) then
       error stop 'operation_finite_element: the values at the corners do not determine the polynomials'
    end if

    ! the reference quadrature: Gauss-Legendre of three points on each
    ! axis of the unit box, collapsed onto the simplex
    do l = 1, 3
       call gauss_node_weight(l, 3, 1.0_dp, gauss(l), gauss_weight(l))
    end do
    nq = 3**d
    allocate(xi(d, nq), w(nq))
    q = 0
    do l = 1, 3
       do m = 1, 3
          do u = 1, merge(3, 1, d == 3)
             q = q + 1
             if (d == 2) then
                xi(:, q) = [gauss(l), gauss(m)]
                w(q)     = gauss_weight(l) * gauss_weight(m)
                if (simplex) then
                   xi(2, q) = gauss(m) * (1.0_dp - gauss(l))
                   w(q)     = w(q) * (1.0_dp - gauss(l))
                end if
             else
                xi(:, q) = [gauss(l), gauss(m), gauss(u)]
                w(q)     = gauss_weight(l) * gauss_weight(m) * gauss_weight(u)
                if (simplex) then
                   xi(2, q) = gauss(m) * (1.0_dp - gauss(l))
                   xi(3, q) = gauss(u) * (1.0_dp - gauss(l)) * (1.0_dp - gauss(m))
                   w(q)     = w(q) * (1.0_dp - gauss(l))**2 * (1.0_dp - gauss(m))
                end if
             end if
          end do
       end do
    end do

    ! the basis and its derivatives along xi at the reference points
    allocate(value(nn, nq), slope(d, nn, nq))
    do i = 1, nn
       unit    = 0.0_dp
       unit(i) = 1.0_dp
       call factors % substitute(unit, coefficients, transposed=.false.)
       do q = 1, nq
          call shape % values(xi(:, q), spread(0.0_dp, 1, d), member)
          value(i, q) = dot_product(coefficients, member(active))
          do r = 1, d
             call shape % derivatives(xi(:, q), spread(0.0_dp, 1, d), unit_vector(r, d), member)
             slope(r, i, q) = dot_product(coefficients, member(active))
          end do
       end do
    end do

    ! each cell: its corners, and at each quadrature point the weight,
    ! the position and the basis's jet along x
    set % dimension   = d
    set % order       = 0
    set % derivatives = derivatives
    allocate(set % node(nn, ne), set % weight(nq, ne), set % position(d, nq, ne))
    allocate(set % basis(0:d * derivatives, 0:0, nn, nq, ne), source=0.0_dp)
    do e = 1, ne
       set % node(:, e) = cells % cell_corners(1:nn, e)
       do q = 1, nq
          jacobian = 0.0_dp
          x        = 0.0_dp
          do i = 1, nn
             do r = 1, d
                jacobian(1:d, r) = jacobian(1:d, r) + cells % corners(1:d, set % node(i, e)) * slope(r, i, q)
             end do
             x(1:d) = x(1:d) + cells % corners(1:d, set % node(i, e)) * value(i, q)
          end do
          call inverted(jacobian, d, inverse, determinant)
          if (abs(determinant) <= epsilon(1.0_dp) * maxval(abs(jacobian(1:d, 1:d)))**d) then
             error stop 'operation_finite_element: a cell of the mesh is degenerate'
          end if
          set % weight(q, e)      = w(q) * abs(determinant)
          set % position(:, q, e) = x(1:d)
          do i = 1, nn
             set % basis(0, 0, i, q, e) = value(i, q)
             if (derivatives < 1) cycle
             do r = 1, d
                set % basis(r, 0, i, q, e) = dot_product(inverse(1:d, r), slope(1:d, i, q))
             end do
          end do
       end do
    end do
    call values_of(cells % cell_centres, centre)
    call tested(set, centre)

  end function corner_elements

  ! the unit vector of the axis r among d
  pure function unit_vector(r, d) result(e)
    integer, intent(in) :: r, d
    integer :: e(d)
    e    = 0
    e(r) = 1
  end function unit_vector

  ! the inverse of a jacobian of two or three coordinates, inverse(r, i)
  ! the derivative of xi_r along x_i, and its determinant
  pure subroutine inverted(jacobian, d, inverse, determinant)
    real(dp), intent(in)  :: jacobian(3, 3)
    integer , intent(in)  :: d
    real(dp), intent(out) :: inverse(3, 3), determinant
    real(dp) :: j(3, 3)
    j = jacobian
    inverse = 0.0_dp
    if (d == 2) then
       determinant = j(1, 1) * j(2, 2) - j(1, 2) * j(2, 1)
       inverse(1, 1) =  j(2, 2) / determinant
       inverse(1, 2) = -j(1, 2) / determinant
       inverse(2, 1) = -j(2, 1) / determinant
       inverse(2, 2) =  j(1, 1) / determinant
    else
       determinant = j(1, 1) * (j(2, 2) * j(3, 3) - j(2, 3) * j(3, 2)) &
            &      - j(1, 2) * (j(2, 1) * j(3, 3) - j(2, 3) * j(3, 1)) &
            &      + j(1, 3) * (j(2, 1) * j(3, 2) - j(2, 2) * j(3, 1))
       inverse(1, 1) = (j(2, 2) * j(3, 3) - j(2, 3) * j(3, 2)) / determinant
       inverse(1, 2) = (j(1, 3) * j(3, 2) - j(1, 2) * j(3, 3)) / determinant
       inverse(1, 3) = (j(1, 2) * j(2, 3) - j(1, 3) * j(2, 2)) / determinant
       inverse(2, 1) = (j(2, 3) * j(3, 1) - j(2, 1) * j(3, 3)) / determinant
       inverse(2, 2) = (j(1, 1) * j(3, 3) - j(1, 3) * j(3, 1)) / determinant
       inverse(2, 3) = (j(1, 3) * j(2, 1) - j(1, 1) * j(2, 3)) / determinant
       inverse(3, 1) = (j(2, 1) * j(3, 2) - j(2, 2) * j(3, 1)) / determinant
       inverse(3, 2) = (j(1, 2) * j(3, 1) - j(1, 1) * j(3, 2)) / determinant
       inverse(3, 3) = (j(1, 1) * j(2, 2) - j(1, 2) * j(2, 1)) / determinant
    end if
  end subroutine inverted

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
