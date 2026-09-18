!=====================================================================!
! THE FINITE-VOLUME METHOD ON A MESH: a derivative at a cell as the
! integral of the quantity over the cell's faces divided by the
! cell's volume, by the divergence theorem,
!
!      d u / d x_a  at c   =  (1 / V_c) sum over the faces f of c of
!                             n_a(f) A_f u_f,
!
! with u_f the value at the face, the face's unit normal n outward
! from c and its area A. A first derivative reads the value at the
! face centre interpolated between the two cells the face separates,
! by the mesh's inverse-distance weights (view_mesh: face_weights,
! the tail cell's share); a second derivative along the axis reads,
! in place of u_f, the derivative along that axis at the face centre
! from the polynomial of degree one fitted over the neighbourhood of
! the two cells (operation_fitting), which on a skewed mesh is the
! derivative along the axis and not along the line of centres. Both
! are exact on the polynomials of degree two at the face centre,
! and on the rectangular box reduce to the central differences.
!
! THE ORDER. The quantities of this framework are values at the
! cell centres, and the rules of a residual read them there. The
! face integral of a reconstruction above degree one returns the
! mean of the derivative over the cell, which differs from its value
! at the centre by O(h^2) whatever the reconstruction, and a
! reconstruction of a higher degree is exact on its polynomials only
! when fitted to cell means, not to centre values. The finite-volume
! derivative of an order above two is therefore not defined on these
! quantities, and is refused with this statement: an order above two
! on centre values is what finite_difference provides. Order one has
! no derivative without a direction of bias, and is refused.
!
! A boundary face, one without a head, reads the tail cell's own
! value; a periodic face reads its head in the tail's frame through
! the neighbourhood's offsets.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_finite_volume

  use util_precision                    , only : dp
  use field_calculus                    , only : field
  use field_forms                       , only : polynomial_form
  use field_stored                      , only : stored_field
  use view_directed_stored              , only : stored_directed_graph
  use view_mesh                         , only : mesh, values_of
  use operation_action                  , only : binding
  use operation_stencil                 , only : stencil, triple_list, combine_triples
  use operation_fitting                 , only : fit
  use operation_derivative_approximation, only : derivative_approximation

  implicit none

  private
  public :: finite_volume

  type, extends(derivative_approximation) :: finite_volume

   contains

     procedure :: derivative_rows

  end type finite_volume

  interface finite_volume
     module procedure create
  end interface finite_volume

contains

  function create(order) result(this)

    integer, intent(in) :: order
    type(finite_volume) :: this

    character(len=250) :: message

    if (order /= 2) then
       write(message,'(a,i0,a)') 'operation_finite_volume: the finite-volume derivative on values at the &
            &cell centres is of order two; an order above two returns the mean of the derivative over &
            &the cell, not its value at the centre, and is what finite_difference provides; order = ', &
            & order, ' is refused'
       error stop trim(message)
    end if
    this % order = order

  end function create

  !===================================================================!
  ! The rows: for each face, the tail cell gains n_a A u_f / V_tail
  ! and the head cell loses n_a A u_f / V_head, the normal being
  ! outward from the tail; u_f is the interpolated value for a first
  ! derivative and the fitted derivative along the axis for a second.
  !===================================================================!

  function derivative_rows(this, cells, axis, derivative) result(rows)

    class(finite_volume), intent(in) :: this
    type(mesh)          , intent(in) :: cells
    integer             , intent(in) :: axis, derivative
    type(stencil) :: rows

    real(dp), allocatable :: areas(:), normals(:), fcentres(:), centres(:), volumes(:), shares(:)
    real(dp), allocatable :: w(:), pts(:), offsets(:,:)
    integer , allocatable :: members(:), orders(:), r(:), c(:), rows_of(:), columns_of(:)
    real(dp), allocatable :: weights(:), weights_of(:)
    type(triple_list) :: triples
    type(polynomial_form) :: shape
    type(fit) :: fitting
    type(stored_directed_graph) :: constellation
    type(stored_field) :: positions(1)
    class(field), allocatable :: fitted
    real(dp) :: flux
    integer :: nv, ne, d, e, t, h, j, npts
    character(len=250) :: message

    nv = cells % num_vertices()
    ne = cells % num_edges()
    d  = cells % dimension
    if (axis < 1 .or. axis > d) then
       write(message,'(a,i0,a,i0)') 'operation_finite_volume: the axis must be one of the mesh; &
            &axis = ', axis, ', dimension = ', d
       error stop trim(message)
    end if
    if (derivative < 1 .or. derivative > 2) then
       write(message,'(a,i0)') 'operation_finite_volume: the faces give the first and the second &
            &derivative; derivative = ', derivative
       error stop trim(message)
    end if

    call values_of(cells % face_area()   , areas)
    call values_of(cells % face_normal() , normals)
    call values_of(cells % face_centre() , fcentres)
    call values_of(cells % cell_centre() , centres)
    call values_of(cells % cell_volume() , volumes)
    call values_of(cells % face_weights(), shares)

    if (derivative == 2) then
       shape = polynomial_form(1, d)
       allocate(orders(d), source=0)
       orders(axis) = 1
    end if

    do e = 1, ne
       t = cells % edge_tail(e)
       h = 0
       if (cells % edge_has_head(e)) h = cells % edge_head(e)
       flux = normals(d * (e - 1) + axis) * areas(e)

       if (derivative == 1) then
          ! the interpolated face value
          if (h > 0) then
             call triples % assign(t, t,  flux * shares(e) / volumes(t))
             call triples % assign(t, h,  flux * (1.0_dp - shares(e)) / volumes(t))
             call triples % assign(h, t, -flux * shares(e) / volumes(h))
             call triples % assign(h, h, -flux * (1.0_dp - shares(e)) / volumes(h))
          else
             call triples % assign(t, t,  flux / volumes(t))
          end if
       else
          ! the derivative along the axis at the face centre, fitted
          ! over the neighbourhood of the two cells in the tail's frame
          if (h > 0) then
             call cells % neighbourhood([t, h], 0, shape % num_members(), members, offsets)
          else
             call cells % neighbourhood([t], 0, shape % num_members(), members, offsets)
          end if
          npts = size(members)
          allocate(pts(d * npts))
          do j = 1, npts
             pts(d * (j - 1) + 1:d * j) = centres(d * (members(j) - 1) + 1:d * members(j)) + offsets(:, j)
          end do
          constellation = stored_directed_graph(npts, tails=[integer ::], heads=[integer ::])
          positions(1)  = stored_field('positions', constellation % vertex_set(), npts, num_components=d)
          call positions(1) % set_real_vector(pts)
          fitting = fit(shape, at=fcentres(d * (e - 1) + 1:d * e), direction=[1.0_dp, 0.0_dp, 0.0_dp], &
               & orders=orders)
          call fitting % apply(constellation, fitting % bind(positions), fitted)
          call fitted % real_vector(w)
          do j = 1, npts
             call triples % assign(t, members(j),  flux * w(j) / volumes(t))
             if (h > 0) call triples % assign(h, members(j), -flux * w(j) / volumes(h))
          end do
          deallocate(pts)
       end if
    end do

    ! a cell read by several faces of one row appears once, its
    ! weights summed
    call triples % entries(r, c, weights)
    call combine_triples(nv, nv, r, c, weights, rows_of, columns_of, weights_of)
    rows = stencil(rows_of, columns_of, weights_of, spread(0.0_dp, 1, nv), label='finite volume')

  end function derivative_rows

end module operation_finite_volume
