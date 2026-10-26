!=====================================================================!
! Module that states the fitted balance: exact edge values compiled
! into one operator.
!
! This is level 3, the statements' shared assembly. Level 2 contains
! optimizers only, operations that drive a residual or govern one that
! does, and this module does neither: it composes. Its level is set by
! the role it serves, statement assembly; its vocabulary is
! level-2-neutral by design, so that every statement defined in the
! framework can share this one assembly without any physics entering
! it. Dependency sets only the lower bound (it calls the fit, so its
! level is at least 2); role sets the level. This module does not
! store what the scales mean or where the headless faces' values were
! computed; the statements specify those in their own vocabulary. What
! it defines is the assembly: one value per edge, fitted exactly on the
! edge's neighbourhood, exchanged through incidence, compiled into one
! stencil operator.
!
! The two-point kernel measures the derivative along the centre line,
! and the balance needs it along the edge's normal; on a skewed mesh
! no solver can correct a value measured in the wrong direction. Here
! that approximation is not made: the geometry is used, not assumed.
! This module defines no mathematics of its own. For each edge it
! composes four steps:
!
!      members ····· structure     the edge's two ends and their
!                                  neighbours, and on a headless
!                                  edge its own centre point
!      positions ··· data          centres on that constellation
!      fit ········· algebra       one apply of a fit over the passed
!                                  form, directed along the edge
!                                  normal at the edge centre
!      exchange ···· incidence     the value enters with plus sign on
!                                  the tail and minus sign on the
!                                  head, once: what one end receives
!                                  the other loses
!
!          o     o     o
!            \   |   /                   an interior face e: the fit over
!      o ---- t -|-> h ---- o            t, h and their neighbours, along
!            /   |   \                   the normal n at the face centre;
!          o     o  n  o                 + on row t, - on row h
!
!      o ---- t ---x                     a boundary face: its own point x
!                                        joins the constellation
!
! and assembles the results into one stencil operator, the compiled
! form. The headless edge's own point contributes its weight through
! the constant, multiplied by the known value the caller specified.
! Higher accuracy is the same composition with wider rings and a form
! with more members.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_fitted_balance

  use util_precision    , only : dp
  use field_forms       , only : form
  use view_mesh         , only : mesh, values_of
  use operation_stencil , only : stencil, triple_list
  use operation_fitting , only : fit
  use util_factorisation, only : dense_factorisation
  use util_precision    , only : spacing_at_one

  implicit none

  private

  public :: fitted_balance_stencil
  public :: fitted_derivative_stencil

contains

  !===================================================================!
  ! Returns the assembled balance. The scales are passed per edge, with
  ! the meaning the caller assigns them; the form specifies what the
  ! values are fitted over; and a headless edge's own point is
  ! specified as an affine relation to its tail,
  !
  !      point = (1 - weight)*tail + value
  !
  ! because one number cannot specify every boundary condition. A
  ! weight of one is the fixed-value case, the point equal to the given
  ! value independently of its tail, and that is the default when no
  ! weights are passed, so a caller with a value alone obtains the
  ! previous behaviour. A weight of zero is the other extreme: the
  ! point equals its tail plus the value.
  !
  ! A headless edge whose flux is known needs no fit: its flux, times
  ! the edge's scale, enters the tail's constant directly.
  !===================================================================!

  function fitted_balance_stencil(m, shape, scales, boundary_values, &
       & boundary_weights, rings, flux_known, boundary_flux) result(op)

    type(mesh) , intent(in)           :: m
    class(form), intent(in)           :: shape
    real(dp)   , intent(in)           :: scales(:)
    real(dp)   , intent(in), optional :: boundary_values(:)
    real(dp)   , intent(in), optional :: boundary_weights(:)
    integer    , intent(in), optional :: rings
    logical    , intent(in), optional :: flux_known(:)
    real(dp)   , intent(in), optional :: boundary_flux(:)

    type(stencil) :: op

    type(fit)             :: fitting
    real(dp), allocatable :: areas(:), normals(:), fcentres(:), centres(:), offsets(:,:)
    integer , allocatable :: rows(:), columns(:), cell_neighbourhood(:)
    real(dp), allocatable :: weights(:), pts(:), w(:), constant(:)
    real(dp), allocatable :: xf(:)
    real(dp)              :: vb, wb
    type(triple_list)     :: triples
    integer               :: nv, ne, e, t, h, j, npts, width, d
    character(len=150)    :: message

    ! A fit needs at least as many points as its form has members, so
    ! the neighbourhood grows ring by ring until it contains that many,
    ! and no further: a wider neighbourhood than the fit needs makes a
    ! local fit global. A given ring count overrides.
    width = 0
    if (present(rings)) width = rings
    if (present(rings) .and. width < 1) then
       write(message,'(a,i0)') 'fitted_balance: fitted_balance_stencil requires a neighbourhood &
            &of at least one ring; rings = ', width
       error stop trim(message)
    end if

    nv = m % num_vertices()
    ne = m % num_edges()
    d = m % dimension

    allocate(xf(d))

    call values_of(m % face_area()  , areas)
    call values_of(m % face_normal(), normals)
    call values_of(m % face_centre(), fcentres)
    call values_of(m % cell_centre(), centres)

    ! The triples' capacity grows by doubling: an assembly that appends
    ! one entry at a time to an array copies the array each time, and
    ! is quadratic in the mesh
    allocate(constant(nv))
    constant = 0.0_dp

    do e = 1, ne

       t = m % edge_tail(e)
       h = 0
       if (m % edge_has_head(e)) h = m % edge_head(e)

       ! Structure: the constellation, the headless edge's own point
       ! last

       ! A headless edge whose flux is known needs no fit: the flux
       ! enters the tail's balance directly
       if (h == 0 .and. present(flux_known)) then
          if (flux_known(e)) then
             if (present(boundary_flux)) constant(t) = constant(t) + scales(e) * boundary_flux(e)
             cycle
          end if
       end if

       if (h > 0) then
          call m % neighbourhood([t, h], width, shape % num_members(), cell_neighbourhood, offsets)
       else
          call m % neighbourhood([t], width, shape % num_members(), cell_neighbourhood, offsets)
       end if
       npts = size(cell_neighbourhood)
       if (h == 0) npts = npts + 1

       ! Data: the positions on it, each in the tail's frame
       allocate(pts(d * npts))
       do j = 1, size(cell_neighbourhood)
          pts(d * j - d + 1 : d * j) = centres(d * cell_neighbourhood(j) - d + 1 : d * cell_neighbourhood(j)) &
               & + offsets(:, j)
       end do
       xf = fcentres(d * e - d + 1 : d * e)
       if (h == 0) pts(d * npts - d + 1 : d * npts) = xf

       ! Algebra: one apply, directed along the normal at the face
       fitting = fit(shape, at=xf, &
            & direction=normals(d * e - d + 1 : d * e), &
            & scale=scales(e))
       w = fitting % weights(pts)

       ! Exchange: plus on the tail's row, minus on the head's
       do j = 1, size(cell_neighbourhood)
          call triples % assign(t, cell_neighbourhood(j), w(j))
          if (h > 0) call triples % assign(h, cell_neighbourhood(j), -w(j))
       end do

       ! The headless edge's own point: (1 - weight) on the tail, the
       ! value into the constant
       if (h == 0) then
          vb = 0.0_dp
          wb = 1.0_dp
          if (present(boundary_values))  vb = boundary_values(e)
          if (present(boundary_weights)) wb = boundary_weights(e)
          constant(t) = constant(t) + w(npts) * vb
          if (abs(1.0_dp - wb) > 0.0_dp) then
             call triples % assign(t, t, w(npts) * (1.0_dp - wb))
          end if
       end if
       deallocate(pts)

    end do

    call triples % entries(rows, columns, weights)
    op = stencil(rows, columns, weights, constant, &
         & label='fitted balance')

  end function fitted_balance_stencil

  !===================================================================!
  ! Returns the derivative of a multi-index at every cell centre,
  ! fitted on the cell's neighbourhood over the passed form: one row
  ! per cell, its weights on the neighbourhood, no constant. The
  ! neighbourhood grows ring by ring until the fit has as many points
  ! as the form has members and the members are independent on them,
  ! or until it has as many rings as given:
  !
  !   ring 1           ring 2
  !     o              o o o
  !   o c o    -->   o o c o o     grown while the members are
  !     o              o o o       dependent on the points
  !===================================================================!

  function fitted_derivative_stencil(m, shape, orders, rings) result(op)

    type(mesh) , intent(in)           :: m
    class(form), intent(in)           :: shape
    integer    , intent(in)           :: orders(:)
    integer    , intent(in), optional :: rings

    type(stencil)         :: op
    type(fit)             :: fitting
    real(dp), allocatable :: centres(:), pts(:), w(:), weights(:), offsets(:,:)
    integer , allocatable :: rows(:), columns(:), cell_neighbourhood(:), active(:)
    type(triple_list)     :: triples
    integer               :: nv, c, j, npts, width, d, needed
    character(len=150)    :: message

    width = 0
    if (present(rings)) width = rings
    if (present(rings) .and. width < 1) then
       write(message,'(a,i0)') 'fitted_derivative: fitted_derivative_stencil requires a &
            &neighbourhood of at least one ring; rings = ', width
       error stop trim(message)
    end if
    nv = m % num_vertices()
    d  = m % dimension
    if (size(orders) /= d) then
       write(message,'(a,i0,a,i0)') 'fitted_derivative: fitted_derivative_stencil requires one &
            &order per coordinate; size(orders) = ', size(orders), ', dimension = ', d
       error stop trim(message)
    end if
    call values_of(m % cell_centre(), centres)

    ! The neighbourhood grows until it has as many points as the form
    ! has active members, and then until the members are independent on
    ! its points: a cubic's x y^2 and x^2 y coincide with x y on the
    ! two-ring diamond of a cartesian grid
    call shape % members(active)
    do c = 1, nv
       needed = size(active)
       do
          call m % neighbourhood([c], width, needed, cell_neighbourhood, offsets)
          npts = size(cell_neighbourhood)
          if (allocated(pts)) deallocate(pts)
          allocate(pts(d * npts))
          do j = 1, npts
             pts(d * j - d + 1 : d * j) = centres(d * cell_neighbourhood(j) - d + 1 : d * cell_neighbourhood(j)) &
                  & + offsets(:, j)
          end do
          if (independent_on(shape, centres(d * c - d + 1 : d * c), pts, active)) exit
          if (present(rings) .or. npts < needed) then
             write(message,'(a,i0,a,i0,a,i0)') 'fitted_derivative: the form members are &
                  &dependent on the neighbourhood at cell ', c, ', and it cannot grow further &
                  &(rings fixed, or npts = ', npts, ' already reached needed = ', needed, ')'
             error stop trim(message)
          end if
          needed = npts + 1
       end do
       fitting = fit(shape, at=centres(d * c - d + 1 : d * c), &
            & direction=[1.0_dp, 0.0_dp, 0.0_dp], orders=orders)
       w = fitting % weights(pts)
       do j = 1, npts
          call triples % assign(c, cell_neighbourhood(j), w(j))
       end do
       deallocate(pts)
    end do
    call triples % entries(rows, columns, weights)
    op = stencil(rows, columns, weights, spread(0.0_dp, 1, nv), label='fitted derivative')

  end function fitted_derivative_stencil

  !===================================================================!
  ! Returns whether the active members of the form are independent on
  ! the points: the Gram matrix of their values, the inactive members
  ! replaced by the identity, has no pivot below the arithmetic's
  ! spacing at one times its largest entry,
  !
  !   b(i, j) = member i at point j      (zero for an inactive member)
  !   g       = b b^T, with g(i, i) = 1 for an inactive member
  !===================================================================!

  logical function independent_on(shape, at, pts, active)

    class(form), intent(in) :: shape
    real(dp)   , intent(in) :: at(:), pts(:)
    integer    , intent(in) :: active(:)

    type(dense_factorisation) :: factor
    real(dp), allocatable     :: b(:,:), g(:,:), phi(:)
    logical , allocatable     :: is_active(:)
    integer                   :: nc, npts, d, j, i

    d    = size(at)
    npts = size(pts) / d
    nc   = shape % num_members()
    allocate(b(nc, npts), phi(nc), is_active(nc))
    is_active = .false.
    do i = 1, size(active)
       if (active(i) >= 1 .and. active(i) <= nc) is_active(active(i)) = .true.
    end do
    do j = 1, npts
       call shape % values(pts(d * j - d + 1 : d * j), at, phi)
       b(:, j) = phi
    end do
    do i = 1, nc
       if (.not. is_active(i)) b(i, :) = 0.0_dp
    end do
    g = matmul(b, transpose(b))
    do i = 1, nc
       if (.not. is_active(i)) g(i, i) = 1.0_dp
    end do
    call factor % factorise(g, spacing_at_one * maxval(abs(g)))
    independent_on = .not. factor % singular()

  end function independent_on

end module operation_fitted_balance
