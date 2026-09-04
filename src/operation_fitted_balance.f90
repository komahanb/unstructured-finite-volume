!=====================================================================!
! The fitted balance: exact edge values compiled into one operator.
!
! LEVEL 3 - the statements' shared assembly. Level 2 contains
! optimizers only: operations that drive a residual or govern one
! that does, and this module does neither - it composes. Its level
! is set by the role it serves, statement assembly; its vocabulary
! is level-2-neutral BY DESIGN, so every statement defined in the
! framework can share this one assembly without any physics
! entering it. Dependency sets only the lower bound (it calls the
! fit, so its level is at least 2); role sets the level. This module
! does not store what the scales mean or where the headless faces'
! values were computed - the statements specify those in their own
! vocabulary. What it defines is the assembly: one value per edge,
! fitted exactly on the edge's neighbourhood, exchanged through
! incidence, compiled into one stencil operator.
!
! The two-point kernel measures the derivative ALONG THE CENTRE
! LINE, and the balance needs it along the EDGE'S NORMAL; on a
! skewed mesh no solver can correct a value measured in the wrong
! direction. Here that approximation is not made: the geometry is
! used, not assumed. This module defines NO mathematics of its own.
! Per edge it composes four steps:
!
!      members ····· structure     the edge's two ends and their
!                                  neighbours, and on a headless
!                                  edge its own centre point
!      positions ··· data          centres on that constellation
!      fit ········· algebra       one apply of a fit over the PASSED
!                                  form, directed along the edge
!                                  normal at the edge centre
!      exchange ···· incidence     the value enters with plus sign on
!                                  the tail and minus sign on the
!                                  head, once - what one end receives
!                                  the other loses
!
! and assembles the results into one stencil operator - the compiled
! form. The headless edge's own point contributes its weight through
! the constant, multiplied by the known value the caller specified.
! Higher accuracy is the same composition with wider rings and a
! form with more members.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_fitted_balance

  use util_precision  , only : dp
  use view_directed, only : directed_graph
  use field_calculus, only : field
  use view_directed_stored        , only : stored_directed_graph
  use field_forms        , only : form
  use field_stored  , only : stored_field
  use view_mesh   , only : mesh, values_of
  use operation_stencil, only : stencil, triple_list
  use operation_fitting      , only : fit

  implicit none

  private
  public :: fitted_balance_stencil
  public :: fitted_derivative_stencil

contains

  !===================================================================!
  ! The assembly. The scales are passed per edge, with the meaning
  ! the caller assigns them; the form specifies what the values are
  ! fitted over; and a headless edge's own point is specified as an
  ! AFFINE RELATION to its tail,
  !
  !      point = (1 - weight)*tail + value
  !
  ! because one number cannot specify every boundary condition. A
  ! weight of one is the fixed-value case - the point equals the
  ! given value, independent of its tail - and that is the default
  ! when no weights are passed, so a caller with a value alone
  ! obtains the previous behaviour. A weight of zero is the other
  ! extreme: the point equals its tail plus the value.
  !===================================================================!

  function fitted_balance_stencil(m, shape, scales, boundary_values, &
       & boundary_weights, rings, flux_known, boundary_flux) result(op)

    type(mesh) , intent(in) :: m
    class(form), intent(in) :: shape
    real(dp)   , intent(in) :: scales(:)
    real(dp)   , intent(in), optional :: boundary_values(:)
    real(dp)   , intent(in), optional :: boundary_weights(:)
    integer    , intent(in), optional :: rings
    logical    , intent(in), optional :: flux_known(:)
    real(dp)   , intent(in), optional :: boundary_flux(:)

    type(stencil) :: op

    type(fit) :: fitting
    type(stored_directed_graph) :: constellation
    type(stored_field)   :: positions
    class(field), allocatable :: fitted
    real(dp), allocatable :: areas(:), normals(:), fcentres(:), centres(:)
    integer , allocatable :: rows(:), columns(:), cell_neighbourhood(:)
    real(dp), allocatable :: weights(:), pts(:), w(:), constant(:)
    real(dp), allocatable :: xf(:)
    real(dp) :: vb, wb
    type(triple_list) :: triples
    integer :: nv, ne, e, t, h, j, npts, width, d

    ! A fit needs at least as many points as its form has members,
    ! so the neighbourhood grows ring by ring until it contains that
    ! many, and no further: a wider neighbourhood than the fit needs
    ! makes a local fit global. A given ring count overrides.
    width = 0
    if (present(rings)) width = rings
    if (present(rings) .and. width < 1) then
       error stop 'fitted_balance: a neighbourhood is at least one ring'
    end if

    nv = m % num_vertices()
    ne = m % num_edges()
    d = m % dimension

    allocate(xf(d))

    call values_of(m % face_area()  , areas)
    call values_of(m % face_normal(), normals)
    call values_of(m % face_centre(), fcentres)
    call values_of(m % cell_centre(), centres)

    ! capacity for the triples grows by doubling: an assembly that
    ! appends one entry at a time to an array copies the array each
    ! time, and is quadratic in the mesh
    allocate(constant(nv))
    constant = 0.0_dp

    do e = 1, ne

       t = m % edge_tail(e)
       h = 0
       if (m % edge_has_head(e)) h = m % edge_head(e)

       ! Structure: the constellation, the headless edge's own
       ! point last.

       ! a headless edge whose flux is known needs no fit: the flux
       ! enters the tail's balance directly
       if (h == 0 .and. present(flux_known)) then
          if (flux_known(e)) then
             if (present(boundary_flux)) constant(t) = constant(t) + scales(e) * boundary_flux(e)
             cycle
          end if
       end if

       call neighbourhood_of(m, e, cell_neighbourhood, width, shape % num_members())
       npts = size(cell_neighbourhood)
       if (h == 0) npts = npts + 1

       ! Data: the positions on it.
       allocate(pts(d * npts))
       do j = 1, size(cell_neighbourhood)
          pts(d * j - d + 1 : d * j) = centres(d * cell_neighbourhood(j) - d + 1 : d * cell_neighbourhood(j))
       end do
       xf = fcentres(d * e - d + 1 : d * e)
       if (h == 0) pts(d * npts - d + 1 : d * npts) = xf

       constellation = stored_directed_graph(npts, tails=[integer ::], heads=[integer ::])
       positions = stored_field('positions', constellation % vertex_set(), &
            & constellation % num_vertices(), num_components=d)
       call positions % set_real_vector(pts)

       ! Algebra: one apply, directed along the normal at the face.
       fitting = fit(shape, at=xf, &
            & direction=normals(d * e - d + 1 : d * e), &
            & scale=scales(e))
       call fitting % apply(constellation, fitting % bind([positions]), fitted)
       call fitted % real_vector(w)
       do j = 1, size(cell_neighbourhood)
          call triples % assign(t, cell_neighbourhood(j), w(j))
          if (h > 0) call triples % assign(h, cell_neighbourhood(j), -w(j))
       end do
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
  ! The derivative of a multi-index at every cell centre, fitted on
  ! the cell's neighbourhood over the passed form: one row per cell,
  ! its weights on the neighbourhood, no constant. The neighbourhood
  ! grows ring by ring until the fit has as many points as the form
  ! has members, or as many rings as given.
  !===================================================================!

  function fitted_derivative_stencil(m, shape, orders, rings) result(op)

    type(mesh) , intent(in) :: m
    class(form), intent(in) :: shape
    integer    , intent(in) :: orders(:)
    integer    , intent(in), optional :: rings

    type(stencil) :: op
    type(fit) :: fitting
    type(stored_directed_graph) :: constellation
    type(stored_field)   :: positions
    class(field), allocatable :: fitted
    real(dp), allocatable :: centres(:), pts(:), w(:), weights(:)
    integer , allocatable :: rows(:), columns(:), cell_neighbourhood(:)
    type(triple_list) :: triples
    integer :: nv, c, j, npts, width, d

    width = 0
    if (present(rings)) width = rings
    if (present(rings) .and. width < 1) then
       error stop 'fitted_derivative: a neighbourhood is at least one ring'
    end if
    nv = m % num_vertices()
    d  = m % dimension
    if (size(orders) /= d) then
       error stop 'fitted_derivative: one order per coordinate of the mesh'
    end if
    call values_of(m % cell_centre(), centres)
    do c = 1, nv
       call cell_neighbourhood_of(m, c, cell_neighbourhood, width, shape % num_members())
       npts = size(cell_neighbourhood)
       allocate(pts(d * npts))
       do j = 1, npts
          pts(d * j - d + 1 : d * j) = centres(d * cell_neighbourhood(j) - d + 1 : d * cell_neighbourhood(j))
       end do
       constellation = stored_directed_graph(npts, tails=[integer ::], heads=[integer ::])
       positions = stored_field('positions', constellation % vertex_set(), &
            & constellation % num_vertices(), num_components=d)
       call positions % set_real_vector(pts)
       fitting = fit(shape, at=centres(d * c - d + 1 : d * c), &
            & direction=[1.0_dp, 0.0_dp, 0.0_dp], orders=orders)
       call fitting % apply(constellation, fitting % bind([positions]), fitted)
       call fitted % real_vector(w)
       do j = 1, npts
          call triples % assign(c, cell_neighbourhood(j), w(j))
       end do
       deallocate(pts)
    end do
    call triples % entries(rows, columns, weights)
    op = stencil(rows, columns, weights, spread(0.0_dp, 1, nv), label='fitted derivative')

  end function fitted_derivative_stencil

  !===================================================================!
  ! A cell's neighbourhood: itself and its neighbours ring by ring,
  ! each once, as many rings as given or until at least the count
  ! required is reached.
  !===================================================================!

  subroutine cell_neighbourhood_of(m, c, cell_neighbourhood, rings, at_least)

    type(mesh), intent(in) :: m
    integer   , intent(in) :: c, rings, at_least
    integer, allocatable, intent(out) :: cell_neighbourhood(:)

    integer, allocatable :: near(:), frontier(:)
    integer :: j, r, k, before

    cell_neighbourhood = [c]
    r = 0
    do
       if (rings > 0) then
          if (r >= rings) exit
       else
          if (r >= 1 .and. size(cell_neighbourhood) >= at_least) exit
       end if
       before   = size(cell_neighbourhood)
       frontier = cell_neighbourhood
       do k = 1, size(frontier)
          call m % adjacent_vertices(frontier(k), near)
          do j = 1, size(near)
             call extend(cell_neighbourhood, near(j))
          end do
       end do
       r = r + 1
       if (size(cell_neighbourhood) == before) exit
    end do

  end subroutine cell_neighbourhood_of

  !===================================================================!
  ! The face's neighbourhood: its two cells and their neighbours,
  ! each once.
  !===================================================================!

  subroutine neighbourhood_of(m, e, cell_neighbourhood, rings, at_least)

    type(mesh), intent(in) :: m
    integer   , intent(in) :: e, rings, at_least
    integer, allocatable, intent(out) :: cell_neighbourhood(:)

    integer, allocatable :: near(:), frontier(:)
    integer :: t, h, j, r, k, before

    t = m % edge_tail(e)
    h = 0
    if (m % edge_has_head(e)) h = m % edge_head(e)

    cell_neighbourhood = [t]
    if (h > 0) call extend(cell_neighbourhood, h)

    ! each ring adds the neighbours of every member of the previous
    ! ring: as many rings as given, or until the neighbourhood
    ! contains at least the count required, or until no new member
    ! is found
    r = 0
    do
       if (rings > 0) then
          if (r >= rings) exit
       else
          if (r >= 1 .and. size(cell_neighbourhood) >= at_least) exit
       end if
       before   = size(cell_neighbourhood)
       frontier = cell_neighbourhood
       do k = 1, size(frontier)
          call m % adjacent_vertices(frontier(k), near)
          do j = 1, size(near)
             call extend(cell_neighbourhood, near(j))
          end do
       end do
       r = r + 1
       if (size(cell_neighbourhood) == before) exit
    end do

  end subroutine neighbourhood_of

  pure subroutine extend(cell_neighbourhood, member)

    integer, allocatable, intent(inout) :: cell_neighbourhood(:)
    integer, intent(in) :: member

    if (any(cell_neighbourhood == member)) return
    cell_neighbourhood = [cell_neighbourhood, member]

  end subroutine extend

end module operation_fitted_balance
