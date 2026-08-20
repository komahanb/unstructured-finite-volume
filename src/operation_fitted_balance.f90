!=====================================================================!
! The fitted balance: exact edge values compiled into one operator.
!
! LEVEL 3 - the statements' shared machinery. Level 2 holds
! optimizers only: things that drive a residual or govern one that
! does, and this module does neither - it wires. Its home is the
! role it serves, statement-making; its words stay level-2-neutral
! ON PURPOSE, so every statement the constitution learns to speak
! can share this one act without any physics leaking between them.
! Dependency sets only the floor (it employs the fit, so it cannot
! sink below 2); role sets the home. This module does not know what
! the scales mean or where the headless faces' values came from -
! the statements say those things in their own vocabulary. What it knows is the act: one value per edge, fitted
! exactly on the edge's neighbourhood, exchanged through incidence,
! compiled into one stencil operator.
!
! The two-point kernel measures the derivative ALONG THE CENTER
! LINE, and the balance needs it along the EDGE'S NORMAL; on a
! skewed mesh no solver can repair a value measured in the wrong
! direction. Here the assumption is never made: the geometry is
! used, not assumed. This module owns NO mathematics. Per edge it
! composes four moves:
!
!      members ····· structure     the edge's two ends and their
!                                  neighbours, and on a headless
!                                  edge its own centre point
!      positions ··· data          centres on that constellation
!      fit ········· algebra       one apply of a fit over the HANDED
!                                  form, aimed along the edge normal
!                                  at the edge centre
!      exchange ···· incidence     the value lands plus on the tail
!                                  and minus on the head, once -
!                                  what one end receives the other
!                                  gave up
!
! and assembles the answers into one stencil operator - the compiled
! form. The headless edge's own point leaves its weight through the
! constant, carrying the known value the caller stated. Higher
! accuracy is the same composition with wider rings and a richer
! form.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_fitted_balance

  use iso_fortran_env    , only : dp => REAL64
  use view_directed, only : directed_graph
  use field_calculus, only : field
  use view_directed_stored        , only : stored_directed_graph
  use field_forms        , only : form
  use field_stored  , only : stored_field
  use view_mesh   , only : mesh
  use operation_stencil, only : stencil
  use operation_fitting      , only : fit

  implicit none

  private
  public :: fitted_balance_stencil

contains

  !===================================================================!
  ! The assembly. The scales arrive per edge, already meaning
  ! whatever the caller means by them; the form says what the values
  ! are fitted over; and a headless edge's own point is stated as an
  ! AFFINE RELATION to its tail,
  !
  !      point = (1 - weight)*tail + value
  !
  ! because one number cannot say everything a boundary says. A
  ! weight of one is the pinned case - the point stands at the given
  ! value, knowing nothing of its tail - and that is the default when
  ! no weights arrive, so a caller with a value alone gets what it
  ! always got. A weight of zero is the other end: the point follows
  ! its tail, offset by the value.
  !===================================================================!

  function fitted_balance_stencil(m, shape, scales, boundary_values, &
       & boundary_weights) result(op)

    type(mesh) , intent(in) :: m
    class(form), intent(in) :: shape
    real(dp)   , intent(in) :: scales(:)
    real(dp)   , intent(in), optional :: boundary_values(:)
    real(dp)   , intent(in), optional :: boundary_weights(:)

    type(stencil) :: op

    type(stored_field) :: fa, fn, fc, fcc
    type(fit) :: fitting
    type(stored_directed_graph) :: constellation
    type(stored_field)   :: positions
    class(field), allocatable :: answer
    real(dp), allocatable :: areas(:), normals(:), fcentres(:), centres(:)
    integer , allocatable :: rows(:), columns(:), hood(:)
    real(dp), allocatable :: weights(:), pts(:), w(:), constant(:)
    real(dp) :: xf(3), vb, wb
    integer :: nv, ne, e, t, h, j, npts

    nv = m % num_vertices()
    ne = m % num_edges()

    fa = m % face_area()
    call fa % real_vector(areas)
    fn = m % face_normal()
    call fn % real_vector(normals)
    fc = m % face_centre()
    call fc % real_vector(fcentres)
    fcc = m % cell_centre()
    call fcc % real_vector(centres)

    allocate(rows(0), columns(0), weights(0))
    allocate(constant(nv))
    constant = 0.0_dp

    do e = 1, ne

       t = m % edge_tail(e)
       h = 0
       if (m % edge_has_head(e)) h = m % edge_head(e)

       ! Structure: the constellation, the headless edge's own
       ! point last.
       call neighbourhood_of(m, e, hood)
       npts = size(hood)
       if (h == 0) npts = npts + 1

       ! Data: the positions on it.
       allocate(pts(3 * npts))
       do j = 1, size(hood)
          pts(3 * j - 2 : 3 * j) = centres(3 * hood(j) - 2 : 3 * hood(j))
       end do
       xf = fcentres(3 * e - 2 : 3 * e)
       if (h == 0) pts(3 * npts - 2 : 3 * npts) = xf

       constellation = stored_directed_graph(npts, tails=[integer ::], heads=[integer ::])
       positions = stored_field('positions', constellation % vertex_set(), &
            & constellation % num_vertices(), num_components=3)
       call positions % set_real_vector(pts)

       ! Algebra: one apply, aimed along the normal at the face.
       fitting = fit(shape, at=xf, &
            & direction=normals(3 * e - 2 : 3 * e), &
            & scale=scales(e))
       call fitting % apply(constellation, [positions], answer)
       call answer % real_vector(w)

       ! The exchange: the value lands once, plus to the tail,
       ! minus to the head; the edge's own point leaves its share
       ! as a constant.
       do j = 1, size(hood)
          rows    = [rows   , t]
          columns = [columns, hood(j)]
          weights = [weights, w(j)]
          if (h > 0) then
             rows    = [rows   , h]
             columns = [columns, hood(j)]
             weights = [weights, -w(j)]
          end if
       end do

       if (h == 0) then

          vb = 0.0_dp
          wb = 1.0_dp
          if (present(boundary_values))  vb = boundary_values(e)
          if (present(boundary_weights)) wb = boundary_weights(e)

          ! The stated part leaves through the constant; the part
          ! that follows the tail stays in the matrix, on the tail's
          ! own diagonal.
          constant(t) = constant(t) + w(npts) * vb

          if (abs(1.0_dp - wb) > 0.0_dp) then
             rows    = [rows   , t]
             columns = [columns, t]
             weights = [weights, w(npts) * (1.0_dp - wb)]
          end if

       end if

       deallocate(pts)

    end do

    op = stencil(rows, columns, weights, constant, &
         & label='fitted balance')

  end function fitted_balance_stencil

  !===================================================================!
  ! The face's neighbourhood: its two cells and their neighbours,
  ! each once.
  !===================================================================!

  subroutine neighbourhood_of(m, e, hood)

    type(mesh), intent(in) :: m
    integer   , intent(in) :: e
    integer, allocatable, intent(out) :: hood(:)

    integer, allocatable :: near(:)
    integer :: t, h, j

    t = m % edge_tail(e)
    h = 0
    if (m % edge_has_head(e)) h = m % edge_head(e)

    hood = [t]
    if (h > 0) call grow(hood, h)

    call m % adjacent_vertices(t, near)
    do j = 1, size(near)
       call grow(hood, near(j))
    end do

    if (h > 0) then
       call m % adjacent_vertices(h, near)
       do j = 1, size(near)
          call grow(hood, near(j))
       end do
    end if

  end subroutine neighbourhood_of

  pure subroutine grow(hood, member)

    integer, allocatable, intent(inout) :: hood(:)
    integer, intent(in) :: member

    if (any(hood == member)) return
    hood = [hood, member]

  end subroutine grow

end module operation_fitted_balance
