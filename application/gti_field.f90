! THE FIELD BLOCK: one block of instants over a spatial mesh.
!
! An instant no longer holds one node's components but every node's,
! and the unknown of instant k, node i, degree d lies at
!
!      ((k - 1) nodes + (i - 1)) degrees + d + 1
!
! so the components of one point stay together, points lie node by
! node within an instant, and instants follow one another, exactly as
! a single node's block does with nodes = 1.
!
! The block is the one gti_block already solves. Its derived rows are
! the family's, replicated once per node - each node's history is its
! own - and its physics is nodal, applied at every (instant, node) as
! at every instant before. The one thing added is the level below: a
! stencil over the same unknowns carrying, on the row the physics
! sits on, minus the framework's diffusion operator over the values
! of one instant, each row divided by its cell's area so that the
! flux balance becomes kappa times the laplacian. It is linear in
! the state and knows nothing of the design, so it enters the apply
! and the tangent and nothing else. Its order is the form's, and the
! form's degree is given.
!
! The functional over a field is the integral over the domain and the
! duration, so its measure at a point is the step times the cell's
! area: one weight per point, in the order the points lie.
module gti_field

  use util_precision          , only : dp
  use operation_family        , only : family
  use operation_coupling      , only : weights_of
  use operation_weight        , only : scheme_weight
  use operation_scheme_stencil, only : derived_constraints
  use operation_stencil       , only : stencil
  use view_directed_stored    , only : stored_directed_graph
  use field_calculus          , only : field
  use field_stored            , only : stored_field
  use physics_integrand       , only : nodal_integrand
  use gti_block               , only : block_residual
  use gti_expansion           , only : block_reach
  use gti_march               , only : unknown, block_of, solved
  use gti_space               , only : room, spatial_operator, coarse_cells

  implicit none

  private
  public :: spatial_rows, field_measure
  public :: field_startup, field_aggregates

contains

  !-------------------------------------------------------------------!
  ! The aggregates multigrid coarsens a field block by: the mesh's
  ! coarse cells, every instant and degree its own, laid out as the
  ! unknowns are.
  !-------------------------------------------------------------------!

  function field_aggregates(space, n, degrees) result(aggregate)

    type(room), intent(in) :: space
    integer   , intent(in) :: n, degrees
    integer, allocatable :: aggregate(:)

    integer, allocatable :: cell(:)
    integer :: k, i, d, nodes, coarse

    cell   = coarse_cells(space)
    coarse = maxval(cell)
    nodes  = space % num_cells
    allocate(aggregate(n * nodes * degrees))

    do k = 1, n
       do i = 1, nodes
          do d = 0, degrees - 1
             aggregate(unknown(k, d, degrees, i, nodes)) = &
                  & ((k - 1) * coarse + (cell(i) - 1)) * degrees + d + 1
          end do
       end do
    end do

  end function field_aggregates

  !-------------------------------------------------------------------!
  ! The spatial rows: -kappa times the laplacian of the values of one
  ! instant, on the row of every node the physics sits on, for every
  ! instant.
  !-------------------------------------------------------------------!

  function spatial_rows(space, kappa, degree, degrees, primary, n) result(rows)

    type(room), intent(in) :: space
    real(dp)  , intent(in) :: kappa
    integer   , intent(in) :: degree, degrees, primary, n
    type(stencil) :: rows

    type(stencil) :: op
    integer , allocatable :: r(:), c(:)
    real(dp), allocatable :: lw(:), w(:), held(:)
    integer :: k, e, m, nodes, row_cell, column_cell

    nodes = space % num_cells
    op    = spatial_operator(space, kappa, degree)
    m     = op % pattern % num_edges()
    call op % weights % real_vector(lw)
    call op % constants % real_vector(held)

    ! a wall holding a value would enter as a source, which the block
    ! has no place for yet; the wall here holds no flux
    if (any(abs(held) > 0.0_dp)) then
       error stop 'gti_field: the spatial operator carries no constant'
    end if

    allocate(r(m * n), c(m * n), w(m * n))
    do k = 1, n
       do e = 1, m
          row_cell    = op % pattern % edge_head(e)
          column_cell = op % pattern % edge_tail(e)
          r((k - 1) * m + e) = unknown(k, primary, degrees, row_cell, nodes)
          c((k - 1) * m + e) = unknown(k, 0, degrees, column_cell, nodes)
          w((k - 1) * m + e) = -lw(e) / space % volume(row_cell)
       end do
    end do

    rows = stencil(r, c, w, spread(0.0_dp, 1, n * nodes * degrees), 'spatial rows')

  end function spatial_rows

  !-------------------------------------------------------------------!
  ! The measure at every point, in the order the points lie: the step
  ! ending at the instant times the cell's area.
  !-------------------------------------------------------------------!

  pure function field_measure(dt, space) result(m)

    real(dp)  , intent(in) :: dt(:)
    type(room), intent(in) :: space
    real(dp) :: m(size(dt) * space % num_cells)

    integer :: k, i

    m = [((dt(k) * space % volume(i), i = 1, space % num_cells), k = 1, size(dt))]

  end function field_measure

  !-------------------------------------------------------------------!
  ! The first h instants of a march, from a field at the first: a
  ! one-history scheme marched over the first h - 1 steps, each split
  ! r ways, and sampled at the instants. The steps given are the
  ! march's own, dt(k) ending at instant k with dt(1) zero.
  !-------------------------------------------------------------------!

  function field_startup(starter, physics, degrees, h, r, dt, space, kappa, degree, design, q0) &
       & result(held)

    class(family)         , intent(in) :: starter
    class(nodal_integrand), intent(in) :: physics
    integer               , intent(in) :: degrees, h, r, degree
    real(dp)              , intent(in) :: dt(:), kappa, design, q0(:)
    type(room)            , intent(in) :: space
    real(dp), allocatable :: held(:)

    type(block_residual) :: rows
    real(dp), allocatable :: fine(:), q(:)
    real(dp) :: achieved
    integer :: k, width, n

    width = space % num_cells * degrees

    if (h == 1) then
       held = q0
       return
    end if
    if (starter % history_depth(degrees - 1) /= 1) then
       error stop 'gti_field: a starter reads one instant back'
    end if

    n    = 1 + (h - 1) * r
    fine = [0.0_dp, (dt(1 + (k - 1) / r + 1) / real(r, dp), k = 1, (h - 1) * r)]

    rows = block_of(starter, physics, degrees, n, fine, q0, nodes=space % num_cells, &
         & spatial=spatial_rows(space, kappa, degree, degrees, starter % primary_degree(degrees - 1), n))
    call solved(rows, design, q, achieved)

    allocate(held(h * width))
    do k = 1, h
       held((k - 1) * width + 1:k * width) = &
            & q((1 + (k - 1) * r - 1) * width + 1:(1 + (k - 1) * r) * width)
    end do

  end function field_startup

end module gti_field
