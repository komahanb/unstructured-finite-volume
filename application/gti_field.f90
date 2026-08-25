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
  use gti_march               , only : solved
  use gti_space               , only : room, spatial_operator, coarse_cells

  implicit none

  private
  public :: field_unknown, field_block_of, field_measure
  public :: consistent_field, field_startup, field_aggregates

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
             aggregate(field_unknown(k, i, d, nodes, degrees)) = &
                  & ((k - 1) * coarse + (cell(i) - 1)) * degrees + d + 1
          end do
       end do
    end do

  end function field_aggregates

  pure integer function field_unknown(instant, node, degree, nodes, degrees) result(at)

    integer, intent(in) :: instant, node, degree, nodes, degrees

    at = ((instant - 1) * nodes + (node - 1)) * degrees + degree + 1

  end function field_unknown

  !-------------------------------------------------------------------!
  ! The family's derived rows over n instants, once per node.
  !-------------------------------------------------------------------!

  function time_rows(scheme, degrees, n, dt, nodes) result(rows)

    class(family), intent(in) :: scheme
    integer      , intent(in) :: degrees, n, nodes
    real(dp)     , intent(in) :: dt(:)
    type(stencil) :: rows

    integer , allocatable :: tails(:), heads(:), source_degree(:), determines(:)
    integer , allocatable :: determined(:), source(:)
    real(dp), allocatable :: w(:), replicated(:)
    integer :: e, i, ne

    call block_reach(scheme, degrees, n, tails, heads, source_degree, determines)
    ne = size(tails)
    call weights_of(scheme_weight(scheme), n, tails, heads, dt, source_degree, determines, w)

    allocate(determined(ne * nodes), source(ne * nodes), replicated(ne * nodes))
    do i = 1, nodes
       do e = 1, ne
          determined((i - 1) * ne + e) = field_unknown(heads(e), i, determines(e), nodes, degrees)
          source((i - 1) * ne + e)     = field_unknown(tails(e), i, source_degree(e), nodes, degrees)
          replicated((i - 1) * ne + e) = w(e)
       end do
    end do

    rows = derived_constraints(determined, source, replicated, n * nodes * degrees, &
         & 'derived rows')

  end function time_rows

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
          r((k - 1) * m + e) = field_unknown(k, row_cell, primary, nodes, degrees)
          c((k - 1) * m + e) = field_unknown(k, column_cell, 0, nodes, degrees)
          w((k - 1) * m + e) = -lw(e) / space % volume(row_cell)
       end do
    end do

    rows = stencil(r, c, w, spread(0.0_dp, 1, n * nodes * degrees), 'spatial rows')

  end function spatial_rows

  !-------------------------------------------------------------------!
  ! The block over n instants: held holds the first history_depth
  ! instants in unknown order. A held vector of the wrong extent stops
  ! the program.
  !-------------------------------------------------------------------!

  function field_block_of(scheme, physics, degrees, n, dt, space, kappa, degree, held) &
       & result(rows)

    class(family)         , intent(in) :: scheme
    class(nodal_integrand), intent(in) :: physics
    integer               , intent(in) :: degrees, n, degree
    real(dp)              , intent(in) :: dt(:), kappa, held(:)
    type(room)            , intent(in) :: space
    type(block_residual) :: rows

    integer, allocatable :: carried(:), at(:)
    integer :: h, k, i, d, nodes, primary

    nodes   = space % num_cells
    h       = scheme % history_depth(degrees - 1)
    primary = scheme % primary_degree(degrees - 1)

    carried = [(((field_unknown(k, i, d, nodes, degrees), d = 0, degrees - 1), &
         &        i = 1, nodes), k = 1, h)]
    at      = [((field_unknown(k, i, 0, nodes, degrees) - 1, i = 1, nodes), k = 1, n)]

    if (size(held) /= size(carried)) then
       error stop 'gti_field: one value per carried component'
    end if

    rows = block_residual(time_rows(scheme, degrees, n, dt, nodes), physics, at, &
         & n * nodes * degrees, degrees, primary, carried, held, &
         & spatial=spatial_rows(space, kappa, degree, degrees, primary, n))

  end function field_block_of

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
  ! The consistent field at one instant: the components below the
  ! highest given at every node, lower(d + 1, i), and the highest at
  ! every node what the physics and the laplacian say it is there.
  ! The same one-instant block as consistent_state, with the level
  ! below attached.
  !-------------------------------------------------------------------!

  function consistent_field(physics, degrees, lower, space, kappa, degree, design) result(q)

    class(nodal_integrand), intent(in) :: physics
    integer               , intent(in) :: degrees, degree
    real(dp)              , intent(in) :: lower(:,:), kappa, design
    type(room)            , intent(in) :: space
    real(dp), allocatable :: q(:)

    type(block_residual) :: rows
    type(stencil) :: none
    integer , allocatable :: carried(:), at(:)
    real(dp), allocatable :: held(:)
    real(dp) :: achieved
    integer :: i, d, nodes

    nodes = space % num_cells

    if (size(lower, 1) /= degrees - 1 .or. size(lower, 2) /= nodes) then
       error stop 'gti_field: the components below the highest are given at every node'
    end if

    none = stencil([integer ::], [integer ::], [real(dp) ::], &
         & spread(0.0_dp, 1, nodes * degrees), 'none')

    carried = [((field_unknown(1, i, d, nodes, degrees), d = 0, degrees - 2), i = 1, nodes)]
    held    = [((lower(d + 1, i), d = 0, degrees - 2), i = 1, nodes)]
    at      = [(field_unknown(1, i, 0, nodes, degrees) - 1, i = 1, nodes)]

    rows = block_residual(none, physics, at, nodes * degrees, degrees, degrees - 1, &
         & carried, held, spatial=spatial_rows(space, kappa, degree, degrees, degrees - 1, 1))

    call solved(rows, design, q, achieved)

  end function consistent_field

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

    rows = field_block_of(starter, physics, degrees, n, fine, space, kappa, degree, q0)
    call solved(rows, design, q, achieved)

    allocate(held(h * width))
    do k = 1, h
       held((k - 1) * width + 1:k * width) = &
            & q((1 + (k - 1) * r - 1) * width + 1:(1 + (k - 1) * r) * width)
    end do

  end function field_startup

end module gti_field
