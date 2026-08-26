! THE FIELD: what a field adds to a block over a spatial mesh.
!
! A block over a field holds every node's components at every moment,
! node by node within a moment, laid out by the constructors in
! gti_march and gti_stage exactly as a single node's block is with
! nodes = 1. What the field adds is the level below: a stencil over
! the nodes carrying minus the framework's diffusion operator, each
! row divided by its cell's area so that the flux balance becomes
! kappa times the laplacian. A block lays it on every moment it
! evaluates its physics at, where it adds to the physics' row. It is
! linear in the state and knows nothing of the design, so it enters
! the apply and the tangent and nothing else. Its order is the form's,
! and the form's degree is given.
!
! The functional over a field is the integral over the domain and the
! duration, so its measure at a point is the step times the cell's
! area: one weight per point, in the order the points lie.
module gti_field

  use util_precision   , only : dp
  use operation_stencil, only : stencil
  use gti_space        , only : room, spatial_operator

  implicit none

  private
  public :: node_operator, field_measure

contains

  !-------------------------------------------------------------------!
  ! The level below as a stencil over the nodes: -kappa times the
  ! laplacian, one row per cell. A wall holding a value would enter
  ! as a source, which a block has no place for; the wall here holds
  ! no flux, and the block refuses a constant if one arrives.
  !-------------------------------------------------------------------!

  function node_operator(space, kappa, degree) result(op)

    type(room), intent(in) :: space
    real(dp)  , intent(in) :: kappa
    integer   , intent(in) :: degree
    type(stencil) :: op

    type(stencil) :: balance
    integer , allocatable :: r(:), c(:)
    real(dp), allocatable :: lw(:), w(:), held(:)
    integer :: e, m

    balance = spatial_operator(space, kappa, degree)
    m       = balance % pattern % num_edges()
    call balance % weights % real_vector(lw)
    call balance % constants % real_vector(held)

    allocate(r(m), c(m), w(m))
    do e = 1, m
       r(e) = balance % pattern % edge_head(e)
       c(e) = balance % pattern % edge_tail(e)
       w(e) = -lw(e) / space % volume(r(e))
    end do

    op = stencil(r, c, w, held, 'level below')

  end function node_operator

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

end module gti_field
