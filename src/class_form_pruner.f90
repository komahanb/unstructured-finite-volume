!=====================================================================!
! The pruner: the form family's first concretion.
!
! The simplest form decision, made a citizen: a basis member the
! points cannot see - one whose column of values vanishes on the
! whole constellation - is struck from the roster before any fit
! runs. What used to hide inside a solve as a pivot trick is now a
! form change, owned by the sector that changes forms, at the slow
! cadence form changes keep.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_form_pruner

  use iso_fortran_env    , only : dp => REAL64
  use graph_calculus     , only : GRAPH_SIDE_VERTEX
  use class_graph_support, only : support
  use graph_forms        , only : form
  use graph_fitting      , only : form_optimizer

  implicit none

  private
  public :: pruner

  type, extends(form_optimizer) :: pruner

     real(dp) :: threshold = 1.0d-12

   contains

     procedure :: adapt

  end type pruner

contains

  !===================================================================!
  ! Strike what the points cannot see. The constant member always
  ! survives: every point sees it.
  !===================================================================!

  subroutine adapt(this, shape, positions)

    class(pruner), intent(in) :: this
    class(form), intent(inout) :: shape
    real(dp), intent(in) :: positions(:)

    real(dp), allocatable :: phi(:), seen(:)
    real(dp) :: middle(3)
    integer :: nc, npts, j, m

    nc   = shape % size_of()
    npts = size(positions) / 3

    ! The members are read about the constellation's own middle, so
    ! what the points cannot see does not depend on where they sit.
    middle = 0.0_dp
    do j = 1, npts
       middle = middle + positions(3 * j - 2 : 3 * j)
    end do
    middle = middle / real(max(npts, 1), dp)

    allocate(phi(nc), seen(nc))
    seen = 0.0_dp

    do j = 1, npts
       call shape % values(positions(3 * j - 2 : 3 * j), middle, phi)
       do m = 1, nc
          seen(m) = seen(m) + phi(m) * phi(m)
       end do
    end do

    ! Membership is the roster: the pruned form's member set is the
    ! kept table entries. A set needs no second list to say who
    ! belongs.
    shape % support = support(GRAPH_SIDE_VERTEX, &
         & pack([(m, m = 1, nc)], seen > this % threshold))

  end subroutine adapt

end module class_form_pruner
