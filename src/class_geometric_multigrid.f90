!=====================================================================!
! Geometric multigrid. A concrete multigrid whose entire contribution
! is its squint: which fine vertices huddle into which coarse part.
! The answer here is by coordinates - nearby cells belong together,
! and the matrix is never consulted.
!
! Each level keeps a graph and the centroids of its vertices. To
! coarsen, the level's graph is partitioned by recursive coordinate
! bisection into compact parts; the parts are the aggregates. The
! next level is then the quotient - the same kind of animal, its
! vertices the parts, its centroids the part averages - so the squint
! composes all the way down. Everything after the squint - smoothed
! prolongation, galerkin coarse operators, the V-cycle, the coarse
! LU - is the multigrid mechanism, inherited untouched.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_geometric_multigrid

  use iso_fortran_env    , only : dp => REAL64
  use interface_multigrid, only : multigrid
  use interface_graph    , only : graph
  use class_stored_graph , only : stored_graph
  use class_mesh         , only : mesh

  implicit none

  private
  public :: geometric_multigrid

  !-------------------------------------------------------------------!
  ! One rung of the geometric hierarchy: a graph and where its
  ! vertices sit
  !-------------------------------------------------------------------!

  type :: mesh_level
     type(stored_graph)    :: g
     real(dp), allocatable :: centroids(:,:)   ! (ndim, num_vertices)
  end type mesh_level

  type, extends(multigrid) :: geometric_multigrid

     type(mesh_level), allocatable :: mesh_levels(:)       ! 1=finest, grown by coarsen
     integer                       :: cells_per_part = 4   ! about this many fine cells huddle into one part

   contains

     procedure :: coarsen

  end type geometric_multigrid

  interface geometric_multigrid
     module procedure create
     module procedure create_from_mesh
  end interface geometric_multigrid

contains

  !===================================================================!
  ! Constructor: the finest rung is the given graph (its edge list
  ! copied into a stored graph) and the centroids of its vertices -
  ! for a mesh graph, the cell centres.
  !===================================================================!

  type(geometric_multigrid) function create(g, centroids) result(this)

    class(graph), intent(in) :: g
    real(dp)    , intent(in) :: centroids(:,:)

    if (.not. allocated(g % edges)) then
       error stop "geometric multigrid: the fine graph must carry an edge list"
    end if

    allocate(this % mesh_levels(this % max_levels))
    this % mesh_levels(1) % g         = stored_graph(g % num_vertices, &
         &                                     g % edges % tail, g % edges % head)
    this % mesh_levels(1) % centroids = centroids

  end function create

  !===================================================================!
  ! From a mesh alone: a mesh is a graph that knows where its
  ! vertices sit, so nothing else needs to be passed.
  !===================================================================!

  type(geometric_multigrid) function create_from_mesh(m) result(this)

    class(mesh), intent(in) :: m

    this = create(m, m % cell_centers)

  end function create_from_mesh

  !===================================================================!
  ! The squint, by coordinates. Bisect the level's graph into compact
  ! parts of about cells_per_part vertices each; the parts are the
  ! aggregates. Then build the next rung - the quotient graph with
  ! part-averaged centroids - so the level below can squint the same
  ! way.
  !===================================================================!

  subroutine coarsen(this, lev, agg, naggr)

    class(geometric_multigrid), intent(inout) :: this
    integer                   , intent(in)    :: lev
    integer, allocatable      , intent(out)   :: agg(:)
    integer                   , intent(out)   :: naggr

    integer, allocatable :: members(:)
    integer :: n, v, k, ndim

    associate(fine => this % mesh_levels(lev))

      n     = fine % g % num_vertices
      naggr = max(1, n/this % cells_per_part)

      call fine % g % partition_rcb(fine % centroids, naggr)

      allocate(agg(n))
      do v = 1, n
         agg(v) = fine % g % part_of(v)
      end do

      ! the next rung: the quotient, centred at its part averages.
      ! a second setup on the same object rebuilds the hierarchy, so a
      ! rung from the previous one may still hold centroids - drop them
      ! (the re-entrancy contract the shared setup driver promises)
      ndim = size(fine % centroids, 1)
      this % mesh_levels(lev+1) % g = stored_graph(fine % g)
      if (allocated(this % mesh_levels(lev+1) % centroids)) then
         deallocate(this % mesh_levels(lev+1) % centroids)
      end if
      allocate(this % mesh_levels(lev+1) % centroids(ndim, naggr))
      do k = 1, naggr
         members = fine % g % owned(k)
         do v = 1, ndim
            this % mesh_levels(lev+1) % centroids(v, k) = &
                 & sum(fine % centroids(v, members))/real(size(members), dp)
         end do
      end do

    end associate

  end subroutine coarsen

end module class_geometric_multigrid
