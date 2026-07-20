!=====================================================================!
! Smoothed-aggregation algebraic multigrid (SA-AMG). A concrete
! multigrid whose entire contribution is its squint: which fine
! vertices huddle into which coarse part. The answer here is by
! strength - the matrix says which neighbours matter.
!
! The division of labour is strict. Strength is matrix work, done
! here: the test |a_ij| >= theta sqrt(|a_ii a_jj|) picks the edges,
! and those edges build a strength graph. The clustering is graph
! work, delegated: partition_aggregate (Vanek's three passes over the
! neighbour queries) discovers the aggregates, and they come back as
! the graph's partition. Everything after the squint - smoothed
! prolongation, galerkin coarse operators, the V-cycle, the coarse
! LU - is the multigrid mechanism, inherited untouched.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_algebraic_multigrid

  use iso_fortran_env    , only : dp => REAL64
  use interface_multigrid, only : multigrid
  use class_stored_graph , only : stored_graph

  implicit none

  private
  public :: algebraic_multigrid

  type, extends(multigrid) :: algebraic_multigrid

     real(dp) :: theta = 0.08_dp   ! strength threshold

   contains

     procedure :: coarsen

  end type algebraic_multigrid

contains

  !===================================================================!
  ! The squint, by strength. The strong off-diagonals of the level's
  ! operator become the edges of a strength graph (a strong entry in
  ! either row joins the pair - the stored adjacency is symmetric and
  ! tolerates the repeat when both rows agree); the graph then
  ! discovers its own aggregates, and they come back as its partition.
  !===================================================================!

  subroutine coarsen(this, lev, agg, naggr)

    class(algebraic_multigrid), intent(inout) :: this
    integer                   , intent(in)    :: lev
    integer, allocatable      , intent(out)   :: agg(:)
    integer                   , intent(out)   :: naggr

    type(stored_graph)    :: strength_graph
    real(dp), allocatable :: d(:)
    integer , allocatable :: tails(:), heads(:)
    integer :: n, i, k, j, nstrong, pass

    ! strength chooses the edges (matrix work): count, then collect
    associate(A => this % levels(lev) % A, theta => this % theta)

      n = A % num_vertices
      d = A % get_diagonal()

      do pass = 1, 2
         nstrong = 0
         do i = 1, n
            do k = A % out_xadj(i), A % out_xadj(i+1) - 1
               j = A % out_adj(k)
               if (j .ne. i .and. strong(A % vals(k), d(i), d(j), theta)) then
                  nstrong = nstrong + 1
                  if (pass .eq. 2) then
                     tails(nstrong) = i
                     heads(nstrong) = j
                  end if
               end if
            end do
         end do
         if (pass .eq. 1) allocate(tails(nstrong), heads(nstrong))
      end do

    end associate

    ! the graph does the clustering (graph work, delegated)
    strength_graph = stored_graph(n, tails, heads)
    call strength_graph % partition_aggregate()

    naggr = strength_graph % nparts
    allocate(agg(n))
    do i = 1, n
       agg(i) = strength_graph % part_of(i)
    end do

  end subroutine coarsen

  !===================================================================!
  ! strength test
  !===================================================================!

  elemental logical function strong(aij, aii, ajj, theta)
    real(dp), intent(in) :: aij, aii, ajj, theta
    strong = abs(aij) .ge. theta*sqrt(abs(aii*ajj))
  end function strong

end module class_algebraic_multigrid
