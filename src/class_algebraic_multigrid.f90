!=====================================================================!
! This module implements smoothed-aggregation algebraic multigrid
! (SA-AMG), a concrete multigrid whose entire contribution is its
! squint: which fine vertices huddle into which coarse part. The
! answer here comes from strength - the matrix says which neighbours
! matter.
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

     real(dp) :: theta = 0.08_dp   ! The strength threshold.

   contains

     procedure :: coarsen

  end type algebraic_multigrid

contains

  !===================================================================!
  ! This is the squint, by strength. The strong off-diagonals of the
  ! level's operator become the edges of a strength graph (a strong
  ! entry in either row joins the pair - the stored adjacency is
  ! symmetric and tolerates the repeat when both rows agree); the
  ! graph then discovers its own aggregates, and they come back as its
  ! partition.
  !
  !     (operator) --strength--> (strength graph) --> (aggregates)
  !===================================================================!

  subroutine coarsen(this, lev, agg, naggr)

    class(algebraic_multigrid), intent(inout) :: this
    integer                   , intent(in)    :: lev
    integer, allocatable      , intent(out)   :: agg(:)
    integer                   , intent(out)   :: naggr

    type(stored_graph)    :: strength_graph
    real(dp), allocatable :: d(:)
    integer , allocatable :: tails(:), heads(:)
    integer :: n, i

    !-----------------------------------------------------------------!
    ! Strength chooses the edges (matrix work); the inherited harvest
    ! walks the rows and collects them.
    !-----------------------------------------------------------------!

    associate(A => this % levels(lev) % A)
      n = A % num_vertices
      d = A % get_diagonal()
      call A % harvest_edges(strong_row, tails, heads, directed = .true.)
    end associate

    ! The graph does the clustering (graph work, delegated).
    strength_graph = stored_graph(n, tails, heads)
    call strength_graph % partition_aggregate()

    naggr = strength_graph % nparts
    allocate(agg(n))
    do i = 1, n
       agg(i) = strength_graph % part_of(i)
    end do

  contains

    !=================================================================!
    ! Collect the strong entries of row i. Both directions are
    ! recorded, and the stored adjacency tolerates the repeat when the
    ! two rows agree.
    !=================================================================!

    pure function strong_row(i) result(cands)
      integer, intent(in)  :: i
      integer, allocatable :: cands(:)
      integer :: e, j, m
      associate(A => this % levels(lev) % A)
        allocate(cands(A % out_xadj(i+1) - A % out_xadj(i)))
        m = 0
        do e = A % out_xadj(i), A % out_xadj(i+1) - 1
           j = A % out_adj(e)
           if (j .ne. i .and. strong(A % vals(e), d(i), d(j), this % theta)) then
              m = m + 1
              cands(m) = j
           end if
        end do
        cands = cands(1:m)
      end associate
    end function strong_row

  end subroutine coarsen

  !===================================================================!
  ! This is the strength test: an entry passes when
  ! |a_ij| >= theta sqrt(|a_ii a_jj|).
  !===================================================================!

  elemental logical function strong(aij, aii, ajj, theta)
    real(dp), intent(in) :: aij, aii, ajj, theta
    strong = abs(aij) .ge. theta*sqrt(abs(aii*ajj))
  end function strong

end module class_algebraic_multigrid
