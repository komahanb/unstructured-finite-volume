!=====================================================================!
! Two-grid multigrid: a minimizer that changes resolution within one
! iteration.
!
! One family, one iteration structure - and this member's policy is
! the coarse correction. A smoother removes the high-frequency part
! of the error at low cost; the smooth remainder, high cost on the
! fine graph, is low cost one level down, so the residual is
! restricted to the blocks, solved there, and the correction is
! prolonged back:
!
!      smooth . restrict . solve coarse . prolong . correct . smooth
!
! GOVERNANCE, twice. The smoother is a stored minimizer sweeping the
! fine statement; the coarse correction is a stored minimizer solving
! the block statement. Multigrid does not iterate on its own - it
! schedules the two minimizers it governs.
!
! THE GALERKIN CONSTRUCTION. The coarse operator is not re-derived: it is
! the fine stencil operator READ THROUGH THE AGGREGATES - each
! dependency (row, column, weight) becomes (block of row, block of
! column, weight), and the dependencies mapped to one block
! pair are combined into one entry (combine_triples), so the
! coarse matrix has one entry per pair. That is R A P with summing
! restriction and injected prolongation, and it satisfies the
! commutation identity the test suite checks:
!
!      solve_coarse( R(A(P e)) ) = e        the coarsened statement
!                                           returns what the fine one
!                                           would, on any vector the
!                                           blocks can express
!
! The explicit path is required: state a stencil operator. The
! interpreted path would re-derive per level; that is a different
! member, not implemented.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_multigrid

  use util_precision  , only : dp
  use view_directed, only : directed_graph
  use view_directed_stored, only : stored_directed_graph
  use operation_stencil, only : stencil, combine_triples
  use operation_minimization , only : minimizer, state, solve_result, SOLVE_INNER_FAILED
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use operation_action       , only : operation
  use graph_fractal          , only : graph
  use field_stored           , only : stored_field
  use util_tally, only : tally_record, linear_solves
  use transform_structure, only : through_blocks

  implicit none

  private

  ! a version number for every coarse statement made
  integer, save :: statements_made = 0
  public :: multigrid

  type, extends(minimizer) :: multigrid

     class(minimizer), allocatable :: smoother
     class(minimizer), allocatable :: coarse

     integer, allocatable :: aggregates(:)
     integer :: nblocks = 0

   contains

     procedure :: name  => multigrid_name
     procedure :: setup
     procedure :: state => multigrid_state
     procedure :: solve

  end type multigrid

contains

  pure function multigrid_name(this) result(name)

    class(multigrid), intent(in) :: this
    character(len=:), allocatable :: name

    associate (u1 => this); end associate

    name = 'two-grid multigrid'

  end function multigrid_name

  !===================================================================!
  ! After state: read the fine stencil through the aggregates and
  ! pass the governed pair their statements. The smoother sweeps the
  ! fine one; the coarse minimizer solves the block one.
  !===================================================================!

  !===================================================================!
  ! State, and where aggregates are already stored, set the two levels
  ! up at once: a governing minimizer that re-states its inner
  ! minimizer at every iteration need not reference the coarse level.
  !===================================================================!

  subroutine multigrid_state(this, action, context, unknown_domain, num_unknowns, &
       & num_components, coupling, stored_inputs)

    class(multigrid)     , intent(inout)        :: this
    class(operation)     , intent(in)           :: action
    class(directed_graph), intent(in)           :: context
    type(graph)          , intent(in)           :: unknown_domain
    integer              , intent(in)           :: num_unknowns
    integer              , intent(in), optional :: num_components
    class(directed_graph), intent(in), optional :: coupling
    type(stored_field)   , intent(in), optional :: stored_inputs(:)

    integer, allocatable :: aggregate_labels(:)

    call state(this, action, context, unknown_domain, num_unknowns, &
         & num_components, coupling, stored_inputs)

    if (allocated(this % aggregates)) then
       aggregate_labels = this % aggregates
       call this % setup(aggregate_labels)
    end if

  end subroutine multigrid_state

  subroutine setup(this, aggregates)

    class(multigrid), intent(inout) :: this
    integer, intent(in) :: aggregates(:)

    type(stencil) :: block_statement
    integer , allocatable :: rows(:), columns(:)
    real(dp), allocatable :: weights(:), zeros(:)

    this % aggregates = aggregates
    this % nblocks    = maxval(aggregates)

    ! The Galerkin construction: every dependency mapped to a block
    ! pair, the dependencies of one pair summed into one entry.
    call through_map(this % action, aggregates, this % nblocks, rows, columns, weights)

    allocate(zeros(this % nblocks))
    zeros = 0.0_dp

    block_statement = stencil(rows, columns, weights, zeros, label='block statement')
    ! versioned, so a direct coarse solver factorises it once and
    ! not once per cycle
    statements_made = statements_made + 1
    call block_statement % versioned(statements_made)

    ! The smoother is a STRUCTURED one - jacobi, gauss-seidel - so it
    ! is passed the dependent-variable coupling explicitly. On this
    ! path the mesh the action executes over IS the coupling of its
    ! unknowns, and stating so here makes that a caller's statement
    ! rather than the minimizer's assumption.
    ! the smoother sweeps a block at a time where the unknowns come in
    ! blocks, and colours the coupling between blocks: the fine
    ! pattern read through the blocks, one vertex each
    this % smoother % block_width = this % block_width
    if (this % block_width > 1) then
       call this % smoother % state(this % action, this % graph, &
            & this % unknown_domain, this % num_unknowns, &
            & coupling = read_through(this % action, size(this % affine), this % block_width))
    else
       call this % smoother % state(this % action, this % graph, &
            & this % unknown_domain, this % num_unknowns, coupling = this % graph)
    end if

    ! The coarse statement stores its own stencil, and that stencil
    ! is exactly the coupling of the coarse unknowns.
    call this % coarse % state(block_statement, block_statement % pattern, &
         & block_statement % pattern % vertex_set(), &
         & block_statement % pattern % num_vertices(), &
         & coupling = block_statement % pattern)

  end subroutine setup

  !===================================================================!
  ! The cycle: smooth, restrict the residual, prolong the correction,
  ! smooth again.
  !===================================================================!

  !===================================================================!
  ! A stencil's dependencies read through a map of its vertices onto
  ! nb blocks: (row, column, weight) becomes (block of row, block of
  ! column, weight), and the dependencies of one block pair are
  ! summed into one entry, in the order combine_triples defines. An
  ! action that is not a stencil stops the program.
  !===================================================================!

  subroutine through_map(fine, block_of, nb, rows, columns, weights)

    class(operation), intent(in) :: fine
    integer         , intent(in) :: block_of(:), nb
    integer , allocatable, intent(out) :: rows(:), columns(:)
    real(dp), allocatable, intent(out) :: weights(:)

    integer , allocatable :: r(:), c(:)
    real(dp), allocatable :: w(:)
    integer :: e, ne

    select type (fine)
    type is (stencil)
       ne = fine % pattern % num_edges()
       allocate(r(ne), c(ne))
       call fine % weights % real_vector(w)
       do e = 1, ne
          r(e) = block_of(fine % pattern % edge_head(e))
          c(e) = block_of(fine % pattern % edge_tail(e))
       end do
       call combine_triples(nb, nb, r, c, w, rows, columns, weights)
    class default
       error stop 'multigrid: state an explicit (stencil) operator'
    end select

  end subroutine through_map

  !===================================================================!
  ! A stencil's pattern read through blocks of consecutive unknowns:
  ! the graph over the blocks with an edge where any unknown of one
  ! reads any unknown of the other, self-edges removed. The coupling
  ! a block smoother colours. num_vertices is the stencil's vertex
  ! count: the extent of the affine part evaluated from it when the
  ! minimizer was stated.
  !===================================================================!

  function read_through(action, num_vertices, width) result(coupling)

    class(operation), intent(in) :: action
    integer         , intent(in) :: num_vertices, width
    type(stored_directed_graph)  :: coupling

    integer , allocatable :: crows(:), ccolumns(:)
    real(dp), allocatable :: cweights(:)
    logical , allocatable :: off_diagonal(:)
    integer :: v, nb

    nb = num_vertices / width
    call through_map(action, [((v - 1) / width + 1, v = 1, num_vertices)], nb, &
         & crows, ccolumns, cweights)
    off_diagonal = crows /= ccolumns
    coupling = stored_directed_graph(nb, tails=pack(ccolumns, off_diagonal), heads=pack(crows, off_diagonal))

  end function read_through

  subroutine solve(this, rhs, x, achieved)

    class(multigrid), intent(inout) :: this
    real(dp), intent(in)    :: rhs(:)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out)   :: achieved

    real(dp), allocatable :: r(:), rc(:), ec(:), e(:)
    real(dp) :: smoothed, coarse_residual
    integer :: it
    type(solve_result) :: outcome

    call tally_record(linear_solves)

    allocate(ec(this % nblocks))

    call this % initialize_residual_history()
    call this % imbalance(rhs, x, r)
    achieved = this % norm(r)
    if (this % terminated(achieved, 0)) return

    do it = 1, this % max_iterations

       call this % smoother % solve(rhs, x, smoothed)

       call this % imbalance(rhs, x, r)
       outcome = this % smoother % result()
       if (outcome % failed() .or. .not. ieee_is_finite(smoothed)) then
          achieved = this % norm(r)
          call this % record_result(achieved, it - 1, SOLVE_INNER_FAILED)
          return
       end if

       ! Restriction: the residual restricted onto the blocks;
       ! prolongation: the correction computed there, prolonged by the
       ! transpose.
       call through_blocks(this % aggregates, this % nblocks, 1, r, rc, transposed=.false.)
       ec = 0.0_dp
       call this % coarse % solve(rc, ec, coarse_residual)
       outcome = this % coarse % result()
       if (outcome % failed() .or. .not. ieee_is_finite(coarse_residual)) then
          achieved = this % norm(r)
          call this % record_result(achieved, it - 1, SOLVE_INNER_FAILED)
          return
       end if
       call through_blocks(this % aggregates, this % nblocks, 1, e, ec, transposed=.true.)
       x = x + e

       call this % smoother % solve(rhs, x, smoothed)

       call this % imbalance(rhs, x, r)
       achieved = this % norm(r)
       outcome = this % smoother % result()
       if (outcome % failed() .or. .not. ieee_is_finite(smoothed)) then
          call this % record_result(achieved, it - 1, SOLVE_INNER_FAILED)
          return
       end if
       if (this % terminated(achieved, it)) return

    end do

  end subroutine solve


end module operation_multigrid
