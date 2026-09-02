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
! The compiled path is required: attach a stencil operator. The
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
  use operation_minimization , only : minimizer, attach
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
     procedure :: attach => multigrid_attach
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
  ! After attach: read the fine stencil through the aggregates and
  ! pass the governed pair their statements. The smoother sweeps the
  ! fine one; the coarse minimizer solves the block one.
  !===================================================================!

  !===================================================================!
  ! Attach, and where aggregates are already stored, set the two levels
  ! up at once: a governing minimizer that re-attaches its inner
  ! minimizer at every iteration need not reference the coarse level.
  !===================================================================!

  subroutine multigrid_attach(this, action, on, unknown_domain, num_unknowns, &
       & num_components, coupling, stored_inputs)

    class(multigrid)     , intent(inout)        :: this
    class(operation)     , intent(in)           :: action
    class(directed_graph), intent(in)           :: on
    type(graph)          , intent(in)           :: unknown_domain
    integer              , intent(in)           :: num_unknowns
    integer              , intent(in), optional :: num_components
    class(directed_graph), intent(in), optional :: coupling
    type(stored_field)   , intent(in), optional :: stored_inputs(:)

    integer, allocatable :: kept(:)

    call attach(this, action, on, unknown_domain, num_unknowns, &
         & num_components, coupling, stored_inputs)

    if (allocated(this % aggregates)) then
       kept = this % aggregates
       call this % setup(kept)
    end if

  end subroutine multigrid_attach

  subroutine setup(this, aggregates)

    class(multigrid), intent(inout) :: this
    integer, intent(in) :: aggregates(:)

    type(stencil) :: block_statement
    integer , allocatable :: rows(:), columns(:)
    real(dp), allocatable :: weights(:), zeros(:)
    integer :: e, ne, b

    this % aggregates = aggregates
    this % nblocks    = maxval(aggregates)

    ! The Galerkin construction: every dependency is mapped to a block pair.
    select type (fine => this % action)

    type is (stencil)

       ne = fine % pattern % num_edges()
       allocate(rows(ne), columns(ne), weights(ne))

       call fine % weights % real_vector(weights)
       do e = 1, ne
          rows(e)    = aggregates(fine % pattern % edge_head(e))
          columns(e) = aggregates(fine % pattern % edge_tail(e))
       end do

       allocate(zeros(this % nblocks))
       zeros = 0.0_dp

       ! many fine dependencies map to one block pair; the
       ! coarse matrix stores their sum as one entry
       block
         integer , allocatable :: crows(:), ccolumns(:)
         real(dp), allocatable :: cweights(:)
         call combine_triples(this % nblocks, this % nblocks, &
              & rows, columns, weights, crows, ccolumns, cweights)
         block_statement = stencil(crows, ccolumns, cweights, &
              & zeros, label='block statement')
         ! versioned, so a direct coarse solver factorises it once and
         ! not once per cycle
         statements_made = statements_made + 1
         call block_statement % versioned(statements_made)
       end block

    class default

       error stop 'multigrid: attach a compiled (stencil) operator'

    end select

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
       call this % smoother % attach(this % action, this % on, &
            & this % unknown_domain, this % num_unknowns, &
            & coupling = read_through(this % action, this % block_width))
    else
       call this % smoother % attach(this % action, this % on, &
            & this % unknown_domain, this % num_unknowns, coupling = this % on)
    end if

    ! The coarse statement stores its own stencil, and that stencil
    ! is exactly the coupling of the coarse unknowns.
    call this % coarse % attach(block_statement, block_statement % pattern, &
         & block_statement % pattern % vertex_set(), &
         & block_statement % pattern % num_vertices(), &
         & coupling = block_statement % pattern)

  end subroutine setup

  !===================================================================!
  ! The cycle: smooth, restrict the residual, prolong the correction,
  ! smooth again.
  !===================================================================!

  !===================================================================!
  ! A stencil's pattern read through blocks of consecutive unknowns:
  ! the graph over the blocks with an edge where any unknown of one
  ! reads any unknown of the other, self-edges removed. The coupling
  ! a block smoother colours.
  !===================================================================!

  function read_through(action, width) result(coupling)

    class(operation), intent(in) :: action
    integer         , intent(in) :: width
    type(stored_directed_graph)  :: coupling

    integer , allocatable :: rows(:), columns(:), crows(:), ccolumns(:)
    real(dp), allocatable :: ones(:), cweights(:)
    logical , allocatable :: kept(:)
    integer :: e, ne, nb

    select type (fine => action)
    type is (stencil)
       ne = fine % pattern % num_edges()
       nb = fine % pattern % num_vertices() / width
       allocate(rows(ne), columns(ne), ones(ne))
       do e = 1, ne
          rows(e)    = (fine % pattern % edge_head(e) - 1) / width + 1
          columns(e) = (fine % pattern % edge_tail(e) - 1) / width + 1
       end do
       ones = 1.0_dp
       call combine_triples(nb, nb, rows, columns, ones, crows, ccolumns, cweights)
       kept = crows /= ccolumns
       coupling = stored_directed_graph(nb, tails=pack(ccolumns, kept), heads=pack(crows, kept))
    class default
       error stop 'multigrid: attach a compiled (stencil) operator'
    end select

  end function read_through

  subroutine solve(this, rhs, x, achieved)

    class(multigrid), intent(inout) :: this
    real(dp), intent(in)    :: rhs(:)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out)   :: achieved

    real(dp), allocatable :: r(:), rc(:), ec(:), e(:)
    real(dp) :: smoothed, coarse_residual
    integer :: it

    call tally_record(linear_solves)

    allocate(ec(this % nblocks))

    call this % begin_imbalance()

    do it = 1, this % max_iterations

       call this % smoother % solve(rhs, x, smoothed)

       call this % imbalance(rhs, x, r)

       achieved = this % norm(r)
       if (this % halted(achieved, it)) return

       ! Restriction: the residual restricted onto the blocks;
       ! prolongation: the correction computed there, prolonged by the
       ! transpose.
       call through_blocks(this % aggregates, this % nblocks, 1, r, rc, transposed=.false.)
       ec = 0.0_dp
       call this % coarse % solve(rc, ec, coarse_residual)
       call through_blocks(this % aggregates, this % nblocks, 1, e, ec, transposed=.true.)
       x = x + e

       call this % smoother % solve(rhs, x, smoothed)

    end do

    call this % imbalance(rhs, x, r)
    achieved = this % norm(r)

  end subroutine solve


end module operation_multigrid
