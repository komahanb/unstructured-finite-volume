!=====================================================================!
! Two-grid multigrid: a minimizer that changes resolution mid-fight.
!
! One family, one story - and this member's policy is the coarse
! detour. A smoother kills the rough part of the error cheaply; the
! smooth remainder, expensive on the fine graph, is cheap one level
! down, so the residual travels to the blocks, is answered there,
! and the answer travels back:
!
!      smooth . restrict . solve coarse . prolong . correct . smooth
!
! GOVERNANCE, twice. The smoother is a held minimizer sweeping the
! fine statement; the coarse answer is a held minimizer solving the
! block statement. Multigrid does not iterate on its own - it
! schedules the two it governs.
!
! THE GALERKIN ROAD. The coarse operator is not re-derived: it is
! the fine stencil operator READ THROUGH THE AGGREGATES - each
! dependency (row, column, weight) becomes (block of row, block of
! column, weight), and the dependencies landing between one block
! pair are combined into one entry (combine_triples), so the
! coarse matrix has one entry per pair. That is R A P with summing
! restriction and injected prolongation, and it obeys the
! commutation square the suite holds:
!
!      solve_coarse( R(A(P e)) ) = e        the coarsened statement
!                                           answers as the fine one
!                                           would have, on anything
!                                           the blocks can express
!
! The compiled road is required: attach a stencil operator. The
! interpreted road would re-derive per level; that is a different
! member for another day.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_multigrid

  use iso_fortran_env    , only : dp => REAL64
  use view_directed, only : directed_graph
  use operation_stencil, only : stencil, combine_triples
  use operation_minimization , only : minimizer

  implicit none

  private
  public :: multigrid

  type, extends(minimizer) :: multigrid

     class(minimizer), allocatable :: smoother
     class(minimizer), allocatable :: coarse

     integer, allocatable :: aggregates(:)
     integer :: nblocks = 0

   contains

     procedure :: name  => multigrid_name
     procedure :: setup
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
  ! hand the governed pair their statements. The smoother sweeps the
  ! fine one; the coarse minimizer answers the block one.
  !===================================================================!

  subroutine setup(this, aggregates)

    class(multigrid), intent(inout) :: this
    integer, intent(in) :: aggregates(:)

    type(stencil) :: block_statement
    integer , allocatable :: rows(:), columns(:)
    real(dp), allocatable :: weights(:), zeros(:)
    integer :: e, ne, b

    this % aggregates = aggregates
    this % nblocks    = maxval(aggregates)

    ! The Galerkin road: every dependency lands between blocks.
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

       ! many fine dependencies land between one block pair; the
       ! coarse matrix carries their sum as one entry
       block
         integer , allocatable :: crows(:), ccolumns(:)
         real(dp), allocatable :: cweights(:)
         call combine_triples(this % nblocks, this % nblocks, &
              & rows, columns, weights, crows, ccolumns, cweights)
         block_statement = stencil(crows, ccolumns, cweights, &
              & zeros, label='block statement')
       end block

    class default

       error stop 'multigrid: attach a compiled (stencil) operator'

    end select

    ! The smoother is a STRUCTURED one - jacobi, gauss-seidel - so it
    ! is handed the dependent-variable coupling explicitly. On this
    ! path the mesh the action executes over IS the coupling of its
    ! unknowns, and saying so here makes that a caller's statement
    ! rather than the minimizer's assumption.
    call this % smoother % attach(this % action, this % on, &
         & this % unknown_domain, this % num_unknowns, coupling = this % on)

    ! The coarse statement carries its own stencil, and that stencil
    ! is exactly the coupling of the coarse unknowns.
    call this % coarse % attach(block_statement, block_statement % pattern, &
         & block_statement % pattern % vertex_set(), &
         & block_statement % pattern % num_vertices(), &
         & coupling = block_statement % pattern)

  end subroutine setup

  !===================================================================!
  ! The cycle: smooth, send the remainder down, bring the answer
  ! back, smooth again.
  !===================================================================!

  subroutine solve(this, rhs, x, achieved)

    class(multigrid), intent(inout) :: this
    real(dp), intent(in)    :: rhs(:)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out)   :: achieved

    real(dp), allocatable :: y(:), r(:), rc(:), ec(:)
    real(dp) :: goal, smoothed, answered
    integer :: it, v

    allocate(rc(this % nblocks), ec(this % nblocks))

    goal = this % tolerance * (1.0_dp + this % norm(rhs))

    do it = 1, this % max_iterations

       call this % smoother % solve(rhs, x, smoothed)

       call this % matvec(x, y)
       r = rhs - y

       achieved = this % norm(r)
       if (achieved < goal) return

       ! Down: the residual gathered onto the blocks.
       rc = 0.0_dp
       do v = 1, size(r)
          rc(this % aggregates(v)) = rc(this % aggregates(v)) + r(v)
       end do

       ! Answered there.
       ec = 0.0_dp
       call this % coarse % solve(rc, ec, answered)

       ! Back: each cell takes its block's correction.
       do v = 1, size(x)
          x(v) = x(v) + ec(this % aggregates(v))
       end do

       call this % smoother % solve(rhs, x, smoothed)

    end do

    call this % matvec(x, y)
    achieved = this % norm(rhs - y)

  end subroutine solve

end module operation_multigrid
