!=====================================================================!
! Abstract multigrid method. Build a coarse-grid hierarchy from a fine
! operator (setup), then use it either as a preconditioner (apply: one
! cycle, z = M^-1 r - the deferred apply inherited from preconditioner)
! or as a standalone solver (solve: cycle until converged).
!
! Concrete kinds extend this: algebraic multigrid (smoothed aggregation,
! class_algebraic_multigrid) builds the hierarchy from the matrix graph;
! a geometric multigrid (later) would build it from a mesh hierarchy.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module interface_multigrid

  use iso_fortran_env        , only : dp => REAL64
  use class_csr              , only : csr_matrix
  use interface_linear_solver, only : preconditioner

  implicit none

  private
  public :: multigrid

  !-------------------------------------------------------------------!
  ! A multigrid is a preconditioner (apply = one cycle) that can also
  ! solve standalone.
  !-------------------------------------------------------------------!

  type, abstract, extends(preconditioner) :: multigrid

   contains

     procedure(setup_interface)     , deferred :: setup       ! build the hierarchy from A
     procedure(solve_interface)     , deferred :: solve       ! A x = b, cycling to tol
     procedure(num_levels_interface), deferred :: num_levels  ! levels in the hierarchy
     ! apply(this, r, z) = one cycle (z = M^-1 r) is inherited deferred
     ! from preconditioner.

  end type multigrid

  !-------------------------------------------------------------------!
  ! Deferred interfaces
  !-------------------------------------------------------------------!

  abstract interface

     ! Build the multigrid hierarchy for the fine operator A.
     subroutine setup_interface(this, A)
       import :: multigrid, csr_matrix
       class(multigrid), intent(inout) :: this
       type(csr_matrix), intent(in)    :: A
     end subroutine setup_interface

     ! Solve A x = b by cycling from the initial x until the relative
     ! residual is below max_tol or max_it cycles are spent; iters returns
     ! the number of cycles taken.
     subroutine solve_interface(this, b, x, max_tol, max_it, iters)
       import :: multigrid, dp
       class(multigrid), intent(in)    :: this
       real(dp)        , intent(in)    :: b(:)
       real(dp)        , intent(inout) :: x(:)
       real(dp)        , intent(in)    :: max_tol
       integer         , intent(in)    :: max_it
       integer         , intent(out)   :: iters
     end subroutine solve_interface

     ! Number of levels in the hierarchy (1 = a single grid).
     pure integer function num_levels_interface(this)
       import :: multigrid
       class(multigrid), intent(in) :: this
     end function num_levels_interface

  end interface

contains

end module interface_multigrid
