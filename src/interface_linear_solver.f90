!=====================================================================!
! Interface for linear solvers, plus the preconditioner contract they
! consume. A linear solver solves A x = b for x; a preconditioner applies
! an approximate inverse z = M^-1 r once per iteration. Both abstract
! types live here so a solver and its preconditioner share one home.
!
! linear_solver extends the common algebraic_solver base.
!
! Author : Komahan Boopathy
!=====================================================================!

module interface_linear_solver

  use iso_fortran_env           , only : dp => REAL64
  use interface_algebraic_solver, only : algebraic_solver
  use interface_assembler       , only : assembler

  implicit none

  private
  public :: linear_solver
  public :: preconditioner

  !===================================================================!
  ! Linear solver datatype
  !===================================================================!

  type, abstract, extends(algebraic_solver) :: linear_solver

     real(dp) :: max_tol
     integer  :: max_it

   contains

     ! type bound procedures
     procedure(solve_interface), deferred :: solve

     ! Convergence monitor: the corrected-system residual r = (b+s) - A x of
     ! the system being solved, and a tabulated iteration history. Shared by
     ! every linear solver so they all report identically.
     procedure :: residual_norm
     procedure :: monitor_step
     procedure, private :: corrected_residual
     procedure, private :: monitor_header

  end type linear_solver

  !===================================================================!
  ! Abstract preconditioner: applies an approximate inverse z = M^-1 r.
  ! A linear solver (e.g. CG) calls apply once per iteration; concrete
  ! preconditioners (algebraic multigrid, jacobi, ...) extend this.
  !===================================================================!

  type, abstract :: preconditioner
   contains
     procedure(apply_interface), deferred :: apply
  end type preconditioner

  !===================================================================!
  ! Deferred interfaces
  !===================================================================!

  interface

     subroutine solve_interface(this, system, x, mode)
       import linear_solver
       import assembler
       import dp
       class(linear_solver)  , intent(in)            :: this
       class(assembler)      , intent(in)            :: system
       real(dp), allocatable , intent(out)           :: x(:)
       integer               , intent(in) , optional :: mode  ! FORWARD (default) / REVERSE
     end subroutine solve_interface

     ! z = M^-1 r  (the approximate-inverse action)
     subroutine apply_interface(this, r, z)
       import preconditioner
       import dp
       class(preconditioner), intent(in)  :: this
       real(dp)             , intent(in)  :: r(:)
       real(dp)             , intent(out) :: z(:)
     end subroutine apply_interface

  end interface

contains

  !===================================================================!
  ! Norm of the system's corrected residual at x - used to seed the drop
  ! with the initial-guess residual before a solver's outer loop starts.
  !===================================================================!

  pure function residual_norm(this, sys, x) result(rnorm)

    class(linear_solver), intent(in) :: this
    class(assembler)    , intent(in) :: sys
    real(dp)            , intent(in) :: x(:)

    real(dp)              :: rnorm
    real(dp), allocatable :: r(:)
    real(dp)              :: bnorm, snorm

    call this % corrected_residual(sys, x, r, bnorm, snorm)

    rnorm = norm2(r)

  end function residual_norm

  !===================================================================!
  ! Corrected residual  r = (b + s) - A x of the system sys, returning
  ! also the source and skew norms ||b|| and ||s|| for the diagnostics.
  !===================================================================!

  pure subroutine corrected_residual(this, sys, x, r, bnorm, snorm)

    class(linear_solver) , intent(in)  :: this
    class(assembler)     , intent(in)  :: sys
    real(dp)             , intent(in)  :: x(:)
    real(dp), allocatable, intent(out) :: r(:)
    real(dp)             , intent(out) :: bnorm
    real(dp)             , intent(out) :: snorm

    real(dp), allocatable :: b(:), s(:), ax(:)

    allocate(b, s, ax, mold = x)

    call sys % get_source(b)
    call sys % get_skew_source(s, x)
    call sys % get_jacobian_vector_product(ax, x)

    r     = (b + s) - ax
    bnorm = norm2(b)
    snorm = norm2(s)

  end subroutine corrected_residual

  !===================================================================!
  ! Column header for the convergence table
  !===================================================================!

  impure subroutine monitor_header(this)

    class(linear_solver), intent(in) :: this

    write(*,'(2x,a5,2x,a5,5(2x,a13))') &
         & "iter", "inner", "||r||", "||r||/||r0||", "||dx||", &
         & "||dx||/||x||", "||s||/||b||"

    write(*,'(2x,a5,2x,a5,5(2x,a13))') &
         & "----", "-----", "-------------", "-------------", &
         & "-------------", "-------------", "-------------"

  end subroutine monitor_header

  !===================================================================!
  ! One row of the convergence table: residual (absolute and as a drop
  ! from the initial residual rnorm0), solution update (absolute and
  ! relative), and the skew fraction. Prints the header on iteration 1.
  !===================================================================!

  impure subroutine monitor_step(this, sys, iter, inner, x, xold, rnorm0)

    class(linear_solver), intent(in) :: this
    class(assembler)    , intent(in) :: sys
    integer             , intent(in) :: iter
    integer             , intent(in) :: inner
    real(dp)            , intent(in) :: x(:)
    real(dp)            , intent(in) :: xold(:)
    real(dp)            , intent(in) :: rnorm0

    real(dp), allocatable :: r(:)
    real(dp)              :: r_abs, r_drop, du_abs, du_rel, skew
    real(dp)              :: bnorm, snorm, xnorm

    call this % corrected_residual(sys, x, r, bnorm, snorm)

    r_abs  = norm2(r)
    du_abs = norm2(x - xold)
    xnorm  = norm2(x)

    ! Residual drop from the initial-guess residual
    r_drop = r_abs
    if (rnorm0 .gt. epsilon(1.0_dp)) r_drop = r_abs/rnorm0

    ! Relative solution update
    du_rel = du_abs
    if (xnorm .gt. epsilon(1.0_dp)) du_rel = du_abs/xnorm

    ! Fraction of the rhs carried by the non-orthogonal skew correction
    skew = snorm
    if (bnorm .gt. epsilon(1.0_dp)) skew = snorm/bnorm

    if (iter .eq. 1) call this % monitor_header()

    write(*,'(2x,i5,2x,i5,5(2x,es13.5))') &
         & iter, inner, r_abs, r_drop, du_abs, du_rel, skew

  end subroutine monitor_step

end module interface_linear_solver
