!=====================================================================!
! Convergence monitor for the linear solvers. Prints a tabulated
! history of the outer (defect-correction) iterations - the usual cfd
! picture of residual and solution-update norms:
!
!   iter  inner    ||r||  ||r||/||r0||   ||dx||  ||dx||/||x||  ||s||/||b||
!   ----  -----  -------  ------------  -------  ------------  -----------
!      1     52  2.2E-12       1.2E-14  1.9E+02       1.1E+00      0.0E+00
!
! The residual is that of the corrected system  r = (b + s) - A x, where
! s is the non-orthogonal skew source. Two extra columns help diagnose:
!   ||r||/||r0|| : drop from the initial-guess residual (orders gained)
!   ||s||/||b||  : how much of the rhs is the skew correction - a direct
!                  fingerprint of mesh non-orthogonality (0 if orthogonal)
!
! Shared by cg / sor / gauss-seidel / gauss-jacobi so they all report
! identically. Seed the drop with rnorm0 = residual_norm(sys, x0).
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module module_solver_monitor

  use iso_fortran_env , only : dp => real64
  use class_assembler , only : assembler

  implicit none

  private
  public :: monitor_header, monitor_step, residual_norm

contains

  !===================================================================!
  ! Column header for the convergence table
  !===================================================================!

  subroutine monitor_header()

    write(*,'(2x,a5,2x,a5,5(2x,a13))') &
         & "iter", "inner", "||r||", "||r||/||r0||", "||dx||", &
         & "||dx||/||x||", "||s||/||b||"

    write(*,'(2x,a5,2x,a5,5(2x,a13))') &
         & "----", "-----", "-------------", "-------------", &
         & "-------------", "-------------", "-------------"

  end subroutine monitor_header

  !===================================================================!
  ! One row of the table: residual (absolute and as a drop from the
  ! initial residual rnorm0), solution update (absolute and relative),
  ! and the skew fraction. Prints the header on the first iteration.
  !===================================================================!

  subroutine monitor_step(sys, iter, inner, x, xold, rnorm0)

    class(assembler), intent(in) :: sys
    integer         , intent(in) :: iter
    integer         , intent(in) :: inner
    real(dp)        , intent(in) :: x(:)
    real(dp)        , intent(in) :: xold(:)
    real(dp)        , intent(in) :: rnorm0

    real(dp), allocatable :: r(:)
    real(dp)              :: r_abs, r_drop, du_abs, du_rel, skew
    real(dp)              :: bnorm, snorm, xnorm

    call corrected_residual(sys, x, r, bnorm, snorm)

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

    if (iter .eq. 1) call monitor_header()

    write(*,'(2x,i5,2x,i5,5(2x,es13.5))') &
         & iter, inner, r_abs, r_drop, du_abs, du_rel, skew

  end subroutine monitor_step

  !===================================================================!
  ! Norm of the corrected residual at x - used to seed the drop with
  ! the initial-guess residual before the outer loop starts.
  !===================================================================!

  function residual_norm(sys, x) result(rnorm)

    class(assembler), intent(in) :: sys
    real(dp)        , intent(in) :: x(:)

    real(dp)              :: rnorm
    real(dp), allocatable :: r(:)
    real(dp)              :: bnorm, snorm

    call corrected_residual(sys, x, r, bnorm, snorm)

    rnorm = norm2(r)

  end function residual_norm

  !===================================================================!
  ! Corrected residual  r = (b + s) - A x, returning also the source
  ! and skew norms ||b|| and ||s|| for the diagnostics (module-private)
  !===================================================================!

  subroutine corrected_residual(sys, x, r, bnorm, snorm)

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

end module module_solver_monitor
