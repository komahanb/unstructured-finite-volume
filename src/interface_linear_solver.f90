!=====================================================================!
! Interface for linear solvers, plus the preconditioner contract they
! consume. A linear solver solves A x = b for x; a preconditioner applies
! an approximate inverse z = M^-1 r once per iteration. Both abstract
! types live here so a solver and its preconditioner share one home.
!
! The default solve is the deferred-correction outer loop shared by
! every solver built on a correction sweep: start from x = b, run the
! sweep (`iterate`, the hook a correction solver overrides), refresh the
! non-orthogonal skew source at the new solution, and repeat until the
! solution stops moving. cg / gauss-seidel / jacobi / sor differ only in
! their sweep; solvers with a different loop structure (gmres, normal
! cg, distributed cg) override solve itself.
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
     integer  :: print_level = 0

     ! iteration-history file written when print_level == -1
     character(len=:), allocatable :: res_file

   contains

     ! Default solve = the deferred-correction outer loop around the
     ! correction sweep `iterate`. correction_solve is also bound by name
     ! so an override of solve can pre-process and then delegate to it.
     procedure :: solve => correction_solve
     procedure :: correction_solve
     procedure :: iterate

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
  ! Default solve: the deferred-correction outer loop. The initial guess
  ! is the source b; each pass corrects the right-hand side with the
  ! skew source at the current solution (zero on the first pass), runs
  ! the solver's sweep, and stops when the solution stops moving. On an
  ! orthogonal mesh the skew source is identically zero and the loop
  ! exits after one corrected pass.
  !===================================================================!

  impure subroutine correction_solve(this, system, x, mode)

    class(linear_solver)  , intent(in)           :: this
    class(assembler)      , intent(in)           :: system
    real(dp), allocatable , intent(out)          :: x(:)
    integer               , intent(in), optional :: mode  ! FORWARD (default) / REVERSE

    ! Locals
    real(dp), allocatable :: xold(:), ss(:)
    real(dp) :: update_norm, rnorm0
    integer  :: iter, num_inner_iters

    ! Initial guess vector for the subspace is "b"
    allocate(x(system % num_state_vars))
    call system % get_source(x)
    if (norm2(x) .lt. epsilon(1.0_dp)) then
       write(*,*) 'warning: zero right-hand side; solution is x = 0'
    end if

    allocate(xold, ss, mold = x)

    if (this % print_level .eq. -1) then
       open(13, file = history_file(this), action = 'write')
       write(13,*) "iteration ", " residual"
    end if

    rnorm0 = 0.0_dp
    if (this % print_level .gt. 0) rnorm0 = this % residual_norm(system, x)

    update_norm = huge(1.0d0)
    iter = 1
    outer_iterations: do while (update_norm .gt. this % max_tol .and. iter .le. this % max_it)

       xold = x

       ! Inner iterations with the solver's correction sweep
       if (iter .eq. 1) then
          ss = 0.0d0
       else
          call system % get_skew_source(ss, x)
       end if
       call this % iterate(system, x, ss, num_inner_iters)

       update_norm = norm2(x - xold)
       if (this % print_level .eq. -1) then
          write(13, *) iter, update_norm
       end if
       if (this % print_level .gt. 0) then
          call this % monitor_step(system, iter, num_inner_iters, x, xold, rnorm0)
       end if
       iter = iter + 1

    end do outer_iterations

    if (this % print_level .eq. -1) then
       close(13)
    end if

    ! Safety net: report if the outer (deferred-correction) loop was
    ! capped at max_it without meeting the tolerance.
    if (update_norm .gt. this % max_tol) then
       write(*,*) "warning: outer correction loop hit max_it without convergence; update_norm =", update_norm
    end if

    deallocate(xold)
    deallocate(ss)

  end subroutine correction_solve

  !===================================================================!
  ! Correction sweep consumed by the default solve: drive the orthogonal
  ! part of the system at the given skew-corrected right-hand side.
  ! Correction solvers (cg, gauss-seidel, jacobi, sor) override this;
  ! solvers that override solve itself never reach it.
  !===================================================================!

  impure subroutine iterate(this, system, x, ss, iter)

    class(linear_solver), intent(in)    :: this
    class(assembler)    , intent(in)    :: system
    real(dp)            , intent(inout) :: x(:)
    real(dp)            , intent(in)    :: ss(:)
    integer             , intent(out)   :: iter

    iter = 0
    error stop 'linear_solver: no correction sweep (iterate) for this solver'

  end subroutine iterate

  !===================================================================!
  ! Name of the iteration-history file (print_level == -1)
  !===================================================================!

  pure function history_file(this) result(fname)

    class(linear_solver), intent(in) :: this
    character(len=:), allocatable    :: fname

    if (allocated(this % res_file)) then
       fname = this % res_file
    else
       fname = 'solver.res'
    end if

  end function history_file

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
