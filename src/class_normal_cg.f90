!=====================================================================!
! Conjugate gradient on the normal equations - nonsymmetric solvers that
! reuse the CG recurrence through the two directions of the ONE product
! (system % get_jacobian_residual_product, FORWARD and REVERSE).
!
!   cgnr - CG on A^T A x = A^T b. minimizes ||b - A x||_2 (= LSQR).
!   cgne - CG on A A^T y = b, x = A^T y (Craig). minimizes ||x - x*||_2.
!
! Both are robust (the normal-equations operator is SPD for any nonsingular
! A, so they never break down) and need only the forward and transpose
! products. BUT they operate on the normal equations, so their convergence is
! governed by kappa(A)^2, not kappa(A): on advection-dominated (highly
! non-normal) operators that squaring makes them converge far slower than
! GMRES. Use them as a cheap, robust baseline for mildly nonsymmetric /
! diffusion-dominated problems; reach for class_gmres when advection bites.
! See test/krylov for the kappa^2 crossover.
!
! Each iteration costs one forward product (A p) and one transpose
! product (A^T r) - assembled entries are never touched (assembled only
! to build preconditioners, never to iterate). The convergence test is on
! the true residual ||b - A x||/||b||. cgnr/cgne are methods on the
! normal_cg solver: iterate dispatches to the kernel; the standalone
! tests wrap a raw operator in the csr_system fixture and call the
! method directly.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_normal_cg

  use iso_fortran_env        , only : dp => REAL64
  use interface_linear_solver, only : linear_solver
  use interface_assembler    , only : assembler
  use module_solve_mode      , only : REVERSE

  implicit none

  private
  public :: normal_cg, CGNR_METHOD, CGNE_METHOD

  ! normal-equations method selector for the linear_solver wrapper
  integer, parameter :: CGNR_METHOD = 1   ! CG on A^T A  (residual-minimizing)
  integer, parameter :: CGNE_METHOD = 2   ! CG on A A^T  (Craig, error-minimizing)

  !-------------------------------------------------------------------!
  ! linear_solver wrapper so the config-driven driver can pick a
  ! normal-equations CG ("cgnr"/"cgne") the way it picks "cg"/"gmres".
  ! Assembles the operator + rhs once and runs the kernel below. Robust
  ! on (mildly) nonsymmetric operators; converges on kappa(A)^2 - reach
  ! for gmres_solver when advection dominates.
  !-------------------------------------------------------------------!

  type, extends(linear_solver) :: normal_cg

     integer :: method = CGNR_METHOD

   contains

     ! the sweep consumed by the inherited outer iteration
     procedure :: iterate
     procedure :: cgnr
     procedure :: cgne

  end type normal_cg

  interface normal_cg
     module procedure construct
  end interface normal_cg

contains

  !===================================================================!
  ! Constructor for the normal_cg linear-solver wrapper
  !===================================================================!

  pure type(normal_cg) function construct(max_it, max_tol, method, &
       & print_level) result(this)

    integer        , intent(in)           :: max_it
    real(dp)       , intent(in)           :: max_tol
    integer        , intent(in), optional :: method
    integer        , intent(in), optional :: print_level

    this % max_it  = max_it
    this % max_tol = max_tol

    if (present(method))      this % method      = method
    if (present(print_level)) this % print_level = print_level

  end function construct

  !===================================================================!
  ! The sweep: assemble the operator and run CGNR or CGNE on the
  ! correction equation A dx = r from dx = 0 (both directions of the
  ! one product: the normal equations are the slope-zero condition of
  ! minimizing the residual's size).
  !===================================================================!

  impure subroutine iterate(this, system, r, dx, iter)

    class(normal_cg)     , intent(inout) :: this
    class(assembler)     , intent(in)    :: system
    real(dp)             , intent(in)    :: r(:)
    real(dp)             , intent(out)   :: dx(:)
    integer              , intent(out)   :: iter

    if (this % method .eq. CGNE_METHOD) then
       call this % cgne(system, r, dx, iter)
    else
       call this % cgnr(system, r, dx, iter)
    end if

  end subroutine iterate

  !===================================================================!
  ! CGNR: CG on A^T A x = A^T b (residual-minimizing / LSQR).
  !===================================================================!

  impure subroutine cgnr(this, system, b, x, iters)

    class(normal_cg), intent(in)  :: this
    class(assembler), intent(in)  :: system
    real(dp)        , intent(in)  :: b(:)
    real(dp)        , intent(out) :: x(:)
    integer         , intent(out) :: iters

    real(dp), allocatable :: r(:), z(:), p(:), w(:)
    real(dp)              :: alpha, beta, zz, zz_new, bnorm, tol
    integer               :: m, n, iter

    m = system % num_state_vars
    n = system % num_state_vars
    allocate(r(m), w(m), z(n), p(n))

    x = 0.0_dp
    r = b                              ! r = b - A x0,  x0 = 0
    call system % get_jacobian_residual_product(z, r, mode = REVERSE)   ! z = A^T r
    p     = z
    zz    = dot_product(z, z)
    bnorm = norm2(b)
    if (bnorm .eq. 0.0_dp) bnorm = 1.0_dp
    tol   = norm2(r)/bnorm

    iter = 0
    do while (tol .gt. this % max_tol .and. iter .lt. this % max_it)

       call system % get_jacobian_residual_product(w, p)                ! w = A p
       alpha = zz / dot_product(w, w)
       x = x + alpha*p
       r = r - alpha*w

       call system % get_jacobian_residual_product(z, r, mode = REVERSE)   ! z = A^T r
       zz_new = dot_product(z, z)
       beta   = zz_new / zz
       p  = z + beta*p
       zz = zz_new

       tol  = norm2(r)/bnorm
       iter = iter + 1
       if (this % print_level .gt. 1) write(*,'(2x,a,i6,a,es12.5)') "cgnr ", iter, "  rel res ", tol

    end do

    iters = iter
    if (this % print_level .gt. 0) write(*,'(1x,a,i0,a,es12.5)') &
         & "cgnr: ", iter, " iters, rel res ", tol

  end subroutine cgnr

  !===================================================================!
  ! CGNE (Craig): CG on A A^T y = b, x = A^T y (error-minimizing).
  !===================================================================!

  impure subroutine cgne(this, system, b, x, iters)

    class(normal_cg), intent(in)  :: this
    class(assembler), intent(in)  :: system
    real(dp)        , intent(in)  :: b(:)
    real(dp)        , intent(out) :: x(:)
    integer         , intent(out) :: iters

    real(dp), allocatable :: r(:), p(:), w(:), atr(:)
    real(dp)              :: alpha, beta, rr, rr_new, bnorm, tol
    integer               :: m, n, iter

    m = system % num_state_vars
    n = system % num_state_vars
    allocate(r(m), w(m), p(n), atr(n))

    x = 0.0_dp
    r = b                              ! r = b - A x0
    call system % get_jacobian_residual_product(p, r, mode = REVERSE)   ! p = A^T r
    rr    = dot_product(r, r)
    bnorm = norm2(b)
    if (bnorm .eq. 0.0_dp) bnorm = 1.0_dp
    tol   = sqrt(rr)/bnorm

    iter = 0
    do while (tol .gt. this % max_tol .and. iter .lt. this % max_it)

       call system % get_jacobian_residual_product(w, p)                ! w = A p
       alpha = rr / dot_product(p, p)
       x = x + alpha*p
       r = r - alpha*w
       rr_new = dot_product(r, r)
       beta   = rr_new / rr

       call system % get_jacobian_residual_product(atr, r, mode = REVERSE)   ! atr = A^T r
       p  = atr + beta*p
       rr = rr_new

       tol  = sqrt(rr)/bnorm
       iter = iter + 1
       if (this % print_level .gt. 1) write(*,'(2x,a,i6,a,es12.5)') "cgne ", iter, "  rel res ", tol

    end do

    iters = iter
    if (this % print_level .gt. 0) write(*,'(1x,a,i0,a,es12.5)') &
         & "cgne: ", iter, " iters, rel res ", tol

  end subroutine cgne

end module class_normal_cg
