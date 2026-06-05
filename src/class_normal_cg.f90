!=====================================================================!
! Conjugate gradient on the normal equations - nonsymmetric solvers that
! reuse the CG recurrence through a transpose matvec (csr_matrix has both
! matvec and matvec_transpose).
!
!   cgnr - CG on A^T A x = A^T b. minimizes ||b - A x||_2 (= LSQR).
!   cgne - CG on A A^T y = b, x = A^T y (Craig). minimizes ||x - x*||_2.
!
! Both are robust (the normal-equations operator is SPD for any nonsingular
! A, so they never break down) and need only A and A^T - which we already
! have. BUT they operate on the normal equations, so their convergence is
! governed by kappa(A)^2, not kappa(A): on advection-dominated (highly
! non-normal) operators that squaring makes them converge far slower than
! GMRES. Use them as a cheap, robust baseline for mildly nonsymmetric /
! diffusion-dominated problems; reach for class_gmres when advection bites.
! See test/krylov for the kappa^2 crossover.
!
! Each iteration costs one matvec (A p) and one transpose matvec (A^T r).
! The convergence test is on the true residual ||b - A x||/||b||.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_normal_cg

  use iso_fortran_env, only : dp => REAL64
  use class_csr      , only : csr_matrix

  implicit none

  private
  public :: cgnr, cgne
  public :: cgnr_last_iters, cgne_last_iters

  integer :: cgnr_last_iters = 0
  integer :: cgne_last_iters = 0

contains

  !===================================================================!
  ! CGNR: CG on A^T A x = A^T b (residual-minimizing / LSQR).
  !===================================================================!

  subroutine cgnr(A, b, x, max_it, max_tol, print_level)

    type(csr_matrix), intent(in)  :: A
    real(dp)        , intent(in)  :: b(:)
    real(dp)        , intent(out) :: x(:)
    integer         , intent(in)  :: max_it
    real(dp)        , intent(in)  :: max_tol
    integer         , intent(in)  :: print_level

    real(dp), allocatable :: r(:), z(:), p(:), w(:)
    real(dp) :: alpha, beta, zz, zz_new, bnorm, tol
    integer  :: m, n, iter

    m = A % nrows; n = A % ncols
    allocate(r(m), w(m), z(n), p(n))

    x = 0.0_dp
    r = b                              ! r = b - A x0,  x0 = 0
    call A % matvec_transpose(r, z)    ! z = A^T r
    p  = z
    zz = dot_product(z, z)
    bnorm = norm2(b); if (bnorm .eq. 0.0_dp) bnorm = 1.0_dp
    tol   = norm2(r)/bnorm

    iter = 0
    do while (tol .gt. max_tol .and. iter .lt. max_it)
       call A % matvec(p, w)           ! w = A p
       alpha = zz / dot_product(w, w)
       x = x + alpha*p
       r = r - alpha*w
       call A % matvec_transpose(r, z) ! z = A^T r
       zz_new = dot_product(z, z)
       beta   = zz_new / zz
       p  = z + beta*p
       zz = zz_new
       tol = norm2(r)/bnorm
       iter = iter + 1
       if (print_level .gt. 1) write(*,'(2x,a,i6,a,es12.5)') "cgnr ", iter, "  rel res ", tol
    end do

    cgnr_last_iters = iter
    if (print_level .gt. 0) write(*,'(1x,a,i0,a,es12.5)') &
         & "cgnr: ", iter, " iters, rel res ", tol

  end subroutine cgnr

  !===================================================================!
  ! CGNE (Craig): CG on A A^T y = b, x = A^T y (error-minimizing).
  !===================================================================!

  subroutine cgne(A, b, x, max_it, max_tol, print_level)

    type(csr_matrix), intent(in)  :: A
    real(dp)        , intent(in)  :: b(:)
    real(dp)        , intent(out) :: x(:)
    integer         , intent(in)  :: max_it
    real(dp)        , intent(in)  :: max_tol
    integer         , intent(in)  :: print_level

    real(dp), allocatable :: r(:), p(:), w(:), atr(:)
    real(dp) :: alpha, beta, rr, rr_new, bnorm, tol
    integer  :: m, n, iter

    m = A % nrows; n = A % ncols
    allocate(r(m), w(m), p(n), atr(n))

    x = 0.0_dp
    r = b                              ! r = b - A x0
    call A % matvec_transpose(r, p)    ! p = A^T r
    rr = dot_product(r, r)
    bnorm = norm2(b); if (bnorm .eq. 0.0_dp) bnorm = 1.0_dp
    tol   = sqrt(rr)/bnorm

    iter = 0
    do while (tol .gt. max_tol .and. iter .lt. max_it)
       call A % matvec(p, w)           ! w = A p
       alpha = rr / dot_product(p, p)
       x = x + alpha*p
       r = r - alpha*w
       rr_new = dot_product(r, r)
       beta   = rr_new / rr
       call A % matvec_transpose(r, atr)   ! atr = A^T r
       p  = atr + beta*p
       rr = rr_new
       tol = sqrt(rr)/bnorm
       iter = iter + 1
       if (print_level .gt. 1) write(*,'(2x,a,i6,a,es12.5)') "cgne ", iter, "  rel res ", tol
    end do

    cgne_last_iters = iter
    if (print_level .gt. 0) write(*,'(1x,a,i0,a,es12.5)') &
         & "cgne: ", iter, " iters, rel res ", tol

  end subroutine cgne

end module class_normal_cg
