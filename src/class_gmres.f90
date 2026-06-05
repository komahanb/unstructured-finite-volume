!=====================================================================!
! Restarted GMRES(m) - the general-purpose nonsymmetric Krylov solver.
!
! Unlike the normal-equations solvers (class_normal_cg) GMRES works on A
! directly, so its convergence is governed by kappa(A) rather than
! kappa(A)^2 - decisive on advection-dominated (non-normal) operators. It
! is also TRANSPOSE-FREE (only A*v), and preconditions naturally. The cost
! is memory/orthogonalization that grow with the Krylov dimension, hence
! the restart parameter m: build an m-dimensional Krylov space, form the
! minimum-residual solution there, restart from it until converged.
!
! Right preconditioning (solve A M^-1 u = b, x = M^-1 u): keeps the Arnoldi
! residual equal to the TRUE residual ||b - A x||, so the Givens estimate
! drives the stopping test directly, and class_amg (or any preconditioner)
! plugs in through the optional argument.
!
! Method: modified Gram-Schmidt Arnoldi -> upper Hessenberg H, reduced to
! upper triangular by Givens rotations as columns arrive (so the residual
! norm |g(j+1)| is available every step), then back-substitution for y and
! x <- x + M^-1 (V y).
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_gmres

  use iso_fortran_env         , only : dp => REAL64
  use class_csr               , only : csr_matrix
  use interface_preconditioner, only : preconditioner

  implicit none

  private
  public :: gmres
  public :: gmres_last_iters          ! total inner iterations of the last solve

  integer :: gmres_last_iters = 0

contains

  !===================================================================!
  ! Solve A x = b with restarted, optionally right-preconditioned GMRES.
  !   max_it  - cap on TOTAL inner iterations (across restarts)
  !   restart - Krylov dimension m before restarting
  !===================================================================!

  subroutine gmres(A, b, x, max_it, restart, max_tol, print_level, precond)

    type(csr_matrix)     , intent(in)           :: A
    real(dp)             , intent(in)           :: b(:)
    real(dp)             , intent(out)          :: x(:)
    integer              , intent(in)           :: max_it
    integer              , intent(in)           :: restart
    real(dp)             , intent(in)           :: max_tol
    integer              , intent(in)           :: print_level
    class(preconditioner), intent(in), optional :: precond

    real(dp), allocatable :: V(:,:), H(:,:), cs(:), sn(:), g(:), y(:)
    real(dp), allocatable :: w(:), z(:), r(:), u(:), du(:)
    real(dp) :: bnorm, beta, tol, tmp
    integer  :: n, m, i, j, kk, total
    logical  :: converged

    n = A % ncols
    m = min(restart, max_it)
    allocate(V(n, m+1), H(m+1, m), cs(m), sn(m), g(m+1), y(m))
    allocate(w(n), z(n), r(n), u(n), du(n))

    x = 0.0_dp
    bnorm = norm2(b); if (bnorm .eq. 0.0_dp) bnorm = 1.0_dp
    total     = 0
    converged = .false.
    tol       = huge(1.0_dp)

    outer: do while (total .lt. max_it .and. .not. converged)

       ! restart residual r = b - A x
       call A % matvec(x, w)
       r    = b - w
       beta = norm2(r)
       tol  = beta/bnorm
       if (tol .le. max_tol) then; converged = .true.; exit outer; end if

       V(:,1) = r/beta
       g      = 0.0_dp
       g(1)   = beta

       j = 0
       inner: do while (j .lt. m .and. total .lt. max_it)
          j     = j + 1
          total = total + 1

          ! w = A M^-1 v_j  (right preconditioning; w = A v_j if none)
          if (present(precond)) then
             call precond % apply(V(:,j), z)
             call A % matvec(z, w)
          else
             call A % matvec(V(:,j), w)
          end if

          ! modified Gram-Schmidt against the existing basis
          do i = 1, j
             H(i,j) = dot_product(w, V(:,i))
             w = w - H(i,j)*V(:,i)
          end do
          H(j+1,j) = norm2(w)
          if (H(j+1,j) .gt. 0.0_dp) V(:,j+1) = w / H(j+1,j)

          ! apply the previous Givens rotations to the new column
          do i = 1, j-1
             tmp      =  cs(i)*H(i,j) + sn(i)*H(i+1,j)
             H(i+1,j) = -sn(i)*H(i,j) + cs(i)*H(i+1,j)
             H(i,j)   = tmp
          end do

          ! new rotation that zeros H(j+1,j), then apply it to H and g
          call givens(H(j,j), H(j+1,j), cs(j), sn(j))
          H(j,j)   = cs(j)*H(j,j) + sn(j)*H(j+1,j)
          H(j+1,j) = 0.0_dp
          tmp    =  cs(j)*g(j)
          g(j+1) = -sn(j)*g(j)
          g(j)   = tmp

          tol = abs(g(j+1))/bnorm
          if (print_level .gt. 1) write(*,'(2x,a,i6,a,es12.5)') "gmres ", total, "  rel res ", tol
          if (tol .le. max_tol) then; converged = .true.; exit inner; end if
       end do inner

       ! back-substitution: H(1:j,1:j) y = g(1:j)
       do i = j, 1, -1
          y(i) = g(i)
          do kk = i+1, j
             y(i) = y(i) - H(i,kk)*y(kk)
          end do
          y(i) = y(i) / H(i,i)
       end do

       ! correction u = V(:,1:j) y, then x <- x + M^-1 u
       u = 0.0_dp
       do i = 1, j
          u = u + y(i)*V(:,i)
       end do
       if (present(precond)) then
          call precond % apply(u, du)
          x = x + du
       else
          x = x + u
       end if

    end do outer

    gmres_last_iters = total
    if (print_level .gt. 0) write(*,'(1x,a,i0,a,i0,a,es12.5)') &
         & "gmres(", m, "): ", total, " iters, rel res ", tol

  end subroutine gmres

  !===================================================================!
  ! Givens rotation [c s; -s c] [a; b] = [r; 0]  (drotg-style, robust)
  !===================================================================!

  subroutine givens(a, b, c, s)
    real(dp), intent(in)  :: a, b
    real(dp), intent(out) :: c, s
    real(dp) :: t, rr
    if (b .eq. 0.0_dp) then
       c = 1.0_dp; s = 0.0_dp
    else if (abs(b) .gt. abs(a)) then
       t = a/b; rr = sqrt(1.0_dp + t*t); s = 1.0_dp/rr; c = s*t
    else
       t = b/a; rr = sqrt(1.0_dp + t*t); c = 1.0_dp/rr; s = c*t
    end if
  end subroutine givens

end module class_gmres
