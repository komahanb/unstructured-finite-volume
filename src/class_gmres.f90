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
! drives the stopping test directly, and class_algebraic_multigrid (or any
! preconditioner) plugs in through the optional member.
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
  use interface_linear_solver , only : linear_solver, preconditioner
  use interface_assembler     , only : assembler

  implicit none

  private
  public :: gmres_solver
  public :: gmres_last_iters          ! total inner iterations of the last solve

  ! Total inner iterations of the most recent solve. Written by the gmres
  ! kernel (which takes `this` as intent(in), so this cannot live on the
  ! object); read by tests comparing iteration counts across operators.
  integer :: gmres_last_iters = 0

  !-------------------------------------------------------------------!
  ! linear_solver wrapper around restarted GMRES so the config-driven
  ! driver can pick "gmres" the way it picks "cg". The gmres kernel below
  ! is a method that runs on any operator it is handed; solve assembles
  ! the csv via get_operator_csr / get_source first. Use it for
  ! nonsymmetric operators (advection) where CG does not apply.
  !-------------------------------------------------------------------!

  type, extends(linear_solver) :: gmres_solver

     integer :: restart = 200

   contains

     ! the sweep consumed by the inherited march (the kernel below runs
     ! on any operator it is handed)
     procedure :: iterate
     procedure :: gmres

  end type gmres_solver

  interface gmres_solver
     module procedure construct
  end interface gmres_solver

contains

  !===================================================================!
  ! Constructor for the gmres_solver linear-solver wrapper. FVAssembler is
  ! optional: omit it to drive the gmres kernel on a raw csr operator.
  !===================================================================!

  pure type(gmres_solver) function construct(max_it, max_tol, restart, &
       & print_level, precond) result(this)

    integer              , intent(in)           :: max_it
    real(dp)             , intent(in)           :: max_tol
    integer              , intent(in), optional :: restart
    integer              , intent(in), optional :: print_level
    class(preconditioner), intent(in), optional :: precond

    this % max_it  = max_it
    this % max_tol = max_tol

    this % res_file = 'gmres.res'

    if (present(restart))     this % restart     = restart
    if (present(print_level)) this % print_level = print_level
    ! the right preconditioner fills the post-operator slot (applied per
    ! arnoldi vector, back-mapped inside the kernel)
    if (present(precond))     allocate(this % post_conditioner, source = precond)

  end function construct

  !===================================================================!
  ! The sweep: assemble the operator and run restarted GMRES on the
  ! correction equation A dx = r from dx = 0.
  !===================================================================!

  impure subroutine iterate(this, system, r, dx, iter)

    class(gmres_solver)  , intent(in)  :: this
    class(assembler)     , intent(in)  :: system
    real(dp)             , intent(in)  :: r(:)
    real(dp)             , intent(out) :: dx(:)
    integer              , intent(out) :: iter

    type(csr_matrix) :: A

    call system % get_operator_csr(A)

    call this % gmres(A, r, dx)
    iter = gmres_last_iters

  end subroutine iterate

  !===================================================================!
  ! Solve A x = b with restarted, optionally right-preconditioned GMRES.
  ! Caps, Krylov dimension and verbosity are read from `this`; the
  ! preconditioner member (if allocated) supplies z = M^-1 r.
  !===================================================================!

  impure subroutine gmres(this, A, b, x)

    class(gmres_solver), intent(in)  :: this
    type(csr_matrix)   , intent(in)  :: A
    real(dp)           , intent(in)  :: b(:)
    real(dp)           , intent(out) :: x(:)

    real(dp), allocatable :: V(:,:), H(:,:), cs(:), sn(:), g(:), y(:)
    real(dp), allocatable :: w(:), z(:), r(:), u(:), du(:)
    real(dp)              :: bnorm, beta, tol, tmp
    integer               :: n, m, i, j, kk, total
    logical               :: converged

    n = A % ncols
    m = min(this % restart, this % max_it)

    allocate(V(n, m+1), H(m+1, m), cs(m), sn(m), g(m+1), y(m))
    allocate(w(n), z(n), r(n), u(n), du(n))

    x = 0.0_dp
    bnorm = norm2(b)
    if (bnorm .eq. 0.0_dp) bnorm = 1.0_dp

    total     = 0
    converged = .false.
    tol       = huge(1.0_dp)

    outer: do while (total .lt. this % max_it .and. .not. converged)

       ! restart residual r = b - A x
       call A % matvec(x, w)
       r    = b - w
       beta = norm2(r)
       tol  = beta/bnorm
       if (tol .le. this % max_tol) then
          converged = .true.
          exit outer
       end if

       V(:,1) = r/beta
       g      = 0.0_dp
       g(1)   = beta

       j = 0
       inner: do while (j .lt. m .and. total .lt. this % max_it)

          j     = j + 1
          total = total + 1

          ! w = A M^-1 v_j  (right preconditioning; w = A v_j if none)
          if (allocated(this % post_conditioner)) then
             call this % post_conditioner % apply(V(:,j), z)
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
          if (this % print_level .gt. 1) then
             write(*,'(2x,a,i6,a,es12.5)') "gmres ", total, "  rel res ", tol
          end if
          if (tol .le. this % max_tol) then
             converged = .true.
             exit inner
          end if

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
       if (allocated(this % post_conditioner)) then
          call this % post_conditioner % apply(u, du)
          x = x + du
       else
          x = x + u
       end if

    end do outer

    gmres_last_iters = total
    if (this % print_level .gt. 0) then
       write(*,'(1x,a,i0,a,i0,a,es12.5)') &
            & "gmres(", m, "): ", total, " iters, rel res ", tol
    end if

  end subroutine gmres

  !===================================================================!
  ! Givens rotation [c s; -s c] [a; b] = [r; 0]  (drotg-style, robust)
  !===================================================================!

  elemental subroutine givens(a, b, c, s)

    real(dp), intent(in)  :: a, b
    real(dp), intent(out) :: c, s

    real(dp) :: t, rr

    if (b .eq. 0.0_dp) then
       c = 1.0_dp
       s = 0.0_dp
    else if (abs(b) .gt. abs(a)) then
       t = a/b
       rr = sqrt(1.0_dp + t*t)
       s = 1.0_dp/rr
       c = s*t
    else
       t = b/a
       rr = sqrt(1.0_dp + t*t)
       c = 1.0_dp/rr
       s = c*t
    end if

  end subroutine givens

end module class_gmres
