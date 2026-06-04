!=====================================================================!
! Coupled-diffusion adjoint gradient (the problem of
! src/coupled_diffusion_adjoint_fv_implementation_spec.md, implemented
! directly - not the spec's 5-layer framework).
!
! Two coupled fields q = [u, v] on a 1D box, steady state:
!
!     R_i = - sum_f A_f F_f  -  V_i S(q_i, p) = 0,
!     F_f = D (q_j - q_i)/d_ij          (centred diffusive flux)
!     S   = C q,  C = [[-alpha, alpha], [beta, -beta]]
!
! design p = [D, alpha, beta]; dirichlet on the two ends, no-flux else
! (1D -> only the two ends are boundaries). Tracking objective
!
!     J = integral_Omega [ 1/2 (q-q_d)^T W (q-q_d)
!                          + 1/2 sum_k gamma_k (p_k - p0_k)^2 ] dOmega.
!
! Adjoint:  R_q^T lambda = -J_q ,  dJ/dp = J_p + R_p^T lambda,
! verified against central finite differences (re-solving the primal at
! each perturbed design). Local physics/function partials are also
! checked against the spec's acceptance values (sections 23.1-23.8).
!
! A nonzero exit (error stop) means a check failed.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program coupled_diffusion_adjoint

  use iso_fortran_env, only : dp => real64

  implicit none

  external :: dgesv

  ! ---- discretisation / problem sizes -----------------------------!
  integer , parameter :: nc = 50          ! cells
  integer , parameter :: ns = 2           ! state components (u, v)
  integer , parameter :: npar = 3         ! design parameters (D, alpha, beta)
  integer , parameter :: ndof = nc*ns

  real(dp), parameter :: Lx = 1.0_dp
  real(dp), parameter :: dx   = Lx/real(nc, dp)
  real(dp), parameter :: vol  = dx          ! 1D cell volume (unit cross-section)
  real(dp), parameter :: area = 1.0_dp      ! face area
  real(dp), parameter :: dint = dx          ! interior cell-centre distance
  real(dp), parameter :: dbnd = 0.5_dp*dx   ! cell-centre to boundary distance

  ! ---- problem data (spec section 22) -----------------------------!
  real(dp), parameter :: qL(ns) = [1.0_dp, 0.0_dp]   ! left dirichlet
  real(dp), parameter :: qR(ns) = [0.0_dp, 1.0_dp]   ! right dirichlet
  real(dp), parameter :: qd(ns) = [0.5_dp, 0.5_dp]   ! tracking target
  real(dp), parameter :: Wm(ns) = [1.0_dp, 1.0_dp]   ! diag state weights (W = I)
  real(dp), parameter :: p0(npar)  = [0.05_dp, 1.0_dp, 0.5_dp]   ! design reference
  real(dp), parameter :: gam(npar) = [1.0e-4_dp, 1.0e-4_dp, 1.0e-4_dp]

  ! ---- working storage --------------------------------------------!
  real(dp) :: p(npar)
  real(dp) :: q(ndof), Rq(ndof,ndof), Rp(ndof,npar)
  real(dp) :: Jval, Jq(ndof), Jp(npar), lambda(ndof)
  real(dp) :: g(npar), gfd(npar)
  integer  :: nfail

  nfail = 0

  call local_law_checks(nfail)

  ! design at which the gradient is evaluated
  p = [0.05_dp, 1.0_dp, 0.5_dp]

  ! primal solve, functional + all partials at the converged state
  call primal_solve(p, q)
  call assemble_state_jacobian(q, p, Rq)
  call assemble_design_partials(q, p, Rp)
  call eval_functional(q, p, Jval, Jq, Jp)

  ! adjoint:  R_q^T lambda = -J_q
  call solve_transpose(Rq, -Jq, lambda)

  ! total gradient  dJ/dp = J_p + R_p^T lambda
  block
    integer :: k
    do k = 1, npar
       g(k) = Jp(k) + dot_product(Rp(:,k), lambda)
    end do
  end block

  ! finite-difference reference (re-solve the primal at each perturbation)
  call fd_gradient(p, gfd)

  ! ---- report ------------------------------------------------------!
  write(*,'(a)')        " coupled-diffusion adjoint gradient dJ/dp,  p = [D, alpha, beta]"
  write(*,'(a,es14.6)') "   J            = ", Jval
  write(*,'(a,3es14.6)')"   adjoint  g   = ", g
  write(*,'(a,3es14.6)')"   finite-diff  = ", gfd
  block
    integer  :: k
    real(dp) :: relerr
    do k = 1, npar
       relerr = abs(g(k) - gfd(k))/max(abs(gfd(k)), 1.0e-30_dp)
       write(*,'(a,i0,a,es12.4)') "   rel err [", k, "] = ", relerr
       if (relerr .gt. 1.0e-5_dp) nfail = nfail + 1
    end do
  end block

  write(*,*) "============================================="
  if (nfail .eq. 0) then
     write(*,*) "coupled-diffusion adjoint gradient matches finite difference"
  else
     write(*,*) nfail, " coupled-diffusion checks FAILED"
     error stop
  end if

contains

  !===================================================================!
  ! Flat dof index for component c (1..ns) of cell i (1..nc)
  !===================================================================!
  pure integer function idx(i, c)
    integer, intent(in) :: i, c
    idx = (i-1)*ns + c
  end function idx

  !===================================================================!
  ! Local physics laws (spec section 15.1) - pointwise, no mesh
  !===================================================================!

  ! source  S = C q :  S_u = alpha(v-u),  S_v = beta(u-v)
  pure function source(qi, pp) result(S)
    real(dp), intent(in) :: qi(ns), pp(npar)
    real(dp)             :: S(ns)
    S(1) = pp(2)*(qi(2) - qi(1))
    S(2) = pp(3)*(qi(1) - qi(2))
  end function source

  ! diffusive normal flux  F = D (qr - ql)/d  (area applied by caller)
  pure function flux(ql, qr, dd, pp) result(F)
    real(dp), intent(in) :: ql(ns), qr(ns), dd, pp(npar)
    real(dp)             :: F(ns)
    F = pp(1)*(qr - ql)/dd
  end function flux

  ! functional state partial  phi_q = W (q - q_d)
  pure function dphi_dq(qi) result(dq)
    real(dp), intent(in) :: qi(ns)
    real(dp)             :: dq(ns)
    dq = Wm*(qi - qd)
  end function dphi_dq

  ! functional design partial  phi_p = gamma_k (p_k - p0_k)
  pure function dphi_dp(pp) result(dpd)
    real(dp), intent(in) :: pp(npar)
    real(dp)             :: dpd(npar)
    dpd = gam*(pp - p0)
  end function dphi_dp

  !===================================================================!
  ! Steady residual  R(q, p)  (spec sign convention, section 24)
  !===================================================================!
  subroutine assemble_residual(qin, pp, R)
    real(dp), intent(in)  :: qin(ndof), pp(npar)
    real(dp), intent(out) :: R(ndof)
    real(dp) :: qi(ns), qj(ns), F(ns), S(ns)
    integer  :: i, c

    R = 0.0_dp

    ! cell source:  R_i -= V S(q_i)
    do i = 1, nc
       qi = qin(idx(i,1):idx(i,ns))
       S  = source(qi, pp)
       do c = 1, ns
          R(idx(i,c)) = R(idx(i,c)) - vol*S(c)
       end do
    end do

    ! interior faces (i, i+1):  R_i -= A F,  R_{i+1} += A F
    do i = 1, nc-1
       qi = qin(idx(i  ,1):idx(i  ,ns))
       qj = qin(idx(i+1,1):idx(i+1,ns))
       F  = flux(qi, qj, dint, pp)
       do c = 1, ns
          R(idx(i  ,c)) = R(idx(i  ,c)) - area*F(c)
          R(idx(i+1,c)) = R(idx(i+1,c)) + area*F(c)
       end do
    end do

    ! dirichlet ends:  R_i -= A F,  ghost = qL / qR, distance dbnd
    block
      real(dp) :: Fb(ns)
      Fb = flux(qin(idx(1,1):idx(1,ns)), qL, dbnd, pp)
      do c = 1, ns
         R(idx(1,c)) = R(idx(1,c)) - area*Fb(c)
      end do
      Fb = flux(qin(idx(nc,1):idx(nc,ns)), qR, dbnd, pp)
      do c = 1, ns
         R(idx(nc,c)) = R(idx(nc,c)) - area*Fb(c)
      end do
    end block
    ! (the other 1D "faces" are no-flux: zero contribution)

  end subroutine assemble_residual

  !===================================================================!
  ! State Jacobian  R_q = dR/dq  (dense)
  !===================================================================!
  subroutine assemble_state_jacobian(qin, pp, A)
    real(dp), intent(in)  :: qin(ndof), pp(npar)
    real(dp), intent(out) :: A(ndof,ndof)
    real(dp) :: kf, kb
    integer  :: i, c

    A  = 0.0_dp
    kf = pp(1)/dint          ! D/d interior
    kb = pp(1)/dbnd          ! D/d boundary

    ! source:  R_{q,ii} -= V C
    do i = 1, nc
       A(idx(i,1),idx(i,1)) = A(idx(i,1),idx(i,1)) - vol*(-pp(2))
       A(idx(i,1),idx(i,2)) = A(idx(i,1),idx(i,2)) - vol*( pp(2))
       A(idx(i,2),idx(i,1)) = A(idx(i,2),idx(i,1)) - vol*( pp(3))
       A(idx(i,2),idx(i,2)) = A(idx(i,2),idx(i,2)) - vol*(-pp(3))
    end do

    ! interior faces: F_qi = -D/d I, F_qj = +D/d I; assembler applies signs
    do i = 1, nc-1
       do c = 1, ns
          A(idx(i  ,c),idx(i  ,c)) = A(idx(i  ,c),idx(i  ,c)) + area*kf
          A(idx(i  ,c),idx(i+1,c)) = A(idx(i  ,c),idx(i+1,c)) - area*kf
          A(idx(i+1,c),idx(i  ,c)) = A(idx(i+1,c),idx(i  ,c)) - area*kf
          A(idx(i+1,c),idx(i+1,c)) = A(idx(i+1,c),idx(i+1,c)) + area*kf
       end do
    end do

    ! dirichlet ends: only the owner diagonal (ghost value is fixed)
    do c = 1, ns
       A(idx(1 ,c),idx(1 ,c)) = A(idx(1 ,c),idx(1 ,c)) + area*kb
       A(idx(nc,c),idx(nc,c)) = A(idx(nc,c),idx(nc,c)) + area*kb
    end do

  end subroutine assemble_state_jacobian

  !===================================================================!
  ! Residual design partials  R_p = dR/dp,  columns [D, alpha, beta]
  !===================================================================!
  subroutine assemble_design_partials(qin, pp, Rpout)
    real(dp), intent(in)  :: qin(ndof), pp(npar)
    real(dp), intent(out) :: Rpout(ndof,npar)
    real(dp) :: qi(ns), qj(ns), FpD(ns)
    integer  :: i, c

    Rpout = 0.0_dp

    ! --- D column (face flux):  F_D = (qr - ql)/d ; R_i -= A F_D, R_j += ---
    do i = 1, nc-1
       qi  = qin(idx(i  ,1):idx(i  ,ns))
       qj  = qin(idx(i+1,1):idx(i+1,ns))
       FpD = (qj - qi)/dint
       do c = 1, ns
          Rpout(idx(i  ,c),1) = Rpout(idx(i  ,c),1) - area*FpD(c)
          Rpout(idx(i+1,c),1) = Rpout(idx(i+1,c),1) + area*FpD(c)
       end do
    end do
    qi  = qin(idx(1,1):idx(1,ns)); FpD = (qL - qi)/dbnd
    do c = 1, ns
       Rpout(idx(1,c),1) = Rpout(idx(1,c),1) - area*FpD(c)
    end do
    qi  = qin(idx(nc,1):idx(nc,ns)); FpD = (qR - qi)/dbnd
    do c = 1, ns
       Rpout(idx(nc,c),1) = Rpout(idx(nc,c),1) - area*FpD(c)
    end do

    ! --- alpha, beta columns (source):  R_p -= V S_p ----------------!
    ! S_alpha = [v-u, 0],  S_beta = [0, u-v]
    do i = 1, nc
       qi = qin(idx(i,1):idx(i,ns))
       Rpout(idx(i,1),2) = Rpout(idx(i,1),2) - vol*(qi(2) - qi(1))
       Rpout(idx(i,2),3) = Rpout(idx(i,2),3) - vol*(qi(1) - qi(2))
    end do

  end subroutine assemble_design_partials

  !===================================================================!
  ! Functional value + partials  J, J_q, J_p
  !===================================================================!
  subroutine eval_functional(qin, pp, J, Jqout, Jpout)
    real(dp), intent(in)  :: qin(ndof), pp(npar)
    real(dp), intent(out) :: J, Jqout(ndof), Jpout(npar)
    real(dp) :: qi(ns), dq(ns), phi_design
    integer  :: i, c

    J      = 0.0_dp
    Jqout  = 0.0_dp
    Jpout  = 0.0_dp

    phi_design = 0.5_dp*sum(gam*(pp - p0)**2)

    do i = 1, nc
       qi = qin(idx(i,1):idx(i,ns))
       dq = qi - qd
       ! J += V (1/2 (q-qd)^T W (q-qd) + phi_design)
       J = J + vol*(0.5_dp*sum(Wm*dq**2) + phi_design)
       ! J_q += V phi_q
       do c = 1, ns
          Jqout(idx(i,c)) = Jqout(idx(i,c)) + vol*Wm(c)*dq(c)
       end do
       ! J_p += V phi_p
       Jpout = Jpout + vol*dphi_dp(pp)
    end do

  end subroutine eval_functional

  !===================================================================!
  ! Primal Newton solve:  R(q, p) = 0  (linear here -> one step)
  !===================================================================!
  subroutine primal_solve(pp, qout)
    real(dp), intent(in)  :: pp(npar)
    real(dp), intent(out) :: qout(ndof)
    real(dp) :: R(ndof), A(ndof,ndof), dq(ndof)
    integer  :: ipiv(ndof), info, it

    qout = 0.0_dp
    do it = 1, 20
       call assemble_residual(qout, pp, R)
       if (sqrt(sum(R**2)) .lt. 1.0e-13_dp) exit
       call assemble_state_jacobian(qout, pp, A)
       dq = -R
       call dgesv(ndof, 1, A, ndof, ipiv, dq, ndof, info)
       if (info .ne. 0) error stop "primal dgesv failed"
       qout = qout + dq
    end do

  end subroutine primal_solve

  !===================================================================!
  ! Solve the transpose system  A^T x = b
  !===================================================================!
  subroutine solve_transpose(A, b, x)
    real(dp), intent(in)  :: A(ndof,ndof), b(ndof)
    real(dp), intent(out) :: x(ndof)
    real(dp) :: At(ndof,ndof)
    integer  :: ipiv(ndof), info

    At = transpose(A)
    x  = b
    call dgesv(ndof, 1, At, ndof, ipiv, x, ndof, info)
    if (info .ne. 0) error stop "adjoint dgesv failed"

  end subroutine solve_transpose

  !===================================================================!
  ! Central finite-difference total gradient (re-solve the primal)
  !===================================================================!
  subroutine fd_gradient(pp, gout)
    real(dp), intent(in)  :: pp(npar)
    real(dp), intent(out) :: gout(npar)
    real(dp) :: pm(npar), qm(ndof), jp, jm, dummy_jq(ndof), dummy_jp(npar), eps
    integer  :: k

    do k = 1, npar
       eps = 1.0e-6_dp*max(1.0_dp, abs(pp(k)))

       pm = pp; pm(k) = pp(k) + eps
       call primal_solve(pm, qm)
       call eval_functional(qm, pm, jp, dummy_jq, dummy_jp)

       pm = pp; pm(k) = pp(k) - eps
       call primal_solve(pm, qm)
       call eval_functional(qm, pm, jm, dummy_jq, dummy_jp)

       gout(k) = (jp - jm)/(2.0_dp*eps)
    end do

  end subroutine fd_gradient

  !===================================================================!
  ! Local-law unit checks against the spec's acceptance values
  ! (sections 23.1, 23.4, 23.7, 23.8)
  !===================================================================!
  subroutine local_law_checks(nf)
    integer, intent(inout) :: nf
    real(dp) :: S(ns), F(ns), dq(ns), dpd(npar)

    ! 23.1 source:  q=[2,5], p=[0.1,3,4] -> S=[9,-12]
    S = source([2.0_dp, 5.0_dp], [0.1_dp, 3.0_dp, 4.0_dp])
    call check_close("source", S, [9.0_dp, -12.0_dp], nf)

    ! 23.4 flux:  qL=[1,2], qR=[3,5], D=0.1, d=0.5 -> [0.4,0.6]
    F = flux([1.0_dp,2.0_dp], [3.0_dp,5.0_dp], 0.5_dp, [0.1_dp,3.0_dp,4.0_dp])
    call check_close("flux", F, [0.4_dp, 0.6_dp], nf)

    ! 23.7 functional state partial: phi_q = W(q-qd); q=[2,5], qd=[.5,.5]
    dq = dphi_dq([2.0_dp, 5.0_dp])
    call check_close("phi_q", dq, [1.5_dp, 4.5_dp], nf)

    ! 23.8 functional design partial: phi_p = gamma (p-p0), offset [.01,.1,-.1]
    dpd = dphi_dp(p0 + [0.01_dp, 0.1_dp, -0.1_dp])
    call check_close("phi_p", dpd, gam*[0.01_dp, 0.1_dp, -0.1_dp], nf)

  end subroutine local_law_checks

  !===================================================================!
  ! Assert two vectors agree (local-law unit checks)
  !===================================================================!
  subroutine check_close(name, got, want, nf)
    character(*), intent(in)    :: name
    real(dp)    , intent(in)    :: got(:), want(:)
    integer     , intent(inout) :: nf
    if (maxval(abs(got - want)) .gt. 1.0e-12_dp) then
       write(*,'(a,a)') " local check FAILED: ", name
       nf = nf + 1
    end if
  end subroutine check_close

end program coupled_diffusion_adjoint
