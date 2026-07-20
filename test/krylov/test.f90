!=====================================================================!
! Nonsymmetric Krylov solvers: CGNR / CGNE (normal equations) vs GMRES.
!
! The test operator is a 2D advection-diffusion 5-point stencil on an NxN
! interior grid (Dirichlet, so it is nonsingular): symmetric Laplacian
! (diag 4, off-diag -1) plus a SKEW advection term of magnitude gamma
! (east/north get -1+gamma, west/south -1-gamma). gamma = 0 is the
! symmetric Laplacian; growing gamma makes it increasingly nonsymmetric and
! non-normal - a tunable stand-in for the Peclet number.
!
! Checks:
!   1. correctness  - cgnr, cgne and gmres all recover the known solution
!                     (b is built as A x_exact) to the solver tolerance.
!   2. kappa^2      - sweep gamma and tabulate CGNR vs GMRES iterations.
!                     CGNR works on A^T A (condition number SQUARED), so it
!                     needs far more iterations than GMRES (which works on A)
!                     - and the gap is exactly why GMRES is needed. assert
!                     gmres_iters < cgnr_iters across the sweep.
!   3. restart      - GMRES(m) with a small restart still converges.
!   4. precond      - AMG (built on the symmetric part) as a right
!                     preconditioner cuts GMRES iterations sharply.
!
! A nonzero exit (error stop) means a check failed.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_krylov

  use iso_fortran_env , only : dp => real64
  use class_csr       , only : csr_matrix
  use class_normal_cg , only : normal_cg, CGNR_METHOD, CGNE_METHOD
  use class_gmres     , only : gmres_solver
  use class_csr_system, only : csr_system
  use class_algebraic_multigrid, only : algebraic_multigrid

  implicit none

  integer, parameter :: N = 24                  ! NxN interior grid
  real(dp), parameter :: TOL = 1.0e-8_dp
  integer :: nfail

  nfail = 0

  call check_correctness(nfail)
  call check_peclet_sweep(nfail)
  call check_restart(nfail)
  call check_amg_precond(nfail)

  write(*,'(a)') " ============================================="
  if (nfail .eq. 0) then
     write(*,'(a)') "  all krylov checks passed"
  else
     write(*,'(1x,i0,a)') nfail, "  krylov checks FAILED"
     error stop
  end if

contains

  !===================================================================!
  ! 2d advection-diffusion 5-point operator, skewness gamma, as CSR
  !===================================================================!
  function advdiff_csr(gamma) result(A)
    real(dp), intent(in) :: gamma
    type(csr_matrix) :: A
    integer, allocatable :: row_ptr(:), col_idx(:)
    real(dp), allocatable :: vals(:)
    integer :: nn, i, j, k, pos, cnt

    nn = N*N
    allocate(row_ptr(nn+1)); row_ptr(1) = 1
    do j = 1, N
       do i = 1, N
          k = (j-1)*N + i
          cnt = 1
          if (j .gt. 1) cnt = cnt + 1
          if (i .gt. 1) cnt = cnt + 1
          if (i .lt. N) cnt = cnt + 1
          if (j .lt. N) cnt = cnt + 1
          row_ptr(k+1) = row_ptr(k) + cnt
       end do
    end do
    allocate(col_idx(row_ptr(nn+1)-1), vals(row_ptr(nn+1)-1))

    pos = 1
    do j = 1, N
       do i = 1, N
          k = (j-1)*N + i
          if (j .gt. 1) then; col_idx(pos) = k-N; vals(pos) = -1.0_dp - gamma; pos = pos+1; end if ! south
          if (i .gt. 1) then; col_idx(pos) = k-1; vals(pos) = -1.0_dp - gamma; pos = pos+1; end if ! west
          col_idx(pos) = k; vals(pos) = 4.0_dp; pos = pos+1                                         ! diag
          if (i .lt. N) then; col_idx(pos) = k+1; vals(pos) = -1.0_dp + gamma; pos = pos+1; end if ! east
          if (j .lt. N) then; col_idx(pos) = k+N; vals(pos) = -1.0_dp + gamma; pos = pos+1; end if ! north
       end do
    end do

    A = csr_matrix(nn, nn, row_ptr, col_idx, vals)
  end function advdiff_csr

  !===================================================================!
  ! known solution and matching rhs b = A x_exact
  !===================================================================!
  subroutine make_problem(A, xex, b)
    type(csr_matrix), intent(in)  :: A
    real(dp), allocatable, intent(out) :: xex(:), b(:)
    integer :: k
    allocate(xex(A % ncols), b(A % num_vertices))
    do k = 1, A % ncols
       xex(k) = sin(0.1_dp*real(k,dp)) + 0.5_dp
    end do
    call A % matvec(xex, b)
  end subroutine make_problem

  ! relative true residual ||A x - b|| / ||b||
  real(dp) function relresid(A, x, b)
    type(csr_matrix), intent(in) :: A
    real(dp)        , intent(in) :: x(:), b(:)
    real(dp), allocatable :: r(:)
    allocate(r(A % num_vertices))
    call A % matvec(x, r)
    relresid = norm2(r - b)/max(norm2(b), tiny(1.0_dp))
  end function relresid

  ! relative solution error ||x - xex|| / ||xex||
  real(dp) function solerr(x, xex)
    real(dp), intent(in) :: x(:), xex(:)
    solerr = norm2(x - xex)/max(norm2(xex), tiny(1.0_dp))
  end function solerr

  !===================================================================!
  ! 1. correctness: all three solvers recover x_exact at a moderate gamma
  !===================================================================!
  subroutine check_correctness(nf)
    integer, intent(inout) :: nf
    integer :: it_nr, it_ne, it_gm
    type(csr_system) :: sys
    type(csr_matrix) :: A
    real(dp), allocatable :: xex(:), b(:), x(:)
    real(dp) :: rr, se
    type(normal_cg)    :: cgnr_solver, cgne_solver
    type(gmres_solver) :: gs

    A = advdiff_csr(0.6_dp)
    call make_problem(A, xex, b)
    allocate(x(A % ncols))

    ! the kernels run on the system contract, never on assembled
    ! entries: wrap the manufactured operator as a system (its
    ! transpose is genuine, so no symmetry claim is needed)
    sys = csr_system(A)

    cgnr_solver = normal_cg(max_it=50000, max_tol=TOL, method=CGNR_METHOD, print_level=0)
    cgne_solver = normal_cg(max_it=50000, max_tol=TOL, method=CGNE_METHOD, print_level=0)
    gs          = gmres_solver(max_it=50000, restart=400, max_tol=TOL, print_level=0)

    write(*,'(a)') " ---- correctness (gamma = 0.6) ----"

    call cgnr_solver % cgnr(sys, b, x, it_nr)
    rr = relresid(A, x, b); se = solerr(x, xex)
    write(*,'(2x,a,es11.3,a,es11.3,a,i0)') "cgnr  resid=", rr, "  err=", se, "  iters=", it_nr
    if (rr .gt. 1.0e-6_dp .or. se .gt. 1.0e-3_dp) nf = nf + 1

    call cgne_solver % cgne(sys, b, x, it_ne)
    rr = relresid(A, x, b); se = solerr(x, xex)
    write(*,'(2x,a,es11.3,a,es11.3,a,i0)') "cgne  resid=", rr, "  err=", se, "  iters=", it_ne
    if (rr .gt. 1.0e-6_dp .or. se .gt. 1.0e-3_dp) nf = nf + 1

    call gs % gmres(sys, b, x, it_gm)
    rr = relresid(A, x, b); se = solerr(x, xex)
    write(*,'(2x,a,es11.3,a,es11.3,a,i0)') "gmres resid=", rr, "  err=", se, "  iters=", it_gm
    if (rr .gt. 1.0e-6_dp .or. se .gt. 1.0e-3_dp) nf = nf + 1
  end subroutine check_correctness

  !===================================================================!
  ! 2. the headline: CGNR (kappa^2) vs GMRES (kappa) iteration counts
  !===================================================================!
  subroutine check_peclet_sweep(nf)
    integer, intent(inout) :: nf
    integer :: it_cgnr
    type(csr_system) :: sys
    real(dp), parameter :: gammas(4) = [0.0_dp, 0.3_dp, 0.6_dp, 0.9_dp]
    type(csr_matrix) :: A
    real(dp), allocatable :: xex(:), b(:), x(:)
    integer :: t
    type(normal_cg)    :: cgnr_solver
    type(gmres_solver) :: gs

    cgnr_solver = normal_cg(max_it=50000, max_tol=TOL, method=CGNR_METHOD, print_level=0)
    gs          = gmres_solver(max_it=50000, restart=400, max_tol=TOL, print_level=0)

    write(*,'(a)') " ---- Peclet sweep: CGNR (A^T A, kappa^2) vs GMRES (A, kappa) ----"
    write(*,'(2x,a8,2x,a10,2x,a10,2x,a8)') "gamma", "cgnr_its", "gmres_its", "ratio"
    do t = 1, size(gammas)
       A = advdiff_csr(gammas(t))
       call make_problem(A, xex, b)
       if (allocated(x)) deallocate(x); allocate(x(A % ncols))
       sys = csr_system(A)

       call cgnr_solver % cgnr(sys, b, x, it_cgnr)
       block
         real(dp) :: rr_c
         integer  :: it_gmres
         rr_c = relresid(A, x, b)
         call gs % gmres(sys, b, x, it_gmres)
         write(*,'(2x,f8.2,2x,i10,2x,i10,2x,f8.1)') gammas(t), it_cgnr, it_gmres, &
              & real(it_cgnr,dp)/real(max(it_gmres,1),dp)
         ! both must converge, and GMRES must take strictly fewer iterations
         if (rr_c .gt. 1.0e-6_dp)              nf = nf + 1
         if (relresid(A, x, b) .gt. 1.0e-6_dp) nf = nf + 1
         if (it_gmres .ge. it_cgnr)            nf = nf + 1
       end block
       deallocate(xex, b)
    end do
  end subroutine check_peclet_sweep

  !===================================================================!
  ! 3. restarted GMRES(m) with a small m still converges
  !===================================================================!
  subroutine check_restart(nf)
    integer, intent(inout) :: nf
    integer :: it_gm
    type(csr_system) :: sys
    type(csr_matrix) :: A
    real(dp), allocatable :: xex(:), b(:), x(:)
    real(dp) :: rr
    type(gmres_solver) :: gs

    A = advdiff_csr(0.6_dp)
    call make_problem(A, xex, b)
    allocate(x(A % ncols))

    gs = gmres_solver(max_it=50000, restart=20, max_tol=TOL, print_level=0)   ! restart every 20
    sys = csr_system(A)
    call gs % gmres(sys, b, x, it_gm)
    rr = relresid(A, x, b)
    write(*,'(a,es11.3,a,i0)') " ---- GMRES(20) restart: resid=", rr, "  iters=", it_gm
    if (rr .gt. 1.0e-6_dp) nf = nf + 1
  end subroutine check_restart

  !===================================================================!
  ! 4. AMG (on the symmetric part) as a right preconditioner for GMRES
  !===================================================================!
  subroutine check_amg_precond(nf)
    integer, intent(inout) :: nf
    type(csr_matrix) :: A, Asym
    type(algebraic_multigrid)        :: M
    real(dp), allocatable :: xex(:), b(:), x(:)
    integer :: it_plain, it_prec
    real(dp) :: rr
    type(gmres_solver) :: gs, gs_prec
    type(csr_system) :: sys

    A    = advdiff_csr(0.6_dp)
    Asym = advdiff_csr(0.0_dp)        ! symmetric part = the SPD laplacian
    call M % setup(Asym)
    call make_problem(A, xex, b)
    allocate(x(A % ncols))

    sys = csr_system(A)
    gs = gmres_solver(max_it=50000, restart=400, max_tol=TOL, print_level=0)
    call gs % gmres(sys, b, x, it_plain)

    gs_prec = gmres_solver(max_it=50000, restart=400, max_tol=TOL, print_level=0, precond=M)
    call gs_prec % gmres(sys, b, x, it_prec)
    rr = relresid(A, x, b)

    write(*,'(a,i0,a,i0,a,es11.3)') " ---- AMG-right-preconditioned GMRES: plain=", &
         & it_plain, "  amg=", it_prec, "  resid=", rr
    if (rr .gt. 1.0e-6_dp)       nf = nf + 1
    if (it_prec .ge. it_plain)   nf = nf + 1     ! preconditioning must help
  end subroutine check_amg_precond

end program test_krylov
