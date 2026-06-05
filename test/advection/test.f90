!=====================================================================!
! Advection (and advection-diffusion) verification.
!
! The advective flux F = v q is integrated through the same flux seam as
! diffusion, but it depends on the face VALUE of q (not its gradient), so
! central differencing gives a skew-symmetric contribution and the operator
! becomes NON-symmetric - the first such operator in the framework, and the
! reason class_gmres exists.
!
! Checks:
!   1. nonsymmetry  - the assembled operator is symmetric with v = 0 and
!                     genuinely nonsymmetric with v /= 0 (A x /= A^T x).
!   2. fidelity     - the assembled csr equals the matrix-free operator with
!                     advection on (the get_operator_csr / get_jvp gate).
!   3. exact 1D     - steady advection-diffusion  v q' - kappa q'' = 0 on the
!                     unit square (q=0 left, q=1 right, zero-flux top/bottom)
!                     has the exact boundary-layer solution
!                       q(x) = (e^{a x} - 1)/(e^{a} - 1),  a = v/kappa.
!                     solved with GMRES, it converges at 2nd order - which
!                     validates the advection sign/magnitude AND the advective
!                     boundary closure (the right wall has q=1, vn /= 0).
!   4. gmres vs cg  - GMRES solves the nonsymmetric system in a handful of
!                     iterations; plain CG (correct only for SPD) is wildly
!                     inefficient here (and not guaranteed to converge at all
!                     for stronger advection).
!   5. upwinding    - at high cell-Peclet central differencing oscillates (the
!                     solution overshoots the physical [0,1] range); the upwind
!                     scheme stays monotone and bounded.
!
! A nonzero exit (error stop) means a check failed.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_advection

  use iso_fortran_env        , only : dp => real64
  use class_gmsh_loader      , only : gmsh_loader
  use class_mesh             , only : mesh
  use class_assembler        , only : assembler, CONVECTION_UPWIND
  use class_advection_flux   , only : advection_diffusion_flux
  use class_diffusion_flux   , only : constant_source
  use class_csr              , only : csr_matrix
  use class_gmres            , only : gmres, gmres_last_iters
  use class_conjugate_gradient, only : conjugate_gradient, cg_last_iters

  implicit none

  real(dp), parameter :: VX = 2.0_dp, KAPPA = 1.0_dp   ! global Peclet = 2
  integer :: nfail
  nfail = 0

  call check_nonsymmetry(nfail)
  call check_fidelity(nfail)
  call check_exact_1d(nfail)
  call check_gmres_vs_cg(nfail)
  call check_upwind(nfail)

  write(*,'(a)') " ============================================="
  if (nfail .eq. 0) then
     write(*,'(a)') "  all advection checks passed"
  else
     write(*,'(1x,i0,a)') nfail, "  advection checks FAILED"
     error stop
  end if

contains

  !===================================================================!
  ! advection-diffusion assembler on square-n: v=(vx,0), kappa; q=0 on the
  ! left wall, q=1 on the right, zero-flux (neumann) top/bottom -> the
  ! solution is 1d in x.
  !===================================================================!
  subroutine make_advdiff(n, vx, kappa, fvm)
    integer , intent(in) :: n
    real(dp), intent(in) :: vx, kappa
    class(assembler), allocatable, intent(out) :: fvm
    class(gmsh_loader), allocatable :: gl
    class(mesh)       , allocatable :: grid
    character(len=64) :: meshfile
    write(meshfile, '(a,i0,a)') "square-", n, ".msh"
    allocate(gl  , source = gmsh_loader(trim(meshfile)))
    allocate(grid, source = mesh(gl))
    allocate(fvm , source = assembler(grid))
    call fvm % set_equation(advection_diffusion_flux([vx, 0.0_dp, 0.0_dp], kappa), &
         &                  constant_source(0.0_dp))
    call fvm % set_dirichlet("BoundaryLeft"  , 0.0_dp)
    call fvm % set_dirichlet("BoundaryRight" , 1.0_dp)
    call fvm % set_neumann  ("BoundaryTop"   , 0.0_dp)
    call fvm % set_neumann  ("BoundaryBottom", 0.0_dp)
  end subroutine make_advdiff

  ! domain x-extent from the boundary face centres (robust to the mesh box)
  subroutine x_extent(grid, xlo, xhi)
    class(mesh), intent(in)  :: grid
    real(dp)   , intent(out) :: xlo, xhi
    integer :: f
    xlo = huge(1.0_dp); xhi = -huge(1.0_dp)
    do f = 1, grid % num_faces
       if (grid % num_face_cells(f) .eq. 1) then
          xlo = min(xlo, grid % face_centers(1, f))
          xhi = max(xhi, grid % face_centers(1, f))
       end if
    end do
  end subroutine x_extent

  !===================================================================!
  ! 1. operator nonsymmetry (v/=0) vs symmetry (v=0)
  !===================================================================!
  subroutine check_nonsymmetry(nf)
    integer, intent(inout) :: nf
    class(assembler), allocatable :: fvm
    type(csr_matrix) :: A
    real(dp), allocatable :: w(:), Aw(:), Atw(:)
    real(dp) :: s_adv, s_diff
    integer  :: i, m

    ! advection on
    call make_advdiff(20, VX, KAPPA, fvm)
    call fvm % get_operator_csr(A)
    m = A % nrows
    allocate(w(m), Aw(m), Atw(m))
    do i = 1, m
       w(i) = sin(0.3_dp*real(i,dp)) + 0.2_dp
    end do
    call A % matvec(w, Aw); call A % matvec_transpose(w, Atw)
    s_adv = maxval(abs(Aw - Atw))/max(maxval(abs(Aw)), tiny(1.0_dp))
    deallocate(fvm)

    ! advection off (v=0) -> pure diffusion -> symmetric
    call make_advdiff(20, 0.0_dp, KAPPA, fvm)
    call fvm % get_operator_csr(A)
    call A % matvec(w, Aw); call A % matvec_transpose(w, Atw)
    s_diff = maxval(abs(Aw - Atw))/max(maxval(abs(Aw)), tiny(1.0_dp))

    write(*,'(a)') " ---- operator symmetry ----"
    write(*,'(2x,a,es11.3)') "||A w - A^T w||/||A w||,  v=0 : ", s_diff
    write(*,'(2x,a,es11.3)') "||A w - A^T w||/||A w||,  v/=0: ", s_adv
    if (s_diff .gt. 1.0e-10_dp) nf = nf + 1   ! diffusion must be symmetric
    if (s_adv  .lt. 1.0e-3_dp ) nf = nf + 1   ! advection must break symmetry
  end subroutine check_nonsymmetry

  !===================================================================!
  ! 2. assembled csr == matrix-free operator, with advection on
  !===================================================================!
  subroutine check_fidelity(nf)
    integer, intent(inout) :: nf
    class(assembler), allocatable :: fvm
    type(csr_matrix) :: A
    real(dp), allocatable :: q(:), y_csr(:), y_mf(:)
    real(dp) :: e
    integer  :: i, m

    call make_advdiff(20, VX, KAPPA, fvm)
    call fvm % get_operator_csr(A)
    m = A % nrows
    allocate(q(m), y_csr(m), y_mf(m))
    do i = 1, m
       q(i) = cos(0.5_dp*real(i,dp)) + 0.4_dp*real(i,dp)
    end do
    call A % matvec(q, y_csr)
    call fvm % get_jacobian_vector_product(y_mf, q)
    e = maxval(abs(y_csr - y_mf))/max(maxval(abs(y_mf)), tiny(1.0_dp))
    write(*,'(a,es11.3)') " ---- csr vs matrix-free (advection): ", e
    if (e .gt. 1.0e-10_dp) nf = nf + 1
  end subroutine check_fidelity

  !===================================================================!
  ! 3. exact 1d advection-diffusion boundary layer, 2nd-order convergence
  !===================================================================!
  subroutine check_exact_1d(nf)
    integer, intent(inout) :: nf
    integer, parameter :: ns(3) = [10, 20, 40]
    class(assembler), allocatable :: fvm
    type(csr_matrix) :: A
    real(dp), allocatable :: b(:), x(:), qex(:)
    real(dp) :: xlo, xhi, a_pe, xc, err(3), order
    integer  :: k, m, c

    write(*,'(a)') " ---- exact 1d advection-diffusion (Pe=2), GMRES ----"
    write(*,'(2x,a6,2x,a12,2x,a8)') "n", "L2 error", "order"
    do k = 1, 3
       call make_advdiff(ns(k), VX, KAPPA, fvm)
       call fvm % get_operator_csr(A)
       m = A % nrows
       allocate(b(m), x(m), qex(m))
       call fvm % get_source(b)
       call gmres(A, b, x, 20000, 200, 1.0e-12_dp, 0)

       ! exact q(x) = (e^{a(x-xlo)} - 1)/(e^{a(xhi-xlo)} - 1),  a = vx/kappa
       call x_extent(fvm % grid, xlo, xhi)
       a_pe = VX/KAPPA
       do c = 1, fvm % grid % num_cells
          xc = fvm % grid % cell_centers(1, c)
          qex(c) = (exp(a_pe*(xc - xlo)) - 1.0_dp)/(exp(a_pe*(xhi - xlo)) - 1.0_dp)
       end do

       err(k) = sqrt(sum(fvm % grid % cell_volumes*(x - qex)**2) &
            &        / sum(fvm % grid % cell_volumes))
       if (k .eq. 1) then
          write(*,'(2x,i6,2x,es12.4,2x,a8)') ns(k), err(k), "  -"
       else
          order = log(err(k-1)/err(k))/log(2.0_dp)
          write(*,'(2x,i6,2x,es12.4,2x,f8.3)') ns(k), err(k), order
       end if
       deallocate(b, x, qex, fvm)
    end do

    order = log(err(2)/err(3))/log(2.0_dp)
    if (err(3) .gt. 5.0e-3_dp) nf = nf + 1            ! converged to the exact
    if (order  .lt. 1.7_dp   ) nf = nf + 1            ! ~2nd order
  end subroutine check_exact_1d

  !===================================================================!
  ! 4. GMRES solves the nonsymmetric system; plain CG does not
  !===================================================================!
  subroutine check_gmres_vs_cg(nf)
    integer, intent(inout) :: nf
    class(assembler)         , allocatable :: fvm
    class(conjugate_gradient), allocatable :: cg
    type(csr_matrix) :: A
    real(dp), allocatable :: b(:), x_g(:), x_c(:), r(:)
    real(dp) :: res_g, res_c, bnorm
    integer  :: m

    call make_advdiff(20, VX, KAPPA, fvm)
    call fvm % get_operator_csr(A)
    m = A % nrows
    allocate(b(m), x_g(m), x_c(m), r(m))
    call fvm % get_source(b)
    bnorm = max(norm2(b), tiny(1.0_dp))

    ! GMRES on the assembled operator
    call gmres(A, b, x_g, 20000, 200, 1.0e-10_dp, 0)
    call A % matvec(x_g, r); res_g = norm2(r - b)/bnorm

    ! plain CG (matrix-free) on the same nonsymmetric operator
    allocate(cg, source = conjugate_gradient(fvm, 5000, 1.0e-10_dp, 0))
    call cg % solve(x_c)
    call A % matvec(x_c, r); res_c = norm2(r - b)/bnorm

    write(*,'(a)') " ---- GMRES vs CG efficiency on the nonsymmetric operator ----"
    write(*,'(2x,a,es11.3,a,i0)') "gmres rel res = ", res_g, "   iters=", gmres_last_iters
    write(*,'(2x,a,es11.3,a,i0)') "cg    rel res = ", res_c, "   iters=", cg_last_iters
    if (res_g .gt. 1.0e-6_dp) nf = nf + 1                         ! GMRES must solve it
    if (cg_last_iters .lt. 10*gmres_last_iters) nf = nf + 1       ! CG far less efficient
  end subroutine check_gmres_vs_cg

  !===================================================================!
  ! 5. upwinding: at high cell-Peclet central differencing oscillates out of
  ! the physical [0,1] range; upwind stays monotone/bounded
  !===================================================================!
  subroutine check_upwind(nf)
    integer, intent(inout) :: nf
    class(assembler), allocatable :: fvm
    type(csr_matrix) :: A
    real(dp), allocatable :: b(:), x(:)
    real(dp) :: minc, maxc, minu, maxu
    integer  :: m

    ! cell-Peclet ~ vx*h/kappa = 80/20 = 4 on square-20 (central oscillates)
    call make_advdiff(20, 80.0_dp, 1.0_dp, fvm)

    ! central (default)
    call fvm % get_operator_csr(A)
    m = A % nrows
    allocate(b(m), x(m))
    call fvm % get_source(b)
    call gmres(A, b, x, 50000, 300, 1.0e-9_dp, 0)
    minc = minval(x); maxc = maxval(x)

    ! upwind
    call fvm % set_convection_scheme(CONVECTION_UPWIND)
    call fvm % get_operator_csr(A)
    call fvm % get_source(b)
    call gmres(A, b, x, 50000, 300, 1.0e-9_dp, 0)
    minu = minval(x); maxu = maxval(x)

    write(*,'(a)') " ---- high-Peclet (cell-Pe~4): central vs upwind range ----"
    write(*,'(2x,a,es11.3,a,es11.3)') "central  min=", minc, "  max=", maxc
    write(*,'(2x,a,es11.3,a,es11.3)') "upwind   min=", minu, "  max=", maxu
    if (minc .gt. -1.0e-3_dp .and. maxc .lt. 1.0_dp + 1.0e-3_dp) nf = nf + 1  ! central must oscillate
    if (minu .lt. -1.0e-6_dp .or.  maxu .gt. 1.0_dp + 1.0e-6_dp) nf = nf + 1  ! upwind must stay bounded
  end subroutine check_upwind

end program test_advection
