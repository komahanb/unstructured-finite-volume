!=====================================================================!
! AMG verification suite (phase 1: serial SA-AMG preconditioned CG).
!
!   1. CSR fidelity      - the assembled sparse operator equals the
!                          matrix-free operator (matvec) and the dense
!                          jacobian (entrywise). THE gate.
!   2. M^-1 SPD          - one V-cycle is symmetric positive definite
!                          (required for the CG preconditioner).
!   3. same solution     - PCG-AMG and plain CG reach the same answer.
!   4. fewer iterations  - PCG-AMG takes markedly fewer CG iterations.
!   5. h-independence    - on square-10/20/40/80, plain-CG iterations grow
!                          with the mesh while AMG iterations stay ~flat
!                          (the point of multigrid).
!
! A nonzero exit (error stop) means a check failed.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_amg

  use iso_fortran_env       , only : dp => real64
  use class_gmsh_loader     , only : gmsh_loader
  use class_mesh            , only : mesh
  use class_assembler       , only : assembler
  use class_diffusion_flux  , only : diffusion_flux, constant_source
  use class_csr             , only : csr_matrix
  use class_algebraic_multigrid, only : algebraic_multigrid
  use class_conjugate_gradient, only : conjugate_gradient

  implicit none

  integer :: nfail
  nfail = 0

  call check_csr_fidelity(nfail)
  call check_amg_spd(nfail)
  call check_same_solution_and_iters(nfail)
  call check_h_independence(nfail)

  write(*,*) "============================================="
  if (nfail .eq. 0) then
     write(*,*) "all amg checks passed"
  else
     write(*,*) nfail, " amg checks FAILED"
     error stop
  end if

contains

  !===================================================================!
  ! Mixed-bc diffusion assembler on box-36 (for the fidelity gate)
  !===================================================================!
  subroutine make_box(fvm)
    class(assembler), allocatable, intent(out) :: fvm
    class(gmsh_loader), allocatable :: gl
    class(mesh)       , allocatable :: grid
    allocate(gl  , source = gmsh_loader("box-36.msh"))
    allocate(grid, source = mesh(gl))
    allocate(fvm , source = assembler(grid))
    call fvm % set_equation(diffusion_flux(2.0_dp), constant_source(1.0_dp))
    call fvm % set_dirichlet("front" , 5.0_dp)
    call fvm % set_dirichlet("back"  , 0.0_dp)
    call fvm % set_neumann  ("top"   , 0.0_dp)
    call fvm % set_neumann  ("bottom", 0.0_dp)
    call fvm % set_dirichlet("left"  , 1.0_dp)
    call fvm % set_dirichlet("right" , 3.0_dp)
  end subroutine make_box

  !===================================================================!
  ! 2d poisson on square-n (homogeneous dirichlet, unit source) - the
  ! canonical SPD elliptic problem for the AMG checks
  !===================================================================!
  subroutine make_square(n, fvm)
    integer, intent(in) :: n
    class(assembler), allocatable, intent(out) :: fvm
    class(gmsh_loader), allocatable :: gl
    class(mesh)       , allocatable :: grid
    character(len=64) :: meshfile
    write(meshfile, '(a,i0,a)') "square-", n, ".msh"
    allocate(gl  , source = gmsh_loader(trim(meshfile)))
    allocate(grid, source = mesh(gl))
    allocate(fvm , source = assembler(grid))
    call fvm % set_equation(diffusion_flux(1.0_dp), constant_source(-1.0_dp))
    call fvm % set_dirichlet("BoundaryLeft"  , 0.0_dp)
    call fvm % set_dirichlet("BoundaryRight" , 0.0_dp)
    call fvm % set_dirichlet("BoundaryTop"   , 0.0_dp)
    call fvm % set_dirichlet("BoundaryBottom", 0.0_dp)
  end subroutine make_square

  !===================================================================!
  ! 1. CSR fidelity vs the matrix-free operator and the dense jacobian
  !===================================================================!
  subroutine check_csr_fidelity(nf)
    integer, intent(inout) :: nf
    class(assembler), allocatable :: fvm
    type(csr_matrix)              :: A
    real(dp), allocatable :: q(:), y_csr(:), y_mf(:), Adense(:,:)
    real(dp) :: e
    integer  :: n, i, j

    call make_box(fvm)
    call fvm % get_operator_csr(A)
    n = A % nrows
    allocate(q(n), y_csr(n), y_mf(n))
    do i = 1, n
       q(i) = sin(real(i,dp)*0.7_dp) + 0.3_dp*real(i,dp)
    end do

    call A % matvec(q, y_csr)
    call fvm % get_jacobian_vector_product(y_mf, q)
    e = maxval(abs(y_csr - y_mf))
    write(*,'(a,es12.4)') " csr matvec vs matrix-free : ", e
    if (e .gt. 1.0e-10_dp) nf = nf + 1

    call fvm % get_jacobian(Adense)
    e = 0.0_dp
    do i = 1, n
       do j = 1, n
          e = max(e, abs(A % get_entry(i,j) - Adense(i,j)))
       end do
    end do
    write(*,'(a,es12.4)') " csr vs dense entrywise    : ", e
    if (e .gt. 1.0e-10_dp) nf = nf + 1
  end subroutine check_csr_fidelity

  !===================================================================!
  ! 2. M^-1 (one V-cycle) is symmetric positive definite
  !===================================================================!
  subroutine check_amg_spd(nf)
    integer, intent(inout) :: nf
    class(assembler), allocatable :: fvm
    type(csr_matrix)              :: A
    type(algebraic_multigrid)                     :: M
    real(dp), allocatable :: u(:), v(:), Mu(:), Mv(:)
    real(dp) :: sym, pd_u, pd_v
    integer  :: n, i

    call make_square(40, fvm)
    call fvm % get_operator_csr(A)
    call M % setup(A)
    n = A % nrows
    allocate(u(n), v(n), Mu(n), Mv(n))
    do i = 1, n
       u(i) = sin(real(i,dp)*0.3_dp)
       v(i) = cos(real(i,dp)*0.5_dp) + 0.2_dp
    end do
    call M % apply(u, Mu)
    call M % apply(v, Mv)
    ! M^-1 must be SYMMETRIC (so PCG is a valid CG) and DEFINITE with a
    ! consistent sign. The diffusion operator A is negative definite
    ! (Aq = + discrete laplacian), so M^-1 ~ A^-1 is negative definite -
    ! that is correct, not a defect; only consistent definiteness matters.
    sym  = abs(dot_product(v, Mu) - dot_product(u, Mv))/max(abs(dot_product(u, Mv)), 1.0e-30_dp)
    pd_u = dot_product(u, Mu)
    pd_v = dot_product(v, Mv)
    write(*,'(a,es12.4,a,es12.4)') " amg M^-1 symmetry rel err : ", sym, "   v^T M^-1 v = ", pd_v
    if (sym .gt. 1.0e-8_dp)    nf = nf + 1   ! symmetric
    if (pd_u*pd_v .le. 0.0_dp) nf = nf + 1   ! same sign => definite
  end subroutine check_amg_spd

  !===================================================================!
  ! 3 + 4. PCG-AMG reaches the same solution as plain CG, in fewer iters
  !===================================================================!
  subroutine check_same_solution_and_iters(nf)
    integer, intent(inout) :: nf
    class(assembler)         , allocatable :: fvm
    class(conjugate_gradient), allocatable :: cg
    type(csr_matrix)                       :: A
    type(algebraic_multigrid)                              :: M
    real(dp), allocatable :: x_cg(:), x_amg(:)
    integer :: iters_cg, iters_amg
    real(dp) :: e

    call make_square(40, fvm)

    ! plain CG
    allocate(cg, source = conjugate_gradient(5000, 1.0e-10_dp, 0))
    call cg % solve(fvm, x_cg)
    iters_cg = cg % last_inner_iters
    deallocate(cg)

    ! PCG-AMG (same operator, same tolerance)
    call fvm % get_operator_csr(A)
    call M % setup(A)
    allocate(cg, source = conjugate_gradient(5000, 1.0e-10_dp, 0, precond = M))
    call cg % solve(fvm, x_amg)
    iters_amg = cg % last_inner_iters
    deallocate(cg)

    e = maxval(abs(x_cg - x_amg))/max(maxval(abs(x_cg)), 1.0e-30_dp)
    write(*,'(a,es12.4)')   " pcg-amg vs plain cg sol   : ", e
    write(*,'(a,i6,a,i6)')  " iterations   plain cg = ", iters_cg, "    pcg-amg = ", iters_amg
    if (e .gt. 1.0e-7_dp)            nf = nf + 1
    if (iters_amg .ge. iters_cg)     nf = nf + 1
  end subroutine check_same_solution_and_iters

  !===================================================================!
  ! 5. h-independence: cg iterations grow with the mesh, amg stays flat
  !===================================================================!
  subroutine check_h_independence(nf)
    integer, intent(inout) :: nf
    integer, parameter :: ns(4) = [10, 20, 40, 80]
    class(assembler)         , allocatable :: fvm
    class(conjugate_gradient), allocatable :: cg
    type(csr_matrix)                       :: A
    type(algebraic_multigrid)                              :: M
    real(dp), allocatable :: x(:)
    integer :: ic(4), ia(4), k

    write(*,'(a)') " h-independence (2d poisson):"
    write(*,'(2x,a6,2x,a10,2x,a10)') "n", "cg_iters", "amg_iters"
    do k = 1, 4
       call make_square(ns(k), fvm)

       allocate(cg, source = conjugate_gradient(20000, 1.0e-8_dp, 0))
       call cg % solve(fvm, x)
       ic(k) = cg % last_inner_iters
       deallocate(cg)

       call fvm % get_operator_csr(A)
       call M % setup(A)
       allocate(cg, source = conjugate_gradient(20000, 1.0e-8_dp, 0, precond = M))
       call cg % solve(fvm, x)
       ia(k) = cg % last_inner_iters
       deallocate(cg)

       write(*,'(2x,i6,2x,i10,2x,i10)') ns(k), ic(k), ia(k)
       deallocate(fvm)
    end do

    ! h-independence: plain CG grows strongly with refinement while AMG
    ! stays nearly flat (sub-doubling over an 8x mesh refinement) and
    ! clearly beats CG at the finest mesh.
    if (ic(4) .le. 3*ic(1))    nf = nf + 1   ! plain cg grows a lot
    if (ia(4) .ge. 2*ia(1))    nf = nf + 1   ! amg iterations stay ~flat
    if (maxval(ia) .ge. ic(4)) nf = nf + 1   ! amg beats cg at the finest mesh
  end subroutine check_h_independence

end program test_amg
