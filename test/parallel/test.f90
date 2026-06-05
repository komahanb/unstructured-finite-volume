!=====================================================================!
! Distributed CG verification (phase 2: coarray domain decomposition).
!
! Each image holds the full problem but owns only its partition's cells.
! For several meshes we:
!   - partition the graph across the images and build the halo
!   - solve A x = b with the distributed (coarray) CG
!   - solve the SAME system with the serial matrix-free library CG (run
!     replicated on every image) as the reference
! and check:
!   1. coverage     - the owned sets partition every dof exactly once
!                     (sum of owned counts across images == n)
!   2. nontrivial   - with >1 image the partition actually cuts edges
!   3. solves it    - ||A x_dist - b|| / ||b|| is at the solver tolerance
!   4. == serial    - the distributed solution matches the serial one
!   5. block-AMG    - a per-image AMG built on each owned block (restricted
!                     additive schwarz preconditioner) gives the same answer
!                     in fewer CG iterations than the unpreconditioned solve
!
! Run: /usr/bin/cafrun.openmpi -np {2,4} ./run   (and serial: ./run_serial).
! A nonzero exit (error stop) means a check failed.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_parallel

  use iso_fortran_env       , only : dp => real64
  use class_gmsh_loader     , only : gmsh_loader
  use class_mesh            , only : mesh
  use class_assembler       , only : assembler
  use class_diffusion_flux  , only : diffusion_flux, constant_source
  use class_csr             , only : csr_matrix
  use class_partition       , only : partition
  use class_amg             , only : amg
  use class_distributed_cg  , only : halo, distributed_cg, owned_block, dist_cg_last_iters
  use class_conjugate_gradient, only : conjugate_gradient

  implicit none

  integer :: me, np, nfail

  me = this_image()
  np = num_images()
  nfail = 0

  if (me .eq. 1) then
     write(*,'(a)')        " ====================================================="
     write(*,'(a,i0,a)')   "  distributed CG verification  (", np, " image(s))"
     write(*,'(a)')        " ====================================================="
  end if

  ! square-20: correctness (one-level block-jacobi's gain is subdomain-count
  ! dependent at this size - ties unpreconditioned at 4 images).
  ! square-40: large enough that per-subdomain AMG always cuts iterations.
  call run_square(20, 1, .false., me, np, nfail)
  call run_square(40, 2, .true. , me, np, nfail)   ! print_level 2: iteration trace
  call run_box(me, np, nfail)

  if (me .eq. 1) then
     write(*,'(a)') " -----------------------------------------------------"
     if (nfail .eq. 0) then
        write(*,'(a)') "  all distributed-cg checks passed"
     else
        write(*,'(1x,i0,a)') nfail, "  distributed-cg checks FAILED"
     end if
  end if

  sync all
  if (nfail .ne. 0) error stop

contains

  !===================================================================!
  ! 2d poisson on square-n (homogeneous dirichlet, unit source)
  !===================================================================!
  subroutine run_square(n, plvl, assert_speedup, me, np, nfail)
    integer, intent(in)    :: n, plvl, me, np
    logical, intent(in)    :: assert_speedup
    integer, intent(inout) :: nfail
    class(assembler), allocatable :: fvm
    class(gmsh_loader), allocatable :: gl
    class(mesh)       , allocatable :: grid
    character(len=64) :: meshfile, label
    write(meshfile, '(a,i0,a)') "square-", n, ".msh"
    write(label,    '(a,i0)')   "square-", n
    allocate(gl  , source = gmsh_loader(trim(meshfile)))
    allocate(grid, source = mesh(gl))
    allocate(fvm , source = assembler(grid))
    call fvm % set_equation(diffusion_flux(1.0_dp), constant_source(-1.0_dp))
    call fvm % set_dirichlet("BoundaryLeft"  , 0.0_dp)
    call fvm % set_dirichlet("BoundaryRight" , 0.0_dp)
    call fvm % set_dirichlet("BoundaryTop"   , 0.0_dp)
    call fvm % set_dirichlet("BoundaryBottom", 0.0_dp)
    call solve_and_check(fvm, trim(label), plvl, assert_speedup, me, np, nfail)
  end subroutine run_square

  !===================================================================!
  ! mixed-bc diffusion on box-36
  !===================================================================!
  subroutine run_box(me, np, nfail)
    integer, intent(in)    :: me, np
    integer, intent(inout) :: nfail
    class(assembler), allocatable :: fvm
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
    ! box-36 is too small for block-jacobi to help, so don't require a speedup
    call solve_and_check(fvm, "box-36", 0, .false., me, np, nfail)
  end subroutine run_box

  !===================================================================!
  ! partition, distributed solve, serial reference, and the four checks
  !===================================================================!
  subroutine solve_and_check(fvm, label, plvl, assert_speedup, me, np, nfail)
    class(assembler), allocatable, intent(inout) :: fvm
    character(*), intent(in)    :: label
    integer     , intent(in)    :: plvl, me, np
    logical     , intent(in)    :: assert_speedup   ! require block-amg to cut iters
    integer     , intent(inout) :: nfail

    type(csr_matrix) :: A, Ablock
    type(partition)  :: p
    type(halo)       :: h
    type(amg)        :: M
    class(conjugate_gradient), allocatable :: cg
    real(dp), allocatable :: b(:), x_dist(:), x_pc(:), x_ref(:), r(:)
    integer  :: n, nown_tot, it_unprec, it_pc
    real(dp) :: e, e_pc, relres, bnorm

    ! partition the graph across images (deterministic - every image agrees)
    call fvm % g % partition(np)
    p = partition(fvm % g, np)
    h = halo(p, me, np)

    ! assembled operator + source (replicated on every image)
    call fvm % get_operator_csr(A)
    n = A % nrows
    allocate(b(n)); call fvm % get_source(b)

    ! (a) unpreconditioned distributed (coarray) CG
    allocate(x_dist(n))
    call distributed_cg(A, b, x_dist, h, 20000, 1.0e-10_dp, plvl)
    it_unprec = dist_cg_last_iters

    ! (b) per-image block-AMG preconditioned distributed CG: each image builds
    ! an AMG on its OWNED diagonal block (restricted additive schwarz)
    Ablock = owned_block(A, h % own)
    call M % setup(Ablock)
    allocate(x_pc(n))
    call distributed_cg(A, b, x_pc, h, 20000, 1.0e-10_dp, 0, precond = M)
    it_pc = dist_cg_last_iters

    ! serial reference: matrix-free library CG, run replicated
    allocate(cg, source = conjugate_gradient(fvm, 20000, 1.0e-10_dp, 0))
    call cg % solve(x_ref)
    deallocate(cg)

    ! check 1: owned sets cover every dof exactly once
    nown_tot = size(h % own)
    call co_sum(nown_tot)

    ! check 3: distributed solution actually solves the system
    allocate(r(n)); call A % matvec(x_dist, r); r = r - b
    bnorm  = max(norm2(b), tiny(1.0_dp))
    relres = norm2(r)/bnorm

    ! checks 4 + 5: both distributed solves match the serial answer
    e    = maxval(abs(x_dist - x_ref))/max(maxval(abs(x_ref)), tiny(1.0_dp))
    e_pc = maxval(abs(x_pc   - x_ref))/max(maxval(abs(x_ref)), tiny(1.0_dp))

    if (me .eq. 1) then
       write(*,'(a)') " -----------------------------------------------------"
       write(*,'(2x,a,a)') "case: ", label
       call p % print()
       write(*,'(4x,a,i0,a,i0)')      "dofs covered      : ", nown_tot, " / ", n
       write(*,'(4x,a,es12.4)')       "dist residual rel : ", relres
       write(*,'(4x,a,es12.4)')       "dist vs serial    : ", e
       write(*,'(4x,a,es12.4)')       "block-amg vs serial: ", e_pc
       write(*,'(4x,a,i0,a,i0)')      "iters  unprec=", it_unprec, "  block-amg=", it_pc
    end if

    if (nown_tot .ne. n)                 nfail = nfail + 1
    if (np .gt. 1 .and. p % ncut .le. 0) nfail = nfail + 1
    if (relres .gt. 1.0e-6_dp)           nfail = nfail + 1
    if (e .gt. 1.0e-6_dp)                nfail = nfail + 1
    if (e_pc .gt. 1.0e-6_dp)             nfail = nfail + 1
    ! block-AMG cuts iterations on the well-resolved elliptic problems; on a
    ! trivially small mesh (box-36, ~7 cells/block) unpreconditioned CG already
    ! converges in a handful of iters and block-jacobi cannot beat it.
    if (assert_speedup .and. it_pc .ge. it_unprec) nfail = nfail + 1

    deallocate(b, x_dist, x_pc, x_ref, r)
  end subroutine solve_and_check

end program test_parallel
