!=====================================================================!
! Partitioned-system verification (coarray domain decomposition).
!
! The solver is the ordinary conjugate gradient - the SAME class the
! serial reference uses. Parallelism lives on the system side: the
! partitioned assembler sums inner products over its owned dofs and
! reduces across images, and computes only its owned rows of each
! jacobian-vector product before the result is assembled.
!
! For several meshes we solve with cg on the partitioned system and
! with cg on an identical serial system (run replicated on every image)
! and check:
!   1. coverage     - the owned sets partition every dof exactly once
!                     (sum of owned counts across images == n)
!   2. nontrivial   - with >1 image the partition actually cuts edges
!   3. solves it    - ||A x - b|| / ||b|| is at the solver tolerance
!   4. == serial    - the partitioned solution matches the serial one
!   5. block-AMG    - a per-image AMG built on each owned block
!                     (additive Schwarz preconditioner) gives the same
!                     answer in fewer CG iterations than unpreconditioned
!   6. RCB quality  - recursive coordinate bisection cuts fewer edges (a
!                     smaller halo) than the BFS partitioner on the squares
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
  use class_partitioned_assembler, only : partitioned_assembler, block_preconditioner
  use class_diffusion_flux  , only : diffusion_flux, constant_source
  use class_csr             , only : csr_matrix
  use class_algebraic_multigrid, only : algebraic_multigrid
  use class_conjugate_gradient, only : conjugate_gradient

  implicit none

  integer :: me, np, nfail

  me = this_image()
  np = num_images()
  nfail = 0

  if (me .eq. 1) then
     write(*,'(a)')        " ====================================================="
     write(*,'(a,i0,a)')   "  partitioned-system cg verification  (", np, " image(s))"
     write(*,'(a)')        " ====================================================="
  end if

  ! square-20: correctness (one-level block-jacobi's gain is subdomain-count
  ! dependent at this size - ties unpreconditioned at 4 images).
  ! square-40: large enough that per-subdomain AMG always cuts iterations.
  call run_square(20, 1, .false., .true., me, np, nfail)
  call run_square(40, 2, .true. , .true., me, np, nfail)   ! print_level 2: iteration trace
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
  subroutine run_square(n, plvl, assert_speedup, assert_cut, me, np, nfail)
    integer, intent(in)    :: n, plvl, me, np
    logical, intent(in)    :: assert_speedup, assert_cut
    integer, intent(inout) :: nfail
    type(partitioned_assembler), allocatable :: fvmp
    class(assembler)           , allocatable :: fvms
    class(gmsh_loader), allocatable :: gl
    class(mesh)       , allocatable :: grid
    character(len=64) :: meshfile, label
    write(meshfile, '(a,i0,a)') "square-", n, ".msh"
    write(label,    '(a,i0)')   "square-", n
    allocate(gl  , source = gmsh_loader(trim(meshfile)))
    allocate(grid, source = mesh(gl))
    allocate(fvmp, source = partitioned_assembler(grid))
    allocate(fvms, source = assembler(grid))
    call square_bc(fvmp); call square_bc(fvms)
    call fvmp % setup_partition()
    call solve_and_check(fvmp, fvms, trim(label), plvl, assert_speedup, assert_cut, me, np, nfail)
  end subroutine run_square

  subroutine square_bc(fvm)
    class(assembler), intent(inout) :: fvm
    call fvm % set_equation(diffusion_flux(1.0_dp), constant_source(-1.0_dp))
    call fvm % set_dirichlet("BoundaryLeft"  , 0.0_dp)
    call fvm % set_dirichlet("BoundaryRight" , 0.0_dp)
    call fvm % set_dirichlet("BoundaryTop"   , 0.0_dp)
    call fvm % set_dirichlet("BoundaryBottom", 0.0_dp)
  end subroutine square_bc

  !===================================================================!
  ! mixed-bc diffusion on box-36
  !===================================================================!
  subroutine run_box(me, np, nfail)
    integer, intent(in)    :: me, np
    integer, intent(inout) :: nfail
    type(partitioned_assembler), allocatable :: fvmp
    class(assembler)           , allocatable :: fvms
    class(gmsh_loader), allocatable :: gl
    class(mesh)       , allocatable :: grid
    allocate(gl  , source = gmsh_loader("box-36.msh"))
    allocate(grid, source = mesh(gl))
    allocate(fvmp, source = partitioned_assembler(grid))
    allocate(fvms, source = assembler(grid))
    call box_bc(fvmp); call box_bc(fvms)
    call fvmp % setup_partition()
    ! box-36 is too small for block-jacobi to help or for a reliable cut win,
    ! so don't assert speedup or the edge-cut comparison - just correctness
    call solve_and_check(fvmp, fvms, "box-36", 0, .false., .false., me, np, nfail)
  end subroutine run_box

  subroutine box_bc(fvm)
    class(assembler), intent(inout) :: fvm
    call fvm % set_equation(diffusion_flux(2.0_dp), constant_source(1.0_dp))
    call fvm % set_dirichlet("front" , 5.0_dp)
    call fvm % set_dirichlet("back"  , 0.0_dp)
    call fvm % set_neumann  ("top"   , 0.0_dp)
    call fvm % set_neumann  ("bottom", 0.0_dp)
    call fvm % set_dirichlet("left"  , 1.0_dp)
    call fvm % set_dirichlet("right" , 3.0_dp)
  end subroutine box_bc

  !===================================================================!
  ! solve on the partitioned system, solve on the serial reference,
  ! and run the six checks
  !===================================================================!
  subroutine solve_and_check(fvmp, fvms, label, plvl, assert_speedup, assert_cut, me, np, nfail)
    type(partitioned_assembler), allocatable, intent(inout) :: fvmp
    class(assembler)           , allocatable, intent(inout) :: fvms
    character(*), intent(in)    :: label
    integer     , intent(in)    :: plvl, me, np
    logical     , intent(in)    :: assert_speedup   ! require block-amg to cut iters
    logical     , intent(in)    :: assert_cut       ! require RCB cut <= BFS cut
    integer     , intent(inout) :: nfail

    type(csr_matrix) :: A, Ablock
    type(mesh) :: gp
    type(algebraic_multigrid)              :: M
    class(conjugate_gradient), allocatable :: cg
    real(dp), allocatable :: b(:), x_dist(:), x_pc(:), x_ref(:), r(:)
    integer  :: n, nown_tot, it_unprec, it_pc, bfs_cut, rcb_cut
    real(dp) :: e, e_pc, relres, bnorm

    ! partition the graph: BFS (placeholder) vs RCB (geometric); compare the
    ! edge cut. both are deterministic, so every image computes the identical
    ! partition and agrees without communication.
    gp      = fvmp % grid
    call gp % partition(np)
    bfs_cut = gp % edge_cut()
    gp      = fvmp % grid
    call gp % partition_rcb(fvmp % grid % cell_centers, np)
    rcb_cut = gp % edge_cut()

    ! assembled operator + source (replicated on every image)
    call fvmp % get_operator_csr(A)
    n = A % nrows
    allocate(b(n)); call fvmp % get_source(b)

    ! (a) unpreconditioned cg on the partitioned system
    allocate(cg, source = conjugate_gradient(20000, 1.0e-10_dp, plvl))
    call cg % solve(fvmp, x_dist)
    it_unprec = cg % last_inner_iters
    deallocate(cg)

    ! (b) per-image block-AMG preconditioned cg: each image builds an AMG
    ! on its OWNED diagonal block (additive Schwarz preconditioner)
    Ablock = A % principal_submatrix(fvmp % own)
    call M % setup(Ablock)
    allocate(cg, source = conjugate_gradient(20000, 1.0e-10_dp, 0, &
         & precond = block_preconditioner(M, fvmp % own)))
    call cg % solve(fvmp, x_pc)
    it_pc = cg % last_inner_iters
    deallocate(cg)

    ! serial reference: the same cg on the serial system, run replicated
    allocate(cg, source = conjugate_gradient(20000, 1.0e-10_dp, 0))
    call cg % solve(fvms, x_ref)
    deallocate(cg)

    ! check 1: owned sets cover every dof exactly once
    nown_tot = size(fvmp % own)
    call co_sum(nown_tot)

    ! check 3: partitioned solution actually solves the system
    allocate(r(n)); call A % matvec(x_dist, r); r = r - b
    bnorm  = max(norm2(b), tiny(1.0_dp))
    relres = norm2(r)/bnorm

    ! checks 4 + 5: both partitioned solves match the serial answer
    e    = maxval(abs(x_dist - x_ref))/max(maxval(abs(x_ref)), tiny(1.0_dp))
    e_pc = maxval(abs(x_pc   - x_ref))/max(maxval(abs(x_ref)), tiny(1.0_dp))

    if (me .eq. 1) then
       write(*,'(a)') " -----------------------------------------------------"
       write(*,'(2x,a,a)') "case: ", label
       call gp % print_partition()
       write(*,'(4x,a,i0,a,i0)')      "dofs covered      : ", nown_tot, " / ", n
       write(*,'(4x,a,i0,a,i0)')      "edge cut  bfs=", bfs_cut, "  rcb=", rcb_cut
       write(*,'(4x,a,es12.4)')       "dist residual rel : ", relres
       write(*,'(4x,a,es12.4)')       "dist vs serial    : ", e
       write(*,'(4x,a,es12.4)')       "block-amg vs serial: ", e_pc
       write(*,'(4x,a,i0,a,i0)')      "iters  unprec=", it_unprec, "  block-amg=", it_pc
    end if

    if (nown_tot .ne. n)                 nfail = nfail + 1
    if (np .gt. 1 .and. gp % ncut .le. 0) nfail = nfail + 1
    if (relres .gt. 1.0e-6_dp)           nfail = nfail + 1
    if (e .gt. 1.0e-6_dp)                nfail = nfail + 1
    if (e_pc .gt. 1.0e-6_dp)             nfail = nfail + 1
    ! block-AMG cuts iterations on the well-resolved elliptic problems; on a
    ! trivially small mesh (box-36, ~7 cells/block) unpreconditioned CG already
    ! converges in a handful of iters and block-jacobi cannot beat it.
    if (assert_speedup .and. it_pc .ge. it_unprec) nfail = nfail + 1
    ! RCB should cut no more edges (smaller/equal halo) than BFS
    if (np .gt. 1 .and. assert_cut .and. rcb_cut .gt. bfs_cut) nfail = nfail + 1

    deallocate(b, x_dist, x_pc, x_ref, r)
  end subroutine solve_and_check

end program test_parallel
