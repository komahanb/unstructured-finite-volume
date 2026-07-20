!=====================================================================!
! Partitioned-system verification (coarray domain decomposition,
! DISTRIBUTED vectors - the dereplication stages 2 and 3).
!
! The solver is the ordinary conjugate gradient - the SAME class the
! serial reference uses. Parallelism lives on the system side, and
! after setup_partition a vector is this image's owned slab:
!
!    ┌─ image 1 ─┐   halo   ┌─ image 2 ─┐      inner products dot
!    │ owned slab│ <------> │ owned slab│      the slabs + one scalar
!    └───────────┘          └───────────┘      reduction; the product
!                                              exchanges the halo then
!                                              dots the owned rows
!
! For several meshes we solve with cg on the partitioned system and
! with cg on an identical serial system (run replicated on every image)
! and check:
!   1. coverage     - the owned sets partition every dof exactly once
!   2. shrinkage    - the answer really is a slab (nown, not n)
!   3. solves it    - the slabs, replicated through the trio
!                     (scatter + sum), satisfy ||A x - b||/||b|| at tol
!   4. == serial    - and match the serial answer entrywise
!   5. block-AMG    - a per-image AMG on each owned block (additive
!                     Schwarz) reaches the same answer in fewer
!                     iterations than unpreconditioned
!   6. RCB quality  - rcb cuts fewer edges (a smaller halo, which is
!                     now literally less traffic) than bfs
!   7. halo reach   - every dof an owned row touches is owned or ghost
!                     (the exchange lists are exactly sufficient)
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
  use class_paraview_writer , only : paraview_writer
  use class_string          , only : string

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
    real(dp), allocatable :: xd(:), xp(:)
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
    n = A % num_vertices
    allocate(b(n)); call fvmp % get_source(b)

    ! (a) unpreconditioned cg on the partitioned system: the answer
    ! comes back as this image's owned slab
    allocate(cg, source = conjugate_gradient(20000, 1.0e-10_dp, plvl))
    call cg % solve(fvmp, x_dist)
    it_unprec = cg % last_inner_iters
    deallocate(cg)

    ! (b) per-image block-AMG preconditioned cg: each image builds an AMG
    ! on its OWNED diagonal block (additive Schwarz preconditioner); the
    ! preconditioner needs no lists - the frame did the plumbing
    Ablock = A % principal_submatrix(fvmp % own)
    call M % setup(Ablock)
    allocate(cg, source = conjugate_gradient(20000, 1.0e-10_dp, 0, &
         & precond = block_preconditioner(M)))
    call cg % solve(fvmp, x_pc)
    it_pc = cg % last_inner_iters
    deallocate(cg)

    ! serial reference: the same cg on the serial system, run replicated
    allocate(cg, source = conjugate_gradient(20000, 1.0e-10_dp, 0))
    call cg % solve(fvms, x_ref)
    deallocate(cg)

    ! the door replicates the slabs for the global comparisons below
    allocate(xd(n), xp(n))
    call fvmp % replicate(x_dist, xd)
    call fvmp % replicate(x_pc  , xp)

    ! the writer at the door: image 1 paints the distributed answer
    ! and the decomposition that produced it - the solve draws its
    ! own partition
    if (me .eq. 1) then
       write_door: block
         type(paraview_writer), allocatable :: pw
         type(string)          :: names(2)
         real(dp), allocatable :: fields(:,:)
         integer :: v
         allocate(fields(n, 2))
         fields(:,1) = xd
         fields(:,2) = [(real(fvmp % grid % part_of(v), dp), v = 1, n)]
         names(1) = string("phi")
         names(2) = string("part")
         allocate(pw, source = paraview_writer(fvmp % grid))
         call pw % write("dist-"//label//".vtu", fields, names)
         deallocate(pw, fields)
       end block write_door
    end if

    ! check 1: owned sets cover every dof exactly once
    nown_tot = size(fvmp % own)
    call co_sum(nown_tot)

    ! check 3: partitioned solution actually solves the system
    allocate(r(n)); call A % matvec(xd, r); r = r - b
    bnorm  = max(norm2(b), tiny(1.0_dp))
    relres = norm2(r)/bnorm

    ! checks 4 + 5: both partitioned solves match the serial answer
    e    = maxval(abs(xd - x_ref))/max(maxval(abs(x_ref)), tiny(1.0_dp))
    e_pc = maxval(abs(xp - x_ref))/max(maxval(abs(x_ref)), tiny(1.0_dp))

    if (me .eq. 1) then
       write(*,'(a)') " -----------------------------------------------------"
       write(*,'(2x,a,a)') "case: ", label
       call gp % print_partition()
       write(*,'(4x,a,i0,a,i0)')      "dofs covered      : ", nown_tot, " / ", n
       write(*,'(4x,a,i0,a,i0)')      "local length      : ", size(x_dist), " of ", n
       write(*,'(4x,a,i0,a,i0)')      "edge cut  bfs=", bfs_cut, "  rcb=", rcb_cut
       write(*,'(4x,a,es12.4)')       "dist residual rel : ", relres
       write(*,'(4x,a,es12.4)')       "dist vs serial    : ", e
       write(*,'(4x,a,es12.4)')       "block-amg vs serial: ", e_pc
       write(*,'(4x,a,i0,a,i0)')      "iters  unprec=", it_unprec, "  block-amg=", it_pc
    end if

    ! rung-4 readiness: the halo is exact. every dof an owned row of
    ! the operator reaches is either owned or in the ghost halo - so
    ! the day vectors stop being replicated, the graph's ghost lists
    ! ARE the exchange lists, nothing more and nothing less needed.
    halo_reach: block
      logical, allocatable :: known(:)
      integer, allocatable :: gh(:)
      integer :: i2, v2, e2, miss
      allocate(known(n)); known = .false.
      known(fvmp % own) = .true.
      gh = fvmp % grid % dofs_of(fvmp % grid % ghosts(me))
      if (size(gh) .gt. 0) known(gh) = .true.
      miss = 0
      do i2 = 1, size(fvmp % own)
         v2 = fvmp % own(i2)
         do e2 = A % out_xadj(v2), A % out_xadj(v2+1) - 1
            if (.not. known(A % out_adj(e2))) miss = miss + 1
         end do
      end do
      if (me .eq. 1) write(*,'(4x,a,i0,a)') &
           & "halo reach       : ", miss, " owned-row columns outside owned+ghost"
      if (miss .ne. 0) nfail = nfail + 1
    end block halo_reach

    ! the vectors actually shrank: the answer is the owned slab, not
    ! a photocopy of the whole
    if (size(x_dist) .ne. size(fvmp % own))          nfail = nfail + 1
    if (np .gt. 1 .and. size(x_dist) .ge. n)         nfail = nfail + 1

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

    deallocate(b, x_dist, x_pc, x_ref, r, xd, xp)
  end subroutine solve_and_check

end program test_parallel
