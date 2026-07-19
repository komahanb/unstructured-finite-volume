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
!   6. geometric         - the geometric squint (rcb on the mesh graph,
!                          quotients down the hierarchy) reaches the same
!                          solution in markedly fewer iterations than
!                          plain CG, staying ~flat over refinement.
!   7. injected cycle    - the cycle is a graph handed in: the W built
!                          as stations and arrows converges; repeated
!                          stations smooth again; a negative station
!                          dips into a refined level below the mesh
!                          (the zoom, built by refine). differential
!                          probes prove each feature steers the trip.
!   8. unstructured      - nothing assumes structure: the same problem
!                          on a crooked triangle mesh, both squints
!                          still beat plain CG to the same answer.
!   9. hierarchy output  - every coarsening level painted onto the
!                          mesh as a cell field through the existing
!                          writer: a quotient is a partition of the
!                          mesh, and a partition is a picture.
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
  use class_geometric_multigrid, only : geometric_multigrid
  use class_stored_graph       , only : stored_digraph
  use class_conjugate_gradient, only : conjugate_gradient
  use class_paraview_writer   , only : paraview_writer
  use class_string            , only : string

  implicit none

  integer :: nfail
  nfail = 0

  call check_csr_fidelity(nfail)
  call check_amg_spd(nfail)
  call check_same_solution_and_iters(nfail)
  call check_h_independence(nfail)
  call check_geometric(nfail)
  call check_cycle_injection(nfail)
  call check_unstructured(nfail)
  call check_hierarchy_output(nfail)

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

  !===================================================================!
  ! 6. geometric multigrid: the spatial squint. Same solution as plain
  ! cg, markedly fewer iterations, ~flat over refinement - the matrix
  ! never chose the aggregates, the coordinates did.
  !===================================================================!

  subroutine check_geometric(nf)
    integer, intent(inout) :: nf
    integer, parameter :: ns(4) = [10, 20, 40, 80]
    class(assembler)         , allocatable :: fvm
    class(conjugate_gradient), allocatable :: cg
    type(csr_matrix)                       :: A
    type(geometric_multigrid)              :: M
    real(dp), allocatable :: x_cg(:), x_geo(:)
    integer  :: ic(4), ig(4), k
    real(dp) :: e

    write(*,'(a)') " geometric multigrid (2d poisson):"
    write(*,'(2x,a6,2x,a10,2x,a10)') "n", "cg_iters", "geo_iters"
    do k = 1, 4
       call make_square(ns(k), fvm)

       allocate(cg, source = conjugate_gradient(20000, 1.0e-8_dp, 0))
       call cg % solve(fvm, x_cg)
       ic(k) = cg % last_inner_iters
       deallocate(cg)

       M = geometric_multigrid(fvm % grid)
       call fvm % get_operator_csr(A)
       call M % setup(A)
       ! re-setup on the same built object must not die - the
       ! re-entrancy contract the shared setup driver promises
       if (k .eq. 1) call M % setup(A)
       allocate(cg, source = conjugate_gradient(20000, 1.0e-8_dp, 0, precond = M))
       call cg % solve(fvm, x_geo)
       ig(k) = cg % last_inner_iters
       deallocate(cg)

       write(*,'(2x,i6,2x,i10,2x,i10)') ns(k), ic(k), ig(k)

       if (k .eq. 3) then
          e = maxval(abs(x_cg - x_geo))/max(maxval(abs(x_cg)), 1.0e-30_dp)
          write(*,'(a,es12.4)') " pcg-geometric vs plain cg : ", e
          if (e .gt. 1.0e-7_dp) nf = nf + 1
       end if

       deallocate(fvm)
    end do

    if (ig(4) .ge. 2*ig(1))    nf = nf + 1   ! geometric iterations stay ~flat
    if (maxval(ig) .ge. ic(4)) nf = nf + 1   ! and beat cg at the finest mesh
  end subroutine check_geometric

  !===================================================================!
  ! 7. the cycle is a graph handed in. Build the W - the deep levels
  ! revisited twice - as stations and arrows, inject it through
  ! apply_cycle, and richardson iteration with one W-trip per step
  ! must reach a manufactured solution.
  !===================================================================!

  subroutine check_cycle_injection(nf)
    integer, intent(inout) :: nf
    class(assembler), allocatable :: fvm
    type(csr_matrix)              :: A
    type(algebraic_multigrid)     :: M
    type(stored_digraph)  :: w_cycle, doubled, refined_trip
    integer , allocatable :: level_visits(:)
    real(dp), allocatable :: x_ref(:), b(:), r(:), z(:)
    integer               :: n, i
    real(dp)              :: bnorm, steer

    call make_square(40, fvm)
    call fvm % get_operator_csr(A)
    call M % setup(A)

    ! the W over the whole hierarchy, written as stations and arrows
    level_visits = w_level_visits(0, M % num_levels() - 1)
    w_cycle = stored_digraph(size(level_visits), &
         & tails   = [(i, i = 1, size(level_visits)-1)], &
         & heads   = [(i+1, i = 1, size(level_visits)-1)], &
         & numbers = level_visits)

    n = A % nrows
    allocate(x_ref(n), b(n), r(n), z(n))
    do i = 1, n
       x_ref(i) = sin(real(i,dp)*0.3_dp) + 0.2_dp
    end do
    call A % matvec(x_ref, b)
    bnorm = norm2(b)

    ! the injected schedule must steer the machinery: one W-trip and
    ! one V-trip from the same residual give different corrections (a
    ! bug that quietly walks the stored V passes every convergence
    ! check, so convergence alone proves nothing about injection)
    call M % apply(b, z)
    call M % apply_cycle(w_cycle, b, r)
    steer = norm2(r - z)/max(norm2(z), 1.0e-30_dp)
    write(*,'(a,es12.4)') " injected w vs built-in v  : ", steer
    if (steer .le. 1.0e-12_dp) nf = nf + 1

    call solve_by_schedule(M, A, w_cycle, b, x_ref, bnorm, &
         & " injected w-cycle  ", nf)

    ! repeated stations smooth again: [0,0,1,1,2,1,0] doubles the
    ! smoothing on the way down and must still converge - and it must
    ! differ from the plain v, which visits the same levels once, so
    ! only the extra smoothing separates them (a bug that collapses
    ! repeats degenerates this schedule to the v and converges anyway)
    doubled = path_schedule(doubled_level_visits(0, M % num_levels() - 1))
    call solve_by_schedule(M, A, doubled, b, x_ref, bnorm, &
         & " doubled smoothing ", nf)
    call M % apply(b, z)
    call M % apply_cycle(doubled, b, r)
    steer = norm2(r - z)/max(norm2(z), 1.0e-30_dp)
    write(*,'(a,es12.4)') " doubled vs plain v        : ", steer
    if (steer .le. 1.0e-12_dp) nf = nf + 1

    ! negative stations are refined levels - the zoom. build one
    ! level below the mesh, and take a trip that dips into it before
    ! running the v: [1, -1, 1, 2, ..., 2, 1]
    call M % refine(1)
    refined_trip = path_schedule([0, -1, v_level_visits(M % num_levels() - 1)])
    call solve_by_schedule(M, A, refined_trip, b, x_ref, bnorm, &
         & " refined excursion ", nf)

    ! and the excursion must matter for what happens AT the refined
    ! level. the control is [0, 0, v...]: same trip, the dip replaced
    ! by staying home - a dead zoom (lift, nothing, project undoes it,
    ! since the interpolation's transpose is its inverse that way)
    ! degenerates the excursion to exactly this control, so only live
    ! refined smoothing separates them
    call M % apply_cycle(path_schedule([0, 0, v_level_visits(M % num_levels() - 1)]), b, z)
    call M % apply_cycle(refined_trip, b, r)
    steer = norm2(r - z)/max(norm2(z), 1.0e-30_dp)
    write(*,'(a,es12.4)') " refined vs stay-home ctrl : ", steer
    if (steer .le. 1.0e-12_dp) nf = nf + 1

  end subroutine check_cycle_injection

  ! the levels the v visits: down the hierarchy and back (levels
  ! count squints from the mesh, 0 the mesh itself)
  pure function v_level_visits(deepest) result(st)
    integer, intent(in)  :: deepest
    integer, allocatable :: st(:)
    integer :: l
    st = [(l, l = 0, deepest), (l, l = deepest-1, 0, -1)]
  end function v_level_visits

  ! richardson iteration, one trip of the given schedule per step,
  ! against a manufactured solution
  subroutine solve_by_schedule(M, A, schedule, b, x_ref, bnorm, label, nf)
    type(algebraic_multigrid), intent(in)    :: M
    type(csr_matrix)         , intent(in)    :: A
    type(stored_digraph)     , intent(in)    :: schedule
    real(dp)                 , intent(in)    :: b(:), x_ref(:), bnorm
    character(len=*)         , intent(in)    :: label
    integer                  , intent(inout) :: nf

    real(dp), allocatable :: x(:), r(:), z(:), Ax(:)
    integer  :: it, ncycles
    real(dp) :: e

    allocate(x(size(b)), r(size(b)), z(size(b)), Ax(size(b)))
    x       = 0.0_dp
    ncycles = 0
    do it = 1, 100
       call A % matvec(x, Ax)
       r = b - Ax
       if (norm2(r)/bnorm .le. 1.0e-10_dp) exit
       ncycles = ncycles + 1
       call M % apply_cycle(schedule, r, z)
       x = x + z
    end do

    e = maxval(abs(x - x_ref))/maxval(abs(x_ref))
    write(*,'(a,a,i0,a,i0,a,es12.4)') label, ": ", &
         & schedule % num_vertices, " stations, ", ncycles, " cycles, error ", e
    if (e .gt. 1.0e-7_dp) nf = nf + 1
    if (ncycles .ge. 100) nf = nf + 1
  end subroutine solve_by_schedule

  ! a straight line of stations from a list of level numbers
  pure function path_schedule(numbers) result(s)
    integer, intent(in)  :: numbers(:)
    type(stored_digraph) :: s
    integer :: i, ns
    ns = size(numbers)
    s = stored_digraph(ns, tails = [(i, i = 1, ns-1)], &
         &                 heads = [(i+1, i = 1, ns-1)], numbers = numbers)
  end function path_schedule

  ! the levels the W visits, written recursively: visit the level, go
  ! deep, come back, go deep again, come back (the l >= nlevel base
  ! case also refuses an unbuilt hierarchy rather than recursing off
  ! the stack)
  pure recursive function w_level_visits(l, nlevel) result(st)
    integer, intent(in) :: l, nlevel
    integer, allocatable :: st(:)
    if (l .ge. nlevel) then
       st = [l]
    else
       st = [l, w_level_visits(l+1, nlevel), l, w_level_visits(l+1, nlevel), l]
    end if
  end function w_level_visits

  ! a v with every downward station doubled - [1,1,2,2,3,2,1] on
  ! three levels: smooth twice on the way down, once on the way up
  pure recursive function doubled_level_visits(l, nlevel) result(st)
    integer, intent(in) :: l, nlevel
    integer, allocatable :: st(:)
    if (l .ge. nlevel) then
       st = [l]
    else
       st = [l, l, doubled_level_visits(l+1, nlevel), l]
    end if
  end function doubled_level_visits


  !===================================================================!
  ! 8. nothing assumes structure. The same poisson problem on a
  ! genuinely crooked mesh - unstructured triangles, no two alike -
  ! and both squints still beat plain cg to the same answer.
  !===================================================================!

  subroutine check_unstructured(nf)
    integer, intent(inout) :: nf
    class(assembler)         , allocatable :: fvm
    class(conjugate_gradient), allocatable :: cg
    type(csr_matrix)                       :: A
    type(algebraic_multigrid)              :: Ma
    type(geometric_multigrid)              :: Mg
    real(dp), allocatable :: x_cg(:), x_a(:), x_g(:)
    integer  :: ic, ia, ig
    real(dp) :: ea, eg

    call make_square_tri(fvm)

    allocate(cg, source = conjugate_gradient(20000, 1.0e-8_dp, 0))
    call cg % solve(fvm, x_cg)
    ic = cg % last_inner_iters
    deallocate(cg)

    call fvm % get_operator_csr(A)

    call Ma % setup(A)
    allocate(cg, source = conjugate_gradient(20000, 1.0e-8_dp, 0, precond = Ma))
    call cg % solve(fvm, x_a)
    ia = cg % last_inner_iters
    deallocate(cg)

    Mg = geometric_multigrid(fvm % grid)
    call Mg % setup(A)
    allocate(cg, source = conjugate_gradient(20000, 1.0e-8_dp, 0, precond = Mg))
    call cg % solve(fvm, x_g)
    ig = cg % last_inner_iters
    deallocate(cg)

    ea = maxval(abs(x_cg - x_a))/max(maxval(abs(x_cg)), 1.0e-30_dp)
    eg = maxval(abs(x_cg - x_g))/max(maxval(abs(x_cg)), 1.0e-30_dp)
    write(*,'(a,i0,a,i0,a,i0)') " unstructured triangles: cg=", ic, &
         & "  amg=", ia, "  geometric=", ig
    if (ea .gt. 1.0e-7_dp) nf = nf + 1
    if (eg .gt. 1.0e-7_dp) nf = nf + 1
    if (ia .ge. ic)        nf = nf + 1
    if (ig .ge. ic)        nf = nf + 1
  end subroutine check_unstructured

  !===================================================================!
  ! 9. the hierarchy written out. Every coarsening level is a
  ! partition of the mesh - each cell's ancestor after l squints -
  ! so the existing writer paints all of them in one call, one field
  ! per level. No coarse mesh file is invented. The painted part
  ! count at each level must equal that level's operator rows: the
  ! picture and the algebra must agree.
  !===================================================================!

  subroutine check_hierarchy_output(nf)

    integer, intent(inout) :: nf

    class(assembler)      , allocatable :: fvm
    class(paraview_writer), allocatable :: painter
    type(string)          , allocatable :: labels(:)
    real(dp)              , allocatable :: fields(:,:)

    type(csr_matrix)          :: A
    type(geometric_multigrid) :: M

    type(algebraic_multigrid) :: M_unconfigured

    integer, allocatable :: fine_ancestors(:)

    character(len=16) :: buf
    logical           :: on_disk, agrees
    integer           :: nl, l, unit, ierr

    call make_square_tri(fvm)

    M = geometric_multigrid(fvm % grid)

    call fvm % get_operator_csr(A)

    ! export is configured at problem setup: name the levels wanted,
    ! positive or negative, and only their ancestry is stored
    call M % export_levels([(l, l = 1, 8), -1])
    call M % setup(A)
    call M % refine(1)

    ! one field per coarsening level: each cell painted by its
    ! ancestor, and the painted part count must equal that level's
    ! operator rows - the picture and the algebra must agree
    nl = M % num_levels() - 1
    allocate(fields(A % nrows, nl), labels(nl))

    agrees = .true.
    do l = 1, nl
       fields(:, l) = real(M % ancestors(l), dp)
       write(buf, '(a,i0)') "level_", l
       labels(l) = string(trim(buf))
       if (nint(maxval(fields(:, l))) .ne. M % levels(l+1) % A % nrows) then
          agrees = .false.
       end if
    end do

    ! a stale picture would satisfy the on-disk check - burn it first
    open(newunit = unit, file = 'hierarchy.vtu', status = 'old', iostat = ierr)
    if (ierr .eq. 0) close(unit, status = 'delete')

    allocate(painter, source = paraview_writer(fvm % grid))
    call painter % write("hierarchy.vtu", fields, labels)

    inquire(file = "hierarchy.vtu", exist = on_disk)

    write(*,'(a,i0,a,i0,a)') " hierarchy painted: ", nl, &
         & " coarsening levels onto ", A % nrows, " cells (hierarchy.vtu)"

    if (.not. agrees)  nf = nf + 1
    if (.not. on_disk) nf = nf + 1

    ! the refined side: each vertex of level -1 knows its mesh cell,
    ! and every cell has exactly its share of children
    fine_ancestors = M % ancestors(-1)
    if (size(fine_ancestors) .ne. M % children_per_vertex*A % nrows) nf = nf + 1
    if (minval(fine_ancestors) .ne. 1 .or. &
         & maxval(fine_ancestors) .ne. A % nrows)                    nf = nf + 1
    if (count(fine_ancestors .eq. 1) .ne. M % children_per_vertex)   nf = nf + 1

    ! storage is on demand: an unconfigured hierarchy keeps nothing
    call M_unconfigured % setup(A)
    if (allocated(M_unconfigured % levels(1) % parts)) nf = nf + 1

  end subroutine check_hierarchy_output

  !===================================================================!
  ! 2d poisson on the unstructured triangle square (same boundary
  ! names as the structured square, so the same problem runs on a
  ! crooked mesh)
  !===================================================================!

  subroutine make_square_tri(fvm)
    class(assembler), allocatable, intent(out) :: fvm
    class(gmsh_loader), allocatable :: gl
    class(mesh)       , allocatable :: grid
    allocate(gl  , source = gmsh_loader("square-tri-40.msh"))
    allocate(grid, source = mesh(gl))
    allocate(fvm , source = assembler(grid))
    call fvm % set_equation(diffusion_flux(1.0_dp), constant_source(-1.0_dp))
    call fvm % set_dirichlet("BoundaryLeft"  , 0.0_dp)
    call fvm % set_dirichlet("BoundaryRight" , 0.0_dp)
    call fvm % set_dirichlet("BoundaryTop"   , 0.0_dp)
    call fvm % set_dirichlet("BoundaryBottom", 0.0_dp)
  end subroutine make_square_tri

end program test_amg
