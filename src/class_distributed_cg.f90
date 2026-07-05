!=====================================================================!
! Distributed-memory conjugate gradient over a coarray halo.
!
! Domain decomposition, "replicated matrix / partitioned work" flavour:
! every image holds the full assembled operator (csr) and right-hand side
! but is responsible only for the rows of the cells it OWNS (class_partition).
! The CG algorithm is unchanged from the serial one; only two things become
! collective:
!
!   matvec  w = A p : each image needs p on its ghost cells (the off-image
!                     columns its owned rows touch). A coarray halo exchange
!                     pulls just that boundary layer from the owning images -
!                     the interface, not the whole vector.
!   dot/norm        : a local sum over owned entries, then co_sum.
!
! With one image (serial -fcoarray=single build) there are no ghosts, the
! exchange is a no-op and co_sum is the identity, so this reduces exactly to
! the serial CG. Real runs use /usr/bin/cafrun.openmpi -np {2,4}.
!
! The matvec operator A is matrix-free in the serial Krylov path; here we use
! the assembled csr (class_assembler % get_operator_csr) so a partitioned row
! loop is trivial. (True per-subdomain assembly is a later, memory-distributed
! refinement; the halo + reductions are the parallel essence and are real.)
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_distributed_cg

  use iso_fortran_env         , only : dp => REAL64
  use class_csr               , only : csr_matrix
  use class_graph             , only : graph
  use interface_assembler     , only : assembler
  ! the concrete spatial assembler is needed (via select type) for the mesh
  ! geometry the RCB partition uses; the rest goes through the abstract base
  use class_assembler         , only : spatial_assembler => assembler
  use interface_linear_solver , only : preconditioner, linear_solver

  implicit none

  private
  public :: halo
  public :: distributed_cg_solver
  public :: dist_cg_last_iters        ! iterations of the last distributed solve

  ! Inner iterations of the most recent distributed_cg (diagnostics / tests).
  integer :: dist_cg_last_iters = 0

  !===================================================================!
  ! Per-image halo: which cells I own, and the off-image cells I must
  ! pull each matvec, grouped by the image that owns them (so the
  ! exchange is one coarray reference per neighbour).
  !===================================================================!

  type :: halo
     integer :: n  = 0                  ! global system size (dofs = cells)
     integer :: me = 1                  ! this image
     integer :: np = 1                  ! number of images
     integer, allocatable :: own(:)     ! global indices I own (my rows)
     integer, allocatable :: recv_ptr(:)! (np+1) csr over source images
     integer, allocatable :: recv_idx(:)! global indices I receive, by owner
   contains
     procedure :: exchange              ! fill ghost slots of a coarray vector
     procedure :: ddot                  ! distributed dot over owned entries
  end type halo

  interface halo
     module procedure make_halo
  end interface halo

  !-------------------------------------------------------------------!
  ! linear_solver wrapper: partitions the assembler's graph across the
  ! coarray images, builds the halo, assembles the operator + rhs, and
  ! runs distributed_cg. With one image (serial -fcoarray=single build)
  ! the partition is trivial and this reduces to ordinary CG. SPD
  ! operators only (diffusion), like the serial conjugate_gradient.
  !-------------------------------------------------------------------!

  type, extends(linear_solver) :: distributed_cg_solver
     ! stateless w.r.t. the system: solve takes the assembler as an argument
     class(preconditioner), allocatable :: precond     ! optional per-image block precond
   contains
     procedure :: solve
     procedure :: distributed_cg
  end type distributed_cg_solver

  interface distributed_cg_solver
     module procedure construct_dist
  end interface distributed_cg_solver

contains

  !===================================================================!
  ! Build this image's halo from the (replicated) partitioned graph.
  !===================================================================!

  pure type(halo) function make_halo(g, me, np) result(h)

    type(graph), intent(in) :: g
    integer    , intent(in) :: me, np

    integer, allocatable :: gh(:), cnt(:), ptr(:)
    integer :: i, o, gc

    h % n  = g % num_vertices
    h % me = me
    h % np = np

    ! my owned rows and my ghosts, sliced straight from the graph's csr
    ! (explicit allocate avoids a spurious realloc-on-assign warning)
    allocate(h % own(g % n_owned(me)))
    h % own = g % own_list(g % own_ptr(me) : g % own_ptr(me+1) - 1)

    ! ghosts grouped by owning image (counting sort over part_of)
    allocate(gh(g % n_ghosts(me)))
    gh = g % gh_list(g % gh_ptr(me) : g % gh_ptr(me+1) - 1)
    allocate(cnt(np)); cnt = 0
    do i = 1, size(gh)
       o = g % part_of(gh(i))
       cnt(o) = cnt(o) + 1
    end do
    allocate(h % recv_ptr(np+1))
    h % recv_ptr(1) = 1
    do o = 1, np
       h % recv_ptr(o+1) = h % recv_ptr(o) + cnt(o)
    end do
    allocate(h % recv_idx(size(gh)))
    allocate(ptr(np)); ptr = h % recv_ptr(1:np)
    do i = 1, size(gh)
       gc = gh(i)
       o = g % part_of(gc)
       h % recv_idx(ptr(o)) = gc
       ptr(o) = ptr(o) + 1
    end do

  end function make_halo

  !===================================================================!
  ! Halo exchange: pull my ghost cells' values from the images that own
  ! them. v is a coarray; on each owning image v(g) holds the authoritative
  ! value for the cells it owns. One strided remote GET per neighbour.
  !===================================================================!

  impure subroutine exchange(this, v)

    class(halo), intent(in)    :: this
    real(dp)   , intent(inout) :: v(:)[*]

    integer :: o, j, g

    sync all                       ! every image has written its owned slots
    do o = 1, this % np
       if (o .eq. this % me) cycle
       ! scalar coindexed GETs (the most robust coarray op; vector-subscripted
       ! remote refs mis-transfer under this OpenCoarrays/OpenMPI)
       do j = this % recv_ptr(o), this % recv_ptr(o+1) - 1
          g = this % recv_idx(j)
          v(g) = v(g)[o]
       end do
    end do
    sync all                       ! all GETs done before owners overwrite
  end subroutine exchange

  !===================================================================!
  ! Distributed dot product: sum over MY owned entries only (so each dof is
  ! counted once across images), then co_sum to the global value.
  !===================================================================!

  impure real(dp) function ddot(this, a, b)

    class(halo), intent(in) :: this
    real(dp)   , intent(in) :: a(:), b(:)

    integer :: i, k

    ddot = 0.0_dp
    do i = 1, size(this % own)
       k = this % own(i)
       ddot = ddot + a(k)*b(k)
    end do
    call co_sum(ddot)
  end function ddot

  !===================================================================!
  ! Distributed CG: solve A x = b with A symmetric (negative-)definite,
  ! partitioned by h % own. x is returned assembled (correct on every image).
  !===================================================================!

  impure subroutine distributed_cg(this, A, b, x, h)

    class(distributed_cg_solver), intent(in)  :: this
    type(csr_matrix)            , intent(in)  :: A
    real(dp)                    , intent(in)  :: b(:)
    real(dp)                    , intent(out) :: x(:)
    type(halo)                  , intent(in)  :: h

    real(dp), allocatable :: p(:)[:]            ! search dir (needs ghosts)
    real(dp), allocatable :: r(:), w(:), z(:)   ! full-length work vectors
    real(dp), allocatable :: r_loc(:), z_loc(:) ! owned-block scratch for precond
    real(dp) :: rho, rho_old, pAp, alpha, beta, bnorm, rnorm, tol
    integer  :: n, nown, iter, i, k

    n    = h % n
    nown = size(h % own)
    allocate(p(n)[*]); allocate(r(n), w(n), z(n), r_loc(nown), z_loc(nown))
    p = 0.0_dp; r = 0.0_dp; w = 0.0_dp; z = 0.0_dp; x = 0.0_dp

    ! initial residual r = b - A x with x = 0  ->  r = b (on owned rows)
    do i = 1, size(h % own)
       k = h % own(i); r(k) = b(k)
    end do

    bnorm = sqrt(h % ddot(b, b))
    if (bnorm .eq. 0.0_dp) bnorm = 1.0_dp
    rho   = h % ddot(r, r)
    rnorm = sqrt(rho)
    tol   = rnorm/bnorm

    call apply_precond(r, z)                   ! z(owned) = M^-1 r(owned) (or r)
    rho = h % ddot(r, z)

    if (this % print_level .gt. 0 .and. h % me .eq. 1) &
         & write(*,'(1x,a,i0,a,es12.5)') "dist-cg: images=", h % np, "  tol0=", tol

    iter = 1
    do while (tol .gt. this % max_tol .and. iter .lt. this % max_it)

       ! search direction p = z + beta p   (owned entries)
       if (iter .eq. 1) then
          do i = 1, size(h % own); k = h % own(i); p(k) = z(k); end do
       else
          beta = rho/rho_old
          do i = 1, size(h % own); k = h % own(i); p(k) = z(k) + beta*p(k); end do
       end if

       ! w = A p  (halo-exchange p, then local rows)
       call h % exchange(p)
       call A % matvec_rows(p, w, h % own)

       pAp   = h % ddot(p, w)
       alpha = rho/pAp

       do i = 1, size(h % own)
          k = h % own(i)
          x(k) = x(k) + alpha*p(k)
          r(k) = r(k) - alpha*w(k)
       end do

       rnorm = sqrt(h % ddot(r, r))
       tol   = rnorm/bnorm

       call apply_precond(r, z)
       rho_old = rho
       rho     = h % ddot(r, z)

       if (this % print_level .gt. 1 .and. h % me .eq. 1) &
            & write(*,'(1x,a,i5,a,es12.5)') "  iter ", iter, "  rel res ", tol

       iter = iter + 1
    end do
    dist_cg_last_iters = iter - 1

    if (this % print_level .gt. 0 .and. h % me .eq. 1) &
         & write(*,'(1x,a,i0,a,es12.5)') "dist-cg: converged in ", iter-1, &
         & " iters, rel res ", tol

    ! assemble the full solution on every image: zero non-owned, co_sum
    block
      real(dp), allocatable :: xfull(:)
      logical, allocatable  :: ismine(:)
      allocate(xfull(n), ismine(n)); xfull = 0.0_dp; ismine = .false.
      do i = 1, size(h % own); ismine(h % own(i)) = .true.; end do
      where (ismine) xfull = x
      call co_sum(xfull)
      x = xfull
    end block

  contains

    !-----------------------------------------------------------------!
    ! z(owned) = M^-1 r(owned) via this image's block preconditioner,
    ! or z(owned) = r(owned) when none is attached (plain CG). The block
    ! solve couples only owned dofs; the global krylov couples subdomains
    ! through the matvec halo (restricted additive schwarz / block jacobi).
    !-----------------------------------------------------------------!
    impure subroutine apply_precond(rr, zz)
      real(dp), intent(in)    :: rr(:)
      real(dp), intent(inout) :: zz(:)
      integer :: ii, kk
      if (allocated(this % precond)) then
         do ii = 1, nown; r_loc(ii) = rr(h % own(ii)); end do
         call this % precond % apply(r_loc, z_loc)
         do ii = 1, nown; zz(h % own(ii)) = z_loc(ii); end do
      else
         do ii = 1, nown; kk = h % own(ii); zz(kk) = rr(kk); end do
      end if
    end subroutine apply_precond

  end subroutine distributed_cg

  !===================================================================!
  ! Constructor for the distributed_cg linear-solver wrapper
  !===================================================================!

  pure type(distributed_cg_solver) function construct_dist(max_it, max_tol, &
       & print_level, precond) result(this)

    integer              , intent(in)           :: max_it
    real(dp)             , intent(in)           :: max_tol
    integer              , intent(in), optional :: print_level
    class(preconditioner), intent(in), optional :: precond

    this % max_it  = max_it
    this % max_tol = max_tol

    if (present(print_level)) this % print_level = print_level
    if (present(precond))     allocate(this % precond, source = precond)

  end function construct_dist

  !===================================================================!
  ! Partition across images, build the halo, assemble, run distributed CG.
  !===================================================================!

  impure subroutine solve(this, system, x, mode)

    class(distributed_cg_solver), intent(in)       :: this
    class(assembler)            , intent(in)       :: system
    real(dp), allocatable       , intent(out)      :: x(:)
    integer              , intent(in), optional :: mode  ! FORWARD (default) / REVERSE

    type(csr_matrix)      :: A
    type(graph)           :: g
    type(halo)            :: h
    real(dp), allocatable :: b(:)
    integer               :: n, me, np

    me = this_image()
    np = num_images()

    ! geometric (RCB) partition of the graph (needs the mesh geometry, so
    ! reach the concrete spatial assembler), then this image's halo
    select type (system)
    type is (spatial_assembler)
       g = system % g % partition_rcb(system % grid % cell_centers, np)
    class default
       error stop "distributed_cg: needs a spatial assembler (mesh geometry)"
    end select
    h = halo(g, me, np)

    call system % get_operator_csr(A)
    n = A % nrows
    allocate(b(n), x(n))
    call system % get_source(b)

    call this % distributed_cg(A, b, x, h)

  end subroutine solve

end module class_distributed_cg
