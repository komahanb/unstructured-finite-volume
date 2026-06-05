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
  use class_partition         , only : partition
  use class_graph             , only : graph
  use class_assembler         , only : assembler
  use interface_linear_solver , only : preconditioner, linear_solver

  implicit none

  private
  public :: halo, distributed_cg, owned_block
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
     class(assembler)     , allocatable :: FVAssembler
     class(preconditioner), allocatable :: precond     ! optional per-image block precond
     integer                            :: print_level = 0
   contains
     procedure :: solve => distributed_cg_solve
  end type distributed_cg_solver

  interface distributed_cg_solver
     module procedure construct_dist
  end interface distributed_cg_solver

contains

  !===================================================================!
  ! Build this image's halo from the (replicated) partition.
  !===================================================================!

  type(halo) function make_halo(p, me, np) result(h)

    type(partition), intent(in) :: p
    integer        , intent(in) :: me, np

    integer, allocatable :: gh(:), cnt(:), ptr(:)
    integer :: i, o, g

    h % n  = p % ncells
    h % me = me
    h % np = np

    ! my owned rows and my ghosts, sliced straight from the partition csr
    ! (explicit allocate avoids a spurious realloc-on-assign warning)
    allocate(h % own(p % n_owned(me)))
    h % own = p % own_list(p % own_ptr(me) : p % own_ptr(me+1) - 1)

    ! ghosts grouped by owning image (counting sort over part_of)
    allocate(gh(p % n_ghosts(me)))
    gh = p % gh_list(p % gh_ptr(me) : p % gh_ptr(me+1) - 1)
    allocate(cnt(np)); cnt = 0
    do i = 1, size(gh)
       o = p % part_of(gh(i))
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
       g = gh(i)
       o = p % part_of(g)
       h % recv_idx(ptr(o)) = g
       ptr(o) = ptr(o) + 1
    end do

  end function make_halo

  !===================================================================!
  ! Halo exchange: pull my ghost cells' values from the images that own
  ! them. v is a coarray; on each owning image v(g) holds the authoritative
  ! value for the cells it owns. One strided remote GET per neighbour.
  !===================================================================!

  subroutine exchange(this, v)

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

  real(dp) function ddot(this, a, b)

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

  subroutine distributed_cg(A, b, x, h, max_it, max_tol, print_level, precond)

    type(csr_matrix)     , intent(in)           :: A
    real(dp)             , intent(in)           :: b(:)
    real(dp)             , intent(out)          :: x(:)
    type(halo)           , intent(in)           :: h
    integer              , intent(in)           :: max_it
    real(dp)             , intent(in)           :: max_tol
    integer              , intent(in)           :: print_level
    class(preconditioner), intent(in), optional :: precond  ! per-image block precond

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

    if (print_level .gt. 0 .and. h % me .eq. 1) &
         & write(*,'(1x,a,i0,a,es12.5)') "dist-cg: images=", h % np, "  tol0=", tol

    iter = 1
    do while (tol .gt. max_tol .and. iter .lt. max_it)

       ! search direction p = z + beta p   (owned entries)
       if (iter .eq. 1) then
          do i = 1, size(h % own); k = h % own(i); p(k) = z(k); end do
       else
          beta = rho/rho_old
          do i = 1, size(h % own); k = h % own(i); p(k) = z(k) + beta*p(k); end do
       end if

       ! w = A p  (halo-exchange p, then local rows)
       call h % exchange(p)
       call local_matvec(A, p, w, h % own)

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

       if (print_level .gt. 1 .and. h % me .eq. 1) &
            & write(*,'(1x,a,i5,a,es12.5)') "  iter ", iter, "  rel res ", tol

       iter = iter + 1
    end do
    dist_cg_last_iters = iter - 1

    if (print_level .gt. 0 .and. h % me .eq. 1) &
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
    subroutine apply_precond(rr, zz)
      real(dp), intent(in)    :: rr(:)
      real(dp), intent(inout) :: zz(:)
      integer :: ii, kk
      if (present(precond)) then
         do ii = 1, nown; r_loc(ii) = rr(h % own(ii)); end do
         call precond % apply(r_loc, z_loc)
         do ii = 1, nown; zz(h % own(ii)) = z_loc(ii); end do
      else
         do ii = 1, nown; kk = h % own(ii); zz(kk) = rr(kk); end do
      end if
    end subroutine apply_precond

  end subroutine distributed_cg

  !===================================================================!
  ! Local matvec: w(row) = sum_j A(row,j) v(j) for the OWNED rows only.
  ! v must already carry valid owned + ghost values (post halo exchange).
  !===================================================================!

  subroutine local_matvec(A, v, w, own)

    type(csr_matrix), intent(in)    :: A
    real(dp)        , intent(in)    :: v(:)
    real(dp)        , intent(inout) :: w(:)
    integer         , intent(in)    :: own(:)

    integer  :: i, row, jj
    real(dp) :: s

    do i = 1, size(own)
       row = own(i)
       s = 0.0_dp
       do jj = A % row_ptr(row), A % row_ptr(row+1)-1
          s = s + A % vals(jj) * v(A % col_idx(jj))
       end do
       w(row) = s
    end do
  end subroutine local_matvec

  !===================================================================!
  ! Owned-owned diagonal block of A in local numbering (1..n_own),
  ! dropping the off-image (ghost) columns. This is the local subdomain
  ! operator on which each image builds its block preconditioner (e.g.
  ! class_amg). With one image it is the whole matrix, so the block
  ! preconditioner becomes the global one.
  !===================================================================!

  type(csr_matrix) function owned_block(A, own) result(B)

    type(csr_matrix), intent(in) :: A
    integer         , intent(in) :: own(:)

    integer , allocatable :: loc(:), row_ptr(:), col_idx(:)
    real(dp), allocatable :: vals(:)
    integer :: nown, il, row, jj, jl, pos, nnz

    nown = size(own)
    allocate(loc(A % nrows)); loc = 0
    do il = 1, nown
       loc(own(il)) = il
    end do

    ! symbolic: count owned-owned entries per local row
    allocate(row_ptr(nown+1)); row_ptr(1) = 1
    do il = 1, nown
       row = own(il)
       pos = 0
       do jj = A % row_ptr(row), A % row_ptr(row+1)-1
          if (loc(A % col_idx(jj)) .gt. 0) pos = pos + 1
       end do
       row_ptr(il+1) = row_ptr(il) + pos
    end do
    nnz = row_ptr(nown+1) - 1
    allocate(col_idx(nnz), vals(nnz))

    ! numeric: copy entries, remapping global columns to local indices
    pos = 1
    do il = 1, nown
       row = own(il)
       do jj = A % row_ptr(row), A % row_ptr(row+1)-1
          jl = loc(A % col_idx(jj))
          if (jl .gt. 0) then
             col_idx(pos) = jl
             vals(pos)    = A % vals(jj)
             pos = pos + 1
          end if
       end do
    end do

    B = csr_matrix(nown, nown, row_ptr, col_idx, vals)
  end function owned_block

  !===================================================================!
  ! Constructor for the distributed_cg linear-solver wrapper
  !===================================================================!

  type(distributed_cg_solver) function construct_dist(FVAssembler, max_it, max_tol, print_level, precond) result(this)
    type(assembler)      , intent(in)           :: FVAssembler
    integer              , intent(in)           :: max_it
    real(dp)             , intent(in)           :: max_tol
    integer              , intent(in), optional :: print_level
    class(preconditioner), intent(in), optional :: precond
    allocate(this % FVAssembler, source = FVAssembler)
    this % max_it  = max_it
    this % max_tol = max_tol
    if (present(print_level)) this % print_level = print_level
    if (present(precond))     allocate(this % precond, source = precond)
  end function construct_dist

  !===================================================================!
  ! Partition across images, build the halo, assemble, run distributed CG.
  !===================================================================!

  subroutine distributed_cg_solve(this, x)

    class(distributed_cg_solver), intent(in)  :: this
    real(dp), allocatable       , intent(out) :: x(:)

    type(csr_matrix)      :: A
    type(graph)           :: g
    type(partition)       :: p
    type(halo)            :: h
    real(dp), allocatable :: b(:)
    integer               :: n, me, np

    me = this_image()
    np = num_images()

    ! local copy of the graph (partition_rcb stamps vertex % part)
    g = this % FVAssembler % g
    call g % partition_rcb(this % FVAssembler % grid % cell_centers, np)
    p = partition(g, np)
    h = halo(p, me, np)

    call this % FVAssembler % get_operator_csr(A)
    n = A % nrows
    allocate(b(n), x(n))
    call this % FVAssembler % get_source(b)

    if (allocated(this % precond)) then
       call distributed_cg(A, b, x, h, this % max_it, this % max_tol, &
            &              this % print_level, precond = this % precond)
    else
       call distributed_cg(A, b, x, h, this % max_it, this % max_tol, this % print_level)
    end if

  end subroutine distributed_cg_solve

end module class_distributed_cg
