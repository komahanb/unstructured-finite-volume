!=====================================================================!
! Abstract multigrid method: the squint made mechanical. Build a coarse
! hierarchy from a fine operator (setup), then use it either as a
! preconditioner (apply: one V-cycle, z = M^-1 r - the deferred apply
! inherited from preconditioner) or as a standalone solver (solve:
! cycle until converged).
!
! Everything mechanical lives here, written once: the level storage,
! the setup driver (coarsen - tentative prolongation over the parts -
! smoothed prolongation - galerkin coarse operator, recursing until the
! coarsest level is small, then a dense LU there), the V-cycle
! (pre-smooth with weighted jacobi, restrict the residual, recurse,
! prolong and correct, post-smooth), and the coarse direct solve.
!
! A concrete kind answers exactly one deferred question: how does this
! level squint - which fine vertices huddle into which coarse part.
! Algebraic multigrid answers by strength (the strength graph's own
! aggregation); geometric multigrid answers by coordinates (recursive
! bisection of the mesh graph). Everything after the squint is shared.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module interface_multigrid

  use iso_fortran_env        , only : dp => REAL64
  use class_csr              , only : csr_matrix
  use interface_linear_solver, only : preconditioner

  implicit none

  private
  public :: multigrid

  !-------------------------------------------------------------------!
  ! One multigrid level
  !-------------------------------------------------------------------!

  type :: multigrid_level
     type(csr_matrix)      :: A           ! operator on this level
     type(csr_matrix)      :: P           ! prolongation (fine x coarse)
     type(csr_matrix)      :: R           ! restriction  = P^T
     real(dp), allocatable :: Dinv(:)     ! 1/diag(A) (smoother)
  end type multigrid_level

  !-------------------------------------------------------------------!
  ! A multigrid is a preconditioner (apply = one cycle) that can also
  ! solve standalone. The hierarchy and every mechanical part live on
  ! this abstract type; a concrete kind only answers how to squint.
  !-------------------------------------------------------------------!

  type, abstract, extends(preconditioner) :: multigrid

     type(multigrid_level), allocatable :: levels(:)   ! 1=finest .. nlevel=coarsest
     integer                            :: nlevel = 0

     ! coarsest level: dense LU factor
     real(dp), allocatable :: Ac_dense(:,:)
     integer , allocatable :: ipiv(:)
     integer               :: ncoarse = 0

     ! parameters of the mechanism
     real(dp) :: jacobi_w    = 0.6667_dp   ! smoother weight (~2/3)
     integer  :: max_levels  = 25
     integer  :: coarse_size = 50
     integer  :: npre        = 1
     integer  :: npost       = 1

   contains

     ! the one deferred question: which fine vertex huddles into which
     ! coarse part on level lev
     procedure(coarsen_interface), deferred :: coarsen

     ! the mechanism, shared
     procedure :: setup       ! build the hierarchy from A
     procedure :: apply       ! z = M^-1 r (one V-cycle)
     procedure :: solve       ! A x = b, cycling to tol
     procedure :: num_levels  ! levels in the hierarchy

     procedure, private :: vcycle
     procedure, private :: smooth
     procedure, private :: spectral_radius

  end type multigrid

  !-------------------------------------------------------------------!
  ! Deferred interface of the squint
  !-------------------------------------------------------------------!

  abstract interface

     ! Stamp each fine vertex of level lev with its coarse part
     ! (agg(v) in 1..naggr). The level's operator is available as
     ! this % levels(lev) % A.
     subroutine coarsen_interface(this, lev, agg, naggr)
       import :: multigrid
       class(multigrid)    , intent(inout) :: this
       integer             , intent(in)    :: lev
       integer, allocatable, intent(out)   :: agg(:)
       integer             , intent(out)   :: naggr
     end subroutine coarsen_interface

  end interface

contains

  !===================================================================!
  ! Build the hierarchy from the fine operator A: at each level ask
  ! the concrete kind for its squint (coarsen), then build the
  ! tentative prolongation over the parts, smooth it, and form the
  ! galerkin coarse operator Ac = R A P - recursing until the coarsest
  ! level is small, where a dense LU waits.
  !===================================================================!

  impure subroutine setup(this, A)

    class(multigrid), intent(inout) :: this
    type(csr_matrix), intent(in)    :: A

    integer               :: lev, n, naggr, info, i
    integer , allocatable :: agg(:)
    real(dp), allocatable :: dinv(:)
    real(dp)              :: omega, rho
    type(csr_matrix)      :: P0, AP0

    ! re-entrant: drop any previously-built hierarchy
    if (allocated(this % levels))   deallocate(this % levels)
    if (allocated(this % Ac_dense)) deallocate(this % Ac_dense)
    if (allocated(this % ipiv))     deallocate(this % ipiv)

    allocate(this % levels(this % max_levels))
    this % levels(1) % A = A

    lev = 1
    coarsening: do

       n = this % levels(lev) % A % nrows
       if (n .le. this % coarse_size .or. lev .ge. this % max_levels) exit coarsening

       ! smoother data on this level
       dinv = 1.0_dp/this % levels(lev) % A % get_diagonal()
       this % levels(lev) % Dinv = dinv

       ! the concrete kind's squint, then the tentative prolongation
       call this % coarsen(lev, agg, naggr)
       if (naggr .ge. n .or. naggr .le. 1) exit coarsening   ! no useful coarsening

       P0 = tentative_prolongation(agg, naggr)

       rho   = this % spectral_radius(lev)
       omega = (4.0_dp/3.0_dp)/rho

       ! smoothed prolongation  P = (I - omega Dinv A) P0 = P0 - omega Dinv (A P0).
       ! the recipe is multigrid's; it is composed from the general csr operations.
       AP0 = this % levels(lev) % A % matmat(P0)
       call AP0 % scale_rows(dinv)
       this % levels(lev) % P = P0 % add(1.0_dp, -omega, AP0)
       this % levels(lev) % R = this % levels(lev) % P % transpose()

       ! galerkin coarse operator  Ac = R A P
       this % levels(lev+1) % A = this % levels(lev) % R % matmat( &
            & this % levels(lev) % A % matmat(this % levels(lev) % P))

       lev = lev + 1

    end do coarsening

    this % nlevel  = lev
    this % ncoarse = this % levels(lev) % A % nrows

    ! coarsest: dense LU (with a tiny diagonal shift if singular, e.g. a
    ! pure-neumann constant null space). self-contained LU so the library
    ! carries no LAPACK dependency.
    call this % levels(lev) % A % to_dense(this % Ac_dense)
    allocate(this % ipiv(this % ncoarse))
    call dense_lu_factor(this % Ac_dense, this % ncoarse, this % ipiv, info)
    if (info .ne. 0) then
       call this % levels(lev) % A % to_dense(this % Ac_dense)
       do i = 1, this % ncoarse
          this % Ac_dense(i,i) = this % Ac_dense(i,i) + 1.0e-10_dp
       end do
       call dense_lu_factor(this % Ac_dense, this % ncoarse, this % ipiv, info)
       if (info .ne. 0) error stop "multigrid: coarse factorization failed"
    end if

  end subroutine setup

  !===================================================================!
  ! z = M^-1 r  (one V-cycle from the finest level)
  !===================================================================!

  pure subroutine apply(this, r, z)
    class(multigrid), intent(in)  :: this
    real(dp)        , intent(in)  :: r(:)
    real(dp)        , intent(out) :: z(:)
    z = 0.0_dp
    call this % vcycle(1, r, z)
  end subroutine apply

  !===================================================================!
  ! Standalone solve A x = b by cycling until ||b - A x||/||b|| < max_tol
  ! or max_it cycles. iters returns the cycles taken.
  !===================================================================!

  pure subroutine solve(this, b, x, max_tol, max_it, iters)

    class(multigrid), intent(in)    :: this
    real(dp)        , intent(in)    :: b(:)
    real(dp)        , intent(inout) :: x(:)
    real(dp)        , intent(in)    :: max_tol
    integer         , intent(in)    :: max_it
    integer         , intent(out)   :: iters

    real(dp), allocatable :: r(:), z(:), Ax(:)
    real(dp)              :: bnorm

    associate(A => this % levels(1) % A)

      allocate(r(A % nrows), z(A % nrows), Ax(A % nrows))

      bnorm = norm2(b); if (bnorm .eq. 0.0_dp) bnorm = 1.0_dp

      do iters = 1, max_it
         call A % matvec(x, Ax)
         r = b - Ax
         if (norm2(r)/bnorm .le. max_tol) exit
         call this % apply(r, z)
         x = x + z
      end do

    end associate

  end subroutine solve

  !===================================================================!
  ! Number of levels in the hierarchy (1 = a single grid)
  !===================================================================!

  pure integer function num_levels(this)
    class(multigrid), intent(in) :: this
    num_levels = this % nlevel
  end function num_levels

  !===================================================================!
  ! Recursive V-cycle worker
  !===================================================================!

  pure recursive subroutine vcycle(this, lev, b, x)

    class(multigrid), intent(in)    :: this
    integer         , intent(in)    :: lev
    real(dp)        , intent(in)    :: b(:)
    real(dp)        , intent(inout) :: x(:)

    real(dp), allocatable :: Ax(:), res(:), rc(:), ec(:), ef(:)
    integer               :: n, nc

    if (lev .eq. this % nlevel) then
       ! coarsest: direct solve  x = Ac^-1 b
       x = b
       call dense_lu_solve(this % Ac_dense, this % ncoarse, this % ipiv, x)
       return
    end if

    associate(L => this % levels(lev))

      n  = L % A % nrows
      nc = L % R % nrows
      allocate(Ax(n), res(n), rc(nc), ec(nc), ef(n))

      ! pre-smooth
      call this % smooth(lev, b, x, this % npre)

      ! residual and restriction to the coarse level
      call L % A % matvec(x, Ax)
      res = b - Ax
      call L % R % matvec(res, rc)

      ! coarse-grid correction (recurse)
      ec = 0.0_dp
      call this % vcycle(lev+1, rc, ec)
      call L % P % matvec(ec, ef)
      x = x + ef

      ! post-smooth
      call this % smooth(lev, b, x, this % npost)

    end associate

  end subroutine vcycle

  !===================================================================!
  ! Weighted-jacobi smoother on level lev (symmetric -> keeps M^-1 SPD):
  !   x <- x + w Dinv (b - A x)
  !===================================================================!

  pure subroutine smooth(this, lev, b, x, nsweep)
    class(multigrid), intent(in)    :: this
    integer         , intent(in)    :: lev
    real(dp)        , intent(in)    :: b(:)
    real(dp)        , intent(inout) :: x(:)
    integer         , intent(in)    :: nsweep
    real(dp), allocatable :: Ax(:)
    integer :: s
    associate(A => this % levels(lev) % A, Dinv => this % levels(lev) % Dinv, w => this % jacobi_w)
      allocate(Ax(A % nrows))
      do s = 1, nsweep
         call A % matvec(x, Ax)
         x = x + w*Dinv*(b - Ax)
      end do
    end associate
  end subroutine smooth

  !===================================================================!
  ! Spectral radius of (Dinv A) on level lev by bounded power iteration,
  ! floored by the gershgorin bound max_i sum_j |Dinv A|_ij - used for the
  ! smoothed-prolongation weight omega = (4/3)/rho
  !===================================================================!

  pure real(dp) function spectral_radius(this, lev) result(rho)
    class(multigrid), intent(in) :: this
    integer         , intent(in) :: lev
    real(dp), allocatable :: v(:), w(:)
    real(dp) :: nv, gers, rowsum
    integer  :: n, it, i, k

    associate(A => this % levels(lev) % A, Dinv => this % levels(lev) % Dinv)

      n = A % nrows
      allocate(v(n), w(n))

      ! non-constant deterministic start (the constant is near-null)
      do i = 1, n
         v(i) = sin(real(i,dp)*0.9_dp) + 0.1_dp
      end do
      v = v/norm2(v)

      rho = 0.0_dp
      do it = 1, 15
         call A % matvec(v, w)
         w  = Dinv*w
         nv = norm2(w)
         if (nv .le. 0.0_dp) exit
         rho = nv
         v   = w/nv
      end do

      ! gershgorin floor
      gers = 0.0_dp
      do i = 1, n
         rowsum = 0.0_dp
         do k = A % row_ptr(i), A % row_ptr(i+1) - 1
            rowsum = rowsum + abs(Dinv(i)*A % vals(k))
         end do
         gers = max(gers, rowsum)
      end do
      rho = max(rho, gers)
      if (rho .le. 0.0_dp) rho = 1.0_dp

    end associate

  end function spectral_radius

  !===================================================================!
  ! Tentative prolongation P0 (n x naggr): piecewise constant over the
  ! parts, column-normalized so the constant near-null space is
  ! represented exactly (orthonormal columns). Built from the stamped
  ! parts alone - no input matrix - so it is a plain builder.
  !===================================================================!

  pure function tentative_prolongation(agg, naggr) result(P0)
    integer         , intent(in) :: agg(:), naggr
    type(csr_matrix)             :: P0
    integer , allocatable :: cnt(:), row_ptr(:), col_idx(:)
    real(dp), allocatable :: vals(:)
    integer :: n, i

    n = size(agg)
    allocate(cnt(naggr)); cnt = 0
    do i = 1, n
       cnt(agg(i)) = cnt(agg(i)) + 1
    end do

    allocate(row_ptr(n+1), col_idx(n), vals(n))
    do i = 1, n
       row_ptr(i) = i
       col_idx(i) = agg(i)
       vals(i)    = 1.0_dp/sqrt(real(cnt(agg(i)), dp))
    end do
    row_ptr(n+1) = n + 1

    P0 = csr_matrix(n, naggr, row_ptr, col_idx, vals)

  end function tentative_prolongation

  !===================================================================!
  ! Self-contained dense LU with partial pivoting (coarse direct solve;
  ! the coarsest level is small, so a plain O(n^3) factor is fine and
  ! keeps the library free of any LAPACK dependency).
  !===================================================================!

  pure subroutine dense_lu_factor(A, n, ipiv, info)
    integer , intent(in)    :: n
    real(dp), intent(inout) :: A(n,n)
    integer , intent(out)   :: ipiv(n)
    integer , intent(out)   :: info
    integer  :: i, j, k, p
    real(dp) :: amax, t
    info = 0
    do k = 1, n
       p    = k
       amax = abs(A(k,k))
       do i = k+1, n
          if (abs(A(i,k)) .gt. amax) then
             amax = abs(A(i,k))
             p = i
          end if
       end do
       ipiv(k) = p
       if (amax .eq. 0.0_dp) then
          info = k
          return
       end if
       if (p .ne. k) then
          do j = 1, n
             t = A(k,j); A(k,j) = A(p,j); A(p,j) = t
          end do
       end if
       do i = k+1, n
          A(i,k) = A(i,k)/A(k,k)
          do j = k+1, n
             A(i,j) = A(i,j) - A(i,k)*A(k,j)
          end do
       end do
    end do
  end subroutine dense_lu_factor

  pure subroutine dense_lu_solve(A, n, ipiv, b)
    integer , intent(in)    :: n
    real(dp), intent(in)    :: A(n,n)
    integer , intent(in)    :: ipiv(n)
    real(dp), intent(inout) :: b(n)
    integer  :: i, k, p
    real(dp) :: t
    ! apply the recorded row swaps
    do k = 1, n
       p = ipiv(k)
       if (p .ne. k) then
          t = b(k); b(k) = b(p); b(p) = t
       end if
    end do
    ! forward solve (unit lower triangle)
    do i = 2, n
       do k = 1, i-1
          b(i) = b(i) - A(i,k)*b(k)
       end do
    end do
    ! back solve (upper triangle)
    do i = n, 1, -1
       do k = i+1, n
          b(i) = b(i) - A(i,k)*b(k)
       end do
       b(i) = b(i)/A(i,i)
    end do
  end subroutine dense_lu_solve

end module interface_multigrid
