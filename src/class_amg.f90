!=====================================================================!
! Smoothed-aggregation algebraic multigrid (SA-AMG), used as a
! preconditioner (one V-cycle = M^-1 r). Built once from an assembled
! sparse CSR operator; the krylov operator A*p stays matrix-free.
!
! setup builds a hierarchy of levels by, at each level:
!   - strength-of-connection  |a_ij| >= theta sqrt(|a_ii a_jj|)
!   - greedy aggregation (Vanek 3-pass): cells -> aggregates
!   - tentative prolongation P0 (piecewise-constant over aggregates,
!     column-normalized; near-null space = the constant vector)
!   - smoothed prolongation   P = (I - omega Dinv A) P0,  omega=(4/3)/rho
!   - galerkin coarse operator Ac = P^T A P
! recursing until the coarsest level is small, then a dense LU there.
!
! apply runs one recursive V-cycle: pre-smooth (weighted jacobi, which is
! symmetric so M^-1 stays SPD for CG), restrict the residual via R=P^T,
! coarse-solve (recurse; dense solve at the bottom), prolong+correct,
! post-smooth. Everything is CSR + vector ops - no mesh, no globals - so
! the hierarchy is reusable per-subdomain when the solver goes parallel.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_amg

  use iso_fortran_env        , only : dp => REAL64
  use class_csr              , only : csr_matrix, csr_transpose, csr_galerkin, &
       &                              csr_jacobi_smoothed_prolongation
  use interface_preconditioner, only : preconditioner

  implicit none

  private
  public :: amg

  !-------------------------------------------------------------------!
  ! One multigrid level
  !-------------------------------------------------------------------!

  type :: amg_level
     type(csr_matrix)      :: A           ! operator on this level
     type(csr_matrix)      :: P           ! prolongation (fine x coarse)
     type(csr_matrix)      :: R           ! restriction  = P^T
     real(dp), allocatable :: Dinv(:)     ! 1/diag(A) (smoother)
  end type amg_level

  !-------------------------------------------------------------------!
  ! The AMG hierarchy / preconditioner
  !-------------------------------------------------------------------!

  type, extends(preconditioner) :: amg

     type(amg_level), allocatable :: levels(:)   ! 1=finest .. nlevel=coarsest
     integer                      :: nlevel = 0

     ! coarsest level: dense LU factor (LAPACK)
     real(dp), allocatable :: Ac_dense(:,:)
     integer , allocatable :: ipiv(:)
     integer               :: ncoarse = 0

     ! parameters
     real(dp) :: theta       = 0.08_dp     ! strength threshold
     real(dp) :: jacobi_w    = 0.6667_dp   ! smoother weight (~2/3)
     integer  :: max_levels  = 25
     integer  :: coarse_size = 50
     integer  :: npre        = 1
     integer  :: npost       = 1

   contains

     procedure :: setup
     procedure :: apply
     procedure :: vcycle

  end type amg

contains

  !===================================================================!
  ! Build the hierarchy from the fine operator A_fine
  !===================================================================!

  subroutine setup(this, A_fine)

    class(amg)      , intent(inout) :: this
    type(csr_matrix), intent(in)    :: A_fine

    integer               :: lev, n, naggr, info, i
    integer , allocatable :: agg(:)
    real(dp), allocatable :: dinv(:)
    real(dp)              :: omega, rho
    type(csr_matrix)      :: P0

    ! re-entrant: drop any previously-built hierarchy
    if (allocated(this % levels))   deallocate(this % levels)
    if (allocated(this % Ac_dense)) deallocate(this % Ac_dense)
    if (allocated(this % ipiv))     deallocate(this % ipiv)

    allocate(this % levels(this % max_levels))
    this % levels(1) % A = A_fine

    lev = 1
    coarsen: do

       n = this % levels(lev) % A % nrows
       if (n .le. this % coarse_size .or. lev .ge. this % max_levels) exit coarsen

       ! smoother data on this level
       dinv = 1.0_dp/this % levels(lev) % A % get_diagonal()
       this % levels(lev) % Dinv = dinv

       ! aggregate, then build tentative + smoothed prolongation
       call aggregate(this % levels(lev) % A, this % theta, agg, naggr)
       if (naggr .ge. n .or. naggr .le. 1) exit coarsen   ! no useful coarsening

       call tentative_prolongation(agg, naggr, P0)

       rho   = spectral_radius(this % levels(lev) % A, dinv)
       omega = (4.0_dp/3.0_dp)/rho

       this % levels(lev) % P = csr_jacobi_smoothed_prolongation( &
            & this % levels(lev) % A, P0, dinv, omega)
       this % levels(lev) % R = csr_transpose(this % levels(lev) % P)

       ! galerkin coarse operator
       this % levels(lev+1) % A = csr_galerkin( &
            & this % levels(lev) % R, this % levels(lev) % A, this % levels(lev) % P)

       lev = lev + 1

    end do coarsen

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
       if (info .ne. 0) error stop "amg: coarse factorization failed"
    end if

  end subroutine setup

  !===================================================================!
  ! z = M^-1 r  (one V-cycle from the finest level)
  !===================================================================!

  subroutine apply(this, r, z)
    class(amg), intent(in)  :: this
    real(dp)  , intent(in)  :: r(:)
    real(dp)  , intent(out) :: z(:)
    z = 0.0_dp
    call this % vcycle(1, r, z)
  end subroutine apply

  !===================================================================!
  ! Recursive V-cycle worker
  !===================================================================!

  recursive subroutine vcycle(this, lev, b, x)

    class(amg), intent(in)    :: this
    integer   , intent(in)    :: lev
    real(dp)  , intent(in)    :: b(:)
    real(dp)  , intent(inout) :: x(:)

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
      call smooth(L % A, L % Dinv, this % jacobi_w, b, x, this % npre)

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
      call smooth(L % A, L % Dinv, this % jacobi_w, b, x, this % npost)

    end associate

  end subroutine vcycle

  !===================================================================!
  ! Weighted-jacobi smoother (symmetric -> keeps M^-1 SPD):
  !   x <- x + w Dinv (b - A x)
  !===================================================================!

  subroutine smooth(A, Dinv, w, b, x, nsweep)
    type(csr_matrix), intent(in)    :: A
    real(dp)        , intent(in)    :: Dinv(:), w, b(:)
    real(dp)        , intent(inout) :: x(:)
    integer         , intent(in)    :: nsweep
    real(dp), allocatable :: Ax(:)
    integer :: s
    allocate(Ax(A % nrows))
    do s = 1, nsweep
       call A % matvec(x, Ax)
       x = x + w*Dinv*(b - Ax)
    end do
  end subroutine smooth

  !===================================================================!
  ! Greedy aggregation (Vanek 3-pass) on the strength graph
  !   strong(i,j):  |a_ij| >= theta sqrt(|a_ii a_jj|),  i /= j
  !===================================================================!

  subroutine aggregate(A, theta, agg, naggr)
    type(csr_matrix)    , intent(in)  :: A
    real(dp)            , intent(in)  :: theta
    integer, allocatable, intent(out) :: agg(:)
    integer             , intent(out) :: naggr
    real(dp), allocatable :: d(:)
    integer :: n, i, k, j
    logical :: isolated

    n = A % nrows
    allocate(agg(n)); agg = 0
    naggr = 0
    d = A % get_diagonal()

    ! pass 1: a node whose strong neighbourhood is entirely free seeds a
    ! new aggregate = the node + its strong neighbours
    do i = 1, n
       if (agg(i) .ne. 0) cycle
       isolated = .true.
       do k = A % row_ptr(i), A % row_ptr(i+1) - 1
          j = A % col_idx(k)
          if (j .ne. i .and. strong(A % vals(k), d(i), d(j), theta)) then
             if (agg(j) .ne. 0) then
                isolated = .false.
                exit
             end if
          end if
       end do
       if (.not. isolated) cycle
       naggr = naggr + 1
       agg(i) = naggr
       do k = A % row_ptr(i), A % row_ptr(i+1) - 1
          j = A % col_idx(k)
          if (j .ne. i .and. strong(A % vals(k), d(i), d(j), theta)) agg(j) = naggr
       end do
    end do

    ! pass 2: sweep each remaining node into a strong, already-aggregated nbr
    do i = 1, n
       if (agg(i) .ne. 0) cycle
       do k = A % row_ptr(i), A % row_ptr(i+1) - 1
          j = A % col_idx(k)
          if (j .ne. i .and. agg(j) .ne. 0 .and. strong(A % vals(k), d(i), d(j), theta)) then
             agg(i) = agg(j)
             exit
          end if
       end do
    end do

    ! pass 3: anything left seeds its own (singleton) aggregate
    do i = 1, n
       if (agg(i) .eq. 0) then
          naggr = naggr + 1
          agg(i) = naggr
       end if
    end do

  end subroutine aggregate

  !===================================================================!
  ! strength test
  !===================================================================!

  pure logical function strong(aij, aii, ajj, theta)
    real(dp), intent(in) :: aij, aii, ajj, theta
    strong = abs(aij) .ge. theta*sqrt(abs(aii*ajj))
  end function strong

  !===================================================================!
  ! Tentative prolongation P0 (n x naggr): piecewise constant over
  ! aggregates, column-normalized so the constant near-null space is
  ! represented exactly (orthonormal columns)
  !===================================================================!

  subroutine tentative_prolongation(agg, naggr, P0)
    integer         , intent(in)  :: agg(:), naggr
    type(csr_matrix), intent(out) :: P0
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

  end subroutine tentative_prolongation

  !===================================================================!
  ! Spectral radius of (Dinv A) by bounded power iteration, floored by
  ! the gershgorin bound max_i sum_j |Dinv A|_ij - used for the smoothed
  ! prolongation weight omega = (4/3)/rho
  !===================================================================!

  real(dp) function spectral_radius(A, Dinv) result(rho)
    type(csr_matrix), intent(in) :: A
    real(dp)        , intent(in) :: Dinv(:)
    real(dp), allocatable :: v(:), w(:)
    real(dp) :: nv, gers, rowsum
    integer  :: n, it, i, k

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

  end function spectral_radius

  !===================================================================!
  ! Self-contained dense LU with partial pivoting (coarse direct solve;
  ! the coarsest level is small, so a plain O(n^3) factor is fine and
  ! keeps the library free of any LAPACK dependency).
  !===================================================================!

  subroutine dense_lu_factor(A, n, ipiv, info)
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

  subroutine dense_lu_solve(A, n, ipiv, b)
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

end module class_amg
