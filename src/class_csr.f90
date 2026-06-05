!=====================================================================!
! Sparse matrix in compressed-sparse-row (CSR) form, plus the handful of
! sparse operations the smoothed-aggregation AMG needs: matrix-vector
! products, transpose, sparse matrix-matrix product (Gustavson), sparse
! sum, the Galerkin triple product R*A*P, and the Jacobi-smoothed
! prolongation (I - omega*Dinv*A)*P0.
!
! Everything operates on a CSR matrix + plain real(dp) vectors - no mesh,
! no globals - so the AMG hierarchy built on it is reusable per-subdomain
! when the solver goes parallel.
!
! Rectangular matrices are supported (nrows /= ncols) so prolongation P
! (fine x coarse) and restriction R = P^T fit the same type.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_csr

  use iso_fortran_env, only : dp => REAL64

  implicit none

  private
  public :: csr_matrix
  public :: csr_transpose, csr_matmat, csr_add, csr_galerkin
  public :: csr_jacobi_smoothed_prolongation

  !-------------------------------------------------------------------!
  ! CSR matrix (1-based row_ptr/col_idx)
  !-------------------------------------------------------------------!

  type :: csr_matrix
     integer               :: nrows = 0
     integer               :: ncols = 0
     integer               :: nnz   = 0
     integer , allocatable :: row_ptr(:)   ! (nrows+1)
     integer , allocatable :: col_idx(:)   ! (nnz)
     real(dp), allocatable :: vals(:)      ! (nnz)
   contains
     procedure :: matvec
     procedure :: matvec_transpose
     procedure :: get_diagonal
     procedure :: get_entry
     procedure :: add_entry
     procedure :: scale_rows
     procedure :: is_symmetric
     procedure :: to_dense
  end type csr_matrix

  interface csr_matrix
     module procedure csr_from_arrays
     module procedure csr_from_pattern
  end interface csr_matrix

contains

  !===================================================================!
  ! Construct from full (row_ptr, col_idx, vals) arrays
  !===================================================================!

  type(csr_matrix) function csr_from_arrays(nrows, ncols, row_ptr, col_idx, vals) result(A)
    integer , intent(in) :: nrows, ncols
    integer , intent(in) :: row_ptr(:), col_idx(:)
    real(dp), intent(in) :: vals(:)
    A % nrows = nrows
    A % ncols = ncols
    A % nnz   = size(col_idx)
    allocate(A % row_ptr(size(row_ptr)))
    allocate(A % col_idx(size(col_idx)))
    allocate(A % vals(size(vals)))
    A % row_ptr = row_ptr
    A % col_idx = col_idx
    A % vals    = vals
  end function csr_from_arrays

  !===================================================================!
  ! Construct from a known sparsity pattern (row_ptr, col_idx) with the
  ! values zeroed - then accumulate with add_entry (used by assembly)
  !===================================================================!

  type(csr_matrix) function csr_from_pattern(nrows, ncols, row_ptr, col_idx) result(A)
    integer, intent(in) :: nrows, ncols
    integer, intent(in) :: row_ptr(:), col_idx(:)
    A % nrows = nrows
    A % ncols = ncols
    A % nnz   = size(col_idx)
    allocate(A % row_ptr(size(row_ptr)))
    allocate(A % col_idx(size(col_idx)))
    A % row_ptr = row_ptr
    A % col_idx = col_idx
    allocate(A % vals(A % nnz))
    A % vals = 0.0_dp
  end function csr_from_pattern

  !===================================================================!
  ! y = A x
  !===================================================================!

  subroutine matvec(this, x, y)
    class(csr_matrix), intent(in)  :: this
    real(dp)         , intent(in)  :: x(:)
    real(dp)         , intent(out) :: y(:)
    integer :: i, k
    do i = 1, this % nrows
       y(i) = 0.0_dp
       do k = this % row_ptr(i), this % row_ptr(i+1) - 1
          y(i) = y(i) + this % vals(k)*x(this % col_idx(k))
       end do
    end do
  end subroutine matvec

  !===================================================================!
  ! y = A^T x
  !===================================================================!

  subroutine matvec_transpose(this, x, y)
    class(csr_matrix), intent(in)  :: this
    real(dp)         , intent(in)  :: x(:)
    real(dp)         , intent(out) :: y(:)
    integer :: i, k
    y(1:this % ncols) = 0.0_dp
    do i = 1, this % nrows
       do k = this % row_ptr(i), this % row_ptr(i+1) - 1
          y(this % col_idx(k)) = y(this % col_idx(k)) + this % vals(k)*x(i)
       end do
    end do
  end subroutine matvec_transpose

  !===================================================================!
  ! Extract the diagonal (square matrices)
  !===================================================================!

  function get_diagonal(this) result(d)
    class(csr_matrix), intent(in) :: this
    real(dp) :: d(this % nrows)
    integer  :: i, k
    d = 0.0_dp
    do i = 1, this % nrows
       do k = this % row_ptr(i), this % row_ptr(i+1) - 1
          if (this % col_idx(k) .eq. i) d(i) = this % vals(k)
       end do
    end do
  end function get_diagonal

  !===================================================================!
  ! a_ij (linear scan within the row; 0 if absent)
  !===================================================================!

  real(dp) function get_entry(this, i, j) result(aij)
    class(csr_matrix), intent(in) :: this
    integer          , intent(in) :: i, j
    integer :: k
    aij = 0.0_dp
    do k = this % row_ptr(i), this % row_ptr(i+1) - 1
       if (this % col_idx(k) .eq. j) then
          aij = this % vals(k)
          return
       end if
    end do
  end function get_entry

  !===================================================================!
  ! Accumulate v into entry (i,j); (i,j) must be in the pattern
  !===================================================================!

  subroutine add_entry(this, i, j, v)
    class(csr_matrix), intent(inout) :: this
    integer          , intent(in)    :: i, j
    real(dp)         , intent(in)    :: v
    integer :: k
    do k = this % row_ptr(i), this % row_ptr(i+1) - 1
       if (this % col_idx(k) .eq. j) then
          this % vals(k) = this % vals(k) + v
          return
       end if
    end do
    error stop "csr add_entry: (i,j) not in sparsity pattern"
  end subroutine add_entry

  !===================================================================!
  ! Scale each row i by s(i) in place
  !===================================================================!

  subroutine scale_rows(this, s)
    class(csr_matrix), intent(inout) :: this
    real(dp)         , intent(in)    :: s(:)
    integer :: i, k
    do i = 1, this % nrows
       do k = this % row_ptr(i), this % row_ptr(i+1) - 1
          this % vals(k) = this % vals(k)*s(i)
       end do
    end do
  end subroutine scale_rows

  !===================================================================!
  ! max |a_ij - a_ji| (square matrices); diagnostic
  !===================================================================!

  real(dp) function is_symmetric(this) result(asym)
    class(csr_matrix), intent(in) :: this
    integer :: i, k, j
    asym = 0.0_dp
    do i = 1, this % nrows
       do k = this % row_ptr(i), this % row_ptr(i+1) - 1
          j = this % col_idx(k)
          asym = max(asym, abs(this % vals(k) - this % get_entry(j, i)))
       end do
    end do
  end function is_symmetric

  !===================================================================!
  ! Dense copy (small coarse operators only)
  !===================================================================!

  subroutine to_dense(this, D)
    class(csr_matrix)    , intent(in)  :: this
    real(dp), allocatable, intent(out) :: D(:,:)
    integer :: i, k
    allocate(D(this % nrows, this % ncols))
    D = 0.0_dp
    do i = 1, this % nrows
       do k = this % row_ptr(i), this % row_ptr(i+1) - 1
          D(i, this % col_idx(k)) = this % vals(k)
       end do
    end do
  end subroutine to_dense

  !===================================================================!
  ! A^T as a fresh CSR (count columns, prefix-sum, scatter)
  !===================================================================!

  function csr_transpose(A) result(At)
    type(csr_matrix), intent(in) :: A
    type(csr_matrix)             :: At
    integer, allocatable :: next(:)
    integer :: i, k, col, dest

    At % nrows = A % ncols
    At % ncols = A % nrows
    At % nnz   = A % nnz
    allocate(At % row_ptr(At % nrows + 1)); At % row_ptr = 0
    allocate(At % col_idx(A % nnz), At % vals(A % nnz))

    ! count entries per column of A (= per row of At)
    do k = 1, A % nnz
       At % row_ptr(A % col_idx(k) + 1) = At % row_ptr(A % col_idx(k) + 1) + 1
    end do
    At % row_ptr(1) = 1
    do i = 1, At % nrows
       At % row_ptr(i+1) = At % row_ptr(i+1) + At % row_ptr(i)
    end do

    allocate(next(At % nrows)); next = At % row_ptr(1:At % nrows)
    do i = 1, A % nrows
       do k = A % row_ptr(i), A % row_ptr(i+1) - 1
          col  = A % col_idx(k)
          dest = next(col)
          At % col_idx(dest) = i
          At % vals(dest)    = A % vals(k)
          next(col)          = dest + 1
       end do
    end do
  end function csr_transpose

  !===================================================================!
  ! C = A B  (Gustavson: symbolic count, then numeric accumulate)
  !===================================================================!

  function csr_matmat(A, B) result(C)
    type(csr_matrix), intent(in) :: A, B
    type(csr_matrix)             :: C
    integer, allocatable :: marker(:)
    integer :: i, ka, kb, j, col, nnz, rstart, length
    real(dp) :: v

    C % nrows = A % nrows
    C % ncols = B % ncols
    allocate(C % row_ptr(C % nrows + 1))
    allocate(marker(B % ncols)); marker = 0

    ! symbolic: distinct columns per row
    nnz = 0
    do i = 1, A % nrows
       C % row_ptr(i) = nnz + 1
       do ka = A % row_ptr(i), A % row_ptr(i+1) - 1
          j = A % col_idx(ka)
          do kb = B % row_ptr(j), B % row_ptr(j+1) - 1
             col = B % col_idx(kb)
             if (marker(col) .ne. i) then
                marker(col) = i
                nnz = nnz + 1
             end if
          end do
       end do
    end do
    C % row_ptr(C % nrows + 1) = nnz + 1
    C % nnz = nnz
    allocate(C % col_idx(nnz), C % vals(nnz))

    ! numeric: marker(col) holds the absolute fill slot for the current row
    marker = 0
    do i = 1, A % nrows
       rstart = C % row_ptr(i)
       length = 0
       do ka = A % row_ptr(i), A % row_ptr(i+1) - 1
          j = A % col_idx(ka)
          v = A % vals(ka)
          do kb = B % row_ptr(j), B % row_ptr(j+1) - 1
             col = B % col_idx(kb)
             if (marker(col) .lt. rstart) then
                marker(col)          = rstart + length
                C % col_idx(marker(col)) = col
                C % vals(marker(col))    = v*B % vals(kb)
                length = length + 1
             else
                C % vals(marker(col)) = C % vals(marker(col)) + v*B % vals(kb)
             end if
          end do
       end do
    end do
  end function csr_matmat

  !===================================================================!
  ! C = alpha A + beta B  (same shape; union pattern)
  !===================================================================!

  function csr_add(alpha, A, beta, B) result(C)
    real(dp)        , intent(in) :: alpha, beta
    type(csr_matrix), intent(in) :: A, B
    type(csr_matrix)             :: C
    integer, allocatable :: marker(:)
    integer :: i, k, col, nnz, rstart, length

    C % nrows = A % nrows
    C % ncols = A % ncols
    allocate(C % row_ptr(C % nrows + 1))
    allocate(marker(C % ncols)); marker = 0

    ! symbolic
    nnz = 0
    do i = 1, C % nrows
       C % row_ptr(i) = nnz + 1
       do k = A % row_ptr(i), A % row_ptr(i+1) - 1
          col = A % col_idx(k)
          if (marker(col) .ne. i) then; marker(col) = i; nnz = nnz + 1; end if
       end do
       do k = B % row_ptr(i), B % row_ptr(i+1) - 1
          col = B % col_idx(k)
          if (marker(col) .ne. i) then; marker(col) = i; nnz = nnz + 1; end if
       end do
    end do
    C % row_ptr(C % nrows + 1) = nnz + 1
    C % nnz = nnz
    allocate(C % col_idx(nnz), C % vals(nnz))

    ! numeric
    marker = 0
    do i = 1, C % nrows
       rstart = C % row_ptr(i)
       length = 0
       do k = A % row_ptr(i), A % row_ptr(i+1) - 1
          col = A % col_idx(k)
          if (marker(col) .lt. rstart) then
             marker(col) = rstart + length
             C % col_idx(marker(col)) = col
             C % vals(marker(col))    = alpha*A % vals(k)
             length = length + 1
          else
             C % vals(marker(col)) = C % vals(marker(col)) + alpha*A % vals(k)
          end if
       end do
       do k = B % row_ptr(i), B % row_ptr(i+1) - 1
          col = B % col_idx(k)
          if (marker(col) .lt. rstart) then
             marker(col) = rstart + length
             C % col_idx(marker(col)) = col
             C % vals(marker(col))    = beta*B % vals(k)
             length = length + 1
          else
             C % vals(marker(col)) = C % vals(marker(col)) + beta*B % vals(k)
          end if
       end do
    end do
  end function csr_add

  !===================================================================!
  ! Galerkin coarse operator  Ac = R A P  (R = P^T)
  !===================================================================!

  function csr_galerkin(R, A, P) result(Ac)
    type(csr_matrix), intent(in) :: R, A, P
    type(csr_matrix)             :: Ac
    type(csr_matrix)             :: AP
    AP = csr_matmat(A, P)
    Ac = csr_matmat(R, AP)
  end function csr_galerkin

  !===================================================================!
  ! Smoothed prolongation  P = (I - omega Dinv A) P0
  !                           = P0 - omega * Dinv .* (A P0)
  !===================================================================!

  function csr_jacobi_smoothed_prolongation(A, P0, Dinv, omega) result(P)
    type(csr_matrix), intent(in) :: A, P0
    real(dp)        , intent(in) :: Dinv(:)
    real(dp)        , intent(in) :: omega
    type(csr_matrix)             :: P
    type(csr_matrix)             :: AP0
    AP0 = csr_matmat(A, P0)        ! (n x ncoarse)
    call AP0 % scale_rows(Dinv)    ! Dinv .* (A P0)
    P = csr_add(1.0_dp, P0, -omega, AP0)
  end function csr_jacobi_smoothed_prolongation

end module class_csr
