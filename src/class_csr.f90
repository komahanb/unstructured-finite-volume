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

  use iso_fortran_env   , only : dp => REAL64
  use class_stored_graph, only : stored_graph

  implicit none

  private
  public :: csr_matrix

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
     ! sparse algebra (was the free csr_* module functions)
     procedure :: transpose
     procedure :: matmat
     procedure :: add
     procedure :: matvec_rows
     procedure :: principal_submatrix
     ! the matrix's own graph: rows are vertices, off-diagonals edges
     procedure :: adjacency_graph
  end type csr_matrix

  interface csr_matrix
     module procedure csr_from_arrays
     module procedure csr_from_pattern
  end interface csr_matrix

contains

  !===================================================================!
  ! Construct from full (row_ptr, col_idx, vals) arrays
  !===================================================================!

  pure type(csr_matrix) function csr_from_arrays(nrows, ncols, row_ptr, col_idx, vals) result(A)
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

  pure type(csr_matrix) function csr_from_pattern(nrows, ncols, row_ptr, col_idx) result(A)
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

  pure subroutine matvec(this, x, y)
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

  pure subroutine matvec_transpose(this, x, y)
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

  pure function get_diagonal(this) result(d)
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

  pure real(dp) function get_entry(this, i, j) result(aij)
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

  impure subroutine add_entry(this, i, j, v)
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

  pure subroutine scale_rows(this, s)
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

  pure real(dp) function is_symmetric(this) result(asym)
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

  pure subroutine to_dense(this, D)
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

  pure function transpose(A) result(At)
    class(csr_matrix), intent(in) :: A
    type(csr_matrix)              :: At
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
  end function transpose

  !===================================================================!
  ! C = A B  (Gustavson: symbolic count, then numeric accumulate)
  !===================================================================!

  pure function matmat(A, B) result(C)
    class(csr_matrix), intent(in) :: A
    type(csr_matrix) , intent(in) :: B
    type(csr_matrix)              :: C
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
  end function matmat

  !===================================================================!
  ! C = alpha A + beta B  (same shape; union pattern)
  !===================================================================!

  pure function add(A, alpha, beta, B) result(C)
    class(csr_matrix), intent(in) :: A
    real(dp)         , intent(in) :: alpha, beta
    type(csr_matrix) , intent(in) :: B
    type(csr_matrix)              :: C
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
  end function add

  !===================================================================!
  ! y(row) = sum_j A(row,j) x(j) for the listed rows only; other entries
  ! of y are left untouched. x must already carry valid values on every
  ! column those rows touch (e.g. owned + ghost after a halo exchange).
  !===================================================================!

  pure subroutine matvec_rows(this, x, y, rows)

    class(csr_matrix), intent(in)    :: this
    real(dp)         , intent(in)    :: x(:)
    real(dp)         , intent(inout) :: y(:)
    integer          , intent(in)    :: rows(:)

    integer  :: i, row, k
    real(dp) :: s

    do i = 1, size(rows)
       row = rows(i)
       s = 0.0_dp
       do k = this % row_ptr(row), this % row_ptr(row+1) - 1
          s = s + this % vals(k) * x(this % col_idx(k))
       end do
       y(row) = s
    end do

  end subroutine matvec_rows

  !===================================================================!
  ! Principal submatrix on the index set idx, renumbered to 1..size(idx):
  ! keep only rows and columns in idx (rows == cols), dropping the rest.
  ! Used to extract a subdomain's owned-owned block from the global
  ! operator (each image then builds a block preconditioner on it).
  !===================================================================!

  pure function principal_submatrix(this, idx) result(B)

    class(csr_matrix), intent(in) :: this
    integer          , intent(in) :: idx(:)
    type(csr_matrix)              :: B

    integer , allocatable :: loc(:), row_ptr(:), col_idx(:)
    real(dp), allocatable :: vals(:)
    integer :: n, il, row, k, jl, pos, nnz

    n = size(idx)

    ! global -> local map (0 = not in the set)
    allocate(loc(this % nrows)); loc = 0
    do il = 1, n
       loc(idx(il)) = il
    end do

    ! symbolic: count kept entries per local row
    allocate(row_ptr(n+1)); row_ptr(1) = 1
    do il = 1, n
       row = idx(il)
       pos = 0
       do k = this % row_ptr(row), this % row_ptr(row+1) - 1
          if (loc(this % col_idx(k)) .gt. 0) pos = pos + 1
       end do
       row_ptr(il+1) = row_ptr(il) + pos
    end do
    nnz = row_ptr(n+1) - 1
    allocate(col_idx(nnz), vals(nnz))

    ! numeric: copy entries, remapping global columns to local indices
    pos = 1
    do il = 1, n
       row = idx(il)
       do k = this % row_ptr(row), this % row_ptr(row+1) - 1
          jl = loc(this % col_idx(k))
          if (jl .gt. 0) then
             col_idx(pos) = jl
             vals(pos)    = this % vals(k)
             pos = pos + 1
          end if
       end do
    end do

    B = csr_matrix(n, n, row_ptr, col_idx, vals)

  end function principal_submatrix

  !===================================================================!
  ! The matrix's own graph: one vertex per row, one edge per coupled
  ! pair of rows - whichever triangle carries the coupling. An entry
  ! above the diagonal always records its edge; one below records it
  ! only when its mirror is structurally absent, so an unsymmetric
  ! pattern (an upwinded operator, say) loses nothing and a symmetric
  ! one records nothing twice. The sparsity pattern is graph-shaped
  ! data - whoever needs to traverse, partition, or refine the
  ! matrix's structure asks for this instead of walking row_ptr by
  ! hand.
  !===================================================================!

  pure type(stored_graph) function adjacency_graph(this) result(g)

    class(csr_matrix), intent(in) :: this

    integer, allocatable :: tails(:), heads(:)
    integer              :: i, k, j, ne, pass

    do pass = 1, 2
       ne = 0
       do i = 1, this % nrows
          do k = this % row_ptr(i), this % row_ptr(i+1) - 1
             j = this % col_idx(k)
             if (j .gt. i .or. (j .lt. i .and. .not. has_entry(this, j, i))) then
                ne = ne + 1
                if (pass .eq. 2) then
                   tails(ne) = min(i, j)
                   heads(ne) = max(i, j)
                end if
             end if
          end do
       end do
       if (pass .eq. 1) allocate(tails(ne), heads(ne))
    end do

    g = stored_graph(this % nrows, tails, heads)

  end function adjacency_graph

  !===================================================================!
  ! Whether the pattern holds an entry at (row, col) - structure
  ! only, the value may be anything
  !===================================================================!

  pure logical function has_entry(this, row, col)

    class(csr_matrix), intent(in) :: this
    integer          , intent(in) :: row, col

    integer :: k

    has_entry = .false.
    do k = this % row_ptr(row), this % row_ptr(row+1) - 1
       if (this % col_idx(k) .eq. col) then
          has_entry = .true.
          return
       end if
    end do

  end function has_entry

end module class_csr
