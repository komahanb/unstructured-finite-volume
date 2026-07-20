!=====================================================================!
! Sparse matrix in compressed-sparse-row form - which IS a weighted
! directed graph, stored. The compressed rows are the out-edge lists;
! this class extends the digraph and keeps NO structure of its own,
! only the weights. One object, two vocabularies, one storage:
!
!    matrix speak          graph speak
!    --------------------  --------------------------------
!    row i                 vertex v            num_vertices
!    column index j        edge v --> j        out_adj
!    row pointer           v's out-edge slice  out_xadj
!    stored entry a_vj     the edge's weight   vals
!    nnz                   num_edges (a diagonal entry is a
!                          self-loop: an edge home to itself)
!
!    [ a b . ]         a          b
!    [ . c . ]   ==   (1)<-.  (1)---->(2)<-.        the matrix,
!    [ d . e ]         ^   |          |    | c      drawn as the
!                      |   '----------'----'        graph it is
!                      | d         e
!                     (3)<--------------.
!                      '----------------'
!
! matvec: every vertex dots its out-edge weights with the values at
! the far ends - the per-vertex inner product the whole solver world
! is built from. A rectangular matrix (prolongation P, fine x coarse)
! is the bipartite case: the vertices are the rows, and their
! out-edges point at column labels 1..ncols on the far side.
!
! Everything the graph interface owns is therefore available ON the
! matrix - it does not carry a graph, it IS one.
!
! The handful of sparse operations the smoothed-aggregation AMG needs
! ride on top: products, transpose, Gustavson matmat, sparse sum, and
! the Galerkin triple product's ingredients.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_csr

  use iso_fortran_env   , only : dp => REAL64
  use interface_graph   , only : digraph, counting_sort
  use class_stored_graph, only : stored_graph

  implicit none

  private
  public :: csr_matrix

  !-------------------------------------------------------------------!
  ! The weighted digraph: structure inherited (num_vertices, out_xadj,
  ! out_adj, num_edges), weights and the bipartite far side here
  !-------------------------------------------------------------------!

  type, extends(digraph) :: csr_matrix
     integer               :: ncols = 0    ! the far side's label count
     real(dp), allocatable :: vals(:)      ! one weight per out-edge
   contains
     ! the directed contract, answered by the stored out-lists (the
     ! in-lists arrive the day a consumer needs them - via transpose)
     procedure :: out_neighbours
     procedure :: in_neighbours
     ! weighted actions along the edges
     procedure :: matvec
     procedure :: matvec_transpose
     procedure :: get_diagonal
     procedure :: get_entry
     procedure :: add_entry
     procedure :: scale_rows
     procedure :: is_symmetric
     procedure :: to_dense
     ! sparse algebra
     procedure :: transpose
     procedure :: matmat
     procedure :: add
     procedure :: matvec_rows
     procedure :: principal_submatrix
     procedure :: local_block
     ! the underlying simple graph: self-loops dropped, direction
     ! forgotten, mirrored pairs recorded once
     procedure :: simple_graph
  end type csr_matrix

  interface csr_matrix
     module procedure csr_from_arrays
     module procedure csr_from_pattern
  end interface csr_matrix

contains

  !===================================================================!
  ! Construct from full (row_ptr, col_idx, vals) arrays - the caller
  ! speaks matrix, the object stores graph
  !===================================================================!

  pure type(csr_matrix) function csr_from_arrays(nrows, ncols, row_ptr, col_idx, vals) result(A)
    integer , intent(in) :: nrows, ncols
    integer , intent(in) :: row_ptr(:), col_idx(:)
    real(dp), intent(in) :: vals(:)
    A % num_vertices = nrows
    A % ncols        = ncols
    A % num_edges    = size(col_idx)
    A % out_xadj     = row_ptr
    A % out_adj      = col_idx
    A % vals         = vals
  end function csr_from_arrays

  !===================================================================!
  ! Construct from a known sparsity pattern (row_ptr, col_idx) with the
  ! weights zeroed - then accumulate with add_entry (used by assembly)
  !===================================================================!

  pure type(csr_matrix) function csr_from_pattern(nrows, ncols, row_ptr, col_idx) result(A)
    integer, intent(in) :: nrows, ncols
    integer, intent(in) :: row_ptr(:), col_idx(:)
    A % num_vertices = nrows
    A % ncols        = ncols
    A % num_edges    = size(col_idx)
    A % out_xadj     = row_ptr
    A % out_adj      = col_idx
    allocate(A % vals(A % num_edges))
    A % vals = 0.0_dp
  end function csr_from_pattern

  !===================================================================!
  ! The directed contract, read off the stored out-lists
  !===================================================================!

  pure function out_neighbours(this, v) result(nbrs)
    class(csr_matrix), intent(in) :: this
    integer          , intent(in) :: v
    integer, allocatable :: nbrs(:)
    nbrs = this % stored_out_neighbours(v)
  end function out_neighbours

  pure function in_neighbours(this, v) result(nbrs)
    class(csr_matrix), intent(in) :: this
    integer          , intent(in) :: v
    integer, allocatable :: nbrs(:)
    nbrs = this % stored_in_neighbours(v)
  end function in_neighbours

  !===================================================================!
  ! y = A x: every vertex dots its out-edges -
  !
  !            w1
  !    (v) ---------> (j1)      y(v) = w1*x(j1) + w2*x(j2) + ...
  !      \ ---------> (j2)      one inner product per vertex,
  !            w2               over its out-edge weights
  !===================================================================!

  pure subroutine matvec(this, x, y)
    class(csr_matrix), intent(in)  :: this
    real(dp)         , intent(in)  :: x(:)
    real(dp)         , intent(out) :: y(:)
    integer :: v, e
    do v = 1, this % num_vertices
       y(v) = 0.0_dp
       do e = this % out_xadj(v), this % out_xadj(v+1) - 1
          y(v) = y(v) + this % vals(e)*x(this % out_adj(e))
       end do
    end do
  end subroutine matvec

  !===================================================================!
  ! y = A^T x: the same edges walked against their arrows - every
  ! vertex PUSHES its value out along its edges instead of pulling
  !===================================================================!

  pure subroutine matvec_transpose(this, x, y)
    class(csr_matrix), intent(in)  :: this
    real(dp)         , intent(in)  :: x(:)
    real(dp)         , intent(out) :: y(:)
    integer :: v, e
    y(1:this % ncols) = 0.0_dp
    do v = 1, this % num_vertices
       do e = this % out_xadj(v), this % out_xadj(v+1) - 1
          y(this % out_adj(e)) = y(this % out_adj(e)) + this % vals(e)*x(v)
       end do
    end do
  end subroutine matvec_transpose

  !===================================================================!
  ! The self-loop weights (square matrices): the diagonal
  !===================================================================!

  pure function get_diagonal(this) result(d)
    class(csr_matrix), intent(in) :: this
    real(dp) :: d(this % num_vertices)
    integer  :: v, e
    d = 0.0_dp
    do v = 1, this % num_vertices
       do e = this % out_xadj(v), this % out_xadj(v+1) - 1
          if (this % out_adj(e) .eq. v) d(v) = this % vals(e)
       end do
    end do
  end function get_diagonal

  !===================================================================!
  ! a_ij: the weight of edge i --> j (0 if the edge is absent)
  !===================================================================!

  pure real(dp) function get_entry(this, i, j) result(aij)
    class(csr_matrix), intent(in) :: this
    integer          , intent(in) :: i, j
    integer :: e
    aij = 0.0_dp
    do e = this % out_xadj(i), this % out_xadj(i+1) - 1
       if (this % out_adj(e) .eq. j) then
          aij = this % vals(e)
          return
       end if
    end do
  end function get_entry

  !===================================================================!
  ! Accumulate v into the weight of edge (i,j); the edge must exist
  !===================================================================!

  impure subroutine add_entry(this, i, j, v)
    class(csr_matrix), intent(inout) :: this
    integer          , intent(in)    :: i, j
    real(dp)         , intent(in)    :: v
    integer :: e
    do e = this % out_xadj(i), this % out_xadj(i+1) - 1
       if (this % out_adj(e) .eq. j) then
          this % vals(e) = this % vals(e) + v
          return
       end if
    end do
    error stop "csr add_entry: (i,j) not in sparsity pattern"
  end subroutine add_entry

  !===================================================================!
  ! Scale every out-edge of vertex i by s(i), in place
  !===================================================================!

  pure subroutine scale_rows(this, s)
    class(csr_matrix), intent(inout) :: this
    real(dp)         , intent(in)    :: s(:)
    integer :: v, e
    do v = 1, this % num_vertices
       do e = this % out_xadj(v), this % out_xadj(v+1) - 1
          this % vals(e) = this % vals(e)*s(v)
       end do
    end do
  end subroutine scale_rows

  !===================================================================!
  ! max |a_ij - a_ji|: how far the weights are from riding both ways
  ! equally (square matrices); diagnostic
  !===================================================================!

  pure real(dp) function is_symmetric(this) result(asym)
    class(csr_matrix), intent(in) :: this
    integer :: v, e, j
    asym = 0.0_dp
    do v = 1, this % num_vertices
       do e = this % out_xadj(v), this % out_xadj(v+1) - 1
          j = this % out_adj(e)
          asym = max(asym, abs(this % vals(e) - this % get_entry(j, v)))
       end do
    end do
  end function is_symmetric

  !===================================================================!
  ! Dense copy (small coarse operators only)
  !===================================================================!

  pure subroutine to_dense(this, D)
    class(csr_matrix)    , intent(in)  :: this
    real(dp), allocatable, intent(out) :: D(:,:)
    integer :: v, e
    allocate(D(this % num_vertices, this % ncols))
    D = 0.0_dp
    do v = 1, this % num_vertices
       do e = this % out_xadj(v), this % out_xadj(v+1) - 1
          D(v, this % out_adj(e)) = this % vals(e)
       end do
    end do
  end subroutine to_dense

  !===================================================================!
  ! A^T as a fresh matrix: every arrow reversed. Grouping the edges
  ! by their far end is the graph's counting kernel; it returns the
  ! new out-lists and the edge permutation, and the transpose is two
  ! gathers through it.
  !===================================================================!

  pure function transpose(A) result(At)
    class(csr_matrix), intent(in) :: A
    type(csr_matrix)              :: At
    integer, allocatable :: rows(:), perm(:)
    integer :: v, e

    ! the tail of each stored edge, in storage order
    allocate(rows(A % num_edges))
    do v = 1, A % num_vertices
       rows(A % out_xadj(v) : A % out_xadj(v+1) - 1) = v
    end do

    call counting_sort(A % ncols, A % out_adj, [(e, e = 1, A % num_edges)], &
         & At % out_xadj, perm)

    At % num_vertices = A % ncols
    At % ncols        = A % num_vertices
    At % num_edges    = A % num_edges
    At % out_adj      = rows(perm)
    At % vals         = A % vals(perm)
  end function transpose

  !===================================================================!
  ! C = A B (Gustavson): C's edges are the two-step paths - vertex v
  ! reaches j through any middle vertex - with path weights
  ! multiplied and parallel paths summed. Symbolic count, then
  ! numeric accumulate through a marker.
  !===================================================================!

  pure function matmat(A, B) result(C)
    class(csr_matrix), intent(in) :: A
    type(csr_matrix) , intent(in) :: B
    type(csr_matrix)              :: C
    integer, allocatable :: marker(:)
    integer :: v, ea, eb, j, col, ne, rstart, length
    real(dp) :: w

    C % num_vertices = A % num_vertices
    C % ncols        = B % ncols
    allocate(C % out_xadj(C % num_vertices + 1))
    allocate(marker(B % ncols)); marker = 0

    ! symbolic: distinct two-step destinations per vertex
    ne = 0
    do v = 1, A % num_vertices
       C % out_xadj(v) = ne + 1
       do ea = A % out_xadj(v), A % out_xadj(v+1) - 1
          j = A % out_adj(ea)
          do eb = B % out_xadj(j), B % out_xadj(j+1) - 1
             col = B % out_adj(eb)
             if (marker(col) .ne. v) then
                marker(col) = v
                ne = ne + 1
             end if
          end do
       end do
    end do
    C % out_xadj(C % num_vertices + 1) = ne + 1
    C % num_edges = ne
    allocate(C % out_adj(ne), C % vals(ne))

    ! numeric: marker(col) holds the absolute fill slot for the vertex
    marker = 0
    do v = 1, A % num_vertices
       rstart = C % out_xadj(v)
       length = 0
       do ea = A % out_xadj(v), A % out_xadj(v+1) - 1
          j = A % out_adj(ea)
          w = A % vals(ea)
          do eb = B % out_xadj(j), B % out_xadj(j+1) - 1
             col = B % out_adj(eb)
             if (marker(col) .lt. rstart) then
                marker(col)              = rstart + length
                C % out_adj(marker(col)) = col
                C % vals(marker(col))    = w*B % vals(eb)
                length = length + 1
             else
                C % vals(marker(col)) = C % vals(marker(col)) + w*B % vals(eb)
             end if
          end do
       end do
    end do
  end function matmat

  !===================================================================!
  ! C = alpha A + beta B (same shape): the union of the two edge
  ! sets, weights summed where the edges coincide
  !===================================================================!

  pure function add(A, alpha, beta, B) result(C)
    class(csr_matrix), intent(in) :: A
    real(dp)         , intent(in) :: alpha, beta
    type(csr_matrix) , intent(in) :: B
    type(csr_matrix)              :: C
    integer, allocatable :: marker(:)
    integer :: v, e, col, ne, rstart, length

    C % num_vertices = A % num_vertices
    C % ncols        = A % ncols
    allocate(C % out_xadj(C % num_vertices + 1))
    allocate(marker(C % ncols)); marker = 0

    ! symbolic
    ne = 0
    do v = 1, C % num_vertices
       C % out_xadj(v) = ne + 1
       do e = A % out_xadj(v), A % out_xadj(v+1) - 1
          col = A % out_adj(e)
          if (marker(col) .ne. v) then; marker(col) = v; ne = ne + 1; end if
       end do
       do e = B % out_xadj(v), B % out_xadj(v+1) - 1
          col = B % out_adj(e)
          if (marker(col) .ne. v) then; marker(col) = v; ne = ne + 1; end if
       end do
    end do
    C % out_xadj(C % num_vertices + 1) = ne + 1
    C % num_edges = ne
    allocate(C % out_adj(ne), C % vals(ne))

    ! numeric
    marker = 0
    do v = 1, C % num_vertices
       rstart = C % out_xadj(v)
       length = 0
       do e = A % out_xadj(v), A % out_xadj(v+1) - 1
          col = A % out_adj(e)
          if (marker(col) .lt. rstart) then
             marker(col)              = rstart + length
             C % out_adj(marker(col)) = col
             C % vals(marker(col))    = alpha*A % vals(e)
             length = length + 1
          else
             C % vals(marker(col)) = C % vals(marker(col)) + alpha*A % vals(e)
          end if
       end do
       do e = B % out_xadj(v), B % out_xadj(v+1) - 1
          col = B % out_adj(e)
          if (marker(col) .lt. rstart) then
             marker(col)              = rstart + length
             C % out_adj(marker(col)) = col
             C % vals(marker(col))    = beta*B % vals(e)
             length = length + 1
          else
             C % vals(marker(col)) = C % vals(marker(col)) + beta*B % vals(e)
          end if
       end do
    end do
  end function add

  !===================================================================!
  ! y(row) = the listed vertices' dots only; other entries of y are
  ! left untouched. x must already carry valid values on every far
  ! end those vertices reach (e.g. owned + ghost after a halo
  ! exchange).
  !===================================================================!

  pure subroutine matvec_rows(this, x, y, rows)

    class(csr_matrix), intent(in)    :: this
    real(dp)         , intent(in)    :: x(:)
    real(dp)         , intent(inout) :: y(:)
    integer          , intent(in)    :: rows(:)

    integer  :: i, v, e
    real(dp) :: s

    do i = 1, size(rows)
       v = rows(i)
       s = 0.0_dp
       do e = this % out_xadj(v), this % out_xadj(v+1) - 1
          s = s + this % vals(e) * x(this % out_adj(e))
       end do
       y(v) = s
    end do

  end subroutine matvec_rows

  !===================================================================!
  ! The induced subgraph on the index set idx, renumbered 1..size(idx):
  ! keep only the vertices in idx and the edges between them, dropping
  ! the rest. Used to extract a subdomain's owned-owned block from the
  ! global operator (each image then builds a block preconditioner on
  ! it).
  !===================================================================!

  pure function principal_submatrix(this, idx) result(B)

    class(csr_matrix), intent(in) :: this
    integer          , intent(in) :: idx(:)
    type(csr_matrix)              :: B

    integer , allocatable :: loc(:), row_ptr(:), col_idx(:)
    real(dp), allocatable :: vals(:)
    integer :: n, il, v, e, jl, pos, ne

    n = size(idx)

    ! global -> local map (0 = not in the set)
    allocate(loc(this % num_vertices)); loc = 0
    do il = 1, n
       loc(idx(il)) = il
    end do

    ! symbolic: count kept edges per local vertex
    allocate(row_ptr(n+1)); row_ptr(1) = 1
    do il = 1, n
       v = idx(il)
       pos = 0
       do e = this % out_xadj(v), this % out_xadj(v+1) - 1
          if (loc(this % out_adj(e)) .gt. 0) pos = pos + 1
       end do
       row_ptr(il+1) = row_ptr(il) + pos
    end do
    ne = row_ptr(n+1) - 1
    allocate(col_idx(ne), vals(ne))

    ! numeric: copy edges, remapping far ends to local labels
    pos = 1
    do il = 1, n
       v = idx(il)
       do e = this % out_xadj(v), this % out_xadj(v+1) - 1
          jl = loc(this % out_adj(e))
          if (jl .gt. 0) then
             col_idx(pos) = jl
             vals(pos)    = this % vals(e)
             pos = pos + 1
          end if
       end do
    end do

    B = csr_matrix(n, n, row_ptr, col_idx, vals)

  end function principal_submatrix

  !===================================================================!
  ! The listed rows, read in a frame: the rows kept whole - every
  ! out-edge survives - with every far end renumbered through loc,
  ! the frame's global -> local map:
  !
  !    global:   r1 r2 ...  --edges--> anywhere the rows reach
  !    local :   1  2  ...  --edges--> loc(far end), within 1..nloc
  !
  ! The rows arrive in their local order, so row l of the block IS
  ! frame position l, and B % matvec on a frame-ordered vector
  ! computes exactly the owned rows of the global product. A far
  ! end the frame cannot see refuses loudly - the halo-reach
  ! invariant (everything a row touches is owned or ghost) is the
  ! caller's contract.
  !
  ! (principal_submatrix is the kin that DROPS foreign far ends;
  ! this one insists the frame covers them.)
  !===================================================================!

  pure function local_block(this, rows, loc, nloc) result(B)

    class(csr_matrix), intent(in) :: this
    integer          , intent(in) :: rows(:)
    integer          , intent(in) :: loc(:)
    integer          , intent(in) :: nloc
    type(csr_matrix)              :: B

    integer :: nr, il, v, e, pos

    nr = size(rows)

    allocate(B % out_xadj(nr + 1))
    B % out_xadj(1) = 1
    do il = 1, nr
       v = rows(il)
       B % out_xadj(il+1) = B % out_xadj(il) &
            &             + (this % out_xadj(v+1) - this % out_xadj(v))
    end do

    B % num_vertices = nr
    B % ncols        = nloc
    B % num_edges    = B % out_xadj(nr+1) - 1
    allocate(B % out_adj(B % num_edges), B % vals(B % num_edges))

    pos = 0
    do il = 1, nr
       v = rows(il)
       do e = this % out_xadj(v), this % out_xadj(v+1) - 1
          if (loc(this % out_adj(e)) .eq. 0) then
             error stop "csr local_block: an edge leaves the frame - " // &
                  & "the halo is not the reach"
          end if
          pos = pos + 1
          B % out_adj(pos) = loc(this % out_adj(e))
          B % vals(pos)    = this % vals(e)
       end do
    end do

  end function local_block

  !===================================================================!
  ! The underlying simple graph: self-loops dropped, direction
  ! forgotten, a mirrored pair of edges recorded once. An edge above
  ! the diagonal always records; one below records only when its
  ! mirror is structurally absent, so an unsymmetric pattern (an
  ! upwinded operator, say) loses nothing and a symmetric one records
  ! nothing twice. This is the view the coarsening and refinement
  ! machinery traverses.
  !===================================================================!

  pure type(stored_graph) function simple_graph(this) result(g)

    class(csr_matrix), intent(in) :: this

    integer, allocatable :: tails(:), heads(:)
    integer              :: v, e, j, ne, pass

    do pass = 1, 2
       ne = 0
       do v = 1, this % num_vertices
          do e = this % out_xadj(v), this % out_xadj(v+1) - 1
             j = this % out_adj(e)
             if (j .gt. v .or. (j .lt. v .and. .not. has_entry(this, j, v))) then
                ne = ne + 1
                if (pass .eq. 2) then
                   tails(ne) = min(v, j)
                   heads(ne) = max(v, j)
                end if
             end if
          end do
       end do
       if (pass .eq. 1) allocate(tails(ne), heads(ne))
    end do

    g = stored_graph(this % num_vertices, tails, heads)

  end function simple_graph

  !===================================================================!
  ! Whether the edge (row, col) exists - structure only, the weight
  ! may be anything
  !===================================================================!

  pure logical function has_entry(this, row, col)

    class(csr_matrix), intent(in) :: this
    integer          , intent(in) :: row, col

    integer :: e

    has_entry = .false.
    do e = this % out_xadj(row), this % out_xadj(row+1) - 1
       if (this % out_adj(e) .eq. col) then
          has_entry = .true.
          return
       end if
    end do

  end function has_entry

end module class_csr
