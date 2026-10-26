!=====================================================================!
! Module that states a factorisation retained: P A = L U with partial
! pivoting, formed once from a square matrix and applied to every
! right side passed to it, against A or against its transpose.
!
!      P A  =  L U        | 1          |   | u  u  u  u |
!                         | l  1       |   |    u  u  u |
!                         | l  l  1    |   |       u  u |
!                         | l  l  l  1 |   |          u |
!
!      stored in place, one array: the l below the diagonal, the u on
!      and above it, the unit diagonal of L implied
!
! The purpose of retaining it is the operation count. A tangent solve,
! an adjoint solve and every order of an expansion read the same
! matrix frozen at the same state, so one factorisation serves them
! all and each costs a substitution. Formed again per solve, each would
! cost the factorisation again. Each factorisation is recorded by its
! caller, so the count is measured rather than modelled.
!
! Solving against the transpose uses the same factors. P A = L U gives
! A^T P^T = U^T L^T, so A^T x = b is U^T y = b, then L^T z = y, then
! x = P^T z, which is the row exchanges of the factorisation applied
! in the reverse of the order they were made.
!
! A pivot at or below the threshold given leaves the factorisation
! singular: no further elimination is done, and a substitution against
! it stops the program. Whether a singular pivot is an error or an
! expected result is the caller's decision, read from singular().
!
! Equilibrated. The matrix factorised is D_r A D_c, its rows and then
! its columns scaled to a largest entry of one, and a substitution
! scales its right side and its solution by them:
!
!      x  =  D_c (D_r A D_c)^-1 D_r b
!      x  =  D_r (D_c A^T D_r)^-1 D_c b       against the transpose
!
! A system of rows of different units, positions, rotations and
! reactions, is solved to the digits its scaled matrix allows.
!
! Block triangular. The matrix is permuted, rows and columns alike, to
! the block lower triangular form its zero pattern allows, and only
! its diagonal blocks are factorised:
!
!      | A_11             |     x_1 = A_11^-1  b_1
!      | A_21 A_22        |     x_2 = A_22^-1 (b_2 - A_21 x_1)
!      | A_31 A_32 A_33   |     x_3 = A_33^-1 (b_3 - A_31 x_1 - A_32 x_2)
!
! The blocks are the strongly connected components of the graph of the
! pattern, an edge from i to j where a_ij is not zero, found by
! Tarjan's algorithm. It emits a component after every component it
! reads, which is the order of a forward substitution:
! x_b = A_bb^-1 (b_b - sum over c before b of A_bc x_c). Against the
! transpose the blocks are solved in reverse. The matrix is nonsingular
! exactly where every diagonal block is. The step of a family of s
! stages is such a matrix, each stage reading the stages before it:
! its factorisation costs s blocks of one stage, about s^2 times less
! than the whole. A row is the law of the unknown of the same index, so
! the permutation is the same on both sides.
!
! From its entries. The block triangular factorisation of a matrix
! stated by its nonzero entries forms the diagonal blocks alone, and
! stores the entries below them by rows: no array of the order of the
! matrix squared is formed, and a substitution reads the entries below
! the blocks once each. Its scaling, its blocks and their factors are
! those of the matrix the entries state.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module util_factorisation

  use util_precision, only : dp

  implicit none

  private

  public :: dense_factorisation

  !-------------------------------------------------------------------!
  ! Type that models the factorisation retained. One of three forms
  ! is stored:
  !
  !      dense             lu, pivot_rows
  !      block triangular  permutation, block_first, diagonal(:), and
  !                        lu the permuted matrix
  !      from its entries  permutation, block_first, diagonal(:), and
  !                        the entries below the blocks by rows
  !
  ! row_scale and column_scale are stored when equilibrated.
  !-------------------------------------------------------------------!

  type :: dense_factorisation

     real(dp), allocatable, private :: lu(:,:)
     real(dp), allocatable, private :: row_scale(:), column_scale(:)
     integer , allocatable, private :: pivot_rows(:)
     integer              , private :: n = 0
     logical              , private :: is_singular = .false.

     ! Block triangular: the unknowns in the order of the blocks, the
     ! first of each block in that order, the permuted matrix, and the
     ! factorisation of each diagonal block
     integer                  , allocatable, private :: permutation(:), block_first(:)
     type(dense_factorisation), allocatable, private :: diagonal(:)

     ! From its entries: the entries below the diagonal blocks, by the
     ! rows of the permuted matrix, their columns permuted alike
     integer , allocatable, private :: lower_start(:), lower_column(:)
     real(dp), allocatable, private :: lower_value(:)

   contains

     ! The factorisation
     procedure :: factorise
     procedure :: factorise_entries

     ! The substitutions against it
     procedure :: substitute
     procedure :: substituted

     ! What it is
     procedure :: order
     procedure :: singular

  end type dense_factorisation

contains

  !===================================================================!
  ! Factorises a square matrix, L U in place in a copy. A matrix that
  ! is not square stops the program. A pivot at or below the threshold
  ! leaves the factorisation singular rather than dividing by it; given
  ! the scale of each column, at or below the threshold times its
  ! column's scale, which equilibrated is one.
  !===================================================================!

  recursive subroutine factorise(this, a, threshold, scales, equilibrated, block_triangular)

    class(dense_factorisation), intent(inout)        :: this
    real(dp)                  , intent(in)           :: a(:,:)
    real(dp)                  , intent(in)           :: threshold
    real(dp)                  , intent(in), optional :: scales(:)
    logical                   , intent(in), optional :: equilibrated, block_triangular

    real(dp), allocatable :: row(:)
    real(dp)              :: factor
    integer               :: n, k, p, j, b
    character(len=150)    :: message

    n = size(a, 1)
    if (size(a, 2) /= n) then
       write(message,'(a,i0,a,i0)') 'util_factorisation: factorise requires a square matrix; &
            &size(a,1) = ', n, ', size(a,2) = ', size(a, 2)
       error stop trim(message)
    end if

    this % n           = n
    this % lu          = a
    this % is_singular = .false.
    if (allocated(this % row_scale)) deallocate(this % row_scale, this % column_scale)
    if (present(equilibrated)) then
       if (equilibrated) then
          ! The rows' largest entries and their scaling, by sweeps of
          ! the columns, the array's memory order
          allocate(this % row_scale(n), source=0.0_dp)
          do k = 1, n
             this % row_scale = max(this % row_scale, abs(a(:, k)))
          end do
          this % row_scale = 1.0_dp / max(tiny(1.0_dp), this % row_scale)
          do k = 1, n
             this % lu(:, k) = this % lu(:, k) * this % row_scale
          end do
          this % column_scale = [(1.0_dp / max(tiny(1.0_dp), maxval(abs(this % lu(:, k)))), k = 1, n)]
          do k = 1, n
             this % lu(:, k) = this % lu(:, k) * this % column_scale(k)
          end do
       end if
    end if
    if (allocated(this % permutation)) deallocate(this % permutation, this % block_first, this % diagonal)
    if (allocated(this % lower_start)) deallocate(this % lower_start, this % lower_column, this % lower_value)
    if (present(block_triangular)) then
       if (block_triangular) then
          call strong_components(this % lu, this % permutation, this % block_first)
          this % lu = this % lu(this % permutation, this % permutation)
          allocate(this % diagonal(size(this % block_first) - 1))
          do b = 1, size(this % diagonal)
             associate (first => this % block_first(b), last => this % block_first(b + 1) - 1)
               if (allocated(this % row_scale) .or. .not. present(scales)) then
                  call this % diagonal(b) % factorise(this % lu(first:last, first:last), threshold)
               else
                  call this % diagonal(b) % factorise(this % lu(first:last, first:last), threshold, &
                       & scales(this % permutation(first:last)))
               end if
               if (this % diagonal(b) % is_singular) then
                  this % is_singular = .true.
                  return
               end if
             end associate
          end do
          return
       end if
    end if
    if (allocated(this % pivot_rows)) deallocate(this % pivot_rows)
    allocate(this % pivot_rows(n), source=0)
    allocate(row(n))

    ! Right-looking elimination, column by column: the multipliers of
    ! column k are formed at once, and each later column is updated as
    ! a whole, so every inner sweep runs down a column, which is the
    ! array's memory order.
    !
    !      column k      k+1 .. n
    !         u_kk   |   u_kj  ...        row k, after the exchange
    !         l_ik   |   a_ij - l_ik u_kj  every later column, at once
    do k = 1, n

       p = k - 1 + maxloc(abs(this % lu(k:n, k)), dim=1)
       this % pivot_rows(k) = p

       if (p /= k) then
          row             = this % lu(k, :)
          this % lu(k, :) = this % lu(p, :)
          this % lu(p, :) = row
       end if

       if (allocated(this % row_scale)) then
          if (abs(this % lu(k, k)) <= threshold) then
             this % is_singular = .true.
             return
          end if
       else if (present(scales)) then
          if (abs(this % lu(k, k)) <= threshold * scales(k)) then
             this % is_singular = .true.
             return
          end if
       else if (abs(this % lu(k, k)) <= threshold) then
          this % is_singular = .true.
          return
       end if

       this % lu(k+1:n, k) = this % lu(k+1:n, k) / this % lu(k, k)
       do j = k + 1, n
          factor = this % lu(k, j)
          if (factor == 0.0_dp) cycle
          this % lu(k+1:n, j) = this % lu(k+1:n, j) - factor * this % lu(k+1:n, k)
       end do

    end do

  end subroutine factorise

  !===================================================================!
  ! Factorises a matrix of order n stated by its entries, block
  ! triangular: weights(e) at row rows(e) and column columns(e), the
  ! weights of an entry stated more than once summed in the order
  ! given. Its scaling, its blocks and the factors of each are those of
  ! factorise on the matrix the entries state, with block_triangular
  ! and, without equilibrated, the largest entry of each column the
  ! scale of its pivots. An entry outside 1..n stops the program.
  !
  !      the entries by rows, compressed
  !
  !      start(i) .. start(i + 1) - 1     the entries of row i
  !      column(e), value(e)              their columns, ascending, and
  !                                       values, duplicates summed and
  !                                       zeros dropped
  !===================================================================!

  subroutine factorise_entries(this, n, rows, columns, weights, threshold, equilibrated)

    class(dense_factorisation), intent(inout) :: this
    integer                   , intent(in)    :: n, rows(:), columns(:)
    real(dp)                  , intent(in)    :: weights(:), threshold
    logical                   , intent(in)    :: equilibrated

    integer , allocatable :: start(:), column(:), place(:), edge_start(:), edge(:), position(:)
    real(dp), allocatable :: value(:), scales(:), block(:,:)
    integer               :: e, i, j, k, b, m, filled
    real(dp)              :: summed

    if (any(rows < 1) .or. any(rows > n) .or. any(columns < 1) .or. any(columns > n)) then
       error stop 'util_factorisation: an entry of the matrix lies outside its order'
    end if
    this % n           = n
    this % is_singular = .false.
    if (allocated(this % lu)) deallocate(this % lu)
    if (allocated(this % pivot_rows)) deallocate(this % pivot_rows)
    if (allocated(this % row_scale)) deallocate(this % row_scale, this % column_scale)
    if (allocated(this % permutation)) deallocate(this % permutation, this % block_first, this % diagonal)
    if (allocated(this % lower_start)) deallocate(this % lower_start, this % lower_column, this % lower_value)

    ! The entries by rows, each row in the order given, then by columns
    ! within it, an entry stated more than once summed from zero in that
    ! order and an entry of zero dropped
    allocate(start(n + 1), source=0)
    do e = 1, size(rows)
       start(rows(e) + 1) = start(rows(e) + 1) + 1
    end do
    start(1) = 1
    do i = 1, n
       start(i + 1) = start(i) + start(i + 1)
    end do
    place = start(1:n)
    allocate(column(size(rows)), value(size(rows)))
    do e = 1, size(rows)
       column(place(rows(e))) = columns(e)
       value(place(rows(e)))  = weights(e)
       place(rows(e)) = place(rows(e)) + 1
    end do
    filled = 0
    do i = 1, n
       call ordered(column(start(i):start(i + 1) - 1), value(start(i):start(i + 1) - 1))
       e = start(i)
       start(i) = filled + 1
       do while (e < start(i + 1))
          j      = column(e)
          summed = 0.0_dp
          do while (e < start(i + 1))
             if (column(e) /= j) exit
             summed = summed + value(e)
             e = e + 1
          end do
          if (summed == 0.0_dp) cycle
          filled = filled + 1
          column(filled) = j
          value(filled)  = summed
       end do
    end do
    start(n + 1) = filled + 1

    ! The rows and then the columns scaled to a largest entry of one; or
    ! the largest entry of each column, the scale of its pivots
    if (equilibrated) then
       allocate(this % row_scale(n), this % column_scale(n), source=0.0_dp)
       do i = 1, n
          do e = start(i), start(i + 1) - 1
             this % row_scale(i) = max(this % row_scale(i), abs(value(e)))
          end do
       end do
       this % row_scale = 1.0_dp / max(tiny(1.0_dp), this % row_scale)
       do i = 1, n
          value(start(i):start(i + 1) - 1) = value(start(i):start(i + 1) - 1) * this % row_scale(i)
       end do
       do e = 1, filled
          this % column_scale(column(e)) = max(this % column_scale(column(e)), abs(value(e)))
       end do
       this % column_scale = 1.0_dp / max(tiny(1.0_dp), this % column_scale)
       do e = 1, filled
          value(e) = value(e) * this % column_scale(column(e))
       end do
    else
       allocate(scales(n), source=0.0_dp)
       do e = 1, filled
          scales(column(e)) = max(scales(column(e)), abs(value(e)))
       end do
    end if

    ! The blocks: the components of the pattern, an edge from i to j
    ! where the scaled entry is not zero and i /= j
    allocate(edge_start(n + 1))
    allocate(edge(filled))
    k = 0
    do i = 1, n
       edge_start(i) = k + 1
       do e = start(i), start(i + 1) - 1
          if (column(e) == i .or. value(e) == 0.0_dp) cycle
          k = k + 1
          edge(k) = column(e)
       end do
    end do
    edge_start(n + 1) = k + 1
    call components_of_edges(edge_start, edge(1:k), this % permutation, this % block_first)
    allocate(position(n))
    position(this % permutation) = [(i, i = 1, n)]

    ! Each diagonal block formed and factorised; the entries below it
    ! stored by its rows
    allocate(this % lower_start(n + 1), this % lower_column(filled), this % lower_value(filled))
    allocate(this % diagonal(size(this % block_first) - 1))
    k = 0
    do b = 1, size(this % diagonal)
       associate (first => this % block_first(b), last => this % block_first(b + 1) - 1)
         m = last - first + 1
         allocate(block(m, m), source=0.0_dp)
         do i = first, last
            this % lower_start(i) = k + 1
            do e = start(this % permutation(i)), start(this % permutation(i) + 1) - 1
               ! A scaled entry of zero is no edge, and may lie above the
               ! blocks
               if (value(e) == 0.0_dp) cycle
               j = position(column(e))
               if (j >= first) then
                  block(i - first + 1, j - first + 1) = value(e)
               else
                  k = k + 1
                  this % lower_column(k) = j
                  this % lower_value(k)  = value(e)
               end if
            end do
         end do
         if (equilibrated) then
            call this % diagonal(b) % factorise(block, threshold)
         else
            call this % diagonal(b) % factorise(block, threshold, scales(this % permutation(first:last)))
         end if
         deallocate(block)
         if (this % diagonal(b) % is_singular) then
            this % is_singular = .true.
            return
         end if
       end associate
    end do
    this % lower_start(n + 1) = k + 1

  contains

    !-----------------------------------------------------------------!
    ! Sorts the entries of a row by their columns, by insertion, so
    ! that the entries of one column keep the order given
    !-----------------------------------------------------------------!

    pure subroutine ordered(c, v)
      integer , intent(inout) :: c(:)
      real(dp), intent(inout) :: v(:)
      integer                 :: p, q, inserted_column
      real(dp)                :: inserted_value
      do p = 2, size(c)
         inserted_column = c(p)
         inserted_value  = v(p)
         q = p - 1
         do while (q >= 1)
            if (c(q) <= inserted_column) exit
            c(q + 1) = c(q)
            v(q + 1) = v(q)
            q = q - 1
         end do
         c(q + 1) = inserted_column
         v(q + 1) = inserted_value
      end do
    end subroutine ordered

  end subroutine factorise_entries

  !===================================================================!
  ! Returns x with A x = b, or with A^T x = b where transposed. A right
  ! side of the wrong extent, or a factorisation left singular or never
  ! formed, stops the program.
  !===================================================================!

  subroutine substitute(this, b, x, transposed)

    class(dense_factorisation), intent(in)  :: this
    real(dp)                  , intent(in)  :: b(:)
    real(dp), allocatable     , intent(out) :: x(:)
    logical                   , intent(in)  :: transposed

    allocate(x(size(b)))
    call this % substituted(b, x, transposed)

  end subroutine substitute

  !===================================================================!
  ! Returns the same substitution into a vector already of the order
  ! of the matrix: a relaxation that solves one block per unknown block
  ! many times allocates nothing per solve.
  !
  ! In block triangular form the blocks are solved in order, or in
  ! reverse against the transpose, each right side less the blocks
  ! already solved:
  !
  !      forward      b = 1 .. nb     y_b = A_bb^-1 (y_b - sum over c < b of A_bc y_c)
  !      transpose    b = nb .. 1     y_b = A_bb^-T (y_b - sum over c > b of A_cb^T y_c)
  !
  ! From its entries, the sum against the transpose is scattered from
  ! each block as it is solved, into the right sides of the blocks
  ! before it.
  !===================================================================!

  recursive subroutine substituted(this, b, x, transposed)

    class(dense_factorisation), intent(in)  :: this
    real(dp)                  , intent(in)  :: b(:)
    real(dp)                  , intent(out) :: x(:)
    logical                   , intent(in)  :: transposed

    real(dp), allocatable :: y(:), part(:)
    real(dp)              :: stored
    integer               :: n, i, j, k, nb
    character(len=150)    :: message

    n = this % n

    if (n == 0) then
       error stop 'util_factorisation: substitute requires a prior factorise, but this % n = 0'
    end if
    if (this % is_singular) then
       error stop 'util_factorisation: substitute requires a non-singular factorisation, but &
            &this factorisation is singular'
    end if
    if (size(b) /= n .or. size(x) /= n) then
       write(message,'(a,i0,a,i0,a,i0)') 'util_factorisation: substitute requires the right side and &
            &the solution to match the matrix; size(b) = ', size(b), ', size(x) = ', size(x), ', n = ', n
       error stop trim(message)
    end if

    x = b
    if (allocated(this % row_scale)) x = b * merge(this % column_scale, this % row_scale, transposed)

    if (allocated(this % permutation)) then
       ! The blocks in order, or in reverse against the transpose, each
       ! right side less the blocks already solved
       nb = size(this % diagonal)
       y  = x(this % permutation)
       allocate(part(n))
       do k = merge(nb, 1, transposed), merge(1, nb, transposed), merge(-1, 1, transposed)
          associate (first => this % block_first(k), last => this % block_first(k + 1) - 1)
            if (allocated(this % lower_start)) then
               if (.not. transposed) then
                  do i = first, last
                     do j = this % lower_start(i), this % lower_start(i + 1) - 1
                        y(i) = y(i) - this % lower_value(j) * y(this % lower_column(j))
                     end do
                  end do
               end if
            else if (transposed) then
               if (last < n) y(first:last) = y(first:last) - matmul(y(last + 1:n), this % lu(last + 1:n, first:last))
            else
               if (first > 1) y(first:last) = y(first:last) - matmul(this % lu(first:last, 1:first - 1), y(1:first - 1))
            end if
            call this % diagonal(k) % substituted(y(first:last), part(first:last), transposed)
            y(first:last) = part(first:last)
            ! Against the transpose, from its entries, the block solved
            ! less from the right sides of the blocks before it
            if (allocated(this % lower_start) .and. transposed) then
               do i = first, last
                  do j = this % lower_start(i), this % lower_start(i + 1) - 1
                     y(this % lower_column(j)) = y(this % lower_column(j)) - this % lower_value(j) * y(i)
                  end do
               end do
            end if
          end associate
       end do
       x(this % permutation) = y
       if (allocated(this % row_scale)) x = x * merge(this % row_scale, this % column_scale, transposed)
       return
    end if

    if (.not. transposed) then

       do k = 1, n
          if (this % pivot_rows(k) /= k) then
             stored                  = x(k)
             x(k)                    = x(this % pivot_rows(k))
             x(this % pivot_rows(k)) = stored
          end if
       end do

       ! L y = b, then U x = y, each a sweep of columns: once x(j) is
       ! computed its column is subtracted from every entry below (or
       ! above)
       do j = 1, n - 1
          x(j+1:n) = x(j+1:n) - x(j) * this % lu(j+1:n, j)
       end do

       do j = n, 1, -1
          x(j) = x(j) / this % lu(j, j)
          x(1:j-1) = x(1:j-1) - x(j) * this % lu(1:j-1, j)
       end do

    else

       ! U^T y = b, then L^T z = y: the transposes run down columns as
       ! dot products, which are the same contiguous sweeps
       do i = 1, n
          x(i) = (x(i) - dot_product(this % lu(1:i-1, i), x(1:i-1))) / this % lu(i, i)
       end do

       do i = n - 1, 1, -1
          x(i) = x(i) - dot_product(this % lu(i+1:n, i), x(i+1:n))
       end do

       do k = n, 1, -1
          if (this % pivot_rows(k) /= k) then
             stored                  = x(k)
             x(k)                    = x(this % pivot_rows(k))
             x(this % pivot_rows(k)) = stored
          end if
       end do

    end if
    if (allocated(this % row_scale)) x = x * merge(this % row_scale, this % column_scale, transposed)

  end subroutine substituted

  !===================================================================!
  ! Returns the strongly connected components of the graph of the
  ! pattern of a, an edge from i to j where a_ij is not zero and
  ! i /= j, by Tarjan's algorithm: the vertices in the order the
  ! components are emitted, each component after every component it
  ! reaches, and the first of each component in that order, with n + 1
  ! after the last.
  !
  !      a pattern          its graph               the permutation
  !      | x . . |          1 <---- 2               [ 1 | 3 2 ]
  !      | x x x |                 ^ |              block_first = 1 2 4
  !      | . x x |                 | v
  !                                 3
  !===================================================================!

  subroutine strong_components(a, permutation, block_first)

    real(dp)            , intent(in)  :: a(:,:)
    integer, allocatable, intent(out) :: permutation(:), block_first(:)

    integer, allocatable :: count(:), start(:), edge(:)
    integer              :: n, i, j

    n = size(a, 1)
    ! The edges of each vertex, in compressed rows, from a sweep of the
    ! columns
    allocate(count(n), source=0)
    do j = 1, n
       do i = 1, n
          if (i /= j .and. a(i, j) /= 0.0_dp) count(i) = count(i) + 1
       end do
    end do
    allocate(start(n + 1))
    start(1) = 1
    do i = 1, n
       start(i + 1) = start(i) + count(i)
    end do
    allocate(edge(start(n + 1) - 1))
    count = 0
    do j = 1, n
       do i = 1, n
          if (i /= j .and. a(i, j) /= 0.0_dp) then
             edge(start(i) + count(i)) = j
             count(i) = count(i) + 1
          end if
       end do
    end do
    call components_of_edges(start, edge, permutation, block_first)

  end subroutine strong_components

  !===================================================================!
  ! Returns the strongly connected components of the graph whose edges
  ! from vertex i are edge(start(i):start(i + 1) - 1), as
  ! strong_components states them. The recursion of Tarjan's algorithm
  ! runs on a stack of its own:
  !
  !      visit(1:depth)    the vertices of the search, root first
  !      next(depth)       the next edge of visit(depth) to follow
  !      path(1:top)       the vertices not yet in a component, the
  !                        stack of the algorithm
  !      index, low        the number of each vertex in the visit, and
  !                        the least number it reaches on the path
  !
  ! A vertex whose low is its own index closes a component: the path
  ! is popped down to it.
  !===================================================================!

  subroutine components_of_edges(start, edge, permutation, block_first)

    integer             , intent(in)  :: start(:), edge(:)
    integer, allocatable, intent(out) :: permutation(:), block_first(:)

    integer, allocatable :: index(:), low(:), path(:), next(:), visit(:)
    logical, allocatable :: on_path(:)
    integer              :: n, v, w, top, depth, number, filled, root, blocks

    n = size(start) - 1
    allocate(index(n), low(n), source=0)
    allocate(on_path(n), source=.false.)
    allocate(path(n), visit(n), next(n))
    allocate(permutation(n), block_first(n + 1))
    blocks = 0
    top    = 0
    depth  = 0
    number = 0
    filled = 0
    do root = 1, n
       if (index(root) /= 0) cycle
       v     = root
       depth = 1
       visit(depth) = v
       next(depth)  = start(v)
       number     = number + 1
       index(v)   = number
       low(v)     = number
       top        = top + 1
       path(top)  = v
       on_path(v) = .true.
       do while (depth > 0)
          v = visit(depth)
          if (next(depth) < start(v + 1)) then
             w = edge(next(depth))
             next(depth) = next(depth) + 1
             if (index(w) == 0) then
                depth = depth + 1
                visit(depth) = w
                next(depth)  = start(w)
                number     = number + 1
                index(w)   = number
                low(w)     = number
                top        = top + 1
                path(top)  = w
                on_path(w) = .true.
             else if (on_path(w)) then
                low(v) = min(low(v), index(w))
             end if
             cycle
          end if
          ! Every edge of v visited: v closes a component where it is
          ! its own root
          if (low(v) == index(v)) then
             blocks = blocks + 1
             block_first(blocks) = filled + 1
             do
                w = path(top)
                top = top - 1
                on_path(w) = .false.
                filled = filled + 1
                permutation(filled) = w
                if (w == v) exit
             end do
          end if
          depth = depth - 1
          if (depth > 0) low(visit(depth)) = min(low(visit(depth)), low(v))
       end do
    end do
    block_first(blocks + 1) = n + 1
    block_first = block_first(1:blocks + 1)

  end subroutine components_of_edges

  !===================================================================!
  ! Returns the order of the matrix factorised; zero before a
  ! factorisation
  !===================================================================!

  pure integer function order(this) result(n)

    class(dense_factorisation), intent(in) :: this

    n = this % n

  end function order

  !===================================================================!
  ! Returns whether the factorisation was left singular
  !===================================================================!

  pure logical function singular(this) result(is_rank_deficient)

    class(dense_factorisation), intent(in) :: this

    is_rank_deficient = this % is_singular

  end function singular

end module util_factorisation
