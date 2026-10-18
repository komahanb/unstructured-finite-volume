!=====================================================================!
! A factorisation retained: P A = L U with partial pivoting, formed
! once from a square matrix and applied to every right side passed to
! it, against A or against its transpose.
!
! The purpose of retaining it is the operation count. A tangent solve,
! an adjoint solve and every order of an expansion read the same
! matrix frozen at the same state, so one factorisation serves them
! all and each costs a substitution. Formed again per solve, each
! would cost the factorisation again. Each factorisation is recorded
! by its caller, so the count is measured rather than modelled.
!
! Solving against the transpose uses the same factors. P A = L U gives
! A^T P^T = U^T L^T, so A^T x = b is U^T y = b, then L^T z = y, then
! x = P^T z, which is the row exchanges of the factorisation applied
! in the reverse of the order they were made.
!
! A pivot at or below the threshold given leaves the factorisation
! singular: no further elimination is done, and a substitution
! against it stops the program. Whether a singular pivot is an error
! or an expected result is the caller's decision, read from singular().
!
! EQUILIBRATED, the matrix factorised is D_r A D_c, its rows and then
! its columns scaled to a largest entry of one, and a substitution
! scales its right side and its solution by them: x = D_c (D_r A
! D_c)^-1 D_r b, and against the transpose x = D_r (D_c A^T D_r)^-1
! D_c b. A system of rows of different units - positions, rotations
! and reactions - is solved to the digits its scaled matrix allows.
!
! BLOCK TRIANGULAR, the matrix is permuted, rows and columns alike, to
! the block lower triangular form its zero pattern allows, and only
! its diagonal blocks are factorised. The blocks are the strongly
! connected components of the graph of the pattern, an edge from i to
! j where a_ij is not zero, found by Tarjan's algorithm. It emits a
! component after every component it reads, which is the order of a
! forward substitution: x_b = A_bb^-1 (b_b - sum over c before b of
! A_bc x_c). Against the transpose the blocks are solved in reverse.
! The matrix is nonsingular exactly where every diagonal block is.
! The step of a family of s stages is such a matrix, each stage
! reading the stages before it: its factorisation costs s blocks of
! one stage, about s^2 times less than the whole. A row is the law of
! the unknown of the same index, so the permutation is the same on
! both sides.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module util_factorisation

  use util_precision  , only : dp

  implicit none

  private
  public :: dense_factorisation

  type :: dense_factorisation

     real(dp), allocatable, private :: lu(:,:)
     real(dp), allocatable, private :: row_scale(:), column_scale(:)
     integer , allocatable, private :: pivot_rows(:)
     integer , private :: n = 0
     logical , private :: is_singular = .false.
     ! block triangular: the unknowns in the order of the blocks, the
     ! first of each block in that order, the permuted matrix, and the
     ! factorisation of each diagonal block
     integer , allocatable, private :: permutation(:), block_first(:)
     type(dense_factorisation), allocatable, private :: diagonal(:)

   contains

     procedure :: factorise
     procedure :: substitute
     procedure :: substituted
     procedure :: order
     procedure :: singular

  end type dense_factorisation

contains

  !===================================================================!
  ! L U of a square matrix, in place in a copy. A matrix that is not
  ! square stops the program. A pivot at or below the threshold leaves
  ! the factorisation singular rather than dividing by it; given the
  ! scale of each column, at or below the threshold times its column's
  ! scale, which equilibrated is one.
  !===================================================================!

  recursive subroutine factorise(this, a, threshold, scales, equilibrated, block_triangular)

    class(dense_factorisation), intent(inout) :: this
    real(dp)                  , intent(in)    :: a(:,:)
    real(dp)                  , intent(in)    :: threshold
    real(dp)                  , intent(in), optional :: scales(:)
    logical                   , intent(in), optional :: equilibrated, block_triangular

    real(dp), allocatable :: row(:)
    real(dp) :: factor
    integer  :: n, k, p, j, b
    character(len=150) :: message

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
          this % row_scale = [(1.0_dp / max(tiny(1.0_dp), maxval(abs(a(k, :)))), k = 1, n)]
          do k = 1, n
             this % lu(k, :) = this % lu(k, :) * this % row_scale(k)
          end do
          this % column_scale = [(1.0_dp / max(tiny(1.0_dp), maxval(abs(this % lu(:, k)))), k = 1, n)]
          do k = 1, n
             this % lu(:, k) = this % lu(:, k) * this % column_scale(k)
          end do
       end if
    end if
    if (allocated(this % permutation)) deallocate(this % permutation, this % block_first, this % diagonal)
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
    ! column k are formed at once, and each later column is updated
    ! as a whole, so every inner sweep runs down a column, which is
    ! the array's memory order.
    do k = 1, n

       p = k - 1 + maxloc(abs(this % lu(k:n, k)), dim=1)
       this % pivot_rows(k) = p

       if (p /= k) then
          row              = this % lu(k, :)
          this % lu(k, :)  = this % lu(p, :)
          this % lu(p, :)  = row
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
  ! x with A x = b, or with A^T x = b where transposed. A right side of
  ! the wrong extent, or a factorisation left singular or never formed,
  ! stops the program.
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
  ! The same substitution into a vector already of the order of the
  ! matrix: a relaxation that solves one block per unknown block many
  ! times allocates nothing per solve.
  !===================================================================!

  recursive subroutine substituted(this, b, x, transposed)

    class(dense_factorisation), intent(in)  :: this
    real(dp)                  , intent(in)  :: b(:)
    real(dp)                  , intent(out) :: x(:)
    logical                   , intent(in)  :: transposed

    real(dp), allocatable :: y(:), part(:)
    real(dp) :: stored
    integer  :: n, i, j, k, nb
    character(len=150) :: message

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
       ! the blocks in order, or in reverse against the transpose, each
       ! right side less the blocks already solved
       nb = size(this % diagonal)
       y  = x(this % permutation)
       allocate(part(n))
       do k = merge(nb, 1, transposed), merge(1, nb, transposed), merge(-1, 1, transposed)
          associate (first => this % block_first(k), last => this % block_first(k + 1) - 1)
            if (transposed) then
               if (last < n) y(first:last) = y(first:last) - matmul(y(last + 1:n), this % lu(last + 1:n, first:last))
            else
               if (first > 1) y(first:last) = y(first:last) - matmul(this % lu(first:last, 1:first - 1), y(1:first - 1))
            end if
            call this % diagonal(k) % substituted(y(first:last), part(first:last), transposed)
            y(first:last) = part(first:last)
          end associate
       end do
       x(this % permutation) = y
       if (allocated(this % row_scale)) x = x * merge(this % row_scale, this % column_scale, transposed)
       return
    end if

    if (.not. transposed) then

       do k = 1, n
          if (this % pivot_rows(k) /= k) then
             stored                    = x(k)
             x(k)                    = x(this % pivot_rows(k))
             x(this % pivot_rows(k))  = stored
          end if
       end do

       ! L y = b, then U x = y, each a sweep of columns: once x(j) is
       ! computed its column is subtracted from every entry below (or above).
       do j = 1, n - 1
          x(j+1:n) = x(j+1:n) - x(j) * this % lu(j+1:n, j)
       end do

       do j = n, 1, -1
          x(j) = x(j) / this % lu(j, j)
          x(1:j-1) = x(1:j-1) - x(j) * this % lu(1:j-1, j)
       end do

    else

       ! U^T y = b, then L^T z = y: the transposes run down columns
       ! as dot products, which are the same contiguous sweeps.
       do i = 1, n
          x(i) = (x(i) - dot_product(this % lu(1:i-1, i), x(1:i-1))) / this % lu(i, i)
       end do

       do i = n - 1, 1, -1
          x(i) = x(i) - dot_product(this % lu(i+1:n, i), x(i+1:n))
       end do

       do k = n, 1, -1
          if (this % pivot_rows(k) /= k) then
             stored                    = x(k)
             x(k)                    = x(this % pivot_rows(k))
             x(this % pivot_rows(k))  = stored
          end if
       end do

    end if
    if (allocated(this % row_scale)) x = x * merge(this % row_scale, this % column_scale, transposed)

  end subroutine substituted

  !===================================================================!
  ! THE STRONGLY CONNECTED COMPONENTS of the graph of the pattern of a,
  ! an edge from i to j where a_ij is not zero and i /= j, by Tarjan's
  ! algorithm, its recursion on a stack of its own. The vertices
  ! in the order the components are emitted, each component after
  ! every component it reaches, and the first of each component in
  ! that order, with n + 1 after the last.
  !===================================================================!

  subroutine strong_components(a, permutation, block_first)

    real(dp)            , intent(in)  :: a(:,:)
    integer, allocatable, intent(out) :: permutation(:), block_first(:)

    integer, allocatable :: count(:), start(:), edge(:), index(:), low(:), path(:), next(:), visit(:)
    logical, allocatable :: on_path(:)
    integer :: n, i, j, v, w, top, depth, number, filled, root

    n = size(a, 1)
    ! the edges of each vertex, in compressed rows, from a sweep of the
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

    allocate(index(n), low(n), source=0)
    allocate(on_path(n), source=.false.)
    allocate(path(n), visit(n), next(n))
    allocate(permutation(n), block_first(0))
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
       number   = number + 1
       index(v) = number
       low(v)   = number
       top      = top + 1
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
                number   = number + 1
                index(w) = number
                low(w)   = number
                top      = top + 1
                path(top)  = w
                on_path(w) = .true.
             else if (on_path(w)) then
                low(v) = min(low(v), index(w))
             end if
             cycle
          end if
          ! every edge of v visited: v closes a component where it is
          ! its own root
          if (low(v) == index(v)) then
             block_first = [block_first, filled + 1]
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
    block_first = [block_first, n + 1]

  end subroutine strong_components

  pure integer function order(this) result(n)

    class(dense_factorisation), intent(in) :: this

    n = this % n

  end function order

  pure logical function singular(this) result(is_rank_deficient)

    class(dense_factorisation), intent(in) :: this

    is_rank_deficient = this % is_singular

  end function singular

end module util_factorisation
