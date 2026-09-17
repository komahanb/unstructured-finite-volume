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
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module util_factorisation

  use util_precision  , only : dp

  implicit none

  private
  public :: dense_factorisation

  type :: dense_factorisation

     real(dp), allocatable, private :: lu(:,:)
     integer , allocatable, private :: pivot_rows(:)
     integer , private :: n = 0
     logical , private :: is_singular = .false.

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
  ! the factorisation singular rather than dividing by it.
  !===================================================================!

  subroutine factorise(this, a, threshold)

    class(dense_factorisation), intent(inout) :: this
    real(dp)                  , intent(in)    :: a(:,:)
    real(dp)                  , intent(in)    :: threshold

    real(dp), allocatable :: row(:)
    real(dp) :: factor
    integer  :: n, k, p, j
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

       if (abs(this % lu(k, k)) <= threshold) then
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

  subroutine substituted(this, b, x, transposed)

    class(dense_factorisation), intent(in)  :: this
    real(dp)                  , intent(in)  :: b(:)
    real(dp)                  , intent(out) :: x(:)
    logical                   , intent(in)  :: transposed

    real(dp) :: stored
    integer  :: n, i, j, k
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

  end subroutine substituted

  pure integer function order(this) result(n)

    class(dense_factorisation), intent(in) :: this

    n = this % n

  end function order

  pure logical function singular(this) result(is_rank_deficient)

    class(dense_factorisation), intent(in) :: this

    is_rank_deficient = this % is_singular

  end function singular

end module util_factorisation
