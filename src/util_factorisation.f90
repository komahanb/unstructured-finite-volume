!=====================================================================!
! A factorisation kept: P A = L U with partial pivoting, formed once
! from a square matrix and applied to as many right sides as are
! brought to it, against A or against its transpose.
!
! The point of keeping it is the count. A tangent solve, an adjoint
! solve and every order of an expansion read the same matrix frozen
! at the same state, so one factorisation serves them all and each
! costs a substitution. Formed afresh per solve, each would cost the
! factorisation over again.
!
!=====================================================================!
!
!                        THE COST MODEL
!
! Every procedure that costs carries its own count, so a caller
! choosing between routes can add rather than guess:
!
!      factorise      n^3 / 3      multiplications
!      substitute     n^2          multiplications, either way round
!
! and the count of each is filed with util_tally, so what a run did
! can be read beside what the model said it would.
!
!=====================================================================!
!
! Solving against the transpose uses the same factors. P A = L U gives
! A^T P^T = U^T L^T, so A^T x = b is U^T y = b, then L^T z = y, then
! x = P^T z, which is the row exchanges of the factorisation applied
! in the reverse of the order they were made.
!
! A pivot at or below the threshold given leaves the factorisation
! singular: no further elimination is done, and a substitution
! against it stops the program. Whether a singular pivot is a fault or
! a fact is the caller's to decide, from singular().
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module util_factorisation

  use iso_fortran_env, only : dp => REAL64
  use util_tally     , only : tally_record, factorisations, linear_solves

  implicit none

  private
  public :: dense_factorisation

  type :: dense_factorisation

     real(dp), allocatable, private :: lu(:,:)
     integer , allocatable, private :: exchanged(:)
     integer , private :: n = 0
     logical , private :: is_singular = .false.
     real(dp), private :: least_pivot = 0.0_dp

   contains

     procedure :: factorise
     procedure :: substitute
     procedure :: order
     procedure :: singular
     procedure :: smallest_pivot
     procedure, nopass :: factorise_cost
     procedure, nopass :: substitute_cost

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
    integer  :: n, k, p, i

    n = size(a, 1)
    if (size(a, 2) /= n) then
       error stop 'util_factorisation: a factorised matrix is square'
    end if

    call tally_record(factorisations)

    this % n           = n
    this % lu          = a
    this % is_singular = .false.
    this % least_pivot = huge(1.0_dp)
    allocate(this % exchanged(n), source=0)
    allocate(row(n))

    do k = 1, n

       p = k - 1 + maxloc(abs(this % lu(k:n, k)), dim=1)
       this % exchanged(k) = p

       if (p /= k) then
          row              = this % lu(k, :)
          this % lu(k, :)  = this % lu(p, :)
          this % lu(p, :)  = row
       end if

       this % least_pivot = min(this % least_pivot, abs(this % lu(k, k)))

       if (abs(this % lu(k, k)) <= threshold) then
          this % is_singular = .true.
          return
       end if

       do i = k + 1, n
          factor                 = this % lu(i, k) / this % lu(k, k)
          this % lu(i, k)        = factor
          this % lu(i, k+1:n)    = this % lu(i, k+1:n) - factor * this % lu(k, k+1:n)
       end do

    end do

  end subroutine factorise

  !===================================================================!
  ! x with A x = b, or with A^T x = b where transposed. A right side of
  ! the wrong extent, or a factorisation left singular or never made,
  ! stops the program.
  !===================================================================!

  subroutine substitute(this, b, x, transposed)

    class(dense_factorisation), intent(in)  :: this
    real(dp)                  , intent(in)  :: b(:)
    real(dp), allocatable     , intent(out) :: x(:)
    logical                   , intent(in)  :: transposed

    real(dp) :: held
    integer  :: n, i, k

    n = this % n

    if (n == 0) then
       error stop 'util_factorisation: a substitution follows a factorisation'
    end if
    if (this % is_singular) then
       error stop 'util_factorisation: a singular factorisation is not substituted against'
    end if
    if (size(b) /= n) then
       error stop 'util_factorisation: the right side matches the matrix'
    end if

    call tally_record(linear_solves)

    x = b

    if (.not. transposed) then

       do k = 1, n
          if (this % exchanged(k) /= k) then
             held                    = x(k)
             x(k)                    = x(this % exchanged(k))
             x(this % exchanged(k))  = held
          end if
       end do

       do i = 2, n
          x(i) = x(i) - dot_product(this % lu(i, 1:i-1), x(1:i-1))
       end do

       do i = n, 1, -1
          x(i) = (x(i) - dot_product(this % lu(i, i+1:n), x(i+1:n))) / this % lu(i, i)
       end do

    else

       do i = 1, n
          x(i) = (x(i) - dot_product(this % lu(1:i-1, i), x(1:i-1))) / this % lu(i, i)
       end do

       do i = n - 1, 1, -1
          x(i) = x(i) - dot_product(this % lu(i+1:n, i), x(i+1:n))
       end do

       do k = n, 1, -1
          if (this % exchanged(k) /= k) then
             held                    = x(k)
             x(k)                    = x(this % exchanged(k))
             x(this % exchanged(k))  = held
          end if
       end do

    end if

  end subroutine substitute

  pure integer function order(this) result(n)

    class(dense_factorisation), intent(in) :: this

    n = this % n

  end function order

  pure logical function singular(this) result(yes)

    class(dense_factorisation), intent(in) :: this

    yes = this % is_singular

  end function singular

  pure real(dp) function smallest_pivot(this) result(least)

    class(dense_factorisation), intent(in) :: this

    least = this % least_pivot

  end function smallest_pivot

  !===================================================================!
  ! THE COST MODEL, in multiplications.
  !===================================================================!

  pure real(dp) function factorise_cost(n) result(cost)

    integer, intent(in) :: n

    cost = real(n, dp) ** 3 / 3.0_dp

  end function factorise_cost

  pure real(dp) function substitute_cost(n) result(cost)

    integer, intent(in) :: n

    cost = real(n, dp) ** 2

  end function substitute_cost

end module util_factorisation
