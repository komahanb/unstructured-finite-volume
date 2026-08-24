!=====================================================================!
! A number together with every mixed derivative of it along n given
! directions.
!
! A quantity built from the four arithmetic operations alone carries
! its derivatives with it if each operation is given the rule for
! them. That is what this type holds: one coefficient per subset of
! the n directions, 2^n in all, the empty subset holding the value
! and the full subset holding the n-th mixed partial. Nothing is
! perturbed and nothing is truncated - each coefficient is the exact
! derivative.
!
!             THE SUBSET INDEX
!
! A subset is a bit mask, and the mask m is held at index m + 1.
! Allocatable assignment does not carry a zero lower bound in the
! compilers this is built with, so the masks are stored from one and
! the two are never confused.
!
!             SYMMETRIC SEEDING
!
! Set every subset of the same size to the same number and the type
! computes a Taylor composition rather than a mixed partial: seed the
! subsets of size k with the k-th derivative of a quantity along one
! parameter, and the coefficient of the full subset comes out as the
! n-th derivative of whatever was built from it. That is the product
! rule read on subsets - splitting a set of size k into two parts
! counts each split once, which is the binomial coefficient Leibniz
! asks for - so nothing here changes; only what the numbers are taken
! to mean does.
!
! The convention is derivatives, not derivatives over factorials: a
! subset of size k holds the k-th derivative itself.
!
!             THE FOUR RULES
!
!      sum        coefficient by coefficient
!      product    the coefficient of a subset m is the sum, over the
!                 ways of splitting m into two disjoint parts, of the
!                 product of the factors' coefficients on those parts
!      quotient   r = a / b read from r b = a subset by subset in
!                 increasing order, each one solved for r
!      power      repeated product, and for a negative exponent the
!                 reciprocal of the positive power
!
! A product costs up to 3^n operations and a quotient the same, which
! is what a mixed partial of degree n costs; the storage is 2^n
! numbers per quantity. The degree is bounded by the width a default
! integer can carry a mask in, which max_subset_width reports.
!
!             WHAT IS REFUSED
!
! A quotient whose divisor has a zero value, and a direction index
! outside one to n. Each stops the program: neither has a derivative
! to report.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module util_derivative_terms

  use iso_fortran_env, only : dp => REAL64

  implicit none

  private
  public :: derivative_terms, value, mixed_partial, max_subset_width
  public :: integer_power
  public :: operator(+), operator(-), operator(*), operator(/)

  type :: derivative_terms

     integer               , private              :: directions = 0
     real(dp)              , private, allocatable :: terms(:)

   contains

     procedure :: set_direction
     procedure :: set_symmetric
     procedure :: num_directions

  end type derivative_terms

  interface derivative_terms
     module procedure create_constant
     module procedure create_like
  end interface derivative_terms

  interface operator(+)
     module procedure terms_plus
  end interface operator(+)

  interface operator(-)
     module procedure terms_minus
     module procedure terms_negated
  end interface operator(-)

  interface operator(*)
     module procedure terms_times
     module procedure real_times_terms
  end interface operator(*)

  interface operator(/)
     module procedure terms_over
  end interface operator(/)

contains

  !===================================================================!
  ! The widest mask a default integer indexes: the bit width less the
  ! sign bit and the bit that 2 raised to the width would carry.
  !===================================================================!

  pure integer function max_subset_width()

    max_subset_width = bit_size(0) - 2

  end function max_subset_width

  !===================================================================!
  ! A constant: its value, and zero in every derivative. A negative
  ! count of directions, or one past the mask width, stops the
  ! program.
  !===================================================================!

  pure function create_constant(x, num_directions) result(this)

    real(dp), intent(in) :: x
    integer , intent(in) :: num_directions
    type(derivative_terms) :: this

    if (num_directions < 0 .or. num_directions > max_subset_width()) then
       error stop 'util_derivative_terms: the count of directions is within the mask width'
    end if

    this % directions = num_directions
    allocate(this % terms(2**num_directions), source=0.0_dp)
    this % terms(1) = x

  end function create_constant

  pure function create_like(x, other) result(this)

    real(dp)              , intent(in) :: x
    type(derivative_terms), intent(in) :: other
    type(derivative_terms) :: this

    this = create_constant(x, other % directions)

  end function create_like

  pure integer function num_directions(this)

    class(derivative_terms), intent(in) :: this

    num_directions = this % directions

  end function num_directions

  !===================================================================!
  ! The derivative of this quantity along direction i, which is the
  ! coefficient of the subset holding i alone. An index outside one
  ! to n stops the program.
  !===================================================================!

  pure subroutine set_direction(this, i, x)

    class(derivative_terms), intent(inout) :: this
    integer                , intent(in)    :: i
    real(dp)               , intent(in)    :: x

    if (i < 1 .or. i > this % directions) then
       error stop 'util_derivative_terms: the direction is one of those declared'
    end if

    this % terms(2**(i - 1) + 1) = x

  end subroutine set_direction

  !===================================================================!
  ! Seed every subset of one size with one number. An order below
  ! zero or past the directions declared stops the program: there is
  ! no such coefficient to set.
  !===================================================================!

  pure subroutine set_symmetric(this, order, x)

    class(derivative_terms), intent(inout) :: this
    integer                , intent(in)    :: order
    real(dp)               , intent(in)    :: x

    integer :: m

    if (order < 0 .or. order > this % directions) then
       error stop 'util_derivative_terms: the order is one the directions carry'
    end if

    do m = 0, size(this % terms) - 1
       if (popcnt(m) == order) this % terms(m + 1) = x
    end do

  end subroutine set_symmetric

  pure real(dp) function value(x)

    type(derivative_terms), intent(in) :: x

    value = x % terms(1)

  end function value

  pure real(dp) function mixed_partial(x)

    type(derivative_terms), intent(in) :: x

    mixed_partial = x % terms(size(x % terms))

  end function mixed_partial

  !===================================================================!
  ! THE FOUR RULES.
  !===================================================================!

  pure function terms_plus(a, b) result(c)

    type(derivative_terms), intent(in) :: a, b
    type(derivative_terms) :: c

    c % directions = a % directions
    c % terms = a % terms + b % terms

  end function terms_plus

  pure function terms_minus(a, b) result(c)

    type(derivative_terms), intent(in) :: a, b
    type(derivative_terms) :: c

    c % directions = a % directions
    c % terms = a % terms - b % terms

  end function terms_minus

  pure function terms_negated(a) result(c)

    type(derivative_terms), intent(in) :: a
    type(derivative_terms) :: c

    c % directions = a % directions
    c % terms = -a % terms

  end function terms_negated

  pure function real_times_terms(x, a) result(c)

    real(dp)              , intent(in) :: x
    type(derivative_terms), intent(in) :: a
    type(derivative_terms) :: c

    c % directions = a % directions
    c % terms = x * a % terms

  end function real_times_terms

  pure function terms_times(a, b) result(c)

    type(derivative_terms), intent(in) :: a, b
    type(derivative_terms) :: c

    integer :: m, s

    c = create_like(0.0_dp, a)

    do m = 0, size(a % terms) - 1
       s = m
       do
          c % terms(m + 1) = c % terms(m + 1) &
               & + a % terms(s + 1) * b % terms(ieor(m, s) + 1)
          if (s == 0) exit
          s = iand(s - 1, m)
       end do
    end do

  end function terms_times

  pure function terms_over(a, b) result(r)

    type(derivative_terms), intent(in) :: a, b
    type(derivative_terms) :: r

    real(dp) :: accumulated
    integer  :: m, s

    if (b % terms(1) == 0.0_dp) then
       error stop 'util_derivative_terms: a quotient divides by a nonzero value'
    end if

    r = create_like(0.0_dp, a)

    do m = 0, size(a % terms) - 1
       accumulated = a % terms(m + 1)
       s = m
       do while (s /= 0)
          accumulated = accumulated - b % terms(s + 1) * r % terms(ieor(m, s) + 1)
          s = iand(s - 1, m)
       end do
       r % terms(m + 1) = accumulated / b % terms(1)
    end do

  end function terms_over

  !===================================================================!
  ! x raised to an integer exponent, of either sign. A zero exponent
  ! is one, whatever the value; a negative exponent divides, and a
  ! zero value then stops the program inside the quotient.
  !===================================================================!

  pure function integer_power(x, n) result(p)

    type(derivative_terms), intent(in) :: x
    integer               , intent(in) :: n
    type(derivative_terms) :: p

    integer :: i

    p = create_like(1.0_dp, x)

    do i = 1, abs(n)
       p = p * x
    end do

    if (n < 0) p = create_like(1.0_dp, x) / p

  end function integer_power

end module util_derivative_terms
