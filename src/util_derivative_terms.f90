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
!             COMPOSITION WITH A FUNCTION
!
! For g = f(a) the coefficient of a subset m is the sum over the set
! partitions of m into blocks B_1 .. B_p of
!
!      f^(p)(a_0)  a_{B_1} ... a_{B_p}
!
! computed without listing partitions. With i the lowest element of m
! and G(k, m) the same sum with f^(p+k) in place of f^(p), the block
! holding i is S + i for some S within m - i and the rest is
! partitioned with one more derivative taken:
!
!      G(k, empty)  =  f^(k)(a_0)
!      G(k, m)      =  sum over S within m - i  of  a_{S+i} G(k+1, m - i - S)
!
! filled in increasing m, for k up to n - |m|; the coefficient of g on
! m is G(0, m). Each entry is a subset sum, so the cost is the 3^n a
! product pays. A function supplies its derivatives at the value of a
! for k = 0 .. n: exp is its own; sin and cos cycle through four; log
! and a real power recur by one division of the value each, so no
! factorial is formed. sqrt is the power one half, and an integer
! power stays a repeated product, which a zero value admits.
!
! A product costs up to 3^n operations and a quotient the same, which
! is what a mixed partial of degree n costs; the storage is 2^n
! numbers per quantity. The degree is bounded by the width a default
! integer can carry a mask in, which max_subset_width reports.
!
!             WHAT IS REFUSED
!
! A quotient whose divisor has a zero value, and a direction index
! outside one to n. log, sqrt and a real power at a value that is not
! positive, and a derivative table shorter than n + 1. Each stops the
! program: there is no derivative to report.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module util_derivative_terms

  use util_precision  , only : dp

  implicit none

  private
  public :: derivative_terms, value, mixed_partial, coefficient, max_subset_width
  public :: integer_power, composed
  public :: operator(+), operator(-), operator(*), operator(/), operator(**)
  public :: sin, cos, exp, log, sqrt

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

  interface operator(**)
     module procedure terms_to_integer
     module procedure terms_to_real
  end interface operator(**)

  interface sin
     module procedure terms_sin
  end interface sin

  interface cos
     module procedure terms_cos
  end interface cos

  interface exp
     module procedure terms_exp
  end interface exp

  interface log
     module procedure terms_log
  end interface log

  interface sqrt
     module procedure terms_sqrt
  end interface sqrt

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
  ! The coefficient of the subset with mask m: the mixed partial along
  ! the directions m holds. A mask outside 0 to 2^n - 1 stops the
  ! program.
  !===================================================================!

  pure real(dp) function coefficient(x, m)

    type(derivative_terms), intent(in) :: x
    integer               , intent(in) :: m

    if (m < 0 .or. m >= size(x % terms)) then
       error stop 'util_derivative_terms: the mask names a subset of the directions'
    end if

    coefficient = x % terms(m + 1)

  end function coefficient

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

  pure function terms_to_integer(x, n) result(p)

    type(derivative_terms), intent(in) :: x
    integer               , intent(in) :: n
    type(derivative_terms) :: p

    p = integer_power(x, n)

  end function terms_to_integer

  !===================================================================!
  ! COMPOSITION WITH A FUNCTION: g = f(a) from the derivatives f(0:n)
  ! of f at the value of a. A table shorter than n + 1 stops the
  ! program: the highest subset has no derivative to draw on.
  !===================================================================!

  pure function composed(a, f) result(g)

    type(derivative_terms), intent(in) :: a
    real(dp)              , intent(in) :: f(0:)
    type(derivative_terms) :: g

    real(dp), allocatable :: table(:,:)
    real(dp) :: accumulated
    integer  :: n, m, k, i, rest, s

    n = a % directions

    if (ubound(f, 1) < n) then
       error stop 'util_derivative_terms: a function supplies n + 1 derivatives'
    end if

    allocate(table(0:n, 0:2**n - 1))
    table(:, 0) = f(0:n)

    do m = 1, 2**n - 1
       i    = ibset(0, trailz(m))
       rest = ieor(m, i)
       do k = 0, n - popcnt(m)
          accumulated = 0.0_dp
          s = rest
          do
             accumulated = accumulated + a % terms(ior(s, i) + 1) * table(k + 1, ieor(rest, s))
             if (s == 0) exit
             s = iand(s - 1, rest)
          end do
          table(k, m) = accumulated
       end do
    end do

    g = create_like(0.0_dp, a)
    do m = 0, 2**n - 1
       g % terms(m + 1) = table(0, m)
    end do

  end function composed

  !===================================================================!
  ! THE DERIVATIVE TABLES, each f^(k) at x for k = 0 .. n.
  !===================================================================!

  pure function terms_exp(a) result(g)

    type(derivative_terms), intent(in) :: a
    type(derivative_terms) :: g

    real(dp) :: f(0:a % directions)

    f = exp(a % terms(1))
    g = composed(a, f)

  end function terms_exp

  pure function terms_sin(a) result(g)

    type(derivative_terms), intent(in) :: a
    type(derivative_terms) :: g

    real(dp) :: f(0:a % directions), cycle(0:3)
    integer  :: k

    cycle = [sin(a % terms(1)), cos(a % terms(1)), -sin(a % terms(1)), -cos(a % terms(1))]
    do k = 0, a % directions
       f(k) = cycle(mod(k, 4))
    end do
    g = composed(a, f)

  end function terms_sin

  pure function terms_cos(a) result(g)

    type(derivative_terms), intent(in) :: a
    type(derivative_terms) :: g

    real(dp) :: f(0:a % directions), cycle(0:3)
    integer  :: k

    cycle = [cos(a % terms(1)), -sin(a % terms(1)), -cos(a % terms(1)), sin(a % terms(1))]
    do k = 0, a % directions
       f(k) = cycle(mod(k, 4))
    end do
    g = composed(a, f)

  end function terms_cos

  !-------------------------------------------------------------------!
  ! log: f^(k) = (-1)^(k-1) (k-1)! / x^k, formed as f^(k) = -f^(k-1)
  ! (k-1) / x from f^(1) = 1 / x. A value that is not positive stops
  ! the program.
  !-------------------------------------------------------------------!

  pure function terms_log(a) result(g)

    type(derivative_terms), intent(in) :: a
    type(derivative_terms) :: g

    real(dp) :: f(0:a % directions), x
    integer  :: k

    x = a % terms(1)
    if (x <= 0.0_dp) then
       error stop 'util_derivative_terms: log is taken at a positive value'
    end if

    f(0) = log(x)
    if (a % directions >= 1) f(1) = 1.0_dp / x
    do k = 2, a % directions
       f(k) = -f(k - 1) * real(k - 1, dp) / x
    end do
    g = composed(a, f)

  end function terms_log

  !-------------------------------------------------------------------!
  ! A real power: f^(k) = p (p-1) .. (p-k+1) x^(p-k), formed as
  ! f^(k) = f^(k-1) (p-k+1) / x from f^(0) = x^p. A value that is not
  ! positive stops the program; an integer exponent is a repeated
  ! product and admits any value.
  !-------------------------------------------------------------------!

  pure function terms_to_real(a, p) result(g)

    type(derivative_terms), intent(in) :: a
    real(dp)              , intent(in) :: p
    type(derivative_terms) :: g

    real(dp) :: f(0:a % directions), x
    integer  :: k

    x = a % terms(1)
    if (x <= 0.0_dp) then
       error stop 'util_derivative_terms: a real power is taken at a positive value'
    end if

    f(0) = x**p
    do k = 1, a % directions
       f(k) = f(k - 1) * (p - real(k - 1, dp)) / x
    end do
    g = composed(a, f)

  end function terms_to_real

  pure function terms_sqrt(a) result(g)

    type(derivative_terms), intent(in) :: a
    type(derivative_terms) :: g

    g = terms_to_real(a, 0.5_dp)

  end function terms_sqrt

end module util_derivative_terms
