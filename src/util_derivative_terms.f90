!=====================================================================!
! Module that states a number together with every mixed derivative of
! it along n given directions.
!
! A quantity built from the four arithmetic operations alone stores its
! derivatives beside its value if each operation is given the rule for
! them. That is what this type stores: one coefficient per subset of
! the n directions, 2^n in all, the empty subset storing the value and
! the full subset storing the n-th mixed partial. Nothing is perturbed
! and nothing is truncated: each coefficient is the exact derivative.
!
! The subset index. A subset is a bit mask, and the mask m is stored at
! index m + 1. Allocatable assignment does not preserve a zero lower
! bound in the compilers this is built with, so the masks are stored
! from one and the two are never confused.
!
!      n = 2     mask m    subset     stored at    coefficient
!                00        { }        1            a
!                01        {1}        2            d a / d s_1
!                10        {2}        3            d a / d s_2
!                11        {1, 2}     4            d^2 a / d s_1 d s_2
!
! Symmetric seeding. Set every subset of the same size to the same
! number and the type computes a Taylor composition rather than a mixed
! partial: seed the subsets of size k with the k-th derivative of a
! quantity along one parameter, and the coefficient of the full subset
! comes out as the n-th derivative of whatever was built from it. That
! is the product rule read on subsets: splitting a set of size k into
! two parts counts each split once, which is the binomial coefficient
! the Leibniz rule requires. So nothing here changes; only the
! interpretation of the numbers does.
!
!      c_{12} = a b_{12} + a_1 b_2 + a_2 b_1 + a_{12} b             subsets
!      c''    = a b''    + 2 a' b'           + a'' b               symmetric
!
! The convention is derivatives, not derivatives over factorials: a
! subset of size k stores the k-th derivative itself.
!
! Subset seeding. Set the coefficient of one subset to one number. A
! quantity whose total derivative along every subset of n design
! directions is known (the state of a system solved for the designs,
! whose tangents are those derivatives) is seeded subset by subset, and
! the coefficient of a subset of anything built from the quantity is
! then the total derivative along that subset: the sum over the set
! partitions of the subset that the chain rule requires is what the
! product rule on subsets computes, so no partition is listed.
! Symmetric seeding is the case of one direction repeated.
!
! The four rules:
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
! Composition with a function. For g = f(a) the coefficient of a
! subset m is the sum over the set partitions of m into blocks
! B_1 .. B_p of
!
!      f^(p)(a_0)  a_{B_1} ... a_{B_p}
!
!      g_{12} = f'(a) a_{12} + f''(a) a_1 a_2          the partitions
!                                                      {12} and {1}{2}
!
! computed without listing partitions. With i the lowest element of m
! and G(k, m) the same sum with f^(p+k) in place of f^(p), the block
! containing i is S + i for some S within m - i and the rest is
! partitioned with one more derivative taken:
!
!      G(k, empty)  =  f^(k)(a_0)
!      G(k, m)      =  sum over S within m - i  of  a_{S+i} G(k+1, m - i - S)
!
! filled in increasing m, for k up to n - |m|; the coefficient of g on
! m is G(0, m). Each entry is a subset sum, so the cost is the 3^n
! operations of a product. A function supplies its derivatives at the
! value of a for k = 0 .. n: exp is its own; sin and cos cycle through
! four; log and a real power recur by one division of the value each,
! so no factorial is formed. sqrt is the power one half, and an integer
! power stays a repeated product, which a zero value admits.
!
! A product costs up to 3^n operations and a quotient the same, which
! is what a mixed partial of degree n costs; the storage is 2^n numbers
! per quantity. The degree is bounded by the width of mask a default
! integer can store, which max_subset_width reports.
!
! The rules on the coefficients. Each rule is also stated on arrays of
! the 2^n coefficients, indexed by mask from zero (product_terms,
! quotient_terms and the others), so that an evaluation over many
! quantities stores their coefficients in one array and allocates none
! per operation. The operators of the type are these rules on its own
! coefficients, so the two forms perform the same arithmetic in the
! same order.
!
! What is rejected: a quotient whose divisor has a zero value, and a
! direction index outside one to n; log, sqrt and a real power at a
! value that is not positive, and a derivative table shorter than
! n + 1. Each stops the program: there is no derivative to report.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module util_derivative_terms

  use util_precision, only : dp

  implicit none

  private

  public :: derivative_terms, value, mixed_partial, coefficient, max_subset_width, symmetric_terms
  public :: extend_directions, partial
  public :: integer_power, composed, leibniz_parts, inner_product
  public :: extended_terms, product_terms, quotient_terms, integer_power_terms, real_power_terms
  public :: sine_terms, cosine_terms, exponential_terms, logarithm_terms
  public :: operator(+), operator(-), operator(*), operator(/), operator(**)
  public :: sin, cos, exp, log, sqrt

  !-------------------------------------------------------------------!
  ! Type that models a quantity and its mixed derivatives: the count of
  ! directions, and one coefficient per subset of them, terms(m + 1)
  ! the coefficient of the mask m
  !-------------------------------------------------------------------!

  type :: derivative_terms

     integer               , private              :: directions = 0
     real(dp)              , private, allocatable :: terms(:)

   contains

     procedure :: set_direction
     procedure :: set_symmetric
     procedure :: set_coefficient
     procedure :: num_directions

  end type derivative_terms

  interface derivative_terms
     module procedure create_constant
     module procedure create_like
     module procedure create_from_coefficients
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
  ! Returns the widest mask a default integer indexes: the bit width
  ! less the sign bit and the bit that 2 raised to the width would
  ! occupy, 30 for a 32-bit integer
  !===================================================================!

  pure integer function max_subset_width()

    max_subset_width = bit_size(0) - 2

  end function max_subset_width

  !===================================================================!
  ! Creates a constant: its value, and zero in every derivative. A
  ! negative count of directions, or one past the mask width, stops the
  ! program.
  !===================================================================!

  pure function create_constant(x, num_directions) result(this)

    real(dp), intent(in)   :: x
    integer , intent(in)   :: num_directions
    type(derivative_terms) :: this

    if (num_directions < 0 .or. num_directions > max_subset_width()) then
       error stop 'util_derivative_terms: num_directions is negative or exceeds the mask &
            &width max_subset_width()'
    end if

    this % directions = num_directions
    allocate(this % terms(2**num_directions), source=0.0_dp)
    this % terms(1) = x

  end function create_constant

  !===================================================================!
  ! Creates a constant over the directions of another quantity
  !===================================================================!

  pure function create_like(x, other) result(this)

    real(dp)              , intent(in) :: x
    type(derivative_terms), intent(in) :: other
    type(derivative_terms)             :: this

    this = create_constant(x, other % directions)

  end function create_like

  !===================================================================!
  ! Creates a quantity from its coefficients, one per subset, indexed
  ! by the mask from zero. A count that is not a power of two within
  ! the mask width names no set of directions and stops the program.
  !===================================================================!

  pure function create_from_coefficients(c) result(this)

    real(dp), intent(in)   :: c(0:)
    type(derivative_terms) :: this

    integer :: n

    n = trailz(max(size(c), 1))
    if (size(c) < 1 .or. 2**n /= size(c) .or. n > max_subset_width()) then
       error stop 'util_derivative_terms: size(c) is empty, not a power of two, or exceeds &
            &the mask width - one coefficient is required per subset'
    end if

    this % directions = n
    allocate(this % terms(2**n))
    this % terms(:) = c

  end function create_from_coefficients

  !===================================================================!
  ! Returns the same quantity over more directions, the added ones
  ! being the highest: a subset of the old directions has the same
  ! mask, so the coefficients are copied in place and every subset that
  ! contains an added direction stores zero. A count below the present
  ! one stops the program: a direction is not discarded.
  !===================================================================!

  pure function extend_directions(x, num_directions) result(this)

    type(derivative_terms), intent(in) :: x
    integer               , intent(in) :: num_directions
    type(derivative_terms)             :: this

    if (num_directions < x % directions) then
       error stop 'util_derivative_terms: extend_directions was called with num_directions &
            &less than x % directions - a direction is never discarded'
    end if
    this = create_constant(0.0_dp, num_directions)
    this % terms(1:size(x % terms)) = x % terms

  end function extend_directions

  !===================================================================!
  ! Returns the same, into the coefficients c over the directions of c
  !===================================================================!

  pure subroutine extended_terms(x, c)

    type(derivative_terms), intent(in)  :: x
    real(dp)              , intent(out) :: c(0:)

    if (size(c) < size(x % terms)) then
       error stop 'util_derivative_terms: extended_terms was given fewer coefficients than x stores - a &
            &direction is never discarded'
    end if
    c = 0.0_dp
    c(0:size(x % terms) - 1) = x % terms

  end subroutine extended_terms

  !===================================================================!
  ! Returns the derivative along direction i, as a quantity over the
  ! other directions: the coefficient of a subset m of those is the
  ! coefficient of m with i added. A direction outside those declared
  ! stops the program.
  !
  !      n = 2, i = 1      new mask    0      1
  !                        old mask    01     11      i's bit inserted
  !===================================================================!

  pure function partial(x, i) result(this)

    type(derivative_terms), intent(in) :: x
    integer               , intent(in) :: i
    type(derivative_terms)             :: this

    integer :: m, low, high, bit

    if (i < 1 .or. i > x % directions) then
       error stop 'util_derivative_terms: i is outside the declared directions 1..x % directions'
    end if
    this = create_constant(0.0_dp, x % directions - 1)
    bit = 2**(i - 1)
    do m = 0, size(this % terms) - 1
       ! The mask over the remaining directions, with i's bit inserted
       low  = iand(m, bit - 1)
       high = m - low
       this % terms(m + 1) = x % terms(low + 2 * high + bit + 1)
    end do

  end function partial

  !===================================================================!
  ! Returns the number of directions
  !===================================================================!

  pure integer function num_directions(this)

    class(derivative_terms), intent(in) :: this

    num_directions = this % directions

  end function num_directions

  !===================================================================!
  ! Sets the derivative of this quantity along direction i, which is
  ! the coefficient of the subset holding i alone. An index outside one
  ! to n stops the program.
  !===================================================================!

  pure subroutine set_direction(this, i, x)

    class(derivative_terms), intent(inout) :: this
    integer                , intent(in)    :: i
    real(dp)               , intent(in)    :: x

    if (i < 1 .or. i > this % directions) then
       error stop 'util_derivative_terms: i is outside the declared directions 1..this % directions'
    end if

    this % terms(2**(i - 1) + 1) = x

  end subroutine set_direction

  !===================================================================!
  ! Seeds every subset of one size with one number. An order below zero
  ! or past the directions declared stops the program: there is no such
  ! coefficient to set.
  !===================================================================!

  pure subroutine set_symmetric(this, order, x)

    class(derivative_terms), intent(inout) :: this
    integer                , intent(in)    :: order
    real(dp)               , intent(in)    :: x

    integer :: m

    if (order < 0 .or. order > this % directions) then
       error stop 'util_derivative_terms: order is negative or exceeds this % directions'
    end if

    do m = 0, size(this % terms) - 1
       if (popcnt(m) == order) this % terms(m + 1) = x
    end do

  end subroutine set_symmetric

  !===================================================================!
  ! Seeds the coefficient of the subset with mask m: the total
  ! derivative along the directions m contains. A mask outside 1 to
  ! 2^n - 1 stops the program: the value is not a derivative, and no
  ! other subset exists.
  !===================================================================!

  pure subroutine set_coefficient(this, m, x)

    class(derivative_terms), intent(inout) :: this
    integer                , intent(in)    :: m
    real(dp)               , intent(in)    :: x

    if (m < 1 .or. m >= size(this % terms)) then
       error stop 'util_derivative_terms: m is zero or outside 1..size(this % terms)-1 - the &
            &mask must name a nonempty subset of the directions'
    end if

    this % terms(m + 1) = x

  end subroutine set_coefficient

  !===================================================================!
  ! Returns the quantity whose k-th derivative along one parameter is
  ! d(k), k = 1 .. size(d), seeded on every subset of the first size(d)
  ! directions of the width given, symmetric seeding over the
  ! parameter; every subset containing a later direction stores zero,
  ! so that a later direction can be seeded on one quantity alone as a
  ! perturbation whose mixed coefficients with the parameter's subsets
  ! are the derivatives along the parameter of the partial derivative in
  ! that quantity. More derivatives than directions stops the program.
  !===================================================================!

  pure function symmetric_terms(x, d, width) result(this)

    real(dp), intent(in)   :: x, d(:)
    integer , intent(in)   :: width
    type(derivative_terms) :: this

    integer :: m

    if (size(d) > width) then
       error stop 'util_derivative_terms: symmetric_terms was given more derivatives than directions'
    end if
    this = create_constant(x, width)
    do m = 1, 2**size(d) - 1
       this % terms(m + 1) = d(popcnt(m))
    end do

  end function symmetric_terms

  !===================================================================!
  ! Returns the value, the coefficient of the empty subset
  !===================================================================!

  pure real(dp) function value(x)

    type(derivative_terms), intent(in) :: x

    value = x % terms(1)

  end function value

  !===================================================================!
  ! Returns the mixed partial along every direction, the coefficient of
  ! the full subset
  !===================================================================!

  pure real(dp) function mixed_partial(x)

    type(derivative_terms), intent(in) :: x

    mixed_partial = x % terms(size(x % terms))

  end function mixed_partial

  !===================================================================!
  ! Returns the coefficient of the subset with mask m: the mixed
  ! partial along the directions m contains. A mask outside 0 to
  ! 2^n - 1 stops the program.
  !===================================================================!

  pure real(dp) function coefficient(x, m)

    type(derivative_terms), intent(in) :: x
    integer               , intent(in) :: m

    if (m < 0 .or. m >= size(x % terms)) then
       error stop 'util_derivative_terms: m is outside 0..size(x % terms)-1'
    end if

    coefficient = x % terms(m + 1)

  end function coefficient

  !===================================================================!
  ! Returns the sum, coefficient by coefficient, the first of the four
  ! rules
  !===================================================================!

  pure function terms_plus(a, b) result(c)

    type(derivative_terms), intent(in) :: a, b
    type(derivative_terms)             :: c

    c % directions = a % directions
    c % terms = a % terms + b % terms

  end function terms_plus

  !===================================================================!
  ! Returns the difference, coefficient by coefficient
  !===================================================================!

  pure function terms_minus(a, b) result(c)

    type(derivative_terms), intent(in) :: a, b
    type(derivative_terms)             :: c

    c % directions = a % directions
    c % terms = a % terms - b % terms

  end function terms_minus

  !===================================================================!
  ! Returns the negation, coefficient by coefficient
  !===================================================================!

  pure function terms_negated(a) result(c)

    type(derivative_terms), intent(in) :: a
    type(derivative_terms)             :: c

    c % directions = a % directions
    c % terms = -a % terms

  end function terms_negated

  !===================================================================!
  ! Returns a real number times a quantity, coefficient by coefficient
  !===================================================================!

  pure function real_times_terms(x, a) result(c)

    real(dp)              , intent(in) :: x
    type(derivative_terms), intent(in) :: a
    type(derivative_terms)             :: c

    c % directions = a % directions
    c % terms = x * a % terms

  end function real_times_terms

  !===================================================================!
  ! Returns the product of two quantities
  !===================================================================!

  pure function terms_times(a, b) result(c)

    type(derivative_terms), intent(in) :: a, b
    type(derivative_terms)             :: c

    c = create_like(0.0_dp, a)
    call product_terms(a % terms, b % terms, c % terms)

  end function terms_times

  !===================================================================!
  ! Returns the quotient of two quantities
  !===================================================================!

  pure function terms_over(a, b) result(r)

    type(derivative_terms), intent(in) :: a, b
    type(derivative_terms)             :: r

    r = create_like(0.0_dp, a)
    call quotient_terms(a % terms, b % terms, r % terms)

  end function terms_over

  !===================================================================!
  ! Returns the product on the coefficients: c(m) the sum over the
  ! splittings of m of a on one part times b on the other. The parts are
  ! the submasks s of m, visited as s = m, (s - 1) and m, ..., 0:
  !
  !      m = 11      s = 11, 10, 01, 00
  !      c_11 = a_11 b_00 + a_10 b_01 + a_01 b_10 + a_00 b_11
  !===================================================================!

  pure subroutine product_terms(a, b, c)

    real(dp), intent(in)  :: a(0:), b(0:)
    real(dp), intent(out) :: c(0:)

    integer :: m, s

    do m = 0, size(a) - 1
       c(m) = 0.0_dp
       s = m
       do
          c(m) = c(m) + a(s) * b(ieor(m, s))
          if (s == 0) exit
          s = iand(s - 1, m)
       end do
    end do

  end subroutine product_terms

  !===================================================================!
  ! Returns the quotient on the coefficients: r b = a solved subset by
  ! subset in increasing order,
  !
  !      r_m = (a_m - sum over nonempty s within m of b_s r_(m - s)) / b_0
  !
  ! every r on the right of a smaller mask, so already solved
  !===================================================================!

  pure subroutine quotient_terms(a, b, r)

    real(dp), intent(in)  :: a(0:), b(0:)
    real(dp), intent(out) :: r(0:)

    real(dp) :: accumulated
    integer  :: m, s

    if (b(0) == 0.0_dp) then
       error stop 'util_derivative_terms: b % terms(1) is zero - a quotient cannot divide by zero'
    end if

    do m = 0, size(a) - 1
       accumulated = a(m)
       s = m
       do while (s /= 0)
          accumulated = accumulated - b(s) * r(ieor(m, s))
          s = iand(s - 1, m)
       end do
       r(m) = accumulated / b(0)
    end do

  end subroutine quotient_terms

  !===================================================================!
  ! Returns the product rule on the full subset, read by the size of
  ! one factor's part: parts(k) is the sum over the subsets S of size k
  ! of a_S b_(full - S), so the parts add to the n-th mixed partial of
  ! a b. Under symmetric seeding parts(k) is C(n, k) a^(k) b^(n-k), the
  ! k-th term of Leibniz. Factors over different counts of directions
  ! stop the program.
  !===================================================================!

  pure function leibniz_parts(a, b) result(parts)

    type(derivative_terms), intent(in) :: a, b
    real(dp)                           :: parts(0:a % directions)

    integer :: full, s

    if (a % directions /= b % directions) then
       error stop 'util_derivative_terms: a % directions and b % directions differ - the &
            &factors of a product must share the same count of directions'
    end if

    full  = size(a % terms) - 1
    parts = 0.0_dp
    do s = 0, full
       parts(popcnt(s)) = parts(popcnt(s)) + a % terms(s + 1) * b % terms(ieor(full, s) + 1)
    end do

  end function leibniz_parts

  !===================================================================!
  ! Returns the inner product sum_i a(i) b(i), the entries an active
  ! mask excludes taken as zero: the same reduction real arrays take
  ! under dot_product, over the terms this module's own product rule
  ! multiplies. An active mask shorter than the factors stops the
  ! program.
  !===================================================================!

  pure function inner_product(a, b, active) result(c)

    type(derivative_terms), intent(in)           :: a(:), b(:)
    logical               , intent(in), optional :: active(:)
    type(derivative_terms)                       :: c

    integer :: i

    if (size(a) /= size(b)) then
       error stop 'util_derivative_terms: size(a) and size(b) differ - an inner product pairs &
            &factors of the same extent'
    end if
    if (present(active)) then
       if (size(active) < size(a)) then
          error stop 'util_derivative_terms: size(active) is smaller than size(a) - the active &
               &mask must cover every entry'
       end if
    end if

    c = create_like(0.0_dp, a(1))
    do i = 1, size(a)
       if (present(active)) then
          if (.not. active(i)) cycle
       end if
       c = c + a(i) * b(i)
    end do

  end function inner_product

  !===================================================================!
  ! Returns x raised to an integer exponent, of either sign. A zero
  ! exponent is one, whatever the value; a negative exponent divides,
  ! and a zero value then stops the program inside the quotient.
  !===================================================================!

  pure function integer_power(x, n) result(p)

    type(derivative_terms), intent(in) :: x
    integer               , intent(in) :: n
    type(derivative_terms)             :: p

    p = create_like(0.0_dp, x)
    call integer_power_terms(x % terms, n, p % terms)

  end function integer_power

  !===================================================================!
  ! Returns the integer power on the coefficients: |n| products from
  ! one, then the reciprocal for a negative n
  !===================================================================!

  pure subroutine integer_power_terms(x, n, p)

    real(dp), intent(in)  :: x(0:)
    integer , intent(in)  :: n
    real(dp), intent(out) :: p(0:)

    real(dp) :: previous(0:size(x) - 1), unit(0:size(x) - 1)
    integer  :: i

    unit    = 0.0_dp
    unit(0) = 1.0_dp
    p = unit
    do i = 1, abs(n)
       previous = p
       call product_terms(previous, x, p)
    end do
    if (n < 0) then
       previous = p
       call quotient_terms(unit, previous, p)
    end if

  end subroutine integer_power_terms

  !===================================================================!
  ! Returns x ** n for an integer n
  !===================================================================!

  pure function terms_to_integer(x, n) result(p)

    type(derivative_terms), intent(in) :: x
    integer               , intent(in) :: n
    type(derivative_terms)             :: p

    p = integer_power(x, n)

  end function terms_to_integer

  !===================================================================!
  ! Returns the composition with a function: g = f(a) from the
  ! derivatives f(0:n) of f at the value of a. A table shorter than
  ! n + 1 stops the program: the full subset has no derivative to read.
  !===================================================================!

  pure function composed(a, f) result(g)

    type(derivative_terms), intent(in) :: a
    real(dp)              , intent(in) :: f(0:)
    type(derivative_terms)             :: g

    g = create_like(0.0_dp, a)
    call composed_terms(a % terms, f, g % terms)

  end function composed

  !===================================================================!
  ! Returns the composition on the coefficients, by the table G(k, m)
  ! of the module header: G(k, 0) = f^(k), each later column the sum
  ! over the submasks of m less its lowest element i
  !===================================================================!

  pure subroutine composed_terms(a, f, g)

    real(dp), intent(in)  :: a(0:), f(0:)
    real(dp), intent(out) :: g(0:)

    real(dp), allocatable :: table(:,:)
    real(dp)              :: accumulated
    integer               :: n, m, k, i, remaining_subset, s

    n = trailz(size(a))

    if (ubound(f, 1) < n) then
       error stop 'util_derivative_terms: ubound(f,1) is below n - the composed function must &
            &supply n + 1 derivatives'
    end if

    allocate(table(0:n, 0:2**n - 1))
    table(:, 0) = f(0:n)

    do m = 1, 2**n - 1
       i    = ibset(0, trailz(m))
       remaining_subset = ieor(m, i)
       do k = 0, n - popcnt(m)
          accumulated = 0.0_dp
          s = remaining_subset
          do
             accumulated = accumulated + a(ior(s, i)) * table(k + 1, ieor(remaining_subset, s))
             if (s == 0) exit
             s = iand(s - 1, remaining_subset)
          end do
          table(k, m) = accumulated
       end do
    end do

    do m = 0, 2**n - 1
       g(m) = table(0, m)
    end do

  end subroutine composed_terms

  !===================================================================!
  ! Returns exp(a). Each function below supplies its derivative table,
  ! f^(k) at the value for k = 0 .. n; exp's is its own value
  !===================================================================!

  pure function terms_exp(a) result(g)

    type(derivative_terms), intent(in) :: a
    type(derivative_terms)             :: g

    g = create_like(0.0_dp, a)
    call exponential_terms(a % terms, g % terms)

  end function terms_exp

  !===================================================================!
  ! Returns exp on the coefficients
  !===================================================================!

  pure subroutine exponential_terms(a, g)

    real(dp), intent(in)  :: a(0:)
    real(dp), intent(out) :: g(0:)

    real(dp) :: f(0:trailz(size(a)))

    f = exp(a(0))
    call composed_terms(a, f, g)

  end subroutine exponential_terms

  !===================================================================!
  ! Returns sin(a)
  !===================================================================!

  pure function terms_sin(a) result(g)

    type(derivative_terms), intent(in) :: a
    type(derivative_terms)             :: g

    g = create_like(0.0_dp, a)
    call sine_terms(a % terms, g % terms)

  end function terms_sin

  !===================================================================!
  ! Returns sin on the coefficients, its derivatives the cycle
  ! sin, cos, -sin, -cos
  !===================================================================!

  pure subroutine sine_terms(a, g)

    real(dp), intent(in)  :: a(0:)
    real(dp), intent(out) :: g(0:)

    real(dp) :: f(0:trailz(size(a))), cycle(0:3)
    integer  :: k

    cycle = [sin(a(0)), cos(a(0)), -sin(a(0)), -cos(a(0))]
    do k = 0, trailz(size(a))
       f(k) = cycle(mod(k, 4))
    end do
    call composed_terms(a, f, g)

  end subroutine sine_terms

  !===================================================================!
  ! Returns cos(a)
  !===================================================================!

  pure function terms_cos(a) result(g)

    type(derivative_terms), intent(in) :: a
    type(derivative_terms)             :: g

    g = create_like(0.0_dp, a)
    call cosine_terms(a % terms, g % terms)

  end function terms_cos

  !===================================================================!
  ! Returns cos on the coefficients, its derivatives the cycle
  ! cos, -sin, -cos, sin
  !===================================================================!

  pure subroutine cosine_terms(a, g)

    real(dp), intent(in)  :: a(0:)
    real(dp), intent(out) :: g(0:)

    real(dp) :: f(0:trailz(size(a))), cycle(0:3)
    integer  :: k

    cycle = [cos(a(0)), -sin(a(0)), -cos(a(0)), sin(a(0))]
    do k = 0, trailz(size(a))
       f(k) = cycle(mod(k, 4))
    end do
    call composed_terms(a, f, g)

  end subroutine cosine_terms

  !===================================================================!
  ! Returns log(a): f^(k) = (-1)^(k-1) (k-1)! / x^k, formed as
  ! f^(k) = -f^(k-1) (k-1) / x from f^(1) = 1 / x. A value that is not
  ! positive stops the program.
  !===================================================================!

  pure function terms_log(a) result(g)

    type(derivative_terms), intent(in) :: a
    type(derivative_terms)             :: g

    g = create_like(0.0_dp, a)
    call logarithm_terms(a % terms, g % terms)

  end function terms_log

  !===================================================================!
  ! Returns log on the coefficients
  !===================================================================!

  pure subroutine logarithm_terms(a, g)

    real(dp), intent(in)  :: a(0:)
    real(dp), intent(out) :: g(0:)

    real(dp) :: f(0:trailz(size(a))), x
    integer  :: k, n

    n = trailz(size(a))
    x = a(0)
    if (x <= 0.0_dp) then
       error stop 'util_derivative_terms: a % terms(1) is not positive - log is only &
            &differentiated at a positive value'
    end if

    f(0) = log(x)
    if (n >= 1) f(1) = 1.0_dp / x
    do k = 2, n
       f(k) = -f(k - 1) * real(k - 1, dp) / x
    end do
    call composed_terms(a, f, g)

  end subroutine logarithm_terms

  !===================================================================!
  ! Returns a real power: f^(k) = p (p-1) .. (p-k+1) x^(p-k), formed as
  ! f^(k) = f^(k-1) (p-k+1) / x from f^(0) = x^p. A value that is not
  ! positive stops the program; an integer exponent is a repeated
  ! product and admits any value.
  !===================================================================!

  pure function terms_to_real(a, p) result(g)

    type(derivative_terms), intent(in) :: a
    real(dp)              , intent(in) :: p
    type(derivative_terms)             :: g

    g = create_like(0.0_dp, a)
    call real_power_terms(a % terms, p, g % terms)

  end function terms_to_real

  !===================================================================!
  ! Returns the real power on the coefficients
  !===================================================================!

  pure subroutine real_power_terms(a, p, g)

    real(dp), intent(in)  :: a(0:), p
    real(dp), intent(out) :: g(0:)

    real(dp) :: f(0:trailz(size(a))), x
    integer  :: k

    x = a(0)
    if (x <= 0.0_dp) then
       error stop 'util_derivative_terms: a % terms(1) is not positive - a real power is only &
            &differentiated at a positive value'
    end if

    f(0) = x**p
    do k = 1, trailz(size(a))
       f(k) = f(k - 1) * (p - real(k - 1, dp)) / x
    end do
    call composed_terms(a, f, g)

  end subroutine real_power_terms

  !===================================================================!
  ! Returns sqrt(a), the power one half
  !===================================================================!

  pure function terms_sqrt(a) result(g)

    type(derivative_terms), intent(in) :: a
    type(derivative_terms)             :: g

    g = terms_to_real(a, 0.5_dp)

  end function terms_sqrt

end module util_derivative_terms
