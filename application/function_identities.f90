! The elementary functions over derivative_terms, checked by
! identities that hold coefficient by coefficient when the composition
! rule is right, and by one closed form.
!
! The quantity a is nonlinear in the directions - a power of a seeded
! linear term - so every subset of every size carries a nonzero
! coefficient and the set-partition sum is exercised in full. The
! identities:
!
!      exp(a + b)  =  exp(a) exp(b)
!      sin(a)^2 + cos(a)^2  =  1
!      sin(a + b)  =  sin a cos b + cos a sin b
!      log(exp(a))  =  a          exp(log(a))  =  a
!      sqrt(a) sqrt(a)  =  a      a^(1/2)  =  sqrt(a)      a^3 (real)  =  a a a
!
! and for a linear in the directions the coefficient of the full
! subset is f^(n)(a_0) times the product of the direction values.
!
! The floor each difference is held to is the arithmetic's: 3^n
! operations per coefficient, each within one epsilon of the largest
! coefficient that entered them - for an identity that cancels, the
! operand cancelled, whose rounding is proportional to its own size
! and not to the result's. Nothing here is a chosen number.
program function_identities

  use util_precision       , only : dp
  use util_derivative_terms, only : derivative_terms, integer_power, mixed_partial, coefficient, &
       & operator(+), operator(-), operator(*), operator(/), operator(**), &
       & sin, cos, exp, log, sqrt

  implicit none

  integer :: n, failures

  failures = 0

  do n = 1, 5
     call identities(n, failures)
  end do

  if (failures > 0) then
     write(*,'(a,i0,a)') ' FAIL : ', failures, ' identities exceeded the floor'
     error stop
  end if
  write(*,'(a)') ' PASS : the elementary functions compose exactly to five directions'

contains

  subroutine identities(n, failures)

    integer, intent(in)    :: n
    integer, intent(inout) :: failures

    type(derivative_terms) :: a, b, one
    real(dp) :: v(n), w(n), closed
    integer  :: i

    do i = 1, n
       v(i) = 0.3_dp + 0.1_dp * real(i, dp)
       w(i) = 0.7_dp - 0.1_dp * real(i, dp)
    end do

    ! nonlinear positive quantities: every subset nonzero
    a = seeded(1.2_dp, v)
    a = integer_power(a, n) + seeded(0.5_dp, w)
    b = seeded(0.8_dp, w)
    b = integer_power(b, n)
    one = derivative_terms(1.0_dp, n)

    call held(n, 'exp(a+b) = exp(a) exp(b)', exp(a + b), exp(a) * exp(b), exp(a + b), failures)
    call held(n, 'sin^2 + cos^2 = 1', sin(a) * sin(a) + cos(a) * cos(a), one, sin(a) * sin(a), failures)
    call held(n, 'sin(a+b) addition', sin(a + b), sin(a) * cos(b) + cos(a) * sin(b), sin(a) * cos(b), failures)
    call held(n, 'log(exp(a)) = a', log(exp(a)), a, exp(a), failures)
    call held(n, 'exp(log(a)) = a', exp(log(a)), a, a, failures)
    call held(n, 'sqrt(a) sqrt(a) = a', sqrt(a) * sqrt(a), a, a, failures)
    call held(n, 'a**0.5 = sqrt(a)', a ** 0.5_dp, sqrt(a), sqrt(a), failures)
    call held(n, 'a**3.0 = a a a', a ** 3.0_dp, a * a * a, a * a * a, failures)
    call held(n, 'a**(-1.0) = 1/a', a ** (-1.0_dp), one / a, one / a, failures)

    ! closed form on a linear quantity: d^n exp / prod dv = exp(x0) prod v
    a = seeded(0.4_dp, v)
    closed = exp(0.4_dp) * product(v)
    call held_scalar(n, 'full partial of exp', mixed_partial(exp(a)), closed, failures)
    closed = sin(0.4_dp + real(n, dp) * acos(-1.0_dp) / 2.0_dp) * product(v)
    call held_scalar(n, 'full partial of sin', mixed_partial(sin(a)), closed, failures)

  end subroutine identities

  function seeded(x, v) result(a)

    real(dp), intent(in) :: x, v(:)
    type(derivative_terms) :: a

    integer :: i

    a = derivative_terms(x, size(v))
    do i = 1, size(v)
       call a % set_direction(i, v(i))
    end do

  end function seeded

  !-------------------------------------------------------------------!
  ! Two quantities agree when every coefficient's difference is within
  ! the floor: 3^n operations, each within one epsilon of the largest
  ! coefficient of the operand that entered them.
  !-------------------------------------------------------------------!

  subroutine held(n, label, got, reference, operand, failures)

    integer               , intent(in)    :: n
    character(len=*)      , intent(in)    :: label
    type(derivative_terms), intent(in)    :: got, reference, operand
    integer               , intent(inout) :: failures

    real(dp) :: worst, scale, floor
    integer  :: m

    worst = 0.0_dp
    scale = 0.0_dp
    do m = 0, 2**n - 1
       worst = max(worst, abs(coefficient(got, m) - coefficient(reference, m)))
       scale = max(scale, abs(coefficient(operand, m)))
    end do
    floor = real(3**n, dp) * epsilon(1.0_dp) * scale

    call reported(n, label, worst, floor, failures)

  end subroutine held

  subroutine held_scalar(n, label, got, reference, failures)

    integer         , intent(in)    :: n
    character(len=*), intent(in)    :: label
    real(dp)        , intent(in)    :: got, reference
    integer         , intent(inout) :: failures

    call reported(n, label, abs(got - reference), &
         & real(3**n, dp) * epsilon(1.0_dp) * abs(reference), failures)

  end subroutine held_scalar

  subroutine reported(n, label, worst, floor, failures)

    integer         , intent(in)    :: n
    character(len=*), intent(in)    :: label
    real(dp)        , intent(in)    :: worst, floor
    integer         , intent(inout) :: failures

    write(*,'(a,i0,a,a28,a,es9.2,a,es9.2)') '   n = ', n, '  ', label, &
         & '  difference ', worst, '  floor ', floor
    if (worst > floor) failures = failures + 1

  end subroutine reported

end program function_identities
