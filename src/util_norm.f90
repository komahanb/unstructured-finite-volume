!=====================================================================!
! Module that computes the Euclidean norm of a real vector without
! overflow or underflow in its squares. The residual reductions and
! the numerical directional differences share it.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module util_norm

  use util_precision, only : dp
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_positive_inf

  implicit none

  private

  public :: euclidean_norm

  ! An integer kind as wide as real(dp), to hold the bits of a subnormal
  integer, parameter :: bits_kind = selected_int_kind(precision(1.0_dp) + 1)

contains

  !===================================================================!
  ! Returns the Euclidean norm of u, sqrt(u(1)**2 + u(2)**2 + ...).
  !
  ! Each entry is a fraction times a power of two, u(i) = f(i) 2**e(i)
  ! with 1/2 <= |f(i)| < 1. The exponent e of the largest entry is
  ! taken out before squaring, so that every scaled entry lies in
  ! (-1, 1) and no square overflows:
  !
  !   exponent        e - digits                 e
  !                       |                      |
  !   u(4)        x       |                      |   skipped
  !   u(3)        ........|..........x           |   2**(e(3) - e) f(3)
  !   u(1)        ........|......................x   the largest
  !
  !   length = sqrt(sum of the scaled squares) * 2**e
  !
  ! An entry below e - digits is skipped: its square cannot change the
  ! sum at working precision. The intrinsic norm2 may still form
  ! unscaled squares in an optimized build.
  !===================================================================!

  pure real(dp) function euclidean_norm(u) result(length)

    real(dp), intent(in) :: u(:)
    real(dp)             :: largest, component, squares, magnitude
    integer              :: i, common_exponent, difference, result_exponent

    ! The largest entry; an entry that is not finite is the norm
    largest = 0.0_dp
    do i = 1, size(u)
       if (.not. ieee_is_finite(u(i))) then
          length = abs(u(i))
          return
       end if
       largest = max(largest, abs(u(i)))
    end do
    length = largest
    if (largest == 0.0_dp) return
    common_exponent = exponent(largest)

    ! The sum of the scaled squares
    squares = 0.0_dp
    do i = 1, size(u)
       if (u(i) == 0.0_dp) cycle
       difference = exponent(u(i)) - common_exponent
       ! Omitting these squares also prevents underflow for mixed
       ! magnitudes
       if (difference < -digits(largest)) cycle
       component = scale(fraction(u(i)), difference)
       squares = squares + component * component
    end do

    ! The scale put back. With sqrt(squares) = g 2**k, 1/2 <= g < 1,
    ! the norm is g 2**(e + k):
    !
    !   e + k above maxexponent, or at it with g > fraction(huge)   +inf
    !   e + k below minexponent                                    a subnormal
    !   otherwise                                                  scale(g, e + k)
    magnitude = sqrt(squares)
    result_exponent = common_exponent + exponent(magnitude)
    if (result_exponent > maxexponent(length)) then
       length = ieee_value(0.0_dp, ieee_positive_inf)
    else if (result_exponent == maxexponent(length) .and. &
         & fraction(magnitude) > fraction(huge(length))) then
       length = ieee_value(0.0_dp, ieee_positive_inf)
    else if (result_exponent < minexponent(length)) then
       ! A positive subnormal is a whole number n of units of the least
       ! subnormal, 2**(minexponent - digits). In the supported IEEE
       ! binary64 and binary128 formats its bit pattern is the integer
       ! n:
       !
       !   0       1        2        3                 2**minexponent
       !   |-------|--------|--------|----- ... -------|
       !    1 unit                                       the least normal
       !
       ! n is rounded and its bits are built directly: even an exact
       ! floating-point product can trap when its result is subnormal.
       difference = minexponent(length) - digits(length)
       component = anint(scale(fraction(magnitude), result_exponent - difference))
       length = transfer(int(component, kind=bits_kind), length)
    else
       length = scale(fraction(magnitude), result_exponent)
    end if

  end function euclidean_norm

end module util_norm
