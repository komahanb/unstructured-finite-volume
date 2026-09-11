! The Euclidean norm of a real vector, with scaling shared by residual
! reductions and numerical directional differences.
module util_norm

  use util_precision, only : dp
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_positive_inf

  implicit none
  private
  public :: euclidean_norm
  integer, parameter :: bits_kind = selected_int_kind(precision(1.0_dp) + 1)

contains

  pure real(dp) function euclidean_norm(u) result(length)

    real(dp), intent(in) :: u(:)
    real(dp) :: largest, component, squares, magnitude
    integer :: i, common_exponent, difference, result_exponent

    ! Remove the common binary exponent before squaring. The intrinsic
    ! norm2 may still form unscaled squares in an optimized build.
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
    squares = 0.0_dp
    do i = 1, size(u)
       if (u(i) == 0.0_dp) cycle
       difference = exponent(u(i)) - common_exponent
       ! These squares cannot affect the sum at working precision;
       ! omitting them also prevents underflow for mixed magnitudes.
       if (difference < -digits(largest)) cycle
       component = scale(fraction(u(i)), difference)
       squares = squares + component * component
    end do
    magnitude = sqrt(squares)
    result_exponent = common_exponent + exponent(magnitude)
    if (result_exponent > maxexponent(length)) then
       length = ieee_value(0.0_dp, ieee_positive_inf)
    else if (result_exponent == maxexponent(length) .and. &
         & fraction(magnitude) > fraction(huge(length))) then
       length = ieee_value(0.0_dp, ieee_positive_inf)
    else if (result_exponent < minexponent(length)) then
       ! In the supported IEEE binary64 and binary128 formats, a positive
       ! subnormal's bit pattern is its integer count of minimum-subnormal
       ! units. Round that count and construct its bits: even an exact
       ! floating-point product can trap when its result is subnormal.
       difference = minexponent(length) - digits(length)
       component = anint(scale(fraction(magnitude), result_exponent - difference))
       length = transfer(int(component, kind=bits_kind), length)
    else
       length = scale(fraction(magnitude), result_exponent)
    end if

  end function euclidean_norm

end module util_norm
