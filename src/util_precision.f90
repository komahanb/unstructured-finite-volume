!=====================================================================!
! The working precision of the whole tower, in one place.
!
! Every module takes its real kind from here and nowhere else, so the
! precision is a property of a build and not of a file. PRECISION=quad
! at build time selects real128; nothing else selects real64. Both
! builds keep their own library and binaries, so they coexist.
!
! The three kinds' spacings at one are given by name, so a caller can
! say which kind a target needs from the floor eps ||A|| ||q|| a
! march can reach, without naming a number.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module util_precision

  use iso_fortran_env, only : real32, real64, real128

  implicit none

  private
  public :: dp, precision_named
  public :: single_spacing, double_spacing, quadruple_spacing, least_kind_for

#ifdef PRECISION_QUAD
  integer, parameter :: dp = real128
#else
  integer, parameter :: dp = real64
#endif

  real(real128), parameter :: single_spacing    = real(epsilon(1.0_real32), real128)
  real(real128), parameter :: double_spacing    = real(epsilon(1.0_real64), real128)
  real(real128), parameter :: quadruple_spacing = epsilon(1.0_real128)

contains

  pure function precision_named() result(name)

    character(len=:), allocatable :: name

    if (dp == real128) then
       name = 'quadruple'
    else
       name = 'double'
    end if

  end function precision_named

  !===================================================================!
  ! The least of the three kinds whose spacing is under the one
  ! required, or none where even the widest is not.
  !===================================================================!

  pure function least_kind_for(required) result(name)

    real(real128), intent(in) :: required
    character(len=:), allocatable :: name

    if (single_spacing < required) then
       name = 'single'
    else if (double_spacing < required) then
       name = 'double'
    else if (quadruple_spacing < required) then
       name = 'quadruple'
    else
       name = 'none'
    end if

  end function least_kind_for

end module util_precision
