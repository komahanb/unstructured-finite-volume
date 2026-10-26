!=====================================================================!
! Module that states the working precision of the whole tower, in one
! place.
!
! Every module takes its real kind from here and nowhere else, so the
! precision is a property of a build and not of a file. PRECISION=quad
! at build time selects real128; otherwise the kind is real64. Each
! build keeps its own library and binaries, so the two coexist.
!
! The spacings at one of the three kinds are given by name, so that a
! caller can determine which kind a target requires from the lower
! bound eps ||A|| ||q|| a march can reach, without stating a number:
!
!   log10 of the spacing at one
!
!   -34                       -16               -7          0
!    |                         |                 |          |
!    quadruple                 double            single
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module util_precision

  use iso_fortran_env, only : real32, real64, real128

  implicit none

  private

  public :: dp, precision_named
  public :: single_spacing, double_spacing, quadruple_spacing, least_kind_for
  public :: spacing_at_one, half_digits

#ifdef PRECISION_QUAD
  integer, parameter :: dp = real128
#else
  integer, parameter :: dp = real64
#endif

  real(real128), parameter :: single_spacing    = real(epsilon(1.0_real32), real128)
  real(real128), parameter :: double_spacing    = real(epsilon(1.0_real64), real128)
  real(real128), parameter :: quadruple_spacing = epsilon(1.0_real128)

  ! The lower bounds of this build's own kind, from which every
  ! tolerance in the tower is measured. The spacing at one is the size
  ! below which a residual cannot be told from zero. Its square root
  ! is the step at which a first difference loses half its digits to
  ! rounding and half to truncation, and so the reduction past which a
  ! first-order iterate gives a caller no further information:
  !
  !   error of a first difference
  !
  !     |  \                                  /
  !     |    \  rounding, eps / h           /  truncation, h
  !     |      \_                        _/
  !     |         \___              ___/
  !     |             \____________/
  !     +---------------------+------------------------ h
  !                       sqrt(eps)
  real(dp), parameter :: spacing_at_one = epsilon(1.0_dp)
  real(dp), parameter :: half_digits    = sqrt(epsilon(1.0_dp))

contains

  !===================================================================!
  ! Returns the name of this build's precision: 'quadruple' or
  ! 'double'
  !===================================================================!

  pure function precision_named() result(name)

    character(len=:), allocatable :: name

    if (dp == real128) then
       name = 'quadruple'
    else
       name = 'double'
    end if

  end function precision_named

  !===================================================================!
  ! Returns the least of the three kinds whose spacing is below the one
  ! required, or 'none' where even the widest is not:
  !
  !   required   above single   between single   between double   below quadruple
  !                             and double       and quadruple
  !
  !   returns    'single'       'double'         'quadruple'      'none'
  !===================================================================!

  pure function least_kind_for(required) result(name)

    real(real128)   , intent(in)  :: required
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
