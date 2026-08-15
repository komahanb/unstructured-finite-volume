!=====================================================================!
! THE SPECIMEN'S COEFFICIENTS - earned at LEVEL 5, and the first place
! in this tower where a number exists at all.
!
!      w1 : E1 -> reals =  [ 2, -1,  0,  3,  4 ]
!      w2 : E2 -> reals =  [ 1,  5, -2,  2 ]
!      w3 : E3 -> reals =  [ 3, -1,  4 ]
!
!                    WHY THE VALUES LIVE ON E, NOT ON X
!
! A coefficient belongs to a dependency OCCURRENCE, not to the member
! it reads and not to the member it writes. e12 and e13 both read b;
! they carry -1 and 0. A field on X0 could not hold both, and a field
! on X1 could not either. E_k is where the two ends meet, and it is
! the only carrier of the three that can seat the number.
!
! This is why Level 0 gave the occurrences first-class identity five
! levels before anything needed a value: a coefficient needs somewhere
! to live that survives being zero.
!
!                        THE ZERO IS THE EXPERIMENT
!
!      e13 : b -> q        w1(e13) = 0
!
! Structurally b -> q is in D1, and Level 5 exists to prove that the
! zero does not remove it. Only w1 carries a zero; the other two
! fields are ordinary, because a second zero would add nothing and
! would blur which assertion is load-bearing.
!
!                     THE ALTERNATE IS A PROBE, NOT PHYSICS
!
! w1_alt = [9, 8, 7, 6, 5] is a second field on the SAME E1, used once
! to show that D1 and w1 determine nothing about each other. It is not
! part of the specimen and nothing above Level 5 should carry it
! forward.
!
!                       THESE ARE ORACLES
!
! COEFF_A1, COEFF_A2, COEFF_A3 and COEFF_A1_ALT are written here as
! plain parameters and are compared against by the level above. No
! level may obtain them from the machinery it is testing.
!
! They are named COEFF_* rather than W1..W3 because Fortran does not
! distinguish case, and a local field variable called w1 would
! silently shadow a parameter called W1 - the same trap the Level-4
! renderer's stub width fell into.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module visualization_values_fixture

  use iso_fortran_env  , only : dp => REAL64
  use graph_carrier    , only : member_set
  use class_graph_field, only : field

  implicit none

  private
  public :: coefficients_of_a1, coefficients_of_a2, coefficients_of_a3
  public :: alternate_coefficients_of_a1
  public :: COEFF_A1, COEFF_A2, COEFF_A3, COEFF_A1_ALT, TOL

  !-------------------------------------------------------------------!
  ! The twelve coefficients, in each occurrence carrier's own
  ! declaration order.
  !
  !      w1(e11) =  2      e11 : a -> p
  !      w1(e12) = -1      e12 : b -> p
  !      w1(e13) =  0      e13 : b -> q     <-- the witness
  !      w1(e14) =  3      e14 : c -> q
  !      w1(e15) =  4      e15 : d -> r
  !-------------------------------------------------------------------!

  real(dp), parameter :: COEFF_A1(5) = [ 2.0_dp, -1.0_dp,  0.0_dp,  3.0_dp, 4.0_dp]
  real(dp), parameter :: COEFF_A2(4) = [ 1.0_dp,  5.0_dp, -2.0_dp,  2.0_dp]
  real(dp), parameter :: COEFF_A3(3) = [ 3.0_dp, -1.0_dp,  4.0_dp]

  real(dp), parameter :: COEFF_A1_ALT(5) = [9.0_dp, 8.0_dp, 7.0_dp, 6.0_dp, 5.0_dp]

  !-------------------------------------------------------------------!
  ! Every value here is a whole number set once and read back once,
  ! so the comparison could be exact. It is written with a tolerance
  ! anyway, because a level that only ever compares reals exactly
  ! teaches the next level the wrong habit.
  !-------------------------------------------------------------------!

  real(dp), parameter :: TOL = 1.0e-14_dp

contains

  !===================================================================!
  ! One coefficient field per operator, on that operator's own
  ! occurrence carrier. The field takes its entry count from the
  ! domain, so a value vector of the wrong length is refused by the
  ! nucleus rather than by this file.
  !===================================================================!

  type(field) function coefficients_of_a1(e1) result(w)

    class(member_set), intent(in) :: e1

    w = field('w1', e1, ncomp=1)
    call w % set_real_vector(COEFF_A1)

  end function coefficients_of_a1

  type(field) function coefficients_of_a2(e2) result(w)

    class(member_set), intent(in) :: e2

    w = field('w2', e2, ncomp=1)
    call w % set_real_vector(COEFF_A2)

  end function coefficients_of_a2

  type(field) function coefficients_of_a3(e3) result(w)

    class(member_set), intent(in) :: e3

    w = field('w3', e3, ncomp=1)
    call w % set_real_vector(COEFF_A3)

  end function coefficients_of_a3

  !===================================================================!
  ! The independence probe: a second field on the SAME E1, with no
  ! value in common with the first and no zero among them.
  !===================================================================!

  type(field) function alternate_coefficients_of_a1(e1) result(w)

    class(member_set), intent(in) :: e1

    w = field('w1-alt', e1, ncomp=1)
    call w % set_real_vector(COEFF_A1_ALT)

  end function alternate_coefficients_of_a1

end module visualization_values_fixture
