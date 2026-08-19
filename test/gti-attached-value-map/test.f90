!=====================================================================!
! The graph-attached value map suite: value_map : graph identity
! -> value status x value buffer, proven over the closed status
! vocabulary
!
!      UNATTACHED -> attach -> UNKNOWN -> mark_known -> KNOWN
!      KNOWN -> mark_unknown -> UNKNOWN        (trust withdrawn)
!      any attached -> detach -> UNATTACHED    (the seat leaves)
!
! with rows keyed on copied identity tokens and never on position:
! detaching an early row re-aims nothing, a copied graph IS its
! identity, and the map outlives the variable that built it.
! Structural knownness lives in fractal_graph. Value knownness
! lives in an attached value map. Reversible changes update maps,
! but do not define the map ontology.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_gti_attached_value_map

  use iso_fortran_env      , only : dp => REAL64
  use fractal_graph        , only : graph
  use gti_value_buffers    , only : gti_value_buffer
  use gti_attached_value_maps, only : gti_attached_value_map, &
       & GTI_VALUE_STATUS_UNATTACHED, GTI_VALUE_STATUS_UNKNOWN, &
       & GTI_VALUE_STATUS_KNOWN

  implicit none

  type(gti_attached_value_map) :: map
  type(graph) :: a, b, c
  type(graph) :: b_copy, keeper

  type(gti_value_buffer) :: values, stored
  real(dp), allocatable  :: rv(:)
  integer :: nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "gti graph-attached value map suite"
  write(*,'(1x,a)') "============================================="

  call a % declare()
  call b % declare()
  call c % declare()

  !-------------------------------------------------------------------!
  ! The status walk of one seat.
  !-------------------------------------------------------------------!

  call report(.not. map % attached(a) .and. &
       & map % status_of(a) == GTI_VALUE_STATUS_UNATTACHED, &
       & "a fresh map holds nothing: unattached is the lawful answer", nfail)

  call map % attach_unknown(a)
  call report(map % attached(a) .and. &
       & map % status_of(a) == GTI_VALUE_STATUS_UNKNOWN, &
       & "attach opens the seat: attached and UNKNOWN", nfail)

  call values % set_real([2.5_dp])
  call map % mark_known(a, values)
  call map % value_of(a, stored)
  call stored % get_real(rv)
  call report(map % status_of(a) == GTI_VALUE_STATUS_KNOWN .and. &
       & size(rv) == 1 .and. abs(rv(1) - 2.5_dp) < 1.0e-14_dp, &
       & "mark_known vouches: KNOWN, and the value reads back", nfail)

  call values % set_real([3.5_dp, 4.5_dp])
  call map % mark_known(a, values)
  call map % value_of(a, stored)
  call stored % get_real(rv)
  call report(size(rv) == 2 .and. abs(rv(1) - 3.5_dp) < 1.0e-14_dp .and. &
       & abs(rv(2) - 4.5_dp) < 1.0e-14_dp, &
       & "an update overwrites: values move, identity does not", nfail)

  call map % mark_unknown(a)
  call report(map % attached(a) .and. &
       & map % status_of(a) == GTI_VALUE_STATUS_UNKNOWN, &
       & "mark_unknown withdraws trust: the seat stays, UNKNOWN", nfail)

  call values % set_real([7.0_dp])
  call map % mark_known(a, values)
  call report(map % status_of(a) == GTI_VALUE_STATUS_KNOWN, &
       & "trust can return: the seat survives its cycles", nfail)

  call map % detach(a)
  call report(.not. map % attached(a) .and. &
       & map % status_of(a) == GTI_VALUE_STATUS_UNATTACHED, &
       & "detach closes the seat: unattached again", nfail)

  !-------------------------------------------------------------------!
  ! Many seats: independence, equal values, and the position law.
  !-------------------------------------------------------------------!

  call map % attach_unknown(a)
  call map % attach_unknown(b)
  call map % attach_unknown(c)

  call values % set_real([9.0_dp])
  call map % mark_known(b, values)

  call report(map % status_of(a) == GTI_VALUE_STATUS_UNKNOWN .and. &
       & map % status_of(b) == GTI_VALUE_STATUS_KNOWN .and. &
       & map % status_of(c) == GTI_VALUE_STATUS_UNKNOWN, &
       & "seats are independent: one vouching moves one status", nfail)

  call map % mark_known(c, values)
  call map % detach(a)
  call map % value_of(b, stored)
  call stored % get_real(rv)
  call report(map % attached(b) .and. map % attached(c) .and. &
       & size(rv) == 1 .and. abs(rv(1) - 9.0_dp) < 1.0e-14_dp .and. &
       & map % status_of(c) == GTI_VALUE_STATUS_KNOWN, &
       & "detaching an early row re-aims nothing: identity, not position", nfail)

  call report(.not. map % attached(a), &
       & "equal values are not identity: a left while c stands", nfail)

  !-------------------------------------------------------------------!
  ! Copies and lifetimes: a copy IS the identity, and the map
  ! outlives the variable that built its rows.
  !-------------------------------------------------------------------!

  b_copy = b
  call map % value_of(b_copy, stored)
  call stored % get_real(rv)
  call report(map % attached(b_copy) .and. &
       & size(rv) == 1 .and. abs(rv(1) - 9.0_dp) < 1.0e-14_dp, &
       & "a copied graph carries the identity: the copy reads b's seat", nfail)

  block
    type(graph) :: temporary
    call temporary % declare()
    keeper = temporary
    call map % attach_unknown(temporary)
  end block
  call report(map % attached(keeper) .and. &
       & map % status_of(keeper) == GTI_VALUE_STATUS_UNKNOWN, &
       & "the map outlives its builder: the token was copied, not borrowed", nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all attached value map checks passed"
  else
     error stop
  end if

contains

  subroutine report(ok, label, nfail)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: label
    integer, intent(inout) :: nfail
    if (ok) then
       write(*,'(1x,a,a)') "PASS : ", label
    else
       write(*,'(1x,a,a)') "FAIL : ", label
       nfail = nfail + 1
    end if
  end subroutine report

end program test_gti_attached_value_map
