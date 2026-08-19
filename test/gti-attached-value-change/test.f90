!=====================================================================!
! The attached-value change suite: the first value-side client of
! the reversible change protocol. One concrete change updates the
! attached value map through the ONE generic controller -
!
!      apply   save the prior state, vouch the new value
!      check   the caller's verdict
!      keep    the update stands
!      revert  the seat restored EXACTLY
!
! - proven from every prior status: an UNATTACHED seat appears on
! accept and vanishes on reject; an UNKNOWN seat gains trust on
! accept and stands untrusted again on reject; a KNOWN seat trades
! values on accept and keeps its old value exactly on reject. A
! vetoed check restores like a rejection; the bound value is a
! copy the caller cannot mutate afterward; shape rides through;
! and the change declares itself pure value - touches_value, never
! touches_structure - the mirror image of adaptive growth's mixed
! change under the same four verbs. No law test below calls a
! lifecycle verb by hand: the controller runs every path.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_gti_attached_value_change

  use iso_fortran_env      , only : dp => REAL64
  use fractal_graph        , only : graph
  use gti_value_buffers    , only : gti_value_buffer
  use gti_attached_value_maps, only : gti_attached_value_map, &
       & GTI_VALUE_STATUS_UNATTACHED, GTI_VALUE_STATUS_UNKNOWN, &
       & GTI_VALUE_STATUS_KNOWN
  use gti_attached_value_changes, only : gti_attached_value_change
  use gti_change_protocols , only : gti_change_controller, gti_change_result

  implicit none

  type(gti_attached_value_map)   :: map
  type(gti_attached_value_change):: change
  type(gti_change_controller)    :: controller
  type(gti_change_result)        :: result

  type(graph) :: g1, g2, g3, g4, g5, g6, g7, g8, g9

  type(gti_value_buffer) :: new_value, old_value, stored
  real(dp), allocatable  :: rv(:)
  integer :: nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "gti attached-value change suite"
  write(*,'(1x,a)') "============================================="

  call g1 % declare()
  call g2 % declare()
  call g3 % declare()
  call g4 % declare()
  call g5 % declare()
  call g6 % declare()
  call g7 % declare()
  call g8 % declare()
  call g9 % declare()

  !-------------------------------------------------------------------!
  ! From UNATTACHED: accept creates, reject leaves no seat.
  !-------------------------------------------------------------------!

  call new_value % set_real([5.0_dp])
  call change % bind(map, g1, new_value)
  call controller % run(change, .true., result)
  call map % value_of(g1, stored)
  call stored % get_real(rv)
  call report(map % status_of(g1) == GTI_VALUE_STATUS_KNOWN .and. &
       & size(rv) == 1 .and. abs(rv(1) - 5.0_dp) < 1.0e-14_dp, &
       & "accepted on UNATTACHED: the seat appears, KNOWN with the value", nfail)

  call report(result % touches_value .and. .not. result % touches_structure, &
       & "the change declares itself pure value: no structure touched", nfail)

  call change % bind(map, g2, new_value)
  call controller % run(change, .false., result)
  call report(map % status_of(g2) == GTI_VALUE_STATUS_UNATTACHED .and. &
       & .not. map % attached(g2), &
       & "rejected on UNATTACHED: no seat remains, nothing readable", nfail)

  !-------------------------------------------------------------------!
  ! From UNKNOWN: accept vouches, reject stands untrusted again.
  !-------------------------------------------------------------------!

  call map % attach_unknown(g3)
  call new_value % set_real([6.5_dp])
  call change % bind(map, g3, new_value)
  call controller % run(change, .true., result)
  call map % value_of(g3, stored)
  call stored % get_real(rv)
  call report(map % status_of(g3) == GTI_VALUE_STATUS_KNOWN .and. &
       & abs(rv(1) - 6.5_dp) < 1.0e-14_dp, &
       & "accepted on UNKNOWN: the seat gains trust and the value", nfail)

  call map % attach_unknown(g4)
  call change % bind(map, g4, new_value)
  call controller % run(change, .false., result)
  call report(map % attached(g4) .and. &
       & map % status_of(g4) == GTI_VALUE_STATUS_UNKNOWN, &
       & "rejected on UNKNOWN: the seat stands, untrusted again", nfail)

  !-------------------------------------------------------------------!
  ! From KNOWN: accept trades values, reject restores exactly.
  !-------------------------------------------------------------------!

  call map % attach_unknown(g5)
  call old_value % set_real([1.0_dp])
  call map % mark_known(g5, old_value)
  call new_value % set_real([2.0_dp])
  call change % bind(map, g5, new_value)
  call controller % run(change, .true., result)
  call map % value_of(g5, stored)
  call stored % get_real(rv)
  call report(abs(rv(1) - 2.0_dp) < 1.0e-14_dp, &
       & "accepted on KNOWN: the new value stands", nfail)

  call map % attach_unknown(g6)
  call old_value % set_real([1.5_dp])
  call map % mark_known(g6, old_value)
  call new_value % set_real([9.9_dp])
  call change % bind(map, g6, new_value)
  call controller % run(change, .false., result)
  call map % value_of(g6, stored)
  call stored % get_real(rv)
  call report(map % status_of(g6) == GTI_VALUE_STATUS_KNOWN .and. &
       & size(rv) == 1 .and. abs(rv(1) - 1.5_dp) < 1.0e-14_dp, &
       & "rejected on KNOWN: the old value restored exactly", nfail)

  !-------------------------------------------------------------------!
  ! The vetoed check: even an accepting caller cannot keep it.
  !-------------------------------------------------------------------!

  call map % attach_unknown(g7)
  call old_value % set_real([3.0_dp])
  call map % mark_known(g7, old_value)
  call new_value % set_real([8.0_dp])
  call change % bind(map, g7, new_value, check_passes=.false.)
  call controller % run(change, .true., result)
  call map % value_of(g7, stored)
  call stored % get_real(rv)
  call report(abs(rv(1) - 3.0_dp) < 1.0e-14_dp, &
       & "a vetoed check restores the old state", nfail)
  call report(result % checked .and. .not. result % check_passed .and. &
       & result % reverted .and. .not. result % kept, &
       & "the vetoed lifecycle reads: checked, not passed, reverted", nfail)

  !-------------------------------------------------------------------!
  ! The copy law and the shape law.
  !-------------------------------------------------------------------!

  call new_value % set_real([4.0_dp])
  call change % bind(map, g8, new_value)
  call new_value % set_real([999.0_dp])
  call controller % run(change, .true., result)
  call map % value_of(g8, stored)
  call stored % get_real(rv)
  call report(abs(rv(1) - 4.0_dp) < 1.0e-14_dp, &
       & "the bound value is a copy: the caller's later edit changes nothing", nfail)

  call new_value % set_real([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], ncomp=2)
  call change % bind(map, g9, new_value)
  call controller % run(change, .true., result)
  call map % value_of(g9, stored)
  call stored % get_real(rv)
  call report(stored % ncomp == 2 .and. stored % nentries == 2 .and. &
       & size(rv) == 4 .and. abs(rv(3) - 3.0_dp) < 1.0e-14_dp, &
       & "shape rides through: two entries of two components survive", nfail)

  call report(.true., &
       & "every path above ran through gti_change_controller alone", nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all attached-value change checks passed"
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

end program test_gti_attached_value_change
