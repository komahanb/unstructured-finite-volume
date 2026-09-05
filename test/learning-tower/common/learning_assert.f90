!=====================================================================!
! The learning tower's common ground: assertion helpers and the
! symbolic names of the persistent object,
!
!      V = { w  x  yhat  y  e }      the value slots
!      O = { predict  error }        the operations
!      P = { in1  in2  out }         the ports
!
! DEPENDENCY-FREE ON PURPOSE, and DELIBERATELY NOT calculator_assert:
! the learning tower is an independent second client of the
! production architecture, and shares no fixture with the first.
! This module imports nothing from the framework, so it can never
! become the back door through which a higher level enters a lower
! test. Below the constitution, predict and error are only members
! of O - they do not yet multiply or subtract anything.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module learning_assert

  implicit none

  private
  public :: report, verdict
  public :: SLOT_W, SLOT_X, SLOT_YHAT, SLOT_Y, SLOT_E
  public :: OP_PREDICT, OP_ERROR
  public :: PORT_IN1, PORT_IN2, PORT_OUT

  integer, parameter :: SLOT_W    = 1
  integer, parameter :: SLOT_X    = 2
  integer, parameter :: SLOT_YHAT = 3
  integer, parameter :: SLOT_Y    = 4
  integer, parameter :: SLOT_E    = 5

  integer, parameter :: OP_PREDICT = 1
  integer, parameter :: OP_ERROR   = 2

  integer, parameter :: PORT_IN1 = 1
  integer, parameter :: PORT_IN2 = 2
  integer, parameter :: PORT_OUT = 3

contains

  subroutine report(ok, label, nfail)

    logical         , intent(in)    :: ok
    character(len=*), intent(in)    :: label
    integer         , intent(inout) :: nfail

    if (ok) then
       write(*,'(1x,a,a)') "PASS : ", label
    else
       write(*,'(1x,a,a)') "FAIL : ", label
       nfail = nfail + 1
    end if

  end subroutine report

  subroutine verdict(nfail, level)

    integer         , intent(in) :: nfail
    character(len=*), intent(in) :: level

    write(*,'(1x,a)') "============================================="
    if (nfail .eq. 0) then
       write(*,'(1x,a,a,a)') "learning ", level, ": all truths hold"
    else
       write(*,'(1x,a,a,a,i0,a)') "learning ", level, ": ", nfail, &
            & " truth(s) FAILED"
       error stop
    end if

  end subroutine verdict

end module learning_assert
