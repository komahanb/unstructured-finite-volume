!=====================================================================!
! The calculator tower's common ground: assertion helpers and the
! symbolic names of the persistent object,
!
!      X = { a b c d e }      O = { +  x }      P = { in1 in2 out }
!
! DEPENDENCY-FREE ON PURPOSE (CALCULATOR.md 5, 22): this module
! imports nothing from the framework, so it can never become the
! back door through which a higher level enters a lower test. The
! constants are only names for integers; every level constructs its
! own framework objects.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module calculator_assert

  implicit none

  private
  public :: report, verdict
  public :: SLOT_A, SLOT_B, SLOT_C, SLOT_D, SLOT_E
  public :: OP_PLUS, OP_TIMES
  public :: PORT_IN1, PORT_IN2, PORT_OUT

  !===================================================================!
  ! The value slots of X, the operations of O, the ports of P -
  ! members named for the reader, stored as the small integers the
  ! carriers count.
  !===================================================================!

  integer, parameter :: SLOT_A = 1
  integer, parameter :: SLOT_B = 2
  integer, parameter :: SLOT_C = 3
  integer, parameter :: SLOT_D = 4
  integer, parameter :: SLOT_E = 5

  integer, parameter :: OP_PLUS  = 1
  integer, parameter :: OP_TIMES = 2

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
       write(*,'(1x,a,a,a)') "calculator ", level, ": all truths hold"
    else
       write(*,'(1x,a,a,a,i0,a)') "calculator ", level, ": ", nfail, &
            & " truth(s) FAILED"
       error stop
    end if

  end subroutine verdict

end module calculator_assert
