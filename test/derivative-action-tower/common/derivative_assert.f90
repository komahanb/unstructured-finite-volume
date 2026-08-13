!=====================================================================!
! The derivative action tower's common ground: assertion helpers and
! the symbolic names of the persistent specimen,
!
!      V = { x  y  u  z }            the value slots
!      O = { product  sum }          the operations
!      P = { in1  in2  out }         the ports
!
! for the computation u = product(x, y), z = sum(u, y) - chosen
! because y reaches z along TWO structural routes while x reaches z
! along one. DEPENDENCY-FREE ON PURPOSE, and deliberately neither
! calculator_assert nor learning_assert: the derivative tower is an
! independent third client of the production architecture and shares
! no fixture with the first two. This module imports nothing from
! the framework, so it can never become the back door through which
! a higher level enters a lower test. No numerical value and no
! derivative coefficient lives here: below constitution, product and
! sum are only members of O, and Gate A never differentiates
! anything.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module derivative_assert

  implicit none

  private
  public :: report, verdict
  public :: SLOT_X, SLOT_Y, SLOT_U, SLOT_Z
  public :: OP_PRODUCT, OP_SUM
  public :: PORT_IN1, PORT_IN2, PORT_OUT

  integer, parameter :: SLOT_X = 1
  integer, parameter :: SLOT_Y = 2
  integer, parameter :: SLOT_U = 3
  integer, parameter :: SLOT_Z = 4

  integer, parameter :: OP_PRODUCT = 1
  integer, parameter :: OP_SUM     = 2

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
       write(*,'(1x,a,a,a)') "derivative ", level, ": all truths hold"
    else
       write(*,'(1x,a,a,a,i0,a)') "derivative ", level, ": ", nfail, &
            & " truth(s) FAILED"
       error stop
    end if

  end subroutine verdict

end module derivative_assert
