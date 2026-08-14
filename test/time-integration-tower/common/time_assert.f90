!=====================================================================!
! The time integration tower's common ground: assertion helpers and
! the symbolic names of the persistent structural specimen,
!
!      Q = { x  y }                      state coordinates
!      T = { t0 t1 t2 t3 t4 }            time instants
!      E = { e1 e2 e3 e4 }               time steps
!
! DEPENDENCY-FREE ON PURPOSE, and deliberately none of the five
! sibling towers' assert modules: this is an independent sixth
! client. It imports nothing from the framework, so it can never
! become the back door through which a higher level enters a lower
! test.
!
! THE NAMES EARN THEIR KEEP HERE. The carriers all count from one,
! so the instant written t0 is member 1, and t4 is member 5 - an
! off-by-one waiting to happen every time a test writes a literal.
! Nothing above may write a bare integer for a member; it writes
! T0..T4, E1..E4, C_X, C_Y, and the numbering lives in exactly one
! place.
!
! NO ORACLES LIVE HERE YET. Gate A has no numbers: no state value,
! no step size, no scheme coefficient, no trajectory. When Level 5
! earns them they will arrive beside these names, and not before.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module time_assert

  implicit none

  private
  public :: report, verdict
  public :: NQ, NT, NE
  public :: C_X, C_Y
  public :: T0, T1, T2, T3, T4
  public :: E1, E2, E3, E4

  ! The cardinalities of the three carriers.
  integer, parameter :: NQ = 2
  integer, parameter :: NT = 5
  integer, parameter :: NE = 4

  !===================================================================!
  ! The members, named for the reader and stored as the small
  ! integers the carriers count. Note where the collision bites: the
  ! integer 1 below is C_X, T0 and E1 at once, and those three facts
  ! are about three different sets.
  !===================================================================!

  ! Q - the state coordinates.
  integer, parameter :: C_X = 1
  integer, parameter :: C_Y = 2

  ! T - the time instants. t0 is member ONE, not member zero.
  integer, parameter :: T0 = 1
  integer, parameter :: T1 = 2
  integer, parameter :: T2 = 3
  integer, parameter :: T3 = 4
  integer, parameter :: T4 = 5

  ! E - the time steps.
  integer, parameter :: E1 = 1
  integer, parameter :: E2 = 2
  integer, parameter :: E3 = 3
  integer, parameter :: E4 = 4

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
       write(*,'(1x,a,a)') "time integration ", level // ": all truths hold"
    else
       write(*,'(1x,a,a)') "time integration ", level // ": FAILED"
       error stop 1
    end if

  end subroutine verdict

end module time_assert
