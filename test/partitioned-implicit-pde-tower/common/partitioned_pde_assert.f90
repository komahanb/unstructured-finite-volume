!=====================================================================!
! The partitioned implicit PDE tower's common ground: assertion
! helpers and the persistent specimen's oracles.
!
!      G:  1 -e1-> 2 -e2-> 3 -e3-> 4 -e4-> 5 -e5-> 6
!
!      A = 2I - L,      A q = b
!
!      q* = [1, 2, 4, 7, 11, 16]
!      b  = [1, 3, 7, 13, 21, 37]
!
! DEPENDENCY-FREE ON PURPOSE, and deliberately none of the four
! sibling towers' assert modules: this is an independent fifth
! client, and the first outside the derivative family. It imports
! nothing from the framework, so it can never become the back door
! through which a higher gate enters a lower test.
!
! The oracles live here because every gate checks against them, but
! they are ORACLES: no gate is permitted to obtain b from a matrix,
! and nothing here knows what a Laplacian is.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module partitioned_pde_assert

  use iso_fortran_env, only : dp => REAL64

  implicit none

  private
  public :: report, verdict
  public :: NV, NE, Q_EXACT, B_EXACT, L_EXACT

  ! The global chain: six vertices, five edges.
  integer, parameter :: NV = 6
  integer, parameter :: NE = 5

  ! The oracles, in GLOBAL vertex order.
  real(dp), parameter :: Q_EXACT(NV) = &
       & [1.0_dp, 2.0_dp, 4.0_dp, 7.0_dp, 11.0_dp, 16.0_dp]
  real(dp), parameter :: L_EXACT(NV) = &
       & [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, -5.0_dp]
  real(dp), parameter :: B_EXACT(NV) = &
       & [1.0_dp, 3.0_dp, 7.0_dp, 13.0_dp, 21.0_dp, 37.0_dp]

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

  subroutine verdict(nfail, gate)

    integer         , intent(in) :: nfail
    character(len=*), intent(in) :: gate

    write(*,'(1x,a)') "============================================="
    if (nfail .eq. 0) then
       write(*,'(1x,a,a,a)') "partitioned pde ", gate, ": all truths hold"
    else
       write(*,'(1x,a,a,a,i0,a)') "partitioned pde ", gate, ": ", &
            & nfail, " truth(s) FAILED"
       error stop
    end if

  end subroutine verdict

end module partitioned_pde_assert
