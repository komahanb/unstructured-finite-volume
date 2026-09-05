!=====================================================================!
! The adjoint tower's common ground: assertion helpers and the
! symbolic names of the persistent specimen,
!
!      V = { p  u  v }        the variable carrier
!      T = { r1  r2  f }      the target carrier
!
! carrying the four roles as SUBOBJECTS - P = {p}, Q = {u,v} in V;
! Y = {r1,r2}, Z = {f} in T - for the implicit problem
!
!      R(q,p) = 0,      f(q,p),      df/dp = f_p - lambda^T R_p.
!
! DEPENDENCY-FREE ON PURPOSE, and deliberately none of
! calculator_assert, learning_assert or derivative_assert: the
! adjoint tower is an independent fourth client of the production
! architecture and shares no fixture with its siblings. This module
! imports nothing from the framework, so it can never become the
! back door through which a higher level enters a lower test.
!
! No coefficient and no numerical value lives here. Below Level 8
! the residual rows and the response are only members of T, and
! nothing in this tower knows what 2, 1, 3, 4, -4 or -11 mean.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module adjoint_assert

  implicit none

  private
  public :: report, verdict
  public :: VAR_P, VAR_U, VAR_V
  public :: TGT_R1, TGT_R2, TGT_F

  ! The variable carrier V: the parameter, then the state slots.
  integer, parameter :: VAR_P = 1
  integer, parameter :: VAR_U = 2
  integer, parameter :: VAR_V = 3

  ! The target carrier T: the residual rows, then the response.
  integer, parameter :: TGT_R1 = 1
  integer, parameter :: TGT_R2 = 2
  integer, parameter :: TGT_F  = 3

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
       write(*,'(1x,a,a,a)') "adjoint ", level, ": all truths hold"
    else
       write(*,'(1x,a,a,a,i0,a)') "adjoint ", level, ": ", nfail, &
            & " truth(s) FAILED"
       error stop
    end if

  end subroutine verdict

end module adjoint_assert
