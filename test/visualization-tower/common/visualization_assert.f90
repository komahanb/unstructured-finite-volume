!=====================================================================!
! The visualization tower's common ground: assertion helpers and the
! symbolic names of the persistent structural specimen - a chain of
! three mathematical operators,
!
!      X0 --A1--> X1 --A2--> X2 --A3--> X3
!
! over seven structurally distinct sets,
!
!      X0 = { a  b  c  d }              state / domain sets
!      X1 = { p  q  r }
!      X2 = { u  v  w }
!      X3 = { m  n }
!
!      E1 = { e11 e12 e13 e14 e15 }     dependency OCCURRENCES
!      E2 = { e21 e22 e23 e24 }
!      E3 = { e31 e32 e33 }
!
! DEPENDENCY-FREE ON PURPOSE, and deliberately none of the six
! sibling towers' assert modules: this is an independent seventh
! client. It imports nothing from the framework, so it can never
! become the back door through which a higher level enters a lower
! test.
!
!                     NOT ONE NUMBER LIVES HERE
!
! The six towers before this one all reached a level where an oracle
! arrived - a step size, an initial state, a trajectory. This module
! holds NO real constant of any kind, and that absence is the
! tower's opening claim rather than an omission.
!
! The operators A1, A2, A3 have no coefficients at Gate A. They are
! not instantiated numerically, they are never evaluated, and no
! vector is ever pushed through them. What DOES exist is their
! structural dependency,
!
!      D1 : X0 -> X1,   D2 : X1 -> X2,   D3 : X2 -> X3
!
! which the levels above derive from primitive incidence and render.
! The hypothesis under test is that the structural picture can
! precede the numbers entirely - so the day a coefficient appears in
! this file, Gate A's central claim has been abandoned.
!
!                       THE NAMES EARN THEIR KEEP
!
! Every set counts from one, so the member written a is 1 and so
! is p, and so is u, and so is m, and so is e11 - five different
! mathematical objects wearing one integer. Nothing above may write
! a bare integer for a member; it writes X0_A, X1_P, X2_U, X3_M,
! E1_1 and their kin, and the numbering lives in exactly one place.
!
! That the same raw integer inhabits several sets is not an
! accident of the specimen to be tidied away. It is Level 0's
! sharpest law: identity is DECLARED, never inferred from what a
! member happens to be called or from how many there are.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module visualization_assert

  implicit none

  private
  public :: report, verdict
  public :: NX0, NX1, NX2, NX3
  public :: NE1, NE2, NE3
  public :: X0_A, X0_B, X0_C, X0_D
  public :: X1_P, X1_Q, X1_R
  public :: X2_U, X2_V, X2_W
  public :: X3_M, X3_N
  public :: E1_1, E1_2, E1_3, E1_4, E1_5
  public :: E2_1, E2_2, E2_3, E2_4
  public :: E3_1, E3_2, E3_3
  public :: ND1, ND2, ND3, ND21, ND31

  !-------------------------------------------------------------------!
  ! The seven cardinalities. |X1| = |X2| = 3 on purpose: two sets
  ! of one size that are not one set is the first thing Level 0
  ! must be able to say.
  !-------------------------------------------------------------------!

  integer, parameter :: NX0 = 4
  integer, parameter :: NX1 = 3
  integer, parameter :: NX2 = 3
  integer, parameter :: NX3 = 2

  integer, parameter :: NE1 = 5
  integer, parameter :: NE2 = 4
  integer, parameter :: NE3 = 3

  !-------------------------------------------------------------------!
  ! The members, by name. The overlap is the point: X0_A, X1_P, X2_U,
  ! X3_M and E1_1 are all the integer 1, and all five are different.
  !-------------------------------------------------------------------!

  integer, parameter :: X0_A = 1, X0_B = 2, X0_C = 3, X0_D = 4
  integer, parameter :: X1_P = 1, X1_Q = 2, X1_R = 3
  integer, parameter :: X2_U = 1, X2_V = 2, X2_W = 3
  integer, parameter :: X3_M = 1, X3_N = 2

  integer, parameter :: E1_1 = 1, E1_2 = 2, E1_3 = 3, E1_4 = 4, E1_5 = 5
  integer, parameter :: E2_1 = 1, E2_2 = 2, E2_3 = 3, E2_4 = 4
  integer, parameter :: E3_1 = 1, E3_2 = 2, E3_3 = 3

  !-------------------------------------------------------------------!
  ! The tuple counts the derived dependencies must have. These are
  ! ORACLES, worked out on paper from the twelve occurrences and
  ! written here so no level may read them off the machinery it is
  ! testing.
  !
  !      |D1|   = 5     one per occurrence of E1 - nothing collapses
  !      |D2|   = 4     one per occurrence of E2
  !      |D3|   = 3     one per occurrence of E3
  !
  !      |D2:1| = 6     NOT 7. The seven source-to-sink walks
  !                     through X1 include b->p->u and b->q->u, two
  !                     witnesses of the one dependency b->u.
  !
  !      |D3:1| = 6     the skeleton of the whole chain
  !-------------------------------------------------------------------!

  integer, parameter :: ND1  = 5
  integer, parameter :: ND2  = 4
  integer, parameter :: ND3  = 3
  integer, parameter :: ND21 = 6
  integer, parameter :: ND31 = 6

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
       write(*,'(1x,a,a)') "visualization ", level // ": all truths hold"
    else
       write(*,'(1x,a,a)') "visualization ", level // ": FAILED"
       error stop 1
    end if

  end subroutine verdict

end module visualization_assert
