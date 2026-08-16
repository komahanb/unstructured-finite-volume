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
! THE ORACLES ARRIVE AT LEVEL 5, and not one instant before. Gate A
! had no numbers at all. From Level 5 the specimen acquires them,
! and they live here because every level above checks against them:
!
!      h    = 1/2                  the step size
!      time = [0, 1/2, 1, 3/2, 2]  the instants' coordinates
!      q0   = [2, 0]               the initial state
!
!      S(q) = [ x, y - x ]         the dynamical action
!      qdot = -S(q)                the house convention
!
! and the three discrete answers the levels above must reach:
!
!      q_FE,1   = [1, 1]                forward euler, by hand
!      q_BE,1   = [4/3, 4/9]            backward euler
!      q_BDF2,2 = [5/6, 47/72]          bdf-2, started from q_BE,1
!
! These are ORACLES. No level is permitted to obtain them from the
! machinery it is testing: forward euler is computed here from
! ordinary arithmetic, and the two implicit answers are the exact
! rational solutions of their own residual equations, verified by
! substitution rather than by a solver.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module time_assert

  use iso_fortran_env, only : dp => REAL64

  implicit none

  private
  public :: report, verdict
  public :: NQ, NT, NE
  public :: C_X, C_Y
  public :: T0, T1, T2, T3, T4
  public :: E1, E2, E3, E4
  public :: H_STEP, TIME_COORD, Q0, Q_FE1, Q_BE1, Q_BDF2, TOL
  public :: NSTEPS, TOL_MARCH
  public :: FE_TRAJECTORY, BE_TRAJECTORY, BDF2_TRAJECTORY
  public :: action_of

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

  !===================================================================!
  ! The numerical specimen, earned at Level 5.
  !===================================================================!

  real(dp), parameter :: H_STEP = 0.5_dp

  ! The instants' coordinates, in T's declaration order.
  real(dp), parameter :: TIME_COORD(NT) = &
       & [0.0_dp, 0.5_dp, 1.0_dp, 1.5_dp, 2.0_dp]

  ! The initial state, in Q's declaration order: x then y.
  real(dp), parameter :: Q0(NQ) = [2.0_dp, 0.0_dp]

  !===================================================================!
  ! The three discrete answers, as exact rationals.
  !
  !   forward euler   q1 = q0 - h S(q0) = [2,0] - (1/2)[2,-2]
  !                      = [1, 1]
  !
  !   backward euler  q1 - q0 + h S(q1) = 0, i.e. (I + hM) q1 = q0
  !                   with M = [[1,0],[-1,1]]
  !                      = [4/3, 4/9]
  !
  !   bdf-2           (3/2)q2 - 2 q1 + (1/2)q0 + (1/2) S(q2) = 0
  !                   started from q0 and the backward-euler q1
  !                      = [5/6, 47/72]
  !===================================================================!

  real(dp), parameter :: Q_FE1(NQ) = [1.0_dp, 1.0_dp]

  real(dp), parameter :: Q_BE1(NQ) = &
       & [4.0_dp / 3.0_dp, 4.0_dp / 9.0_dp]

  real(dp), parameter :: Q_BDF2(NQ) = &
       & [5.0_dp / 6.0_dp, 47.0_dp / 72.0_dp]

  real(dp), parameter :: TOL = 1.0e-10_dp

  !===================================================================!
  ! The full four-step marches, earned at Level 8. Each column is one
  ! instant, t0 through t4, and each row a state coordinate:
  !
  !      row 1 = x        row 2 = y
  !
  ! FORWARD EULER    q_{n+1} = q_n - h S(q_n)
  !
  !      x : 2, 1, 1/2, 1/4, 1/8
  !      y : 0, 1, 1,   3/4, 1/2
  !
  ! BACKWARD EULER   (I + hM) q_{n+1} = q_n,  M = [[1,0],[-1,1]]
  !
  !      x : 2, 4/3, 8/9,   16/27, 32/81
  !      y : 0, 4/9, 16/27, 16/27, 128/243
  !
  ! BDF-2, started with one backward-euler step, as production does:
  !
  !      x_{n+1} = x_n - x_{n-1}/4
  !      y_{n+1} = y_n - y_{n-1}/4 + x_{n+1}/4
  !
  !      x : 2, 4/3, 5/6,   1/2, 7/24
  !      y : 0, 4/9, 47/72, 2/3, 83/144
  !
  ! ORACLES, every one. They are the exact rational solutions of the
  ! recurrences above, and no level may obtain them from the
  ! machinery it is testing.
  !===================================================================!

  integer, parameter :: NSTEPS = 4

  real(dp), parameter :: FE_TRAJECTORY(NQ, 0:NSTEPS) = reshape( &
       & [2.0_dp        , 0.0_dp        , &
       &  1.0_dp        , 1.0_dp        , &
       &  0.5_dp        , 1.0_dp        , &
       &  0.25_dp       , 0.75_dp       , &
       &  0.125_dp      , 0.5_dp        ], [NQ, NSTEPS + 1])

  real(dp), parameter :: BE_TRAJECTORY(NQ, 0:NSTEPS) = reshape( &
       & [2.0_dp             , 0.0_dp               , &
       &  4.0_dp/3.0_dp      , 4.0_dp/9.0_dp        , &
       &  8.0_dp/9.0_dp      , 16.0_dp/27.0_dp      , &
       &  16.0_dp/27.0_dp    , 16.0_dp/27.0_dp      , &
       &  32.0_dp/81.0_dp    , 128.0_dp/243.0_dp    ], [NQ, NSTEPS + 1])

  real(dp), parameter :: BDF2_TRAJECTORY(NQ, 0:NSTEPS) = reshape( &
       & [2.0_dp             , 0.0_dp               , &
       &  4.0_dp/3.0_dp      , 4.0_dp/9.0_dp        , &
       &  5.0_dp/6.0_dp      , 47.0_dp/72.0_dp      , &
       &  0.5_dp             , 2.0_dp/3.0_dp        , &
       &  7.0_dp/24.0_dp     , 83.0_dp/144.0_dp     ], [NQ, NSTEPS + 1])

  !===================================================================!
  ! A SEPARATE tolerance for the marched answers, and it is separate
  ! on purpose. Levels 5-7 substitute exact states into residuals and
  ! meet TOL. Level 8 walks four steps of a finite-difference Newton
  ! stack, whose error is a MEASURED quantity, not a guessed one.
  !
  ! MEASURED, on this specimen: the worst deviation from the exact
  ! rational trajectory is 2.22e-16 for backward euler and the same
  ! for bdf-2 - one unit in the last place. The stack converges
  ! essentially exactly here because the residual is affine and the
  ! problem is 2x2.
  !
  ! The constant below therefore sits four orders ABOVE what was
  ! observed: loose enough that a different platform's last bits
  ! cannot fail it, tight enough that it remains a real constraint.
  ! A tolerance chosen to make a test pass without knowing the true
  ! error is how a genuine discrepancy hides, and both marches
  ! print their worst error so the margin stays visible.
  !===================================================================!

  real(dp), parameter :: TOL_MARCH = 1.0e-12_dp

contains

  !===================================================================!
  ! The dynamical action, in the plainest arithmetic there is:
  !
  !      S([x, y]) = [ x, y - x ]
  !
  ! An ORACLE, and deliberately a duplicate of what the Level-6
  ! fixture will compute through production types. A test that
  ! checked production against production would check nothing.
  !===================================================================!

  pure function action_of(q) result(s)

    real(dp), intent(in) :: q(NQ)
    real(dp)             :: s(NQ)

    s(1) = q(1)
    s(2) = q(2) - q(1)

  end function action_of

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
