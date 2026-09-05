!=====================================================================!
! TIME INTEGRATION TOWER . LEVEL 0 . CARRIER
!
! The level answers one question: WHAT SETS EXIST BEFORE TIME HAS
! DIRECTION OR VALUES.
!
!      Q = { x  y }                 state coordinates
!      T = { t0 t1 t2 t3 t4 }       time instants
!      E = { e1 e2 e3 e4 }          time steps
!
! Three carriers, and nothing joins them. No step knows an instant,
! no instant carries a value, and the words BEFORE and AFTER have no
! meaning at this level - direction arrives one rung up, as relation
! structure.
!
! WHY Q IS HERE AT THE BOTTOM. Time integration has two
! conceptually independent axes, the state coordinate and the time
! instant, and the whole tower is an argument that they must not be
! collapsed. That argument has to start before there is anything to
! collapse them WITH: no value, no field, no graph. So Q is
! declared here, beside T, answering nothing.
!
! The rung's hazard, met by a second independent specimen after the
! partitioned tower's V/E/K: all three carriers enumerate from one,
! so the integer 1 is simultaneously x, t0 and e1. No count and no
! numeral distinguishes these sets - only identity does.
!
! The level imports carriers and the assert module, and NOTHING
! relational. Its fixture is earned here; the relation fixture one
! rung up is out of reach by construction, and the import gate
! proves that mechanically rather than by promise.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program time_level_0

  use time_assert          , only : report, verdict
  use time_assert          , only : NQ, NT, NE
  use time_assert          , only : C_X, T0, T4, E1, E4
  use graph_fractal        , only : graph
  use map_set        , only : set_map
  use time_carriers_fixture, only : time_carriers

  implicit none

  type(graph) :: q, t, e
  type(set_map) :: sets
  integer           :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "time integration tower . level 0 . carrier"
  write(*,'(1x,a)') "============================================="

  call time_carriers(sets, q, t, e)

  call check_cardinalities(nfail)
  call check_identities(nfail)
  call check_enumeration(nfail)
  call check_boundaries(nfail)

  call verdict(nfail, "level 0")

contains

  subroutine check_cardinalities(nfail)

    integer, intent(inout) :: nfail

    call report(sets % num_members_of(q) .eq. NQ, "Q counts two state coordinates", nfail)
    call report(sets % num_members_of(t) .eq. NT, "T counts five time instants", nfail)
    call report(sets % num_members_of(e) .eq. NE, "E counts four time steps", nfail)

  end subroutine check_cardinalities

  !===================================================================!
  ! THE rung's truth, and the one the whole tower is built to keep:
  ! three different sets, and the numerals cannot tell them apart.
  !===================================================================!

  subroutine check_identities(nfail)

    integer, intent(inout) :: nfail

    call report(.not. q % same_as(t), &
         & "Q is not T: a state coordinate is not a time instant", nfail)
    call report(.not. q % same_as(e), &
         & "Q is not E: a state coordinate is not a time step", nfail)
    call report(.not. t % same_as(e), &
         & "T is not E: an instant is not a step between instants", nfail)

    call report(sets % has(q, 1) .and. sets % has(t, 1) .and. sets % has(e, 1), &
         & "and the integer 1 is a member of ALL THREE - it is x, " // &
         & "t0 and e1 at once: identity does this work, never the " // &
         & "numerals", nfail)

    call report(C_X .eq. T0 .and. T0 .eq. E1, &
         & "the names agree with that hazard rather than hiding " // &
         & "it: C_X, T0 and E1 are the same integer", nfail)

  end subroutine check_identities

  subroutine check_enumeration(nfail)

    integer, intent(inout) :: nfail

    call report(round_trips(sets, q) .and. round_trips(sets, t) .and. &
         &      round_trips(sets, e), &
         & "member and local_index invert on every carrier", nfail)

  end subroutine check_enumeration

  !===================================================================!
  ! Each carrier refuses at its own edge, and the edges differ -
  ! which is one more thing a numeral cannot tell you.
  !===================================================================!

  subroutine check_boundaries(nfail)

    integer, intent(inout) :: nfail

    call report(.not. sets % has(q, 0) .and. .not. sets % has(q, NQ + 1), &
         & "an outsider is rejected by Q", nfail)
    call report(.not. sets % has(t, 0) .and. .not. sets % has(t, NT + 1), &
         & "and by T, at its own edge", nfail)
    call report(.not. sets % has(e, 0) .and. .not. sets % has(e, NE + 1), &
         & "and by E, at another", nfail)

    call report(sets % has(t, T4) .and. .not. sets % has(e, T4) .and. &
         &      .not. sets % has(q, T4), &
         & "the last instant t4 is a member of T alone - E and Q " // &
         & "are not that long", nfail)
    call report(sets % has(e, E4) .and. sets % has(t, E4) .and. .not. sets % has(q, E4), &
         & "while e4's numeral IS a member of T, meaning t3, and " // &
         & "means nothing at all in Q", nfail)

  end subroutine check_boundaries

  !===================================================================!
  ! The two enumeration laws, asked of the map that answers them. The
  ! map is a dummy, not a local: the identity handed in is described
  ! THERE and nowhere else, and a second map would describe nothing.
  !===================================================================!

  logical function round_trips(sets, s)

    type(set_map)  , intent(in) :: sets
    type(graph), intent(in) :: s

    integer :: i, m

    round_trips = .true.
    do i = 1, sets % num_members_of(s)
       m = sets % member_of(s, i)
       round_trips = round_trips .and. &
            & (sets % member_of(s, sets % index_in(s, m)) .eq. m)
       round_trips = round_trips .and. (sets % index_in(s, m) .eq. i)
    end do

  end function round_trips

end program time_level_0
