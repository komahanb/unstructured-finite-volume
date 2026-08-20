!=====================================================================!
! PARTITIONED IMPLICIT PDE TOWER . LEVEL 0 . CARRIER
!
! The level answers one question: WHAT ARE THE SETS, before there is
! incidence, graph topology, partitioning, fields or a PDE.
!
!      V = { 1 2 3 4 5 6 }         global vertices
!      E = { e1 e2 e3 e4 e5 }      global edges
!      K = { part1 part2 }         partition labels
!
! Before the chain is a graph, its members and its partition labels
! are merely carriers. Nothing joins them yet: no edge knows a
! vertex, no part owns anything, and the word "borrowed" has no
! meaning at this level.
!
! The rung's hazard is worth naming at once. All three carriers
! enumerate from one, so the integer 1 is simultaneously a vertex,
! an edge and a part. No count and no numeral distinguishes these
! sets - only identity does, and every level above depends on that
! being true here first.
!
! The level imports carriers and the assert module, and NOTHING
! relational. Its fixture, chain_carriers_fixture, is earned here;
! the relation fixture one rung up is out of reach by construction,
! and the import gate proves that mechanically rather than by
! promise.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program partitioned_pde_level_0

  use partitioned_pde_assert , only : report, verdict
  use graph_fractal        , only : graph
  use map_set        , only : set_map
  use chain_carriers_fixture , only : chain_carriers

  implicit none

  type(graph) :: v, e, k
  type(set_map) :: sets
  integer           :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "partitioned pde tower . level 0 . carrier"
  write(*,'(1x,a)') "============================================="

  call chain_carriers(sets, v, e, k)

  call check_cardinalities(nfail)
  call check_identities(nfail)
  call check_enumeration(nfail)
  call check_boundaries(nfail)

  call verdict(nfail, "level 0")

contains

  subroutine check_cardinalities(nfail)

    integer, intent(inout) :: nfail

    call report(sets % num_members_of(v) .eq. 6, "V counts six global vertices", nfail)
    call report(sets % num_members_of(e) .eq. 5, "E counts five global edges", nfail)
    call report(sets % num_members_of(k) .eq. 2, "K counts two partition labels", nfail)

  end subroutine check_cardinalities

  !===================================================================!
  ! THE rung's truth: three different sets, and the numerals cannot
  ! tell them apart.
  !===================================================================!

  subroutine check_identities(nfail)

    integer, intent(inout) :: nfail

    call report(.not. v % same_as(e), "V is not E", nfail)
    call report(.not. v % same_as(k), "V is not K", nfail)
    call report(.not. e % same_as(k), "E is not K", nfail)

    call report(sets % has(v, 1) .and. sets % has(e, 1) .and. sets % has(k, 1), &
         & "and the integer 1 is a member of ALL THREE: identity " // &
         & "does this work, never the numerals", nfail)

  end subroutine check_identities

  subroutine check_enumeration(nfail)

    integer, intent(inout) :: nfail

    call report(round_trips(v) .and. round_trips(e) .and. &
         &      round_trips(k), &
         & "member and local_index invert on every carrier", nfail)

  end subroutine check_enumeration

  subroutine check_boundaries(nfail)

    integer, intent(inout) :: nfail

    call report(.not. sets % has(v, 7) .and. .not. sets % has(v, 0), &
         & "an outsider is rejected by V", nfail)
    call report(.not. sets % has(e, 6) .and. .not. sets % has(k, 3), &
         & "and by E and K, each at its own edge", nfail)

  end subroutine check_boundaries

  logical function round_trips(s)

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

end program partitioned_pde_level_0
