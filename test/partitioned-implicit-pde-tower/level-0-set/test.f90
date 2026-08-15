!=====================================================================!
! PARTITIONED IMPLICIT PDE TOWER . LEVEL 0 . SET
!
! The level answers one question: WHAT ARE THE SETS, before there is
! incidence, graph topology, partitioning, fields or a PDE.
!
!      V = { 1 2 3 4 5 6 }         global vertices
!      E = { e1 e2 e3 e4 e5 }      global edges
!      K = { part1 part2 }         partition labels
!
! Before the chain is a graph, its members and its partition labels
! are merely sets. Nothing joins them yet: no edge knows a
! vertex, no part owns anything, and the word "borrowed" has no
! meaning at this level.
!
! The rung's hazard is worth naming at once. All three sets
! enumerate from one, so the integer 1 is simultaneously a vertex,
! an edge and a part. No count and no numeral distinguishes these
! sets - only identity does, and every level above depends on that
! being true here first.
!
! The level imports sets and the assert module, and NOTHING
! relational. Its fixture, chain_sets_fixture, is earned here;
! the relation fixture one rung up is out of reach by construction,
! and the import gate proves that mechanically rather than by
! promise.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program partitioned_pde_level_0

  use partitioned_pde_assert , only : report, verdict
  use graph_set          , only : index_set, set
  use chain_sets_fixture , only : chain_sets

  implicit none

  type(index_set) :: v, e, k
  integer           :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "partitioned pde tower . level 0 . set"
  write(*,'(1x,a)') "============================================="

  call chain_sets(v, e, k)

  call check_cardinalities(nfail)
  call check_identities(nfail)
  call check_enumeration(nfail)
  call check_boundaries(nfail)

  call verdict(nfail, "level 0")

contains

  subroutine check_cardinalities(nfail)

    integer, intent(inout) :: nfail

    call report(v % size() .eq. 6, "V counts six global vertices", nfail)
    call report(e % size() .eq. 5, "E counts five global edges", nfail)
    call report(k % size() .eq. 2, "K counts two partition labels", nfail)

  end subroutine check_cardinalities

  !===================================================================!
  ! THE rung's truth: three different sets, and the numerals cannot
  ! tell them apart.
  !===================================================================!

  subroutine check_identities(nfail)

    integer, intent(inout) :: nfail

    call report(.not. v % equals(e), "V is not E", nfail)
    call report(.not. v % equals(k), "V is not K", nfail)
    call report(.not. e % equals(k), "E is not K", nfail)

    call report(v % has(1) .and. e % has(1) .and. k % has(1), &
         & "and the integer 1 is a member of ALL THREE: identity " // &
         & "does this work, never the numerals", nfail)

  end subroutine check_identities

  subroutine check_enumeration(nfail)

    integer, intent(inout) :: nfail

    call report(round_trips(v) .and. round_trips(e) .and. &
         &      round_trips(k), &
         & "member and local_index invert on every set", nfail)

  end subroutine check_enumeration

  subroutine check_boundaries(nfail)

    integer, intent(inout) :: nfail

    call report(.not. v % has(7) .and. .not. v % has(0), &
         & "an outsider is rejected by V", nfail)
    call report(.not. e % has(6) .and. .not. k % has(3), &
         & "and by E and K, each at its own edge", nfail)

  end subroutine check_boundaries

  logical function round_trips(s)

    class(set), intent(in) :: s

    integer :: i, m

    round_trips = .true.
    do i = 1, s % size()
       m = s % member(i)
       round_trips = round_trips .and. &
            & (s % member(s % local_index(m)) .eq. m)
       round_trips = round_trips .and. (s % local_index(m) .eq. i)
    end do

  end function round_trips

end program partitioned_pde_level_0
