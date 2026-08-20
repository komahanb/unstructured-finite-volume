!=====================================================================!
! CALCULATOR TOWER . LEVEL 0 . CARRIERS
!
! The level answers one question: WHAT MEMBERS MAY EXIST. Three
! independent domains are declared,
!
!      X = { a b c d e }        the value slots
!      O = { +  x }             the operations
!      P = { in1 in2 out }      the ports
!
! and nothing else is true yet. No relation, no graph, no value, no
! arithmetic: the symbol + is only a member of O here. The imports
! of this file ARE the negative truth of the level - the set modules
! and nothing above it (CALCULATOR.md 7).
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program calculator_level_0

  use calculator_assert, only : report, verdict
  use calculator_assert, only : SLOT_A, SLOT_E, OP_PLUS, PORT_OUT
  use graph_fractal        , only : set_graph => graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set        , only : set_map

  implicit none

  type(set_graph) :: x, o, p
  integer           :: nfail
  type(set_map)     :: sets

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "calculator tower . level 0 . carriers"
  write(*,'(1x,a)') "============================================="

  ! The three declarations. This is everything Level 0 may do.
  call x % declare()
  call sets % bind(x, counted_set_representation(5))
  call o % declare()
  call sets % bind(o, counted_set_representation(2))
  call p % declare()
  call sets % bind(p, counted_set_representation(3))

  call check_cardinalities(nfail)
  call check_distinct_identities(nfail)
  call check_enumeration_round_trips(nfail)
  call check_membership_boundary(nfail)

  call verdict(nfail, "level 0")

contains

  !===================================================================!
  ! |X| = 5, |O| = 2, |P| = 3 - exactly.
  !===================================================================!

  subroutine check_cardinalities(nfail)

    integer, intent(inout) :: nfail

    call report(sets % size_of(x) .eq. 5, &
         & "X counts its five value slots", nfail)
    call report(sets % size_of(o) .eq. 2, &
         & "O counts its two operations", nfail)
    call report(sets % size_of(p) .eq. 3, &
         & "P counts its three ports", nfail)

  end subroutine check_cardinalities

  !===================================================================!
  ! Three declarations, three worlds: no pair is the same domain,
  ! whatever the integers inside may suggest.
  !===================================================================!

  subroutine check_distinct_identities(nfail)

    integer, intent(inout) :: nfail

    call report(.not. x % same_as(o), &
         & "X is not O", nfail)
    call report(.not. x % same_as(p), &
         & "X is not P", nfail)
    call report(.not. o % same_as(p), &
         & "O is not P", nfail)

  end subroutine check_distinct_identities

  !===================================================================!
  ! The two enumeration laws, on every member of every carrier:
  !
  !      member(local_index(m)) = m
  !      local_index(member(i)) = i
  !
  ! which together make enumeration injective.
  !===================================================================!

  subroutine check_enumeration_round_trips(nfail)

    integer, intent(inout) :: nfail

    logical :: ok
    integer :: i

    ok = .true.
    do i = 1, sets % size_of(x)
       ok = ok .and. (sets % member_of(x, sets % index_in(x, sets % member_of(x, i))) &
            &         .eq. sets % member_of(x, i))
       ok = ok .and. (sets % index_in(x, sets % member_of(x, i)) .eq. i)
    end do
    call report(ok, "member and local_index invert on X", nfail)

    ok = .true.
    do i = 1, sets % size_of(o)
       ok = ok .and. (sets % member_of(o, sets % index_in(o, sets % member_of(o, i))) &
            &         .eq. sets % member_of(o, i))
       ok = ok .and. (sets % index_in(o, sets % member_of(o, i)) .eq. i)
    end do
    call report(ok, "member and local_index invert on O", nfail)

    ok = .true.
    do i = 1, sets % size_of(p)
       ok = ok .and. (sets % member_of(p, sets % index_in(p, sets % member_of(p, i))) &
            &         .eq. sets % member_of(p, i))
       ok = ok .and. (sets % index_in(p, sets % member_of(p, i)) .eq. i)
    end do
    call report(ok, "member and local_index invert on P", nfail)

  end subroutine check_enumeration_round_trips

  !===================================================================!
  ! Membership knows its own: a, e and + and out belong where they
  ! were declared; the outsider belongs nowhere.
  !===================================================================!

  subroutine check_membership_boundary(nfail)

    integer, intent(inout) :: nfail

    call report(sets % has_in(x, SLOT_A) .and. sets % has_in(x, SLOT_E), &
         & "a and e are value slots", nfail)
    call report(sets % has_in(o, OP_PLUS) .and. sets % has_in(p, PORT_OUT), &
         & "+ is an operation and out is a port", nfail)

    call report(.not. sets % has_in(x, 6) .and. .not. sets % has_in(x, 0), &
         & "an outsider is rejected by X", nfail)
    call report(.not. sets % has_in(o, 3), &
         & "and by O", nfail)
    call report(.not. sets % has_in(p, 4), &
         & "and by P", nfail)

  end subroutine check_membership_boundary

end program calculator_level_0
