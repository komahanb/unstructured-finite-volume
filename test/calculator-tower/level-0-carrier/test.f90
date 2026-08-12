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
! of this file ARE the negative truth of the level - graph_carrier
! and nothing above it (CALCULATOR.md 7).
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program calculator_level_0

  use calculator_assert, only : report, verdict
  use calculator_assert, only : SLOT_A, SLOT_E, OP_PLUS, PORT_OUT
  use graph_carrier    , only : counted_set

  implicit none

  type(counted_set) :: x, o, p
  integer           :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "calculator tower . level 0 . carriers"
  write(*,'(1x,a)') "============================================="

  ! The three declarations. This is everything Level 0 may do.
  x = counted_set('value-slots', 5)
  o = counted_set('operations' , 2)
  p = counted_set('ports'      , 3)

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

    call report(x % size() .eq. 5, &
         & "X counts its five value slots", nfail)
    call report(o % size() .eq. 2, &
         & "O counts its two operations", nfail)
    call report(p % size() .eq. 3, &
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
    do i = 1, x % size()
       ok = ok .and. (x % member(x % local_index(x % member(i))) &
            &         .eq. x % member(i))
       ok = ok .and. (x % local_index(x % member(i)) .eq. i)
    end do
    call report(ok, "member and local_index invert on X", nfail)

    ok = .true.
    do i = 1, o % size()
       ok = ok .and. (o % member(o % local_index(o % member(i))) &
            &         .eq. o % member(i))
       ok = ok .and. (o % local_index(o % member(i)) .eq. i)
    end do
    call report(ok, "member and local_index invert on O", nfail)

    ok = .true.
    do i = 1, p % size()
       ok = ok .and. (p % member(p % local_index(p % member(i))) &
            &         .eq. p % member(i))
       ok = ok .and. (p % local_index(p % member(i)) .eq. i)
    end do
    call report(ok, "member and local_index invert on P", nfail)

  end subroutine check_enumeration_round_trips

  !===================================================================!
  ! Membership knows its own: a, e and + and out belong where they
  ! were declared; the outsider belongs nowhere.
  !===================================================================!

  subroutine check_membership_boundary(nfail)

    integer, intent(inout) :: nfail

    call report(x % has(SLOT_A) .and. x % has(SLOT_E), &
         & "a and e are value slots", nfail)
    call report(o % has(OP_PLUS) .and. p % has(PORT_OUT), &
         & "+ is an operation and out is a port", nfail)

    call report(.not. x % has(6) .and. .not. x % has(0), &
         & "an outsider is rejected by X", nfail)
    call report(.not. o % has(3), &
         & "and by O", nfail)
    call report(.not. p % has(4), &
         & "and by P", nfail)

  end subroutine check_membership_boundary

end program calculator_level_0
