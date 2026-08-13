!=====================================================================!
! LEARNING TOWER . LEVEL 0 . CARRIERS
!
! The level answers one question: WHAT SYMBOLIC KINDS OF MEMBERS
! EXIST. Three independent domains are declared,
!
!      V = { w  x  yhat  y  e }      the value slots
!      O = { predict  error }        the operations
!      P = { in1  in2  out }         the ports
!
! and nothing else is true yet. No relation, no graph, no field, no
! minimizer, no model, no neuron: predict is only a member of O
! here. The imports of this file ARE the negative truth -
! learning_assert, graph_carrier, and nothing above.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program learning_level_0

  use learning_assert, only : report, verdict
  use learning_assert, only : SLOT_W, SLOT_E, OP_PREDICT, PORT_OUT
  use graph_carrier  , only : counted_set

  implicit none

  type(counted_set) :: v, o, p
  integer           :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "learning tower . level 0 . carriers"
  write(*,'(1x,a)') "============================================="

  ! The three declarations. This is everything Level 0 may do.
  v = counted_set('value-slots', 5)
  o = counted_set('operations' , 2)
  p = counted_set('ports'      , 3)

  call check_cardinalities(nfail)
  call check_distinct_identities(nfail)
  call check_enumeration_round_trips(nfail)
  call check_membership_boundary(nfail)

  call verdict(nfail, "level 0")

contains

  subroutine check_cardinalities(nfail)

    integer, intent(inout) :: nfail

    call report(v % size() .eq. 5, &
         & "V counts its five value slots", nfail)
    call report(o % size() .eq. 2, &
         & "O counts its two operations", nfail)
    call report(p % size() .eq. 3, &
         & "P counts its three ports", nfail)

  end subroutine check_cardinalities

  subroutine check_distinct_identities(nfail)

    integer, intent(inout) :: nfail

    call report(.not. v % same_as(o), &
         & "V is not O", nfail)
    call report(.not. v % same_as(p), &
         & "V is not P", nfail)
    call report(.not. o % same_as(p), &
         & "O is not P", nfail)

  end subroutine check_distinct_identities

  !===================================================================!
  ! The two enumeration laws on every member of every carrier:
  ! member(local_index(m)) = m and local_index(member(i)) = i.
  !===================================================================!

  subroutine check_enumeration_round_trips(nfail)

    integer, intent(inout) :: nfail

    logical :: ok
    integer :: i

    ok = .true.
    do i = 1, v % size()
       ok = ok .and. (v % member(v % local_index(v % member(i))) &
            &         .eq. v % member(i))
       ok = ok .and. (v % local_index(v % member(i)) .eq. i)
    end do
    call report(ok, "member and local_index invert on V", nfail)

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

  subroutine check_membership_boundary(nfail)

    integer, intent(inout) :: nfail

    call report(v % has(SLOT_W) .and. v % has(SLOT_E), &
         & "w and e are value slots", nfail)
    call report(o % has(OP_PREDICT) .and. p % has(PORT_OUT), &
         & "predict is an operation and out is a port", nfail)

    call report(.not. v % has(6) .and. .not. v % has(0), &
         & "an outsider is rejected by V", nfail)
    call report(.not. o % has(3), &
         & "and by O", nfail)
    call report(.not. p % has(4), &
         & "and by P", nfail)

  end subroutine check_membership_boundary

end program learning_level_0
