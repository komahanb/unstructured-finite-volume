!=====================================================================!
! DERIVATIVE ACTION TOWER . LEVEL 0 . CARRIERS
!
! The level answers one question: WHAT SYMBOLIC KINDS OF MEMBERS
! EXIST. Three independent domains are declared,
!
!      V = { x  y  u  z }            the value slots
!      O = { product  sum }          the operations
!      P = { in1  in2  out }         the ports
!
! and nothing else is true yet. No relation, no graph, no field, no
! derivative, no tangent, no cotangent, no numerical law: product is
! only a member of O here. Nothing about a DERIVATIVE application
! required a new Level-0 ontology - the ordinary carriers suffice,
! and the imports of this file ARE the negative truth:
! derivative_assert, the set modules, and nothing above.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program derivative_level_0

  use derivative_assert, only : report, verdict
  use derivative_assert, only : SLOT_X, SLOT_Z, OP_PRODUCT, PORT_OUT
  use fractal_graph        , only : set_graph => graph
  use graph_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use graph_set_map        , only : set_map

  implicit none

  type(set_graph) :: v, o, p
  integer           :: nfail
  type(set_map)     :: sets

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "derivative action tower . level 0 . carriers"
  write(*,'(1x,a)') "============================================="

  ! The three declarations. This is everything Level 0 may do.
  call v % declare()
  call sets % bind(v, counted_set_representation(4))
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

  subroutine check_cardinalities(nfail)

    integer, intent(inout) :: nfail

    call report(sets % size_of(v) .eq. 4, &
         & "V counts its four value slots", nfail)
    call report(sets % size_of(o) .eq. 2, &
         & "O counts its two operations", nfail)
    call report(sets % size_of(p) .eq. 3, &
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
    do i = 1, sets % size_of(v)
       ok = ok .and. (sets % member_of(v, sets % index_in(v, sets % member_of(v, i))) &
            &         .eq. sets % member_of(v, i))
       ok = ok .and. (sets % index_in(v, sets % member_of(v, i)) .eq. i)
    end do
    call report(ok, "member and local_index invert on V", nfail)

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

  subroutine check_membership_boundary(nfail)

    integer, intent(inout) :: nfail

    call report(sets % has_in(v, SLOT_X) .and. sets % has_in(v, SLOT_Z), &
         & "x and z are value slots", nfail)
    call report(sets % has_in(o, OP_PRODUCT) .and. sets % has_in(p, PORT_OUT), &
         & "product is an operation and out is a port", nfail)

    call report(.not. sets % has_in(v, 5) .and. .not. sets % has_in(v, 0), &
         & "an outsider is rejected by V", nfail)
    call report(.not. sets % has_in(o, 3), &
         & "and by O", nfail)
    call report(.not. sets % has_in(p, 4), &
         & "and by P", nfail)

  end subroutine check_membership_boundary

end program derivative_level_0
