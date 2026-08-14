!=====================================================================!
! TIME INTEGRATION TOWER . LEVEL 1 . RELATION
!
! The level answers one question: WHAT PRIMITIVE FACTS GIVE TIME ITS
! DIRECTION.
!
!      Tail <= E x T     e1->t0  e2->t1  e3->t2  e4->t3
!      Head <= E x T     e1->t1  e2->t2  e3->t3  e4->t4
!
! Two relations, and that is all this level knows. DIRECTION IS
! RELATION STRUCTURE - not the order of a loop index, not the sign
! of a subtraction, and not the fact that 2 is bigger than 1. A step
! knows which end it leaves because Tail says so and Head says
! otherwise.
!
! Q PARTICIPATES IN NO RELATION, and that absence is this level's
! sharpest content rather than an omission. After this rung:
!
!      time has structure
!      state coordinates exist
!      nothing says a state coordinate IS a time instant
!
! Inventing an incidence to attach Q to the time chain would
! manufacture exactly the conflation the tower exists to refuse.
!
! Signatures are pinned by carrier IDENTITY, because the ids collide
! across Q, T and E. And the collision here is sharper than a
! coincidence of small numbers: over this specimen Tail's extension
! is { [1,1] [2,2] [3,3] [4,4] }, which is tuple-for-tuple what a
! six-vertex chain's tail map looks like over entirely different
! carriers. The integers carry none of the meaning. The signature
! carries all of it.
!
! Nothing here knows about q(t), a history field, a derivative, a
! scheme coefficient, Euler, BDF or a marcher.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program time_level_1

  use time_assert           , only : report, verdict
  use time_assert           , only : NT, NE
  use time_assert           , only : T0, T1, T2, T3, T4
  use time_assert           , only : E1, E2, E3, E4
  use graph_carrier         , only : counted_set, member_set
  use graph_binary_relation , only : csr_relation
  use time_carriers_fixture , only : time_carriers
  use time_relations_fixture, only : tail_relation, head_relation

  implicit none

  type(counted_set)  :: q, t, e
  type(csr_relation) :: tail, head
  integer            :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "time integration tower . level 1 . relation"
  write(*,'(1x,a)') "============================================="

  call time_carriers(q, t, e)
  tail = tail_relation(e, t)
  head = head_relation(e, t)

  call check_signatures(nfail)
  call check_incidence_extension(nfail)
  call check_direction_is_structure(nfail)
  call check_state_axis_is_untouched(nfail)

  call verdict(nfail, "level 1")

contains

  !===================================================================!
  ! Every slot answers a DECLARED carrier by identity. Sizes could
  ! not do this alone: Tail and Head share both of theirs exactly.
  !===================================================================!

  subroutine check_signatures(nfail)

    integer, intent(inout) :: nfail

    class(member_set), allocatable :: d

    d = tail % domain(1)
    call report(d % same_as(e), "Tail runs from the steps", nfail)
    d = tail % domain(2)
    call report(d % same_as(t), "into the instants", nfail)

    d = head % domain(1)
    call report(d % same_as(e), "Head runs from the steps", nfail)
    d = head % domain(2)
    call report(d % same_as(t), "into the instants as well", nfail)

  end subroutine check_signatures

  !===================================================================!
  ! The chain, stated as two maps: step e_i leaves instant t_(i-1)
  ! and enters instant t_i.
  !===================================================================!

  subroutine check_incidence_extension(nfail)

    integer, intent(inout) :: nfail

    logical :: ok

    ok = tail % num_tuples() .eq. NE .and. head % num_tuples() .eq. NE
    ok = ok .and. tail % has([E1, T0]) .and. head % has([E1, T1])
    ok = ok .and. tail % has([E2, T1]) .and. head % has([E2, T2])
    ok = ok .and. tail % has([E3, T2]) .and. head % has([E3, T3])
    ok = ok .and. tail % has([E4, T3]) .and. head % has([E4, T4])
    call report(ok, &
         & "Tail and Head hold four facts each: e_i leaves t_(i-1) " // &
         & "and enters t_i", nfail)

    call report(.not. tail % has([E1, T1]) .and. &
         &      .not. head % has([E1, T0]), &
         & "and nothing more - e1 does not leave t1, nor enter t0", &
         & nfail)

  end subroutine check_incidence_extension

  !===================================================================!
  ! THE level's truth. Every step has exactly one tail and exactly
  ! one head, so a step is a directed thing - and the two ends
  ! DISAGREE for every step, which is what makes the direction
  ! real rather than notational. Had Tail and Head agreed anywhere,
  ! that step would be a loop and time would stand still on it.
  !===================================================================!

  subroutine check_direction_is_structure(nfail)

    integer, intent(inout) :: nfail

    integer :: i, m, j, tails, heads
    logical :: ok, disagree

    ok       = .true.
    disagree = .true.
    do i = 1, e % size()
       m = e % member(i)
       tails = 0
       heads = 0
       do j = 1, t % size()
          if (tail % has([m, t % member(j)])) tails = tails + 1
          if (head % has([m, t % member(j)])) heads = heads + 1
       end do
       ok = ok .and. (tails .eq. 1) .and. (heads .eq. 1)

       do j = 1, t % size()
          if (tail % has([m, t % member(j)]) .and. &
              & head % has([m, t % member(j)])) disagree = .false.
       end do
    end do

    call report(ok, &
         & "every step has exactly ONE tail and exactly ONE head", &
         & nfail)
    call report(disagree, &
         & "and no step leaves the instant it enters: DIRECTION IS " // &
         & "RELATION STRUCTURE, not loop index order", nfail)

  end subroutine check_direction_is_structure

  !===================================================================!
  ! The absence that is content. Q was declared at Level 0 and no
  ! fact on this rung mentions it - checked by asking the relations
  ! for their carriers rather than by reading this file's imports.
  !===================================================================!

  subroutine check_state_axis_is_untouched(nfail)

    integer, intent(inout) :: nfail

    class(member_set), allocatable :: d
    integer                        :: k
    logical                        :: mentions_q

    mentions_q = .false.
    do k = 1, tail % arity()
       d = tail % domain(k)
       mentions_q = mentions_q .or. d % same_as(q)
       d = head % domain(k)
       mentions_q = mentions_q .or. d % same_as(q)
    end do

    call report(.not. mentions_q, &
         & "no slot of Tail or Head answers Q: time has structure, " // &
         & "state coordinates exist, and NOTHING says a coordinate " // &
         & "is an instant", nfail)

    call report(q % size() .eq. 2 .and. .not. q % same_as(t), &
         & "Q is still there, still two, still not T - unrelated is " // &
         & "not absent", nfail)

  end subroutine check_state_axis_is_untouched

end program time_level_1
