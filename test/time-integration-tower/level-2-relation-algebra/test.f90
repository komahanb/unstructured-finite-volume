!=====================================================================!
! TIME INTEGRATION TOWER . LEVEL 2 . RELATION ALGEBRA
!
! The level answers one question: WHAT TEMPORAL DEPENDENCIES FOLLOW
! FROM THE PRIMITIVE INCIDENCE. Two things follow, and neither is
! written down as data:
!
!      A1 = Head o Tail^T : T -> T    t0->t1 t1->t2 t2->t3 t3->t4
!      A2 = A1 o A1       : T -> T    t0->t2 t1->t3 t2->t4
!
! The chain comes back through the steps without anyone writing
! t0 -> t1, and the two-step reach comes back without anyone writing
! t0 -> t2. That is the level's central statement:
!
!      ONE-STEP AND TWO-STEP TEMPORAL REACH ARE DERIVED,
!      NOT STORED INDEPENDENTLY.
!
!                    WHAT A2 IS NOT
!
! A2 is NOT BDF2. It is not a scheme, not a stencil and not a
! discretization. It says exactly one thing: an instant two steps
! later is structurally REACHABLE.
!
! The test that keeps that honest is the refusal below: A2 does NOT
! relate t0 to t3. A relation that WERE a two-step scheme's
! dependency would have to see the intervening history from t0; the
! two-step reach lands only on even offsets and sees nothing at t3.
! Whatever BDF2 will later consume - present, one-step history,
! two-step history - is an interpretation a Level-6 scheme lays on
! this structure, and no part of the structure itself.
!
!      temporal reach  !=  temporal discretization scheme
!
! No union appears here and no transitive closure. Writing A1 U A2
! merely because the notation exists would be inventing algebra for
! a reader instead of for a caller.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program time_level_2

  use time_assert           , only : report, verdict
  use time_assert           , only : NT
  use time_assert           , only : T0, T1, T2, T3, T4
  use fractal_graph        , only : set_graph => graph
  use graph_set_map        , only : set_map
  use graph_binary_relation , only : csr_relation
  use time_carriers_fixture , only : time_carriers
  use time_relations_fixture, only : tail_relation, head_relation
  use time_algebra_fixture  , only : derive_one_step_reach, &
       &                             derive_two_step_reach

  implicit none

  type(set_graph)          :: q, t, e
  type(set_map)          :: sets
  type(csr_relation), target :: tail, head, a1
  type(csr_relation)         :: a2
  integer                    :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "time integration tower . level 2 . algebra"
  write(*,'(1x,a)') "============================================="

  call time_carriers(sets, q, t, e)
  tail = tail_relation(e, t, sets)
  head = head_relation(e, t, sets)

  a1 = derive_one_step_reach(tail, head, sets)
  a2 = derive_two_step_reach(a1, sets)

  call check_one_step_reach(nfail)
  call check_two_step_reach(nfail)
  call check_the_two_differ(nfail)
  call check_reach_is_not_a_scheme(nfail)

  call verdict(nfail, "level 2")

contains

  !===================================================================!
  ! A1 = Head o Tail^T : the chain of instants, recovered through
  ! the steps without anyone ever writing t0 -> t1.
  !===================================================================!

  subroutine check_one_step_reach(nfail)

    integer, intent(inout) :: nfail

    type(set_graph) :: d
    logical                        :: ok

    d = a1 % domain(1)
    call report(d % same_as(t), "A1 runs from the instants", nfail)
    d = a1 % domain(2)
    call report(d % same_as(t), "back into the instants", nfail)

    ok = a1 % num_tuples() .eq. 4
    ok = ok .and. a1 % has([T0, T1]) .and. a1 % has([T1, T2])
    ok = ok .and. a1 % has([T2, T3]) .and. a1 % has([T3, T4])
    call report(ok, &
         & "A1 = { t0->t1, t1->t2, t2->t3, t3->t4 } - derived " // &
         & "through the steps, never written", nfail)

    call report(.not. a1 % has([T1, T0]) .and. .not. a1 % has([T0, T2]), &
         & "and it is directed and one-step: no reverse pair, no " // &
         & "shortcut", nfail)

  end subroutine check_one_step_reach

  !===================================================================!
  ! A2 = A1 o A1 : two steps of the same reach, and three facts
  ! rather than four - the chain runs out.
  !===================================================================!

  subroutine check_two_step_reach(nfail)

    integer, intent(inout) :: nfail

    type(set_graph) :: d
    logical                        :: ok

    d = a2 % domain(1)
    call report(d % same_as(t), "A2 runs from the instants too", nfail)
    d = a2 % domain(2)
    call report(d % same_as(t), "and back into the same instants", nfail)

    ok = a2 % num_tuples() .eq. 3
    ok = ok .and. a2 % has([T0, T2]) .and. a2 % has([T1, T3])
    ok = ok .and. a2 % has([T2, T4])
    call report(ok, &
         & "A2 = { t0->t2, t1->t3, t2->t4 } - three facts, not " // &
         & "four: two steps do not fit from t3", nfail)

  end subroutine check_two_step_reach

  !===================================================================!
  ! Same signature, same carrier, different relations - by identity
  ! and by extension both. Identity is the address; a shared
  ! signature has never been one.
  !===================================================================!

  subroutine check_the_two_differ(nfail)

    integer, intent(inout) :: nfail

    call report(.not. a1 % same_as(a2), &
         & "A1 and A2 are different relations by IDENTITY, though " // &
         & "they share a signature", nfail)

    call report(a1 % num_tuples() .ne. a2 % num_tuples() .and. &
         &      a1 % has([T0, T1]) .and. .not. a2 % has([T0, T1]) .and. &
         &      a2 % has([T0, T2]) .and. .not. a1 % has([T0, T2]), &
         & "and by EXTENSION: t0 reaches t1 in one step and t2 in " // &
         & "two, and neither relation holds the other's fact", nfail)

  end subroutine check_the_two_differ

  !===================================================================!
  ! THE level's guard against its own most likely misreading. A2 is
  ! reach, and reach skips: from t0 it lands on t2 and t4 only. A
  ! two-step SCHEME would have to consume the history at t1 and t3;
  ! this relation cannot even see t3 from t0.
  !
  ! Whoever later builds BDF2 must therefore CONSTITUTE it out of
  ! A1 and A2 with a scheme's own semantics. They will not find it
  ! lying here.
  !===================================================================!

  subroutine check_reach_is_not_a_scheme(nfail)

    integer, intent(inout) :: nfail

    integer :: i, m, reached
    logical :: ok

    call report(.not. a2 % has([T0, T3]), &
         & "A2 does NOT relate t0 to t3: two-step reach lands on " // &
         & "even offsets and sees no odd one", nfail)

    call report(.not. a2 % has([T0, T1]) .and. .not. a2 % has([T3, T4]), &
         & "nor does it hold any one-step fact - A2 IS NOT BDF2, " // &
         & "and it is not a stencil of anything", nfail)

    ! Every instant reaches at most one instant under either
    ! relation: these are chains, and reach along a chain is a
    ! partial function.
    ok = .true.
    do i = 1, sets % size_of(t)
       m = sets % member_of(t, i)
       reached = count_images(a2, m)
       ok = ok .and. (reached .le. 1)
    end do
    call report(ok, &
         & "and each instant reaches at most one other under A2 - " // &
         & "structure, with no coefficient anywhere in it", nfail)

  end subroutine check_reach_is_not_a_scheme

  integer function count_images(r, m)

    type(csr_relation), intent(in) :: r
    integer           , intent(in) :: m

    integer :: j

    count_images = 0
    do j = 1, sets % size_of(t)
       if (r % has([m, sets % member_of(t, j)])) count_images = count_images + 1
    end do

  end function count_images

end program time_level_2
