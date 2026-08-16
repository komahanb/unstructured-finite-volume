!=====================================================================!
! VISUALIZATION TOWER . LEVEL 1 . RELATION
!
! The level answers one question: WHAT PRIMITIVE FACTS CONSTITUTE THE
! THREE OPERATORS' STRUCTURE.
!
! Twelve of them, each an OCCURRENCE with two typed ends:
!
!      T_k  <=  E_k x X_(k-1)      one tail per occurrence
!      H_k  <=  E_k x X_k          one head per occurrence
!
! and every end is checked against the carrier that is entitled to
! hold it - not against an integer, and not against a count.
!
!                      THE SIGNATURES ARE RECTANGULAR
!
! T1 relates E1 to X0; H1 relates E1 to X1; and X0 is not X1. That
! single fact is what separates this specimen from every mesh the
! repository has read so far, where the tail and the head of an edge
! land in ONE vertex carrier. Here the two ends of an occurrence live
! in two different worlds on purpose, because an operator's input and
! output domains are two different worlds.
!
! Level 4 will find that this is exactly the shape the ordinary-graph
! profile cannot read, and Level 1 is where the shape is established.
!
!                       WHAT THIS LEVEL REFUSES TO SAY
!
! There is no D1 here, and no D2, and no D3. A dependency of the form
! a -> p is not a primitive fact about this specimen; it is a
! consequence of two primitive facts about e11, and Level 2 draws
! that consequence. If D1 were written out beside T1 and H1, Level 2
! would be comparing two hand-typed tables and the composition law
! would go untested.
!
! There is also no number anywhere. T and H say WHICH members an
! occurrence touches; by how much is a question Gate A never asks.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program visualization_level_1

  use visualization_assert , only : report, verdict
  use visualization_assert , only : NE1, NE2, NE3
  use visualization_assert , only : X0_A, X0_B, X0_C, X0_D
  use visualization_assert , only : X1_P, X1_Q, X1_R
  use visualization_assert , only : X2_U, X2_V, X2_W
  use visualization_assert , only : X3_M, X3_N
  use fractal_graph        , only : set_graph => graph
  use graph_set_map        , only : set_map
  use graph_label_map      , only : label_map
  use graph_binary_relation, only : csr_relation, binary_relation
  use visualization_carriers_fixture, only : structural_carriers, label_for
  use visualization_relations_fixture, only : occurrences_of_a1
  use visualization_relations_fixture, only : occurrences_of_a2
  use visualization_relations_fixture, only : occurrences_of_a3

  implicit none

  type(set_graph)  :: x0, x1, x2, x3, e1, e2, e3
  type(set_map)     :: sets
  type(label_map)     :: labels
  type(csr_relation) :: t1, h1, t2, h2, t3, h3
  integer            :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "visualization tower . level 1 . relation"
  write(*,'(1x,a)') "============================================="

  call structural_carriers(x0, x1, x2, x3, e1, e2, e3, sets, labels)
  call occurrences_of_a1(e1, x0, x1, t1, h1, sets)
  call occurrences_of_a2(e2, x1, x2, t2, h2, sets)
  call occurrences_of_a3(e3, x2, x3, t3, h3, sets)

  call check_the_signatures_are_typed(nfail)
  call check_the_signatures_are_rectangular(nfail)
  call check_every_occurrence_has_two_ends(nfail)
  call check_the_twelve_occurrences(nfail)
  call check_no_dependency_is_stored(nfail)

  call verdict(nfail, "level 1")

contains

  !===================================================================!
  ! Each of the six relations relates the two carriers it is entitled
  ! to relate, and says so by structural identity rather than by
  ! size - which would not distinguish X1 from X2 at all.
  !===================================================================!

  subroutine check_the_signatures_are_typed(nfail)

    integer, intent(inout) :: nfail

    call report(signed(t1, e1, x0) .and. signed(h1, e1, x1), &
         & "T1 <= E1 x X0 and H1 <= E1 x X1 - A1's occurrences read " // &
         & "from X0 and write to X1", nfail)

    call report(signed(t2, e2, x1) .and. signed(h2, e2, x2), &
         & "T2 <= E2 x X1 and H2 <= E2 x X2", nfail)

    call report(signed(t3, e3, x2) .and. signed(h3, e3, x3), &
         & "T3 <= E3 x X2 and H3 <= E3 x X3", nfail)

    call report(t1 % arity() .eq. 2 .and. h3 % arity() .eq. 2, &
         & "and every primitive relation is binary: an occurrence and " // &
         & "one end of it", nfail)

  end subroutine check_the_signatures_are_typed

  !===================================================================!
  ! THE FACT THAT SHAPES THE WHOLE TOWER. The two ends of an
  ! occurrence do not live in one carrier.
  !===================================================================!

  subroutine check_the_signatures_are_rectangular(nfail)

    integer, intent(inout) :: nfail

    type(set_graph) :: tail_end, head_end

    tail_end = t1 % target()
    head_end = h1 % target()

    call report(.not. tail_end % same_as(head_end), &
         & "A1's TAIL END AND HEAD END ARE DIFFERENT CARRIERS - X0 is " // &
         & "not X1, so this is no mesh edge", nfail)

    tail_end = t2 % target()
    head_end = h2 % target()
    call report(.not. tail_end % same_as(head_end), &
         & "and neither are A2's - X1 is not X2, though both hold " // &
         & "three members", nfail)

    tail_end = t3 % target()
    head_end = h3 % target()
    call report(.not. tail_end % same_as(head_end), &
         & "and neither are A3's - X2 is not X3", nfail)

  end subroutine check_the_signatures_are_rectangular

  !===================================================================!
  ! Exactly one tail and exactly one head per occurrence, on every
  ! occurrence of every operator - asked through the fibre views, so
  ! the answer comes from the relation and not from the table that
  ! built it.
  !===================================================================!

  subroutine check_every_occurrence_has_two_ends(nfail)

    integer, intent(inout) :: nfail

    call report(exactly_one_end_each(t1, h1, e1), &
         & "each of A1's five occurrences has ONE tail and ONE head", &
         & nfail)

    call report(exactly_one_end_each(t2, h2, e2), &
         & "each of A2's four occurrences has ONE tail and ONE head", &
         & nfail)

    call report(exactly_one_end_each(t3, h3, e3), &
         & "each of A3's three occurrences has ONE tail and ONE head", &
         & nfail)

    call report(t1 % num_tuples() .eq. NE1 .and. h1 % num_tuples() .eq. NE1 .and. &
         &      t2 % num_tuples() .eq. NE2 .and. h2 % num_tuples() .eq. NE2 .and. &
         &      t3 % num_tuples() .eq. NE3 .and. h3 % num_tuples() .eq. NE3, &
         & "and no occurrence is missing: 5 + 4 + 3 = TWELVE " // &
         & "OCCURRENCES, twelve tails, twelve heads", nfail)

  end subroutine check_every_occurrence_has_two_ends

  !===================================================================!
  ! The twelve, one at a time, by name. This is the specimen's whole
  ! primitive content; everything above is derived from these lines
  ! and nothing else.
  !===================================================================!

  subroutine check_the_twelve_occurrences(nfail)

    integer, intent(inout) :: nfail

    call report(runs(t1, h1, 1, X0_A, X1_P) .and. &
         &      runs(t1, h1, 2, X0_B, X1_P) .and. &
         &      runs(t1, h1, 3, X0_B, X1_Q) .and. &
         &      runs(t1, h1, 4, X0_C, X1_Q) .and. &
         &      runs(t1, h1, 5, X0_D, X1_R), &
         & "A1 : e11 a->p, e12 b->p, e13 b->q, e14 c->q, e15 d->r", nfail)

    call report(runs(t2, h2, 1, X1_P, X2_U) .and. &
         &      runs(t2, h2, 2, X1_Q, X2_U) .and. &
         &      runs(t2, h2, 3, X1_Q, X2_V) .and. &
         &      runs(t2, h2, 4, X1_R, X2_W), &
         & "A2 : e21 p->u, e22 q->u, e23 q->v, e24 r->w", nfail)

    call report(runs(t3, h3, 1, X2_U, X3_M) .and. &
         &      runs(t3, h3, 2, X2_V, X3_N) .and. &
         &      runs(t3, h3, 3, X2_W, X3_N), &
         & "A3 : e31 u->m, e32 v->n, e33 w->n", nfail)

    ! Two occurrences of A1 write to the same member of X1, and two
    ! read from the same member of X0. Both are ordinary; neither is
    ! a duplicate, because the OCCURRENCES are distinct members.
    call report(t1 % has([2, X0_B]) .and. t1 % has([3, X0_B]) .and. &
         &      h1 % has([1, X1_P]) .and. h1 % has([2, X1_P]), &
         & "e12 and e13 both read b, e11 and e12 both write p - " // &
         & "distinct occurrences may share an end", nfail)

  end subroutine check_the_twelve_occurrences

  !===================================================================!
  ! Nothing at this level says a -> p. Every primitive fact is a fact
  ! ABOUT AN OCCURRENCE: its first slot is always an occurrence
  ! carrier, never a state carrier.
  !===================================================================!

  subroutine check_no_dependency_is_stored(nfail)

    integer, intent(inout) :: nfail

    call report(about_an_occurrence(t1, e1) .and. about_an_occurrence(h1, e1) .and. &
         &      about_an_occurrence(t2, e2) .and. about_an_occurrence(h2, e2) .and. &
         &      about_an_occurrence(t3, e3) .and. about_an_occurrence(h3, e3), &
         & "every one of the six primitive relations is ABOUT AN " // &
         & "OCCURRENCE: E_k stands in its first slot", nfail)

    call report(.not. relates(t1, x0, x1) .and. .not. relates(h1, x0, x1) .and. &
         &      .not. relates(t2, x1, x2) .and. .not. relates(h2, x1, x2) .and. &
         &      .not. relates(t3, x2, x3) .and. .not. relates(h3, x2, x3), &
         & "AND NOT ONE OF THEM IS X_(k-1) x X_k - the dependency D_k " // &
         & "is not stored here, it is Level 2's to derive", nfail)

  end subroutine check_no_dependency_is_stored

  !===================================================================!
  ! Helpers. Every one of them asks the relation, never the table
  ! that built it.
  !===================================================================!

  logical function signed(r, first, second)

    class(binary_relation), intent(in) :: r
    type(set_graph)     , intent(in) :: first, second

    type(set_graph) :: d

    d = r % source()
    signed = d % same_as(first)
    d = r % target()
    signed = signed .and. d % same_as(second)

  end function signed

  logical function relates(r, first, second)

    class(binary_relation), intent(in) :: r
    type(set_graph)     , intent(in) :: first, second

    relates = signed(r, first, second)

  end function relates

  logical function about_an_occurrence(r, occurrences)

    class(binary_relation), intent(in) :: r
    type(set_graph)     , intent(in) :: occurrences

    type(set_graph) :: d

    d = r % source()
    about_an_occurrence = d % same_as(occurrences)

  end function about_an_occurrence

  logical function exactly_one_end_each(tail, head, occurrences)

    class(binary_relation), target, intent(in) :: tail, head
    type(set_graph)             , intent(in) :: occurrences

    integer, pointer :: fibre(:)
    integer          :: k, e

    exactly_one_end_each = .true.
    do k = 1, sets % size_of(occurrences)
       e = sets % member_of(occurrences, k)
       fibre => tail % image_view(e)
       exactly_one_end_each = exactly_one_end_each .and. (size(fibre) .eq. 1)
       fibre => head % image_view(e)
       exactly_one_end_each = exactly_one_end_each .and. (size(fibre) .eq. 1)
    end do

  end function exactly_one_end_each

  logical function runs(tail, head, occurrence, from, to)

    class(binary_relation), intent(in) :: tail, head
    integer               , intent(in) :: occurrence, from, to

    runs = tail % has([occurrence, from]) .and. head % has([occurrence, to])

  end function runs

end program visualization_level_1
