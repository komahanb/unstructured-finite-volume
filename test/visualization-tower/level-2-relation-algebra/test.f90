!=====================================================================!
! VISUALIZATION TOWER . LEVEL 2 . RELATION ALGEBRA
!
! The level answers one question: WHAT DEPENDENCY STRUCTURE FOLLOWS
! ALGEBRAICALLY.
!
! From twelve primitive occurrences and nothing else:
!
!      D_k    = H_k o T_k^T                  each operator's own
!                                            dependency
!      D_2:1  = D_2 o D_1     : X0 -> X2     two operators' worth
!      D_3:1  = D_3 o D_2 o D_1 : X0 -> X3   the whole chain's
!                                            structural skeleton
!
! and then, reversing every arrow:
!
!      D_k^T                                 reverse dependency
!      D_rev  = D_1^T o D_2^T o D_3^T : X3 -> X0
!
! with the law this tower exists to state:
!
!      (D_3 o D_2 o D_1)^T  =  D_1^T o D_2^T o D_3^T.
!
!                        THE WITNESS COLLAPSE
!
! Seven walks run from X0 to X2 through X1, and D_2:1 holds six
! tuples. Both b->p->u and b->q->u witness the dependency b->u, and
! a relation is a SET: two witnesses of one fact are one fact.
!
! This is not an optimization and not a de-duplication convenience.
! It is what makes D_2:1 the structure of the composed operator
! rather than an accounting of paths through it - and it is the
! nucleus's constructor that enforces it, so the level checks that
! the nucleus did, rather than doing it itself.
!
!                   STRUCTURAL TRANSPOSE, NOT ADJOINT
!
! D^T reverses the ORIENTATION OF POSSIBLE DEPENDENCE: it says which
! members of X0 could have reached a given member of X3. That is all
! it says. There is no coefficient to transpose, no vector to apply
! it to, no inner product under which anything is dual, and no
! sensitivity being propagated. The single adjoint-shaped statement
! Gate A makes is
!
!      supp(A^T)  =  supp(A)^T
!
! and even that is a statement about supports, made in a tower where
! no A has been given entries.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program visualization_level_2

  use visualization_assert , only : report, verdict
  use visualization_assert , only : ND1, ND2, ND3, ND21, ND31
  use visualization_assert , only : X0_A, X0_B, X0_C, X0_D
  use visualization_assert , only : X1_P, X1_Q, X1_R
  use visualization_assert , only : X2_U, X2_V, X2_W
  use visualization_assert , only : X3_M, X3_N
  use graph_set        , only : index_set, set
  use graph_relation       , only : relation
  use graph_binary_relation, only : csr_relation, binary_relation
  use graph_binary_relation, only : transposed_view, transpose_of
  use visualization_sets_fixture , only : structural_sets
  use visualization_relations_fixture, only : occurrences_of_a1
  use visualization_relations_fixture, only : occurrences_of_a2
  use visualization_relations_fixture, only : occurrences_of_a3
  use visualization_algebra_fixture  , only : derive_dependency
  use visualization_algebra_fixture  , only : derive_composition
  use visualization_algebra_fixture  , only : materialized_transpose
  use visualization_algebra_fixture  , only : named, same_extension

  implicit none

  type(index_set)          :: x0, x1, x2, x3, e1, e2, e3
  type(csr_relation)         :: t1, h1, t2, h2, t3, h3
  type(csr_relation), target :: d1, d2, d3, d21, d31
  type(csr_relation), target :: d1t, d2t, d3t
  type(csr_relation)         :: middle, drev, d31t
  integer                    :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "visualization tower . level 2 . relation algebra"
  write(*,'(1x,a)') "============================================="

  call structural_sets(x0, x1, x2, x3, e1, e2, e3)
  call occurrences_of_a1(e1, x0, x1, t1, h1)
  call occurrences_of_a2(e2, x1, x2, t2, h2)
  call occurrences_of_a3(e3, x2, x3, t3, h3)

  ! ---- FORWARD. D_k = H_k o T_k^T, spelled compose_binary(T_k^T, H_k).
  d1 = derive_dependency('D1', t1, h1)
  d2 = derive_dependency('D2', t2, h2)
  d3 = derive_dependency('D3', t3, h3)

  ! ---- and along the chain. D_2:1 = D_2 o D_1, then D_3 after that.
  d21 = derive_composition('D2:1', d1,  d2)
  d31 = derive_composition('D3:1', d21, d3)

  ! ---- REVERSE. Every arrow turned, then composed the other way.
  d1t = materialized_transpose('D1^T', d1)
  d2t = materialized_transpose('D2^T', d2)
  d3t = materialized_transpose('D3^T', d3)

  middle = derive_composition('D2^T o D3^T', d3t, d2t)
  drev   = derive_composition('Drev'       , middle, d1t)

  d31t = materialized_transpose('D3:1^T', d31)

  call check_the_three_dependencies(nfail)
  call check_the_intermediate_composition(nfail)
  call check_the_witness_collapse(nfail)
  call check_the_full_composition(nfail)
  call check_the_reverse_dependencies(nfail)
  call check_the_transpose_composition_law(nfail)
  call check_naming_changes_nothing(nfail)

  call verdict(nfail, "level 2")

contains

  !===================================================================!
  ! D1, D2, D3 - derived, never written down. Each lands on the two
  ! sets its operator relates, and holds exactly the tuples the
  ! occurrences imply.
  !===================================================================!

  subroutine check_the_three_dependencies(nfail)

    integer, intent(inout) :: nfail

    call report(runs_from_to(d1, x0, x1) .and. &
         &      holds_exactly(d1, reshape([X0_A, X1_P, X0_B, X1_P, &
         &                                 X0_B, X1_Q, X0_C, X1_Q, &
         &                                 X0_D, X1_R], [2, ND1])), &
         & "D1 : X0 -> X1 = { a->p, b->p, b->q, c->q, d->r } - " // &
         & "DERIVED from E1's occurrences, not stored", nfail)

    call report(runs_from_to(d2, x1, x2) .and. &
         &      holds_exactly(d2, reshape([X1_P, X2_U, X1_Q, X2_U, &
         &                                 X1_Q, X2_V, X1_R, X2_W], &
         &                                [2, ND2])), &
         & "D2 : X1 -> X2 = { p->u, q->u, q->v, r->w }", nfail)

    call report(runs_from_to(d3, x2, x3) .and. &
         &      holds_exactly(d3, reshape([X2_U, X3_M, X2_V, X3_N, &
         &                                 X2_W, X3_N], [2, ND3])), &
         & "D3 : X2 -> X3 = { u->m, v->n, w->n }", nfail)

    call report(d1 % num_tuples() .eq. ND1 .and. &
         &      d2 % num_tuples() .eq. ND2 .and. &
         &      d3 % num_tuples() .eq. ND3, &
         & "5, 4 and 3 tuples - one per occurrence, because no two " // &
         & "occurrences of one operator share BOTH ends", nfail)

  end subroutine check_the_three_dependencies

  !===================================================================!
  ! D_2:1 = D_2 o D_1, and the seven walks that make its six tuples.
  !===================================================================!

  subroutine check_the_intermediate_composition(nfail)

    integer, intent(inout) :: nfail

    call report(runs_from_to(d21, x0, x2) .and. &
         &      holds_exactly(d21, reshape([X0_A, X2_U, X0_B, X2_U, &
         &                                  X0_B, X2_V, X0_C, X2_U, &
         &                                  X0_C, X2_V, X0_D, X2_W], &
         &                                 [2, ND21])), &
         & "D2:1 : X0 -> X2 = { a->u, b->u, b->v, c->u, c->v, d->w }", &
         & nfail)

    call report(d21 % num_tuples() .eq. ND21 .and. walks(d1, d2) .eq. 7, &
         & "SEVEN WALKS, SIX TUPLES - the composed structure is not " // &
         & "an accounting of paths", nfail)

    ! X1 has been forgotten. It was the middle domain of the
    ! composition, and composition is what forgets it.
    call report(.not. mentions(d21, x1), &
         & "and X1 appears in neither slot of D2:1 - the intermediate " // &
         & "set is spent by the composition", nfail)

  end subroutine check_the_intermediate_composition

  !===================================================================!
  ! The collapse, named member by named member. b reaches u through
  ! p AND through q; the dependency b->u is one fact.
  !===================================================================!

  subroutine check_the_witness_collapse(nfail)

    integer, intent(inout) :: nfail

    call report(witnesses(d1, d2, X0_B, X2_U) .eq. 2 .and. &
         &      d21 % has([X0_B, X2_U]), &
         & "b->p->u and b->q->u are TWO WITNESSES of b->u", nfail)

    call report(count_tuples(d21, X0_B, X2_U) .eq. 1, &
         & "and D2:1 holds b->u ONCE: two witnesses do not make two " // &
         & "dependency tuples", nfail)

    call report(witnesses(d1, d2, X0_A, X2_U) .eq. 1 .and. &
         &      witnesses(d1, d2, X0_D, X2_W) .eq. 1, &
         & "a->u and d->w have one witness each, and are also held " // &
         & "once - the law is set semantics, not de-duplication", nfail)

    call report(witnesses(d1, d2, X0_A, X2_V) .eq. 0 .and. &
         &      .not. d21 % has([X0_A, X2_V]), &
         & "a->v has NO witness and is absent: composition adds " // &
         & "nothing it cannot walk to", nfail)

  end subroutine check_the_witness_collapse

  !===================================================================!
  ! D_3:1 - the structural skeleton of the complete operator chain.
  ! Boolean composition of relations; no matrix has been multiplied,
  ! because no matrix has entries.
  !===================================================================!

  subroutine check_the_full_composition(nfail)

    integer, intent(inout) :: nfail

    type(csr_relation) :: other_way

    call report(runs_from_to(d31, x0, x3) .and. &
         &      holds_exactly(d31, reshape([X0_A, X3_M, X0_B, X3_M, &
         &                                  X0_B, X3_N, X0_C, X3_M, &
         &                                  X0_C, X3_N, X0_D, X3_N], &
         &                                 [2, ND31])), &
         & "D3:1 : X0 -> X3 = { a->m, b->m, b->n, c->m, c->n, d->n } " // &
         & "- THE CHAIN'S SKELETON", nfail)

    ! Composition is associative, and the tower may as well check it
    ! rather than assume the order it happened to bracket in.
    other_way = derive_composition('D3:1 the other way', &
         &                         d1, derive_composition('D3 o D2', d2, d3))
    call report(same_extension(d31, other_way), &
         & "(D3 o D2) o D1 = D3 o (D2 o D1) - the bracketing of the " // &
         & "chain is not part of its structure", nfail)

    call report(.not. mentions(d31, x1) .and. .not. mentions(d31, x2), &
         & "and both intermediate sets are spent: D3:1 relates " // &
         & "only where the chain starts and where it ends", nfail)

  end subroutine check_the_full_composition

  !===================================================================!
  ! The reverse chain: X3 -> X2 -> X1 -> X0, each arrow a transpose
  ! of the forward one. Orientation is the whole content.
  !===================================================================!

  subroutine check_the_reverse_dependencies(nfail)

    integer, intent(inout) :: nfail

    type(transposed_view) :: live

    call report(runs_from_to(d1t, x1, x0) .and. &
         &      runs_from_to(d2t, x2, x1) .and. &
         &      runs_from_to(d3t, x3, x2), &
         & "D1^T : X1 -> X0, D2^T : X2 -> X1, D3^T : X3 -> X2 - " // &
         & "THE REVERSE STRUCTURAL CHAIN", nfail)

    ! Read every tuple in the signature it belongs to. The obvious
    ! spelling of this check - "D1^T does not hold a->p" - cannot be
    ! written at all, because a and p are both the integer 1 and the
    ! pair is its own mirror. Level 0 said the integers are not the
    ! members; here is what that costs when it is forgotten.
    call report(turned_around(d1, d1t), &
         & "every tuple x->y of D1 appears in D1^T as y->x - THE " // &
         & "ARROW REALLY IS TURNED, on all five", nfail)

    call report(d1 % has([X0_C, X1_Q]) .and. d1t % has([X1_Q, X0_C]) .and. &
         &      .not. d1t % has([X1_R, X0_B]), &
         & "c->q becomes q->c, while r->b - which D1 never held as " // &
         & "b->r - is absent from D1^T too", nfail)

    call report(d1t % num_tuples() .eq. ND1 .and. &
         &      d2t % num_tuples() .eq. ND2 .and. &
         &      d3t % num_tuples() .eq. ND3, &
         & "turning the arrows changes no count - transpose moves " // &
         & "orientation, not content", nfail)

    ! The materialized copy against the nucleus's own live view.
    live = transpose_of(d1)
    call report(same_extension(d1t, live), &
         & "and the owned transpose equals transpose_of(D1) exactly - " // &
         & "materializing a view adds ownership and nothing else", nfail)

    call report(same_extension(materialized_transpose('back again', d1t), d1), &
         & "(D1^T)^T = D1 extensionally - an involution about " // &
         & "content, never about identity stamps", nfail)

  end subroutine check_the_reverse_dependencies

  !===================================================================!
  ! THE LAW THIS LEVEL EXISTS FOR.
  !
  !      (D3 o D2 o D1)^T  =  D1^T o D2^T o D3^T
  !
  ! compared as extensions - domains by equals, count, and
  ! membership both ways - never as a list of tuples in some order.
  !===================================================================!

  subroutine check_the_transpose_composition_law(nfail)

    integer, intent(inout) :: nfail

    call report(runs_from_to(drev, x3, x0) .and. runs_from_to(d31t, x3, x0), &
         & "both reverse compositions run X3 -> X0", nfail)

    call report(same_extension(drev, d31t), &
         & "(D3 o D2 o D1)^T = D1^T o D2^T o D3^T - COMPOSE THEN " // &
         & "REVERSE IS REVERSE THEN COMPOSE, BACKWARDS", nfail)

    call report(drev % num_tuples() .eq. ND31 .and. &
         &      drev % has([X3_M, X0_A]) .and. drev % has([X3_N, X0_D]) .and. &
         &      .not. drev % has([X3_M, X0_D]) .and. &
         &      .not. drev % has([X3_N, X0_A]), &
         & "m is reached from a, b, c; n from b, c, d - and d never " // &
         & "reaches m", nfail)

    ! Orientation is part of equality. D3:1 and its transpose hold
    ! mirror-image tuples and are NOT the same relation.
    call report(.not. same_extension(d31, d31t), &
         & "and D3:1 is NOT equal to D3:1^T - a relation over X0 x X3 " // &
         & "is not one over X3 x X0, however its tuples read", nfail)

  end subroutine check_the_transpose_composition_law

  !===================================================================!
  ! The relabelling that must change nothing. Every derived relation
  ! above passed through `named`; here is the proof that it came out
  ! the other side with the same extension it went in with.
  !===================================================================!

  subroutine check_naming_changes_nothing(nfail)

    integer, intent(inout) :: nfail

    call report(same_extension(named('call it something else', d31), d31) .and. &
         &      same_extension(named('or this', d1), d1), &
         & "naming a relation changes its label and nothing else", nfail)

    call report(d1 % name() .eq. 'D1' .and. d31 % name() .eq. 'D3:1' .and. &
         &      d3t % name() .eq. 'D3^T', &
         & "and the derived relations carry the names the pictures " // &
         & "will print - metadata, carried by the object", nfail)

  end subroutine check_naming_changes_nothing

  !===================================================================!
  ! Helpers.
  !===================================================================!

  logical function runs_from_to(r, first, second)

    class(binary_relation), intent(in) :: r
    class(set)     , intent(in) :: first, second

    class(set), allocatable :: d

    d = r % source()
    runs_from_to = d % equals(first)
    d = r % target()
    runs_from_to = runs_from_to .and. d % equals(second)

  end function runs_from_to

  !-------------------------------------------------------------------!
  ! Every tuple of the forward relation, present reversed in the
  ! transposed one - the whole extension, not a sample.
  !-------------------------------------------------------------------!

  logical function turned_around(forward, reverse)

    class(relation), intent(in) :: forward, reverse

    integer, allocatable :: table(:,:)
    integer              :: j

    call forward % tuples(table)
    turned_around = (forward % num_tuples() .eq. reverse % num_tuples())
    do j = 1, size(table, 2)
       turned_around = turned_around .and. reverse % has([table(2, j), table(1, j)])
    end do

  end function turned_around

  logical function mentions(r, set)

    class(relation)  , intent(in) :: r
    class(set), intent(in) :: set

    class(set), allocatable :: d
    integer                        :: k

    mentions = .false.
    do k = 1, r % arity()
       d = r % domain(k)
       mentions = mentions .or. d % equals(set)
    end do

  end function mentions

  !-------------------------------------------------------------------!
  ! Exactness: every listed tuple is held, and nothing else is.
  !-------------------------------------------------------------------!

  logical function holds_exactly(r, table)

    class(relation), intent(in) :: r
    integer        , intent(in) :: table(:,:)

    integer :: j

    holds_exactly = (r % num_tuples() .eq. size(table, 2))
    do j = 1, size(table, 2)
       holds_exactly = holds_exactly .and. r % has(table(:, j))
    end do

  end function holds_exactly

  !-------------------------------------------------------------------!
  ! How many intermediates carry x to z, and how many walks there are
  ! in all. Computed from the two dependencies directly, so the
  ! composition never gets to answer a question about itself.
  !-------------------------------------------------------------------!

  integer function witnesses(first, second, from, to)

    class(binary_relation), intent(in) :: first, second
    integer               , intent(in) :: from, to

    class(set), allocatable :: middle
    integer                        :: k, y

    middle    = first % target()
    witnesses = 0
    do k = 1, middle % size()
       y = middle % member(k)
       if (first % has([from, y]) .and. second % has([y, to])) then
          witnesses = witnesses + 1
       end if
    end do

  end function witnesses

  integer function walks(first, second)

    class(binary_relation), intent(in) :: first, second

    class(set), allocatable :: start, finish
    integer                        :: i, j

    start  = first % source()
    finish = second % target()
    walks  = 0
    do i = 1, start % size()
       do j = 1, finish % size()
          walks = walks + witnesses(first, second, &
               &                    start % member(i), finish % member(j))
       end do
    end do

  end function walks

  !-------------------------------------------------------------------!
  ! How many times a relation's own tuple table lists one pair. One,
  ! if the set semantics hold.
  !-------------------------------------------------------------------!

  integer function count_tuples(r, from, to)

    class(relation), intent(in) :: r
    integer        , intent(in) :: from, to

    integer, allocatable :: table(:,:)
    integer              :: j

    call r % tuples(table)
    count_tuples = 0
    do j = 1, size(table, 2)
       if (table(1, j) .eq. from .and. table(2, j) .eq. to) then
          count_tuples = count_tuples + 1
       end if
    end do

  end function count_tuples

end program visualization_level_2
