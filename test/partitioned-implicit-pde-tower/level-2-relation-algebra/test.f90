!=====================================================================!
! PARTITIONED IMPLICIT PDE TOWER . LEVEL 2 . RELATION ALGEBRA
!
! The level answers one question: WHAT FOLLOWS ALGEBRAICALLY FROM
! THE PRIMITIVE FACTS. Three things follow, and none is written down
! as data:
!
!      A = Head o Tail^T : V -> V          1->2, 2->3, ... , 5->6
!      TailOwner = Own^T o Tail : E -> K   e1,e2,e3 -> part1
!                                          e4,e5    -> part2
!      HeadOwner = Own^T o Head : E -> K   e1,e2    -> part1
!                                          e3,e4,e5 -> part2
!
! The two ownership maps are why this level exists rather than being
! marked N/A - and why it is careful about a distinction that is
! easy to lose.
!
! DERIVATION is what the algebra does. Name a relational expression
! and its extension follows with no further choices: TailOwner(e3) =
! part1 is a theorem, given TailOwner.
!
! POLICY SELECTION is what the algebra does NOT do. Tail and Head
! are both vertices of an edge, both have owners, and composing
! through either gives a map E -> K that is total and single-valued.
! Both therefore satisfy the reconstruction law - one edge, one
! owner - so uniqueness of the owner does not pick one out. The two
! agree on e1, e2, e4 and e5, and disagree at the crossing edge:
!
!      TailOwner(e3) = part1        HeadOwner(e3) = part2
!
! The earlier gate-shaped tower found tail-ownership OPERATIONALLY -
! it imposed an assembly law on a probe field and read back what the
! production partitioner had done. That was the right way to find
! it. What this level adds is a vocabulary: it states BOTH candidate
! policies as relational mathematics before any partitioner exists,
! so Level 4's reading of production becomes a SELECTION between two
! named alternatives rather than a bare observation.
!
!      RELATION ALGEBRA DERIVES THE CONSEQUENCES OF A POLICY;
!      IT DOES NOT CHOOSE THE POLICY FOR US.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program partitioned_pde_level_2

  use partitioned_pde_assert , only : report, verdict
  use graph_carrier          , only : counted_set, member_set
  use graph_binary_relation  , only : csr_relation
  use chain_carriers_fixture , only : chain_carriers
  use chain_relations_fixture, only : tail_relation, head_relation, &
       &                              own_relation
  use chain_algebra_fixture  , only : derive_adjacency, &
       &                              derive_tail_owner, derive_head_owner

  implicit none

  type(counted_set)          :: v, e, k
  type(csr_relation), target :: tail, head, own
  type(csr_relation)         :: adj, tail_owner, head_owner
  integer                    :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "partitioned pde tower . level 2 . algebra"
  write(*,'(1x,a)') "============================================="

  call chain_carriers(v, e, k)
  tail = tail_relation(e, v)
  head = head_relation(e, v)
  own  = own_relation(k, v)

  adj        = derive_adjacency(tail, head)
  tail_owner = derive_tail_owner(tail, own)
  head_owner = derive_head_owner(head, own)

  call check_adjacency(nfail)
  call check_tail_owner(nfail)
  call check_head_owner(nfail)
  call check_both_are_functions(nfail)
  call check_where_the_policies_part(nfail)

  call verdict(nfail, "level 2")

contains

  !===================================================================!
  ! A = Head o Tail^T : the chain, recovered as vertex-to-vertex
  ! succession without anyone ever writing 1->2.
  !===================================================================!

  subroutine check_adjacency(nfail)

    integer, intent(inout) :: nfail

    class(member_set), allocatable :: d
    integer                        :: i
    logical                        :: ok

    d = adj % domain(1)
    call report(d % same_as(v), "A runs from the vertices", nfail)
    d = adj % domain(2)
    call report(d % same_as(v), "back into the vertices", nfail)

    ok = adj % num_tuples() .eq. 5
    do i = 1, 5
       ok = ok .and. adj % has([i, i + 1])
    end do
    call report(ok, &
         & "A = { 1->2, 2->3, 3->4, 4->5, 5->6 } - derived through " // &
         & "the edges, never written", nfail)

    call report(.not. adj % has([2, 1]) .and. .not. adj % has([1, 3]), &
         & "and it is directed and one-step: no reverse pair, no " // &
         & "shortcut", nfail)

  end subroutine check_adjacency

  !===================================================================!
  ! The FIRST candidate policy: anchor ownership at the vertex the
  ! edge leaves. Its extension is a theorem once the expression is
  ! named; the naming is not.
  !===================================================================!

  subroutine check_tail_owner(nfail)

    integer, intent(inout) :: nfail

    class(member_set), allocatable :: d

    d = tail_owner % domain(1)
    call report(d % same_as(e), "TailOwner runs from the edges", nfail)
    d = tail_owner % domain(2)
    call report(d % same_as(k), "into the partition labels", nfail)

    call report(tail_owner % num_tuples() .eq. 5 .and. &
         &      tail_owner % has([1, 1]) .and. &
         &      tail_owner % has([2, 1]) .and. &
         &      tail_owner % has([3, 1]) .and. &
         &      tail_owner % has([4, 2]) .and. &
         &      tail_owner % has([5, 2]), &
         & "TailOwner = { e1,e2,e3 -> part1 ; e4,e5 -> part2 }", nfail)

  end subroutine check_tail_owner

  !===================================================================!
  ! The SECOND candidate policy, derived by exactly the same two
  ! lines of algebra with Head in place of Tail. Nothing in the
  ! primitive facts prefers one over the other.
  !===================================================================!

  subroutine check_head_owner(nfail)

    integer, intent(inout) :: nfail

    class(member_set), allocatable :: d

    d = head_owner % domain(1)
    call report(d % same_as(e), "HeadOwner runs from the edges too", nfail)
    d = head_owner % domain(2)
    call report(d % same_as(k), "into the same partition labels", nfail)

    call report(head_owner % num_tuples() .eq. 5 .and. &
         &      head_owner % has([1, 1]) .and. &
         &      head_owner % has([2, 1]) .and. &
         &      head_owner % has([3, 2]) .and. &
         &      head_owner % has([4, 2]) .and. &
         &      head_owner % has([5, 2]), &
         & "HeadOwner = { e1,e2 -> part1 ; e3,e4,e5 -> part2 }", nfail)

  end subroutine check_head_owner

  !===================================================================!
  ! THE level's correction to an easy overclaim: BOTH candidates are
  ! total functions E -> K. Each is proved on its own, because each
  ! is its own theorem - and because together they show that "every
  ! edge has exactly one owner" cannot be the thing that selects
  ! tail-ownership. Unique ownership does NOT imply tail ownership.
  !===================================================================!

  subroutine check_both_are_functions(nfail)

    integer, intent(inout) :: nfail

    call report(is_total_function(tail_owner), &
         & "TailOwner is a TOTAL FUNCTION E -> K: every edge one " // &
         & "owner, no edge two", nfail)

    call report(is_total_function(head_owner), &
         & "and so is HeadOwner - so the reconstruction law is " // &
         & "satisfied by BOTH, and cannot choose between them", nfail)

  end subroutine check_both_are_functions

  !===================================================================!
  ! Where the two policies part company, and it is exactly one edge:
  ! the one the cut crosses. Everywhere else the question does not
  ! arise, which is why a specimen without a crossing edge would
  ! have hidden the choice entirely.
  !===================================================================!

  subroutine check_where_the_policies_part(nfail)

    integer, intent(inout) :: nfail

    integer :: ge
    logical :: ok

    ok = .true.
    do ge = 1, 5
       if (ge .eq. 3) cycle
       ok = ok .and. same_owner(ge)
    end do
    call report(ok, &
         & "TailOwner and HeadOwner AGREE on e1, e2, e4 and e5 - " // &
         & "away from the cut the anchor does not matter", nfail)

    call report(.not. same_owner(3), &
         & "and DISAGREE on the crossing edge e3 = 3->4, the only " // &
         & "edge whose two vertices have different owners", nfail)

    call report(tail_owner % has([3, 1]) .and. &
         &      .not. tail_owner % has([3, 2]), &
         & "TailOwner(e3) = part1, the owner of vertex 3", nfail)

    call report(head_owner % has([3, 2]) .and. &
         &      .not. head_owner % has([3, 1]), &
         & "HeadOwner(e3) = part2, the owner of vertex 4 - and the " // &
         & "algebra has no further opinion about which is right", &
         & nfail)

  end subroutine check_where_the_policies_part

  !===================================================================!
  ! A relation is a TOTAL FUNCTION E -> K when it answers the whole
  ! question, and the whole question is four things at once: it runs
  ! between the right two carriers, it holds no fact outside them,
  ! every edge gets an answer, and no edge gets two. The arity check
  ! is what keeps the per-edge count honest - without it a relation
  ! could carry tuples whose first slot is not an edge at all, and
  ! the loop below would never look.
  !===================================================================!

  logical function is_total_function(policy)

    type(csr_relation), intent(in) :: policy

    class(member_set), allocatable :: d
    integer                        :: i, j, m, owners

    d = policy % domain(1)
    is_total_function = d % same_as(e)
    d = policy % domain(2)
    is_total_function = is_total_function .and. d % same_as(k)

    ! No fact lives outside E x K: exactly one tuple per edge, and
    ! the loop below then accounts for every one of them.
    is_total_function = is_total_function .and. &
         & (policy % num_tuples() .eq. e % size())

    do i = 1, e % size()
       m = e % member(i)
       owners = 0
       do j = 1, k % size()
          if (policy % has([m, k % member(j)])) owners = owners + 1
       end do
       is_total_function = is_total_function .and. (owners .eq. 1)
    end do

  end function is_total_function

  logical function same_owner(ge)

    integer, intent(in) :: ge

    integer :: j, m

    same_owner = .true.
    do j = 1, k % size()
       m = k % member(j)
       same_owner = same_owner .and. &
            & (tail_owner % has([ge, m]) .eqv. head_owner % has([ge, m]))
    end do

  end function same_owner

end program partitioned_pde_level_2
