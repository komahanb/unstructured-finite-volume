!=====================================================================!
! PARTITIONED IMPLICIT PDE TOWER . LEVEL 2 . RELATION ALGEBRA
!
! The level answers one question: WHAT FOLLOWS ALGEBRAICALLY FROM
! THE PRIMITIVE FACTS. Two things follow, and neither is written
! down as data:
!
!      A = Head o Tail^T : V -> V          1->2, 2->3, ... , 5->6
!      EdgeOwner = Own^T o Tail : E -> K   e1,e2,e3 -> part1
!                                          e4,e5    -> part2
!
! The second derivation is the important one, and it is why this
! level exists rather than being marked N/A.
!
! The earlier gate-shaped tower discovered tail-ownership
! OPERATIONALLY: it imposed an assembly law on a probe field and
! then read back what the production partitioner had done. That was
! the right way to find it. But it left the law resting on an
! observation of production.
!
! Here the same rule is DERIVED, from Own and Tail alone, before any
! partitioner exists to be observed:
!
!      an edge belongs to whichever part owns the vertex it LEAVES
!
! so that Level 4 can check production against mathematics instead
! of against a previous reading of production. An edge-ownership
! convention that came from nowhere becomes a theorem with a
! two-line proof.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program partitioned_pde_level_2

  use partitioned_pde_assert , only : report, verdict
  use graph_carrier          , only : counted_set, member_set
  use graph_binary_relation  , only : csr_relation
  use chain_relations_fixture, only : chain_carriers, tail_relation, &
       &                              head_relation, own_relation
  use chain_algebra_fixture  , only : derive_adjacency, derive_edge_owner

  implicit none

  type(counted_set)          :: v, e, k
  type(csr_relation), target :: tail, head, own
  type(csr_relation)         :: adj, edge_owner
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
  edge_owner = derive_edge_owner(tail, own)

  call check_adjacency(nfail)
  call check_edge_owner(nfail)
  call check_edge_owner_is_total_and_single(nfail)

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
  ! THE level's theorem: edge ownership follows from VERTEX
  ! ownership through the tail map. No partitioner has been
  ! consulted; this is what Level 4 will hold production to.
  !===================================================================!

  subroutine check_edge_owner(nfail)

    integer, intent(inout) :: nfail

    class(member_set), allocatable :: d

    d = edge_owner % domain(1)
    call report(d % same_as(e), "EdgeOwner runs from the edges", nfail)
    d = edge_owner % domain(2)
    call report(d % same_as(k), "into the partition labels", nfail)

    call report(edge_owner % num_tuples() .eq. 5 .and. &
         &      edge_owner % has([1, 1]) .and. &
         &      edge_owner % has([2, 1]) .and. &
         &      edge_owner % has([3, 1]) .and. &
         &      edge_owner % has([4, 2]) .and. &
         &      edge_owner % has([5, 2]), &
         & "EdgeOwner = { e1,e2,e3 -> part1 ; e4,e5 -> part2 }", nfail)

    call report(edge_owner % has([3, 1]) .and. &
         &      .not. edge_owner % has([3, 2]), &
         & "the CROSSING edge e3 = 3->4 belongs to part1 - the " // &
         & "owner of its TAIL - and this is derived, not observed", &
         & nfail)

  end subroutine check_edge_owner

  !===================================================================!
  ! And it is a function: every edge has exactly one owner. That is
  ! what will make one global entity yield one assembled
  ! contribution when Level 5 reconstructs a field.
  !===================================================================!

  subroutine check_edge_owner_is_total_and_single(nfail)

    integer, intent(inout) :: nfail

    integer :: i, j, m, owners
    logical :: ok

    ok = .true.
    do i = 1, e % size()
       m = e % member(i)
       owners = 0
       do j = 1, k % size()
          if (edge_owner % has([m, k % member(j)])) owners = owners + 1
       end do
       ok = ok .and. (owners .eq. 1)
    end do
    call report(ok, &
         & "every edge has exactly ONE owner: the reconstruction " // &
         & "law, stated relationally before any field exists", nfail)

  end subroutine check_edge_owner_is_total_and_single

end program partitioned_pde_level_2
