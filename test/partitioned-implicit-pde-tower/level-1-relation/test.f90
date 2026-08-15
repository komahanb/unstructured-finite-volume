!=====================================================================!
! PARTITIONED IMPLICIT PDE TOWER . LEVEL 1 . RELATION
!
! The level answers one question: WHAT PRIMITIVE FACTS DESCRIBE THE
! CHAIN AND ITS INTENDED OWNERSHIP, before either is interpreted as
! a graph or a partition.
!
!      Tail <= E x V     e_i -> i
!      Head <= E x V     e_i -> i+1
!      Own  <= K x V     part1 -> 1,2,3    part2 -> 4,5,6
!
! Three relations, and that is all this level knows. There is no
! overlap, no borrowed member, no local numbering, no graph object
! and no partitioner. Ownership here is an INTENTION stated over
! sets - not yet a thing any machinery has realized.
!
! Every signature is pinned by set IDENTITY, because the ids
! collide across V, E and K and orientation could not otherwise be
! established. The collisions are real and specific:
!
!      [1,1]  is a fact of Tail  -  e1 leaves vertex 1
!             and a fact of Own  -  part1 owns vertex 1
!
!      [1,2]  is a fact of Head  -  e1 enters vertex 2
!             and a fact of Own  -  part1 owns vertex 2
!
! Identical integers, unrelated meanings, and nothing in the pair
! itself says which is which. Only the SIGNATURE does - E x V in one
! case, K x V in the other:
!
!      raw integer tuple  /=  typed relational fact.
!
! No pair belongs to all three relations, and none could: Tail and
! Head are disjoint on this chain, and a K-indexed relation admits
! only 1 and 2 in its first slot. That last point is a set
! truth, not an arithmetic one - Own refuses [3,4] outright because
! K has no member 3.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program partitioned_pde_level_1

  use partitioned_pde_assert , only : report, verdict
  use graph_set          , only : index_set, set
  use graph_binary_relation  , only : csr_relation
  use chain_sets_fixture , only : chain_sets
  use chain_relations_fixture, only : tail_relation, head_relation, &
       &                              own_relation

  implicit none

  type(index_set)  :: v, e, k
  type(csr_relation) :: tail, head, own
  integer            :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "partitioned pde tower . level 1 . relation"
  write(*,'(1x,a)') "============================================="

  call chain_sets(v, e, k)
  tail = tail_relation(e, v)
  head = head_relation(e, v)
  own  = own_relation(k, v)

  call check_signatures(nfail)
  call check_incidence_extension(nfail)
  call check_ownership_extension(nfail)
  call check_orientation_is_signature(nfail)

  call verdict(nfail, "level 1")

contains

  !===================================================================!
  ! Every slot answers a DECLARED set by identity. Sizes could
  ! not do this: E and K differ in size from V, but Tail and Head
  ! share both of theirs exactly.
  !===================================================================!

  subroutine check_signatures(nfail)

    integer, intent(inout) :: nfail

    class(set), allocatable :: d

    d = tail % domain(1)
    call report(d % equals(e), "Tail runs from the edges", nfail)
    d = tail % domain(2)
    call report(d % equals(v), "into the vertices", nfail)

    d = head % domain(1)
    call report(d % equals(e), "Head runs from the edges", nfail)
    d = head % domain(2)
    call report(d % equals(v), "into the vertices as well", nfail)

    d = own % domain(1)
    call report(d % equals(k), "Own runs from the partition labels", nfail)
    d = own % domain(2)
    call report(d % equals(v), "into the vertices", nfail)

  end subroutine check_signatures

  !===================================================================!
  ! The chain, stated as two maps: edge i leaves vertex i and enters
  ! vertex i+1.
  !===================================================================!

  subroutine check_incidence_extension(nfail)

    integer, intent(inout) :: nfail

    integer :: i
    logical :: ok

    ok = tail % num_tuples() .eq. 5 .and. head % num_tuples() .eq. 5
    do i = 1, 5
       ok = ok .and. tail % has([i, i])
       ok = ok .and. head % has([i, i + 1])
    end do
    call report(ok, &
         & "Tail and Head hold five facts each: e_i leaves i and " // &
         & "enters i+1", nfail)

    call report(.not. tail % has([1, 2]) .and. .not. head % has([1, 1]), &
         & "and nothing more - e1 does not leave vertex 2, nor " // &
         & "enter vertex 1", nfail)

  end subroutine check_incidence_extension

  !===================================================================!
  ! The INTENDED ownership: a total map from vertices to parts,
  ! stated before anything exists to realize it.
  !===================================================================!

  subroutine check_ownership_extension(nfail)

    integer, intent(inout) :: nfail

    integer :: i, count_owners
    logical :: ok

    call report(own % num_tuples() .eq. 6 .and. &
         &      own % has([1, 1]) .and. own % has([1, 2]) .and. &
         &      own % has([1, 3]) .and. own % has([2, 4]) .and. &
         &      own % has([2, 5]) .and. own % has([2, 6]), &
         & "Own = { part1 -> 1,2,3 ; part2 -> 4,5,6 }", nfail)

    ! Every vertex has exactly one owner - the law that will make
    ! reconstruction possible three levels above this one.
    ok = .true.
    do i = 1, v % size()
       count_owners = 0
       if (own % has([1, v % member(i)])) count_owners = count_owners + 1
       if (own % has([2, v % member(i)])) count_owners = count_owners + 1
       ok = ok .and. (count_owners .eq. 1)
    end do
    call report(ok, &
         & "and every vertex has exactly ONE owner: no vertex is " // &
         & "shared, none is orphaned", nfail)

  end subroutine check_ownership_extension

  !===================================================================!
  ! Why signatures matter here more than anywhere: the ids collide,
  ! so one and the same pair of integers is a fact of two different
  ! relations meaning two unrelated things.
  !===================================================================!

  subroutine check_orientation_is_signature(nfail)

    integer, intent(inout) :: nfail

    call report(tail % has([1, 1]) .and. own % has([1, 1]), &
         & "the pair [1,1] is 'e1 leaves vertex 1' in Tail over " // &
         & "E x V, and 'part1 owns vertex 1' in Own over K x V", &
         & nfail)

    call report(head % has([1, 2]) .and. own % has([1, 2]), &
         & "and [1,2] is 'e1 enters vertex 2' in Head, and 'part1 " // &
         & "owns vertex 2' in Own: the signature is the semantics", &
         & nfail)

    call report(head % has([3, 4]) .and. .not. tail % has([3, 4]) .and. &
         &      .not. own % has([3, 4]), &
         & "no pair is a fact of all three - [3,4] is Head's alone: " // &
         & "Tail says e3 leaves vertex 3, and Own refuses it because " // &
         & "K has no member 3, a set truth, not an arithmetic one", &
         & nfail)

  end subroutine check_orientation_is_signature

end program partitioned_pde_level_1
