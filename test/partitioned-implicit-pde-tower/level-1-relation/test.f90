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
! carriers - not yet a thing any machinery has realized.
!
! Every signature is pinned by carrier IDENTITY, because the ids
! collide across V, E and K and orientation could not otherwise be
! established: the pair [3,4] is a legitimate member of Tail, Head
! and Own alike, meaning three different things each time.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program partitioned_pde_level_1

  use partitioned_pde_assert , only : report, verdict
  use graph_carrier          , only : counted_set, member_set
  use graph_binary_relation  , only : csr_relation
  use chain_relations_fixture, only : chain_carriers, tail_relation, &
       &                              head_relation, own_relation

  implicit none

  type(counted_set)  :: v, e, k
  type(csr_relation) :: tail, head, own
  integer            :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "partitioned pde tower . level 1 . relation"
  write(*,'(1x,a)') "============================================="

  call chain_carriers(v, e, k)
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
  ! Every slot answers a DECLARED carrier by identity. Sizes could
  ! not do this: E and K differ in size from V, but Tail and Head
  ! share both of theirs exactly.
  !===================================================================!

  subroutine check_signatures(nfail)

    integer, intent(inout) :: nfail

    class(member_set), allocatable :: d

    d = tail % domain(1)
    call report(d % same_as(e), "Tail runs from the edges", nfail)
    d = tail % domain(2)
    call report(d % same_as(v), "into the vertices", nfail)

    d = head % domain(1)
    call report(d % same_as(e), "Head runs from the edges", nfail)
    d = head % domain(2)
    call report(d % same_as(v), "into the vertices as well", nfail)

    d = own % domain(1)
    call report(d % same_as(k), "Own runs from the partition labels", nfail)
    d = own % domain(2)
    call report(d % same_as(v), "into the vertices", nfail)

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
  ! so one and the same pair means three different things depending
  ! only on which relation is asked.
  !===================================================================!

  subroutine check_orientation_is_signature(nfail)

    integer, intent(inout) :: nfail

    call report(tail % has([3, 3]) .and. head % has([3, 4]) .and. &
         &      own % has([1, 3]), &
         & "the pair [3,4] reads as 'e3 enters vertex 4' in Head - " // &
         & "and would read as 'part3 owns vertex 4' in a relation " // &
         & "from K, if K had a third member", nfail)
    call report(.not. own % has([3, 4]), &
         & "Own refuses it outright: K has no member 3, and that " // &
         & "is a carrier truth, not an arithmetic one", nfail)

  end subroutine check_orientation_is_signature

end program partitioned_pde_level_1
