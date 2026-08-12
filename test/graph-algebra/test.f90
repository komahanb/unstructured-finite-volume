!=====================================================================!
! The relation algebra suite: the laws of the three earned
! primitives (AGENTS.md 9, level 2) - restriction, projection,
! binary composition - pinned generically, apart from the
! calculator that earned them. Nothing unearned is tested here:
! identity and associativity wait for the callers that need chains.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_graph_algebra

  use graph_carrier        , only : counted_set, subset_set, member_set
  use graph_relation       , only : stored_relation, relation
  use graph_relation_algebra, only : restrict_slot, project_slots, &
       &                             compose_binary

  implicit none

  integer :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "graph relation algebra suite (level 2)"
  write(*,'(1x,a)') "============================================="

  call check_restriction(nfail)
  call check_projection(nfail)
  call check_composition(nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all algebra checks passed"
  else
     write(*,'(1x,a,i0,a)') "FAILED: ", nfail, " algebra check(s)"
     error stop
  end if

contains

  subroutine report(ok, label, nfail)

    logical         , intent(in)    :: ok
    character(len=*), intent(in)    :: label
    integer         , intent(inout) :: nfail

    if (ok) then
       write(*,'(1x,a,a)') "PASS : ", label
    else
       write(*,'(1x,a,a)') "FAIL : ", label
       nfail = nfail + 1
    end if

  end subroutine report

  !===================================================================!
  ! Restriction keeps exactly the admitted tuples, holds the
  ! signature unchanged, and restricting by the full domain itself
  ! is the lawful identity restriction (A embeds in A).
  !===================================================================!

  subroutine check_restriction(nfail)

    integer, intent(inout) :: nfail

    type(counted_set)              :: a, b, c
    type(subset_set)               :: some_b, nobody
    type(stored_relation)          :: r, narrowed
    class(member_set), allocatable :: d
    integer, allocatable           :: rt(:,:)
    integer                        :: j
    logical                        :: ok

    a = counted_set('a-things', 3)
    b = counted_set('b-things', 4)
    c = counted_set('c-things', 2)

    r = stored_relation('r', [a, b, c], &
         & reshape([1,1,1,  1,2,2,  2,2,1,  3,4,2], [3, 4]))

    some_b = subset_set('some-b', b, [2])

    narrowed = restrict_slot(r, 2, some_b)

    call report(narrowed % num_tuples() .eq. 2, &
         & "restriction keeps exactly the admitted tuples", nfail)
    call report(narrowed % has([1, 2, 2]) .and. narrowed % has([2, 2, 1]), &
         & "and they are the right ones", nfail)
    call report(narrowed % arity() .eq. 3, &
         & "the signature's arity stands unchanged", nfail)
    d = narrowed % domain(2)
    call report(d % same_as(b), &
         & "and the restricted slot still answers its full domain", nfail)

    ! Full-domain restriction is the identity, extensionally: equal
    ! count and every original tuple present - for two sets of equal
    ! finite size, that is equality.
    narrowed = restrict_slot(r, 2, b)
    call r % tuples(rt)
    ok = narrowed % num_tuples() .eq. r % num_tuples()
    do j = 1, size(rt, 2)
       ok = ok .and. narrowed % has(rt(:, j))
    end do
    call report(ok, &
         & "restricting by the full domain is the identity, as sets", nfail)

    ! The empty subset admits nothing; the signature stands whole.
    nobody   = subset_set('nobody', b, [integer ::])
    narrowed = restrict_slot(r, 2, nobody)
    call report(narrowed % num_tuples() .eq. 0, &
         & "restriction by the empty subset is the empty relation", nfail)
    call report(narrowed % arity() .eq. 3, &
         & "whose arity is the original's", nfail)
    d = narrowed % domain(1)
    ok = d % same_as(a)
    d = narrowed % domain(2)
    ok = ok .and. d % same_as(b)
    d = narrowed % domain(3)
    ok = ok .and. d % same_as(c)
    call report(ok, &
         & "and whose signature is the original's, slot for slot", nfail)

  end subroutine check_restriction

  !===================================================================!
  ! Projection answers exactly the chosen slots in the chosen
  ! order - [2,1] is the reversed signature - and its image is a
  ! set: tuples that collapse, collapse.
  !===================================================================!

  subroutine check_projection(nfail)

    integer, intent(inout) :: nfail

    type(counted_set)              :: a, b, c
    type(stored_relation)          :: r, image, none
    class(member_set), allocatable :: d

    a = counted_set('a-things', 3)
    b = counted_set('b-things', 4)
    c = counted_set('c-things', 2)

    ! Two tuples agree on (slot1, slot2); they differ only in slot 3.
    r = stored_relation('r', [a, b, c], &
         & reshape([1,1,1,  1,1,2,  2,3,1], [3, 3]))

    image = project_slots(r, [1, 2])
    call report(image % arity() .eq. 2 .and. image % num_tuples() .eq. 2, &
         & "projection collapses what it makes indistinct", nfail)
    call report(image % has([1, 1]) .and. image % has([2, 3]), &
         & "and holds exactly the projected set", nfail)

    image = project_slots(r, [2, 1])
    d = image % domain(1)
    call report(d % same_as(b), &
         & "the chosen order is structural: slot one of [2,1] is B", nfail)
    d = image % domain(2)
    call report(d % same_as(a), &
         & "and slot two is A", nfail)
    call report(image % has([1, 1]) .and. image % has([3, 2]), &
         & "with the tuples reversed to match", nfail)

    image = project_slots(r, [3])
    call report(image % arity() .eq. 1 .and. image % num_tuples() .eq. 2, &
         & "projection to one slot is a unary relation, deduplicated", nfail)

    ! The empty relation projects to the empty relation, carrying
    ! exactly the selected signature.
    none  = stored_relation('none', [a, b, c], &
         & reshape([integer ::], [3, 0]))
    image = project_slots(none, [3, 1])
    call report(image % num_tuples() .eq. 0 .and. image % arity() .eq. 2, &
         & "the empty relation projects to the empty relation", nfail)
    d = image % domain(1)
    call report(d % same_as(c), &
         & "whose first selected slot is C", nfail)
    d = image % domain(2)
    call report(d % same_as(a), &
         & "and whose second is A", nfail)

  end subroutine check_projection

  !===================================================================!
  ! Composition is existential: (a, c) stands wherever SOME b
  ! carries the chain, however many witnesses there are - the
  ! result is a set. The empty chain composes to the empty
  ! relation, which is an answer, not an error.
  !===================================================================!

  subroutine check_composition(nfail)

    integer, intent(inout) :: nfail

    type(counted_set)              :: a, b, c
    type(stored_relation)          :: r_ab, r_bc
    class(relation), allocatable   :: chained
    class(member_set), allocatable :: d

    a = counted_set('a-things', 2)
    b = counted_set('b-things', 3)
    c = counted_set('c-things', 2)

    ! a1 reaches c1 through b1 AND through b2: two witnesses, one
    ! tuple. a2 reaches nothing.
    r_ab = stored_relation('ab', [a, b], &
         & reshape([1,1,  1,2,  2,3], [2, 3]))
    r_bc = stored_relation('bc', [b, c], &
         & reshape([1,1,  2,1,  2,2], [2, 3]))

    chained = compose_binary(r_ab, r_bc)

    call report(chained % num_tuples() .eq. 2, &
         & "two witnesses, one tuple: composition is a set", nfail)
    call report(chained % has([1, 1]) .and. chained % has([1, 2]), &
         & "and holds exactly the chained pairs", nfail)
    call report(.not. chained % has([2, 1]) .and. &
         &      .not. chained % has([2, 2]), &
         & "a member with no chain relates to nothing", nfail)

    d = chained % domain(1)
    call report(d % same_as(a), &
         & "the result runs from the first source", nfail)
    d = chained % domain(2)
    call report(d % same_as(c), &
         & "to the second target", nfail)

    ! No b-chain at all: the empty composition is a relation.
    r_ab = stored_relation('ab', [a, b], reshape([1, 1], [2, 1]))
    r_bc = stored_relation('bc', [b, c], reshape([3, 1], [2, 1]))
    chained = compose_binary(r_ab, r_bc)
    call report(chained % num_tuples() .eq. 0, &
         & "no witness anywhere composes to the empty relation", nfail)

  end subroutine check_composition

end program test_graph_algebra
