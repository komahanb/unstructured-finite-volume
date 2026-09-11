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

  use graph_fractal        , only : graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set        , only : set_map
  use map_inclusion  , only : inclusion_map, declared_subobject
  use relation_finitary       , only : stored_relation, relation
  use relation_algebra, only : restrict_slot, project_slots, &
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

  subroutine report(satisfied, label, nfail)

    logical         , intent(in)    :: satisfied
    character(len=*), intent(in)    :: label
    integer         , intent(inout) :: nfail

    if (satisfied) then
       write(*,'(1x,a,a)') "PASS : ", label
    else
       write(*,'(1x,a,a)') "FAIL : ", label
       nfail = nfail + 1
    end if

  end subroutine report

  !===================================================================!
  ! Restriction retains exactly the admitted tuples, preserves the
  ! signature unchanged, and restricting by the full domain itself
  ! is the lawful identity restriction (A embeds in A).
  !===================================================================!

  subroutine check_restriction(nfail)

    integer, intent(inout) :: nfail

    type(graph)              :: a, b, c
    type(graph)               :: some_b, empty_subset
    type(stored_relation)          :: r, restricted_relation
    type(graph) :: d
    integer, allocatable           :: rt(:,:)
    integer                        :: j
    logical                        :: satisfied
    type(set_map)     :: sets
    type(inclusion_map)     :: inclusions

    call a % declare()
    call sets % bind(a, counted_set_representation(3))
    call b % declare()
    call sets % bind(b, counted_set_representation(4))
    call c % declare()
    call sets % bind(c, counted_set_representation(2))

    r = stored_relation('r', [a, b, c], &
         & reshape([1,1,1,  1,2,2,  2,2,1,  3,4,2], [3, 4]), sets)

    call some_b % declare()
    call sets       % bind(some_b, listed_set_representation([2]))
    call inclusions % include_in(some_b, b)

    restricted_relation = restrict_slot(r, 2, some_b, sets, inclusions)

    call report(restricted_relation % num_tuples() .eq. 2, &
         & "restriction keeps exactly the admitted tuples", nfail)
    call report(restricted_relation % has([1, 2, 2]) .and. restricted_relation % has([2, 2, 1]), &
         & "and they are the right ones", nfail)
    call report(restricted_relation % arity() .eq. 3, &
         & "the signature's arity is unchanged", nfail)
    d = restricted_relation % domain(2)
    call report(d % same_as(b), &
         & "and the restricted position still returns its full domain", nfail)

    ! Full-domain restriction is the identity, extensionally: equal
    ! count and every original tuple present - for two sets of equal
    ! finite size, that is equality.
    restricted_relation = restrict_slot(r, 2, b, sets, inclusions)
    call r % tuples(rt)
    satisfied = restricted_relation % num_tuples() .eq. r % num_tuples()
    do j = 1, size(rt, 2)
       satisfied = satisfied .and. restricted_relation % has(rt(:, j))
    end do
    call report(satisfied, &
         & "restricting by the full domain is the identity, as sets", nfail)

    ! The empty subset admits nothing; the signature is unchanged.
    call empty_subset % declare()
    call sets       % bind(empty_subset, listed_set_representation([integer ::]))
    call inclusions % include_in(empty_subset, b)
    restricted_relation = restrict_slot(r, 2, empty_subset, sets, inclusions)
    call report(restricted_relation % num_tuples() .eq. 0, &
         & "restriction by the empty subset is the empty relation", nfail)
    call report(restricted_relation % arity() .eq. 3, &
         & "whose arity is the original's", nfail)
    d = restricted_relation % domain(1)
    satisfied = d % same_as(a)
    d = restricted_relation % domain(2)
    satisfied = satisfied .and. d % same_as(b)
    d = restricted_relation % domain(3)
    satisfied = satisfied .and. d % same_as(c)
    call report(satisfied, &
         & "and whose signature is the original's, slot for slot", nfail)

  end subroutine check_restriction

  !===================================================================!
  ! Projection returns exactly the selected positions in the chosen
  ! order - [2,1] is the reversed signature - and its image is a
  ! set: tuples that collapse, collapse.
  !===================================================================!

  subroutine check_projection(nfail)

    integer, intent(inout) :: nfail

    type(graph)              :: a, b, c
    type(stored_relation)          :: r, image, none
    type(graph) :: d
    type(set_map)     :: sets

    call a % declare()
    call sets % bind(a, counted_set_representation(3))
    call b % declare()
    call sets % bind(b, counted_set_representation(4))
    call c % declare()
    call sets % bind(c, counted_set_representation(2))

    ! Two tuples agree on (slot1, slot2); they differ only in slot 3.
    r = stored_relation('r', [a, b, c], &
         & reshape([1,1,1,  1,1,2,  2,3,1], [3, 3]), sets)

    image = project_slots(r, [1, 2], sets)
    call report(image % arity() .eq. 2 .and. image % num_tuples() .eq. 2, &
         & "projection collapses what it makes indistinct", nfail)
    call report(image % has([1, 1]) .and. image % has([2, 3]), &
         & "and contains exactly the projected set", nfail)

    image = project_slots(r, [2, 1], sets)
    d = image % domain(1)
    call report(d % same_as(b), &
         & "the chosen order is structural: slot one of [2,1] is B", nfail)
    d = image % domain(2)
    call report(d % same_as(a), &
         & "and slot two is A", nfail)
    call report(image % has([1, 1]) .and. image % has([3, 2]), &
         & "with the tuples reversed to match", nfail)

    image = project_slots(r, [3], sets)
    call report(image % arity() .eq. 1 .and. image % num_tuples() .eq. 2, &
         & "projection to one slot is a unary relation, deduplicated", nfail)

    ! The empty relation projects to the empty relation, preserving
    ! exactly the selected signature.
    none  = stored_relation('none', [a, b, c], &
         & reshape([integer ::], [3, 0]), sets)
    image = project_slots(none, [3, 1], sets)
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
  ! Composition is existential: (a, c) exists wherever SOME b
  ! satisfies both constituent relations, however many witnesses there are - the
  ! result is a set. The empty chain composes to the empty
  ! relation, which is a valid result.
  !===================================================================!

  subroutine check_composition(nfail)

    integer, intent(inout) :: nfail

    type(graph)              :: a, b, c
    type(stored_relation)          :: p_ab, p_bc
    class(relation), allocatable   :: chained
    type(graph) :: d
    type(set_map)     :: sets

    call a % declare()
    call sets % bind(a, counted_set_representation(2))
    call b % declare()
    call sets % bind(b, counted_set_representation(3))
    call c % declare()
    call sets % bind(c, counted_set_representation(2))

    ! a1 reaches c1 through b1 AND through b2: two witnesses, one
    ! tuple. a2 reaches nothing.
    p_ab = stored_relation('ab', [a, b], &
         & reshape([1,1,  1,2,  2,3], [2, 3]), sets)
    p_bc = stored_relation('bc', [b, c], &
         & reshape([1,1,  2,1,  2,2], [2, 3]), sets)

    chained = compose_binary(p_ab, p_bc, sets)

    call report(chained % num_tuples() .eq. 2, &
         & "two witnesses, one tuple: composition is a set", nfail)
    call report(chained % has([1, 1]) .and. chained % has([1, 2]), &
         & "and contains exactly the composed pairs", nfail)
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
    p_ab = stored_relation('ab', [a, b], reshape([1, 1], [2, 1]), sets)
    p_bc = stored_relation('bc', [b, c], reshape([3, 1], [2, 1]), sets)
    chained = compose_binary(p_ab, p_bc, sets)
    call report(chained % num_tuples() .eq. 0, &
         & "no witness anywhere composes to the empty relation", nfail)

  end subroutine check_composition

end program test_graph_algebra
