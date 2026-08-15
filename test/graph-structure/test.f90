!=====================================================================!
! The related graph suite: the laws of the container (AGENTS.md,
! level 3, phase 4A).
!
! G = (S, R), and nothing else: the checks below build a graph out
! of a CALCULATOR - operations, values, ports - precisely so that
! no vertex or edge can sneak in as an assumption. The container
! must hold arbitrary sets and arbitrary-arity relations,
! refuse a relation over unheld domains, let two relations of one
! signature coexist, and - the ownership law made testable - hand
! out references into stable owned storage that a borrower's fibre
! views can ride while the graph lives.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_graph_structure

  use graph_identity       , only : token
  use graph_set        , only : index_set, set, unrelated_graph, declared_set
  use graph_relation       , only : stored_relation, relation, declared_domain
  use graph_binary_relation, only : binary_relation, csr_relation
  use graph_structure      , only : related_graph, declared_set, &
       &                            declared_relation, &
       &                            declare_relations, forget_relations

  implicit none

  integer :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "graph structure suite (AGENTS phase 4A)"
  write(*,'(1x,a)') "============================================="

  call check_generic_container(nfail)
  call check_same_signature_coexists(nfail)
  call check_stable_borrows(nfail)
  call check_graph_identity(nfail)
  call check_transition_maps(nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all structure checks passed"
  else
     write(*,'(1x,a,i0,a)') "FAILED: ", nfail, " structure check(s)"
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
  ! The calculator of AGENTS.md 8.2 as a related graph: three
  ! sets with no graph words on them, one ternary relation, one
  ! binary. The container holds what it is handed and answers what
  ! it holds.
  !===================================================================!

  subroutine check_generic_container(nfail)

    integer, intent(inout) :: nfail

    type(index_set)              :: ops, vals, ports, unheld
    type(stored_relation)          :: flow
    type(csr_relation)             :: dep
    type(related_graph), target :: g
    class(set), pointer     :: sp
    class(relation)  , pointer     :: rp

    ops   = index_set('operations', 3)
    vals  = index_set('values'    , 5)
    ports = index_set('ports'     , 2)     ! 1 = input, 2 = output

    ! (+ reads value 2 on input) ... (* writes value 5 on output)
    flow = stored_relation('flow', &
         & [declared_domain(ops), declared_domain(vals), declared_domain(ports)], &
         & reshape([1,2,1,  1,3,1,  1,5,2,  2,5,1,  2,4,2], [3, 5]))

    dep = csr_relation('feeds', ops, ops, reshape([1, 2], [2, 1]))

    g = related_graph('calculator', &
         & [declared_set(ops), declared_set(vals), declared_set(ports)], &
         & [declared_relation(flow), declared_relation(dep)])

    call report(g % num_sets() .eq. 3 .and. &
         &      g % num_relations() .eq. 2, &
         & "the container counts what it was handed", nfail)

    sp => g % set_at(2)
    call report(sp % equals(vals), &
         & "a seat answers the declared domain, by identity", nfail)

    rp => g % relation_at(1)
    call report(rp % equals(flow) .and. rp % arity() .eq. 3, &
         & "a ternary relation sits beside a binary one, untroubled", &
         & nfail)

    call report(g % holds_set(ports), &
         & "the graph knows its own domains", nfail)
    unheld = index_set('a domain this graph does not hold', 2)
    call report(.not. g % holds_set(unheld), &
         & "and disowns everyone else's", nfail)

  end subroutine check_generic_container

  !===================================================================!
  ! Relation identity is the address, never the signature: two
  ! relations over the same slots coexist in one graph.
  !===================================================================!

  subroutine check_same_signature_coexists(nfail)

    integer, intent(inout) :: nfail

    type(index_set)              :: ops
    type(csr_relation)             :: physical, scheduled
    type(related_graph), target :: g
    class(relation), pointer       :: r1, r2

    ops = index_set('operations', 3)

    physical  = csr_relation('feeds' , ops, ops, reshape([1, 2], [2, 1]))
    scheduled = csr_relation('awaits', ops, ops, reshape([2, 3], [2, 1]))

    g = related_graph('plan', [declared_set(ops)], &
         & [declared_relation(physical), declared_relation(scheduled)])

    r1 => g % relation_at(1)
    r2 => g % relation_at(2)
    call report(.not. r1 % equals(r2), &
         & "same signature, two relations - no collision", nfail)
    call report(r1 % has([1, 2]) .and. .not. r2 % has([1, 2]), &
         & "and each keeps its own tuples", nfail)

  end subroutine check_same_signature_coexists

  !===================================================================!
  ! The ownership law, made testable: a reference obtained from the
  ! graph is returned from the accessor, survives the return, and
  ! its fibre borrows read true while the graph lives. No copies
  ! anywhere on this road.
  !===================================================================!

  subroutine check_stable_borrows(nfail)

    integer, intent(inout) :: nfail

    type(index_set)               :: ops
    type(csr_relation)              :: dep
    type(related_graph), target  :: g
    class(binary_relation), pointer :: br
    integer, pointer                :: f(:)

    ops = index_set('operations', 4)
    dep = csr_relation('feeds', ops, ops, &
         & reshape([1,2,  1,3,  2,4,  3,4], [2, 4]))

    g = related_graph('plan', [declared_set(ops)], [declared_relation(dep)])

    br => borrow_binary(g, 1)
    call report(associated(br), &
         & "a reference returned from an accessor survives the return", &
         & nfail)

    f => br % image_view(1)
    call report(size(f) .eq. 2 .and. all(f .eq. [2, 3]), &
         & "and its fibre borrow reads true while the graph lives", nfail)

    f => br % preimage_view(4)
    call report(size(f) .eq. 2 .and. all(f .eq. [2, 3]), &
         & "both directions, no copy anywhere on the road", nfail)

  end subroutine check_stable_borrows

  !===================================================================!
  ! borrow_binary walks the accessor road a caller would: reference
  ! out of the graph, narrowed to binary, handed back - a borrow
  ! crossing two procedure boundaries before its first use.
  !===================================================================!

  function borrow_binary(g, k) result(br)

    class(related_graph), target, intent(in) :: g
    integer                        , intent(in) :: k
    class(binary_relation), pointer             :: br

    class(relation), pointer :: rp

    br => null()
    rp => g % relation_at(k)
    select type (rp)
    class is (binary_relation)
       br => rp
    end select

  end function borrow_binary

  !===================================================================!
  ! The graph is the third citizen on the identity roll.
  !===================================================================!

  subroutine check_graph_identity(nfail)

    integer, intent(inout) :: nfail

    type(index_set)        :: ops
    type(unrelated_graph)  :: g, h, copy
    type(token)            :: tk

    ops = index_set('operations', 2)

    ! Identically stocked and still two graphs, because identity is
    ! declared and never extensional. No relation is declared here:
    ! the claim is about identity, so the object is G = (S, empty).
    g = unrelated_graph('one', [declared_set(ops)])
    h = unrelated_graph('two', [declared_set(ops)])

    tk = g % id()
    call report(tk % declared(), &
         & "a declared graph carries a token that has signed", nfail)
    call report(.not. g % equals(h), &
         & "two declarations are two graphs", nfail)
    copy = g
    call report(copy % equals(g), &
         & "a copy is the same declared graph", nfail)
    call report(g % name() == 'one', &
         & "and the name is the reader's, as everywhere", nfail)

  end subroutine check_graph_identity


  !===================================================================!
  ! THE TRANSITION MAPS. Declaring relations carries (S, empty) to
  ! (S, R); forgetting them carries it back. S survives both by
  ! IDENTITY - not merely by extension - and the graph identity
  ! survives neither, because each map declares a new object.
  !===================================================================!

  subroutine check_transition_maps(nfail)

    integer, intent(inout) :: nfail

    type(index_set)      , target :: ops, vals
    type(csr_relation)            :: dep
    type(unrelated_graph), target :: g0, back
    type(related_graph)  , target :: g1
    class(set)           , pointer :: p
    class(relation)      , pointer :: rp
    integer                        :: s
    logical                        :: kept

    ops  = index_set('operations', 3)
    vals = index_set('values'    , 5)
    dep  = csr_relation('feeds', ops, ops, reshape([1, 2], [2, 1]))

    g0 = unrelated_graph('start', [declared_set(ops), declared_set(vals)])
    g1 = declare_relations(g0, [declared_relation(dep)], 'declared')

    call report(g0 % num_relations() .eq. 0 .and. g1 % num_relations() .eq. 1, &
         & "declaring relations carries |R| from 0 to what it was handed", nfail)
    call report(g0 % num_relations() .eq. 0, &
         & "and leaves the source untouched - a map, never a mutation", nfail)
    call report(g1 % num_sets() .eq. g0 % num_sets(), &
         & "S is preserved in count", nfail)

    kept = .true.
    do s = 1, g1 % num_sets()
       p => g1 % set_at(s)
       kept = kept .and. p % equals(g0 % set_at(s))
    end do
    call report(kept, &
         & "and preserved BY IDENTITY: every set is the set it came from", nfail)

    rp => g1 % relation_at(1)
    call report(rp % equals(dep), &
         & "the declared relation is the relation that was handed in", nfail)
    call report(.not. g1 % equals(g0), &
         & "and the result is a NEW graph - identity is not carried", nfail)

    back = forget_relations(g1, 'forgotten')

    call report(back % num_relations() .eq. 0, &
         & "forgetting relations empties R", nfail)
    call report(back % num_sets() .eq. g0 % num_sets(), &
         & "and leaves S standing", nfail)

    kept = .true.
    do s = 1, back % num_sets()
       p => back % set_at(s)
       kept = kept .and. p % equals(g0 % set_at(s))
    end do
    call report(kept, &
         & "S survives the round trip by identity, both ways", nfail)
    call report(.not. back % equals(g0) .and. .not. back % equals(g1), &
         & "while the round trip lands on a THIRD graph, never the first", nfail)

  end subroutine check_transition_maps

end program test_graph_structure
