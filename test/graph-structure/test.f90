!=====================================================================!
! The relational graph suite: the laws of the container (AGENTS.md,
! level 3, phase 4A).
!
! G = (S, R), and nothing else: the checks below build a graph out
! of a CALCULATOR - operations, values, ports - precisely so that
! no vertex or edge can sneak in as an assumption. The container
! must hold arbitrary carriers and arbitrary-arity relations,
! refuse a relation over foreign domains, let two relations of one
! signature coexist, and - the ownership law made testable - hand
! out references into stable owned storage that a borrower's fibre
! views can ride while the graph lives.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_graph_structure

  use graph_identity       , only : token
  use graph_carrier        , only : counted_set, member_set
  use graph_relation       , only : stored_relation, relation, slot
  use graph_binary_relation, only : binary_relation, csr_relation
  use graph_structure      , only : relational_graph, held_set, &
       &                            held_relation

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
  ! The calculator of AGENTS.md 8.2 as a relational graph: three
  ! carriers with no graph words on them, one ternary relation, one
  ! binary. The container holds what it is handed and answers what
  ! it holds.
  !===================================================================!

  subroutine check_generic_container(nfail)

    integer, intent(inout) :: nfail

    type(counted_set)              :: ops, vals, ports, foreign
    type(stored_relation)          :: flow
    type(csr_relation)             :: dep
    type(relational_graph), target :: g
    class(member_set), pointer     :: sp
    class(relation)  , pointer     :: rp

    ops   = counted_set('operations', 3)
    vals  = counted_set('values'    , 5)
    ports = counted_set('ports'     , 2)     ! 1 = input, 2 = output

    ! (+ reads value 2 on input) ... (* writes value 5 on output)
    flow = stored_relation('flow', &
         & [slot(ops), slot(vals), slot(ports)], &
         & reshape([1,2,1,  1,3,1,  1,5,2,  2,5,1,  2,4,2], [3, 5]))

    dep = csr_relation('feeds', ops, ops, reshape([1, 2], [2, 1]))

    g = relational_graph('calculator', &
         & [held_set(ops), held_set(vals), held_set(ports)], &
         & [held_relation(flow), held_relation(dep)])

    call report(g % num_member_sets() .eq. 3 .and. &
         &      g % num_relations() .eq. 2, &
         & "the container counts what it was handed", nfail)

    sp => g % member_set_at(2)
    call report(sp % same_as(vals), &
         & "a seat answers the declared domain, by identity", nfail)

    rp => g % relation_at(1)
    call report(rp % same_as(flow) .and. rp % arity() .eq. 3, &
         & "a ternary relation sits beside a binary one, untroubled", &
         & nfail)

    call report(g % holds_set(ports), &
         & "the graph knows its own domains", nfail)
    foreign = counted_set('foreign', 2)
    call report(.not. g % holds_set(foreign), &
         & "and disowns everyone else's", nfail)

  end subroutine check_generic_container

  !===================================================================!
  ! Relation identity is the address, never the signature: two
  ! relations over the same slots coexist in one graph.
  !===================================================================!

  subroutine check_same_signature_coexists(nfail)

    integer, intent(inout) :: nfail

    type(counted_set)              :: ops
    type(csr_relation)             :: physical, scheduled
    type(relational_graph), target :: g
    class(relation), pointer       :: r1, r2

    ops = counted_set('operations', 3)

    physical  = csr_relation('feeds' , ops, ops, reshape([1, 2], [2, 1]))
    scheduled = csr_relation('awaits', ops, ops, reshape([2, 3], [2, 1]))

    g = relational_graph('plan', [held_set(ops)], &
         & [held_relation(physical), held_relation(scheduled)])

    r1 => g % relation_at(1)
    r2 => g % relation_at(2)
    call report(.not. r1 % same_as(r2), &
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

    type(counted_set)               :: ops
    type(csr_relation)              :: dep
    type(relational_graph), target  :: g
    class(binary_relation), pointer :: br
    integer, pointer                :: f(:)

    ops = counted_set('operations', 4)
    dep = csr_relation('feeds', ops, ops, &
         & reshape([1,2,  1,3,  2,4,  3,4], [2, 4]))

    g = relational_graph('plan', [held_set(ops)], [held_relation(dep)])

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

    class(relational_graph), target, intent(in) :: g
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

    type(counted_set)      :: ops
    type(relational_graph) :: g, h, copy
    type(token)            :: tk
    type(held_relation)    :: none(0)

    ops = counted_set('operations', 2)

    g = relational_graph('one', [held_set(ops)], none)
    h = relational_graph('two', [held_set(ops)], none)

    tk = g % id()
    call report(tk % declared(), &
         & "a declared graph carries a token that has signed", nfail)
    call report(.not. g % same_as(h), &
         & "two declarations are two graphs", nfail)
    copy = g
    call report(copy % same_as(g), &
         & "a copy is the same declared graph", nfail)
    call report(g % name() == 'one', &
         & "and the name is the reader's, as everywhere", nfail)

  end subroutine check_graph_identity

end program test_graph_structure
