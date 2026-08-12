!=====================================================================!
! The carrier suite: the laws of the member set (AGENTS.md, level 0,
! phase 1).
!
! A member set answers id, name, size, and members - and its
! identity is STRUCTURAL: the same declared domain, never the same
! numbers. The checks below are the acceptance laws of section 4,
! plus the phase-1 wiring: a stored graph hands out its vertex and
! edge carriers beside the old vocabulary, and hands out the SAME
! domain every time it is asked.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_graph_carrier

  use graph_carrier, only : member_set, counted_set
  use class_graph  , only : stored_graph

  implicit none

  integer :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "graph carrier suite (AGENTS phase 1)"
  write(*,'(1x,a)') "============================================="

  call check_counted_contract(nfail)
  call check_structural_identity(nfail)
  call check_graph_hands_out_carriers(nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all carrier checks passed"
  else
     write(*,'(1x,a,i0,a)') "FAILED: ", nfail, " carrier check(s)"
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
  ! The contract: id, name, size, member, members - each answered
  ! meaningfully, including on the empty domain.
  !===================================================================!

  subroutine check_counted_contract(nfail)

    integer, intent(inout) :: nfail

    type(counted_set)    :: cells, none
    integer, allocatable :: idx(:)
    integer              :: k
    logical              :: ok

    cells = counted_set('cells', 6)

    call report(cells % id() > 0, &
         & "a declared domain carries a nonzero stamp", nfail)
    call report(cells % name() == 'cells', &
         & "and the name it was declared with", nfail)
    call report(cells % size() .eq. 6, &
         & "size counts the members", nfail)
    call report(cells % member(4) .eq. 4, &
         & "the counted set's k-th member is k", nfail)

    call cells % members(idx)
    ok = size(idx) .eq. 6
    do k = 1, size(idx)
       ok = ok .and. (idx(k) .eq. k)
    end do
    call report(ok, &
         & "members enumerates 1 to n, in order", nfail)

    ! Domain validity (AGENTS.md section 61): every member the set
    ! reports belongs to the set.
    ok = .true.
    do k = 1, size(idx)
       ok = ok .and. (idx(k) >= 1) .and. (idx(k) <= cells % size())
    end do
    call report(ok, &
         & "every reported member belongs to the domain", nfail)

    ! Membership is a primitive: one comparison, never an
    ! enumeration and a search. This is what a relation signature
    ! leans on when it validates a tuple.
    call report(cells % has(1) .and. cells % has(4) .and. cells % has(6), &
         & "has says yes to every member", nfail)
    call report(.not. (cells % has(0) .or. cells % has(7) .or. &
         &             cells % has(-3)), &
         & "and no to everything outside the domain", nfail)

    none = counted_set('nothing', 0)
    call none % members(idx)
    call report(none % size() .eq. 0 .and. size(idx) .eq. 0, &
         & "the empty domain is a domain", nfail)
    call report(.not. none % has(1), &
         & "and holds nothing", nfail)

  end subroutine check_counted_contract

  !===================================================================!
  ! Identity is structural. Equal contents prove nothing; a copy is
  ! the same domain; an undeclared set is no domain at all.
  !===================================================================!

  subroutine check_structural_identity(nfail)

    integer, intent(inout) :: nfail

    type(counted_set) :: cells, faces, again, copy, raw

    cells = counted_set('cells', 4)
    faces = counted_set('faces', 4)

    call report(.not. cells % same_as(faces), &
         & "four cells and four faces are different worlds", nfail)

    again = counted_set('cells', 4)
    call report(.not. cells % same_as(again), &
         & "a second declaration is a second domain, same name or not", nfail)

    copy = cells
    call report(copy % same_as(cells) .and. cells % same_as(copy), &
         & "a copy is the same declared domain, both ways round", nfail)

    call report(cells % same_as(cells), &
         & "a domain is itself", nfail)

    call report(.not. raw % same_as(raw), &
         & "an undeclared set equals nothing, itself included", nfail)

    call report(cells % id() /= faces % id() .and. &
         &      cells % id() /= again % id(), &
         & "no two declarations share a stamp", nfail)

  end subroutine check_structural_identity

  !===================================================================!
  ! The phase-1 wiring: a stored graph declares its two carriers at
  ! construction and answers the SAME domain at every asking, while
  ! all the old vocabulary stands untouched beside them.
  !===================================================================!

  subroutine check_graph_hands_out_carriers(nfail)

    integer, intent(inout) :: nfail

    type(stored_graph) :: g, h
    type(counted_set)  :: vs, es, vs2

    g = stored_graph(3, tails=[1, 2], heads=[2, 3])
    h = stored_graph(3, tails=[1, 2], heads=[2, 3])

    vs = g % vertex_set()
    es = g % edge_set()

    call report(vs % size() .eq. g % num_vertices(), &
         & "the vertex carrier counts what num_vertices counts", nfail)
    call report(es % size() .eq. g % num_edges(), &
         & "the edge carrier counts what num_edges counts", nfail)
    call report(vs % name() == 'vertices' .and. es % name() == 'edges', &
         & "the carriers wear the graph's own words", nfail)

    call report(.not. vs % same_as(es), &
         & "a graph's two sides are two domains", nfail)

    vs2 = g % vertex_set()
    call report(vs % same_as(vs2), &
         & "asked twice, the graph answers the same domain", nfail)

    call report(.not. vs % same_as(h % vertex_set()), &
         & "an identical twin graph still owns its own domains", nfail)

  end subroutine check_graph_hands_out_carriers

end program test_graph_carrier
