!=====================================================================!
! The set suite: the laws of the member set (AGENTS.md, level 0,
! phase 1).
!
! A member set answers id, name, size, and members - and its
! identity is STRUCTURAL: the same declared domain, never the same
! numbers. The checks below are the acceptance laws of section 4,
! plus the phase-1 wiring: a stored graph hands out its vertex and
! edge sets beside the old vocabulary, and hands out the SAME
! domain every time it is asked.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_graph_set

  use graph_identity, only : token
  use graph_set , only : set, index_set, subset
  use graph_set , only : unrelated_graph, declared_set
  use class_graph   , only : stored_graph

  implicit none

  integer :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "graph set suite (AGENTS phase 1)"
  write(*,'(1x,a)') "============================================="

  call check_counted_contract(nfail)
  call check_structural_identity(nfail)
  call check_subsets(nfail)
  call check_embedding_order(nfail)
  call check_empty_subset(nfail)
  call check_graph_hands_out_sets(nfail)
  call check_unrelated_graph(nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all set checks passed"
  else
     write(*,'(1x,a,i0,a)') "FAILED: ", nfail, " set check(s)"
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

    type(index_set)    :: cells, none
    type(token)          :: t
    integer, allocatable :: idx(:)
    integer              :: k
    logical              :: ok

    cells = index_set('cells', 6)

    ! id answers the whole opaque token, never a bare local integer.
    t = cells % id()
    call report(t % declared(), &
         & "a declared domain carries a token that has signed", nfail)
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

    ! The inverse enumeration: where does a member stand. Zero for
    ! an outsider, and the two round-trip laws that make the pair
    ! member/local_index a bijection on the domain.
    call report(cells % local_index(4) .eq. 4 .and. &
         &      cells % local_index(0) .eq. 0 .and. &
         &      cells % local_index(7) .eq. 0, &
         & "local_index answers the standing, zero for outsiders", nfail)
    ok = .true.
    do k = 1, cells % size()
       ok = ok .and. (cells % member(cells % local_index(k)) .eq. k)
       ok = ok .and. (cells % local_index(cells % member(k)) .eq. k)
    end do
    call report(ok, &
         & "member and local_index invert each other, both ways", nfail)

    none = index_set('nothing', 0)
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

    type(index_set) :: cells, faces, again, copy, raw
    type(token)       :: tc, tf, ta

    cells = index_set('cells', 4)
    faces = index_set('faces', 4)

    call report(.not. cells % equals(faces), &
         & "four cells and four faces are different worlds", nfail)

    again = index_set('cells', 4)
    call report(.not. cells % equals(again), &
         & "a second declaration is a second domain, same name or not", nfail)

    copy = cells
    call report(copy % equals(cells) .and. cells % equals(copy), &
         & "a copy is the same declared domain, both ways round", nfail)

    call report(cells % equals(cells), &
         & "a domain is itself", nfail)

    call report(.not. raw % equals(raw), &
         & "an undeclared set equals nothing, itself included", nfail)

    tc = cells % id()
    tf = faces % id()
    ta = again % id()
    call report(.not. tc % matches(tf) .and. .not. tc % matches(ta), &
         & "no two declarations share a stamp", nfail)

  end subroutine check_structural_identity

  !===================================================================!
  ! The subobject (AGENTS.md 6 and 37, phase 5): S c--> A. A subset
  ! IS a member set - it answers the whole five-question contract,
  ! signs its own identity, and adds one law: every member belongs
  ! to the ambient, sealed at construction. This is the destination
  ! the old edgeless-graph support was reaching for.
  !===================================================================!

  subroutine check_subsets(nfail)

    integer, intent(inout) :: nfail

    type(index_set)              :: faces
    type(subset)               :: walls, hot
    class(set), allocatable :: amb
    integer, allocatable           :: idx(:)
    integer                        :: k
    logical                        :: ok

    faces = index_set('faces', 6)
    walls = subset('walls', faces, [2, 5, 2])

    call report(walls % size() .eq. 2, &
         & "a member handed in twice is in the subset once", nfail)
    call walls % members(idx)
    call report(all(idx .eq. [2, 5]), &
         & "first appearances stand, in their first order", nfail)

    call report(walls % has(2) .and. walls % has(5) .and. &
         &      .not. walls % has(3) .and. .not. walls % has(9), &
         & "membership answers the chosen family alone", nfail)

    ok = .true.
    do k = 1, walls % size()
       ok = ok .and. (walls % local_index(walls % member(k)) .eq. k)
       ok = ok .and. (walls % member(walls % local_index( &
            &         walls % member(k))) .eq. walls % member(k))
    end do
    call report(ok, &
         & "member and local_index invert each other on the subset", nfail)

    amb = walls % ambient()
    call report(amb % equals(faces), &
         & "the ambient is the domain the subset was carved from", nfail)
    call report(.not. walls % equals(faces), &
         & "a subobject is its own declared domain, not its host", nfail)

    ok = .true.
    do k = 1, walls % size()
       ok = ok .and. amb % has(walls % member(k))
    end do
    call report(ok, &
         & "the inclusion law: every member belongs to the ambient", nfail)

    ! Subobjects nest: a subset of a subset, each layer keeping its
    ! own ambient.
    hot = subset('hot-walls', walls, [5])
    call report(hot % size() .eq. 1 .and. hot % member(1) .eq. 5, &
         & "a subset of a subset is a subset", nfail)
    amb = hot % ambient()
    call report(amb % equals(walls), &
         & "whose ambient is the layer above, not the ground", nfail)
    call report(faces % has(hot % member(1)), &
         & "and whose members reach the ground transitively", nfail)

  end subroutine check_subsets

  !===================================================================!
  ! The embedding order: A precedes A; a subset precedes its whole
  ! ambient chain; nothing else precedes anything. This query is
  ! what replaces the old side flag - a consumer asks where a
  ! domain ultimately lives, and no select type ever answers.
  !===================================================================!

  subroutine check_embedding_order(nfail)

    integer, intent(inout) :: nfail

    type(index_set) :: faces, other
    type(subset)  :: walls, hot, cold

    faces = index_set('faces', 8)
    other = index_set('other', 8)
    walls = subset('walls', faces, [2, 5, 7])
    hot   = subset('hot'  , walls, [5, 7])
    cold  = subset('cold' , faces, [1])

    call report(faces % is_subobject_of(faces), &
         & "a domain is a subobject of itself", nfail)
    call report(.not. faces % is_subobject_of(other), &
         & "and of no unrelated domain, sizes notwithstanding", nfail)

    call report(walls % is_subobject_of(faces), &
         & "a subset is a subobject of its ambient", nfail)
    call report(hot % is_subobject_of(walls) .and. &
         &      hot % is_subobject_of(faces), &
         & "and of the whole chain above it, transitively", nfail)
    call report(hot % is_subobject_of(hot), &
         & "while remaining a subobject of itself", nfail)

    call report(.not. faces % is_subobject_of(walls), &
         & "the order points one way: the ambient never descends", nfail)
    call report(.not. cold % is_subobject_of(walls), &
         & "and siblings of one ambient are unrelated", nfail)
    call report(.not. hot % is_subobject_of(other), &
         & "no chain reaches a domain it was never carved from", nfail)

  end subroutine check_embedding_order

  !===================================================================!
  ! The empty subset is a valid domain: declared, counted at zero,
  ! holding nothing, still embedded where it was carved.
  !===================================================================!

  subroutine check_empty_subset(nfail)

    integer, intent(inout) :: nfail

    type(index_set)              :: faces
    type(subset)               :: none
    class(set), allocatable :: amb
    integer, allocatable           :: idx(:)

    faces = index_set('faces', 6)
    none  = subset('nothing', faces, [integer ::])

    call none % members(idx)
    call report(none % size() .eq. 0 .and. size(idx) .eq. 0, &
         & "the empty subset is a domain", nfail)
    call report(.not. none % has(2) .and. none % local_index(2) .eq. 0, &
         & "that holds nothing and places nothing", nfail)

    amb = none % ambient()
    call report(amb % equals(faces), &
         & "yet remembers what it was carved from", nfail)
    call report(none % is_subobject_of(faces), &
         & "and stands embedded there, empty or not", nfail)

  end subroutine check_empty_subset

  !===================================================================!
  ! The phase-1 wiring: a stored graph declares its two sets at
  ! construction and answers the SAME domain at every asking, while
  ! all the old vocabulary stands untouched beside them.
  !===================================================================!

  subroutine check_graph_hands_out_sets(nfail)

    integer, intent(inout) :: nfail

    type(stored_graph) :: g, h
    type(index_set)  :: vs, es, vs2

    g = stored_graph(3, tails=[1, 2], heads=[2, 3])
    h = stored_graph(3, tails=[1, 2], heads=[2, 3])

    vs = g % vertex_set()
    es = g % edge_set()

    call report(vs % size() .eq. g % num_vertices(), &
         & "the vertex set counts what num_vertices counts", nfail)
    call report(es % size() .eq. g % num_edges(), &
         & "the edge set counts what num_edges counts", nfail)
    call report(vs % name() == 'vertices' .and. es % name() == 'edges', &
         & "the sets wear the graph's own words", nfail)

    call report(.not. vs % equals(es), &
         & "a graph's two sides are two domains", nfail)

    vs2 = g % vertex_set()
    call report(vs % equals(vs2), &
         & "asked twice, the graph answers the same domain", nfail)

    call report(.not. vs % equals(h % vertex_set()), &
         & "an identical twin graph still owns its own domains", nfail)

  end subroutine check_graph_hands_out_sets


  !===================================================================!
  ! THE RELATION-FREE GRAPH, G = (S, empty), for every S.
  !
  ! The empty graph and the three-domain graph are the two cases the
  ! abstract parent could not express; the one-domain case is here
  ! beside them to show that `set` is a SPECIALIZATION of this
  ! object and not a rival to it.
  !===================================================================!

  subroutine check_unrelated_graph(nfail)

    integer, intent(inout) :: nfail

    type(index_set)      , target :: a, b, c
    type(unrelated_graph), target :: empty, one, three
    type(declared_set)                :: nothing(0)
    class(set)           , pointer :: p, q
    type(token)                    :: stamp

    a = index_set('cells', 4)
    b = index_set('faces', 6)
    c = index_set('parts', 2)

    empty = unrelated_graph('the empty graph', nothing)
    one   = unrelated_graph('one domain',      [declared_set(a)])
    three = unrelated_graph('three domains',   [declared_set(a), declared_set(b), declared_set(c)])

    call report(empty % num_sets() .eq. 0, &
         & "the empty graph G = (empty, empty) is a declared object", nfail)
    call report(three % num_sets() .eq. 3, &
         & "and a graph may hold three domains with nothing between them", nfail)
    call report(one % num_sets() .eq. 1, &
         & "one domain needs no special type: it is the same object", nfail)

    call report(empty % num_relations() .eq. 0 .and. &
         &      one   % num_relations() .eq. 0 .and. &
         &      three % num_relations() .eq. 0, &
         & "|R| = 0 is the law, and it does not vary", nfail)

    p => three % set_at(2)
    call report(p % equals(b), &
         & "set_at answers the declared domain, by identity", nfail)
    q => three % set_at(2)
    call report(associated(p, q), &
         & "and BORROWS one owned seat - asked twice, the same storage", nfail)

    call report(three % holds_set(a) .and. three % holds_set(b) .and. &
         &      three % holds_set(c), &
         & "a graph knows every domain it holds", nfail)
    call report(.not. three % holds_set(index_set('cells', 4)), &
         & "and disowns an identical stranger", nfail)

    stamp = empty % id()
    call report(stamp % declared(), &
         & "the empty graph signs like any other", nfail)

    call report(.not. empty % equals(one), &
         & "two declarations are two graphs, empty or not", nfail)

  end subroutine check_unrelated_graph

end program test_graph_set
