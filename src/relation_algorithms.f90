!=====================================================================!
! LEVEL 4 OF THE NEW TOWER . THE GRAPH ALGORITHMS
!
! The level defines WHAT GRAPH-THEORETIC PREDICATES CAN BE
! EVALUATED on an interpretation. The algorithms are defined HERE, as
! free module procedures over a binary relation on one domain - never
! as methods on the relation or on any container:
! traversal acts on structure, it is not structure (AGENTS.md 18,
! CALCULATOR.md 11). This module contains exactly what its first
! caller - the calculator's dependency traversal - requires:
!
!      sources      the members with empty preimage
!      sinks        the members with empty image
!      reachable    whether a directed path exists
!      topological_order   every arrow forward, or rejection
!
! No components, no colouring, no condensation, no visitor
! machinery: each is added when a caller requires it.
!
!                    SUBOBJECTS, NOT INTEGER LISTS
!
! Sources and sinks are returned as declared subobjects of the
! view's own domain - so they have identity, membership, size and
! local_index without additional operations, and their enumeration is
! CANONICAL BY THE DOMAIN'S DECLARATION ORDER: the scan traverses V by
! local index, and a member's numeric value never orders anything. An
! isolated member, with empty image and empty preimage, is both by
! definition.
!
!                    CONVENTIONS, FIXED
!
! reachable(v, v) is TRUE for every member v, by the zero-length
! path; an endpoint outside the domain returns FALSE, never an
! index into invalid storage. A topological order is undefined on a
! cycle: the algorithm REJECTS the input rather than constructing an
! order. A cyclic relation remains a valid structure - a valid
! interpretation is not a valid input to every algorithm.
!
!                    COSTS, AS WRITTEN
!
! Nothing here is optimized before a large caller exists. With n = |V|,
! m = |A|, and the carrier's own lookup cost T_idx:
!
!      sources/sinks       O(n) fibre reads, each costing T_idx -
!                          and then the subobject declaration,
!                          which validates every retained member
!                          against the ambient (T_has each) and
!                          removes duplicates with its current
!                          quadratic upper bound check
!      reachable           breadth-first, O(n + m) fibre reads and
!                          visited marks, each mark one T_idx
!      topological_order   the plain deterministic Kahn algorithm:
!                          each round rescans the domain for the
!                          first member of zero indegree in
!                          declaration order - O(n^2) scanning plus
!                          O(n + m) fibre and T_idx work. A priority
!                          structure is not yet required.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module relation_algorithms

  use graph_fractal           , only : graph
  use relation_finitary          , only : relation
  use relation_binary   , only : binary_relation, transposed_relation, transpose_of
  use map_set           , only : set_map
  use map_set_store     , only : set_store

  implicit none

  private
  public :: sources, sinks, reachable, topological_order

contains

  !===================================================================!
  ! The members with empty preimage, as a subobject of the domain, in
  ! the domain's own order.
  !===================================================================!

  subroutine sources(adjacency, sets, chosen)

    class(relation), target      , intent(in)    :: adjacency
    type(set_store)              , intent(inout) :: sets
    type(graph)              , intent(out)   :: chosen

    class(binary_relation), pointer :: a
    type(graph)      :: dom

    call require_adjacency(adjacency, a, dom)
    call declare_unpointed(a, 'sources', dom, sets, chosen)

  end subroutine sources

  !===================================================================!
  ! The members with empty image: the sources of the converse,
  ! read through the transposed view, so the search is written once.
  !===================================================================!

  subroutine sinks(adjacency, sets, chosen)

    class(relation), target      , intent(in)    :: adjacency
    type(set_store)              , intent(inout) :: sets
    type(graph)              , intent(out)   :: chosen

    class(binary_relation), pointer :: a
    type(transposed_relation), target :: converse
    type(graph)      :: dom

    call require_adjacency(adjacency, a, dom)
    converse = transpose_of(a)
    call declare_unpointed(converse, 'sinks', dom, sets, chosen)

  end subroutine sinks

  !===================================================================!
  ! The members with an empty preimage, declared as a subobject of the
  ! domain in the domain's own order.
  !===================================================================!

  subroutine declare_unpointed(a, label, dom, sets, chosen)

    class(binary_relation), target, intent(in)    :: a
    character(len=*)              , intent(in)    :: label
    type(graph)                   , intent(in)    :: dom
    type(set_store)               , intent(inout) :: sets
    type(graph)                   , intent(out)   :: chosen

    integer, allocatable :: retained(:)
    integer, pointer     :: fibre(:)
    integer              :: i, n, m, size_of_dom

    size_of_dom = sets % num_members_of(dom)

    allocate(retained(size_of_dom))
    n = 0
    do i = 1, size_of_dom
       m = sets % member_of(dom, i)
       fibre => a % preimage_view(m)
       if (size(fibre) == 0) then
          n = n + 1
          retained(n) = m
       end if
    end do

    call sets % declare_subobject(chosen, retained(1:n), label, dom)

  end subroutine declare_unpointed

  !===================================================================!
  ! Whether a directed path exists. Every member reaches itself by the
  ! zero-length path; a non-member reaches nothing and is reached by
  ! nothing. Breadth-first over the successor fibres, visited
  ! marked by local index.
  !===================================================================!

  logical function reachable(adjacency, sets, from, to)

    class(relation), target      , intent(in) :: adjacency
    type(set_map)                , intent(in) :: sets
    integer                      , intent(in) :: from
    integer                      , intent(in) :: to

    class(binary_relation), pointer :: a
    type(graph)      :: dom
    logical, allocatable :: visited(:)
    integer, allocatable :: queue(:)
    integer, pointer     :: fibre(:)
    integer              :: head, tail, v, j, s, n

    reachable = .false.

    call require_adjacency(adjacency, a, dom)
    if (.not. (sets % has(dom, from) .and. sets % has(dom, to))) return

    if (from == to) then
       reachable = .true.
       return
    end if

    n = sets % num_members_of(dom)
    allocate(visited(n), queue(n))
    visited = .false.

    head = 1
    tail = 1
    queue(1) = from
    visited(sets % index_in(dom, from)) = .true.

    do while (head <= tail)
       v = queue(head)
       head = head + 1
       fibre => a % image_view(v)
       do j = 1, size(fibre)
          s = fibre(j)
          if (s == to) then
             reachable = .true.
             return
          end if
          if (.not. visited(sets % index_in(dom, s))) then
             visited(sets % index_in(dom, s)) = .true.
             tail = tail + 1
             queue(tail) = s
          end if
       end do
    end do

  end function reachable

  !===================================================================!
  ! The deterministic Kahn algorithm: n rounds, each taking the FIRST
  ! member of zero indegree in the domain's declaration order -
  ! local_index, never numeric value - so one graph has one order.
  ! Members are returned as member values. A cycle leaves no member of
  ! zero indegree before the n rounds are complete, and the algorithm
  ! rejects the input.
  !===================================================================!

  subroutine topological_order(adjacency, sets, order, acyclic)

    class(relation), target      , intent(in)  :: adjacency
    type(set_map)                , intent(in)  :: sets
    integer, allocatable         , intent(out) :: order(:)
    logical, optional            , intent(out) :: acyclic

    class(binary_relation), pointer :: a
    type(graph)      :: dom
    integer, allocatable :: indegree(:)
    logical, allocatable :: placed(:)
    integer, pointer     :: fibre(:)
    integer              :: n, i, j, round, selected

    call require_adjacency(adjacency, a, dom)
    n   = sets % num_members_of(dom)

    allocate(indegree(n), placed(n), order(n))
    placed = .false.
    if (present(acyclic)) acyclic = .true.
    do i = 1, n
       fibre => a % preimage_view(sets % member_of(dom, i))
       indegree(i) = size(fibre)
    end do

    do round = 1, n
       selected = 0
       do i = 1, n
          if (.not. placed(i) .and. indegree(i) == 0) then
             selected = i
             exit
          end if
       end do
       if (selected == 0) then
          if (present(acyclic)) then
             acyclic = .false.
             order   = order(1:round - 1)
             return
          end if
          error stop 'relation_algorithms: a topological order needs an acyclic graph'
       end if

       placed(selected) = .true.
       order(round) = sets % member_of(dom, selected)

       fibre => a % image_view(sets % member_of(dom, selected))
       do j = 1, size(fibre)
          i = sets % index_in(dom, fibre(j))
          indegree(i) = indegree(i) - 1
       end do
    end do

  end subroutine topological_order

  !===================================================================!
  ! The precondition every algorithm checks: the adjacency must be a
  ! binary relation over one domain (source and target the same set).
  ! Either violation stops the program, because a traversal over two
  ! domains has no single member set to enumerate.
  !===================================================================!

  subroutine require_adjacency(adjacency, a, dom)

    class(relation), target        , intent(in)  :: adjacency
    class(binary_relation), pointer, intent(out) :: a
    type(graph)                , intent(out) :: dom

    type(graph) :: s, t

    select type (adjacency)
    class is (binary_relation)
       a => adjacency
    class default
       error stop 'relation_algorithms: the adjacency is a binary relation'
    end select

    s = a % source()
    t = a % target()
    if (.not. s % same_as(t)) then
       error stop 'relation_algorithms: the adjacency runs over one domain'
    end if

    dom = s

  end subroutine require_adjacency


end module relation_algorithms
