!=====================================================================!
! LEVEL 4 OF THE NEW TOWER . THE GRAPH ALGORITHMS
!
! The level answers one question: WHAT GRAPH-THEORETIC QUESTIONS CAN
! BE ASKED of an interpretation. The algorithms live HERE, as free
! module procedures over the directed adjacency view - never as
! methods on relation, on the container, or on the view itself:
! traversal acts on structure, it is not structure (AGENTS.md 18,
! CALCULATOR.md 11). This module holds exactly what its first real
! caller - the calculator's dependency walk - has earned:
!
!      sources      the members nothing points at
!      sinks        the members that point at nothing
!      reachable    is there a directed path
!      topological_order   every arrow forward, or refusal
!
! No components, no colouring, no condensation, no visitor
! machinery: each waits for the caller that earns it.
!
!                    SUBOBJECTS, NOT INTEGER LISTS
!
! Sources and sinks are answered as subset_set subobjects of the
! view's own domain - so they carry identity, membership, size and
! local_index for free, and their enumeration is CANONICAL BY THE
! DOMAIN'S DECLARATION ORDER: the scan walks V by local index, and a
! member's numeric value never orders anything. An isolated member,
! pointing nowhere and pointed at by nothing, is lawfully both.
!
!                    CONVENTIONS, PINNED
!
! reachable(v, v) is TRUE for every member v, by the zero-length
! path; an endpoint outside the domain answers FALSE, never an
! index into invalid storage. A topological order is undefined on a
! cycle: the walk REFUSES loudly rather than inventing an order.
! The cyclic view itself remains a lawful directed graph - a valid
! interpretation is not a valid input to every algorithm.
!
!                    COSTS, AS WRITTEN
!
! Nothing here is optimized ahead of a large caller. With n = |V|,
! m = |A|, and the carrier's own lookup cost T_idx:
!
!      sources/sinks       O(n) fibre borrows, each paying T_idx
!      reachable           breadth-first, O(n + m) fibre reads and
!                          visited stamps, each stamp one T_idx
!      topological_order   the plain deterministic Kahn: each round
!                          rescans the domain for the first ready
!                          member in declaration order - O(n^2)
!                          scanning plus O(n + m) fibre and T_idx
!                          work. A priority structure is unearned.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module graph_algorithms

  use graph_carrier, only : member_set, subset_set
  use graph_profile, only : directed_adjacency_view

  implicit none

  private
  public :: sources, sinks, reachable, topological_order

contains

  !===================================================================!
  ! The members nothing points at, as a subobject of the domain, in
  ! the domain's own order.
  !===================================================================!

  function sources(view) result(chosen)

    type(directed_adjacency_view), intent(in) :: view
    type(subset_set)                          :: chosen

    class(member_set), allocatable :: dom
    integer, allocatable           :: keep(:)
    integer, pointer               :: fibre(:)
    integer                        :: i, n, m

    dom = view % domain()
    allocate(keep(dom % size()))
    n = 0
    do i = 1, dom % size()
       m = dom % member(i)
       fibre => view % predecessors_view(m)
       if (size(fibre) == 0) then
          n = n + 1
          keep(n) = m
       end if
    end do

    chosen = subset_set('sources', dom, keep(1:n))

  end function sources

  !===================================================================!
  ! The members that point at nothing, likewise.
  !===================================================================!

  function sinks(view) result(chosen)

    type(directed_adjacency_view), intent(in) :: view
    type(subset_set)                          :: chosen

    class(member_set), allocatable :: dom
    integer, allocatable           :: keep(:)
    integer, pointer               :: fibre(:)
    integer                        :: i, n, m

    dom = view % domain()
    allocate(keep(dom % size()))
    n = 0
    do i = 1, dom % size()
       m = dom % member(i)
       fibre => view % successors_view(m)
       if (size(fibre) == 0) then
          n = n + 1
          keep(n) = m
       end if
    end do

    chosen = subset_set('sinks', dom, keep(1:n))

  end function sinks

  !===================================================================!
  ! Is there a directed path. Every member reaches itself by the
  ! zero-length path; an outsider reaches nothing and is reached by
  ! nothing. Breadth-first over the successor fibres, visited
  ! stamped by local index.
  !===================================================================!

  logical function reachable(view, from, to)

    type(directed_adjacency_view), intent(in) :: view
    integer                      , intent(in) :: from
    integer                      , intent(in) :: to

    class(member_set), allocatable :: dom
    logical, allocatable           :: visited(:)
    integer, allocatable           :: queue(:)
    integer, pointer               :: fibre(:)
    integer                        :: head, tail, v, j, s

    reachable = .false.

    dom = view % domain()
    if (.not. (dom % has(from) .and. dom % has(to))) return

    if (from == to) then
       reachable = .true.
       return
    end if

    allocate(visited(dom % size()), queue(dom % size()))
    visited = .false.

    head = 1
    tail = 1
    queue(1) = from
    visited(dom % local_index(from)) = .true.

    do while (head <= tail)
       v = queue(head)
       head = head + 1
       fibre => view % successors_view(v)
       do j = 1, size(fibre)
          s = fibre(j)
          if (s == to) then
             reachable = .true.
             return
          end if
          if (.not. visited(dom % local_index(s))) then
             visited(dom % local_index(s)) = .true.
             tail = tail + 1
             queue(tail) = s
          end if
       end do
    end do

  end function reachable

  !===================================================================!
  ! The deterministic Kahn walk: n rounds, each taking the FIRST
  ! ready member in the domain's declaration order - local_index,
  ! never numeric value - so one graph has one answer. Members come
  ! back as member values. A cycle leaves no ready member before
  ! the walk is done, and the walk refuses.
  !===================================================================!

  subroutine topological_order(view, order)

    type(directed_adjacency_view), intent(in)  :: view
    integer, allocatable         , intent(out) :: order(:)

    class(member_set), allocatable :: dom
    integer, allocatable           :: indegree(:)
    logical, allocatable           :: placed(:)
    integer, pointer               :: fibre(:)
    integer                        :: n, i, j, round, pick

    dom = view % domain()
    n   = dom % size()

    allocate(indegree(n), placed(n), order(n))
    placed = .false.
    do i = 1, n
       fibre => view % predecessors_view(dom % member(i))
       indegree(i) = size(fibre)
    end do

    do round = 1, n
       pick = 0
       do i = 1, n
          if (.not. placed(i) .and. indegree(i) == 0) then
             pick = i
             exit
          end if
       end do
       if (pick == 0) then
          error stop 'graph_algorithms: a topological order needs an acyclic graph'
       end if

       placed(pick) = .true.
       order(round) = dom % member(pick)

       fibre => view % successors_view(dom % member(pick))
       do j = 1, size(fibre)
          i = dom % local_index(fibre(j))
          indegree(i) = indegree(i) - 1
       end do
    end do

  end subroutine topological_order

end module graph_algorithms
