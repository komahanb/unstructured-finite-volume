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
  use relation_binary   , only : integer_fibre, binary_relation, transposed_relation, transpose_of
  use map_set           , only : set_map
  use map_set_store     , only : set_store

  implicit none

  private
  public :: sources, sinks, reachable, topological_order

  ! Per-call traversal state. References remain valid only during the
  ! enclosing algorithm; module readers require no captured procedure.
  type :: reachability_context
     type(set_map), pointer :: sets => null()
     type(graph) :: domain
     logical, pointer :: visited(:) => null()
     integer, pointer :: queue(:) => null()
     integer :: tail = 0
     integer :: target_member = 0
     logical :: is_reachable = .false.
  end type reachability_context

  type :: topological_order_context
     type(set_map), pointer :: sets => null()
     type(graph) :: domain
     integer, pointer :: indegree(:) => null()
     ! the members of zero indegree not yet enumerated: a priority
     ! queue by local index in an array, whose first entry is the minimum
     integer, pointer :: ready(:) => null()
     integer :: num_ready = 0
  end type topological_order_context

contains

  !===================================================================!
  ! The members with empty preimage, as a subobject of the domain, in
  ! the domain's own order.
  !===================================================================!

  subroutine sources(adjacency, sets, source_set)

    class(relation), target      , intent(in)    :: adjacency
    type(set_store)              , intent(inout) :: sets
    type(graph)              , intent(out)   :: source_set

    class(binary_relation), pointer :: a
    type(graph)      :: dom

    call require_adjacency(adjacency, a, dom)
    call declare_unpointed(a, 'sources', dom, sets, source_set)

  end subroutine sources

  !===================================================================!
  ! The members with empty image: the sources of the converse,
  ! read through the transposed view, so the search is written once.
  !===================================================================!

  subroutine sinks(adjacency, sets, sink_set)

    class(relation), target      , intent(in)    :: adjacency
    type(set_store)              , intent(inout) :: sets
    type(graph)              , intent(out)   :: sink_set

    class(binary_relation), pointer :: a
    type(transposed_relation), target :: converse
    type(graph)      :: dom

    call require_adjacency(adjacency, a, dom)
    converse = transpose_of(a)
    call declare_unpointed(converse, 'sinks', dom, sets, sink_set)

  end subroutine sinks

  !===================================================================!
  ! The members with an empty preimage, declared as a subobject of the
  ! domain in the domain's own order.
  !===================================================================!

  subroutine declare_unpointed(a, label, dom, sets, unpointed_set)

    class(binary_relation), target, intent(in)    :: a
    character(len=*)              , intent(in)    :: label
    type(graph)                   , intent(in)    :: dom
    type(set_store)               , intent(inout) :: sets
    type(graph)                   , intent(out)   :: unpointed_set

    integer, allocatable :: retained(:)
    type(integer_fibre) :: fibre
    integer              :: i, n, m, size_of_dom

    size_of_dom = sets % num_members_of(dom)

    allocate(retained(size_of_dom))
    n = 0
    do i = 1, size_of_dom
       m = sets % member_of(dom, i)
       fibre = a % preimage_view(m)
       if (fibre % num_members() == 0) then
          n = n + 1
          retained(n) = m
       end if
    end do

    call sets % declare_subobject(unpointed_set, retained(1:n), label, dom)

  end subroutine declare_unpointed

  !===================================================================!
  ! Whether a directed path exists. Every member reaches itself by the
  ! zero-length path; a non-member reaches nothing and is reached by
  ! nothing. Breadth-first over the successor fibres, visited
  ! marked by local index.
  !===================================================================!

  logical function reachable(adjacency, sets, from, to)

    class(relation), target      , intent(in) :: adjacency
    type(set_map), target        , intent(in) :: sets
    integer                      , intent(in) :: from
    integer                      , intent(in) :: to

    class(binary_relation), pointer :: a
    type(reachability_context) :: context
    logical, allocatable, target :: visited(:)
    integer, allocatable, target :: queue(:)
    type(integer_fibre) :: fibre
    integer              :: head, v, n

    reachable = .false.

    call require_adjacency(adjacency, a, context % domain)
    if (.not. (sets % has(context % domain, from) .and. sets % has(context % domain, to))) return

    if (from == to) then
       reachable = .true.
       return
    end if

    n = sets % num_members_of(context % domain)
    allocate(visited(n), queue(n))
    visited = .false.
    context % sets => sets
    context % visited => visited
    context % queue => queue
    context % target_member = to

    head = 1
    context % tail = 1
    queue(1) = from
    visited(sets % index_in(context % domain, from)) = .true.

    do while (head <= context % tail)
       v = queue(head)
       head = head + 1
       fibre = a % image_view(v)
       call fibre % read(visit_successors, context)
       if (context % is_reachable) then
          reachable = .true.
          return
       end if
    end do

  end function reachable

  subroutine visit_successors(members, context)

    integer, intent(in) :: members(:)
    class(*), intent(inout) :: context
    integer :: j, s, member_index

    select type (context)
    type is (reachability_context)
      do j = 1, size(members)
         s = members(j)
         if (s == context % target_member) then
            context % is_reachable = .true.
            return
         end if
         member_index = context % sets % index_in(context % domain, s)
         if (.not. context % visited(member_index)) then
            context % visited(member_index) = .true.
            context % tail = context % tail + 1
            context % queue(context % tail) = s
         end if
      end do
    class default
       error stop 'relation_algorithms: visit_successors requires a reachability_context, but &
            &the actual context has a different dynamic type'
    end select

  end subroutine visit_successors

  !===================================================================!
  ! The deterministic Kahn algorithm: n rounds, each taking the FIRST
  ! member of zero indegree in the domain's declaration order -
  ! local_index, never numeric value - so one graph has one order.
  ! The members of zero indegree form a priority queue by local index,
  ! so a round costs log n and the order (n + e) log n rather than n^2.
  ! Members are returned as member values. A cycle leaves no member of
  ! zero indegree before the n rounds are complete, and the algorithm
  ! rejects the input.
  !===================================================================!

  subroutine topological_order(adjacency, sets, order, acyclic)

    class(relation), target      , intent(in)  :: adjacency
    type(set_map), target        , intent(in)  :: sets
    integer, allocatable         , intent(out) :: order(:)
    logical, optional            , intent(out) :: acyclic

    class(binary_relation), pointer :: a
    type(topological_order_context) :: context
    integer, allocatable, target :: indegree(:), ready(:)
    type(integer_fibre) :: fibre
    integer              :: n, i, order_index, selected
    character(len=150) :: message

    call require_adjacency(adjacency, a, context % domain)
    n   = sets % num_members_of(context % domain)

    allocate(indegree(n), ready(n), order(n))
    context % sets => sets
    context % indegree => indegree
    context % ready => ready
    if (present(acyclic)) acyclic = .true.
    do i = 1, n
       fibre = a % preimage_view(sets % member_of(context % domain, i))
       indegree(i) = fibre % num_members()
    end do
    do i = 1, n
       if (indegree(i) == 0) call push_ready(context, i)
    end do

    do order_index = 1, n
       if (context % num_ready == 0) then
          if (present(acyclic)) then
             acyclic = .false.
             order   = order(1:order_index - 1)
             return
          end if
          write(message,'(a,i0,a,i0,a)') 'relation_algorithms: topological_order requires an &
               &acyclic graph; no member of zero indegree remains after ', order_index - 1, &
               & ' of ', n, ' rounds'
          error stop trim(message)
       end if

       selected = pop_ready(context)
       order(order_index) = sets % member_of(context % domain, selected)

       fibre = a % image_view(sets % member_of(context % domain, selected))
       call fibre % read(decrease_indegrees, context)
    end do

  end subroutine topological_order

  ! The priority queue of members of zero indegree, a binary tree in an
  ! array: the parent of position at is at / 2, and no parent exceeds
  ! its children.
  subroutine push_ready(context, member_index)

    type(topological_order_context), intent(inout) :: context
    integer, intent(in) :: member_index
    integer :: at, parent, exchanged

    context % num_ready = context % num_ready + 1
    at = context % num_ready
    context % ready(at) = member_index
    do while (at > 1)
       parent = at / 2
       if (context % ready(parent) <= context % ready(at)) exit
       exchanged = context % ready(parent)
       context % ready(parent) = context % ready(at)
       context % ready(at) = exchanged
       at = parent
    end do

  end subroutine push_ready

  integer function pop_ready(context) result(member_index)

    type(topological_order_context), intent(inout) :: context
    integer :: at, child, exchanged

    member_index = context % ready(1)
    context % ready(1) = context % ready(context % num_ready)
    context % num_ready = context % num_ready - 1
    at = 1
    do
       child = 2 * at
       if (child > context % num_ready) exit
       if (child < context % num_ready) then
          if (context % ready(child + 1) < context % ready(child)) child = child + 1
       end if
       if (context % ready(at) <= context % ready(child)) exit
       exchanged = context % ready(at)
       context % ready(at) = context % ready(child)
       context % ready(child) = exchanged
       at = child
    end do

  end function pop_ready

  subroutine decrease_indegrees(members, context)

    integer, intent(in) :: members(:)
    class(*), intent(inout) :: context
    integer :: j, member_index

    select type (context)
    type is (topological_order_context)
      do j = 1, size(members)
         member_index = context % sets % index_in(context % domain, members(j))
         context % indegree(member_index) = context % indegree(member_index) - 1
         if (context % indegree(member_index) == 0) call push_ready(context, member_index)
      end do
    class default
       error stop 'relation_algorithms: decrease_indegrees requires a topological_order_context, &
            &but the actual context has a different dynamic type'
    end select

  end subroutine decrease_indegrees

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
       error stop 'relation_algorithms: require_adjacency requires a binary_relation, but the &
            &actual adjacency has a different dynamic type'
    end select

    s = a % source()
    t = a % target()
    if (.not. s % same_as(t)) then
       error stop 'relation_algorithms: require_adjacency requires source and target to be the &
            &same domain, but they differ'
    end if

    dom = s

  end subroutine require_adjacency


end module relation_algorithms
