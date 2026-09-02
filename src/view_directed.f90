!=====================================================================!
! THE DIRECTED GRAPH, AS A VIEW
!
! One view of the kernel graph:
!
!            D = ( V, E, tail, head )
!
!            V              vertex set identity, and its count
!            E              edge set identity, and its count
!            tail, head     E -> V
!
! and its transpose D^T = ( V, E, head, tail ) is the same view read in
! the reverse orientation: a stored graph stores both orientations and
! switches between them without rebuilding an edge, so that
! (D^T)^T = D exactly.
!
! That is the whole of the role. `directed` is what the structure IS -
! two finite domains and two maps between them - and it is a view over
! the ontology, never a kind of graph.
!
! Named vertex and edge subsets are not extra
! structure: each is a subobject of V or of E, returned as a set graph
! identity and described in the caller's set store. Neighbourhood queries
! are compositions of tail and head, materialized because they are
! called inside loops.
!
! This module was split out of graph_grammar, which had become the
! one module every legacy dependency passed through
! (doc/final-codebase-cutover-plan.md, PR2). The type was renamed from
! `graph`, and then from `ordinary_` to `directed_`, because `ordinary`
! names no mathematical role and `directed` names the role this
! contract has.
!
! WHAT THIS MODULE EXPORTS. One name: the abstract type. graph is
! imported to write the signatures below and is NOT re-exported - a
! consumer that requires the kernel graph imports it from the kernel,
! which defines it. Importing a name to write a declaration is not the
! same as re-exporting it.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!
!
!                      WHAT A GRAPH IS MADE OF
!
! A vertex is a member of V. An edge joins two vertices, tail to head.
!
!                            e
!                     i ----------> j       edge_tail(e) = i
!                                           edge_head(e) = j
!                                           edge_has_head(e) = .true.
!
!                            b
!                     i ----------o         edge_tail(b) = i
!                                           edge_has_head(b) = .false.
!
! The second edge is attached to vertex i alone. That is a boundary
! face, and this is how it is written without introducing a fictitious
! cell beyond the boundary.
!
!
!=====================================================================!
!
!                        CAN A GRAPH CHANGE?
!
! CAN A GRAPH CHANGE? No. Everything a graph stores - structure,
! tags, its relation to the whole it was partitioned from - is set at
! construction, and no procedure below accepts data afterwards.
! When an operation computes a new quantity, the result is returned
! through that operation's output argument. The reason is
! repeatability: the same query on a graph evaluated twice returns the
! same value twice, whatever executed between the two evaluations.
!
!=====================================================================!

module view_directed

  !===================================================================!
  ! THE DOMAIN IS A GRAPH, AND ITS INTERPRETATION IS THE CALLER'S.
  !
  ! graph is the kernel's graph, renamed on import for one reason
  ! only: this module and the kernel both define graph types, and a
  ! reader of a signature must be able to tell which. The COLLISION is
  ! gone - the abstract type below is `directed_graph` now, and no
  ! other module exports a type called `graph` to any module that
  ! uses this one. What remains is a convenience rename, not a
  ! disambiguation.
  !
  ! A domain-producing symbol here returns WHICH set. Where the result
  ! is a set the graph already stores, the identity is all the symbol
  ! returns, and the caller reconstructs the extension from a count it
  ! can already read. Where the result is a subset declared on demand,
  ! the symbol writes its extension, name and embedding into one set
  ! store. The graph does not own those associated data, and callers
  ! no longer pass three maps through every signature.
  !===================================================================!

  use graph_fractal       , only : graph
  use map_set_store , only : set_store

  implicit none

  private

  public :: directed_graph
  public :: forward, reverse

  ! the two orientations a directed view is read in: along its edges,
  ! tail before head, or against them, which is reading its transpose
  integer, parameter :: forward = 1
  integer, parameter :: reverse = 2
  public :: SIDE_VERTEX
  public :: SIDE_EDGE

  !===================================================================!
  ! The two sides of a directed graph an operation's output may be
  ! defined on. Output-side identity only - not field-domain
  ! identity and not a subset: domains are set graph identities.
  !===================================================================!

  integer, parameter :: SIDE_VERTEX = 1
  integer, parameter :: SIDE_EDGE   = 2

  !===================================================================!
  ! GRAPH. The structure query interface.
  !
  ! Fifteen symbols, all queries: identity, counts, incidence, the
  ! tagged edge set, and neighbourhoods. A graph returns values; it
  ! performs no algorithm. Algorithms are applied to it from the
  ! levels above, which is what keeps this contract small.
  !
  ! THE GRAPH STORES NO VALUES. A field references its domain; the
  ! reference never points the other way. What an operation reads is
  ! passed at construction, as a field argument the compiler can
  ! check - a name passed as a string would defer the same binding to
  ! run time. Vocabulary that names particular data (a cell volume,
  ! a face normal) belongs to the level that defines those quantities,
  ! as typed procedures on its concretes, never as string keys here.
  ! The one string below is the tag, and it is data, not a symbol:
  ! it originates outside the code, in the mesh file that named its
  ! boundary groups.
  !
  ! A NAMED SET IS A SET GRAPH. The whole sets are the graph's own
  ! carriers - one stable identity, queried twice, returning the same
  ! identity; the subsets declared on demand return a NEW identity and
  ! bind its extension into the caller's set store,
  !
  !      vertex_set             tagged_edges('edge')
  !      the vertex carrier     a new set { 11 14 19 } c--> edges
  !
  ! and membership, size, order and position are queries on the
  ! representation the caller stores - not on the graph, which records
  ! only which set it named.
  !
  ! A part graph is still a graph. The part graph stores the relation
  ! to the whole - how many parts, which part owns what, and the index
  ! maps both ways - because an assembler must read that relation
  ! rather than construct one.
  !===================================================================!

  !===================================================================!
  ! THE HIERARCHY, STATED ONCE. Three types, three roles, three names:
  !
  !     graph_fractal :: graph                 the ontology, G=(B1,B2)
  !     view_directed :: directed_graph  this contract, D
  !     view_directed_stored :: stored_directed_graph   one stored realization
  !
  ! The migration state that made this type share the ontology's name
  ! is resolved. What remains of it is the module name view_directed_stored,
  ! which states less than the type inside it does.
  !
  ! Do not present this as ontology, and add no NEW concretion: one
  ! realization is what a contract needs to be instantiated.
  !===================================================================!

  type, abstract :: directed_graph

   contains

     ! Identity and size.
     procedure(directed_id_interface)    , deferred :: id
     procedure(directed_count_interface) , deferred :: num_vertices
     procedure(directed_count_interface) , deferred :: num_edges

     ! The carrier map (migration, AGENTS.md 5B): the graph's two
     ! persistent declared domains, for consumers that must determine
     ! which set a field domain is defined on. This root is already
     ! explicitly the directed vertex/edge contract: V and E, by
     ! identity.
     procedure(member_set_interface), deferred :: vertex_set
     procedure(member_set_interface), deferred :: edge_set

     ! Incidence: the two integer edge fields that ARE the structure.
     procedure(directed_edge_end_interface)     , deferred :: edge_tail
     procedure(directed_edge_end_interface)     , deferred :: edge_head
     procedure(directed_edge_has_head_interface), deferred :: edge_has_head

     ! The one named subset, declared on demand: a new set each call,
     ! so each call binds its extension, its label and its declared
     ! embedding into the caller's set store - called twice, it
     ! returns two sets. The whole vertex and edge sets are the
     ! carriers above, vertex_set and edge_set: stable identities, no
     ! binding.
     procedure(directed_tagged_set_interface), deferred :: tagged_edges

     ! Neighbourhoods. Called inside loops, so the results are bare
     ! indices and the procedures are pure; returning a graph here
     ! would allocate three times per neighbour query.
     procedure(directed_from_vertex_interface), deferred :: incident_edges
     procedure(directed_from_vertex_interface), deferred :: adjacent_vertices
     procedure(directed_from_vertex_interface), deferred :: outgoing_edges
     procedure(directed_from_vertex_interface), deferred :: incoming_edges
     procedure(directed_from_vertex_interface), deferred :: outgoing_vertices
     procedure(directed_from_vertex_interface), deferred :: incoming_vertices

     !----------------------------------------------------------------!
     ! THE PARTITION RELATION IS NOT HERE. Its queries - how many
     ! parts, which part owns what, the maps each way, the owned,
     ! halo and overlap subsets - are queries on r <= S_part x S_whole,
     ! not on D = (V, E, tail, head). They are partition_relation: a
     ! value the partitioner writes, a graph stores, and a caller
     ! queries directly.
     !----------------------------------------------------------------!

  end type directed_graph

  abstract interface
     !===============================================================!
     ! Structure: identity, counts, incidence.
     !===============================================================!

     pure integer function directed_id_interface(this)
       import :: directed_graph
       class(directed_graph), intent(in) :: this
     end function directed_id_interface

     pure integer function directed_count_interface(this)
       import :: directed_graph
       class(directed_graph), intent(in) :: this
     end function directed_count_interface

     !---------------------------------------------------------------!
     ! A domain the graph already stores: identity, and nothing else.
     ! The extension is 1..num_vertices() or 1..num_edges(), which the
     ! caller can already read, so a counted representation is one
     ! constructor call, and no map need be returned with the result.
     !---------------------------------------------------------------!

     ! Not pure: a set graph contains a pointer component, so copying
     ! one out of an INTENT(IN) dummy is barred from a pure subprogram
     ! (F2018 C1594). Identity is still returned by value.
     type(graph) function member_set_interface(this)
       import :: directed_graph, graph
       class(directed_graph), intent(in) :: this
     end function member_set_interface

     pure integer function directed_edge_end_interface(this, edge_index)
       import :: directed_graph
       class(directed_graph), intent(in) :: this
       integer, intent(in) :: edge_index
     end function directed_edge_end_interface

     pure logical function directed_edge_has_head_interface(this, edge_index)
       import :: directed_graph
       class(directed_graph), intent(in) :: this
       integer, intent(in) :: edge_index
     end function directed_edge_has_head_interface

     !===============================================================!
     ! THE DECLARED SUBSET. Called once, when an operation begins, so the
     ! cost is incurred per sweep and not per cell.
     !
     ! Each call declares a NEW set - a new identity: a subset
     ! declares its own identity, so two calls to tagged_edges() are
     ! two domains.
     !
     ! What the result needs beyond identity, the call binds in the
     ! caller's set store: listed extension, name and embedding into
     ! the graph's own carrier. The store is referenced for the
     ! duration of the call and never retained.
     !===============================================================!

     subroutine directed_tagged_set_interface(this, tag, sets, members)
       import :: directed_graph, graph, set_store
       class(directed_graph)       , intent(in)    :: this
       character(len=*)   , intent(in)    :: tag
       type(set_store)    , intent(inout) :: sets
       type(graph)    , intent(out)   :: members
     end subroutine directed_tagged_set_interface

     !===============================================================!
     ! Neighbourhoods. Bare indices, pure, loop-safe.
     !===============================================================!

     pure subroutine directed_from_vertex_interface(this, vertex_index, indices)
       import :: directed_graph
       class(directed_graph), intent(in) :: this
       integer, intent(in) :: vertex_index
       integer, allocatable, intent(out) :: indices(:)
     end subroutine directed_from_vertex_interface
  end interface

end module view_directed
