!=====================================================================!
! THE DIRECTED GRAPH, AS A VIEW
!
! One reading of the kernel graph:
!
!            D = ( V, E, tail, head )
!
!            V              vertex set identity, and its count
!            E              edge set identity, and its count
!            tail, head     E -> V
!
! That is the whole of the role. `directed` is what the structure IS -
! two finite domains and two maps between them - and it is a view over
! the ontology, never a kind of graph.
!
! Named, carved and tagged vertex and edge sets are not extra
! structure: each is a subobject of V or of E, answered as a set graph
! identity and described in the caller's maps. Neighbourhood queries
! are compositions of tail and head, materialized because they are
! asked inside loops.
!
! This module was carved out of graph_grammar, which had become the
! one place everything legacy met (doc/final-codebase-cutover-plan.md,
! PR2). The type was renamed from `graph`, and then from `ordinary_`
! to `directed_`, because `ordinary` names no mathematical role and
! `directed` names the one this contract actually holds.
!
! WHAT IT LENDS. One name: the abstract type. graph is imported
! to spell the signatures below and is NOT re-exported - a consumer
! wanting the kernel graph asks the kernel, which mints it. Importing
! a name to write a declaration is not the same act as lending it on.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!
!
!                      WHAT A GRAPH IS MADE OF
!
! A vertex is a thing. An edge joins two of them, tail to head.
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
! face, and this is how it is written without inventing an imaginary
! cell on the far side of the wall.
!
!
!=====================================================================!
!
!                        CAN A GRAPH CHANGE?
!
! CAN A GRAPH CHANGE? No. Everything a graph holds - structure,
! tags, its relation to the whole it came from - goes in at
! construction, and no procedure below accepts data afterwards.
! When an operation computes something new, the result leaves through
! that operation's output argument. The reason is repeatability: ask
! a graph the same question twice and it gives the same answer twice,
! no matter what ran in between.
!
!=====================================================================!

module view_directed

  !===================================================================!
  ! THE DOMAIN IS A GRAPH, AND ITS INTERPRETATION IS THE CALLER'S.
  !
  ! graph is the kernel's graph, renamed on import for one reason
  ! only: this module and the kernel both speak of graphs, and a
  ! reader of a signature must be able to tell which. The COLLISION is
  ! gone - the abstract type below is `directed_graph` now, and no
  ! other module lends a type called `graph` to anyone who reads this
  ! one. What remains is a convenience rename, not a disambiguation.
  !
  ! A domain-producing symbol here answers WHICH set. Where the answer
  ! is a set the graph already holds, that is all it answers, and the
  ! caller reconstructs the extension from a count it can already read.
  ! Where the answer is a set CARVED on demand, the symbol must also
  ! say how its members are stored, what it is called, and what it was
  ! carved from - so it takes the three maps and binds into them. The
  ! graph does not own them. It never learns what a set means.
  !===================================================================!

  use graph_fractal       , only : graph
  use map_set       , only : set_map
  use map_label     , only : label_map
  use map_inclusion , only : inclusion_map

  implicit none

  private

  public :: directed_graph
  public :: SIDE_VERTEX
  public :: SIDE_EDGE

  !===================================================================!
  ! The two sides of a directed graph an operation's output may
  ! land on. Output-landing identity only - not field-domain
  ! identity and not a subset: domains are set graph identities.
  !===================================================================!

  integer, parameter :: SIDE_VERTEX = 1
  integer, parameter :: SIDE_EDGE   = 2

  !===================================================================!
  ! GRAPH. The reader of structure.
  !
  ! Thirty-six symbols, all queries: identity, counts, incidence,
  ! named sets, and neighbourhoods. A graph answers; it
  ! never acts. Algorithms act on it from the levels above, which is
  ! what keeps this contract small.
  !
  ! THE GRAPH CARRIES NO VALUES. A field references its domain; the
  ! reference never points the other way. What an operation reads it
  ! is handed at construction, as a field argument the compiler can
  ! see - a name passed as a string would hide the same binding until
  ! run time. Vocabulary that names particular data (a cell volume,
  ! a face normal) belongs to the level that owns those words, as
  ! typed procedures on its concretes, never as string keys here.
  ! The one string below is the tag, and it is data, not a symbol:
  ! it originates outside the code, in the mesh file that named its
  ! boundary groups.
  !
  ! A NAMED SET IS A SET GRAPH. The full sets answer the graph's own
  ! carrier - one stable identity, asked twice, answering once; the
  ! carved ones answer a FRESH identity and bind what it means into
  ! the caller's maps,
  !
  !      all_vertices           tagged_edges('wall')
  !      the vertex carrier     a new set { 11 14 19 } c--> edges
  !
  ! and membership, size, order and standing are questions for the
  ! representation the caller holds - not for the graph, which knows
  ! only which set it named.
  !
  ! THE FRAME. How a part relates to the whole it was cut from:
  !
  !    owned      this part answers for the value here
  !    borrowed   this part only reads it; someone else owns it
  !    overlap    everything this part must see to finish what it owns
  !
  !            part 1                        part 2
  !       +---------------+            +---------------+
  !       |  o    o    o  |            |  o    o    o  |
  !       |  o    o    o--|------------|--b    o    o  |
  !       +---------------+            +---------------+
  !                    \______________/
  !                       the overlap of part 1
  !
  ! A part graph is still a graph. It holds the relation back to the
  ! whole - how many parts, which part owns what, and the index maps
  ! both ways - because an assembler must read that relation rather
  ! than invent one.
  !===================================================================!

  !===================================================================!
  ! THE HIERARCHY, STATED ONCE. Three types, three roles, three names:
  !
  !     graph_fractal :: graph                 the ontology, G=(B1,B2)
  !     view_directed :: directed_graph  this contract, D
  !     view_directed_stored :: stored_directed_graph   one stored realization
  !
  ! The migration debt that made this type share the ontology's name
  ! is discharged. What remains of it is the module name view_directed_stored,
  ! which says less than the type inside it does.
  !
  ! Do not present this as ontology, and add no NEW concretion: one
  ! realization is what a contract needs to be inhabited.
  !===================================================================!

  type, abstract :: directed_graph

   contains

     ! Identity and size.
     procedure(directed_id_interface)    , deferred :: id
     procedure(directed_count_interface) , deferred :: num_vertices
     procedure(directed_count_interface) , deferred :: num_edges

     ! The carrier bridge (migration, AGENTS.md 5B): the graph's two
     ! persistent declared domains, for consumers that must ask
     ! where a field domain ultimately lives. This root is already
     ! explicitly the directed vertex/edge contract: V and E, by
     ! identity.
     procedure(member_set_interface), deferred :: vertex_set
     procedure(member_set_interface), deferred :: edge_set

     ! Incidence: the two integer edge fields that ARE the structure.
     procedure(directed_edge_end_interface)     , deferred :: edge_tail
     procedure(directed_edge_end_interface)     , deferred :: edge_head
     procedure(directed_edge_has_head_interface), deferred :: edge_has_head

     ! The named sets, split by whether the answer already exists.
     !
     !   all_*        IS the carrier - stable identity, no binding
     !   the rest     CARVED on demand - a fresh set each call, so
     !                each call binds its extension, its label and
     !                its declared embedding into the caller's maps
     !
     ! The split is not cosmetic: the first kind may be asked twice
     ! and answer one set, the second kind answers two.
     procedure(member_set_interface)   , deferred :: all_vertices
     procedure(directed_carved_set_interface), deferred :: interior_vertices
     procedure(directed_carved_set_interface), deferred :: boundary_vertices
     procedure(directed_tagged_set_interface), deferred :: tagged_vertices
     procedure(member_set_interface)   , deferred :: all_edges
     procedure(directed_carved_set_interface), deferred :: interior_edges
     procedure(directed_carved_set_interface), deferred :: boundary_edges
     procedure(directed_tagged_set_interface), deferred :: tagged_edges

     ! Carved by ownership, one part at a time.
     procedure(directed_part_set_interface), deferred :: owned_vertices
     procedure(directed_part_set_interface), deferred :: borrowed_vertices
     procedure(directed_part_set_interface), deferred :: overlap_vertices
     procedure(directed_part_set_interface), deferred :: owned_edges
     procedure(directed_part_set_interface), deferred :: borrowed_edges
     procedure(directed_part_set_interface), deferred :: overlap_edges

     ! Neighbourhoods. Called inside loops, so the answers are bare
     ! indices and the procedures are pure; handing back a graph here
     ! would allocate three times per neighbour query.
     procedure(directed_from_vertex_interface), deferred :: incident_edges
     procedure(directed_from_vertex_interface), deferred :: adjacent_vertices
     procedure(directed_from_vertex_interface), deferred :: outgoing_edges
     procedure(directed_from_vertex_interface), deferred :: incoming_edges
     procedure(directed_from_vertex_interface), deferred :: outgoing_vertices
     procedure(directed_from_vertex_interface), deferred :: incoming_vertices

     !----------------------------------------------------------------!
     ! THE PARTITION RELATION IS GONE FROM HERE, and this note stands
     ! in its place so nobody puts it back.
     !
     ! Eight bindings sat at this point - num_parts,
     ! has_part_relation, the two maps each way, and the two ownership
     ! questions. Not one of them is a question about
     ! D = (V, E, tail, head): they are the tuples of
     ! r <= S_part x S_whole read both ways, r's provenance, and an
     ! integer field on the part's members. They are
     ! partition_relation now - a value the cut writes, a graph
     ! carries, and the four verbs are handed.
     !
     ! The six SETS above stayed. owned, borrowed and overlap carve a
     ! subobject of V or of E and bind what it means into the caller's
     ! maps, which is a view question and always was.
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
     ! A domain the graph already holds: identity, and nothing else.
     ! The extension is 1..num_vertices() or 1..num_edges(), which the
     ! caller can already read, so a counted representation is one
     ! constructor call away and no map need travel with the answer.
     !---------------------------------------------------------------!

     ! Not pure: a set graph carries a pointer component, so copying
     ! one out of an INTENT(IN) dummy is barred from a pure subprogram
     ! (F2018 C1594). Identity is still answered by value.
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
     ! THE CARVED SETS. Called once, when an operation begins, so the
     ! cost is paid per sweep and not per cell.
     !
     ! Each call declares a NEW set - a fresh identity - because that
     ! is what these have always done: the old code built a fresh
     ! subset_set per call, and a subset signs its own identity. Two
     ! calls to boundary_vertices() were never one domain, and are not
     ! one domain now.
     !
     ! What the answer needs beyond identity, it binds:
     !
     !     sets         the listed extension - who belongs
     !     labels       the name the old subset carried
     !     inclusions   the embedding into the graph's own carrier,
     !                  without which subobject questions go silent
     !
     ! All three are the CALLER'S, borrowed for the duration of the
     ! call and never stored. That is the difference between stating a
     ! dependency and hiding one.
     !===============================================================!

     subroutine directed_carved_set_interface(this, sets, labels, &
          & inclusions, members)
       import :: directed_graph, graph, set_map, label_map, inclusion_map
       class(directed_graph)       , intent(in)    :: this
       type(set_map)      , intent(inout) :: sets
       type(label_map)    , intent(inout) :: labels
       type(inclusion_map), intent(inout) :: inclusions
       type(graph)    , intent(out)   :: members
     end subroutine directed_carved_set_interface

     subroutine directed_tagged_set_interface(this, tag, sets, labels, &
          & inclusions, members)
       import :: directed_graph, graph, set_map, label_map, inclusion_map
       class(directed_graph)       , intent(in)    :: this
       character(len=*)   , intent(in)    :: tag
       type(set_map)      , intent(inout) :: sets
       type(label_map)    , intent(inout) :: labels
       type(inclusion_map), intent(inout) :: inclusions
       type(graph)    , intent(out)   :: members
     end subroutine directed_tagged_set_interface

     subroutine directed_part_set_interface(this, part_id, sets, labels, &
          & inclusions, members)
       import :: directed_graph, graph, set_map, label_map, inclusion_map
       class(directed_graph)       , intent(in)    :: this
       integer            , intent(in)    :: part_id
       type(set_map)      , intent(inout) :: sets
       type(label_map)    , intent(inout) :: labels
       type(inclusion_map), intent(inout) :: inclusions
       type(graph)    , intent(out)   :: members
     end subroutine directed_part_set_interface

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
