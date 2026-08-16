!=====================================================================!
! THE ORDINARY GRAPH, AS A VIEW
!
! One reading of the kernel graph: the ordinary binary graph the old
! solvers speak - vertices, edges, incidence, named sets,
! neighbourhoods. It is a VIEW and not an ontology, and the module
! name now says so where the type name still cannot.
!
! This module was carved out of graph_grammar, which had become the
! one place everything legacy met. Nothing here is new: the contract
! is moved verbatim, so that its consumers import the thing they
! actually use and the compatibility module can be drained to
! nothing (doc/final-codebase-cutover-plan.md, PR2).
!
! WHAT IT LENDS. One name: the abstract type. set_graph is imported
! to spell the signatures below and is NOT re-exported - a consumer
! wanting the kernel graph asks the kernel, which mints it. Importing
! a name to write a declaration is not the same act as lending it on.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!
!
!                      WHAT A GRAPH IS MADE OF
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

module graph_ordinary_view

  !===================================================================!
  ! THE DOMAIN IS A GRAPH, AND ITS INTERPRETATION IS THE CALLER'S.
  !
  ! set_graph is the kernel's graph, renamed on import for one reason
  ! only: this module and the kernel both speak of graphs, and a
  ! reader of a signature must be able to tell which. The COLLISION is
  ! gone - the abstract type below is `ordinary_graph` now, and no
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

  use fractal_graph       , only : set_graph => graph
  use graph_set_map       , only : set_map
  use graph_label_map     , only : label_map
  use graph_inclusion_map , only : inclusion_map

  implicit none

  private

  public :: ordinary_graph

  !===================================================================!
  ! GRAPH. The reader of structure.
  !
  ! Thirty-six symbols, all queries: identity, counts, incidence,
  ! named sets, neighbourhoods, and the frame. A graph answers; it
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
  ! MIGRATION DEBT, AND WHAT IS LEFT OF IT. This abstract type is
  ! named graph but is not the graph ontology. The ontology is
  ! G=(B1,B2) in src/fractal_graph.f90; this is a view over it. The
  ! module name has been fixed; the TYPE name has not, and that is the
  ! whole of the remaining debt. Do not present it as ontology.
  !
  ! The line this note replaces said "do not extend it", which read
  ! oddly beside class_graph, which does. It always meant: no NEW
  ! concretion. The existing one is the contract's reason to exist.
  !===================================================================!

  type, abstract :: ordinary_graph

   contains

     ! Identity and size.
     procedure(graph_id_interface)    , deferred :: id
     procedure(graph_count_interface) , deferred :: num_vertices
     procedure(graph_count_interface) , deferred :: num_edges

     ! The carrier bridge (migration, AGENTS.md 5B): the graph's two
     ! persistent declared domains, for consumers that must ask
     ! where a field domain ultimately lives. This root is already
     ! explicitly the ordinary vertex/edge compatibility contract.
     procedure(set_graph_interface), deferred :: vertex_set
     procedure(set_graph_interface), deferred :: edge_set

     ! Incidence: the two integer edge fields that ARE the structure.
     procedure(graph_edge_end_interface)     , deferred :: edge_tail
     procedure(graph_edge_end_interface)     , deferred :: edge_head
     procedure(graph_edge_has_head_interface), deferred :: edge_has_head

     ! The named sets, split by whether the answer already exists.
     !
     !   all_*        IS the carrier - stable identity, no binding
     !   the rest     CARVED on demand - a fresh set each call, so
     !                each call binds its extension, its label and
     !                its declared embedding into the caller's maps
     !
     ! The split is not cosmetic: the first kind may be asked twice
     ! and answer one set, the second kind answers two.
     procedure(set_graph_interface)   , deferred :: all_vertices
     procedure(graph_carved_set_interface), deferred :: interior_vertices
     procedure(graph_carved_set_interface), deferred :: boundary_vertices
     procedure(graph_tagged_set_interface), deferred :: tagged_vertices
     procedure(set_graph_interface)   , deferred :: all_edges
     procedure(graph_carved_set_interface), deferred :: interior_edges
     procedure(graph_carved_set_interface), deferred :: boundary_edges
     procedure(graph_tagged_set_interface), deferred :: tagged_edges

     ! The frame's sets, one part at a time.
     procedure(graph_part_set_interface), deferred :: owned_vertices
     procedure(graph_part_set_interface), deferred :: borrowed_vertices
     procedure(graph_part_set_interface), deferred :: overlap_vertices
     procedure(graph_part_set_interface), deferred :: owned_edges
     procedure(graph_part_set_interface), deferred :: borrowed_edges
     procedure(graph_part_set_interface), deferred :: overlap_edges

     ! Neighbourhoods. Called inside loops, so the answers are bare
     ! indices and the procedures are pure; handing back a graph here
     ! would allocate three times per neighbour query.
     procedure(graph_from_vertex_interface), deferred :: incident_edges
     procedure(graph_from_vertex_interface), deferred :: adjacent_vertices
     procedure(graph_from_vertex_interface), deferred :: outgoing_edges
     procedure(graph_from_vertex_interface), deferred :: incoming_edges
     procedure(graph_from_vertex_interface), deferred :: outgoing_vertices
     procedure(graph_from_vertex_interface), deferred :: incoming_vertices

     !----------------------------------------------------------------!
     ! LEGACY PARTITION FRAME - MOVES IN THE NEXT FRAME PR.
     !
     ! The frame's relations: counts, ownership, and the index maps
     ! both ways between a part and its whole. These eight are the
     ! only bindings on this type that are not a view question. Read
     ! for what they are, they are already spoken for elsewhere:
     !
     !     global_*_index    the EXTENSION of a subobject   set_map
     !     part_*_index      the same map read backwards    set_map
     !     has_part_relation, num_parts    provenance       inclusion_map
     !     *_owner_part      an integer field on the set    a field
     !
     ! They are carried here UNCHANGED, and deliberately: moving them
     ! in this phase would strand class_graph_partitioner and
     ! class_graph_assembler, which are their only two consumers
     ! outside the implementation. The six frame SETS above are NOT
     ! part of this - they carve and bind, which a view may do.
     !----------------------------------------------------------------!
     procedure(graph_count_interface)            , deferred :: num_parts
     procedure(graph_has_part_relation_interface), deferred :: has_part_relation
     procedure(graph_global_id_interface)        , deferred :: global_vertex_index
     procedure(graph_global_id_interface)        , deferred :: global_edge_index
     procedure(graph_part_id_interface)          , deferred :: part_vertex_index
     procedure(graph_part_id_interface)          , deferred :: part_edge_index
     procedure(graph_owner_part_interface)       , deferred :: vertex_owner_part
     procedure(graph_owner_part_interface)       , deferred :: edge_owner_part

  end type ordinary_graph

  abstract interface
     !===============================================================!
     ! Structure: identity, counts, incidence.
     !===============================================================!

     pure integer function graph_id_interface(this)
       import :: ordinary_graph
       class(ordinary_graph), intent(in) :: this
     end function graph_id_interface

     pure integer function graph_count_interface(this)
       import :: ordinary_graph
       class(ordinary_graph), intent(in) :: this
     end function graph_count_interface

     !---------------------------------------------------------------!
     ! A domain the graph already holds: identity, and nothing else.
     ! The extension is 1..num_vertices() or 1..num_edges(), which the
     ! caller can already read, so a counted representation is one
     ! constructor call away and no map need travel with the answer.
     !---------------------------------------------------------------!

     ! Not pure: a set graph carries a pointer component, so copying
     ! one out of an INTENT(IN) dummy is barred from a pure subprogram
     ! (F2018 C1594). Identity is still answered by value.
     type(set_graph) function set_graph_interface(this)
       import :: ordinary_graph, set_graph
       class(ordinary_graph), intent(in) :: this
     end function set_graph_interface

     pure integer function graph_edge_end_interface(this, edge_index)
       import :: ordinary_graph
       class(ordinary_graph), intent(in) :: this
       integer, intent(in) :: edge_index
     end function graph_edge_end_interface

     pure logical function graph_edge_has_head_interface(this, edge_index)
       import :: ordinary_graph
       class(ordinary_graph), intent(in) :: this
       integer, intent(in) :: edge_index
     end function graph_edge_has_head_interface

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

     subroutine graph_carved_set_interface(this, sets, labels, &
          & inclusions, members)
       import :: ordinary_graph, set_graph, set_map, label_map, inclusion_map
       class(ordinary_graph)       , intent(in)    :: this
       type(set_map)      , intent(inout) :: sets
       type(label_map)    , intent(inout) :: labels
       type(inclusion_map), intent(inout) :: inclusions
       type(set_graph)    , intent(out)   :: members
     end subroutine graph_carved_set_interface

     subroutine graph_tagged_set_interface(this, tag, sets, labels, &
          & inclusions, members)
       import :: ordinary_graph, set_graph, set_map, label_map, inclusion_map
       class(ordinary_graph)       , intent(in)    :: this
       character(len=*)   , intent(in)    :: tag
       type(set_map)      , intent(inout) :: sets
       type(label_map)    , intent(inout) :: labels
       type(inclusion_map), intent(inout) :: inclusions
       type(set_graph)    , intent(out)   :: members
     end subroutine graph_tagged_set_interface

     subroutine graph_part_set_interface(this, part_id, sets, labels, &
          & inclusions, members)
       import :: ordinary_graph, set_graph, set_map, label_map, inclusion_map
       class(ordinary_graph)       , intent(in)    :: this
       integer            , intent(in)    :: part_id
       type(set_map)      , intent(inout) :: sets
       type(label_map)    , intent(inout) :: labels
       type(inclusion_map), intent(inout) :: inclusions
       type(set_graph)    , intent(out)   :: members
     end subroutine graph_part_set_interface

     !===============================================================!
     ! Neighbourhoods. Bare indices, pure, loop-safe.
     !===============================================================!

     pure subroutine graph_from_vertex_interface(this, vertex_index, indices)
       import :: ordinary_graph
       class(ordinary_graph), intent(in) :: this
       integer, intent(in) :: vertex_index
       integer, allocatable, intent(out) :: indices(:)
     end subroutine graph_from_vertex_interface

     !===============================================================!
     ! The frame: how a part relates to the whole.
     !===============================================================!

     !---------------------------------------------------------------!
     ! Does this graph hold a partition record?
     !
     !      straight off a mesh file  ->  no
     !      out of a partitioner      ->  yes
     !
     ! An assembler checks this first. Without the relation there is
     ! no way back to whole-graph order, and this check reports that
     ! rather than assuming a map.
     !---------------------------------------------------------------!

     pure logical function graph_has_part_relation_interface(this)
       import :: ordinary_graph
       class(ordinary_graph), intent(in) :: this
     end function graph_has_part_relation_interface

     !---------------------------------------------------------------!
     ! What was this cell called in the whole graph?
     !
     !      whole    1   2   3   4   5   6   7   8
     !                       |   |           |
     !      part 2           1   2           3
     !
     !      global_vertex_index(2) = 4
     !
     ! Read that as: the part numbers this cell 2; the whole graph it
     ! was cut from numbered it 4. For an uncut graph the result
     ! equals the argument, because the graph is its own whole. The
     ! same map on an edgeless member set is what names the members.
     !---------------------------------------------------------------!

     pure integer function graph_global_id_interface(this, index)
       import :: ordinary_graph
       class(ordinary_graph), intent(in) :: this
       integer, intent(in) :: index
     end function graph_global_id_interface

     !---------------------------------------------------------------!
     ! The same map read backwards: where does whole-graph vertex 4
     ! sit inside part 2? The part is named because the partition is
     ! replicated - every image builds every part, so one image
     ! regularly queries a part it does not own.
     !---------------------------------------------------------------!

     pure integer function graph_part_id_interface(this, global_index, part_id)
       import :: ordinary_graph
       class(ordinary_graph), intent(in) :: this
       integer, intent(in) :: global_index
       integer, intent(in) :: part_id
     end function graph_part_id_interface

     !---------------------------------------------------------------!
     ! Whose job is this one?
     !
     !      vertex_owner_part(4) = 2       part 2 answers for it;
     !                                     everyone else only reads it
     !
     ! This is what stops a shared cell being counted twice when the
     ! parts are added back together.
     !---------------------------------------------------------------!

     pure integer function graph_owner_part_interface(this, index)
       import :: ordinary_graph
       class(ordinary_graph), intent(in) :: this
       integer, intent(in) :: index
     end function graph_owner_part_interface
  end interface

end module graph_ordinary_view
