!=====================================================================!
! The graph, and the graph that has no edges.
!
! A graph is a set of vertices and the edges between them. Everything
! structural in this library is one. A support is the case with no
! edges at all: pure membership, which is why it answers one extra
! question - which side of its host graph its members name.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module structure_graph

  implicit none

  private

  public :: graph
  public :: graph_support
  public :: GRAPH_SIDE_VERTEX
  public :: GRAPH_SIDE_EDGE

  ! Which side of its host a support references. A finite axis,
  ! absorbed: it rides as an answer, never as a pair of types.
  integer, parameter :: GRAPH_SIDE_VERTEX = 1
  integer, parameter :: GRAPH_SIDE_EDGE   = 2

  !===================================================================!
  ! GRAPH. The reader of structure.
  !
  ! Thirty-four symbols, all queries: identity, counts, incidence,
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
  ! A NAMED SET IS ITSELF A GRAPH: the edgeless graph of its members,
  ! whose global indices point back into the graph that named them.
  !
  !      all_vertices           tagged_edges('wall')
  !      +-------------+        +-------------+
  !      | 1 2 3 4 5 6 |        |  11  14  19 |     no edges:
  !      +-------------+        +-------------+     pure membership
  !
  ! Level 1 names this degenerate case the support; here it needs no
  ! name, because an edgeless graph already answers every question a
  ! membership list is ever asked: its size is num_vertices, its
  ! members are its global indices.
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

  type, abstract :: graph

   contains

     ! Identity and size.
     procedure(graph_id_interface)    , deferred :: id
     procedure(graph_count_interface) , deferred :: num_vertices
     procedure(graph_count_interface) , deferred :: num_edges

     ! Incidence: the two integer edge fields that ARE the structure.
     procedure(graph_edge_end_interface)     , deferred :: edge_tail
     procedure(graph_edge_end_interface)     , deferred :: edge_head
     procedure(graph_edge_has_head_interface), deferred :: edge_has_head

     ! The named sets. Each answer is itself a graph: the edgeless
     ! graph of the members.
     procedure(graph_member_set_interface), deferred :: all_vertices
     procedure(graph_member_set_interface), deferred :: interior_vertices
     procedure(graph_member_set_interface), deferred :: boundary_vertices
     procedure(graph_tagged_set_interface), deferred :: tagged_vertices
     procedure(graph_member_set_interface), deferred :: all_edges
     procedure(graph_member_set_interface), deferred :: interior_edges
     procedure(graph_member_set_interface), deferred :: boundary_edges
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

     ! The frame's relations: counts, ownership, and the index maps
     ! both ways between a part and its whole.
     procedure(graph_count_interface)            , deferred :: num_parts
     procedure(graph_has_part_relation_interface), deferred :: has_part_relation
     procedure(graph_global_id_interface)        , deferred :: global_vertex_index
     procedure(graph_global_id_interface)        , deferred :: global_edge_index
     procedure(graph_part_id_interface)          , deferred :: part_vertex_index
     procedure(graph_part_id_interface)          , deferred :: part_edge_index
     procedure(graph_owner_part_interface)       , deferred :: vertex_owner_part
     procedure(graph_owner_part_interface)       , deferred :: edge_owner_part

  end type graph

  !===================================================================!
  ! GRAPH_SUPPORT. The grammar's graph at edge count zero: pure
  ! membership, no incidence.
  !
  !      all_vertices           tagged_edges('wall')
  !      +-------------+        +-------------+
  !      | 1 2 3 4 5 6 |        |  11  14  19 |
  !      +-------------+        +-------------+
  !
  ! Almost everything a membership list is asked, the inherited
  ! grammar already answers: its size is num_vertices, its members
  ! are its global indices, its identity is its id. One question has
  ! no inherited answer, and it is the one symbol added here: which
  ! SIDE of the host graph the members reference - its vertices or
  ! its edges. An index list cannot say this about itself, and the
  ! graph no longer carries data that could.
  !===================================================================!

  type, abstract, extends(graph) :: graph_support

   contains

     procedure(support_side_interface), deferred :: side

  end type graph_support

  abstract interface

     pure integer function graph_id_interface(this)
       import :: graph
       class(graph), intent(in) :: this
     end function graph_id_interface
     pure integer function graph_count_interface(this)
       import :: graph
       class(graph), intent(in) :: this
     end function graph_count_interface
     pure integer function graph_edge_end_interface(this, edge_index)
       import :: graph
       class(graph), intent(in) :: this
       integer, intent(in) :: edge_index
     end function graph_edge_end_interface
     pure logical function graph_edge_has_head_interface(this, edge_index)
       import :: graph
       class(graph), intent(in) :: this
       integer, intent(in) :: edge_index
     end function graph_edge_has_head_interface
     subroutine graph_member_set_interface(this, members)
       import :: graph
       class(graph), intent(in) :: this
       class(graph), allocatable, intent(out) :: members
     end subroutine graph_member_set_interface
     subroutine graph_tagged_set_interface(this, tag, members)
       import :: graph
       class(graph), intent(in) :: this
       character(len=*), intent(in) :: tag
       class(graph), allocatable, intent(out) :: members
     end subroutine graph_tagged_set_interface
     subroutine graph_part_set_interface(this, part_id, members)
       import :: graph
       class(graph), intent(in) :: this
       integer, intent(in) :: part_id
       class(graph), allocatable, intent(out) :: members
     end subroutine graph_part_set_interface
     pure subroutine graph_from_vertex_interface(this, vertex_index, indices)
       import :: graph
       class(graph), intent(in) :: this
       integer, intent(in) :: vertex_index
       integer, allocatable, intent(out) :: indices(:)
     end subroutine graph_from_vertex_interface
     pure logical function graph_has_part_relation_interface(this)
       import :: graph
       class(graph), intent(in) :: this
     end function graph_has_part_relation_interface
     pure integer function graph_global_id_interface(this, index)
       import :: graph
       class(graph), intent(in) :: this
       integer, intent(in) :: index
     end function graph_global_id_interface
     pure integer function graph_part_id_interface(this, global_index, part_id)
       import :: graph
       class(graph), intent(in) :: this
       integer, intent(in) :: global_index
       integer, intent(in) :: part_id
     end function graph_part_id_interface
     pure integer function graph_owner_part_interface(this, index)
       import :: graph
       class(graph), intent(in) :: this
       integer, intent(in) :: index
     end function graph_owner_part_interface
     pure integer function support_side_interface(this)
       import :: graph_support
       class(graph_support), intent(in) :: this
     end function support_side_interface

  end interface

end module structure_graph
