!=====================================================================!
! Abstract graph architecture interfaces.
!
! This module contains only abstract types and deferred type-bound
! procedures. Concrete graph storage, partitioning algorithms, field
! layouts, reductions, transforms, and operations must extend these
! contracts in separate class modules.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!
!
!                          THE CENTRAL LAW
!
! Everything here serves one sentence. Split a graph into parts, work
! on the parts, put the answer back together.
!
!        G'  =  P^-1 ( A ( P ( G ) ) )
!
!             G            P(G)           A(P(G))          G'
!             |             |                |              |
!             +-----P------>+-------A------->+----P^-1----->+
!             |             |                |              |
!        whole graph    the parts      worked-on parts  whole again
!
! When data rides along, it goes the same way and by the same maps:
!
!        D'  =  P^-1 ( A ( P ( G, D ) ) )
!
! P is the partitioner and nothing else. P^-1 is the assembler and
! nothing else. A is an operation and does neither.
!
!=====================================================================!
!
!                       WHAT A GRAPH IS MADE OF
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
! The second one hangs off vertex i and stops. That is a boundary
! face, and this is how it is written without inventing an imaginary
! cell on the far side of the wall.
!
!=====================================================================!
!
!                            THE TAXONOMY
!
! Six families, and each answers a different question.
!
!    graph .............. what is joined to what
!    graph_support ...... which ones we are talking about
!    graph_data ......... what values they carry
!    graph_reduction .... how many values become one
!    graph_transform .... how one graph becomes another
!    graph_operation .... how data becomes other data
!
! As a tree, with what each one adds to its parent:
!
!    graph
!
!    graph_support                     a chosen set of ids
!      graph_vertex_support            ... of vertices
!      graph_edge_support              ... of edges
!
!    graph_data                        a name and units
!      graph_field                     one value per support entry
!        graph_vertex_field            ... on vertices
!        graph_edge_field              ... on edges
!      graph_functional                one value, reduced from a field
!
!    graph_reduction                   field  -> functional
!
!    graph_transform                   graph  -> graph
!      graph_partitioner               whole  -> parts
!      graph_assembler                 parts  -> whole
!      graph_coarsener                 fine   -> coarse
!      graph_refiner                   coarse -> fine
!
!    graph_operation                   graph and data -> data
!      graph_vertex_operation
!        graph_vertex_field_operation
!        graph_vertex_functional_operation
!      graph_edge_operation
!        graph_edge_field_operation
!        graph_edge_functional_operation
!
! An operation is placed by two questions and no others. Where does it
! live, and what does it hand back?
!
!                        gives back a field    gives back one value
!                        ------------------    --------------------
!    lives on vertices    vertex_field_op      vertex_functional_op
!    lives on edges       edge_field_op        edge_functional_op
!
! That is the whole taxonomy. No word from fluid dynamics belongs in
! it. A flux, a limiter, a gradient, a Riemann solver, a boundary
! condition, a file writer, a graph colouring - each is a concrete
! class hanging off one of those four leaves.
!
! A sequence of operations is itself an operation, which is why there
! is nothing here that schedules them.
!
!=====================================================================!
!
!                     HOW FINITE VOLUME SITS ON IT
!
!    cell, control volume ......... vertex
!    interior face ................ edge with a head
!    boundary face ................ edge without a head
!    cell volume, cell centre ..... vertex field
!    face area, normal, centroid .. edge field
!    boundary group name .......... character edge field
!    the unknown q ................ vertex field
!    a flux F ..................... edge field
!    a residual ................... vertex field, from a balance
!    an objective J ............... functional
!
! A balance is cell terms plus face terms folded onto the cells they
! touch. Incidence does the folding, and it does it exactly once:
!
!                          z_e
!                   i ------------> j        y_i = y_i - z_e
!                                            y_j = y_j + z_e
!
! The word residual is what a solver calls the balance after it has
! been worked out. It does not name anything here.
!
!=====================================================================!
!
!                          THE PARTITION FRAME
!
!    owned      this part is on the hook for the answer here
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
! A part graph is still a graph. It carries the relation back to the
! whole - how many parts, which part owns what, and the id maps both
! ways - because the assembler must read that relation rather than
! invent one.
!
!=====================================================================!
!
!                        THREE LAWS, NOT PROCEDURES
!
! These replace procedures rather than accompany them.
!
! THE FIELD ORDERING LAW. A field stores its values in the order its
! support lists its ids, and stores components fastest within each
! entry. So a value sits at
!
!        (entry_position - 1) * num_components + component
!
!        support ids     7        7        3        3
!        component       1        2        1        2
!                     +--------+--------+--------+--------+
!        values       |  v(1)  |  v(2)  |  v(3)  |  v(4)  |
!                     +--------+--------+--------+--------+
!
! Because that law holds, degree-of-freedom counting, the local
! numbering of a part's vector, and its inverse are all arithmetic on
! a known layout. None of them are graph structure and none of them
! appear below.
!
! THE READ-ONLY LAW. A graph answers for data but never accepts any.
! Whatever a graph carries was fixed when it was built - geometry,
! tags, and the relation between a part and the whole. Everything
! computed travels as an operation's output and is never stored back.
! Without this law a graph slowly becomes a bag of state, which is
! what happened to the one this file replaces.
!
! THE OVERWRITE LAW. An operation writes its result, it never adds to
! it. The output argument is intent(inout) only so a caller can lend a
! buffer instead of forcing a fresh allocation every call:
!
!        call term % apply(g, data, y)      y is written, not summed
!
! A balance that sums many terms does its own summing.
!
!=====================================================================!

module abstract_graph_types

  use iso_fortran_env, only : dp => REAL64

  implicit none

  private

  public :: graph
  public :: graph_support
  public :: graph_vertex_support
  public :: graph_edge_support

  public :: graph_data
  public :: graph_field
  public :: graph_vertex_field
  public :: graph_edge_field
  public :: graph_functional
  public :: graph_reduction

  public :: graph_transform
  public :: graph_partitioner
  public :: graph_assembler
  public :: graph_coarsener
  public :: graph_refiner

  public :: graph_operation
  public :: graph_vertex_operation
  public :: graph_edge_operation
  public :: graph_vertex_field_operation
  public :: graph_vertex_functional_operation
  public :: graph_edge_field_operation
  public :: graph_edge_functional_operation

  public :: GRAPH_SUPPORT_VERTEX
  public :: GRAPH_SUPPORT_EDGE

  public :: GRAPH_FIELD_INTEGER
  public :: GRAPH_FIELD_REAL
  public :: GRAPH_FIELD_COMPLEX
  public :: GRAPH_FIELD_LOGICAL
  public :: GRAPH_FIELD_CHARACTER

  ! A support holds vertex ids or edge ids, never a mixture.
  integer, parameter :: GRAPH_SUPPORT_VERTEX = 1
  integer, parameter :: GRAPH_SUPPORT_EDGE   = 2

  ! The kind of value a field or a functional carries.
  integer, parameter :: GRAPH_FIELD_INTEGER   = 1
  integer, parameter :: GRAPH_FIELD_REAL      = 2
  integer, parameter :: GRAPH_FIELD_COMPLEX   = 3
  integer, parameter :: GRAPH_FIELD_LOGICAL   = 4
  integer, parameter :: GRAPH_FIELD_CHARACTER = 5

  !===================================================================!
  ! SUPPORTS. A chosen set of vertex ids, or a chosen set of edge ids.
  ! It owns no values and no algorithms.
  !
  !      all_vertices             boundary_edges('wall')
  !      +-------------+          +-------------+
  !      | 1 2 3 4 5 6 |          |  11  14  19 |
  !      +-------------+          +-------------+
  !
  ! These three types carry the only thing that tells a vertex from an
  ! edge anywhere else in this file. Take them away and the vertex and
  ! edge fields, and the vertex and edge operations, have nothing left
  ! to declare.
  !===================================================================!

  type, abstract :: graph_support
   contains

     !----------------------------------------------------------------!
     ! Which kind of id is in here, and how many of them.
     !----------------------------------------------------------------!

     procedure(support_kind_interface), deferred :: kind
     procedure(support_size_interface), deferred :: size

  end type graph_support

  type, abstract, extends(graph_support) :: graph_vertex_support
   contains

     !----------------------------------------------------------------!
     ! Hand over the vertex numbers.
     !----------------------------------------------------------------!

     procedure(vertex_support_ids_interface), deferred :: vertex_ids

  end type graph_vertex_support

  type, abstract, extends(graph_support) :: graph_edge_support
   contains

     !----------------------------------------------------------------!
     ! Hand over the edge numbers.
     !----------------------------------------------------------------!

     procedure(edge_support_ids_interface), deferred :: edge_ids

  end type graph_edge_support

  !===================================================================!
  ! THE GRAPH. Structure and nothing else: which vertices exist, which
  ! edges exist, where each edge goes, which sets can be named, how a
  ! part relates to the whole, and what data came with it.
  !
  ! It owns no geometry, no physics, no solver behaviour, no storage.
  !
  ! Two speeds live here, and they answer differently.
  !
  !   named sets       asked once, when an operation begins
  !                    ---> hand back a support
  !
  !   walking          asked once per vertex, inside a loop, a billion
  !                    times over a run
  !                    ---> hand back plain ids, nothing to allocate
  !===================================================================!

  type, abstract :: graph
   contains

     !----------------------------------------------------------------!
     ! Who am I and how big am I.
     !----------------------------------------------------------------!

     procedure(graph_id_interface)          , deferred :: id
     procedure(graph_num_vertices_interface), deferred :: num_vertices
     procedure(graph_num_edges_interface)   , deferred :: num_edges

     !----------------------------------------------------------------!
     ! Where an edge goes.
     !
     !            e                          b
     !     i ----------> j            i ----------o
     !
     !     tail i, head j             tail i, no head
     !     an interior face           a boundary face
     !
     ! An edge without a head lets a wall hang off one cell without
     ! inventing an imaginary cell on the far side of it.
     !----------------------------------------------------------------!

     procedure(graph_edge_tail_interface)    , deferred :: edge_tail
     procedure(graph_edge_head_interface)    , deferred :: edge_head
     procedure(graph_edge_has_head_interface), deferred :: edge_has_head

     !----------------------------------------------------------------!
     ! The named vertex sets.
     !
     !      all       o   o   o   o   o   o
     !      interior      o   o   o
     !      boundary  o           o   o   o
     !      tagged    .. whichever ones carry the name you ask for
     !
     ! These are named rather than chosen by a number, so that the
     ! compiler checks the promise. With a selector integer every
     ! implementation could read "boundary" its own way and nothing
     ! would catch it.
     !----------------------------------------------------------------!

     procedure(graph_vertex_set_interface)     , deferred :: all_vertices
     procedure(graph_vertex_set_interface)     , deferred :: interior_vertices
     procedure(graph_vertex_set_interface)     , deferred :: boundary_vertices
     procedure(graph_tagged_vertices_interface), deferred :: tagged_vertices

     !----------------------------------------------------------------!
     ! The named edge sets, the same four.
     !
     ! An interior flux runs on interior_edges. A wall condition runs
     ! on tagged_edges('wall'). Neither is a special kind of operation;
     ! they are the same kind of operation pointed at different sets.
     !----------------------------------------------------------------!

     procedure(graph_edge_set_interface)    , deferred :: all_edges
     procedure(graph_edge_set_interface)    , deferred :: interior_edges
     procedure(graph_edge_set_interface)    , deferred :: boundary_edges
     procedure(graph_tagged_edges_interface), deferred :: tagged_edges

     !----------------------------------------------------------------!
     ! The named sets of one part.
     !
     !         part 1                        part 2
     !    +---------------+            +---------------+
     !    |  o    o    o  |            |  o    o    o  |
     !    |  o    o    o--|------------|--b    o    o  |
     !    +---------------+            +---------------+
     !                 \______________/
     !                    overlap of part 1
     !
     ! The part is named rather than assumed, because a partition here
     ! is replicated: every image builds every part, so one image
     ! routinely asks about a part it does not own.
     !----------------------------------------------------------------!

     procedure(graph_part_vertices_interface), deferred :: owned_vertices
     procedure(graph_part_vertices_interface), deferred :: borrowed_vertices
     procedure(graph_part_vertices_interface), deferred :: overlap_vertices
     procedure(graph_part_edges_interface)   , deferred :: owned_edges
     procedure(graph_part_edges_interface)   , deferred :: borrowed_edges
     procedure(graph_part_edges_interface)   , deferred :: overlap_edges

     !----------------------------------------------------------------!
     ! Walking the graph, paying no attention to which way edges point.
     !
     !                      k
     !                      |
     !                 j ---i--- m
     !                      |
     !                      l
     !
     !         adjacent_vertices(i) = j, k, l, m
     !         incident_edges(i)    = the four lines touching i
     !----------------------------------------------------------------!

     procedure(graph_from_vertex_interface), deferred :: incident_edges
     procedure(graph_from_vertex_interface), deferred :: adjacent_vertices

     !----------------------------------------------------------------!
     ! Walking the graph the way its edges point.
     !
     !                   a         b
     !                    \       /
     !                     v     v
     !                       \ /
     !                        i
     !                       / \
     !                      v   v
     !                     c     d
     !
     !         incoming_vertices(i) = a, b     incoming_edges(i) = 2
     !         outgoing_vertices(i) = c, d     outgoing_edges(i) = 2
     !
     ! These four replace the separate directed graph type entirely.
     ! Direction stops being a different kind of graph and becomes
     ! something any graph can be asked about.
     !----------------------------------------------------------------!

     procedure(graph_from_vertex_interface), deferred :: outgoing_edges
     procedure(graph_from_vertex_interface), deferred :: incoming_edges
     procedure(graph_from_vertex_interface), deferred :: outgoing_vertices
     procedure(graph_from_vertex_interface), deferred :: incoming_vertices

     !----------------------------------------------------------------!
     ! How a part relates to the whole.
     !
     !      whole graph        1   2   3   4   5   6   7   8
     !                             |   |           |
     !      part 2                 1   2           3
     !
     !      full_vertex_id(2) = 4        vertex_owner_part(4) = 2
     !      part_vertex_id(4, 2) = 2
     !
     ! This lives on the graph rather than on the partitioner, because
     ! the partitioner has finished by the time anyone asks, and the
     ! assembler must read the relation rather than invent one.
     !----------------------------------------------------------------!

     procedure(graph_num_parts_interface)        , deferred :: num_parts
     procedure(graph_has_part_relation_interface), deferred :: has_part_relation
     procedure(graph_full_id_interface)          , deferred :: full_vertex_id
     procedure(graph_full_id_interface)          , deferred :: full_edge_id
     procedure(graph_part_id_interface)          , deferred :: part_vertex_id
     procedure(graph_part_id_interface)          , deferred :: part_edge_id
     procedure(graph_owner_part_interface)       , deferred :: vertex_owner_part
     procedure(graph_owner_part_interface)       , deferred :: edge_owner_part

     !----------------------------------------------------------------!
     ! The data the graph came with - geometry, tags, the part
     ! relation. Read only, by the law at the top of this file.
     !----------------------------------------------------------------!

     procedure(graph_has_data_interface), deferred :: has_data
     procedure(graph_get_data_interface), deferred :: get_data

  end type graph

  !===================================================================!
  ! DATA. A field carries a value for every entry of a support. A
  ! functional carries one value, reduced from a field.
  !
  !      field on a support              reduction        functional
  !
  !      +----+----+----+----+                              +-----+
  !      | q1 | q2 | q3 | q4 |  ------- sum, norm, ------>   |  J  |
  !      +----+----+----+----+          average, min         +-----+
  !
  ! Both say what kind of value they hold and offer a plain vector or
  ! a plain scalar of that kind. Those are adapters for talking to the
  ! outside world, not the way values are stored.
  !===================================================================!

  type, abstract :: graph_data
   contains

     !----------------------------------------------------------------!
     ! What this data is called, and what it is measured in.
     !----------------------------------------------------------------!

     procedure(data_name_interface) , deferred :: name
     procedure(data_units_interface), deferred :: units

  end type graph_data

  type, abstract, extends(graph_data) :: graph_field
   contains

     !----------------------------------------------------------------!
     ! Where the values sit, and how many of them there are.
     !----------------------------------------------------------------!

     procedure(field_support_interface)       , deferred :: support
     procedure(field_num_components_interface), deferred :: num_components
     procedure(field_num_entries_interface)   , deferred :: num_entries
     procedure(field_value_kind_interface)    , deferred :: value_kind

     !----------------------------------------------------------------!
     ! Plain vectors in and out, one pair per kind of value.
     !
     ! Character fields are real data, not a curiosity: boundary group
     ! names, material names, region labels. A numeric operation may
     ! refuse one, but the contract must not assume every field is a
     ! number.
     !----------------------------------------------------------------!

     procedure(field_get_integer_interface)  , deferred :: get_integer_vector
     procedure(field_set_integer_interface)  , deferred :: set_integer_vector
     procedure(field_get_real_interface)     , deferred :: get_real_vector
     procedure(field_set_real_interface)     , deferred :: set_real_vector
     procedure(field_get_complex_interface)  , deferred :: get_complex_vector
     procedure(field_set_complex_interface)  , deferred :: set_complex_vector
     procedure(field_get_logical_interface)  , deferred :: get_logical_vector
     procedure(field_set_logical_interface)  , deferred :: set_logical_vector
     procedure(field_get_character_interface), deferred :: get_character_vector
     procedure(field_set_character_interface), deferred :: set_character_vector

  end type graph_field

  type, abstract, extends(graph_field) :: graph_vertex_field
   contains

     !----------------------------------------------------------------!
     ! The same support, handed back as a vertex support.
     !----------------------------------------------------------------!

     procedure(vertex_field_support_interface), deferred :: vertex_support

  end type graph_vertex_field

  type, abstract, extends(graph_field) :: graph_edge_field
   contains

     !----------------------------------------------------------------!
     ! The same support, handed back as an edge support.
     !----------------------------------------------------------------!

     procedure(edge_field_support_interface), deferred :: edge_support

  end type graph_edge_field

  !===================================================================!
  ! A FUNCTIONAL. One value, of any kind a field can hold.
  !
  ! Complex is what lets a complex-step derivative live through a
  ! reduction. The derivative is the imaginary part, and a real-only
  ! functional throws it away:
  !
  !      (2.0, 1e-20)  +  (3.0, 3e-20)  =  (5.0, 4e-20)
  !                                              \
  !                                               the answer we wanted
  !
  ! Logical is what lets a question such as "is this graph acyclic"
  ! come back as true or false rather than as a one or a zero.
  !
  ! A functional also serves as the work state of a reduction while it
  ! is still running. That is why a reduction needs no numeric type of
  ! its own.
  !===================================================================!

  type, abstract, extends(graph_data) :: graph_functional
   contains

     !----------------------------------------------------------------!
     ! Which kind of value is in here.
     !----------------------------------------------------------------!

     procedure(functional_value_kind_interface), deferred :: value_kind

     !----------------------------------------------------------------!
     ! One value in and out, per kind, mirroring the field above.
     !----------------------------------------------------------------!

     procedure(functional_get_integer_interface)  , deferred :: get_integer_value
     procedure(functional_set_integer_interface)  , deferred :: set_integer_value
     procedure(functional_get_real_interface)     , deferred :: get_real_value
     procedure(functional_set_real_interface)     , deferred :: set_real_value
     procedure(functional_get_complex_interface)  , deferred :: get_complex_value
     procedure(functional_set_complex_interface)  , deferred :: set_complex_value
     procedure(functional_get_logical_interface)  , deferred :: get_logical_value
     procedure(functional_set_logical_interface)  , deferred :: set_logical_value
     procedure(functional_get_character_interface), deferred :: get_character_value
     procedure(functional_set_character_interface), deferred :: set_character_value

  end type graph_functional

  !===================================================================!
  ! A REDUCTION turns a field on a support into a functional.
  !
  ! Four steps, so that it still works when the graph is in pieces:
  !
  !      part 1 -> identity -> accumulate ---+
  !                                          +--> combine -> finalize
  !      part 2 -> identity -> accumulate ---+
  !
  ! The four steps are what keeps an average honest. Two parts with
  ! means 2 and 7 do not average to 4.5; they average to 4, because
  ! the running sum and the running count travel together and the
  ! division happens once, at the end.
  !
  ! The work state is a functional that is not finished yet, which is
  ! what keeps every procedure here free of any numeric type. A
  ! running average carries a sum and a count; an argument minimum
  ! carries a value and a location. Those live inside the concrete
  ! work state, never in this contract.
  !
  ! The four building steps are pure. The convenience that runs them
  ! all is not, deliberately, because a reduction across images sums
  ! there.
  !
  ! The three that hand back a state are intent(inout) rather than
  ! intent(out) because Fortran forbids a polymorphic intent(out)
  ! argument in a pure procedure. Each must therefore clear a live
  ! state before allocating, which also lets a caller lend one.
  !===================================================================!

  type, abstract :: graph_reduction
   contains

     !----------------------------------------------------------------!
     ! Start empty, fold values in, join two parts, finish once.
     !----------------------------------------------------------------!

     procedure(reduction_identity_interface)  , deferred :: identity
     procedure(reduction_accumulate_interface), deferred :: accumulate
     procedure(reduction_combine_interface)   , deferred :: combine
     procedure(reduction_finalize_interface)  , deferred :: finalize

     !----------------------------------------------------------------!
     ! All four in one call, for when the caller has the whole thing.
     !----------------------------------------------------------------!

     procedure(reduction_reduce_interface), deferred :: reduce

  end type graph_reduction

  !===================================================================!
  ! TRANSFORMS map between graphs and between graph data.
  !
  !      frame                          resolution
  !
  !      whole  --partitioner-->  parts     fine  --coarsener-->  coarse
  !      parts  --assembler---->  whole     coarse --refiner---->  fine
  !
  ! Partitioning and assembly change the frame without changing how
  ! much detail the graph holds. Coarsening and refinement change the
  ! detail without changing the frame. Both are transforms; they are
  ! not the same family.
  !
  ! All four stay this thin deliberately. Everything about how the
  ! pieces relate rides on the graph that comes out, which is why the
  ! graph above answers for the part relation.
  !===================================================================!

  type, abstract :: graph_transform
   contains

     !----------------------------------------------------------------!
     ! Does this transform apply to that graph, and to that data.
     !----------------------------------------------------------------!

     procedure(transform_on_graph_interface), deferred :: defined_on_graph
     procedure(transform_on_data_interface) , deferred :: defined_on_data

  end type graph_transform

  type, abstract, extends(graph_transform) :: graph_partitioner
   contains

     !----------------------------------------------------------------!
     ! P. Cut the graph up, then carry the data across the same cut.
     !----------------------------------------------------------------!

     procedure(partition_graph_interface), deferred :: partition_graph
     procedure(partition_data_interface) , deferred :: partition_data

  end type graph_partitioner

  type, abstract, extends(graph_transform) :: graph_assembler
   contains

     !----------------------------------------------------------------!
     ! P^-1, and only that. Put the graph back in whole-graph order and
     ! bring the data back with it.
     !
     ! Assembler does not mean "everything that makes a residual". No
     ! physics, no boundary conditions, no matrices, no files, no
     ! solver behaviour lives here.
     !----------------------------------------------------------------!

     procedure(assemble_graph_interface), deferred :: assemble_graph
     procedure(assemble_data_interface) , deferred :: assemble_data

  end type graph_assembler

  type, abstract, extends(graph_transform) :: graph_coarsener
   contains

     !----------------------------------------------------------------!
     ! C. Merge vertices into blocks and lift the data onto them.
     !
     !         o o o o                 O   O
     !         o o o o     ------>
     !         o o o o                 O   O
     !
     !         twelve vertices         four blocks, same shape,
     !                                 less detail
     !----------------------------------------------------------------!

     procedure(coarsen_graph_interface), deferred :: coarsen_graph
     procedure(coarsen_data_interface) , deferred :: coarsen_data

  end type graph_coarsener

  type, abstract, extends(graph_transform) :: graph_refiner
   contains

     !----------------------------------------------------------------!
     ! R. Split vertices and carry the data down. The other way round.
     !----------------------------------------------------------------!

     procedure(refine_graph_interface), deferred :: refine_graph
     procedure(refine_data_interface) , deferred :: refine_data

  end type graph_refiner

  !===================================================================!
  ! OPERATIONS. The whole taxonomy is two questions: does it live on
  ! vertices or on edges, and does it hand back a field or one value.
  !
  !                        gives back a field    gives back one value
  !                        ------------------    --------------------
  !    lives on vertices    vertex_field_op      vertex_functional_op
  !    lives on edges       edge_field_op        edge_functional_op
  !
  ! Interior, boundary, tagged, owned, borrowed and overlap are not
  ! kinds of operation. They are the sets an operation is pointed at.
  ! The type says what is being worked out; the support says where.
  !
  ! A sequence of operations is itself an operation, which is why
  ! there is nothing here that schedules them.
  !===================================================================!

  type, abstract :: graph_operation
   contains

     !----------------------------------------------------------------!
     ! What this operation is called.
     !----------------------------------------------------------------!

     procedure(operation_name_interface), deferred :: name

  end type graph_operation

  type, abstract, extends(graph_operation) :: graph_vertex_operation
   contains

     !----------------------------------------------------------------!
     ! Which vertices it works on. Asked once, not per vertex.
     !----------------------------------------------------------------!

     procedure(vertex_operation_support_interface), deferred :: support

  end type graph_vertex_operation

  type, abstract, extends(graph_operation) :: graph_edge_operation
   contains

     !----------------------------------------------------------------!
     ! Which edges it works on. Asked once, not per edge.
     !----------------------------------------------------------------!

     procedure(edge_operation_support_interface), deferred :: support

  end type graph_edge_operation

  type, abstract, extends(graph_vertex_operation) :: graph_vertex_field_operation
   contains

     !----------------------------------------------------------------!
     ! Cell terms: a source, a reaction, a mass term, a gradient, a
     ! limiter, a preconditioner, a Newton step, a time step.
     !----------------------------------------------------------------!

     procedure(vertex_field_apply_interface), deferred :: apply

  end type graph_vertex_field_operation

  type, abstract, extends(graph_vertex_operation) :: graph_vertex_functional_operation
   contains

     !----------------------------------------------------------------!
     ! Cells reduced to one number: an objective, an energy, a norm.
     !----------------------------------------------------------------!

     procedure(vertex_functional_apply_interface), deferred :: apply

  end type graph_vertex_functional_operation

  type, abstract, extends(graph_edge_operation) :: graph_edge_field_operation
   contains

     !----------------------------------------------------------------!
     ! Face terms: a diffusion flux, an advection flux, a Riemann
     ! solver, a wall condition, one edge's share of a matrix.
     !----------------------------------------------------------------!

     procedure(edge_field_apply_interface), deferred :: apply

  end type graph_edge_field_operation

  type, abstract, extends(graph_edge_operation) :: graph_edge_functional_operation
   contains

     !----------------------------------------------------------------!
     ! Faces reduced to one number: a face objective, a jump norm.
     !----------------------------------------------------------------!

     procedure(edge_functional_apply_interface), deferred :: apply

  end type graph_edge_functional_operation

  !===================================================================!
  ! THE DEFERRED PROCEDURES
  !===================================================================!

  abstract interface

     !===============================================================!
     ! Supports.
     !===============================================================!

     pure integer function support_kind_interface(this)
       import :: graph_support
       class(graph_support), intent(in) :: this
     end function support_kind_interface

     pure integer function support_size_interface(this)
       import :: graph_support
       class(graph_support), intent(in) :: this
     end function support_size_interface

     pure subroutine vertex_support_ids_interface(this, ids)
       import :: graph_vertex_support
       class(graph_vertex_support), intent(in) :: this
       integer, allocatable, intent(out) :: ids(:)
     end subroutine vertex_support_ids_interface

     pure subroutine edge_support_ids_interface(this, ids)
       import :: graph_edge_support
       class(graph_edge_support), intent(in) :: this
       integer, allocatable, intent(out) :: ids(:)
     end subroutine edge_support_ids_interface

     !===============================================================!
     ! Graph identity, size, and edge endpoints.
     !===============================================================!

     pure integer function graph_id_interface(this)
       import :: graph
       class(graph), intent(in) :: this
     end function graph_id_interface

     pure integer function graph_num_vertices_interface(this)
       import :: graph
       class(graph), intent(in) :: this
     end function graph_num_vertices_interface

     pure integer function graph_num_edges_interface(this)
       import :: graph
       class(graph), intent(in) :: this
     end function graph_num_edges_interface

     pure integer function graph_edge_tail_interface(this, edge_id)
       import :: graph
       class(graph), intent(in) :: this
       integer, intent(in) :: edge_id
     end function graph_edge_tail_interface

     pure integer function graph_edge_head_interface(this, edge_id)
       import :: graph
       class(graph), intent(in) :: this
       integer, intent(in) :: edge_id
     end function graph_edge_head_interface

     pure logical function graph_edge_has_head_interface(this, edge_id)
       import :: graph
       class(graph), intent(in) :: this
       integer, intent(in) :: edge_id
     end function graph_edge_has_head_interface

     !===============================================================!
     ! The named sets. Asked once, when an operation begins, so the
     ! answer is a support and the cost is spread over a whole sweep.
     !===============================================================!

     subroutine graph_vertex_set_interface(this, support)
       import :: graph, graph_vertex_support
       class(graph), intent(in) :: this
       class(graph_vertex_support), allocatable, intent(out) :: support
     end subroutine graph_vertex_set_interface

     subroutine graph_edge_set_interface(this, support)
       import :: graph, graph_edge_support
       class(graph), intent(in) :: this
       class(graph_edge_support), allocatable, intent(out) :: support
     end subroutine graph_edge_set_interface

     subroutine graph_tagged_vertices_interface(this, tag, support)
       import :: graph, graph_vertex_support
       class(graph), intent(in) :: this
       character(len=*), intent(in) :: tag
       class(graph_vertex_support), allocatable, intent(out) :: support
     end subroutine graph_tagged_vertices_interface

     subroutine graph_tagged_edges_interface(this, tag, support)
       import :: graph, graph_edge_support
       class(graph), intent(in) :: this
       character(len=*), intent(in) :: tag
       class(graph_edge_support), allocatable, intent(out) :: support
     end subroutine graph_tagged_edges_interface

     subroutine graph_part_vertices_interface(this, part_id, support)
       import :: graph, graph_vertex_support
       class(graph), intent(in) :: this
       integer, intent(in) :: part_id
       class(graph_vertex_support), allocatable, intent(out) :: support
     end subroutine graph_part_vertices_interface

     subroutine graph_part_edges_interface(this, part_id, support)
       import :: graph, graph_edge_support
       class(graph), intent(in) :: this
       integer, intent(in) :: part_id
       class(graph_edge_support), allocatable, intent(out) :: support
     end subroutine graph_part_edges_interface

     !===============================================================!
     ! Walking the graph. Asked once per vertex, inside loops, so the
     ! answer is bare ids and the procedure is pure. Handing back a
     ! support here would allocate three times per neighbour query.
     !===============================================================!

     pure subroutine graph_from_vertex_interface(this, vertex_id, ids)
       import :: graph
       class(graph), intent(in) :: this
       integer, intent(in) :: vertex_id
       integer, allocatable, intent(out) :: ids(:)
     end subroutine graph_from_vertex_interface

     !===============================================================!
     ! How a part relates to the whole.
     !===============================================================!

     !---------------------------------------------------------------!
     ! How many pieces this graph was cut into. One, for a graph that
     ! was never cut.
     !---------------------------------------------------------------!

     pure integer function graph_num_parts_interface(this)
       import :: graph
       class(graph), intent(in) :: this
     end function graph_num_parts_interface

     !---------------------------------------------------------------!
     ! Does this graph remember where it came from?
     !
     !      straight off a mesh file  ->  no
     !      out of a partitioner      ->  yes
     !
     ! An assembler asks this first. Without the relation there is no
     ! way back to whole-graph order, and it should say so rather than
     ! guess.
     !---------------------------------------------------------------!

     pure logical function graph_has_part_relation_interface(this)
       import :: graph
       class(graph), intent(in) :: this
     end function graph_has_part_relation_interface

     !---------------------------------------------------------------!
     ! Part numbering out to whole-graph numbering.
     !
     !      whole    1   2   3   4   5   6   7   8
     !                       |   |           |
     !      part 2           1   2           3
     !
     !      full_vertex_id(2) = 4
     !
     ! An assembler walks a part's own numbering and uses this to find
     ! where each value belongs in the whole. A file writer uses it to
     ! print cell numbers a person will recognise.
     !---------------------------------------------------------------!

     pure integer function graph_full_id_interface(this, part_local_id)
       import :: graph
       class(graph), intent(in) :: this
       integer, intent(in) :: part_local_id
     end function graph_full_id_interface

     !---------------------------------------------------------------!
     ! The same map read backwards: where does whole-graph vertex 4
     ! sit inside part 2?
     !
     !      part_vertex_id(4, 2) = 2
     !
     ! The part is named because the partition here is replicated -
     ! every image builds every part, so one image regularly asks
     ! about a part it does not own.
     !---------------------------------------------------------------!

     pure integer function graph_part_id_interface(this, full_id, part_id)
       import :: graph
       class(graph), intent(in) :: this
       integer, intent(in) :: full_id
       integer, intent(in) :: part_id
     end function graph_part_id_interface

     !---------------------------------------------------------------!
     ! Whose job is this one?
     !
     !      vertex_owner_part(4) = 2       part 2 answers for it
     !                                     everyone else only reads it
     !
     ! This is what stops a shared cell being counted twice when the
     ! parts are added back together.
     !---------------------------------------------------------------!

     pure integer function graph_owner_part_interface(this, id)
       import :: graph
       class(graph), intent(in) :: this
       integer, intent(in) :: id
     end function graph_owner_part_interface

     !===============================================================!
     ! The data the graph came with.
     !===============================================================!

     !---------------------------------------------------------------!
     ! Does this graph carry something by that name? An operation asks
     ! before it reaches, so a mesh without face normals says no
     ! rather than crashing.
     !---------------------------------------------------------------!

     pure logical function graph_has_data_interface(this, name)
       import :: graph
       class(graph), intent(in) :: this
       character(len=*), intent(in) :: name
     end function graph_has_data_interface

     !---------------------------------------------------------------!
     ! Everything the mesh knew when it was built, fetched by name.
     !
     !      get_data('cell_volume')     -> a real vertex field
     !      get_data('cell_centre')     -> a real vertex field, 3 wide
     !      get_data('face_area')       -> a real edge field
     !      get_data('face_normal')     -> a real edge field, 3 wide
     !      get_data('boundary_name')   -> a character edge field
     !
     ! This is how geometry reaches an operation without being passed
     ! down through every call. A flux asks the graph for the normal
     ! and the distance; nobody has to hand them over.
     !
     ! Read only. Whatever gets computed leaves as an operation's
     ! output and never comes back here to be stored.
     !---------------------------------------------------------------!

     subroutine graph_get_data_interface(this, name, data)
       import :: graph, graph_data
       class(graph), intent(in) :: this
       character(len=*), intent(in) :: name
       class(graph_data), allocatable, intent(out) :: data
     end subroutine graph_get_data_interface

     !===============================================================!
     ! Data identity.
     !===============================================================!

     pure function data_name_interface(this) result(name)
       import :: graph_data
       class(graph_data), intent(in) :: this
       character(len=:), allocatable :: name
     end function data_name_interface

     pure function data_units_interface(this) result(units)
       import :: graph_data
       class(graph_data), intent(in) :: this
       character(len=:), allocatable :: units
     end function data_units_interface

     !===============================================================!
     ! Fields.
     !===============================================================!

     subroutine field_support_interface(this, support)
       import :: graph_field, graph_support
       class(graph_field), intent(in) :: this
       class(graph_support), allocatable, intent(out) :: support
     end subroutine field_support_interface

     subroutine vertex_field_support_interface(this, support)
       import :: graph_vertex_field, graph_vertex_support
       class(graph_vertex_field), intent(in) :: this
       class(graph_vertex_support), allocatable, intent(out) :: support
     end subroutine vertex_field_support_interface

     subroutine edge_field_support_interface(this, support)
       import :: graph_edge_field, graph_edge_support
       class(graph_edge_field), intent(in) :: this
       class(graph_edge_support), allocatable, intent(out) :: support
     end subroutine edge_field_support_interface

     pure integer function field_num_components_interface(this)
       import :: graph_field
       class(graph_field), intent(in) :: this
     end function field_num_components_interface

     pure integer function field_num_entries_interface(this)
       import :: graph_field
       class(graph_field), intent(in) :: this
     end function field_num_entries_interface

     pure integer function field_value_kind_interface(this)
       import :: graph_field
       class(graph_field), intent(in) :: this
     end function field_value_kind_interface

     !---------------------------------------------------------------!
     ! THE PLAIN-VECTOR ADAPTERS. A way in and out for code that does
     ! not speak graph.
     !
     !      field  ---get--->  [ v1 v2 v3 v4 ]  ---> a Krylov solver,
     !                                                a file writer,
     !                                                an outside library
     !      field  <--set----  [ v1 v2 v3 v4 ]  <--- the answer coming
     !                                                back
     !
     ! Fetch once, work in plain arrays where the arithmetic is free,
     ! write back once. Scaling and adding are not graph theory and
     ! have no procedures here on purpose.
     !
     ! One pair per kind of value. Integer for a colouring, a part
     ! number, a visit order. Real for the ordinary state. Complex for
     ! a complex-step derivative. Logical for a mask. Character for
     ! boundary group names and material names - those are real data
     ! too, and a numeric operation is free to refuse them.
     !---------------------------------------------------------------!

     pure subroutine field_get_integer_interface(this, values)
       import :: graph_field
       class(graph_field), intent(in) :: this
       integer, allocatable, intent(out) :: values(:)
     end subroutine field_get_integer_interface

     pure subroutine field_set_integer_interface(this, values)
       import :: graph_field
       class(graph_field), intent(inout) :: this
       integer, intent(in) :: values(:)
     end subroutine field_set_integer_interface

     pure subroutine field_get_real_interface(this, values)
       import :: graph_field, dp
       class(graph_field), intent(in) :: this
       real(dp), allocatable, intent(out) :: values(:)
     end subroutine field_get_real_interface

     pure subroutine field_set_real_interface(this, values)
       import :: graph_field, dp
       class(graph_field), intent(inout) :: this
       real(dp), intent(in) :: values(:)
     end subroutine field_set_real_interface

     pure subroutine field_get_complex_interface(this, values)
       import :: graph_field, dp
       class(graph_field), intent(in) :: this
       complex(dp), allocatable, intent(out) :: values(:)
     end subroutine field_get_complex_interface

     pure subroutine field_set_complex_interface(this, values)
       import :: graph_field, dp
       class(graph_field), intent(inout) :: this
       complex(dp), intent(in) :: values(:)
     end subroutine field_set_complex_interface

     pure subroutine field_get_logical_interface(this, values)
       import :: graph_field
       class(graph_field), intent(in) :: this
       logical, allocatable, intent(out) :: values(:)
     end subroutine field_get_logical_interface

     pure subroutine field_set_logical_interface(this, values)
       import :: graph_field
       class(graph_field), intent(inout) :: this
       logical, intent(in) :: values(:)
     end subroutine field_set_logical_interface

     pure subroutine field_get_character_interface(this, values)
       import :: graph_field
       class(graph_field), intent(in) :: this
       character(len=:), allocatable, intent(out) :: values(:)
     end subroutine field_get_character_interface

     pure subroutine field_set_character_interface(this, values)
       import :: graph_field
       class(graph_field), intent(inout) :: this
       character(len=*), intent(in) :: values(:)
     end subroutine field_set_character_interface

     !===============================================================!
     ! Functionals.
     !===============================================================!

     pure integer function functional_value_kind_interface(this)
       import :: graph_functional
       class(graph_functional), intent(in) :: this
     end function functional_value_kind_interface

     pure subroutine functional_get_integer_interface(this, value)
       import :: graph_functional
       class(graph_functional), intent(in) :: this
       integer, intent(out) :: value
     end subroutine functional_get_integer_interface

     pure subroutine functional_set_integer_interface(this, value)
       import :: graph_functional
       class(graph_functional), intent(inout) :: this
       integer, intent(in) :: value
     end subroutine functional_set_integer_interface

     pure subroutine functional_get_real_interface(this, value)
       import :: graph_functional, dp
       class(graph_functional), intent(in) :: this
       real(dp), intent(out) :: value
     end subroutine functional_get_real_interface

     pure subroutine functional_set_real_interface(this, value)
       import :: graph_functional, dp
       class(graph_functional), intent(inout) :: this
       real(dp), intent(in) :: value
     end subroutine functional_set_real_interface

     pure subroutine functional_get_complex_interface(this, value)
       import :: graph_functional, dp
       class(graph_functional), intent(in) :: this
       complex(dp), intent(out) :: value
     end subroutine functional_get_complex_interface

     pure subroutine functional_set_complex_interface(this, value)
       import :: graph_functional, dp
       class(graph_functional), intent(inout) :: this
       complex(dp), intent(in) :: value
     end subroutine functional_set_complex_interface

     pure subroutine functional_get_logical_interface(this, value)
       import :: graph_functional
       class(graph_functional), intent(in) :: this
       logical, intent(out) :: value
     end subroutine functional_get_logical_interface

     pure subroutine functional_set_logical_interface(this, value)
       import :: graph_functional
       class(graph_functional), intent(inout) :: this
       logical, intent(in) :: value
     end subroutine functional_set_logical_interface

     pure subroutine functional_get_character_interface(this, value)
       import :: graph_functional
       class(graph_functional), intent(in) :: this
       character(len=:), allocatable, intent(out) :: value
     end subroutine functional_get_character_interface

     pure subroutine functional_set_character_interface(this, value)
       import :: graph_functional
       class(graph_functional), intent(inout) :: this
       character(len=*), intent(in) :: value
     end subroutine functional_set_character_interface

     !===============================================================!
     ! Reductions. The state is an unfinished functional, which is what
     ! keeps every numeric type out of this contract.
     !===============================================================!

     !---------------------------------------------------------------!
     ! The whole four-step dance, and why it is four steps:
     !
     !   part 1  [q q q]  --accumulate-->  (sum 6, count 3) --+
     !                                                        +--> combine
     !   part 2  [q q]    --accumulate-->  (sum 14, count 2) -+        |
     !                                                                 v
     !                                                    (sum 20, count 5)
     !                                                                 |
     !                                                             finalize
     !                                                                 |
     !                                                                 v
     !                                                            J = 4.0
     !
     ! Means of 2 and 7 do not average to 4.5. They average to 4,
     ! because the sum and the count travel together and the division
     ! happens once, at the very end. A reduction that finished early
     ! on each part would get this wrong.
     !
     ! Start empty. For a sum that is zero, for a product one, for a
     ! minimum the largest number there is.
     !---------------------------------------------------------------!

     pure subroutine reduction_identity_interface(this, state)
       import :: graph_reduction, graph_functional
       class(graph_reduction), intent(in) :: this
       class(graph_functional), allocatable, intent(inout) :: state
     end subroutine reduction_identity_interface

     !---------------------------------------------------------------!
     ! Fold one part's values into the running state.
     !
     ! The optional measure is what turns a bare sum into an integral:
     ! weight each cell by its volume, or each face by its area, and
     ! the answer stops depending on how finely the mesh was cut.
     !
     !      sum       J = sum q_i
     !      integral  J = sum q_i V_i          <- measure is the volume
     !      average   J = sum q_i V_i / sum V_i
     !      norm      J = sqrt( sum q_i^2 V_i )
     !---------------------------------------------------------------!

     pure subroutine reduction_accumulate_interface(this, input_graph, field, support, state, measure)
       import :: graph_reduction, graph, graph_field
       import :: graph_support, graph_functional
       class(graph_reduction), intent(in) :: this
       class(graph), intent(in) :: input_graph
       class(graph_field), intent(in) :: field
       class(graph_support), intent(in) :: support
       class(graph_functional), intent(inout) :: state
       class(graph_field), intent(in), optional :: measure
     end subroutine reduction_accumulate_interface

     !---------------------------------------------------------------!
     ! Join two part answers into one. This is what lets the same
     ! reduction serve one image or a thousand, and it must not care
     ! which order the parts arrive in.
     !---------------------------------------------------------------!

     pure subroutine reduction_combine_interface(this, left, right, combined)
       import :: graph_reduction, graph_functional
       class(graph_reduction), intent(in) :: this
       class(graph_functional), intent(in) :: left
       class(graph_functional), intent(in) :: right
       class(graph_functional), allocatable, intent(inout) :: combined
     end subroutine reduction_combine_interface

     !---------------------------------------------------------------!
     ! Finish, once, after every part has been folded in.
     !
     !      a sum        nothing left to do
     !      an average   now divide the sum by the count
     !      a norm       now take the square root
     !
     ! Doing either of those any earlier is how a parallel run quietly
     ! gets a different answer from a serial one.
     !---------------------------------------------------------------!

     pure subroutine reduction_finalize_interface(this, state, functional)
       import :: graph_reduction, graph_functional
       class(graph_reduction), intent(in) :: this
       class(graph_functional), intent(in) :: state
       class(graph_functional), allocatable, intent(inout) :: functional
     end subroutine reduction_finalize_interface

     !---------------------------------------------------------------!
     ! All four in one call, for a caller holding the whole thing.
     !
     ! Not pure, and that is deliberate: a reduction spread across
     ! images has to talk to the other images somewhere, and this is
     ! the one place in the file where that is allowed.
     !---------------------------------------------------------------!
     subroutine reduction_reduce_interface(this, input_graph, field, support, functional, measure)
       import :: graph_reduction, graph, graph_field
       import :: graph_support, graph_functional
       class(graph_reduction), intent(in) :: this
       class(graph), intent(in) :: input_graph
       class(graph_field), intent(in) :: field
       class(graph_support), intent(in) :: support
       class(graph_functional), allocatable, intent(inout) :: functional
       class(graph_field), intent(in), optional :: measure
     end subroutine reduction_reduce_interface

     !===============================================================!
     ! Transforms.
     !===============================================================!

     !---------------------------------------------------------------!
     ! A gate, asked before anything is attempted.
     !
     !      an assembler refuses a graph with no part relation on it
     !      a coarsener refuses a graph already down to one vertex
     !      a geometric partitioner refuses a graph with no coordinates
     !
     ! Better to answer no here than to fail halfway through and leave
     ! a half-built graph behind.
     !---------------------------------------------------------------!

     pure logical function transform_on_graph_interface(this, input_graph)
       import :: graph_transform, graph
       class(graph_transform), intent(in) :: this
       class(graph), intent(in) :: input_graph
     end function transform_on_graph_interface

     !---------------------------------------------------------------!
     ! The same gate for the data. A transform that can move a cell
     ! field may have nothing sensible to say about a face field, or
     ! about names rather than numbers.
     !---------------------------------------------------------------!

     pure logical function transform_on_data_interface(this, input_graph, input_data)
       import :: graph_transform, graph, graph_data
       class(graph_transform), intent(in) :: this
       class(graph), intent(in) :: input_graph
       class(graph_data), intent(in) :: input_data
     end function transform_on_data_interface

     !---------------------------------------------------------------!
     ! P. Cut the graph up.
     !
     !         o---o---o---o---o---o
     !                     :                 cut where few edges cross,
     !         o---o---o   :   o---o---o     so the parts talk little
     !                part 1     part 2
     !
     ! Different ways to choose the cut, all the same contract:
     !
     !      breadth first    grow each part outward from a seed cell
     !      geometric        cut on cell coordinates, recursively
     !      aggregate        glue small clusters together
     !      adopted          take a map somebody else already computed
     !
     ! What comes out is still a graph. It simply also knows how it
     ! relates back to the whole.
     !---------------------------------------------------------------!

     subroutine partition_graph_interface(this, full_graph, part_graph)
       import :: graph_partitioner, graph
       class(graph_partitioner), intent(in) :: this
       class(graph), intent(in) :: full_graph
       class(graph), allocatable, intent(out) :: part_graph
     end subroutine partition_graph_interface

     !---------------------------------------------------------------!
     ! Carry the data across the very same cut. Geometry, tags and the
     ! current state all travel by the map the graph already holds, so
     ! the values cannot drift out of step with the structure.
     !---------------------------------------------------------------!

     subroutine partition_data_interface(this, full_graph, full_data, part_graph, part_data)
       import :: graph_partitioner, graph, graph_data
       class(graph_partitioner), intent(in) :: this
       class(graph), intent(in) :: full_graph
       class(graph_data), intent(in) :: full_data
       class(graph), intent(in) :: part_graph
       class(graph_data), allocatable, intent(out) :: part_data
     end subroutine partition_data_interface

     !---------------------------------------------------------------!
     ! P^-1, and only that. Put the graph back in whole-graph order.
     !
     ! The law this has to satisfy:
     !
     !      assemble( partition( G ) )  ==  G
     !
     ! Assembler does not mean "everything that makes a residual". No
     ! physics, no boundary conditions, no matrices, no files, no
     ! solver behaviour belongs in here.
     !---------------------------------------------------------------!

     subroutine assemble_graph_interface(this, part_graph, full_graph)
       import :: graph_assembler, graph
       class(graph_assembler), intent(in) :: this
       class(graph), intent(in) :: part_graph
       class(graph), allocatable, intent(out) :: full_graph
     end subroutine assemble_graph_interface

     !---------------------------------------------------------------!
     ! Bring the data back with it, and the same law applies:
     !
     !      assemble( partition( G, D ) )  ==  ( G, D )
     !
     ! Only owned values are collected. A borrowed value is somebody
     ! else's copy, and counting it too is how a conserved quantity
     ! quietly stops being conserved.
     !---------------------------------------------------------------!

     subroutine assemble_data_interface(this, part_graph, part_data, full_graph, full_data)
       import :: graph_assembler, graph, graph_data
       class(graph_assembler), intent(in) :: this
       class(graph), intent(in) :: part_graph
       class(graph_data), intent(in) :: part_data
       class(graph), intent(in) :: full_graph
       class(graph_data), allocatable, intent(out) :: full_data
     end subroutine assemble_data_interface

     !---------------------------------------------------------------!
     ! C. Fewer, larger vertices.
     !
     !      o o o o                O   O
     !      o o o o    ------>              a multigrid level, where
     !      o o o o                O   O    the slow, smooth part of
     !                                      the error is cheap to kill
     !
     ! Also how a coarse graph gets built for a first guess, or for a
     ! quick look at a mesh too big to draw.
     !---------------------------------------------------------------!

     subroutine coarsen_graph_interface(this, fine_graph, coarse_graph)
       import :: graph_coarsener, graph
       class(graph_coarsener), intent(in) :: this
       class(graph), intent(in) :: fine_graph
       class(graph), allocatable, intent(out) :: coarse_graph
     end subroutine coarsen_graph_interface

     !---------------------------------------------------------------!
     ! Lift the data up to the coarse graph. Several fine values land
     ! on one coarse vertex, so this has to say how they merge - added
     ! for a residual, averaged for a state, volume-weighted when the
     ! cells are unequal.
     !---------------------------------------------------------------!

     subroutine coarsen_data_interface(this, fine_graph, fine_data, coarse_graph, coarse_data)
       import :: graph_coarsener, graph, graph_data
       class(graph_coarsener), intent(in) :: this
       class(graph), intent(in) :: fine_graph
       class(graph_data), intent(in) :: fine_data
       class(graph), intent(in) :: coarse_graph
       class(graph_data), allocatable, intent(out) :: coarse_data
     end subroutine coarsen_data_interface

     !---------------------------------------------------------------!
     ! R. The other way. One vertex becomes several.
     !
     !         O   O                  o o o o
     !                   ------>      o o o o     more detail where an
     !         O   O                  o o o o     error measure says it
     !                                            is needed
     !
     ! Used to sharpen a mesh around a shock or a boundary layer, and
     ! to carry a coarse multigrid correction back down.
     !---------------------------------------------------------------!

     subroutine refine_graph_interface(this, coarse_graph, fine_graph)
       import :: graph_refiner, graph
       class(graph_refiner), intent(in) :: this
       class(graph), intent(in) :: coarse_graph
       class(graph), allocatable, intent(out) :: fine_graph
     end subroutine refine_graph_interface

     !---------------------------------------------------------------!
     ! Carry the data down. One coarse value feeds several fine ones,
     ! so this says how it spreads - copied straight down, or
     ! interpolated so the result stays smooth across the new cells.
     !---------------------------------------------------------------!

     subroutine refine_data_interface(this, coarse_graph, coarse_data, fine_graph, fine_data)
       import :: graph_refiner, graph, graph_data
       class(graph_refiner), intent(in) :: this
       class(graph), intent(in) :: coarse_graph
       class(graph_data), intent(in) :: coarse_data
       class(graph), intent(in) :: fine_graph
       class(graph_data), allocatable, intent(out) :: fine_data
     end subroutine refine_data_interface

     !===============================================================!
     ! Operations. The output is intent(inout) so a caller may lend a
     ! buffer; by the overwrite law the operation still writes its
     ! result rather than adding to it.
     !===============================================================!

     !---------------------------------------------------------------!
     ! What to call this operation in a log line, a configuration
     ! file, or the message when it refuses to run.
     !
     !      'diffusion flux'   'wall condition'   'cell source'
     !---------------------------------------------------------------!

     pure function operation_name_interface(this) result(name)
       import :: graph_operation
       class(graph_operation), intent(in) :: this
       character(len=:), allocatable :: name
     end function operation_name_interface

     !---------------------------------------------------------------!
     ! Which cells am I responsible for? Asked once; everything after
     ! it is a loop. The type says what is worked out, this says where.
     !
     !      a source everywhere   -> all_vertices
     !          o  o  o  o  o  o      every cell
     !
     !      a heater in a patch   -> tagged_vertices('heater')
     !          .  .  o  o  .  .      only the ones carrying the name
     !
     !      my share of the work  -> owned_vertices(me)
     !          o  o  o  .  .  .      the rest belong to someone else
     !---------------------------------------------------------------!

     subroutine vertex_operation_support_interface(this, input_graph, support)
       import :: graph_vertex_operation, graph, graph_vertex_support
       class(graph_vertex_operation), intent(in) :: this
       class(graph), intent(in) :: input_graph
       class(graph_vertex_support), allocatable, intent(out) :: support
     end subroutine vertex_operation_support_interface

     !---------------------------------------------------------------!
     ! The same question for faces.
     !
     !      an interior flux   -> interior_edges
     !      a wall condition   -> tagged_edges('wall')
     !      an inlet           -> tagged_edges('inlet')
     !
     ! Notice that a boundary condition is not a special kind of
     ! operation here. It is the same edge operation aimed at a
     ! different set of edges.
     !---------------------------------------------------------------!

     subroutine edge_operation_support_interface(this, input_graph, support)
       import :: graph_edge_operation, graph, graph_edge_support
       class(graph_edge_operation), intent(in) :: this
       class(graph), intent(in) :: input_graph
       class(graph_edge_support), allocatable, intent(out) :: support
     end subroutine edge_operation_support_interface

     !---------------------------------------------------------------!
     ! CELL TERMS. A value per cell goes in, a value per cell comes
     ! out.
     !
     !         what each cell holds        what this term contributes
     !
     !            o   o   o                     +   +   +
     !            o   o   o     --apply-->      +   +   +
     !            o   o   o                     +   +   +
     !
     ! Who asks for this:
     !
     !      a source or reaction     S(q) added to every cell
     !      a mass or storage term   how much a cell can hold
     !      a gradient               grad q at the cell centre, so a
     !                               face can reconstruct to second
     !                               order
     !      a limiter                clip that gradient near a shock
     !      a preconditioner         y = M^-1 x, one sweep
     !      a Newton step            u_new from u_old and the residual
     !      a time step              q at the next instant
     !      a balance                every cell term plus every face
     !                               term folded in - what a solver
     !                               calls the residual
     !---------------------------------------------------------------!

     subroutine vertex_field_apply_interface(this, input_graph, input_data, output)
       import :: graph_vertex_field_operation, graph
       import :: graph_data, graph_vertex_field
       class(graph_vertex_field_operation), intent(in) :: this
       class(graph), intent(in) :: input_graph
       class(graph_data), intent(in), optional :: input_data(:)
       class(graph_vertex_field), allocatable, intent(inout) :: output
     end subroutine vertex_field_apply_interface

     !---------------------------------------------------------------!
     ! CELLS IN, ONE NUMBER OUT.
     !
     !            o  o  o  o
     !            o  o  o  o    --apply-->    J
     !            o  o  o  o
     !
     ! Who asks for this:
     !
     !      total energy in the field
     !      an objective a design optimiser is trying to move
     !      the residual norm a solver stops on
     !      an error measure that decides where to refine
     !---------------------------------------------------------------!

     subroutine vertex_functional_apply_interface(this, input_graph, input_data, output)
       import :: graph_vertex_functional_operation, graph
       import :: graph_data, graph_functional
       class(graph_vertex_functional_operation), intent(in) :: this
       class(graph), intent(in) :: input_graph
       class(graph_data), intent(in), optional :: input_data(:)
       class(graph_functional), allocatable, intent(inout) :: output
     end subroutine vertex_functional_apply_interface

     !---------------------------------------------------------------!
     ! FACE TERMS. This is where the physics between two cells lives.
     !
     !                          F_e
     !            (i) --------------------> (j)
     !            q_i                        q_j
     !
     !      diffusion    F = -k (q_j - q_i) / d
     !      advection    F = u q, taken from whichever side is upwind
     !      a Riemann solver decides F when the two sides disagree
     !
     ! At a wall there is no cell on the far side:
     !
     !            (i) ----------------o
     !            q_i                  the tag says what is out there
     !
     ! Either way the balance then folds F onto the cells it touches,
     ! once and only once:
     !
     !            y_i = y_i - F_e        y_j = y_j + F_e
     !
     ! One more customer: a matrix. A sparse matrix is a graph with a
     ! number on every edge, so filling those numbers in is an edge
     ! field operation like any other.
     !---------------------------------------------------------------!

     subroutine edge_field_apply_interface(this, input_graph, input_data, output)
       import :: graph_edge_field_operation, graph
       import :: graph_data, graph_edge_field
       class(graph_edge_field_operation), intent(in) :: this
       class(graph), intent(in) :: input_graph
       class(graph_data), intent(in), optional :: input_data(:)
       class(graph_edge_field), allocatable, intent(inout) :: output
     end subroutine edge_field_apply_interface

     !---------------------------------------------------------------!
     ! FACES IN, ONE NUMBER OUT.
     !
     !      tagged_edges('wall')
     !           |  |  |  |  |      --apply-->    J
     !
     ! Who asks for this:
     !
     !      drag and lift, by adding up what every wall face pushes
     !      the heat passing through a named surface
     !      how much mass an inlet or outlet lets through
     !      the size of the jumps between cells, as a smoothness check
     !---------------------------------------------------------------!

     subroutine edge_functional_apply_interface(this, input_graph, input_data, output)
       import :: graph_edge_functional_operation, graph
       import :: graph_data, graph_functional
       class(graph_edge_functional_operation), intent(in) :: this
       class(graph), intent(in) :: input_graph
       class(graph_data), intent(in), optional :: input_data(:)
       class(graph_functional), allocatable, intent(inout) :: output
     end subroutine edge_functional_apply_interface

  end interface

end module abstract_graph_types
