!=====================================================================!
! LEVEL 0 OF THE STRATIFICATION . THE GRAMMAR
!
! The ground level of the tower. This module answers one question and
! no other: WHAT MAY EXIST. It holds four abstract roles and their
! operation symbols - no mathematics, no algorithms, no physics.
! Those arrive on the levels above, each one a partial concretion of
! what stands here.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!
!
!                         THE STRATIFICATION
!
! The codebase is a tower of levels. Walking down an edge of the
! tower binds a parameter of an abstraction (concretion); walking up
! unbinds it (abstraction); the retired name survives as an argument
! value (absorption). One module per level; the extends edge carries
! the level; no level number ever enters a name.
!
!    level 0   the grammar        what may exist          THIS FILE
!    level 1   the calculus       how quantities relate
!    level 2   the minimization   how answers are found
!    level 3   the constitution   what the material says
!    level 4   the statement      what is asked
!
! Level 1 concretes these roles into the named citizens of graph
! mathematics, one bound parameter per name:
!
!    graph -----[edge count = 0]----------> support
!    field -----[domain size = 1]---------> functional
!    operation -[answer on the one entry]-> reduction
!    operation -[input from the one entry]> broadcast
!    transform -[verb = partition | assemble | coarsen | refine]
!
!=====================================================================!
!
!                         THE ADMISSION LAW
!
! Nothing enters this contract, at any level, except through four
! rules:
!
!    ABSORPTION     a choice from a finite list (vertex or edge, sum
!                   or max, real or complex) rides as an argument or
!                   a constant, never as a type. A type is earned
!                   only when the ROLE of an argument changes.
!
!    GENERATION     a question whose answer composes from other
!                   answers earns no procedure. Only generators
!                   enter; compositions are written at call sites.
!
!    COMMUTATION    two choices are independent axes only if applying
!                   them in either order gives one result. If the
!                   order matters they are one coupled concern and
!                   are handled in one place.
!
!    INHABITATION   no abstract type without a concrete citizen that
!                   answers every symbol meaningfully. Empty types
!                   and dead procedures are structural defects.
!
! Under these rules the minimum is exactly the four roles below. A
! fifth role always either absorbs into one of the four or forces
! dead procedures somewhere. The census of this level: four roles,
! fifty-five operation symbols.
!
!=====================================================================!
!
!                           THE FOUR ROLES
!
! Every citizen of every level is one of these, and each answers a
! different question:
!
!    graph ............ structure     what is joined to what
!    graph_field ...... value         what the members carry
!    graph_operation .. verb within   how data becomes other data
!    graph_transform .. verb between  how one graph becomes another
!
! The definition is self-hosting. A graph is two member sets and two
! integer edge fields, tail and head. A member set is itself a graph
! with no edges. A field's domain is such a graph. Structure is
! data, and the grammar closes on itself.
!
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
!=====================================================================!
!
!                           THE THREE RULES
!
! Three questions come up over and over in the procedures below:
! where does a value sit in a field, can a graph change after it is
! built, and what does apply do to a buffer it is handed? Each has
! one answer, fixed here, and the whole tower is written assuming it.
!
! WHERE A VALUE SITS. Suppose a domain lists three cells in the
! order 7, 3, 11, and a field on that domain holds two components
! per entry. The field keeps its values in that same order, with the
! components of each entry side by side:
!
!        cell            7        7        3        3      ...
!        component       1        2        1        2
!                     +--------+--------+--------+--------+--
!        values       |  v(1)  |  v(2)  |  v(3)  |  v(4)  |
!                     +--------+--------+--------+--------+--
!
! So if a cell is the p-th entry of the domain, its component c is
! the number at position
!
!        (p - 1) * num_components + c
!
! and anyone holding the flat vector finds any value by this formula
! alone. Once the layout is fixed, degree counting and vector
! numbering are one-line formulas, and a formula needs no procedure.
!
! CAN A GRAPH CHANGE? No. Everything a graph holds - structure,
! tags, its relation to the whole it came from - goes in at
! construction, and no procedure below accepts data afterwards.
! When an operation computes something new, the result leaves through
! that operation's output argument. The reason is repeatability: ask
! a graph the same question twice and it gives the same answer twice,
! no matter what ran in between.
!
! WHAT apply DOES TO A LENT BUFFER. It writes the result into it and
! never adds to what was there. The output argument is intent(inout)
! for one reason only: a caller that already holds a buffer of the
! right shape can lend it and save an allocation. Lending changes the
! cost of the call, not its meaning.
!
!=====================================================================!

module graph_grammar

  use iso_fortran_env, only : dp => REAL64

  implicit none

  private

  public :: graph
  public :: graph_field
  public :: graph_operation
  public :: graph_transform

  public :: GRAPH_FIELD_INTEGER
  public :: GRAPH_FIELD_REAL
  public :: GRAPH_FIELD_COMPLEX
  public :: GRAPH_FIELD_LOGICAL
  public :: GRAPH_FIELD_CHARACTER

  ! The kind of value a field holds. Five roads, one field: an
  ! absorbed axis. Integer for a colouring or a part number, real for
  ! the ordinary state, complex for a complex-step derivative,
  ! logical for a mask, character for boundary and material names.
  integer, parameter :: GRAPH_FIELD_INTEGER   = 1
  integer, parameter :: GRAPH_FIELD_REAL      = 2
  integer, parameter :: GRAPH_FIELD_COMPLEX   = 3
  integer, parameter :: GRAPH_FIELD_LOGICAL   = 4
  integer, parameter :: GRAPH_FIELD_CHARACTER = 5

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
  ! GRAPH_FIELD. The carrier of values.
  !
  ! A field is a function: values over a domain, and the domain is a
  ! graph - a member set of some host. One value kind per field,
  ! num_components values per entry, laid out by the formula in the
  ! header.
  !
  !      field  ---get--->  [ v1 v2 v3 v4 ]  ---> a solver, a file
  !                                               writer, an outside
  !      field  <--set----  [ v1 v2 v3 v4 ]  <--- library; the answer
  !                                               coming back
  !
  ! Fetch once, work in plain arrays where the arithmetic is free,
  ! write back once. Scaling and adding are not graph theory and have
  ! no procedures here on purpose. One get/set pair per value kind:
  ! the five roads of one absorbed axis.
  !
  ! A field whose domain has ONE entry is a single number. Level 1
  ! names that case the functional; nothing in this role depends on
  ! the size of the domain, which is why the name can wait.
  !
  ! num_entries repeats the size of the domain. That repetition is a
  ! priced convenience for call sites, recorded here so nobody
  ! mistakes it for a generator.
  !===================================================================!

  type, abstract :: graph_field

   contains

     ! Identity: what it is called and what unit it carries.
     procedure(field_name_interface), deferred :: name
     procedure(field_name_interface), deferred :: units

     ! Where it lives and how it is shaped.
     procedure(field_domain_interface), deferred :: domain
     procedure(field_count_interface) , deferred :: num_components
     procedure(field_count_interface) , deferred :: num_entries
     procedure(field_count_interface) , deferred :: value_kind

     ! The plain-vector adapters, one pair per value kind.
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

  !===================================================================!
  ! GRAPH_OPERATION. The verb within a graph: data in, data out.
  !
  !      input_graph, input_data(:)  --- apply --->  output
  !
  ! Three symbols. name says what it is. domain says where the answer
  ! lives: a member set of the input graph. apply does the work,
  ! under the lent-buffer rule of the header.
  !
  ! A concrete operation is handed the fields it reads when it is
  ! constructed - a coefficient, a measure, a geometry field arrives
  ! as an argument the compiler checks, and apply fetches nothing by
  ! name. This is what keeps the call structure visible from
  ! construction to result.
  !
  ! The concretions of this role, arriving on level 1, split by
  ! order: first-order kernels act on fields (a differential
  ! operator, a balance), and the one higher-order citizen acts on
  ! operators (the walk, a traversal whose kernel is not yet bound).
  ! The maps that touch the one-entry domain - the reduction and the
  ! broadcast - stand beside this role rather than under it, one
  ! deliberate deviation, recorded where they are declared.
  !===================================================================!

  type, abstract :: graph_operation

   contains

     procedure(operation_name_interface)  , deferred :: name
     procedure(operation_domain_interface), deferred :: domain
     procedure(operation_apply_interface) , deferred :: apply

  end type graph_operation

  !===================================================================!
  ! GRAPH_TRANSFORM. The verb between graphs.
  !
  ! Two symbols at this level, both admissibility questions: may this
  ! transform act on that graph, and on that data riding on it. The
  ! maps themselves are verbs, and verbs are level-1 concretions -
  ! partition, assemble, coarsen, refine - each pair judged by its
  ! round-trip law:
  !
  !      exact        assemble(partition(G)) = G      both ways
  !      one-sided    coarsen(refine(G)) = G          one way only
  !
  ! The central law of the whole tower is a sentence about this role:
  ! split a graph into parts, work on the parts, put the answer back
  ! together.
  !
  !        G'  =  P^-1 ( A ( P ( G ) ) )
  !
  !             G            P(G)           A(P(G))          G'
  !             |             |                |              |
  !             +-----P------>+-------A------->+----P^-1----->+
  !             |             |                |              |
  !        whole graph    the parts      worked-on parts  whole again
  !
  ! P is a transform and nothing else. P^-1 is a transform and
  ! nothing else. A is an operation and does neither.
  !===================================================================!

  type, abstract :: graph_transform

   contains

     procedure(transform_on_graph_interface), deferred :: defined_on_graph
     procedure(transform_on_data_interface) , deferred :: defined_on_data

  end type graph_transform

  abstract interface

     !===============================================================!
     ! Structure: identity, counts, incidence.
     !===============================================================!

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

     !===============================================================!
     ! The named sets. Called once, when an operation begins, so each
     ! answer is a whole graph - the edgeless graph of the members -
     ! and the cost is paid once per sweep, not once per cell.
     !===============================================================!

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

     !===============================================================!
     ! Neighbourhoods. Bare indices, pure, loop-safe.
     !===============================================================!

     pure subroutine graph_from_vertex_interface(this, vertex_index, indices)
       import :: graph
       class(graph), intent(in) :: this
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
       import :: graph
       class(graph), intent(in) :: this
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
       import :: graph
       class(graph), intent(in) :: this
       integer, intent(in) :: index
     end function graph_global_id_interface

     !---------------------------------------------------------------!
     ! The same map read backwards: where does whole-graph vertex 4
     ! sit inside part 2? The part is named because the partition is
     ! replicated - every image builds every part, so one image
     ! regularly queries a part it does not own.
     !---------------------------------------------------------------!

     pure integer function graph_part_id_interface(this, global_index, part_id)
       import :: graph
       class(graph), intent(in) :: this
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
       import :: graph
       class(graph), intent(in) :: this
       integer, intent(in) :: index
     end function graph_owner_part_interface

     !===============================================================!
     ! Value: identity, domain, shape.
     !===============================================================!

     pure function field_name_interface(this) result(name)
       import :: graph_field
       class(graph_field), intent(in) :: this
       character(len=:), allocatable :: name
     end function field_name_interface

     subroutine field_domain_interface(this, domain)
       import :: graph_field, graph
       class(graph_field), intent(in) :: this
       class(graph), allocatable, intent(out) :: domain
     end subroutine field_domain_interface

     pure integer function field_count_interface(this)
       import :: graph_field
       class(graph_field), intent(in) :: this
     end function field_count_interface

     !===============================================================!
     ! Value: the plain-vector adapters, one pair per kind.
     !===============================================================!

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
     ! Verb within: name, domain, apply.
     !===============================================================!

     pure function operation_name_interface(this) result(name)
       import :: graph_operation
       class(graph_operation), intent(in) :: this
       character(len=:), allocatable :: name
     end function operation_name_interface

     subroutine operation_domain_interface(this, input_graph, domain)
       import :: graph_operation, graph
       class(graph_operation), intent(in) :: this
       class(graph), intent(in) :: input_graph
       class(graph), allocatable, intent(out) :: domain
     end subroutine operation_domain_interface

     subroutine operation_apply_interface(this, input_graph, input_data, output)
       import :: graph_operation, graph, graph_field
       class(graph_operation), intent(in) :: this
       class(graph), intent(in) :: input_graph
       class(graph_field), intent(in), optional :: input_data(:)
       class(graph_field), allocatable, intent(inout) :: output
     end subroutine operation_apply_interface

     !===============================================================!
     ! Verb between: the two admissibility questions.
     !===============================================================!

     pure logical function transform_on_graph_interface(this, input_graph)
       import :: graph_transform, graph
       class(graph_transform), intent(in) :: this
       class(graph), intent(in) :: input_graph
     end function transform_on_graph_interface

     pure logical function transform_on_data_interface(this, input_graph, input_data)
       import :: graph_transform, graph, graph_field
       class(graph_transform), intent(in) :: this
       class(graph), intent(in) :: input_graph
       class(graph_field), intent(in) :: input_data
     end function transform_on_data_interface

  end interface

end module graph_grammar
