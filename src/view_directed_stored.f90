!=====================================================================!
! A concrete graph that stores its structure.
!
! Hand it a vertex count and an edge list and it works out, once, the
! neighbour lists every later query reads:
!
!            tails  1 1 2 3            e1: 1 -> 2
!            heads  2 3 4 4            e2: 1 -> 3
!                                      e3: 2 -> 4
!                                      e4: 3 -> 4
!
!                        (1)
!                       /   \
!                      v     v
!                    (2)     (3)
!                      \     /
!                       v   v
!                        (4)
!
! An edge whose head is not a real vertex has no head at all. That is
! a boundary face: it is attached to one cell alone, and no imaginary
! cell is invented on the far side of the wall.
!
! Four compressed lists are built at construction and never rebuilt -
! edges touching a vertex, vertices next to it, and the same two split
! by which way the edges point. Every walking query is then a slice of
! an array, which is what lets those queries stay pure and cheap
! enough to sit inside a loop over a million cells.
!
!=====================================================================!
!
!                       WHAT THIS GRAPH DOES NOT DO
!
! It holds no geometry, no physics, no solver state, and no algorithm.
! Colouring, traversal order, partitioning and the rest are operations
! and transforms that read a graph; they are not things a graph does.
! Each convenience procedure added here erodes that separation, one
! procedure at a time.
!
! IT CARRIES NO VALUES. A field references its domain; the reference
! never points the other way. What an operation reads, it is handed
! at construction, as arguments the compiler can see. The one string
! kept is the tag, which is data: the mesh file named its boundary
! groups, and those names flow in from outside the code.
!
! READ ONLY. No procedure puts data on a graph after construction.
! Anything computed leaves as an operation's output. Without this rule
! the graph accumulates state, and its answers come to depend on the
! order of past calls rather than on the mesh it was built from.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module view_directed_stored

  use view_directed, only : directed_graph, forward, reverse
  use relation_algorithms, only : topological_order
  use graph_fractal      , only : graph
  use relation_binary, only : group_by_key, csr_relation
  use relation_partition, only : partition_relation
  use view_directed     , only : SIDE_VERTEX, SIDE_EDGE
  use map_set_representation, only : counted_set_representation
  use map_set      , only : set_map
  use map_label    , only : label_map
  use map_inclusion, only : inclusion_map
  use map_carving  , only : carve

  implicit none

  private

  integer, parameter :: SELECT_INTERIOR = 1, SELECT_BOUNDARY = 2, SELECT_TAGGED = 3
  public :: stored_directed_graph

  !===================================================================!
  ! A graph that keeps its own structure in arrays.
  !===================================================================!

  type, extends(directed_graph) :: stored_directed_graph

     integer :: number = 1
     integer :: nv     = 0
     integer :: ne     = 0

     !----------------------------------------------------------------!
     ! THE ORIENTATION. D = (V, E, tail, head) and its transpose
     ! D^T = (V, E, head, tail) are one stored object read two ways:
     ! both endpoint lists and both compressed directions are kept,
     ! and reversed says which is which. The transpose of the
     ! transpose is the object itself, exactly, and the identity is
     ! the same relation's.
     !----------------------------------------------------------------!
     logical :: reversed = .false.

     !----------------------------------------------------------------!
     ! Edge endpoints. A head of zero means the edge has none.
     !----------------------------------------------------------------!

     integer, allocatable :: tail(:)
     integer, allocatable :: head(:)

     !----------------------------------------------------------------!
     ! The compressed lists, built once. Vertex v's incident edges are
     ! einc(xinc(v) : xinc(v+1)-1), and so on for the other three.
     !----------------------------------------------------------------!

     integer, allocatable :: xinc(:), einc(:)
     integer, allocatable :: xadj(:), vadj(:)
     integer, allocatable :: xout(:), eout(:)
     integer, allocatable :: xin(:) , ein(:)

     !----------------------------------------------------------------!
     ! Names carried by vertices and edges. Blank means untagged.
     !----------------------------------------------------------------!

     character(len=:), allocatable :: vtag(:)
     character(len=:), allocatable :: etag(:)

     !----------------------------------------------------------------!
     ! HOW THIS GRAPH RELATES TO A WHOLE ONE, and it is one component
     ! now rather than six: the relation r <= S_part x S_whole. The
     ! record is a REPRESENTATION - part-local numbering, ownership,
     ! provenance - and none of it is a question about
     ! D = (V, E, tail, head), which is why it stopped being a binding
     ! on the contract and became a value the graph carries.
     !
     ! A graph straight off a mesh file carries the IDENTITY relation:
     ! one part, itself, every member its own name.
     !----------------------------------------------------------------!

     type(partition_relation) :: whole_rel

     !----------------------------------------------------------------!
     ! The graph's two carriers (AGENTS.md, phase 1): its vertices
     ! and its edges as declared set GRAPHS, stamped once at
     ! construction and handed out beside the old vocabulary. Every
     ! call to vertex_set answers the SAME domain, so a relation
     ! signature can hold onto the identity.
     !
     ! Identity only. The extension is 1..nv and 1..ne, which the
     ! graph already answers through num_vertices/num_edges, so
     ! storing a representation here would store a second copy of a
     ! fact - and storing a map would make the graph a registry of
     ! interpretations, which it is not.
     !----------------------------------------------------------------!

     type(graph) :: vset
     type(graph) :: eset

   contains

     !----------------------------------------------------------------!
     ! Identity and size.
     !----------------------------------------------------------------!

     procedure :: id
     procedure :: num_vertices
     procedure :: num_edges

     !----------------------------------------------------------------!
     ! The carriers, beside the old vocabulary (AGENTS.md, phase 1).
     !----------------------------------------------------------------!

     procedure :: vertex_set
     procedure :: edge_set
     procedure :: name_carriers

     !----------------------------------------------------------------!
     ! Where an edge goes.
     !----------------------------------------------------------------!

     procedure :: edge_tail
     procedure :: edge_head
     procedure :: edge_has_head
     procedure :: transpose
     procedure :: transposed
     procedure :: loop

     !----------------------------------------------------------------!
     ! The named vertex sets.
     !----------------------------------------------------------------!

     procedure :: interior_vertices
     procedure :: boundary_vertices
     procedure :: tagged_vertices

     !----------------------------------------------------------------!
     ! The named edge sets.
     !----------------------------------------------------------------!

     procedure :: interior_edges
     procedure :: boundary_edges
     procedure :: tagged_edges

     !----------------------------------------------------------------!
     ! The named sets of one part.
     !----------------------------------------------------------------!

     procedure :: owned_vertices
     procedure :: borrowed_vertices
     procedure :: overlap_vertices
     procedure :: owned_edges
     procedure :: borrowed_edges
     procedure :: overlap_edges

     !----------------------------------------------------------------!
     ! Walking, without regard to direction and with it.
     !----------------------------------------------------------------!

     procedure :: incident_edges
     procedure :: adjacent_vertices
     procedure :: outgoing_edges
     procedure :: incoming_edges
     procedure :: outgoing_vertices
     procedure :: incoming_vertices

     !----------------------------------------------------------------!
     ! How a part relates to the whole: ONE accessor, handing back the
     ! relation by value. The eight questions that used to stand here
     ! are r's, and a caller that needs them takes r and asks it -
     ! which is also what lets the four verbs be handed one explicitly.
     !
     ! whole_relation, not relation: this graph CONTAINS relations -
     ! its incidence and its adjacency are two of them - and this is
     ! not one of those. It is the relation to the whole.
     !----------------------------------------------------------------!

     procedure :: whole_relation

     !----------------------------------------------------------------!
     ! The structure read as relations (AGENTS.md section 16):
     ! T <= E x V (edge to tail) and H <= E x V (edge to head; a
     ! boundary edge is an absence in H). Derived from the stored
     ! table when asked, so a pattern graph or a part graph that
     ! nobody reads relationally never pays for them - the
     ! section-66 benchmark caught the eager version costing every
     ! construction 2.2x.
     !----------------------------------------------------------------!

     procedure :: tail_relation
     procedure :: head_relation

  end type stored_directed_graph

  !===================================================================!
  ! Constructor.
  !===================================================================!

  interface stored_directed_graph
     module procedure create
  end interface stored_directed_graph

contains

  !===================================================================!
  ! Build a graph from a vertex count and an edge list.
  !
  ! A head of zero (or anything outside 1..nv) means the edge has no
  ! head - a boundary face. Tags are optional; an untagged graph
  ! returns an empty set for every tagged query.
  !===================================================================!

  type(stored_directed_graph) function create(nv, tails, heads, vtags, etags, &
       &                             number, vglobal, vowner, eglobal, &
       &                             eowner, num_parts, whole_vertices, whole_edges, &
       &                             num_whole_vertices, num_whole_edges) result(this)

    integer           , intent(in)           :: nv
    integer           , intent(in)           :: tails(:)
    integer           , intent(in)           :: heads(:)
    character(len=*)  , intent(in), optional :: vtags(:)
    character(len=*)  , intent(in), optional :: etags(:)
    integer           , intent(in), optional :: number

    !----------------------------------------------------------------!
    ! The tuples of r, for a graph that is a piece of a larger one:
    ! what each of its own members is called in the whole, and which
    ! part owns it. These arrive HERE or not at all - a graph that
    ! could be told its relation afterwards would answer the same
    ! question two ways in one lifetime, which is the one thing the
    ! grammar says a graph may never do. Present vglobal is what makes
    ! a graph a piece; absent, it is a whole.
    !----------------------------------------------------------------!

    integer           , intent(in), optional :: vglobal(:)
    integer           , intent(in), optional :: vowner(:)
    integer           , intent(in), optional :: eglobal(:)
    integer           , intent(in), optional :: eowner(:)
    integer           , intent(in), optional :: num_parts

    !----------------------------------------------------------------!
    ! The far side of r, by identity and count. A relation that could
    ! not say WHICH sets it relates would let a caller holding two of
    ! them hand the wrong one to the wrong graph and be wrong in
    ! silence.
    !----------------------------------------------------------------!

    type(graph)   , intent(in), optional :: whole_vertices, whole_edges
    integer           , intent(in), optional :: num_whole_vertices, num_whole_edges

    integer :: e

    this % nv = nv
    this % ne = size(tails)

    ! Declare the two domains once, here, so every later answer
    ! carries one identity per side for this graph's whole life.
    call this % vset % declare()
    call this % eset % declare()

    if (present(number)) this % number = number

    allocate(this % tail, source=tails)
    allocate(this % head(this % ne))

    ! Normalize every missing head to zero, so one test answers
    ! everywhere afterwards.
    do e = 1, this % ne
       if (heads(e) >= 1 .and. heads(e) <= nv) then
          this % head(e) = heads(e)
       else
          this % head(e) = 0
       end if
    end do

    !----------------------------------------------------------------!
    ! The relation goes in through the door, or the identity relation
    ! does. A graph that could be told its relation afterwards would
    ! answer one question two ways in one lifetime; present vglobal is
    ! what makes a graph a piece, and absent, it is a whole.
    !----------------------------------------------------------------!

    if (present(vglobal)) then
       this % whole_rel = partition_relation( &
            & part_vertices    = this % vset, num_part_vertices = this % nv, &
            & part_edges    = this % eset, num_part_edges = this % ne, &
            & whole_vertices   = merge_set(whole_vertices, this % vset),   &
            & num_whole_vertices = merge_count(num_whole_vertices, this % nv), &
            & whole_edges   = merge_set(whole_edges, this % eset),   &
            & num_whole_edges = merge_count(num_whole_edges, this % ne), &
            & number  = this % number,                               &
            & num_parts  = merge_count(num_parts, 1),                      &
            & vglobal = vglobal,                                     &
            & vowner  = pick_owner(vowner, size(vglobal), this % number), &
            & eglobal = pick_global(eglobal, this % ne),             &
            & eowner  = pick_owner(eowner, this % ne, this % number))
    else
       this % whole_rel = partition_relation( &
            & this % vset, this % nv, this % eset, this % ne)
    end if

    if (present(vtags)) allocate(this % vtag, source=vtags)
    if (present(etags)) allocate(this % etag, source=etags)

    ! Everything the mesh knew arrives here and never changes again.

    call build_incidence(this % nv, this % tail, this % head, &
         &               this % xinc, this % einc)

    call build_adjacency(this % nv, this % tail, this % head, &
         &               this % xinc, this % einc, this % xadj, this % vadj)

    call build_directed(this % nv, this % tail, this % xout, this % eout)
    call build_directed(this % nv, this % head, this % xin , this % ein )

  end function create

  !===================================================================!
  ! Every edge touches its tail, and its head when it has one. Count
  ! first, then fill: two passes and no growing arrays.
  !===================================================================!

  pure subroutine build_incidence(nv, tail, head, xptr, elist)

    integer             , intent(in)  :: nv
    integer             , intent(in)  :: tail(:), head(:)
    integer, allocatable, intent(out) :: xptr(:), elist(:)

    integer, allocatable :: keys(:), values(:)
    integer :: ne, e

    ! one (endpoint, edge) pair per end, interleaved tail-then-head
    ! so each vertex's fibre keeps the single-pass edge order; a
    ! missing head is key zero and belongs to no vertex
    ne = size(tail)
    allocate(keys(2 * ne), values(2 * ne))
    do e = 1, ne
       keys(2 * e - 1) = tail(e)
       keys(2 * e)     = head(e)
       values(2 * e - 1) = e
       values(2 * e)     = e
    end do

    call group_by_key(nv, keys, values, xptr, elist)

  end subroutine build_incidence

  !===================================================================!
  ! A vertex's neighbours are the far ends of the edges touching it,
  ! each counted once however many edges join the pair.
  !===================================================================!

  pure subroutine build_adjacency(nv, tail, head, xinc, einc, xptr, vlist)

    integer             , intent(in)  :: nv
    integer             , intent(in)  :: tail(:), head(:)
    integer             , intent(in)  :: xinc(:), einc(:)
    integer, allocatable, intent(out) :: xptr(:), vlist(:)

    integer, allocatable :: seen(:), scratch(:)
    integer :: v, k, e, other, ndistinct, total

    allocate(xptr(nv + 1))
    allocate(seen(nv))
    seen = 0

    ! First pass counts the distinct neighbours of every vertex.
    xptr(1) = 1
    do v = 1, nv
       ndistinct = 0
       do k = xinc(v), xinc(v + 1) - 1
          e = einc(k)
          other = far_end(tail(e), head(e), v)
          if (other >= 1 .and. other /= v) then
             if (seen(other) /= v) then
                seen(other) = v
                ndistinct = ndistinct + 1
             end if
          end if
       end do
       xptr(v + 1) = xptr(v) + ndistinct
    end do

    total = xptr(nv + 1) - 1
    allocate(vlist(max(total, 0)))

    ! Second pass writes them, with the marker reset so the same test
    ! can run again.
    seen = 0
    allocate(scratch, source=xptr(1:nv))
    do v = 1, nv
       do k = xinc(v), xinc(v + 1) - 1
          e = einc(k)
          other = far_end(tail(e), head(e), v)
          if (other >= 1 .and. other /= v) then
             if (seen(other) /= v) then
                seen(other) = v
                vlist(scratch(v)) = other
                scratch(v) = scratch(v) + 1
             end if
          end if
       end do
    end do

  end subroutine build_adjacency

  !===================================================================!
  ! Given both ends of an edge and one of them, name the other.
  ! Answers zero when the edge has no head, which is how a boundary
  ! face reports that there is nothing on the far side.
  !===================================================================!

  pure integer function far_end(tail, head, here)

    integer, intent(in) :: tail, head, here

    if (tail == here) then
       far_end = head
    else
       far_end = tail
    end if

  end function far_end

  !===================================================================!
  ! Group the edges by one of their endpoints - tails to get the
  ! outgoing lists, heads to get the incoming ones. An endpoint of
  ! zero belongs to no vertex and is skipped.
  !===================================================================!

  pure subroutine build_directed(nv, endpoint, xptr, elist)

    integer             , intent(in)  :: nv
    integer             , intent(in)  :: endpoint(:)
    integer, allocatable, intent(out) :: xptr(:), elist(:)

    integer, allocatable :: identity(:)
    integer :: e

    ! group the edges by the chosen endpoint; an endpoint of zero
    ! belongs to no vertex and is skipped by the kernel
    allocate(identity(size(endpoint)))
    identity = [(e, e = 1, size(endpoint))]

    call group_by_key(nv, endpoint, identity, xptr, elist)

  end subroutine build_directed

  !===================================================================!
  ! Identity and size.
  !===================================================================!

  pure integer function id(this)

    class(stored_directed_graph), intent(in) :: this

    id = this % number

  end function id

  !===================================================================!
  ! How many vertices.
  !===================================================================!

  pure integer function num_vertices(this)

    class(stored_directed_graph), intent(in) :: this

    num_vertices = this % nv

  end function num_vertices

  !===================================================================!
  ! The two carriers, as declared at construction. Copies of one
  ! stamped domain: every call answers a set that same_as agrees is
  ! the same set (AGENTS.md, phase 1).
  !===================================================================!

  type(graph) function vertex_set(this)

    class(stored_directed_graph), intent(in) :: this

    vertex_set = this % vset

  end function vertex_set

  type(graph) function edge_set(this)

    class(stored_directed_graph), intent(in) :: this

    edge_set = this % eset

  end function edge_set

  !===================================================================!
  ! What this graph calls its own two domains, bound into the CALLER'S
  ! label map. The graph knows the names; it does not keep the map, so
  ! a caller that names nothing never calls this and carries nothing.
  !===================================================================!

  subroutine name_carriers(this, labels)

    class(stored_directed_graph), intent(in)    :: this
    type(label_map)    , intent(inout) :: labels

    call labels % bind(this % vset, 'vertices')
    call labels % bind(this % eset, 'edges')

  end subroutine name_carriers

  !===================================================================!
  ! How many edges.
  !===================================================================!

  pure integer function num_edges(this)

    class(stored_directed_graph), intent(in) :: this

    num_edges = this % ne

  end function num_edges

  !===================================================================!
  ! Where an edge goes.
  !===================================================================!
  ! THE LOOP over the graph: its vertices in an order every edge
  ! respects - a tail before its head in the forward orientation, a
  ! head before its tail in reverse, which is the loop over the
  ! transpose. It exists only where the graph has no cycle; a graph
  ! with one has no loop, and the request stops the program. A march
  ! is this loop forward and its adjoint this loop in reverse, and
  ! neither writes an instant's number down. The order is the one
  ! topological sort in the tree, over the graph's own adjacency.
  !===================================================================!

  function loop(this, orientation) result(order)

    class(stored_directed_graph), intent(in) :: this
    integer, intent(in), optional :: orientation
    integer, allocatable :: order(:)

    type(set_map)      :: sets
    type(csr_relation) :: adjacency
    integer, allocatable :: table(:,:)
    logical :: acyclic
    integer :: way, e, n

    way = forward
    if (present(orientation)) way = orientation
    if (way /= forward .and. way /= reverse) then
       error stop 'stored_directed_graph: a loop runs forward or in reverse'
    end if

    allocate(table(2, this % ne))
    n = 0
    do e = 1, this % ne
       if (.not. this % edge_has_head(e)) cycle
       n = n + 1
       if (way == forward) then
          table(:, n) = [this % edge_tail(e), this % edge_head(e)]
       else
          table(:, n) = [this % edge_head(e), this % edge_tail(e)]
       end if
    end do

    call sets % bind(this % vset, counted_set_representation(this % nv))
    adjacency = csr_relation('adjacency', this % vset, this % vset, table(:, 1:n), sets)
    call topological_order(adjacency, sets, order, acyclic)
    if (.not. acyclic) then
       error stop 'stored_directed_graph: a graph with a cycle has no loop'
    end if

  end function loop

  !===================================================================!
  ! The transpose: the same object read the other way, every edge's
  ! tail its head and head its tail, so that transposing twice gives
  ! back what was there. An edge without a head would become an edge
  ! without a tail, which is not an edge; such a graph has no
  ! transpose and the request stops the program.
  !===================================================================!

  type(stored_directed_graph) function transpose(this) result(turned)

    class(stored_directed_graph), intent(in) :: this

    if (any(this % head < 1)) then
       error stop 'stored_directed_graph: a graph with an edge without a head has no transpose'
    end if

    turned = this
    turned % reversed = .not. this % reversed

  end function transpose

  pure logical function transposed(this)

    class(stored_directed_graph), intent(in) :: this

    transposed = this % reversed

  end function transposed

  !===================================================================!

  pure integer function edge_tail(this, edge_index)

    class(stored_directed_graph), intent(in) :: this
    integer            , intent(in) :: edge_index

    if (this % reversed) then
       edge_tail = this % head(edge_index)
    else
       edge_tail = this % tail(edge_index)
    end if

  end function edge_tail

  !===================================================================!
  ! The vertex an edge enters:  (i) --e--> (j)  answers j. A
  ! boundary edge enters nothing and answers zero.
  !===================================================================!

  pure integer function edge_head(this, edge_index)

    class(stored_directed_graph), intent(in) :: this
    integer            , intent(in) :: edge_index

    if (this % reversed) then
       edge_head = this % tail(edge_index)
    else
       edge_head = this % head(edge_index)
    end if

  end function edge_head

  !===================================================================!
  ! Whether the edge enters a vertex at all; false marks a boundary
  ! edge.
  !===================================================================!

  pure logical function edge_has_head(this, edge_index)

    class(stored_directed_graph), intent(in) :: this
    integer            , intent(in) :: edge_index

    edge_has_head = this % edge_head(edge_index) >= 1

  end function edge_has_head

  !===================================================================!
  ! The named vertex subsets. A boundary vertex is one that touches a
  ! boundary edge; an interior vertex is one that does not.
  !===================================================================!

  !===================================================================!
  ! The vertices that touch no boundary edge.
  !===================================================================!

  subroutine interior_vertices(this, sets, labels, inclusions, members)

    class(stored_directed_graph), intent(in)    :: this
    type(set_map)      , intent(inout) :: sets
    type(label_map)    , intent(inout) :: labels
    type(inclusion_map), intent(inout) :: inclusions
    type(graph)    , intent(out)   :: members

    call carve(members, selected(this, SELECT_INTERIOR, .true.), 'interior_vertices', &
         & this % vset, sets, labels, inclusions)

  end subroutine interior_vertices

  !===================================================================!
  ! The vertices that touch a boundary edge.
  !===================================================================!

  subroutine boundary_vertices(this, sets, labels, inclusions, members)

    class(stored_directed_graph), intent(in)    :: this
    type(set_map)      , intent(inout) :: sets
    type(label_map)    , intent(inout) :: labels
    type(inclusion_map), intent(inout) :: inclusions
    type(graph)    , intent(out)   :: members

    call carve(members, selected(this, SELECT_BOUNDARY, .true.), 'boundary_vertices', &
         & this % vset, sets, labels, inclusions)

  end subroutine boundary_vertices

  !===================================================================!
  ! The vertices carrying this tag - a mesh's named patches arrive
  ! here.
  !===================================================================!

  subroutine tagged_vertices(this, tag, sets, labels, inclusions, members)

    class(stored_directed_graph), intent(in)    :: this
    character(len=*)   , intent(in)    :: tag
    type(set_map)      , intent(inout) :: sets
    type(label_map)    , intent(inout) :: labels
    type(inclusion_map), intent(inout) :: inclusions
    type(graph)    , intent(out)   :: members

    call carve(members, selected(this, SELECT_TAGGED, .true., tag), 'tagged_vertices', &
         & this % vset, sets, labels, inclusions)

  end subroutine tagged_vertices


  !===================================================================!
  ! Does any edge touching this vertex stop here rather than holding
  ! on to another vertex?
  !===================================================================!

  !===================================================================!
  ! The members selected by one predicate, on either side: interior
  ! (a vertex touching no boundary edge; an edge with a head),
  ! boundary (the complement), or tagged with a name. Three
  ! predicates times two sides, written once.
  !===================================================================!

  pure function selected(this, which, on_vertices, tag) result(pick)

    class(stored_directed_graph), intent(in) :: this
    integer                     , intent(in) :: which
    logical                     , intent(in) :: on_vertices
    character(len=*), optional  , intent(in) :: tag
    integer, allocatable :: pick(:)

    integer :: i, n, k
    logical :: keep

    n = merge(this % nv, this % ne, on_vertices)
    allocate(pick(n))
    k = 0
    do i = 1, n
       select case (which)
       case (SELECT_INTERIOR, SELECT_BOUNDARY)
          if (on_vertices) then
             keep = .not. touches_boundary(this, i)
          else
             keep = this % edge_has_head(i)
          end if
          if (which == SELECT_BOUNDARY) keep = .not. keep
       case (SELECT_TAGGED)
          keep = .false.
          if (.not. present(tag)) error stop 'stored_directed_graph: a tagged selection names its tag'
          if (on_vertices) then
             if (allocated(this % vtag)) keep = trim(this % vtag(i)) == tag
          else
             if (allocated(this % etag)) keep = trim(this % etag(i)) == tag
          end if
       case default
          error stop 'stored_directed_graph: a selection is interior, boundary or tagged'
       end select
       if (keep) then
          k = k + 1
          pick(k) = i
       end if
    end do
    pick = pick(1:k)

  end function selected

  pure logical function touches_boundary(this, v)

    class(stored_directed_graph), intent(in) :: this
    integer            , intent(in) :: v

    integer :: k

    touches_boundary = .false.
    do k = this % xinc(v), this % xinc(v + 1) - 1
       if (.not. this % edge_has_head(this % einc(k))) then
          touches_boundary = .true.
          return
       end if
    end do

  end function touches_boundary

  !===================================================================!
  ! The named edge subsets. A boundary edge is one with no head.
  !===================================================================!

  !===================================================================!
  ! The edges with a head: both ends real.
  !===================================================================!

  subroutine interior_edges(this, sets, labels, inclusions, members)

    class(stored_directed_graph), intent(in)    :: this
    type(set_map)      , intent(inout) :: sets
    type(label_map)    , intent(inout) :: labels
    type(inclusion_map), intent(inout) :: inclusions
    type(graph)    , intent(out)   :: members

    call carve(members, selected(this, SELECT_INTERIOR, .false.), 'interior_edges', &
         & this % eset, sets, labels, inclusions)

  end subroutine interior_edges

  !===================================================================!
  ! The edges with no head - the open ends of the graph.
  !===================================================================!

  subroutine boundary_edges(this, sets, labels, inclusions, members)

    class(stored_directed_graph), intent(in)    :: this
    type(set_map)      , intent(inout) :: sets
    type(label_map)    , intent(inout) :: labels
    type(inclusion_map), intent(inout) :: inclusions
    type(graph)    , intent(out)   :: members

    call carve(members, selected(this, SELECT_BOUNDARY, .false.), 'boundary_edges', &
         & this % eset, sets, labels, inclusions)

  end subroutine boundary_edges

  !===================================================================!
  ! The edges carrying this tag - a mesh's named patches arrive
  ! here.
  !===================================================================!

  subroutine tagged_edges(this, tag, sets, labels, inclusions, members)

    class(stored_directed_graph), intent(in)    :: this
    character(len=*)   , intent(in)    :: tag
    type(set_map)      , intent(inout) :: sets
    type(label_map)    , intent(inout) :: labels
    type(inclusion_map), intent(inout) :: inclusions
    type(graph)    , intent(out)   :: members

    call carve(members, selected(this, SELECT_TAGGED, .false., tag), 'tagged_edges', &
         & this % eset, sets, labels, inclusions)

  end subroutine tagged_edges

  !===================================================================!
  ! The named sets of one part.
  !
  ! A graph that was never cut owns everything and borrows nothing,
  ! whichever part the query names. A partitioner fills in the owner
  ! arrays and these answers become real.
  !===================================================================!

  subroutine owned_vertices(this, part_id, sets, labels, inclusions, members)

    class(stored_directed_graph), intent(in)    :: this
    integer            , intent(in)    :: part_id
    type(set_map)      , intent(inout) :: sets
    type(label_map)    , intent(inout) :: labels
    type(inclusion_map), intent(inout) :: inclusions
    type(graph)    , intent(out)   :: members

    call carve(members, owner_matches(this % whole_rel, this % nv, part_id, .true., .true.), 'owned_vertices', &
         & this % vset, sets, labels, inclusions)

  end subroutine owned_vertices

  !===================================================================!
  ! The vertices this part reads but does not own - the neighbours'
  ! cells along the cut.
  !===================================================================!

  subroutine borrowed_vertices(this, part_id, sets, labels, inclusions, members)

    class(stored_directed_graph), intent(in)    :: this
    integer            , intent(in)    :: part_id
    type(set_map)      , intent(inout) :: sets
    type(label_map)    , intent(inout) :: labels
    type(inclusion_map), intent(inout) :: inclusions
    type(graph)    , intent(out)   :: members

    call carve(members, owner_matches(this % whole_rel, this % nv, part_id, .true., .false.), 'borrowed_vertices', &
         & this % vset, sets, labels, inclusions)

  end subroutine borrowed_vertices

  !===================================================================!
  ! The overlap is everything this part must see to finish what it
  ! owns: what it owns, plus what it borrows.
  !===================================================================!

  subroutine overlap_vertices(this, part_id, sets, labels, inclusions, members)

    class(stored_directed_graph), intent(in)    :: this
    integer            , intent(in)    :: part_id
    type(set_map)      , intent(inout) :: sets
    type(label_map)    , intent(inout) :: labels
    type(inclusion_map), intent(inout) :: inclusions
    type(graph)    , intent(out)   :: members

    integer, allocatable :: owned(:), borrowed(:)

    allocate(owned   , source=owner_matches(this % whole_rel, this % nv, part_id, .true., .true.))
    allocate(borrowed, source=owner_matches(this % whole_rel, this % nv, part_id, .true., .false.))

    call carve(members, [owned, borrowed], 'overlap_vertices', &
         & this % vset, sets, labels, inclusions)

  end subroutine overlap_vertices

  !===================================================================!
  ! The edges whose keeper is this part.
  !===================================================================!

  subroutine owned_edges(this, part_id, sets, labels, inclusions, members)

    class(stored_directed_graph), intent(in)    :: this
    integer            , intent(in)    :: part_id
    type(set_map)      , intent(inout) :: sets
    type(label_map)    , intent(inout) :: labels
    type(inclusion_map), intent(inout) :: inclusions
    type(graph)    , intent(out)   :: members

    call carve(members, owner_matches(this % whole_rel, this % ne, part_id, .false., .true.), 'owned_edges', &
         & this % eset, sets, labels, inclusions)

  end subroutine owned_edges

  !===================================================================!
  ! The edges this part reads but does not own.
  !===================================================================!

  subroutine borrowed_edges(this, part_id, sets, labels, inclusions, members)

    class(stored_directed_graph), intent(in)    :: this
    integer            , intent(in)    :: part_id
    type(set_map)      , intent(inout) :: sets
    type(label_map)    , intent(inout) :: labels
    type(inclusion_map), intent(inout) :: inclusions
    type(graph)    , intent(out)   :: members

    call carve(members, owner_matches(this % whole_rel, this % ne, part_id, .false., .false.), 'borrowed_edges', &
         & this % eset, sets, labels, inclusions)

  end subroutine borrowed_edges

  !===================================================================!
  ! Owned and borrowed together: every edge this part can see.
  !===================================================================!

  subroutine overlap_edges(this, part_id, sets, labels, inclusions, members)

    class(stored_directed_graph), intent(in)    :: this
    integer            , intent(in)    :: part_id
    type(set_map)      , intent(inout) :: sets
    type(label_map)    , intent(inout) :: labels
    type(inclusion_map), intent(inout) :: inclusions
    type(graph)    , intent(out)   :: members

    integer, allocatable :: owned(:), borrowed(:)

    allocate(owned   , source=owner_matches(this % whole_rel, this % ne, part_id, .false., .true.))
    allocate(borrowed, source=owner_matches(this % whole_rel, this % ne, part_id, .false., .false.))

    call carve(members, [owned, borrowed], 'overlap_edges', &
         & this % eset, sets, labels, inclusions)

  end subroutine overlap_edges

  !===================================================================!
  ! Collect the indices a part owns, or the ones it does not.
  !
  ! An uncut graph has no ownership record to read. The result: everything
  ! owned, nothing borrowed - correct for a graph that is the whole
  ! of itself.
  !===================================================================!

  !===================================================================!
  ! Optional arguments, defaulted where the relation demands a value.
  ! A piece told its global vertex names but not its whole's identity
  ! answers about itself, which is the honest reading of silence.
  !===================================================================!

  function merge_set(given, fallback) result(s)
    type(graph), intent(in), optional :: given
    type(graph), intent(in)           :: fallback
    type(graph)                       :: s
    if (present(given)) then
       s = given
    else
       s = fallback
    end if
  end function merge_set

  pure integer function merge_count(given, fallback)
    integer, intent(in), optional :: given
    integer, intent(in)           :: fallback
    if (present(given)) then
       merge_count = given
    else
       merge_count = fallback
    end if
  end function merge_count

  pure function pick_global(given, n) result(g)
    integer, intent(in), optional :: given(:)
    integer, intent(in)           :: n
    integer, allocatable          :: g(:)
    integer                       :: i
    if (present(given)) then
       g = given
    else
       g = [(i, i = 1, n)]
    end if
  end function pick_global

  pure function pick_owner(given, n, mine) result(o)
    integer, intent(in), optional :: given(:)
    integer, intent(in)           :: n, mine
    integer, allocatable          :: o(:)
    integer                       :: i
    if (present(given)) then
       o = given
    else
       o = [(mine, i = 1, n)]
    end if
  end function pick_owner

  pure function owner_matches(r, n, part_id, on_vertices, want_owned) result(pick)

    type(partition_relation), intent(in) :: r
    integer                             , intent(in) :: n
    integer                             , intent(in) :: part_id
    logical                             , intent(in) :: on_vertices
    logical                             , intent(in) :: want_owned

    integer, allocatable :: pick(:)
    integer :: i, k, owns

    ! A graph born whole owns everything and borrows nothing. That is
    ! not a special case: it is what the identity relation means.
    if (.not. r % has_part_relation()) then
       if (want_owned) then
          pick = [(i, i = 1, n)]
       else
          allocate(pick(0))
       end if
       return
    end if

    allocate(pick(n))
    k = 0
    do i = 1, n
       owns = r % owner_part(i, on_vertices)
       if ((owns == part_id) .eqv. want_owned) then
          k = k + 1
          pick(k) = i
       end if
    end do
    pick = pick(1:k)

  end function owner_matches

  !===================================================================!
  ! Walking the graph. Each of these is a slice of a list built once
  ! at construction, which is what keeps them pure and cheap enough to
  ! call per vertex.
  !===================================================================!

  pure subroutine incident_edges(this, vertex_index, indices)

    class(stored_directed_graph), intent(in)   :: this
    integer            , intent(in)   :: vertex_index
    integer, allocatable, intent(out) :: indices(:)

    indices = this % einc(this % xinc(vertex_index) : this % xinc(vertex_index + 1) - 1)

  end subroutine incident_edges

  !===================================================================!
  ! The vertices one edge away, either direction.
  !===================================================================!

  pure subroutine adjacent_vertices(this, vertex_index, indices)

    class(stored_directed_graph), intent(in)   :: this
    integer            , intent(in)   :: vertex_index
    integer, allocatable, intent(out) :: indices(:)

    indices = this % vadj(this % xadj(vertex_index) : this % xadj(vertex_index + 1) - 1)

  end subroutine adjacent_vertices

  !===================================================================!
  ! The edges leaving this vertex.
  !===================================================================!

  pure subroutine outgoing_edges(this, vertex_index, indices)

    class(stored_directed_graph), intent(in)   :: this
    integer            , intent(in)   :: vertex_index
    integer, allocatable, intent(out) :: indices(:)

    if (this % reversed) then
       indices = this % ein(this % xin(vertex_index) : this % xin(vertex_index + 1) - 1)
    else
       indices = this % eout(this % xout(vertex_index) : this % xout(vertex_index + 1) - 1)
    end if

  end subroutine outgoing_edges

  !===================================================================!
  ! The edges entering this vertex.
  !===================================================================!

  pure subroutine incoming_edges(this, vertex_index, indices)

    class(stored_directed_graph), intent(in)   :: this
    integer            , intent(in)   :: vertex_index
    integer, allocatable, intent(out) :: indices(:)

    if (this % reversed) then
       indices = this % eout(this % xout(vertex_index) : this % xout(vertex_index + 1) - 1)
    else
       indices = this % ein(this % xin(vertex_index) : this % xin(vertex_index + 1) - 1)
    end if

  end subroutine incoming_edges

  !===================================================================!
  ! Where the outgoing edges land, and where the incoming ones came
  ! from. An edge with no head leads nowhere and is left out.
  !===================================================================!

  pure subroutine outgoing_vertices(this, vertex_index, indices)

    class(stored_directed_graph), intent(in)   :: this
    integer            , intent(in)   :: vertex_index
    integer, allocatable, intent(out) :: indices(:)

    integer, allocatable :: edges(:)
    integer :: k, n

    call this % outgoing_edges(vertex_index, edges)
    allocate(indices(size(edges)))
    n = 0
    do k = 1, size(edges)
       if (this % edge_has_head(edges(k))) then
          n = n + 1
          indices(n) = this % edge_head(edges(k))
       end if
    end do
    indices = indices(1:n)

  end subroutine outgoing_vertices

  !===================================================================!
  ! The vertices whose edges enter this one - the upstream
  ! neighbours.
  !===================================================================!

  pure subroutine incoming_vertices(this, vertex_index, indices)

    class(stored_directed_graph), intent(in)   :: this
    integer            , intent(in)   :: vertex_index
    integer, allocatable, intent(out) :: indices(:)

    integer, allocatable :: edges(:)
    integer :: k

    call this % incoming_edges(vertex_index, edges)
    allocate(indices(size(edges)))
    do k = 1, size(edges)
       indices(k) = this % edge_tail(edges(k))
    end do

  end subroutine incoming_vertices

  !===================================================================!
  ! THE RELATION TO THE WHOLE, HANDED BACK BY VALUE.
  !
  ! Eight questions used to stand here as bindings on the contract:
  ! how many parts, which part owns what, and the maps both ways. Not
  ! one of them is a question about D = (V, E, tail, head). They are
  ! r's - r <= S_part x S_whole - and this graph answers only WHICH
  ! relation it stands in.
  !
  ! By value, and deliberately: the four verbs are HANDED r, so what
  ! they receive must be something that cannot change under them when
  ! the graph it came from goes out of scope.
  !===================================================================!

  type(partition_relation) function whole_relation(this)

    class(stored_directed_graph), intent(in) :: this

    whole_relation = this % whole_rel

  end function whole_relation


  !===================================================================!
  ! The structure read as relations, derived on request from the
  ! stored table over counted coordinates (1..nv, 1..ne), which
  ! keeps every query on the result O(1). A caller holding T and H
  ! may compose, transpose, and query them as relations; the
  ! graph's own answers keep reading the compiled snapshots, and a
  ! graph nobody reads relationally never builds these.
  !===================================================================!

  type(csr_relation) function tail_relation(this)

    class(stored_directed_graph), intent(in) :: this

    type(set_map) :: sets
    integer, allocatable :: table(:,:)
    integer :: k

    call sets % bind(this % vset, counted_set_representation(this % nv))
    call sets % bind(this % eset, counted_set_representation(this % ne))

    allocate(table(2, this % ne))
    do k = 1, this % ne
       table(:, k) = [k, this % edge_tail(k)]
    end do

    tail_relation = csr_relation('edge tails', this % eset, &
         & this % vset, table, sets)

  end function tail_relation

  type(csr_relation) function head_relation(this)

    class(stored_directed_graph), intent(in) :: this

    type(set_map) :: sets
    integer, allocatable :: table(:,:)
    integer :: nh, k

    call sets % bind(this % vset, counted_set_representation(this % nv))
    call sets % bind(this % eset, counted_set_representation(this % ne))

    nh = 0
    do k = 1, this % ne
       if (this % edge_has_head(k)) nh = nh + 1
    end do
    allocate(table(2, nh))
    nh = 0
    do k = 1, this % ne
       if (this % edge_has_head(k)) then
          nh = nh + 1
          table(:, nh) = [k, this % edge_head(k)]
       end if
    end do

    head_relation = csr_relation('edge heads', this % eset, &
         & this % vset, table, sets)

  end function head_relation

end module view_directed_stored
