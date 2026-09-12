!=====================================================================!
! A concrete graph that stores its structure.
!
! Given a vertex count and an edge list, the constructor computes,
! once, the neighbour lists every later query reads:
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
! An edge whose head is not a vertex of the graph has no head. That is
! a boundary face: it is attached to one cell alone, and no fictitious
! cell is created beyond the boundary.
!
! Four compressed lists are built at construction and never rebuilt -
! the edges incident to a vertex, the vertices adjacent to it, and the
! same two split by edge direction. Every traversal query is then a
! slice of an array, which is what lets those queries stay pure and of
! low enough cost to be called inside a loop over a million cells.
!
!=====================================================================!
!
!                       WHAT THIS GRAPH DOES NOT CONTAIN
!
! The graph stores no geometry, no physics, no solver state, and no
! algorithm. Colouring, traversal order, partitioning and the rest are
! operations and transforms that read a graph; they are not procedures
! of the graph. Each convenience procedure added here weakens that
! separation, one procedure at a time.
!
! THE GRAPH STORES NO VALUES. A field references its domain; the
! reference never points the other way. What an operation reads is
! passed to the operation at construction, as explicit arguments. The
! one string stored is the tag, which is data: the mesh file named its
! boundary groups, and those names are read from outside the code.
!
! READ ONLY. No procedure puts data on a graph after construction.
! Anything computed leaves as an operation's output. Without this rule
! the graph accumulates state, and its results come to depend on the
! order of past calls rather than on the mesh it was built from.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module view_directed_stored

  use view_directed, only : directed_graph, forward, reverse
  use relation_algorithms, only : topological_order
  use graph_fractal      , only : graph
  use token_identity, only : token
  use relation_binary, only : group_by_key, csr_relation
  use relation_partition, only : partition_relation
  use map_set_representation, only : counted_set_representation, listed_set_representation
  use map_set      , only : set_map
  use map_set_store, only : set_store

  implicit none

  private

  public :: stored_directed_graph

  !===================================================================!
  ! A graph that stores its own structure in arrays.
  !===================================================================!

  type, extends(directed_graph) :: stored_directed_graph

     private

     integer :: number = 1
     integer :: nv     = 0
     integer :: ne     = 0

     !----------------------------------------------------------------!
     ! THE ORIENTATION. D = (V, E, tail, head) and its transpose
     ! D^T = (V, E, head, tail) have the same carrier identities:
     ! both endpoint lists and both compressed directions are stored,
     ! and reversed records which is which. The transpose of the
     ! transpose restores the original orientation. Each graph value
     ! owns its arrays and every copy of the value copies them, so a
     ! value is valid under every Fortran copy mechanism. reverse
     ! changes the orientation of this value in place at constant
     ! cost; transpose is an owning copy in the other orientation.
     !----------------------------------------------------------------!
     logical :: reversed = .false.

     !----------------------------------------------------------------!
     ! Edge endpoints. A head of zero means the edge has none; the
     ! number of such edges, counted once at construction, decides
     ! whether the orientation can be reversed.
     !----------------------------------------------------------------!

     integer, allocatable :: tail(:)
     integer, allocatable :: head(:)
     integer :: num_without_head = 0

     !----------------------------------------------------------------!
     ! The compressed lists, built once. Vertex v's incident edges are
     ! einc(xinc(v) : xinc(v+1)-1), and so on for the other three.
     !----------------------------------------------------------------!

     integer, allocatable :: xinc(:), einc(:)
     integer, allocatable :: xadj(:), vadj(:)
     integer, allocatable :: xout(:), eout(:)
     integer, allocatable :: xin(:) , ein(:)

     !----------------------------------------------------------------!
     ! Names stored on vertices and edges. Blank means untagged.
     !----------------------------------------------------------------!

     character(len=:), allocatable :: vtag(:)
     character(len=:), allocatable :: etag(:)

     !----------------------------------------------------------------!
     ! HOW THIS GRAPH RELATES TO A WHOLE ONE, and it is one component
     ! now rather than six: the relation r <= S_part x S_whole. The
     ! record is a REPRESENTATION - part-local numbering, ownership,
     ! provenance - and none of it is a predicate on
     ! D = (V, E, tail, head), which is why it stopped being a binding
     ! on the contract and became a value the graph stores.
     !
     ! A graph read directly from a mesh file stores the IDENTITY relation:
     ! one part, itself, every member its own name.
     !----------------------------------------------------------------!

     type(partition_relation) :: whole_rel

     !----------------------------------------------------------------!
     ! The graph's two carriers (AGENTS.md, phase 1): its vertices
     ! and its edges as declared set GRAPHS, declared once at
     ! construction and returned beside the old vocabulary. Every
     ! call to vertex_set returns the SAME domain, so a relation
     ! signature can retain the identity.
     !
     ! Identity only. The extension is 1..nv and 1..ne, which the
     ! graph already returns through num_vertices/num_edges, so
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

     !----------------------------------------------------------------!
     ! The endpoints of an edge.
     !----------------------------------------------------------------!

     procedure :: edge_tail
     procedure :: edge_head
     procedure :: edge_has_head
     procedure :: transpose
     procedure :: reverse => reverse_orientation
     procedure :: transposed
     procedure :: loop

     !----------------------------------------------------------------!
     ! The named edge set.
     !----------------------------------------------------------------!

     procedure :: tagged_edges

     !----------------------------------------------------------------!
     ! Traversal, without regard to direction and with it.
     !----------------------------------------------------------------!

     procedure :: incident_edges
     procedure :: adjacent_vertices
     procedure :: outgoing_edges
     procedure :: incoming_edges
     procedure :: outgoing_vertices
     procedure :: incoming_vertices
     procedure :: read_incoming

     !----------------------------------------------------------------!
     ! How a part relates to the whole: ONE accessor, returning the
     ! relation by value. The eight queries that used to stand here
     ! are r's, and a caller that needs them takes r and queries it -
     ! which is also what lets the four verbs be passed one explicitly.
     !
     ! whole_relation, not relation: this graph CONTAINS relations -
     ! its incidence and its adjacency are two of them - and this is
     ! not one of those. It is the relation to the whole.
     !----------------------------------------------------------------!

     procedure :: whole_relation

  end type stored_directed_graph

  !===================================================================!
  ! Constructor.
  !===================================================================!

  interface stored_directed_graph
     module procedure create
  end interface stored_directed_graph

  abstract interface
     subroutine incoming_reader(offsets, indices, sources)
       integer, intent(in) :: offsets(:), indices(:), sources(:)
     end subroutine incoming_reader
  end interface

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
    ! part owns it. These are passed HERE or not at all - a graph whose
    ! relation could be set afterwards would return two results for
    ! the same query in one lifetime, which the grammar forbids.
    ! Present vglobal is what makes a graph a piece; absent, it is a
    ! whole.
    !----------------------------------------------------------------!

    integer           , intent(in), optional :: vglobal(:)
    integer           , intent(in), optional :: vowner(:)
    integer           , intent(in), optional :: eglobal(:)
    integer           , intent(in), optional :: eowner(:)
    integer           , intent(in), optional :: num_parts

    !----------------------------------------------------------------!
    ! The far side of r, by identity and count. A relation that did
    ! not record WHICH sets it relates would let a caller that stores
    ! two of them pass the wrong one to the wrong graph without any
    ! error.
    !----------------------------------------------------------------!

    type(graph)   , intent(in), optional :: whole_vertices, whole_edges
    integer           , intent(in), optional :: num_whole_vertices, num_whole_edges

    integer :: e, np, part_number, nwv, nwe
    type(token) :: whole_identity

    ! Validate the input arrays before reading endpoints or constructing
    ! any compressed index. Every edge has a tail in the vertex set;
    ! the documented boundary convention applies only to its head.
    if (nv < 0) error stop 'stored_directed_graph: the vertex count is nonnegative'
    if (size(heads) /= size(tails)) then
       error stop 'stored_directed_graph: head and tail arrays have equal extent'
    end if
    if (any(tails < 1) .or. any(tails > nv)) then
       error stop 'stored_directed_graph: every tail belongs to the vertex set'
    end if
    if (present(vtags)) then
       if (size(vtags) /= nv) error stop 'stored_directed_graph: one tag per vertex'
    end if
    if (present(etags)) then
       if (size(etags) /= size(tails)) error stop 'stored_directed_graph: one tag per edge'
    end if

    if (present(vglobal)) then
       if (size(vglobal) /= nv) error stop 'stored_directed_graph: one global index per vertex'
       np = merge_count(num_parts, 1)
       part_number = merge_count(number, 1)
       nwv = merge_count(num_whole_vertices, nv)
       nwe = merge_count(num_whole_edges, size(tails))
       if (present(whole_vertices)) then
          whole_identity = whole_vertices % id()
          if (.not. whole_identity % declared()) then
             error stop 'stored_directed_graph: the whole vertex carrier is declared'
          end if
       else if (nwv /= nv) then
          error stop 'stored_directed_graph: one vertex carrier has one extent'
       end if
       if (present(whole_edges)) then
          whole_identity = whole_edges % id()
          if (.not. whole_identity % declared()) then
             error stop 'stored_directed_graph: the whole edge carrier is declared'
          end if
       else if (nwe /= size(tails)) then
          error stop 'stored_directed_graph: one edge carrier has one extent'
       end if
       if (np < 1 .or. part_number < 1 .or. part_number > np) then
          error stop 'stored_directed_graph: the part number belongs to the set of parts'
       end if
       if (nwv < 0 .or. nwe < 0) then
          error stop 'stored_directed_graph: whole carrier counts are nonnegative'
       end if
       if (any(vglobal < 1) .or. any(vglobal > nwv)) then
          error stop 'stored_directed_graph: global vertices belong to the whole carrier'
       end if
       if (.not. injective(vglobal)) then
          error stop 'stored_directed_graph: global vertex indices are injective'
       end if
       if (present(vowner)) then
          if (size(vowner) /= nv) error stop 'stored_directed_graph: one owner per vertex'
          if (any(vowner < 1) .or. any(vowner > np)) then
             error stop 'stored_directed_graph: vertex owners belong to the set of parts'
          end if
       end if
       if (present(eglobal)) then
          if (size(eglobal) /= size(tails)) error stop 'stored_directed_graph: one global index per edge'
          if (any(eglobal < 1) .or. any(eglobal > nwe)) then
             error stop 'stored_directed_graph: global edges belong to the whole carrier'
          end if
          if (.not. injective(eglobal)) then
             error stop 'stored_directed_graph: global edge indices are injective'
          end if
       else if (size(tails) > nwe) then
          error stop 'stored_directed_graph: global edges belong to the whole carrier'
       end if
       if (present(eowner)) then
          if (size(eowner) /= size(tails)) error stop 'stored_directed_graph: one owner per edge'
          if (any(eowner < 1) .or. any(eowner > np)) then
             error stop 'stored_directed_graph: edge owners belong to the set of parts'
          end if
       end if
    else if (present(vowner) .or. present(eglobal) .or. present(eowner) .or. &
         & present(num_parts) .or. present(whole_vertices) .or. present(whole_edges) .or. &
         & present(num_whole_vertices) .or. present(num_whole_edges)) then
       error stop 'stored_directed_graph: a partition relation requires global vertex indices'
    end if

    this % nv = nv
    this % ne = size(tails)

    ! Declare the two domains once, here, so every later result
    ! has one identity per side for this graph's whole lifetime.
    call this % vset % declare()
    call this % eset % declare()

    if (present(number)) this % number = number

    allocate(this % tail, source=tails)
    allocate(this % head(this % ne))

    ! Normalise every missing head to zero, so one test suffices
    ! everywhere afterwards.
    this % num_without_head = 0
    do e = 1, this % ne
       if (heads(e) >= 1 .and. heads(e) <= nv) then
          this % head(e) = heads(e)
       else
          this % head(e) = 0
          this % num_without_head = this % num_without_head + 1
       end if
    end do

    !----------------------------------------------------------------!
    ! The relation is stored here at construction, or the identity
    ! relation is. A graph whose relation could be set afterwards would
    ! return two results for one query in one lifetime; present vglobal
    ! is what makes a graph a piece, and absent, it is a whole.
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
            & vowner  = select_owner(vowner, size(vglobal), this % number), &
            & eglobal = select_global(eglobal, this % ne),             &
            & eowner  = select_owner(eowner, this % ne, this % number))
    else
       this % whole_rel = partition_relation( &
            & this % vset, this % nv, this % eset, this % ne)
    end if

    if (present(vtags)) allocate(this % vtag, source=vtags)
    if (present(etags)) allocate(this % etag, source=etags)

    ! Everything the mesh recorded is stored here and never changes again.

    call build_incidence(this % nv, this % tail, this % head, &
         &               this % xinc, this % einc)

    call build_adjacency(this % nv, this % tail, this % head, &
         &               this % xinc, this % einc, this % xadj, this % vadj)

    call build_directed(this % nv, this % tail, this % xout, this % eout)
    call build_directed(this % nv, this % head, this % xin , this % ein )

  end function create

  ! Injectivity is a cardinality equality. The listed representation
  ! removes repetitions through its O(n log n) inverse construction,
  ! using storage proportional to this part, not the whole carrier.
  pure logical function injective(values)

    integer, intent(in) :: values(:)
    type(listed_set_representation) :: members

    members = listed_set_representation(values)
    injective = members % num_members() == size(values)

  end function injective

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
    ! so each vertex's fibre retains the single-pass edge order; a
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
  ! each counted once however many edges join the pair. The incident
  ! list is already grouped by vertex, so one pass writes the distinct
  ! far ends; a marker records the vertex each far end was last
  ! written for.
  !===================================================================!

  pure subroutine build_adjacency(nv, tail, head, xinc, einc, xptr, vlist)

    integer             , intent(in)  :: nv
    integer             , intent(in)  :: tail(:), head(:)
    integer             , intent(in)  :: xinc(:), einc(:)
    integer, allocatable, intent(out) :: xptr(:), vlist(:)

    integer, allocatable :: last_marked(:)
    integer :: v, k, e, other, total

    allocate(xptr(nv + 1), vlist(size(einc)), last_marked(nv))
    last_marked = 0

    total   = 0
    xptr(1) = 1
    do v = 1, nv
       do k = xinc(v), xinc(v + 1) - 1
          e = einc(k)
          other = opposite_endpoint(tail(e), head(e), v)
          if (other >= 1 .and. other /= v) then
             if (last_marked(other) /= v) then
                last_marked(other) = v
                total = total + 1
                vlist(total) = other
             end if
          end if
       end do
       xptr(v + 1) = total + 1
    end do
    vlist = vlist(1:total)

  end subroutine build_adjacency

  !===================================================================!
  ! Given both ends of an edge and one of them, return the other.
  ! Returns zero when the edge has no head, which is how a boundary
  ! face records that there is no cell beyond it.
  !===================================================================!

  pure integer function opposite_endpoint(tail, head, endpoint)

    integer, intent(in) :: tail, head, endpoint

    if (tail == endpoint) then
       opposite_endpoint = head
    else
       opposite_endpoint = tail
    end if

  end function opposite_endpoint

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
  ! declared domain: every call returns a set that same_as reports is
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
  ! How many edges.
  !===================================================================!

  pure integer function num_edges(this)

    class(stored_directed_graph), intent(in) :: this

    num_edges = this % ne

  end function num_edges

  !===================================================================!
  ! The endpoints of an edge.
  !===================================================================!
  ! THE LOOP over the graph: its vertices in an order every edge
  ! respects - a tail before its head in the forward orientation, a
  ! head before its tail in reverse, which is the loop over the
  ! transpose. It exists only where the graph has no cycle; a graph
  ! with one has no loop, and the request stops the program. A march
  ! is this loop forward and its adjoint this loop in reverse, and
  ! neither records an instant's index. The order is the one
  ! topological sort in the source tree, over the graph's own adjacency.
  !===================================================================!

  function loop(this, orientation) result(order)

    class(stored_directed_graph), intent(in) :: this
    integer, intent(in), optional :: orientation
    integer, allocatable :: order(:)

    type(set_map)      :: sets
    type(csr_relation) :: adjacency
    integer, allocatable :: table(:,:)
    logical :: acyclic
    integer :: direction, e, n

    direction = forward
    if (present(orientation)) direction = orientation
    if (direction /= forward .and. direction /= reverse) then
       error stop 'stored_directed_graph: a loop runs forward or in reverse'
    end if

    allocate(table(2, this % ne))
    n = 0
    do e = 1, this % ne
       if (.not. this % edge_has_head(e)) cycle
       n = n + 1
       if (direction == forward) then
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
  ! The transpose: an owning copy read in the reverse orientation,
  ! every edge's tail its head and head its tail, so that transposing
  ! twice returns the original. An edge without a head would become an
  ! edge without a tail, which is not an edge; such a graph has no
  ! transpose and the request stops the program. Intrinsic assignment
  ! copies the allocatable storage: the cost is proportional to the
  ! graph, and the copy is valid after this value is finalized or
  ! replaced. A value that is only read in the other orientation from
  ! now on is reversed in place instead.
  !===================================================================!

  type(stored_directed_graph) function transpose(this) result(transposed_graph)

    class(stored_directed_graph), intent(in) :: this

    if (this % num_without_head > 0) then
       error stop 'stored_directed_graph: a graph with an edge without a head has no transpose'
    end if

    transposed_graph = this
    transposed_graph % reversed = .not. this % reversed

  end function transpose

  !===================================================================!
  ! Reverse the orientation of this value in place: the same
  ! immutable arrays read with tail and head exchanged, at constant
  ! cost and without allocation. Reversing twice restores the
  ! orientation. The refusal is the transpose's.
  !===================================================================!

  subroutine reverse_orientation(this)

    class(stored_directed_graph), intent(inout) :: this

    if (this % num_without_head > 0) then
       error stop 'stored_directed_graph: a graph with an edge without a head has no transpose'
    end if

    this % reversed = .not. this % reversed

  end subroutine reverse_orientation

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
  ! The vertex an edge enters:  (i) --e--> (j)  returns j. A
  ! boundary edge enters no vertex and returns zero.
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
  ! The edges with this tag - a mesh's named patches are queried
  ! here. An untagged graph returns an empty set.
  !===================================================================!

  subroutine tagged_edges(this, tag, sets, members)

    class(stored_directed_graph), intent(in)    :: this
    character(len=*)   , intent(in)    :: tag
    type(set_store)    , intent(inout) :: sets
    type(graph)    , intent(out)   :: members

    logical, allocatable :: tagged(:)
    integer :: e

    allocate(tagged(this % ne), source=.false.)
    if (allocated(this % etag)) tagged = [(trim(this % etag(e)) == tag, e = 1, this % ne)]

    call sets % declare_subobject(members, pack([(e, e = 1, this % ne)], tagged), &
         & 'tagged_edges', this % eset)

  end subroutine tagged_edges

  !===================================================================!
  ! Optional arguments, defaulted where the relation requires a value.
  ! A piece given its global vertex names but not its whole's identity
  ! defaults the whole to itself.
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

  pure function select_global(given, n) result(g)
    integer, intent(in), optional :: given(:)
    integer, intent(in)           :: n
    integer, allocatable          :: g(:)
    integer                       :: i
    if (present(given)) then
       g = given
    else
       g = [(i, i = 1, n)]
    end if
  end function select_global

  pure function select_owner(given, n, own_part) result(o)
    integer, intent(in), optional :: given(:)
    integer, intent(in)           :: n, own_part
    integer, allocatable          :: o(:)
    integer                       :: i
    if (present(given)) then
       o = given
    else
       o = [(own_part, i = 1, n)]
    end if
  end function select_owner

  !===================================================================!
  ! Traversing the graph. Each of these is a slice of a list built
  ! once at construction, which is what keeps them pure and of low
  ! enough cost to call per vertex.
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
  ! Read the incoming compressed incidence once for a whole traversal.
  ! For vertex v, indices(offsets(v):offsets(v+1)-1) lists its edges,
  ! and sources(e) is edge e's tail in the current orientation.
  !
  ! The reader receives INTENT(IN) arrays, without TARGET or POINTER.
  ! It can neither alter the topology nor retain a conforming pointer
  ! alias. Storage remains owned by this graph throughout the call;
  ! no neighbourhood allocation or per-edge dispatch is introduced.
  !===================================================================!

  subroutine read_incoming(this, reader)

    class(stored_directed_graph), intent(in) :: this
    procedure(incoming_reader) :: reader

    ! A default value has empty carriers and no allocated indices.
    ! Its incidence has the same empty read as a constructed empty graph.
    if (.not. allocated(this % xin)) then
       call reader([1], [integer ::], [integer ::])
    else if (this % reversed) then
       call reader(this % xout, this % eout, this % head)
    else
       call reader(this % xin, this % ein, this % tail)
    end if

  end subroutine read_incoming

  !===================================================================!
  ! The heads of the outgoing edges, and the tails of the incoming
  ! ones. An edge with no head has no head vertex and is omitted.
  !===================================================================!

  pure subroutine outgoing_vertices(this, vertex_index, indices)

    class(stored_directed_graph), intent(in)   :: this
    integer            , intent(in)   :: vertex_index
    integer, allocatable, intent(out) :: indices(:)

    integer, allocatable :: edges(:)
    integer :: k

    call this % outgoing_edges(vertex_index, edges)
    indices = [(this % edge_head(edges(k)), k = 1, size(edges))]
    indices = pack(indices, indices >= 1)

  end subroutine outgoing_vertices

  !===================================================================!
  ! The vertices whose edges enter this one - the in-neighbours.
  !===================================================================!

  pure subroutine incoming_vertices(this, vertex_index, indices)

    class(stored_directed_graph), intent(in)   :: this
    integer            , intent(in)   :: vertex_index
    integer, allocatable, intent(out) :: indices(:)

    integer, allocatable :: edges(:)
    integer :: k

    call this % incoming_edges(vertex_index, edges)
    indices = [(this % edge_tail(edges(k)), k = 1, size(edges))]

  end subroutine incoming_vertices

  !===================================================================!
  ! THE RELATION TO THE WHOLE, RETURNED BY VALUE.
  !
  ! Eight queries used to stand here as bindings on the contract:
  ! how many parts, which part owns what, and the maps both ways. Not
  ! one of them is a predicate on D = (V, E, tail, head). They are
  ! r's - r <= S_part x S_whole - and this graph returns only WHICH
  ! relation it is in.
  !
  ! By value, and deliberately: the four verbs are PASSED r, so what
  ! they receive must be something that cannot change under them when
  ! the graph it came from goes out of scope.
  !===================================================================!

  type(partition_relation) function whole_relation(this)

    class(stored_directed_graph), intent(in) :: this

    whole_relation = this % whole_rel

  end function whole_relation

end module view_directed_stored
