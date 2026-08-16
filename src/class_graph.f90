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

module class_graph

  use graph_ordinary_view, only : graph
  use fractal_graph      , only : set_graph => graph
  use graph_calculus     , only : GRAPH_SIDE_VERTEX, GRAPH_SIDE_EDGE
  use graph_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use graph_set_map      , only : set_map
  use graph_label_map    , only : label_map
  use graph_inclusion_map, only : inclusion_map

  implicit none

  private
  public :: stored_graph

  !===================================================================!
  ! A graph that keeps its own structure in arrays.
  !===================================================================!

  type, extends(graph) :: stored_graph

     integer :: number = 1
     integer :: nv     = 0
     integer :: ne     = 0

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
     ! How this graph relates to a whole one. A graph straight off a
     ! mesh file has no relation: one part, itself, identity maps. A
     ! partitioner fills these in.
     !----------------------------------------------------------------!

     ! Which part this one is has a name already - the graph's own
     ! number, which the partitioner stamps as it cuts. A second
     ! component saying the same thing is a second thing to keep
     ! true.
     logical :: cut    = .false.
     integer :: nparts = 1

     integer, allocatable :: vowner(:), eowner(:)
     integer, allocatable :: vglobal(:) , eglobal(:)

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

     type(set_graph) :: vset
     type(set_graph) :: eset

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

     !----------------------------------------------------------------!
     ! The named vertex sets.
     !----------------------------------------------------------------!

     procedure :: all_vertices
     procedure :: interior_vertices
     procedure :: boundary_vertices
     procedure :: tagged_vertices

     !----------------------------------------------------------------!
     ! The named edge sets.
     !----------------------------------------------------------------!

     procedure :: all_edges
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
     ! How a part relates to the whole.
     !----------------------------------------------------------------!

     procedure :: num_parts
     procedure :: has_part_relation
     procedure :: global_vertex_index
     procedure :: global_edge_index
     procedure :: part_vertex_index
     procedure :: part_edge_index
     procedure :: vertex_owner_part
     procedure :: edge_owner_part

  end type stored_graph

  !===================================================================!
  ! Constructor.
  !===================================================================!

  interface stored_graph
     module procedure create
  end interface stored_graph

contains

  !===================================================================!
  ! Build a graph from a vertex count and an edge list.
  !
  ! A head of zero (or anything outside 1..nv) means the edge has no
  ! head - a boundary face. Tags are optional; an untagged graph
  ! returns an empty set for every tagged query.
  !===================================================================!

  type(stored_graph) function create(nv, tails, heads, vtags, etags, &
       &                             number, vglobal, vowner, eglobal, &
       &                             eowner, nparts) result(this)

    integer           , intent(in)           :: nv
    integer           , intent(in)           :: tails(:)
    integer           , intent(in)           :: heads(:)
    character(len=*)  , intent(in), optional :: vtags(:)
    character(len=*)  , intent(in), optional :: etags(:)
    integer           , intent(in), optional :: number

    !----------------------------------------------------------------!
    ! The frame, for a graph that is a piece of a larger one: what
    ! each of its own numbers was called in the whole, and which part
    ! owns it. These arrive HERE or not at all - a graph that could
    ! be told its frame afterwards would answer the same question two
    ! ways in one lifetime, which is the one thing the grammar says a
    ! graph may never do. Present vglobal is what makes a graph a
    ! piece; absent, it is a whole.
    !----------------------------------------------------------------!

    integer           , intent(in), optional :: vglobal(:)
    integer           , intent(in), optional :: vowner(:)
    integer           , intent(in), optional :: eglobal(:)
    integer           , intent(in), optional :: eowner(:)
    integer           , intent(in), optional :: nparts

    integer :: e

    this % nv = nv
    this % ne = size(tails)

    ! Declare the two domains once, here, so every later answer
    ! carries one identity per side for this graph's whole life.
    call this % vset % declare()
    call this % eset % declare()

    if (present(number)) this % number = number

    if (present(vglobal)) then
       this % cut = .true.
       allocate(this % vglobal, source=vglobal)
    end if
    if (present(vowner))  allocate(this % vowner , source=vowner)
    if (present(eglobal)) allocate(this % eglobal, source=eglobal)
    if (present(eowner))  allocate(this % eowner , source=eowner)
    if (present(nparts))  this % nparts = nparts

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

    integer :: ne, e, v, n

    ne = size(tail)

    allocate(xptr(nv + 1))
    xptr = 0

    do e = 1, ne
       xptr(tail(e) + 1) = xptr(tail(e) + 1) + 1
       if (head(e) >= 1) xptr(head(e) + 1) = xptr(head(e) + 1) + 1
    end do

    xptr(1) = 1
    do v = 1, nv
       xptr(v + 1) = xptr(v + 1) + xptr(v)
    end do

    n = xptr(nv + 1) - 1
    allocate(elist(max(n, 0)))

    block
      integer, allocatable :: fill(:)
      allocate(fill, source=xptr(1:nv))
      do e = 1, ne
         v = tail(e)
         elist(fill(v)) = e
         fill(v) = fill(v) + 1
         if (head(e) >= 1) then
            v = head(e)
            elist(fill(v)) = e
            fill(v) = fill(v) + 1
         end if
      end do
    end block

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

    integer :: ne, e, v, n

    ne = size(endpoint)

    allocate(xptr(nv + 1))
    xptr = 0

    do e = 1, ne
       if (endpoint(e) >= 1) xptr(endpoint(e) + 1) = xptr(endpoint(e) + 1) + 1
    end do

    xptr(1) = 1
    do v = 1, nv
       xptr(v + 1) = xptr(v + 1) + xptr(v)
    end do

    n = xptr(nv + 1) - 1
    allocate(elist(max(n, 0)))

    block
      integer, allocatable :: fill(:)
      allocate(fill, source=xptr(1:nv))
      do e = 1, ne
         v = endpoint(e)
         if (v >= 1) then
            elist(fill(v)) = e
            fill(v) = fill(v) + 1
         end if
      end do
    end block

  end subroutine build_directed

  !===================================================================!
  ! Identity and size.
  !===================================================================!

  pure integer function id(this)

    class(stored_graph), intent(in) :: this

    id = this % number

  end function id

  !===================================================================!
  ! How many vertices.
  !===================================================================!

  pure integer function num_vertices(this)

    class(stored_graph), intent(in) :: this

    num_vertices = this % nv

  end function num_vertices

  !===================================================================!
  ! The two carriers, as declared at construction. Copies of one
  ! stamped domain: every call answers a set that same_as agrees is
  ! the same set (AGENTS.md, phase 1).
  !===================================================================!

  type(set_graph) function vertex_set(this)

    class(stored_graph), intent(in) :: this

    vertex_set = this % vset

  end function vertex_set

  type(set_graph) function edge_set(this)

    class(stored_graph), intent(in) :: this

    edge_set = this % eset

  end function edge_set

  !===================================================================!
  ! What this graph calls its own two domains, bound into the CALLER'S
  ! label map. The graph knows the names; it does not keep the map, so
  ! a caller that names nothing never calls this and carries nothing.
  !===================================================================!

  subroutine name_carriers(this, labels)

    class(stored_graph), intent(in)    :: this
    type(label_map)    , intent(inout) :: labels

    call labels % bind(this % vset, 'vertices')
    call labels % bind(this % eset, 'edges')

  end subroutine name_carriers

  !===================================================================!
  ! How many edges.
  !===================================================================!

  pure integer function num_edges(this)

    class(stored_graph), intent(in) :: this

    num_edges = this % ne

  end function num_edges

  !===================================================================!
  ! Where an edge goes.
  !===================================================================!

  pure integer function edge_tail(this, edge_index)

    class(stored_graph), intent(in) :: this
    integer            , intent(in) :: edge_index

    edge_tail = this % tail(edge_index)

  end function edge_tail

  !===================================================================!
  ! The vertex an edge enters:  (i) --e--> (j)  answers j. A
  ! boundary edge enters nothing and answers zero.
  !===================================================================!

  pure integer function edge_head(this, edge_index)

    class(stored_graph), intent(in) :: this
    integer            , intent(in) :: edge_index

    edge_head = this % head(edge_index)

  end function edge_head

  !===================================================================!
  ! Whether the edge enters a vertex at all; false marks a boundary
  ! edge.
  !===================================================================!

  pure logical function edge_has_head(this, edge_index)

    class(stored_graph), intent(in) :: this
    integer            , intent(in) :: edge_index

    edge_has_head = this % head(edge_index) >= 1

  end function edge_has_head

  !===================================================================!
  ! The named vertex sets. A boundary vertex is one that touches a
  ! boundary edge; an interior vertex is one that does not.
  !===================================================================!

  type(set_graph) function all_vertices(this) result(members)

    class(stored_graph), intent(in) :: this

    members = this % vset

  end function all_vertices

  !===================================================================!
  ! The vertices that touch no boundary edge.
  !===================================================================!

  subroutine interior_vertices(this, sets, labels, inclusions, members)

    class(stored_graph), intent(in)    :: this
    type(set_map)      , intent(inout) :: sets
    type(label_map)    , intent(inout) :: labels
    type(inclusion_map), intent(inout) :: inclusions
    type(set_graph)    , intent(out)   :: members

    integer, allocatable :: pick(:)
    integer :: v, n

    allocate(pick(this % nv))
    n = 0
    do v = 1, this % nv
       if (.not. touches_boundary(this, v)) then
          n = n + 1
          pick(n) = v
       end if
    end do

    call carve(members, pick(1:n), 'interior_vertices', &
         & this % vset, sets, labels, inclusions)

  end subroutine interior_vertices

  !===================================================================!
  ! The vertices that touch a boundary edge.
  !===================================================================!

  subroutine boundary_vertices(this, sets, labels, inclusions, members)

    class(stored_graph), intent(in)    :: this
    type(set_map)      , intent(inout) :: sets
    type(label_map)    , intent(inout) :: labels
    type(inclusion_map), intent(inout) :: inclusions
    type(set_graph)    , intent(out)   :: members

    integer, allocatable :: pick(:)
    integer :: v, n

    allocate(pick(this % nv))
    n = 0
    do v = 1, this % nv
       if (touches_boundary(this, v)) then
          n = n + 1
          pick(n) = v
       end if
    end do

    call carve(members, pick(1:n), 'boundary_vertices', &
         & this % vset, sets, labels, inclusions)

  end subroutine boundary_vertices

  !===================================================================!
  ! The vertices carrying this tag - a mesh's named patches arrive
  ! here.
  !===================================================================!

  subroutine tagged_vertices(this, tag, sets, labels, inclusions, members)

    class(stored_graph), intent(in)    :: this
    character(len=*)   , intent(in)    :: tag
    type(set_map)      , intent(inout) :: sets
    type(label_map)    , intent(inout) :: labels
    type(inclusion_map), intent(inout) :: inclusions
    type(set_graph)    , intent(out)   :: members

    integer, allocatable :: pick(:)
    integer :: v, n

    allocate(pick(this % nv))
    n = 0
    if (allocated(this % vtag)) then
       do v = 1, this % nv
          if (trim(this % vtag(v)) == tag) then
             n = n + 1
             pick(n) = v
          end if
       end do
    end if

    call carve(members, pick(1:n), 'tagged_vertices', &
         & this % vset, sets, labels, inclusions)

  end subroutine tagged_vertices

  !===================================================================!
  ! CARVE. The one gate every named subset passes through.
  !
  ! A carved set is a NEW set - it signs a fresh identity, exactly as
  ! the subset_set it replaces did - and three things must be said
  ! about it or it is not usable:
  !
  !     its extension     which members, in this order
  !     its label         what the old subset called itself
  !     its embedding     which carrier it was carved from
  !
  ! They are bound together HERE rather than at each of the twelve
  ! call sites, because the third is the one an author forgets: a
  ! missing representation stops the program at the first query, and a
  ! missing label answers '', but a missing inclusion answers FALSE to
  ! is_subobject_of - quietly, and only on a real mesh.
  !
  ! The maps are the caller's. This routine writes into them and keeps
  ! nothing.
  !===================================================================!

  subroutine carve(members, roll, label, ambient, sets, labels, inclusions)

    type(set_graph)    , intent(out)   :: members
    integer            , intent(in)    :: roll(:)
    character(len=*)   , intent(in)    :: label
    type(set_graph)    , intent(in)    :: ambient
    type(set_map)      , intent(inout) :: sets
    type(label_map)    , intent(inout) :: labels
    type(inclusion_map), intent(inout) :: inclusions

    call members % declare()

    call sets       % bind(members, listed_set_representation(roll))
    call labels     % bind(members, label)
    call inclusions % include_in(members, ambient)

  end subroutine carve

  !===================================================================!
  ! Does any edge touching this vertex stop here rather than holding
  ! on to another vertex?
  !===================================================================!

  pure logical function touches_boundary(this, v)

    class(stored_graph), intent(in) :: this
    integer            , intent(in) :: v

    integer :: k

    touches_boundary = .false.
    do k = this % xinc(v), this % xinc(v + 1) - 1
       if (this % head(this % einc(k)) < 1) then
          touches_boundary = .true.
          return
       end if
    end do

  end function touches_boundary

  !===================================================================!
  ! The named edge sets. A boundary edge is one with no head.
  !===================================================================!

  type(set_graph) function all_edges(this) result(members)

    class(stored_graph), intent(in) :: this

    members = this % eset

  end function all_edges

  !===================================================================!
  ! The edges with a head: both ends real.
  !===================================================================!

  subroutine interior_edges(this, sets, labels, inclusions, members)

    class(stored_graph), intent(in)    :: this
    type(set_map)      , intent(inout) :: sets
    type(label_map)    , intent(inout) :: labels
    type(inclusion_map), intent(inout) :: inclusions
    type(set_graph)    , intent(out)   :: members

    integer, allocatable :: pick(:)
    integer :: e, n

    allocate(pick(this % ne))
    n = 0
    do e = 1, this % ne
       if (this % head(e) >= 1) then
          n = n + 1
          pick(n) = e
       end if
    end do

    call carve(members, pick(1:n), 'interior_edges', &
         & this % eset, sets, labels, inclusions)

  end subroutine interior_edges

  !===================================================================!
  ! The edges with no head - the open ends of the graph.
  !===================================================================!

  subroutine boundary_edges(this, sets, labels, inclusions, members)

    class(stored_graph), intent(in)    :: this
    type(set_map)      , intent(inout) :: sets
    type(label_map)    , intent(inout) :: labels
    type(inclusion_map), intent(inout) :: inclusions
    type(set_graph)    , intent(out)   :: members

    integer, allocatable :: pick(:)
    integer :: e, n

    allocate(pick(this % ne))
    n = 0
    do e = 1, this % ne
       if (this % head(e) < 1) then
          n = n + 1
          pick(n) = e
       end if
    end do

    call carve(members, pick(1:n), 'boundary_edges', &
         & this % eset, sets, labels, inclusions)

  end subroutine boundary_edges

  !===================================================================!
  ! The edges carrying this tag - a mesh's named patches arrive
  ! here.
  !===================================================================!

  subroutine tagged_edges(this, tag, sets, labels, inclusions, members)

    class(stored_graph), intent(in)    :: this
    character(len=*)   , intent(in)    :: tag
    type(set_map)      , intent(inout) :: sets
    type(label_map)    , intent(inout) :: labels
    type(inclusion_map), intent(inout) :: inclusions
    type(set_graph)    , intent(out)   :: members

    integer, allocatable :: pick(:)
    integer :: e, n

    allocate(pick(this % ne))
    n = 0
    if (allocated(this % etag)) then
       do e = 1, this % ne
          if (trim(this % etag(e)) == tag) then
             n = n + 1
             pick(n) = e
          end if
       end do
    end if

    call carve(members, pick(1:n), 'tagged_edges', &
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

    class(stored_graph), intent(in)    :: this
    integer            , intent(in)    :: part_id
    type(set_map)      , intent(inout) :: sets
    type(label_map)    , intent(inout) :: labels
    type(inclusion_map), intent(inout) :: inclusions
    type(set_graph)    , intent(out)   :: members

    call carve(members, owner_matches(this % vowner, this % nv, part_id, this % cut, .true.), 'owned_vertices', &
         & this % vset, sets, labels, inclusions)

  end subroutine owned_vertices

  !===================================================================!
  ! The vertices this part reads but does not own - the neighbours'
  ! cells along the cut.
  !===================================================================!

  subroutine borrowed_vertices(this, part_id, sets, labels, inclusions, members)

    class(stored_graph), intent(in)    :: this
    integer            , intent(in)    :: part_id
    type(set_map)      , intent(inout) :: sets
    type(label_map)    , intent(inout) :: labels
    type(inclusion_map), intent(inout) :: inclusions
    type(set_graph)    , intent(out)   :: members

    call carve(members, owner_matches(this % vowner, this % nv, part_id, this % cut, .false.), 'borrowed_vertices', &
         & this % vset, sets, labels, inclusions)

  end subroutine borrowed_vertices

  !===================================================================!
  ! The overlap is everything this part must see to finish what it
  ! owns: what it owns, plus what it borrows.
  !===================================================================!

  subroutine overlap_vertices(this, part_id, sets, labels, inclusions, members)

    class(stored_graph), intent(in)    :: this
    integer            , intent(in)    :: part_id
    type(set_map)      , intent(inout) :: sets
    type(label_map)    , intent(inout) :: labels
    type(inclusion_map), intent(inout) :: inclusions
    type(set_graph)    , intent(out)   :: members

    integer, allocatable :: owned(:), borrowed(:)

    allocate(owned   , source=owner_matches(this % vowner, this % nv, part_id, this % cut, .true.))
    allocate(borrowed, source=owner_matches(this % vowner, this % nv, part_id, this % cut, .false.))

    call carve(members, [owned, borrowed], 'overlap_vertices', &
         & this % vset, sets, labels, inclusions)

  end subroutine overlap_vertices

  !===================================================================!
  ! The edges whose keeper is this part.
  !===================================================================!

  subroutine owned_edges(this, part_id, sets, labels, inclusions, members)

    class(stored_graph), intent(in)    :: this
    integer            , intent(in)    :: part_id
    type(set_map)      , intent(inout) :: sets
    type(label_map)    , intent(inout) :: labels
    type(inclusion_map), intent(inout) :: inclusions
    type(set_graph)    , intent(out)   :: members

    call carve(members, owner_matches(this % eowner, this % ne, part_id, this % cut, .true.), 'owned_edges', &
         & this % eset, sets, labels, inclusions)

  end subroutine owned_edges

  !===================================================================!
  ! The edges this part reads but does not own.
  !===================================================================!

  subroutine borrowed_edges(this, part_id, sets, labels, inclusions, members)

    class(stored_graph), intent(in)    :: this
    integer            , intent(in)    :: part_id
    type(set_map)      , intent(inout) :: sets
    type(label_map)    , intent(inout) :: labels
    type(inclusion_map), intent(inout) :: inclusions
    type(set_graph)    , intent(out)   :: members

    call carve(members, owner_matches(this % eowner, this % ne, part_id, this % cut, .false.), 'borrowed_edges', &
         & this % eset, sets, labels, inclusions)

  end subroutine borrowed_edges

  !===================================================================!
  ! Owned and borrowed together: every edge this part can see.
  !===================================================================!

  subroutine overlap_edges(this, part_id, sets, labels, inclusions, members)

    class(stored_graph), intent(in)    :: this
    integer            , intent(in)    :: part_id
    type(set_map)      , intent(inout) :: sets
    type(label_map)    , intent(inout) :: labels
    type(inclusion_map), intent(inout) :: inclusions
    type(set_graph)    , intent(out)   :: members

    integer, allocatable :: owned(:), borrowed(:)

    allocate(owned   , source=owner_matches(this % eowner, this % ne, part_id, this % cut, .true.))
    allocate(borrowed, source=owner_matches(this % eowner, this % ne, part_id, this % cut, .false.))

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

  pure function owner_matches(owner, n, part_id, cut, want_owned) result(pick)

    integer, allocatable, intent(in) :: owner(:)
    integer             , intent(in) :: n
    integer             , intent(in) :: part_id
    logical             , intent(in) :: cut
    logical             , intent(in) :: want_owned

    integer, allocatable :: pick(:)
    integer :: i, k

    if (.not. cut .or. .not. allocated(owner)) then
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
       if ((owner(i) == part_id) .eqv. want_owned) then
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

    class(stored_graph), intent(in)   :: this
    integer            , intent(in)   :: vertex_index
    integer, allocatable, intent(out) :: indices(:)

    indices = this % einc(this % xinc(vertex_index) : this % xinc(vertex_index + 1) - 1)

  end subroutine incident_edges

  !===================================================================!
  ! The vertices one edge away, either direction.
  !===================================================================!

  pure subroutine adjacent_vertices(this, vertex_index, indices)

    class(stored_graph), intent(in)   :: this
    integer            , intent(in)   :: vertex_index
    integer, allocatable, intent(out) :: indices(:)

    indices = this % vadj(this % xadj(vertex_index) : this % xadj(vertex_index + 1) - 1)

  end subroutine adjacent_vertices

  !===================================================================!
  ! The edges leaving this vertex.
  !===================================================================!

  pure subroutine outgoing_edges(this, vertex_index, indices)

    class(stored_graph), intent(in)   :: this
    integer            , intent(in)   :: vertex_index
    integer, allocatable, intent(out) :: indices(:)

    indices = this % eout(this % xout(vertex_index) : this % xout(vertex_index + 1) - 1)

  end subroutine outgoing_edges

  !===================================================================!
  ! The edges entering this vertex.
  !===================================================================!

  pure subroutine incoming_edges(this, vertex_index, indices)

    class(stored_graph), intent(in)   :: this
    integer            , intent(in)   :: vertex_index
    integer, allocatable, intent(out) :: indices(:)

    indices = this % ein(this % xin(vertex_index) : this % xin(vertex_index + 1) - 1)

  end subroutine incoming_edges

  !===================================================================!
  ! Where the outgoing edges land, and where the incoming ones came
  ! from. An edge with no head leads nowhere and is left out.
  !===================================================================!

  pure subroutine outgoing_vertices(this, vertex_index, indices)

    class(stored_graph), intent(in)   :: this
    integer            , intent(in)   :: vertex_index
    integer, allocatable, intent(out) :: indices(:)

    integer :: k, n, lo, hi

    lo = this % xout(vertex_index)
    hi = this % xout(vertex_index + 1) - 1

    allocate(indices(max(hi - lo + 1, 0)))
    n = 0
    do k = lo, hi
       if (this % head(this % eout(k)) >= 1) then
          n = n + 1
          indices(n) = this % head(this % eout(k))
       end if
    end do
    indices = indices(1:n)

  end subroutine outgoing_vertices

  !===================================================================!
  ! The vertices whose edges enter this one - the upstream
  ! neighbours.
  !===================================================================!

  pure subroutine incoming_vertices(this, vertex_index, indices)

    class(stored_graph), intent(in)   :: this
    integer            , intent(in)   :: vertex_index
    integer, allocatable, intent(out) :: indices(:)

    integer :: k, n, lo, hi

    lo = this % xin(vertex_index)
    hi = this % xin(vertex_index + 1) - 1

    allocate(indices(max(hi - lo + 1, 0)))
    n = 0
    do k = lo, hi
       n = n + 1
       indices(n) = this % tail(this % ein(k))
    end do
    indices = indices(1:n)

  end subroutine incoming_vertices

  !===================================================================!
  ! How a part relates to the whole. An uncut graph is the whole of
  ! itself: one part, identity maps, everything owned by part one.
  !===================================================================!

  pure integer function num_parts(this)

    class(stored_graph), intent(in) :: this

    num_parts = this % nparts

  end function num_parts

  !===================================================================!
  ! Whether this graph knows how it sits in a whole; false for a
  ! graph born whole.
  !===================================================================!

  pure logical function has_part_relation(this)

    class(stored_graph), intent(in) :: this

    has_part_relation = this % cut

  end function has_part_relation

  !===================================================================!
  ! The whole-graph name of this part's vertex. A graph that is not
  ! a part answers the index unchanged.
  !===================================================================!

  pure integer function global_vertex_index(this, index)

    class(stored_graph), intent(in) :: this
    integer            , intent(in) :: index

    if (this % cut .and. allocated(this % vglobal)) then
       global_vertex_index = this % vglobal(index)
    else
       global_vertex_index = index
    end if

  end function global_vertex_index

  !===================================================================!
  ! The same, for an edge.
  !===================================================================!

  pure integer function global_edge_index(this, index)

    class(stored_graph), intent(in) :: this
    integer            , intent(in) :: index

    if (this % cut .and. allocated(this % eglobal)) then
       global_edge_index = this % eglobal(index)
    else
       global_edge_index = index
    end if

  end function global_edge_index

  !===================================================================!
  ! The map read backwards. Zero means the whole-graph index does not
  ! appear in that part at all.
  !===================================================================!

  pure integer function part_vertex_index(this, global_index, part_id)

    class(stored_graph), intent(in) :: this
    integer            , intent(in) :: global_index
    integer            , intent(in) :: part_id

    if (this % cut .and. part_id /= this % number) then
       part_vertex_index = 0
    else
       part_vertex_index = reverse_lookup(this % vglobal, global_index, this % cut)
    end if

  end function part_vertex_index

  !===================================================================!
  ! The part's own index for a whole-graph edge; zero for an edge
  ! this part does not hold.
  !===================================================================!

  pure integer function part_edge_index(this, global_index, part_id)

    class(stored_graph), intent(in) :: this
    integer            , intent(in) :: global_index
    integer            , intent(in) :: part_id

    if (this % cut .and. part_id /= this % number) then
       part_edge_index = 0
    else
       part_edge_index = reverse_lookup(this % eglobal, global_index, this % cut)
    end if

  end function part_edge_index

  !===================================================================!
  ! Find which of the part's own indices holds a given whole-graph
  ! index.
  !===================================================================!

  pure integer function reverse_lookup(global, global_index, cut)

    integer, allocatable, intent(in) :: global(:)
    integer             , intent(in) :: global_index
    logical             , intent(in) :: cut

    integer :: i

    if (.not. cut .or. .not. allocated(global)) then
       reverse_lookup = global_index
       return
    end if

    reverse_lookup = 0
    do i = 1, size(global)
       if (global(i) == global_index) then
          reverse_lookup = i
          return
       end if
    end do

  end function reverse_lookup

  !===================================================================!
  ! Which part keeps this vertex.
  !===================================================================!

  pure integer function vertex_owner_part(this, index)

    class(stored_graph), intent(in) :: this
    integer            , intent(in) :: index

    if (this % cut .and. allocated(this % vowner)) then
       vertex_owner_part = this % vowner(index)
    else
       vertex_owner_part = 1
    end if

  end function vertex_owner_part

  !===================================================================!
  ! Which part keeps this edge.
  !===================================================================!

  pure integer function edge_owner_part(this, index)

    class(stored_graph), intent(in) :: this
    integer            , intent(in) :: index

    if (this % cut .and. allocated(this % eowner)) then
       edge_owner_part = this % eowner(index)
    else
       edge_owner_part = 1
    end if

  end function edge_owner_part


end module class_graph
