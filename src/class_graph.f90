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
! a boundary face: it hangs off one cell and stops, and no imaginary
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
! FETCH GRAPH DATA ONCE. get_data hands back a copy, because the
! contract says the answer is allocatable and freshly made. On a
! million cells that copy is megabytes. An operation fetches what it
! needs at the top of apply, before it starts looping - never inside a
! loop over faces.
!
! READ ONLY. No procedure puts data on a graph after construction.
! Anything computed leaves as an operation's output. Without this rule
! the graph accumulates state, and its answers come to depend on the
! order of past calls rather than on the mesh it was built from.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_graph

  use abstract_graph_types , only : graph, graph_data
  use abstract_graph_types , only : graph_vertex_support, graph_edge_support
  use class_graph_support  , only : vertex_support, edge_support
  use class_graph_field    , only : vertex_field, edge_field

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

     logical :: cut    = .false.
     integer :: nparts = 1
     integer :: me     = 1

     integer, allocatable :: vowner(:), eowner(:)
     integer, allocatable :: vfull(:) , efull(:)

     !----------------------------------------------------------------!
     ! Whatever the mesh knew when it was built. One array per kind,
     ! because a Fortran array holds one type and there is exactly one
     ! concrete vertex field and one concrete edge field.
     !----------------------------------------------------------------!

     type(vertex_field), allocatable :: vdata(:)
     type(edge_field)  , allocatable :: edata(:)

   contains

     !----------------------------------------------------------------!
     ! Who am I and how big am I.
     !----------------------------------------------------------------!

     procedure :: id           => g_id
     procedure :: num_vertices => g_num_vertices
     procedure :: num_edges    => g_num_edges

     !----------------------------------------------------------------!
     ! Where an edge goes.
     !----------------------------------------------------------------!

     procedure :: edge_tail     => g_edge_tail
     procedure :: edge_head     => g_edge_head
     procedure :: edge_has_head => g_edge_has_head

     !----------------------------------------------------------------!
     ! The named vertex sets.
     !----------------------------------------------------------------!

     procedure :: all_vertices      => g_all_vertices
     procedure :: interior_vertices => g_interior_vertices
     procedure :: boundary_vertices => g_boundary_vertices
     procedure :: tagged_vertices   => g_tagged_vertices

     !----------------------------------------------------------------!
     ! The named edge sets.
     !----------------------------------------------------------------!

     procedure :: all_edges      => g_all_edges
     procedure :: interior_edges => g_interior_edges
     procedure :: boundary_edges => g_boundary_edges
     procedure :: tagged_edges   => g_tagged_edges

     !----------------------------------------------------------------!
     ! The named sets of one part.
     !----------------------------------------------------------------!

     procedure :: owned_vertices    => g_owned_vertices
     procedure :: borrowed_vertices => g_borrowed_vertices
     procedure :: overlap_vertices  => g_overlap_vertices
     procedure :: owned_edges       => g_owned_edges
     procedure :: borrowed_edges    => g_borrowed_edges
     procedure :: overlap_edges     => g_overlap_edges

     !----------------------------------------------------------------!
     ! Walking, without regard to direction and with it.
     !----------------------------------------------------------------!

     procedure :: incident_edges    => g_incident_edges
     procedure :: adjacent_vertices => g_adjacent_vertices
     procedure :: outgoing_edges    => g_outgoing_edges
     procedure :: incoming_edges    => g_incoming_edges
     procedure :: outgoing_vertices => g_outgoing_vertices
     procedure :: incoming_vertices => g_incoming_vertices

     !----------------------------------------------------------------!
     ! How a part relates to the whole.
     !----------------------------------------------------------------!

     procedure :: num_parts         => g_num_parts
     procedure :: has_part_relation => g_has_part_relation
     procedure :: full_vertex_index    => g_full_vertex_index
     procedure :: full_edge_index      => g_full_edge_index
     procedure :: part_vertex_index    => g_part_vertex_index
     procedure :: part_edge_index      => g_part_edge_index
     procedure :: vertex_owner_part => g_vertex_owner_part
     procedure :: edge_owner_part   => g_edge_owner_part

     !----------------------------------------------------------------!
     ! The data the graph came with.
     !----------------------------------------------------------------!

     procedure :: has_data => g_has_data
     procedure :: get_data => g_get_data

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
  ! simply answers no to every tagged query.
  !===================================================================!

  type(stored_graph) function create(nv, tails, heads, vtags, etags, &
       &                             vdata, edata, number) result(this)

    integer           , intent(in)           :: nv
    integer           , intent(in)           :: tails(:)
    integer           , intent(in)           :: heads(:)
    character(len=*)  , intent(in), optional :: vtags(:)
    character(len=*)  , intent(in), optional :: etags(:)
    type(vertex_field), intent(in), optional :: vdata(:)
    type(edge_field)  , intent(in), optional :: edata(:)
    integer           , intent(in), optional :: number

    integer :: e

    this % nv = nv
    this % ne = size(tails)

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

    if (present(vtags)) allocate(this % vtag, source=vtags)
    if (present(etags)) allocate(this % etag, source=etags)

    ! Everything the mesh knew arrives here and never changes again.
    if (present(vdata)) allocate(this % vdata, source=vdata)
    if (present(edata)) allocate(this % edata, source=edata)

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

  pure integer function g_id(this)

    class(stored_graph), intent(in) :: this

    g_id = this % number

  end function g_id

  pure integer function g_num_vertices(this)

    class(stored_graph), intent(in) :: this

    g_num_vertices = this % nv

  end function g_num_vertices

  pure integer function g_num_edges(this)

    class(stored_graph), intent(in) :: this

    g_num_edges = this % ne

  end function g_num_edges

  !===================================================================!
  ! Where an edge goes.
  !===================================================================!

  pure integer function g_edge_tail(this, edge_index)

    class(stored_graph), intent(in) :: this
    integer            , intent(in) :: edge_index

    g_edge_tail = this % tail(edge_index)

  end function g_edge_tail

  pure integer function g_edge_head(this, edge_index)

    class(stored_graph), intent(in) :: this
    integer            , intent(in) :: edge_index

    g_edge_head = this % head(edge_index)

  end function g_edge_head

  pure logical function g_edge_has_head(this, edge_index)

    class(stored_graph), intent(in) :: this
    integer            , intent(in) :: edge_index

    g_edge_has_head = this % head(edge_index) >= 1

  end function g_edge_has_head

  !===================================================================!
  ! The named vertex sets. A boundary vertex is one that touches a
  ! boundary edge; an interior vertex is one that does not.
  !===================================================================!

  subroutine g_all_vertices(this, support)

    class(stored_graph), intent(in)                       :: this
    class(graph_vertex_support), allocatable, intent(out) :: support

    integer :: v

    allocate(support, source=vertex_support([(v, v = 1, this % nv)]))

  end subroutine g_all_vertices

  subroutine g_interior_vertices(this, support)

    class(stored_graph), intent(in)                       :: this
    class(graph_vertex_support), allocatable, intent(out) :: support

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

    allocate(support, source=vertex_support(pick(1:n)))

  end subroutine g_interior_vertices

  subroutine g_boundary_vertices(this, support)

    class(stored_graph), intent(in)                       :: this
    class(graph_vertex_support), allocatable, intent(out) :: support

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

    allocate(support, source=vertex_support(pick(1:n)))

  end subroutine g_boundary_vertices

  subroutine g_tagged_vertices(this, tag, support)

    class(stored_graph), intent(in)                       :: this
    character(len=*), intent(in)                          :: tag
    class(graph_vertex_support), allocatable, intent(out) :: support

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

    allocate(support, source=vertex_support(pick(1:n)))

  end subroutine g_tagged_vertices

  !===================================================================!
  ! Does any edge touching this vertex stop here rather than carrying
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

  subroutine g_all_edges(this, support)

    class(stored_graph), intent(in)                     :: this
    class(graph_edge_support), allocatable, intent(out) :: support

    integer :: e

    allocate(support, source=edge_support([(e, e = 1, this % ne)]))

  end subroutine g_all_edges

  subroutine g_interior_edges(this, support)

    class(stored_graph), intent(in)                     :: this
    class(graph_edge_support), allocatable, intent(out) :: support

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

    allocate(support, source=edge_support(pick(1:n)))

  end subroutine g_interior_edges

  subroutine g_boundary_edges(this, support)

    class(stored_graph), intent(in)                     :: this
    class(graph_edge_support), allocatable, intent(out) :: support

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

    allocate(support, source=edge_support(pick(1:n)))

  end subroutine g_boundary_edges

  subroutine g_tagged_edges(this, tag, support)

    class(stored_graph), intent(in)                     :: this
    character(len=*), intent(in)                        :: tag
    class(graph_edge_support), allocatable, intent(out) :: support

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

    allocate(support, source=edge_support(pick(1:n)))

  end subroutine g_tagged_edges

  !===================================================================!
  ! The named sets of one part.
  !
  ! A graph that was never cut owns everything and borrows nothing,
  ! whichever part is asked about. A partitioner fills in the owner
  ! arrays and these answers become real.
  !===================================================================!

  subroutine g_owned_vertices(this, part_id, support)

    class(stored_graph), intent(in)                       :: this
    integer, intent(in)                                   :: part_id
    class(graph_vertex_support), allocatable, intent(out) :: support

    allocate(support, source=vertex_support(owner_matches(this % vowner, this % nv, part_id, this % cut, .true.)))

  end subroutine g_owned_vertices

  subroutine g_borrowed_vertices(this, part_id, support)

    class(stored_graph), intent(in)                       :: this
    integer, intent(in)                                   :: part_id
    class(graph_vertex_support), allocatable, intent(out) :: support

    allocate(support, source=vertex_support(owner_matches(this % vowner, this % nv, part_id, this % cut, .false.)))

  end subroutine g_borrowed_vertices

  !===================================================================!
  ! The overlap is everything this part must see to finish what it
  ! owns: what it owns, plus what it borrows.
  !===================================================================!

  subroutine g_overlap_vertices(this, part_id, support)

    class(stored_graph), intent(in)                       :: this
    integer, intent(in)                                   :: part_id
    class(graph_vertex_support), allocatable, intent(out) :: support

    integer, allocatable :: owned(:), borrowed(:)

    allocate(owned   , source=owner_matches(this % vowner, this % nv, part_id, this % cut, .true.))
    allocate(borrowed, source=owner_matches(this % vowner, this % nv, part_id, this % cut, .false.))

    allocate(support, source=vertex_support([owned, borrowed]))

  end subroutine g_overlap_vertices

  subroutine g_owned_edges(this, part_id, support)

    class(stored_graph), intent(in)                     :: this
    integer, intent(in)                                 :: part_id
    class(graph_edge_support), allocatable, intent(out) :: support

    allocate(support, source=edge_support(owner_matches(this % eowner, this % ne, part_id, this % cut, .true.)))

  end subroutine g_owned_edges

  subroutine g_borrowed_edges(this, part_id, support)

    class(stored_graph), intent(in)                     :: this
    integer, intent(in)                                 :: part_id
    class(graph_edge_support), allocatable, intent(out) :: support

    allocate(support, source=edge_support(owner_matches(this % eowner, this % ne, part_id, this % cut, .false.)))

  end subroutine g_borrowed_edges

  subroutine g_overlap_edges(this, part_id, support)

    class(stored_graph), intent(in)                     :: this
    integer, intent(in)                                 :: part_id
    class(graph_edge_support), allocatable, intent(out) :: support

    integer, allocatable :: owned(:), borrowed(:)

    allocate(owned   , source=owner_matches(this % eowner, this % ne, part_id, this % cut, .true.))
    allocate(borrowed, source=owner_matches(this % eowner, this % ne, part_id, this % cut, .false.))

    allocate(support, source=edge_support([owned, borrowed]))

  end subroutine g_overlap_edges

  !===================================================================!
  ! Collect the indices a part owns, or the ones it does not.
  !
  ! An uncut graph has no owner stamps to read. It answers that it
  ! owns everything and borrows nothing, which is the truthful answer
  ! for a graph that is the whole of itself.
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

  pure subroutine g_incident_edges(this, vertex_index, indices)

    class(stored_graph), intent(in)   :: this
    integer            , intent(in)   :: vertex_index
    integer, allocatable, intent(out) :: indices(:)

    indices = this % einc(this % xinc(vertex_index) : this % xinc(vertex_index + 1) - 1)

  end subroutine g_incident_edges

  pure subroutine g_adjacent_vertices(this, vertex_index, indices)

    class(stored_graph), intent(in)   :: this
    integer            , intent(in)   :: vertex_index
    integer, allocatable, intent(out) :: indices(:)

    indices = this % vadj(this % xadj(vertex_index) : this % xadj(vertex_index + 1) - 1)

  end subroutine g_adjacent_vertices

  pure subroutine g_outgoing_edges(this, vertex_index, indices)

    class(stored_graph), intent(in)   :: this
    integer            , intent(in)   :: vertex_index
    integer, allocatable, intent(out) :: indices(:)

    indices = this % eout(this % xout(vertex_index) : this % xout(vertex_index + 1) - 1)

  end subroutine g_outgoing_edges

  pure subroutine g_incoming_edges(this, vertex_index, indices)

    class(stored_graph), intent(in)   :: this
    integer            , intent(in)   :: vertex_index
    integer, allocatable, intent(out) :: indices(:)

    indices = this % ein(this % xin(vertex_index) : this % xin(vertex_index + 1) - 1)

  end subroutine g_incoming_edges

  !===================================================================!
  ! Where the outgoing edges land, and where the incoming ones came
  ! from. An edge with no head leads nowhere and is left out.
  !===================================================================!

  pure subroutine g_outgoing_vertices(this, vertex_index, indices)

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

  end subroutine g_outgoing_vertices

  pure subroutine g_incoming_vertices(this, vertex_index, indices)

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

  end subroutine g_incoming_vertices

  !===================================================================!
  ! How a part relates to the whole. An uncut graph is the whole of
  ! itself: one part, identity maps, everything owned by part one.
  !===================================================================!

  pure integer function g_num_parts(this)

    class(stored_graph), intent(in) :: this

    g_num_parts = this % nparts

  end function g_num_parts

  pure logical function g_has_part_relation(this)

    class(stored_graph), intent(in) :: this

    g_has_part_relation = this % cut

  end function g_has_part_relation

  pure integer function g_full_vertex_index(this, index)

    class(stored_graph), intent(in) :: this
    integer            , intent(in) :: index

    if (this % cut .and. allocated(this % vfull)) then
       g_full_vertex_index = this % vfull(index)
    else
       g_full_vertex_index = index
    end if

  end function g_full_vertex_index

  pure integer function g_full_edge_index(this, index)

    class(stored_graph), intent(in) :: this
    integer            , intent(in) :: index

    if (this % cut .and. allocated(this % efull)) then
       g_full_edge_index = this % efull(index)
    else
       g_full_edge_index = index
    end if

  end function g_full_edge_index

  !===================================================================!
  ! The map read backwards. Zero means the whole-graph index does not
  ! appear in that part at all.
  !===================================================================!

  pure integer function g_part_vertex_index(this, full_index, part_id)

    class(stored_graph), intent(in) :: this
    integer            , intent(in) :: full_index
    integer            , intent(in) :: part_id

    if (this % cut .and. part_id /= this % me) then
       g_part_vertex_index = 0
    else
       g_part_vertex_index = reverse_lookup(this % vfull, full_index, this % cut)
    end if

  end function g_part_vertex_index

  pure integer function g_part_edge_index(this, full_index, part_id)

    class(stored_graph), intent(in) :: this
    integer            , intent(in) :: full_index
    integer            , intent(in) :: part_id

    if (this % cut .and. part_id /= this % me) then
       g_part_edge_index = 0
    else
       g_part_edge_index = reverse_lookup(this % efull, full_index, this % cut)
    end if

  end function g_part_edge_index

  !===================================================================!
  ! Find which local index carries a given whole-graph index.
  !===================================================================!

  pure integer function reverse_lookup(full, full_index, cut)

    integer, allocatable, intent(in) :: full(:)
    integer             , intent(in) :: full_index
    logical             , intent(in) :: cut

    integer :: i

    if (.not. cut .or. .not. allocated(full)) then
       reverse_lookup = full_index
       return
    end if

    reverse_lookup = 0
    do i = 1, size(full)
       if (full(i) == full_index) then
          reverse_lookup = i
          return
       end if
    end do

  end function reverse_lookup

  pure integer function g_vertex_owner_part(this, index)

    class(stored_graph), intent(in) :: this
    integer            , intent(in) :: index

    if (this % cut .and. allocated(this % vowner)) then
       g_vertex_owner_part = this % vowner(index)
    else
       g_vertex_owner_part = 1
    end if

  end function g_vertex_owner_part

  pure integer function g_edge_owner_part(this, index)

    class(stored_graph), intent(in) :: this
    integer            , intent(in) :: index

    if (this % cut .and. allocated(this % eowner)) then
       g_edge_owner_part = this % eowner(index)
    else
       g_edge_owner_part = 1
    end if

  end function g_edge_owner_part

  !===================================================================!
  ! The data the graph came with, fetched by name.
  !
  ! get_data hands back a copy. Fetch what an operation needs once, at
  ! the top of apply, never inside a loop over faces.
  !===================================================================!

  pure logical function g_has_data(this, name)

    class(stored_graph), intent(in) :: this
    character(len=*)   , intent(in) :: name

    g_has_data = find_vertex_data(this, name) > 0 .or. find_edge_data(this, name) > 0

  end function g_has_data

  subroutine g_get_data(this, name, data)

    class(stored_graph), intent(in)             :: this
    character(len=*)   , intent(in)             :: name
    class(graph_data), allocatable, intent(out) :: data

    integer :: k

    k = find_vertex_data(this, name)
    if (k > 0) then
       allocate(data, source=this % vdata(k))
       return
    end if

    k = find_edge_data(this, name)
    if (k > 0) allocate(data, source=this % edata(k))

  end subroutine g_get_data

  !===================================================================!
  ! Which slot carries that name, or zero if none does.
  !===================================================================!

  pure integer function find_vertex_data(this, name)

    class(stored_graph), intent(in) :: this
    character(len=*)   , intent(in) :: name

    integer :: k

    find_vertex_data = 0
    if (.not. allocated(this % vdata)) return

    do k = 1, size(this % vdata)
       if (this % vdata(k) % name() == name) then
          find_vertex_data = k
          return
       end if
    end do

  end function find_vertex_data

  pure integer function find_edge_data(this, name)

    class(stored_graph), intent(in) :: this
    character(len=*)   , intent(in) :: name

    integer :: k

    find_edge_data = 0
    if (.not. allocated(this % edata)) return

    do k = 1, size(this % edata)
       if (this % edata(k) % name() == name) then
          find_edge_data = k
          return
       end if
    end do

  end function find_edge_data

end module class_graph
