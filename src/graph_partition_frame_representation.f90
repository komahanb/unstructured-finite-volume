!=====================================================================!
! THE PARTITION FRAME, AS A REPRESENTATION
!
! How a part sits inside the whole it was cut from. This is the record
! a partitioner writes and an assembler reads, and it was living, until
! now, as eight deferred bindings on the directed graph contract - where
! it did not belong, because none of the eight is a question about
! D = (V, E, tail, head).
!
! Read for what they are, the eight were never structure:
!
!     global_*_index      the EXTENSION of a subobject, part -> whole
!     part_*_index        the same map read backwards
!     *_owner_part        an integer FIELD on the part's members
!     num_parts, cut      provenance: whether this is a part at all
!
!            part 2                       the whole
!       +---------------+          1  2  3  4  5  6  7  8
!       |  1   2   3    |                |  |        |
!       +---------------+                4  5        7
!         vglobal = [4 5 7]
!
! So global_vertex_index(2) = 5: the part calls that cell 2, and the
! graph it was cut from calls it 5.
!
! WHAT IT CARRIES, AND WHY IT IS MORE THAN COORDINATES. A bare
! numbering could not say WHICH part it numbers, so a caller holding
! two frames could hand the wrong one to the wrong graph and be wrong
! in silence. The frame therefore carries the four set identities and
! their counts beside the arrays, and answers `describes(g)` - the one
! question that makes a captured frame checkable against the graph it
! claims to be about.
!
! WHAT IT MUST NOT CARRY, and does not: no set_map, no label_map, no
! inclusion_map. A representation says HOW members are numbered here.
! What a set MEANS is the caller's, held in the caller's maps, and
! passed at the semantic boundary - never smuggled inside an action.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module graph_partition_frame_representation

  use fractal_graph      , only : set_graph => graph
  use graph_directed_view, only : directed_graph

  implicit none

  private
  public :: partition_frame_representation

  !===================================================================!
  ! One frame: the part's identity within its whole, and the two
  ! coordinate arrays that carry it there and back.
  !
  ! A graph born whole is its own part: `cut` is false, one part,
  ! every index its own, everything owned here. That is not a special
  ! case bolted on - it is what an identity frame IS, and the queries
  ! below answer it without an allocated array anywhere.
  !===================================================================!

  type :: partition_frame_representation

     ! WHICH sets, by identity: the part's two domains, and the two
     ! they were cut from. Identity, never extension - the frame says
     ! which set it numbers, not who belongs to it.
     type(set_graph), private :: part_verts
     type(set_graph), private :: part_edges
     type(set_graph), private :: whole_verts
     type(set_graph), private :: whole_edges

     ! HOW MANY, on both sides. A count is not a membership question.
     integer, private :: n_part_verts  = 0
     integer, private :: n_part_edges  = 0
     integer, private :: n_whole_verts = 0
     integer, private :: n_whole_edges = 0

     ! WHICH part this is, and how many there are.
     integer, private :: number = 1
     integer, private :: nparts = 1
     logical, private :: cut    = .false.

     ! The coordinates: part-local index -> whole-graph index, and the
     ! part that answers for each member.
     integer, allocatable, private :: vglobal(:), eglobal(:)
     integer, allocatable, private :: vowner(:) , eowner(:)

   contains

     procedure :: num_parts
     procedure :: has_part_relation
     procedure :: part_id

     procedure :: global_vertex_index
     procedure :: global_edge_index
     procedure :: part_vertex_index
     procedure :: part_edge_index
     procedure :: vertex_owner_part
     procedure :: edge_owner_part

     procedure :: describes
     procedure :: whole_vertex_set
     procedure :: whole_edge_set
     procedure :: num_whole_vertices
     procedure :: num_whole_edges

  end type partition_frame_representation

  interface partition_frame_representation
     module procedure create_frame
     module procedure create_identity_frame
  end interface partition_frame_representation

contains

  !===================================================================!
  ! The cut frame: a partitioner has decided who owns what, and the
  ! part knows both its own domains and the ones it came from.
  !===================================================================!

  type(partition_frame_representation) function create_frame( &
       & part_verts, n_part_verts, part_edges, n_part_edges, &
       & whole_verts, n_whole_verts, whole_edges, n_whole_edges, &
       & number, nparts, vglobal, vowner, eglobal, eowner) result(this)

    type(set_graph), intent(in) :: part_verts, part_edges
    integer        , intent(in) :: n_part_verts, n_part_edges
    type(set_graph), intent(in) :: whole_verts, whole_edges
    integer        , intent(in) :: n_whole_verts, n_whole_edges
    integer        , intent(in) :: number, nparts
    integer        , intent(in) :: vglobal(:), vowner(:)
    integer        , intent(in) :: eglobal(:), eowner(:)

    this % part_verts    = part_verts
    this % part_edges    = part_edges
    this % whole_verts   = whole_verts
    this % whole_edges   = whole_edges
    this % n_part_verts  = n_part_verts
    this % n_part_edges  = n_part_edges
    this % n_whole_verts = n_whole_verts
    this % n_whole_edges = n_whole_edges
    this % number        = number
    this % nparts        = nparts
    this % cut           = .true.

    allocate(this % vglobal, source=vglobal)
    allocate(this % vowner , source=vowner)
    allocate(this % eglobal, source=eglobal)
    allocate(this % eowner , source=eowner)

  end function create_frame

  !===================================================================!
  ! The identity frame: a graph born whole is the whole of itself.
  ! One part, every index its own, nothing borrowed. No array is
  ! allocated, because none is needed to say `the same'.
  !===================================================================!

  type(partition_frame_representation) function create_identity_frame( &
       & verts, nverts, edges, nedges) result(this)

    type(set_graph), intent(in) :: verts, edges
    integer        , intent(in) :: nverts, nedges

    this % part_verts    = verts
    this % part_edges    = edges
    this % whole_verts   = verts
    this % whole_edges   = edges
    this % n_part_verts  = nverts
    this % n_part_edges  = nedges
    this % n_whole_verts = nverts
    this % n_whole_edges = nedges
    this % number        = 1
    this % nparts        = 1
    this % cut           = .false.

  end function create_identity_frame

  !===================================================================!
  ! Provenance: how many parts, which one this is, and whether the
  ! record is a cut at all.
  !===================================================================!

  pure integer function num_parts(this)
    class(partition_frame_representation), intent(in) :: this
    num_parts = this % nparts
  end function num_parts

  pure logical function has_part_relation(this)
    class(partition_frame_representation), intent(in) :: this
    has_part_relation = this % cut
  end function has_part_relation

  pure integer function part_id(this)
    class(partition_frame_representation), intent(in) :: this
    part_id = this % number
  end function part_id

  !===================================================================!
  ! The map out: what the whole graph calls this part's member. An
  ! uncut frame answers the index unchanged, because it is the whole.
  !===================================================================!

  pure integer function global_vertex_index(this, index)
    class(partition_frame_representation), intent(in) :: this
    integer                              , intent(in) :: index
    global_vertex_index = outward(this % vglobal, index, this % cut)
  end function global_vertex_index

  pure integer function global_edge_index(this, index)
    class(partition_frame_representation), intent(in) :: this
    integer                              , intent(in) :: index
    global_edge_index = outward(this % eglobal, index, this % cut)
  end function global_edge_index

  pure integer function outward(global, index, cut)
    integer, allocatable, intent(in) :: global(:)
    integer             , intent(in) :: index
    logical             , intent(in) :: cut
    if (cut .and. allocated(global)) then
       outward = global(index)
    else
       outward = index
    end if
  end function outward

  !===================================================================!
  ! The map home: where a whole-graph index sits inside a named part.
  ! Zero means it does not sit there at all - the honest answer for a
  ! member this part never held.
  !===================================================================!

  pure integer function part_vertex_index(this, global_index, part_id)
    class(partition_frame_representation), intent(in) :: this
    integer                              , intent(in) :: global_index, part_id
    part_vertex_index = inward(this % vglobal, global_index, part_id, &
         &                     this % number, this % cut)
  end function part_vertex_index

  pure integer function part_edge_index(this, global_index, part_id)
    class(partition_frame_representation), intent(in) :: this
    integer                              , intent(in) :: global_index, part_id
    part_edge_index = inward(this % eglobal, global_index, part_id, &
         &                   this % number, this % cut)
  end function part_edge_index

  pure integer function inward(global, global_index, asked, mine, cut)
    integer, allocatable, intent(in) :: global(:)
    integer             , intent(in) :: global_index, asked, mine
    logical             , intent(in) :: cut
    integer :: k
    inward = 0
    if (cut .and. asked /= mine) return
    if (.not. cut .or. .not. allocated(global)) then
       inward = global_index
       return
    end if
    do k = 1, size(global)
       if (global(k) == global_index) then
          inward = k
          return
       end if
    end do
  end function inward

  !===================================================================!
  ! Ownership: which part answers for a member. This is what stops a
  ! shared cell being counted twice when the parts are added back
  ! together. An uncut frame owns everything itself.
  !===================================================================!

  pure integer function vertex_owner_part(this, index)
    class(partition_frame_representation), intent(in) :: this
    integer                              , intent(in) :: index
    vertex_owner_part = owner(this % vowner, index, this % cut, this % number)
  end function vertex_owner_part

  pure integer function edge_owner_part(this, index)
    class(partition_frame_representation), intent(in) :: this
    integer                              , intent(in) :: index
    edge_owner_part = owner(this % eowner, index, this % cut, this % number)
  end function edge_owner_part

  pure integer function owner(owners, index, cut, mine)
    integer, allocatable, intent(in) :: owners(:)
    integer             , intent(in) :: index, mine
    logical             , intent(in) :: cut
    if (cut .and. allocated(owners)) then
       if (index >= 1 .and. index <= size(owners)) then
          owner = owners(index)
          return
       end if
    end if
    owner = mine
  end function owner

  !===================================================================!
  ! The whole this part came from: identity and count, for a caller
  ! that must size or name the destination of an assembly.
  !===================================================================!

  type(set_graph) function whole_vertex_set(this)
    class(partition_frame_representation), intent(in) :: this
    whole_vertex_set = this % whole_verts
  end function whole_vertex_set

  type(set_graph) function whole_edge_set(this)
    class(partition_frame_representation), intent(in) :: this
    whole_edge_set = this % whole_edges
  end function whole_edge_set

  pure integer function num_whole_vertices(this)
    class(partition_frame_representation), intent(in) :: this
    num_whole_vertices = this % n_whole_verts
  end function num_whole_vertices

  pure integer function num_whole_edges(this)
    class(partition_frame_representation), intent(in) :: this
    num_whole_edges = this % n_whole_edges
  end function num_whole_edges

  !===================================================================!
  ! THE ONE QUESTION A CAPTURED FRAME MUST ANSWER: is this the graph I
  ! am about? Identity first, because two graphs of equal size are not
  ! the same graph; counts second, because a frame written for a part
  ! of six cells cannot address a part of seven.
  !
  ! This is what lets an action HOLD a frame - bound once, checked at
  ! use - instead of asking a graph for a record that was never a
  ! question about its structure.
  !===================================================================!

  logical function describes(this, g)

    class(partition_frame_representation), intent(in) :: this
    class(directed_graph)                , intent(in) :: g

    type(set_graph) :: v, e

    v = g % vertex_set()
    e = g % edge_set()

    describes = v % same_as(this % part_verts)      .and. &
         &      e % same_as(this % part_edges)      .and. &
         &      g % num_vertices() == this % n_part_verts .and. &
         &      g % num_edges()    == this % n_part_edges

  end function describes

end module graph_partition_frame_representation
