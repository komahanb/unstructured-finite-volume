!=====================================================================!
! THE PARTITION RELATION
!
! One object, two opposite verbs. A cut and an assembly are not two
! records; they are one relation read in the two directions,
!
!      r  <=  S_part x S_whole
!
! and the whole point of this module is that the same value drives all
! four motions:
!
!      partition_graph(r, G)     ->  G_p        r is written here
!      partition_data (r, D_G)   ->  D_p        r read forward
!      assemble_graph (r, G_p)   ->  G          r read backward
!      assemble_data  (r, D_p)   ->  D_G        r read backward
!
! Nothing else may write it. A partition that invented one relation
! and an assembly that invented another would agree only by accident,
! and would disagree in exactly the place - a partition boundary -
! where the disagreement is hardest to see.
!
!                        WHAT THE RELATION IS
!
! Two of them, one per carrier, and they have the same shape:
!
!      r_V  <=  V_part x V_whole
!      r_E  <=  E_part x E_whole
!
!            part 2                       the whole
!       +---------------+          1  2  3  4  5  6  7  8
!       |  1   2   3    |                |  |        |
!       +---------------+                4  5        7
!         vglobal = [4 5 7]
!
! so (2, 5) is in r_V: the part calls that member 2, and the graph it
! was cut from calls it 5. Each is TOTAL on the part side - every
! member of the part has a name in the whole - and INJECTIVE - two
! members of one part never share a name. That is what makes the
! backward read a function on the image, and it is the only reason
! `assemble(partition(G)) = G` can hold at all.
!
! A relation stored as one array per carrier, therefore, and read
! forward or backward on demand. The arrays ARE the tuples; they are
! not a numbering that stands beside the mathematics.
!
!                    WHAT RIDES BESIDE THE RELATION
!
! Ownership is not part of r. It is a FIELD on the part's members,
!
!      own : S_part -> K
!
! taking values in the set of parts, and it answers a question r
! cannot: when the parts are added back together, which one speaks
! for a member that several of them hold? Exactly one, or a conserved
! quantity is counted twice. So it travels here, beside r, because
! the two are read together on every backward motion and separating
! them would let a caller hold one without the other.
!
!                   WHY IT CARRIES SET IDENTITIES TOO
!
! A bare pair of arrays could not say WHICH sets it relates, so a
! caller holding two relations could hand the wrong one to the wrong
! part and be wrong in silence. That is not hypothetical: it is the
! defect this design was built after, and it is why the four set
! identities and their counts sit beside the arrays, and why
! `describes(g)` exists - the one question that makes a relation
! checkable against the graph it claims to relate.
!
!      ONE RELATION PER PART. NEVER ONE RELATION ACROSS TWO PARTS.
!
!                       WHAT IT MUST NOT CARRY
!
! No set_map, no label_map, no inclusion_map, and it has none. A
! representation says HOW members are numbered here. What a set MEANS
! is the caller's, held in the caller's maps, and passed at the
! semantic boundary - never smuggled inside an action.
!
! WHY partition_relation AND NOT relation. graph_relation already
! owns `relation` and `stored_relation`, and they are a different
! thing: they have identity - they declare, they sign, they answer
! same_as - and they are built THROUGH a set map that describes their
! domains. This one has no identity, is copied freely into an action,
! and must be constructible where no map exists yet: inside the cut,
! before anything has described the part's carriers. Two public types
! named `relation` would be the ambiguity, not the economy.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module graph_partition_relation

  use fractal_graph      , only : set_graph => graph
  use graph_directed_view, only : directed_graph

  implicit none

  private
  public :: partition_relation

  !===================================================================!
  ! One relation: which sets it relates, how many members on each
  ! side, and the tuples themselves - one array per carrier.
  !
  ! A graph born whole stands in the IDENTITY relation to itself:
  ! one part, every member its own name, everything owned here. That
  ! is not a special case bolted on - it is what the identity relation
  ! IS, and the queries below answer it without an allocated array
  ! anywhere.
  !===================================================================!

  type :: partition_relation

     ! WHICH sets are related, by identity: the part's two carriers,
     ! and the two of the whole. Identity, never extension - it says
     ! which sets it relates, not who belongs to them.
     type(set_graph), private :: part_verts
     type(set_graph), private :: part_edges
     type(set_graph), private :: whole_verts
     type(set_graph), private :: whole_edges

     ! HOW MANY, on both sides. A count is not a membership question.
     integer, private :: n_part_verts  = 0
     integer, private :: n_part_edges  = 0
     integer, private :: n_whole_verts = 0
     integer, private :: n_whole_edges = 0

     ! WHICH part this relation is for, and how many there are.
     integer, private :: number = 1
     integer, private :: nparts = 1
     logical, private :: cut    = .false.

     ! The tuples of r_V and r_E, and the ownership field beside them.
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

  end type partition_relation

  interface partition_relation
     module procedure create_relation
     module procedure create_identity_relation
  end interface partition_relation

contains

  !===================================================================!
  ! The cut relation: a partitioner has decided who holds and who owns
  ! what, and the relation knows both sets it relates and both counts.
  !===================================================================!

  type(partition_relation) function create_relation( &
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

  end function create_relation

  !===================================================================!
  ! The identity relation: a graph born whole is the whole of itself.
  ! One part, every member its own name, nothing borrowed. No array is
  ! allocated, because none is needed to say `the same'.
  !===================================================================!

  type(partition_relation) function &
       & create_identity_relation(verts, nverts, edges, nedges) result(this)

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

  end function create_identity_relation

  !===================================================================!
  ! Provenance: how many parts, which one this relation is for, and
  ! whether there is a part-whole relation here at all.
  !===================================================================!

  pure integer function num_parts(this)
    class(partition_relation), intent(in) :: this
    num_parts = this % nparts
  end function num_parts

  pure logical function has_part_relation(this)
    class(partition_relation), intent(in) :: this
    has_part_relation = this % cut
  end function has_part_relation

  pure integer function part_id(this)
    class(partition_relation), intent(in) :: this
    part_id = this % number
  end function part_id

  !===================================================================!
  ! THE FORWARD READ: what the whole calls this part's member. The
  ! identity relation answers the index unchanged, because it is the
  ! whole.
  !===================================================================!

  pure integer function global_vertex_index(this, index)
    class(partition_relation), intent(in) :: this
    integer                                  , intent(in) :: index
    global_vertex_index = outward(this % vglobal, index, this % cut)
  end function global_vertex_index

  pure integer function global_edge_index(this, index)
    class(partition_relation), intent(in) :: this
    integer                                  , intent(in) :: index
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
  ! THE BACKWARD READ: where a whole-graph member sits inside a named
  ! part. Zero means the tuple is not in r at all - the honest answer
  ! for a member this part never held.
  !===================================================================!

  pure integer function part_vertex_index(this, global_index, part_id)
    class(partition_relation), intent(in) :: this
    integer                                  , intent(in) :: global_index, part_id
    part_vertex_index = inward(this % vglobal, global_index, part_id, &
         &                     this % number, this % cut)
  end function part_vertex_index

  pure integer function part_edge_index(this, global_index, part_id)
    class(partition_relation), intent(in) :: this
    integer                                  , intent(in) :: global_index, part_id
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
  ! Ownership: which part speaks for a member. This is what stops a
  ! shared cell being counted twice when the parts are added back
  ! together. The identity relation owns everything itself.
  !===================================================================!

  pure integer function vertex_owner_part(this, index)
    class(partition_relation), intent(in) :: this
    integer                                  , intent(in) :: index
    vertex_owner_part = owner(this % vowner, index, this % cut, this % number)
  end function vertex_owner_part

  pure integer function edge_owner_part(this, index)
    class(partition_relation), intent(in) :: this
    integer                                  , intent(in) :: index
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
  ! The far side of the relation: identity and count, for a caller
  ! that must size or name the destination of an assembly.
  !===================================================================!

  type(set_graph) function whole_vertex_set(this)
    class(partition_relation), intent(in) :: this
    whole_vertex_set = this % whole_verts
  end function whole_vertex_set

  type(set_graph) function whole_edge_set(this)
    class(partition_relation), intent(in) :: this
    whole_edge_set = this % whole_edges
  end function whole_edge_set

  pure integer function num_whole_vertices(this)
    class(partition_relation), intent(in) :: this
    num_whole_vertices = this % n_whole_verts
  end function num_whole_vertices

  pure integer function num_whole_edges(this)
    class(partition_relation), intent(in) :: this
    num_whole_edges = this % n_whole_edges
  end function num_whole_edges

  !===================================================================!
  ! THE ONE QUESTION A RELATION MUST ANSWER BEFORE IT IS USED: is this
  ! the part I relate? Identity first, because two graphs of equal
  ! size are not the same graph; counts second, because a relation
  ! written over a part of six members cannot address a part of seven.
  !
  ! This is what lets the four verbs take r as an ARGUMENT and still
  ! be safe - the caller may hold two relations, and handing the wrong
  ! one to a part is caught here rather than producing a wrong answer
  ! quietly.
  !===================================================================!

  logical function describes(this, g)

    class(partition_relation), intent(in) :: this
    class(directed_graph)                    , intent(in) :: g

    type(set_graph) :: v, e

    v = g % vertex_set()
    e = g % edge_set()

    describes = v % same_as(this % part_verts)      .and. &
         &      e % same_as(this % part_edges)      .and. &
         &      g % num_vertices() == this % n_part_verts .and. &
         &      g % num_edges()    == this % n_part_edges

  end function describes

end module graph_partition_relation
