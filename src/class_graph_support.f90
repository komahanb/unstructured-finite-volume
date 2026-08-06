!=====================================================================!
! The concrete support: an edgeless graph of members.
!
! A support is a chosen set of indices, and by the tower's reading a
! set of indices IS a graph - one with members and no incidence. This
! type stores the set as a plain integer array plus one fact no index
! list can state about itself: which SIDE of its host the members
! reference.
!
!      all_vertices               tagged_edges('wall')
!      +-------------------+      +--------------+
!      | 1  2  3  4  5  6  |      | 11  14  19   |
!      +-------------------+      +--------------+
!       side = vertex               side = edge
!
! Every structure question inherited from the grammar has a
! degenerate answer here, and each one is the honest truth about an
! edgeless graph: no edges, no boundary, one part, and the members'
! names ARE the global indices.
!
! A support is requested once, when an operation begins, and then
! looped over. It is never built inside a loop over cells or faces.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_graph_support

  use graph_grammar , only : graph, graph_field
  use graph_calculus, only : graph_support
  use graph_calculus, only : GRAPH_SIDE_VERTEX, GRAPH_SIDE_EDGE

  implicit none

  private
  public :: support

  !===================================================================!
  ! One support: the members, their side, and a number for identity.
  !===================================================================!

  type, extends(graph_support) :: support

     integer :: number = 1
     integer :: my_side = GRAPH_SIDE_VERTEX
     integer, allocatable :: members(:)

   contains

     ! The calculus's one question.
     procedure :: side

     ! Identity and size.
     procedure :: id
     procedure :: num_vertices
     procedure :: num_edges

     ! Incidence: an edgeless graph has none.
     procedure :: edge_tail
     procedure :: edge_head
     procedure :: edge_has_head

     ! The named sets.
     procedure :: all_vertices
     procedure :: interior_vertices
     procedure :: boundary_vertices
     procedure :: tagged_vertices
     procedure :: all_edges      => no_edges
     procedure :: interior_edges => no_edges
     procedure :: boundary_edges => no_edges
     procedure :: tagged_edges

     ! The frame's sets.
     procedure :: owned_vertices
     procedure :: borrowed_vertices
     procedure :: overlap_vertices
     procedure :: owned_edges    => no_part_edges
     procedure :: borrowed_edges => no_part_edges
     procedure :: overlap_edges  => no_part_edges

     ! Neighbourhoods: nobody is next to anybody.
     procedure :: incident_edges     => no_neighbours
     procedure :: adjacent_vertices  => no_neighbours
     procedure :: outgoing_edges     => no_neighbours
     procedure :: incoming_edges     => no_neighbours
     procedure :: outgoing_vertices  => no_neighbours
     procedure :: incoming_vertices  => no_neighbours

     ! The frame's relations.
     procedure :: num_parts
     procedure :: has_part_relation
     procedure :: global_vertex_index
     procedure :: global_edge_index
     procedure :: part_vertex_index
     procedure :: part_edge_index
     procedure :: vertex_owner_part
     procedure :: edge_owner_part

     ! A convenience beyond the contract: the whole member list in
     ! one call, with the empty-set promise.
     procedure :: member_indices

  end type support

  interface support
     module procedure create
  end interface support

contains

  !===================================================================!
  ! Build a support from a side and a list of member indices. The
  ! support keeps a copy.
  !===================================================================!

  pure type(support) function create(side, indices, number) result(this)

    integer, intent(in)           :: side
    integer, intent(in)           :: indices(:)
    integer, intent(in), optional :: number

    this % my_side = side
    allocate(this % members, source=indices)
    if (present(number)) this % number = number

  end function create

  !===================================================================!
  ! Which side of the host the members reference.
  !===================================================================!

  pure integer function side(this)

    class(support), intent(in) :: this

    side = this % my_side

  end function side

  !===================================================================!
  ! Identity and size. The vertices of this graph are the members;
  ! edges it has none.
  !===================================================================!

  pure integer function id(this)

    class(support), intent(in) :: this

    id = this % number

  end function id

  pure integer function num_vertices(this)

    class(support), intent(in) :: this

    if (allocated(this % members)) then
       num_vertices = size(this % members)
    else
       num_vertices = 0
    end if

  end function num_vertices

  pure integer function num_edges(this)

    class(support), intent(in) :: this

    associate (u1 => this); end associate

    num_edges = 0

  end function num_edges

  !===================================================================!
  ! Incidence on an edgeless graph: there is no edge to ask about.
  ! Zero is the out-of-range answer, matching a head that is absent.
  !===================================================================!

  pure integer function edge_tail(this, edge_index)

    class(support), intent(in) :: this
    integer, intent(in) :: edge_index

    associate (u1 => this, u2 => edge_index); end associate

    edge_tail = 0

  end function edge_tail

  pure integer function edge_head(this, edge_index)

    class(support), intent(in) :: this
    integer, intent(in) :: edge_index

    associate (u1 => this, u2 => edge_index); end associate

    edge_head = 0

  end function edge_head

  pure logical function edge_has_head(this, edge_index)

    class(support), intent(in) :: this
    integer, intent(in) :: edge_index

    associate (u1 => this, u2 => edge_index); end associate

    edge_has_head = .false.

  end function edge_has_head

  !===================================================================!
  ! The named sets. All vertices is the support itself. With no
  ! edges there is no boundary, so the interior is everything and
  ! the boundary is empty. A support carries no tags, so a tagged
  ! query returns the empty set rather than a guess.
  !===================================================================!

  subroutine all_vertices(this, members)

    class(support), intent(in) :: this
    class(graph), allocatable, intent(out) :: members

    allocate(members, source=this)

  end subroutine all_vertices

  subroutine interior_vertices(this, members)

    class(support), intent(in) :: this
    class(graph), allocatable, intent(out) :: members

    allocate(members, source=this)

  end subroutine interior_vertices

  subroutine boundary_vertices(this, members)

    class(support), intent(in) :: this
    class(graph), allocatable, intent(out) :: members

    call empty_set(this % my_side, members)

  end subroutine boundary_vertices

  subroutine tagged_vertices(this, tag, members)

    class(support), intent(in) :: this
    character(len=*), intent(in) :: tag
    class(graph), allocatable, intent(out) :: members

    associate (unused => tag)
    end associate
    call empty_set(this % my_side, members)

  end subroutine tagged_vertices

  subroutine no_edges(this, members)

    class(support), intent(in) :: this
    class(graph), allocatable, intent(out) :: members

    associate (u1 => this); end associate

    call empty_set(GRAPH_SIDE_EDGE, members)

  end subroutine no_edges

  subroutine tagged_edges(this, tag, members)

    class(support), intent(in) :: this
    character(len=*), intent(in) :: tag
    class(graph), allocatable, intent(out) :: members

    associate (u1 => this, u2 => tag); end associate

    call empty_set(GRAPH_SIDE_EDGE, members)

  end subroutine tagged_edges

  !===================================================================!
  ! The frame's sets. A support is never cut: one part, which owns
  ! everything, borrows nothing, and overlaps exactly itself.
  !===================================================================!

  subroutine owned_vertices(this, part_id, members)

    class(support), intent(in) :: this
    integer, intent(in) :: part_id
    class(graph), allocatable, intent(out) :: members

    if (part_id == 1) then
       allocate(members, source=this)
    else
       call empty_set(this % my_side, members)
    end if

  end subroutine owned_vertices

  subroutine borrowed_vertices(this, part_id, members)

    class(support), intent(in) :: this
    integer, intent(in) :: part_id
    class(graph), allocatable, intent(out) :: members

    associate (unused => part_id)
    end associate
    call empty_set(this % my_side, members)

  end subroutine borrowed_vertices

  subroutine overlap_vertices(this, part_id, members)

    class(support), intent(in) :: this
    integer, intent(in) :: part_id
    class(graph), allocatable, intent(out) :: members

    call this % owned_vertices(part_id, members)

  end subroutine overlap_vertices

  subroutine no_part_edges(this, part_id, members)

    class(support), intent(in) :: this
    integer, intent(in) :: part_id
    class(graph), allocatable, intent(out) :: members

    associate (u1 => this, u2 => part_id); end associate

    call empty_set(GRAPH_SIDE_EDGE, members)

  end subroutine no_part_edges

  !===================================================================!
  ! Neighbourhoods. No edges, so every neighbourhood is empty, and
  ! the empty-set promise holds: a zero-length array, never an
  ! unallocated one.
  !===================================================================!

  pure subroutine no_neighbours(this, vertex_index, indices)

    class(support), intent(in) :: this
    integer, intent(in) :: vertex_index
    integer, allocatable, intent(out) :: indices(:)

    associate (u1 => this, u2 => vertex_index); end associate

    allocate(indices(0))

  end subroutine no_neighbours

  !===================================================================!
  ! The frame's relations. Uncut, so one part and no stored record.
  ! The members' names ARE the global indices: the p-th vertex of
  ! this graph is member p, and the whole graph it references called
  ! that member members(p). This is the self-hosting closure made
  ! literal.
  !===================================================================!

  pure integer function num_parts(this)

    class(support), intent(in) :: this

    associate (u1 => this); end associate

    num_parts = 1

  end function num_parts

  pure logical function has_part_relation(this)

    class(support), intent(in) :: this

    associate (u1 => this); end associate

    has_part_relation = .false.

  end function has_part_relation

  pure integer function global_vertex_index(this, index)

    class(support), intent(in) :: this
    integer, intent(in) :: index

    if (allocated(this % members) .and. &
         & index >= 1 .and. index <= this % num_vertices()) then
       global_vertex_index = this % members(index)
    else
       global_vertex_index = 0
    end if

  end function global_vertex_index

  pure integer function global_edge_index(this, index)

    class(support), intent(in) :: this
    integer, intent(in) :: index

    associate (u1 => this, u2 => index); end associate

    global_edge_index = 0

  end function global_edge_index

  !===================================================================!
  ! The map read backwards: where does host index g sit among the
  ! members? A search, honest about absence.
  !===================================================================!

  pure integer function part_vertex_index(this, global_index, part_id)

    class(support), intent(in) :: this
    integer, intent(in) :: global_index
    integer, intent(in) :: part_id

    integer :: p

    part_vertex_index = 0
    if (part_id /= 1) return
    if (.not. allocated(this % members)) return

    do p = 1, size(this % members)
       if (this % members(p) == global_index) then
          part_vertex_index = p
          return
       end if
    end do

  end function part_vertex_index

  pure integer function part_edge_index(this, global_index, part_id)

    class(support), intent(in) :: this
    integer, intent(in) :: global_index
    integer, intent(in) :: part_id

    associate (u1 => this, u2 => global_index, u3 => part_id); end associate

    part_edge_index = 0

  end function part_edge_index

  pure integer function vertex_owner_part(this, index)

    class(support), intent(in) :: this
    integer, intent(in) :: index

    associate (u1 => this, u2 => index); end associate

    vertex_owner_part = 1

  end function vertex_owner_part

  pure integer function edge_owner_part(this, index)

    class(support), intent(in) :: this
    integer, intent(in) :: index

    associate (u1 => this, u2 => index); end associate

    edge_owner_part = 0

  end function edge_owner_part

  !===================================================================!
  ! The whole member list in one call, with the empty-set promise:
  ! a zero-length array, never an unallocated one.
  !===================================================================!

  pure subroutine member_indices(this, indices)

    class(support), intent(in)        :: this
    integer, allocatable, intent(out) :: indices(:)

    if (allocated(this % members)) then
       indices = this % members
    else
       allocate(indices(0))
    end if

  end subroutine member_indices

  !===================================================================!
  ! The empty set on a chosen side. One home for the answer every
  ! degenerate query shares.
  !===================================================================!

  subroutine empty_set(side, members)

    integer, intent(in) :: side
    class(graph), allocatable, intent(out) :: members

    integer :: none(0)

    allocate(members, source=support(side, none))

  end subroutine empty_set

end module class_graph_support
