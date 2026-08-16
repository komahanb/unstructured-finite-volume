!=====================================================================!
! LEVEL 4 OF THE NEW TOWER . THE ORDINARY GRAPH PROFILE
!
! The container is general; the ordinary directed graph is ONE
! SCHEMA read over it (AGENTS.md 16). The schema this profile
! demands is two binary relations over one pair of domains,
!
!      T  <=  E x V        the tail: every edge, exactly one
!      H  <=  E x V        the head: at most one - a boundary
!                          half-edge is an ABSENCE in H, not an
!                          imaginary member anywhere
!
! and from those two alone the whole ordinary vocabulary is
! derived:
!
!      edge_tail(e)        the one member of T's fibre at e
!      edge_head(e)        H's fibre at e, or zero for the wall
!      outgoing_edges(v)   T-preimage of v
!      incoming_edges(v)   H-preimage of v
!      incident_edges(v)   the two preimages, merged ascending
!      adjacent, out/in vertices ... all readings of T and H
!
! ONE SOURCE OF TRUTH. The profile stores no tail array, no head
! array, no adjacency of its own: every answer is read through the
! two relations' fibre views at the moment of asking. What the old
! stored_graph precomputed, this derives - and the compatibility
! suite holds the two against each other, query for query.
!
! A BORROWER, BY POLICY. The graph owns T and H in stable storage;
! this profile holds them by pointer (the ownership law of
! graph_binary_relation), so the graph must outlive the profile.
! The vertex and edge carriers ride along as copies - a copy IS
! the declared domain, and costs two integers.
!
! ORDER, CANONICAL. A relation is a set: how its tuples were handed
! in is no part of what it IS, so no profile answer may depend on
! it. Every derived list is CANONICALIZED to the edge carrier's own
! enumeration - fibres sorted by the carrier's local_index, never
! by the members' integer values, which a sparse carrier declares
! in whatever order it pleases. On a counted carrier this is the
! ascending order the old stored_graph produced from its builds.
! Hand T and H their tuples shuffled and every answer stands.
!
! TWO SCHEMA LAWS BEYOND THE FIBRES. E and V are distinct declared
! domains - an ordinary graph whose edges ARE its vertices is a
! category error, refused. And a self-loop is adjacent to nothing
! new: the old adjacency excluded the vertex itself (other /= v),
! and so does this one, while incident_edges still counts the loop
! twice, once per end, as the old build stamped it.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module graph_profile

  use graph_carrier        , only : member_set
  use graph_relation       , only : relation
  use graph_binary_relation, only : binary_relation
  use fractal_graph        , only : graph
  use graph_relational_view, only : relational_binding, &
       & num_relations, relation_at

  implicit none

  private
  public :: ordinary_graph_view, directed_adjacency_view

  type :: ordinary_graph_view

     class(member_set), allocatable, private :: verts
     class(member_set), allocatable, private :: edges

     class(binary_relation), pointer, private :: tails => null()
     class(binary_relation), pointer, private :: heads => null()

   contains

     procedure :: num_vertices
     procedure :: num_edges

     procedure :: edge_tail
     procedure :: edge_head
     procedure :: edge_has_head

     procedure :: outgoing_edges
     procedure :: incoming_edges
     procedure :: incident_edges
     procedure :: adjacent_vertices
     procedure :: outgoing_vertices
     procedure :: incoming_vertices

  end type ordinary_graph_view

  interface ordinary_graph_view
     module procedure create_view
  end interface ordinary_graph_view

  !===================================================================!
  ! THE DIRECTED ADJACENCY: the second schema this level reads, and
  ! a smaller one. Where the ordinary view wants E, V and two
  ! endpoint relations, this one interprets a single binary
  ! relation whose two slots are ONE declared domain,
  !
  !      A  <=  V x V        first slot FROM, second slot TO
  !
  ! No edge carrier is manufactured, no tail and head relations are
  ! derived, no orientation metadata rides anywhere: the tuple
  ! order of a same-domain binary relation IS the direction.
  !
  ! OWNERSHIP, AS EVERYWHERE. The constructor takes the graph and a
  ! SELECTOR - an identity used only to find the graph-owned
  ! relation. What the view borrows is g's own stable storage,
  ! never the selector object, so the selector may die the moment
  ! the view is born; the graph must outlive the view, as the
  ! ownership policy has always said.
  !
  ! The view interprets; it does not compute. Sources, sinks,
  ! reachability and ordering act on it from outside
  ! (graph_algorithms) - traversal is not storage.
  !===================================================================!

  type :: directed_adjacency_view

     class(member_set)     , allocatable, private :: over
     class(binary_relation), pointer    , private :: adjacency => null()

   contains

     procedure :: domain            => adjacency_domain
     procedure :: successors_view
     procedure :: predecessors_view

  end type directed_adjacency_view

  interface directed_adjacency_view
     module procedure create_adjacency_view
  end interface directed_adjacency_view

contains

  !===================================================================!
  ! Read the schema off a relational graph: which owned relation
  ! plays T, which plays H. Refusals guard the schema:
  !
  !      not binary            the profile reads fibres both ways
  !      mismatched domains    T and H must share E and V
  !      a tailless edge       every edge has exactly one tail
  !      a two-tailed edge     ditto
  !      a two-headed edge     at most one head
  !===================================================================!

  type(ordinary_graph_view) function create_view(g, binding, tail_at, &
       & head_at) result(this)

    type(graph)             , intent(in) :: g
    type(relational_binding), intent(in) :: binding
    integer                 , intent(in) :: tail_at
    integer                 , intent(in) :: head_at

    class(relation), pointer       :: r
    class(member_set), allocatable :: d
    integer, pointer               :: f(:)
    integer                        :: k, e

    r => relation_at(g, binding, tail_at)
    select type (r)
    class is (binary_relation)
       this % tails => r
    class default
       error stop 'graph_profile: the tail relation must be binary'
    end select

    r => relation_at(g, binding, head_at)
    select type (r)
    class is (binary_relation)
       this % heads => r
    class default
       error stop 'graph_profile: the head relation must be binary'
    end select

    this % edges = this % tails % source()
    this % verts = this % tails % target()

    if (this % edges % same_as(this % verts)) then
       error stop 'graph_profile: edges and vertices are distinct domains'
    end if

    d = this % heads % source()
    if (.not. d % same_as(this % edges)) then
       error stop 'graph_profile: the head relation must share the tail''s domains'
    end if
    d = this % heads % target()
    if (.not. d % same_as(this % verts)) then
       error stop 'graph_profile: the head relation must share the tail''s domains'
    end if

    do k = 1, this % edges % size()
       e = this % edges % member(k)
       f => this % tails % image_view(e)
       if (size(f) /= 1) then
          error stop 'graph_profile: every edge has exactly one tail'
       end if
       f => this % heads % image_view(e)
       if (size(f) > 1) then
          error stop 'graph_profile: no edge has two heads'
       end if
    end do

  end function create_view

  pure integer function num_vertices(this)

    class(ordinary_graph_view), intent(in) :: this

    num_vertices = this % verts % size()

  end function num_vertices

  pure integer function num_edges(this)

    class(ordinary_graph_view), intent(in) :: this

    num_edges = this % edges % size()

  end function num_edges

  !===================================================================!
  ! The edge ends: T's one fibre member, H's one or none. Zero
  ! speaks for the wall, as it always has.
  !===================================================================!

  integer function edge_tail(this, edge_index)

    class(ordinary_graph_view), intent(in) :: this
    integer                   , intent(in) :: edge_index

    integer, pointer :: f(:)

    f => this % tails % image_view(edge_index)
    edge_tail = f(1)

  end function edge_tail

  integer function edge_head(this, edge_index)

    class(ordinary_graph_view), intent(in) :: this
    integer                   , intent(in) :: edge_index

    integer, pointer :: f(:)

    f => this % heads % image_view(edge_index)
    if (size(f) == 0) then
       edge_head = 0
    else
       edge_head = f(1)
    end if

  end function edge_head

  logical function edge_has_head(this, edge_index)

    class(ordinary_graph_view), intent(in) :: this
    integer                   , intent(in) :: edge_index

    integer, pointer :: f(:)

    f => this % heads % image_view(edge_index)
    edge_has_head = size(f) > 0

  end function edge_has_head

  !===================================================================!
  ! The directed fibres, straight off the relations.
  !===================================================================!

  subroutine outgoing_edges(this, vertex_index, indices)

    class(ordinary_graph_view), intent(in)  :: this
    integer                   , intent(in)  :: vertex_index
    integer, allocatable      , intent(out) :: indices(:)

    call this % tails % preimage(vertex_index, indices)
    call canonical(this, indices)

  end subroutine outgoing_edges

  subroutine incoming_edges(this, vertex_index, indices)

    class(ordinary_graph_view), intent(in)  :: this
    integer                   , intent(in)  :: vertex_index
    integer, allocatable      , intent(out) :: indices(:)

    call this % heads % preimage(vertex_index, indices)
    call canonical(this, indices)

  end subroutine incoming_edges

  !===================================================================!
  ! The canonical order: the edge carrier's OWN ENUMERATION, keyed
  ! by local_index - never by the members' integer values, which a
  ! sparse carrier is free to hand out in any order it declared.
  ! For a counted carrier the two coincide; for { 30 10 20 } they
  ! do not, and the declaration wins. Fibres are small; an
  ! insertion sort is honest and allocation-free.
  !===================================================================!

  pure subroutine canonical(this, list)

    class(ordinary_graph_view), intent(in)    :: this
    integer                   , intent(inout) :: list(:)

    integer :: i, j, key, keypos

    do i = 2, size(list)
       key    = list(i)
       keypos = this % edges % local_index(key)
       j      = i - 1
       do while (j >= 1)
          if (this % edges % local_index(list(j)) <= keypos) exit
          list(j + 1) = list(j)
          j = j - 1
       end do
       list(j + 1) = key
    end do

  end subroutine canonical

  !===================================================================!
  ! Incidence is the two fibres merged ascending - the old einc
  ! order, derived instead of stored. An edge in both fibres (a
  ! self-loop) appears twice, exactly as the old build stamped it.
  !===================================================================!

  subroutine incident_edges(this, vertex_index, indices)

    class(ordinary_graph_view), intent(in)  :: this
    integer                   , intent(in)  :: vertex_index
    integer, allocatable      , intent(out) :: indices(:)

    integer, allocatable :: out_f(:), in_f(:)
    integer              :: i, j, n

    call this % outgoing_edges(vertex_index, out_f)
    call this % incoming_edges(vertex_index, in_f)

    allocate(indices(size(out_f) + size(in_f)))

    i = 1
    j = 1
    n = 0
    do while (i <= size(out_f) .or. j <= size(in_f))
       n = n + 1
       if (j > size(in_f)) then
          indices(n) = out_f(i)
          i = i + 1
       else if (i > size(out_f)) then
          indices(n) = in_f(j)
          j = j + 1
       else if (this % edges % local_index(out_f(i)) <= &
            &   this % edges % local_index(in_f(j))) then
          indices(n) = out_f(i)
          i = i + 1
       else
          indices(n) = in_f(j)
          j = j + 1
       end if
    end do

  end subroutine incident_edges

  !===================================================================!
  ! Neighbours: the far end of every incident edge, each named
  ! once, the wall contributing nothing, and - as the old adjacency
  ! always held - the vertex itself never its own neighbour, however
  ! many loops it wears.
  !===================================================================!

  subroutine adjacent_vertices(this, vertex_index, indices)

    class(ordinary_graph_view), intent(in)  :: this
    integer                   , intent(in)  :: vertex_index
    integer, allocatable      , intent(out) :: indices(:)

    integer, allocatable :: inc(:)
    integer              :: k, far, n

    call this % incident_edges(vertex_index, inc)

    allocate(indices(size(inc)))
    n = 0
    do k = 1, size(inc)
       if (this % edge_tail(inc(k)) == vertex_index) then
          far = this % edge_head(inc(k))
       else
          far = this % edge_tail(inc(k))
       end if
       if (far >= 1 .and. far /= vertex_index .and. &
            & .not. any(indices(1:n) == far)) then
          n = n + 1
          indices(n) = far
       end if
    end do
    indices = indices(1:n)

  end subroutine adjacent_vertices

  !===================================================================!
  ! Downstream and upstream: ends of the directed fibres, walls
  ! skipped, duplicates kept - exactly the old reading.
  !===================================================================!

  subroutine outgoing_vertices(this, vertex_index, indices)

    class(ordinary_graph_view), intent(in)  :: this
    integer                   , intent(in)  :: vertex_index
    integer, allocatable      , intent(out) :: indices(:)

    integer, allocatable :: f(:)
    integer              :: k, h, n

    call this % outgoing_edges(vertex_index, f)

    allocate(indices(size(f)))
    n = 0
    do k = 1, size(f)
       h = this % edge_head(f(k))
       if (h >= 1) then
          n = n + 1
          indices(n) = h
       end if
    end do
    indices = indices(1:n)

  end subroutine outgoing_vertices

  subroutine incoming_vertices(this, vertex_index, indices)

    class(ordinary_graph_view), intent(in)  :: this
    integer                   , intent(in)  :: vertex_index
    integer, allocatable      , intent(out) :: indices(:)

    integer, allocatable :: f(:)
    integer              :: k, n

    call this % incoming_edges(vertex_index, f)

    allocate(indices(size(f)))
    n = 0
    do k = 1, size(f)
       n = n + 1
       indices(n) = this % edge_tail(f(k))
    end do
    indices = indices(1:n)

  end subroutine incoming_vertices

  !===================================================================!
  ! Read the directed adjacency off a relational graph: find the
  ! graph-owned relation the selector names, and refuse everything
  ! the interpretation cannot mean. The graph may lawfully CONTAIN
  ! all the refused shapes - only this reading refuses them:
  !
  !      not owned by the graph      no borrow without an owner
  !      not binary                  an adjacency is a pair
  !      two different domains       an adjacency runs over one
  !===================================================================!

  type(directed_adjacency_view) function create_adjacency_view(g, &
       & binding, selector) result(this)

    type(graph)             , intent(in) :: g
    type(relational_binding), intent(in) :: binding
    class(relation)         , intent(in) :: selector

    class(relation), pointer       :: rp
    class(member_set), allocatable :: s, t
    integer                        :: k
    logical                        :: found

    found = .false.
    do k = 1, num_relations(g)
       rp => relation_at(g, binding, k)
       if (rp % same_as(selector)) then
          found = .true.
          exit
       end if
    end do
    if (.not. found) then
       error stop 'graph_profile: the graph does not own the selected relation'
    end if

    select type (rp)
    class is (binary_relation)
       this % adjacency => rp
    class default
       error stop 'graph_profile: a directed adjacency reads a binary relation'
    end select

    s = this % adjacency % source()
    t = this % adjacency % target()
    if (.not. s % same_as(t)) then
       error stop 'graph_profile: a directed adjacency runs over one domain'
    end if

    allocate(this % over, source=s)

  end function create_adjacency_view

  !===================================================================!
  ! The one domain, as a copy - the same declared domain.
  !===================================================================!

  function adjacency_domain(this) result(domain)

    class(directed_adjacency_view), intent(in) :: this
    class(member_set), allocatable             :: domain

    allocate(domain, source=this % over)

  end function adjacency_domain

  !===================================================================!
  ! The fibres of the interpretation: where a member points, and
  ! what points at it - borrows into the graph-owned index, costs
  ! as the binary citizen states them.
  !===================================================================!

  function successors_view(this, member) result(fibre)

    class(directed_adjacency_view), intent(in) :: this
    integer                       , intent(in) :: member
    integer, pointer                           :: fibre(:)

    fibre => this % adjacency % image_view(member)

  end function successors_view

  function predecessors_view(this, member) result(fibre)

    class(directed_adjacency_view), intent(in) :: this
    integer                       , intent(in) :: member
    integer, pointer                           :: fibre(:)

    fibre => this % adjacency % preimage_view(member)

  end function predecessors_view

end module graph_profile
