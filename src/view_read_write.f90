!=====================================================================!
! A BIPARTITE DIGRAPH OVER TWO PARTS.
!
! Let B = (X u Y, A) be a digraph whose vertex set is the disjoint
! union of two PARTS, and whose arcs all cross between them: no arc
! joins two vertices of the same part. Such a digraph is BIPARTITE,
! and the parts are its PARTITE SETS.
!
!             DIRECTION IS THE ONLY INFORMATION NEEDED
!
! Nothing beyond direction is needed to specify which way a value moves.
! An arc y -> x enters x; an arc x -> y leaves it. So for a vertex x
! of the first part:
!
!      N-(x)   its IN-NEIGHBOURHOOD   the y it is entered from
!      N+(x)   its OUT-NEIGHBOURHOOD  the y it points to
!
! and for a vertex y of the second part, the same two neighbourhoods
! read the other way. One relation, two views, and no second kind
! of arc.
!
!                    N-(x1)          N+(x1)
!                      |               |
!         X   . . . ( x1 ) . . . . ( x2 ) . . .        the first part
!                    /  \          /   \
!                   /    \        /     \
!         Y   . ( y1 )  ( y2 ) ( y3 )  ( y4 ) . .      the second part
!
!             THE PROJECTION
!
! Two vertices of one part are joined in the PROJECTION onto that
! part when a directed path of length two runs between them through
! the other:
!
!      x1 -> x2   in the projection onto X
!               iff  there is y with  x1 -> y -> x2
!
! The projection is a digraph on one part alone, and it is DERIVED:
! nothing states it, and it cannot disagree with the arcs it comes
! from. A caller that needs an order over one part requests the
! projection and takes its topological order.
!
!             THE TWO-HOP NEIGHBOURHOODS
!
! A vertex of one part meets another of the same part only through a
! vertex of the other. Following two arcs gives the pair
!
!      previous(u) = { w : w -> y -> u  for some y }
!      next(u)     = { w : u -> y -> w  for some y }
!
! which is the composition of the relation with itself. These are
! the in- and out-neighbourhoods of u in the projection, and they
! are DERIVED: no arc states them.
!
!         previous(u)        u        next(u)
!            ( w1 )                    ( w2 )
!               \                       /
!                \ ->  [ y1 ]  ( u ) ->/  [ y2 ] ->
!
! A caller whose first part is its operations reads them as which
! operations must run before this one and which after. Taking one
! step at a time gives the immediate ones; the transitive closure of
! `next` gives everything downstream, and of `previous` everything
! upstream.
!
!             SOURCES AND SINKS
!
! A vertex of in-degree zero is a SOURCE: nothing enters it. A vertex
! of out-degree zero is a SINK: nothing leaves it. In a bipartite
! digraph these are exactly the boundary - a source of Y is entered
! by no vertex of X, and a sink of Y points at none.
!
!             WHAT A CALLER MAY READ INTO IT
!
! This module fixes no meaning on either part. A caller that takes X
! for its operations and Y for its data reads the neighbourhoods as
! what each operation reads and writes, and the projection onto X as
! the order the operations must run in. A caller that takes them the
! other way reads the same relation the other way. The mathematics
! is the same and this module states only the mathematics.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module view_read_write

  use view_directed        , only : forward, reverse
  use view_directed_stored , only : stored_directed_graph

  implicit none

  private

  public :: bipartite_digraph, crossing
  public :: FIRST_PART, SECOND_PART

  !===================================================================!
  ! WHICH PART A VERTEX BELONGS TO. The parts are numbered, not
  ! named, because naming them would fix a meaning this module does
  ! not have.
  !===================================================================!

  integer, parameter :: FIRST_PART  = 1
  integer, parameter :: SECOND_PART = 2

  !===================================================================!
  ! THE BIPARTITE DIGRAPH. The two parts are stored as counts and the
  ! arcs as one digraph over their disjoint union: a vertex of the
  ! first part is numbered from one, a vertex of the second continues
  ! after it. Callers never read that numbering - every procedure here
  ! takes a part and a vertex within it.
  !
  !      vertex numbering   1 .. p          first part
  !                         p+1 .. p+q      second part
  !===================================================================!

  type :: bipartite_digraph

     private

     type(stored_directed_graph) :: arcs
     integer :: first_order  = 0        ! how many vertices in part one
     integer :: second_order = 0        ! how many in part two

   contains

     ! the counts
     procedure :: order_of_part
     procedure :: size_of_digraph

     ! the neighbourhoods, which is where direction is read
     procedure :: in_neighbourhood
     procedure :: out_neighbourhood
     procedure :: in_degree
     procedure :: out_degree

     ! the boundary
     procedure :: is_source
     procedure :: is_sink

     ! the two-hop neighbourhoods within one part
     procedure :: previous
     procedure :: next

     ! the digraph induced on one part alone
     procedure :: projection

     ! whether two vertices of one part share a neighbour, which is
     ! the predicate an independent set requires of them
     procedure :: share_a_neighbour

     ! the private numbering in this module
     procedure, private :: require_position
     procedure, private :: offset

  end type bipartite_digraph

  interface bipartite_digraph
     module procedure crossing
  end interface bipartite_digraph

contains

  !===================================================================!
  ! THE CROSSING ARCS. Every arc is given as the part and vertex it
  ! leaves and the part and vertex it enters. An arc whose ends lie
  ! in the same part stops the program: a digraph with such an arc is
  ! not bipartite, and this type would be inconsistent with its
  ! definition.
  !===================================================================!

  function crossing(first_order, second_order, from_part, from_vertex, &
       & to_part, to_vertex) result(this)

    integer, intent(in) :: first_order, second_order
    integer, intent(in) :: from_part(:), from_vertex(:)
    integer, intent(in) :: to_part(:), to_vertex(:)
    type(bipartite_digraph) :: this

    integer, allocatable :: tails(:), heads(:)
    integer :: a, n
    character(len=150) :: message

    n = size(from_part)
    if (size(from_vertex) /= n .or. size(to_part) /= n .or. size(to_vertex) /= n) then
       write(message,'(a,i0,a,i0,a,i0,a,i0)') 'view_read_write: crossing requires one end per &
            &arc in each array; size(from_part) = ', n, ', size(from_vertex) = ', &
            & size(from_vertex), ', size(to_part) = ', size(to_part), ', size(to_vertex) = ', &
            & size(to_vertex)
       error stop trim(message)
    end if
    if (first_order < 0 .or. second_order < 0) then
       write(message,'(a,i0,a,i0)') 'view_read_write: crossing requires a non-negative vertex &
            &count per part; first_order = ', first_order, ', second_order = ', second_order
       error stop trim(message)
    end if

    this % first_order  = first_order
    this % second_order = second_order

    allocate(tails(n), heads(n))
    do a = 1, n
       if (from_part(a) == to_part(a)) then
          write(message,'(a,i0,a,i0)') 'view_read_write: crossing requires an arc to cross &
               &parts, but arc ', a, ' has both ends in part ', from_part(a)
          error stop trim(message)
       end if
       call this % require_position(from_part(a), from_vertex(a))
       call this % require_position(to_part(a)  , to_vertex(a))
       tails(a) = this % offset(from_part(a)) + from_vertex(a)
       heads(a) = this % offset(to_part(a))   + to_vertex(a)
    end do

    this % arcs = stored_directed_graph(first_order + second_order, tails=tails, heads=heads)

  end function crossing

  !===================================================================!
  ! THE NUMBERING. A vertex of the second part continues after the
  ! first: its whole number is its number within the part plus the
  ! part's offset. Nothing outside this module reads that.
  !===================================================================!

  pure integer function offset(this, part)
    class(bipartite_digraph), intent(in) :: this
    integer                 , intent(in) :: part
    offset = 0
    if (part == SECOND_PART) offset = this % first_order
  end function offset

  !===================================================================!
  ! A PART IS ONE OF THE TWO, and a vertex is one the part contains.
  ! Neither is defaulted: a label outside the two would otherwise be
  ! treated as the first part and give a structure that passes every
  ! structural check and relates the wrong vertices.
  !===================================================================!

  subroutine require_position(this, part, vertex)
    class(bipartite_digraph), intent(in) :: this
    integer                 , intent(in) :: part, vertex
    character(len=150) :: message
    if (part /= FIRST_PART .and. part /= SECOND_PART) then
       write(message,'(a,i0,a,i0,a,i0)') 'view_read_write: require_position requires part = &
            &FIRST_PART (', FIRST_PART, ') or SECOND_PART (', SECOND_PART, '); part = ', part
       error stop trim(message)
    end if
    if (vertex < 1 .or. vertex > this % order_of_part(part)) then
       write(message,'(a,i0,a,i0)') 'view_read_write: require_position requires a vertex the &
            &part contains; vertex = ', vertex, ', order_of_part = ', this % order_of_part(part)
       error stop trim(message)
    end if
  end subroutine require_position

  pure integer function order_of_part(this, part)
    class(bipartite_digraph), intent(in) :: this
    integer                 , intent(in) :: part
    order_of_part = -1
    if (part == FIRST_PART ) order_of_part = this % first_order
    if (part == SECOND_PART) order_of_part = this % second_order
  end function order_of_part

  pure integer function size_of_digraph(this)
    class(bipartite_digraph), intent(in) :: this
    size_of_digraph = this % arcs % num_edges()
  end function size_of_digraph

  pure integer function other_part(part)
    integer, intent(in) :: part
    other_part = -1
    if (part == FIRST_PART ) other_part = SECOND_PART
    if (part == SECOND_PART) other_part = FIRST_PART
  end function other_part

  !===================================================================!
  ! ONE STEP along the arcs, or against them, from a vertex in the
  ! whole numbering: the stored graph's out- or in-neighbours, in
  ! the order the arcs were given.
  !===================================================================!

  pure subroutine step(this, whole, backwards, vertices)
    class(bipartite_digraph), intent(in)  :: this
    integer                 , intent(in)  :: whole
    logical                 , intent(in)  :: backwards
    integer, allocatable    , intent(out) :: vertices(:)
    if (backwards) then
       call this % arcs % incoming_vertices(whole, vertices)
    else
       call this % arcs % outgoing_vertices(whole, vertices)
    end if
  end subroutine step

  !===================================================================!
  ! THE IN-NEIGHBOURHOOD N-(v): the vertices an arc runs from into v,
  ! and THE OUT-NEIGHBOURHOOD N+(v): the vertices an arc runs from v
  ! into. Every one of them lies in the other part, so no part need
  ! be returned with them: one step in the whole numbering, then the
  ! other part's offset removed.
  !===================================================================!

  subroutine neighbourhood(this, part, vertex, backwards, neighbours)
    class(bipartite_digraph), intent(in)  :: this
    integer                 , intent(in)  :: part, vertex
    logical                 , intent(in)  :: backwards
    integer, allocatable    , intent(out) :: neighbours(:)
    call this % require_position(part, vertex)
    call step(this, this % offset(part) + vertex, backwards, neighbours)
    neighbours = neighbours - this % offset(other_part(part))
  end subroutine neighbourhood

  subroutine in_neighbourhood(this, part, vertex, neighbours)
    class(bipartite_digraph), intent(in)  :: this
    integer                 , intent(in)  :: part, vertex
    integer, allocatable    , intent(out) :: neighbours(:)
    call neighbourhood(this, part, vertex, .true., neighbours)
  end subroutine in_neighbourhood

  subroutine out_neighbourhood(this, part, vertex, neighbours)
    class(bipartite_digraph), intent(in)  :: this
    integer                 , intent(in)  :: part, vertex
    integer, allocatable    , intent(out) :: neighbours(:)
    call neighbourhood(this, part, vertex, .false., neighbours)
  end subroutine out_neighbourhood

  integer function in_degree(this, part, vertex)
    class(bipartite_digraph), intent(in) :: this
    integer                 , intent(in) :: part, vertex
    integer, allocatable :: neighbours(:)
    call this % in_neighbourhood(part, vertex, neighbours)
    in_degree = size(neighbours)
  end function in_degree

  integer function out_degree(this, part, vertex)
    class(bipartite_digraph), intent(in) :: this
    integer                 , intent(in) :: part, vertex
    integer, allocatable :: neighbours(:)
    call this % out_neighbourhood(part, vertex, neighbours)
    out_degree = size(neighbours)
  end function out_degree

  !===================================================================!
  ! THE BOUNDARY. A source is entered by nothing; a sink points at
  ! nothing. Both are ordinary states of a vertex and neither is an
  ! error.
  !===================================================================!

  logical function is_source(this, part, vertex)
    class(bipartite_digraph), intent(in) :: this
    integer                 , intent(in) :: part, vertex
    is_source = this % in_degree(part, vertex) == 0
  end function is_source

  logical function is_sink(this, part, vertex)
    class(bipartite_digraph), intent(in) :: this
    integer                 , intent(in) :: part, vertex
    is_sink = this % out_degree(part, vertex) == 0
  end function is_sink

  !===================================================================!
  ! THE PATHS OF LENGTH TWO from a vertex, one arc into the other part
  ! and one more the same way: the far end of each path, one entry
  ! per path, in the order the arcs are stored - within the vertex's
  ! own part, so the part's offset is removed. A path that returns to
  ! the vertex is omitted. Forward follows the arcs; backwards runs
  ! against them.
  !===================================================================!

  function two_step_endpoints(this, part, vertex, backwards) result(ends)

    class(bipartite_digraph), intent(in) :: this
    integer                 , intent(in) :: part, vertex
    logical                 , intent(in) :: backwards
    integer, allocatable :: ends(:)

    integer, allocatable :: first_neighbours(:), second_neighbours(:)
    integer :: y, u, k, total

    call this % require_position(part, vertex)
    u = this % offset(part) + vertex
    call step(this, u, backwards, first_neighbours)
    ! the paths are counted, then placed: linear in their number
    total = 0
    do y = 1, size(first_neighbours)
       call step(this, first_neighbours(y), backwards, second_neighbours)
       total = total + count(second_neighbours /= u)
    end do
    allocate(ends(total))
    total = 0
    do y = 1, size(first_neighbours)
       call step(this, first_neighbours(y), backwards, second_neighbours)
       do k = 1, size(second_neighbours)
          if (second_neighbours(k) == u) cycle
          total = total + 1
          ends(total) = second_neighbours(k) - this % offset(part)
       end do
    end do

  end function two_step_endpoints

  !===================================================================!
  ! PREVIOUS: the vertices of u's own part that reach it by a
  ! directed path of length two - each written something u reads -
  ! and NEXT: the ones u reaches - each reads something u writes.
  ! Each vertex once, at its first path.
  !===================================================================!

  subroutine previous(this, part, vertex, vertices)
    class(bipartite_digraph), intent(in)  :: this
    integer                 , intent(in)  :: part, vertex
    integer, allocatable    , intent(out) :: vertices(:)
    vertices = distinct(two_step_endpoints(this, part, vertex, .true.))
  end subroutine previous

  subroutine next(this, part, vertex, vertices)
    class(bipartite_digraph), intent(in)  :: this
    integer                 , intent(in)  :: part, vertex
    integer, allocatable    , intent(out) :: vertices(:)
    vertices = distinct(two_step_endpoints(this, part, vertex, .false.))
  end subroutine next

  pure function distinct(list) result(distinct_members)
    integer, intent(in)  :: list(:)
    integer, allocatable :: distinct_members(:)
    integer :: i
    distinct_members = pack(list, [(all(list(1:i-1) /= list(i)), i = 1, size(list))])
  end function distinct

  !===================================================================!
  ! THE PROJECTION ONTO ONE PART: u -> w once for every directed path
  ! of length two from u to w through the other part.
  !
  !      ( u ) ---> [ y ] ---> ( w )      in the digraph
  !      ( u ) -------------> ( w )       in the projection
  !
  ! It is derived, so it cannot disagree with the arcs. Its
  ! topological order is the order the part must be visited in.
  !===================================================================!

  function projection(this, part) result(induced)

    class(bipartite_digraph), intent(in) :: this
    integer                 , intent(in) :: part
    type(stored_directed_graph) :: induced

    integer, allocatable :: tails(:), heads(:), ends(:), first(:)
    integer :: u, n

    ! the paths are counted by vertex, then placed: linear in their
    ! number, where growing one array by concatenation is quadratic
    n = this % order_of_part(part)
    allocate(first(n + 1))
    first(1) = 1
    do u = 1, n
       ends = two_step_endpoints(this, part, u, .false.)
       first(u + 1) = first(u) + size(ends)
    end do
    allocate(tails(first(n + 1) - 1), heads(first(n + 1) - 1))
    do u = 1, n
       ends = two_step_endpoints(this, part, u, .false.)
       tails(first(u):first(u + 1) - 1) = u
       heads(first(u):first(u + 1) - 1) = ends
    end do

    induced = stored_directed_graph(n, tails=tails, heads=heads)

  end function projection

  !===================================================================!
  ! WHETHER TWO VERTICES OF ONE PART MEET THE SAME VERTEX of the
  ! other, in either direction: whether their neighbourhoods in the
  ! stored graph intersect. Two that do not are non-adjacent in every
  ! projection, and a set of pairwise non-adjacent vertices is an
  ! INDEPENDENT SET.
  !===================================================================!

  logical function share_a_neighbour(this, part, one, other)

    class(bipartite_digraph), intent(in) :: this
    integer                 , intent(in) :: part, one, other

    integer, allocatable :: of_one(:), of_other(:)
    integer :: i

    call this % require_position(part, one)
    call this % require_position(part, other)
    call this % arcs % adjacent_vertices(this % offset(part) + one  , of_one)
    call this % arcs % adjacent_vertices(this % offset(part) + other, of_other)

    share_a_neighbour = any([(any(of_other == of_one(i)), i = 1, size(of_one))])

  end function share_a_neighbour

end module view_read_write
