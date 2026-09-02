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

     ! the numbering, kept private to this module
     procedure, private :: require_position
     procedure, private :: whole_of
     procedure, private :: within

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

    n = size(from_part)
    if (size(from_vertex) /= n .or. size(to_part) /= n .or. size(to_vertex) /= n) then
       error stop 'view_read_write: an arc has one end in each part'
    end if
    if (first_order < 0 .or. second_order < 0) then
       error stop 'view_read_write: a part has a non-negative vertex count'
    end if

    this % first_order  = first_order
    this % second_order = second_order

    allocate(tails(n), heads(n))
    do a = 1, n
       if (from_part(a) == to_part(a)) then
          error stop 'view_read_write: an arc of a bipartite digraph crosses its parts'
       end if
       call this % require_position(from_part(a), from_vertex(a))
       call this % require_position(to_part(a)  , to_vertex(a))
       tails(a) = this % whole_of(from_part(a), from_vertex(a))
       heads(a) = this % whole_of(to_part(a)  , to_vertex(a))
    end do

    this % arcs = stored_directed_graph(first_order + second_order, tails=tails, heads=heads)

  end function crossing

  !===================================================================!
  ! THE NUMBERING, both ways. A vertex of the second part continues
  ! after the first, and nothing outside this module reads that.
  !===================================================================!

  !===================================================================!
  ! A PART IS ONE OF THE TWO, and a vertex is one the part contains.
  ! Neither is defaulted: a label outside the two would otherwise be
  ! treated as the first part and give a structure that passes every
  ! structural check and relates the wrong vertices.
  !===================================================================!

  subroutine require_position(this, part, vertex)
    class(bipartite_digraph), intent(in) :: this
    integer                 , intent(in) :: part, vertex
    if (part /= FIRST_PART .and. part /= SECOND_PART) then
       error stop 'view_read_write: a part is the first or the second'
    end if
    if (vertex < 1 .or. vertex > this % order_of_part(part)) then
       error stop 'view_read_write: a vertex is one the part contains'
    end if
  end subroutine require_position

  pure integer function whole_of(this, part, vertex)
    class(bipartite_digraph), intent(in) :: this
    integer                 , intent(in) :: part, vertex
    whole_of = vertex
    if (part == SECOND_PART) whole_of = this % first_order + vertex
  end function whole_of

  pure integer function within(this, whole) result(vertex)
    class(bipartite_digraph), intent(in) :: this
    integer                 , intent(in) :: whole
    vertex = whole
    if (whole > this % first_order) vertex = whole - this % first_order
  end function within

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

  !===================================================================!
  ! THE IN-NEIGHBOURHOOD N-(v): the vertices an arc runs from into v.
  ! Every one of them lies in the other part, so no part need be
  ! returned with them.
  !===================================================================!

  subroutine in_neighbourhood(this, part, vertex, neighbours)

    class(bipartite_digraph), intent(in)  :: this
    integer                 , intent(in)  :: part, vertex
    integer, allocatable    , intent(out) :: neighbours(:)

    integer, allocatable :: incident(:)
    integer :: v, a, kept

    call this % require_position(part, vertex)
    v = this % whole_of(part, vertex)
    call this % arcs % incident_edges(v, incident)
    allocate(neighbours(size(incident)))
    kept = 0
    do a = 1, size(incident)
       if (this % arcs % edge_head(incident(a)) == v) then
          kept = kept + 1
          neighbours(kept) = this % within(this % arcs % edge_tail(incident(a)))
       end if
    end do
    neighbours = neighbours(1:kept)

  end subroutine in_neighbourhood

  !===================================================================!
  ! THE OUT-NEIGHBOURHOOD N+(v): the vertices an arc runs from v into.
  !===================================================================!

  subroutine out_neighbourhood(this, part, vertex, neighbours)

    class(bipartite_digraph), intent(in)  :: this
    integer                 , intent(in)  :: part, vertex
    integer, allocatable    , intent(out) :: neighbours(:)

    integer, allocatable :: incident(:)
    integer :: v, a, kept

    call this % require_position(part, vertex)
    v = this % whole_of(part, vertex)
    call this % arcs % incident_edges(v, incident)
    allocate(neighbours(size(incident)))
    kept = 0
    do a = 1, size(incident)
       if (this % arcs % edge_tail(incident(a)) == v) then
          kept = kept + 1
          neighbours(kept) = this % within(this % arcs % edge_head(incident(a)))
       end if
    end do
    neighbours = neighbours(1:kept)

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
  ! PREVIOUS: the vertices of u's own part that reach it by a
  ! directed path of length two. Where the part is a set of
  ! operations, these are the ones that must run before u, because
  ! each writes something u reads.
  !===================================================================!

  subroutine previous(this, part, vertex, vertices)

    class(bipartite_digraph), intent(in)  :: this
    integer                 , intent(in)  :: part, vertex
    integer, allocatable    , intent(out) :: vertices(:)

    call two_hops(this, part, vertex, .true., vertices)

  end subroutine previous

  !===================================================================!
  ! NEXT: the vertices of u's own part that u reaches by a directed
  ! path of length two - the ones that must run after it, because
  ! each reads something u writes.
  !===================================================================!

  subroutine next(this, part, vertex, vertices)

    class(bipartite_digraph), intent(in)  :: this
    integer                 , intent(in)  :: part, vertex
    integer, allocatable    , intent(out) :: vertices(:)

    call two_hops(this, part, vertex, .false., vertices)

  end subroutine next

  !===================================================================!
  ! Both directions of the same traversal: one arc into the other
  ! part, then one more the same way. A vertex is never its own
  ! previous or next, so u is removed if the traversal returns to it.
  !===================================================================!

  subroutine two_hops(this, part, vertex, backwards, vertices)

    class(bipartite_digraph), intent(in)  :: this
    integer                 , intent(in)  :: part, vertex
    logical                 , intent(in)  :: backwards
    integer, allocatable    , intent(out) :: vertices(:)

    integer, allocatable :: across(:), back(:), keep(:)
    integer :: y, w, kept

    if (backwards) then
       call this % in_neighbourhood(part, vertex, across)
    else
       call this % out_neighbourhood(part, vertex, across)
    end if

    allocate(keep(this % order_of_part(part)))
    kept = 0
    do y = 1, size(across)
       if (backwards) then
          call this % in_neighbourhood(other_part(part), across(y), back)
       else
          call this % out_neighbourhood(other_part(part), across(y), back)
       end if
       do w = 1, size(back)
          if (back(w) == vertex) cycle
          if (kept > 0) then
             if (any(keep(1:kept) == back(w))) cycle
          end if
          kept = kept + 1
          keep(kept) = back(w)
       end do
    end do
    vertices = keep(1:kept)

  end subroutine two_hops

  !===================================================================!
  ! THE PROJECTION ONTO ONE PART: u -> w when a directed path of
  ! length two runs from u to w through the other part.
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

    integer, allocatable :: crossings(:), backs(:), tails(:), heads(:)
    integer :: u, y, w, n, counted, pass

    n = this % order_of_part(part)

    do pass = 1, 2
       counted = 0
       do u = 1, n
          call this % out_neighbourhood(part, u, crossings)
          do y = 1, size(crossings)
             call this % out_neighbourhood(other_part(part), crossings(y), backs)
             do w = 1, size(backs)
                if (backs(w) == u) cycle
                counted = counted + 1
                if (pass == 2) then
                   tails(counted) = u
                   heads(counted) = backs(w)
                end if
             end do
          end do
       end do
       if (pass == 1) allocate(tails(counted), heads(counted))
    end do

    induced = stored_directed_graph(n, tails=tails(1:counted), heads=heads(1:counted))

  end function projection

  pure integer function other_part(part)
    integer, intent(in) :: part
    other_part = -1
    if (part == FIRST_PART ) other_part = SECOND_PART
    if (part == SECOND_PART) other_part = FIRST_PART
  end function other_part

  !===================================================================!
  ! WHETHER TWO VERTICES OF ONE PART MEET THE SAME VERTEX of the
  ! other, in either direction. Two that do not are non-adjacent in
  ! every projection, and a set of pairwise non-adjacent vertices is
  ! an INDEPENDENT SET.
  !===================================================================!

  logical function share_a_neighbour(this, part, one, other)

    class(bipartite_digraph), intent(in) :: this
    integer                 , intent(in) :: part, one, other

    integer, allocatable :: one_in(:), one_out(:), other_in(:), other_out(:)
    integer :: i

    call this % in_neighbourhood (part, one  , one_in )
    call this % out_neighbourhood(part, one  , one_out)
    call this % in_neighbourhood (part, other, other_in )
    call this % out_neighbourhood(part, other, other_out)

    share_a_neighbour = .false.
    do i = 1, size(one_in)
       if (any(other_in == one_in(i)) .or. any(other_out == one_in(i))) share_a_neighbour = .true.
    end do
    do i = 1, size(one_out)
       if (any(other_in == one_out(i)) .or. any(other_out == one_out(i))) share_a_neighbour = .true.
    end do

  end function share_a_neighbour

end module view_read_write
