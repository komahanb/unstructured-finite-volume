!=====================================================================!
! THE PRODUCTION-PATTERN INTERPRETER - earned at LEVEL 6.
!
! It reads an ORDINARY production graph - whatever
! discretization_operator % dependencies() hands back - and renders
! its directed adjacency as a grid, beside the signature that grid
! came from.
!
!                    WHY THIS IS A THIRD RENDERER
!
!      L4  structural_renderer_fixture     typed relations over
!                                          member sets
!      L5  valued_renderer_fixture         a relation and a field on
!                                          its occurrence carrier
!      L6  production_pattern_renderer     an ordinary graph, as
!                                          production returns one
!
! and the arrows point one way only: L5 imports L4, L6 imports L4,
! and L4 imports neither. Level 4 must stay able to render mathematics
! that has no numbers and no production vocabulary in it at all - so
! the ordinary-graph words live here, and never there.
!
!                        THE SIGNATURE IS THE POINT
!
! Every picture this file draws carries its SIGNATURE line, because
! Level 6's whole finding is that a grid alone forgets something:
!
!      D2 : X1 -> X2            two declared carriers, three members
!                               each, and X1 is not X2
!
!      pattern : vertices -> vertices
!                               ONE declared carrier, standing in
!                               both places
!
! Those two objects can render an identical grid. They are not the
! same mathematical object, and a representation that showed only the
! grid would have thrown away exactly the difference.
!
!                      WHAT same_coordinate_pattern MEANS
!
! Given a typed relation R : X -> Y and a production graph P, it
! enumerates
!
!      columns = X's declaration order
!      rows    = Y's declaration order
!
! and compares Boolean occupancy cell by cell against P's directed
! adjacency, read at the SAME local coordinates through P's own
! vertex carrier.
!
! IT PROVES EXACTLY ONE THING:
!
!      the same Boolean occupancy, by local coordinates
!
! and it proves NONE of these:
!
!      same relation
!      same carriers
!      same signature
!      same mathematical object
!
! It is deliberately not called same_structure, and it must never be
! read as if it were. Where the shapes do not even line up - |X| or
! |Y| differing from |V| - it answers false, and coordinate_shapes_fit
! says so separately, so a caller can tell "different pattern" from
! "not comparable at all".
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module production_pattern_renderer_fixture

  use fractal_graph        , only : set_graph => graph
  use graph_set_map        , only : set_map
  use graph_label_map      , only : label_map
  use graph_relation                , only : relation
  use graph_directed_view           , only : directed_graph
  use visualization_carriers_fixture, only : label_for
  use structural_renderer_fixture   , only : picture

  implicit none

  private
  public :: pattern_picture, production_has
  public :: same_coordinate_pattern, coordinate_shapes_fit
  public :: signature_of_pattern, signature_of_relation

  character, parameter :: FILLED   = '#'
  character, parameter :: EMPTY    = '.'
  integer  , parameter :: MIN_STUB = 7

contains

  !===================================================================!
  ! Does the production graph carry a directed edge from one member to
  ! another. Asked of the graph's own edge ends, never of a stored
  ! adjacency this file keeps.
  !
  ! A headless edge - production's wall - is not an arrow to anywhere,
  ! and is skipped rather than silently read as an edge to zero.
  !===================================================================!

  logical function production_has(p, from, to)

    class(directed_graph), intent(in) :: p
    integer     , intent(in) :: from, to

    integer :: e

    production_has = .false.
    do e = 1, p % num_edges()
       if (.not. p % edge_has_head(e)) cycle
       if (p % edge_tail(e) .eq. from .and. p % edge_head(e) .eq. to) then
          production_has = .true.
          return
       end if
    end do

  end function production_has

  !===================================================================!
  ! The signature a production pattern actually has: one vertex
  ! carrier, standing in both places, named as it named itself.
  !===================================================================!

  function signature_of_pattern(p, labels) result(text)

    class(directed_graph), intent(in) :: p
    type(label_map), intent(in) :: labels

    character(len=:), allocatable :: text
    type(set_graph)             :: v

    v    = p % vertex_set()
    text = carrier_name(v, labels) // ' -> ' // carrier_name(v, labels)

  end function signature_of_pattern

  !===================================================================!
  ! The signature a typed relation has: two declared carriers, which
  ! may or may not be the same one.
  !===================================================================!

  function signature_of_relation(r, labels) result(text)

    class(relation), intent(in) :: r
    type(label_map), intent(in) :: labels

    character(len=:), allocatable  :: text
    type(set_graph) :: from, to

    from = r % domain(1)
    to   = r % domain(2)
    text = carrier_name(from, labels) // ' -> ' // carrier_name(to, labels)

  end function signature_of_relation

  !===================================================================!
  ! Can the two even be laid over one another: does the production
  ! graph have as many vertices as the relation has columns, AND as
  ! many as it has rows.
  !
  ! For a square typed relation this can hold. For a rectangular one
  ! it cannot, because |V| is a single number.
  !===================================================================!

  logical function coordinate_shapes_fit(r, p, sets)

    class(relation), intent(in) :: r
    class(directed_graph)   , intent(in) :: p
    type(set_map)  , intent(in) :: sets

    type(set_graph) :: cols, rows

    coordinate_shapes_fit = .false.
    if (r % arity() .ne. 2) return

    cols = r % domain(1)
    rows = r % domain(2)

    coordinate_shapes_fit = (sets % size_of(cols) .eq. p % num_vertices()) .and. &
         &                  (sets % size_of(rows) .eq. p % num_vertices())

  end function coordinate_shapes_fit

  !===================================================================!
  ! Boolean occupancy, by local coordinates. See the header: this
  ! proves ONLY that, and never that the two are the same object.
  !===================================================================!

  logical function same_coordinate_pattern(r, p, sets)

    class(relation), intent(in) :: r
    class(directed_graph)   , intent(in) :: p
    type(set_map)  , intent(in) :: sets

    type(set_graph) :: cols, rows
    type(set_graph)              :: verts
    integer                        :: i, j

    same_coordinate_pattern = coordinate_shapes_fit(r, p, sets)
    if (.not. same_coordinate_pattern) return

    cols  = r % domain(1)
    rows  = r % domain(2)
    verts = p % vertex_set()

    do j = 1, sets % size_of(cols)
       do i = 1, sets % size_of(rows)
          same_coordinate_pattern = same_coordinate_pattern .and. &
               & (r % has([sets % member_of(cols, j), sets % member_of(rows, i)]) .eqv. &
               &  production_has(p, sets % member_of(verts, j), sets % member_of(verts, i)))
       end do
    end do

  end function same_coordinate_pattern

  !===================================================================!
  ! The production pattern as a grid.
  !
  !      rows    = the vertex carrier, in ITS declaration order,
  !                read as the HEAD of an arrow
  !      columns = the same carrier, read as the TAIL
  !
  ! One carrier on both axes, which is precisely what the picture is
  ! here to make visible. The title and the signature stand above the
  ! grid so no reader can mistake this for a typed relation.
  !===================================================================!

  type(picture) function pattern_picture(p, title, sets, labels) result(pic)

    class(directed_graph)    , intent(in) :: p
    character(len=*), intent(in) :: title
    type(set_map)  , intent(in) :: sets
    type(label_map), intent(in) :: labels

    type(set_graph) :: verts
    integer           :: stub, wide, i, j, at, n

    verts = p % vertex_set()
    n     = sets % size_of(verts)

    stub = max(widest(verts, sets, labels), MIN_STUB) + 2
    wide = widest(verts, sets, labels) + 1

    allocate(pic % line(3 + n))
    pic % line = repeat(' ', len(pic % line))

    call put(pic % line(1), 1, title)
    call put(pic % line(2), 1, 'signature: ' // signature_of_pattern(p, labels))

    do j = 1, n
       at = stub + (j - 1) * wide
       call put(pic % line(3), at, label_for(verts, sets % member_of(verts, j), labels))
    end do

    do i = 1, n
       call put(pic % line(3 + i), 1, label_for(verts, sets % member_of(verts, i), labels))
       do j = 1, n
          at = stub + (j - 1) * wide
          if (production_has(p, sets % member_of(verts, j), sets % member_of(verts, i))) then
             call put(pic % line(3 + i), at, FILLED)
          else
             call put(pic % line(3 + i), at, EMPTY)
          end if
       end do
    end do

  end function pattern_picture

  !===================================================================!
  ! Small mechanics, kept local so Level 4 need not export any.
  !===================================================================!

  integer function widest(carrier, sets, labels)

    type(set_graph), intent(in) :: carrier
    type(set_map)  , intent(in) :: sets
    type(label_map), intent(in) :: labels

    integer :: k

    widest = 1
    do k = 1, sets % size_of(carrier)
       widest = max(widest, len(label_for(carrier, sets % member_of(carrier, k), labels)))
    end do

  end function widest

  function carrier_name(carrier, labels) result(text)

    type(set_graph), intent(in) :: carrier
    type(label_map), intent(in) :: labels

    character(len=:), allocatable :: text

    text = labels % label_of(carrier)
    if (len(text) .eq. 0) text = '?'

  end function carrier_name

  subroutine put(line, at, text)

    character(len=*), intent(inout) :: line
    integer         , intent(in)    :: at
    character(len=*), intent(in)    :: text

    if (at + len(text) - 1 .gt. len(line)) then
       error stop 'production_pattern_renderer_fixture: the picture is wider than its page'
    end if

    line(at : at + len(text) - 1) = text

  end subroutine put

end module production_pattern_renderer_fixture
