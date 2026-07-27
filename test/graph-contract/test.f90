!=====================================================================!
! The graph-contract suite exercises the concrete classes that stand
! behind abstract_graph_types, with no mesh and no solver in sight.
!
! It checks:
!   1. Supports: a set of ids goes in and the same set comes back, the
!      size agrees, each kind names itself, and an empty support
!      answers with a zero-length array rather than nothing at all.
!   2. The field ordering law: values run in support order with
!      components fastest inside each entry, so a value sits at
!      (entry - 1) * num_components + component. Nothing else in the
!      library can catch a mistake here - a round trip through
!      partition and assembly would carry a wrong layout out and back
!      unharmed - so it is checked directly.
!   3. The entry count is the support size, not the value count. A
!      three-cell field two components wide has three entries and six
!      values.
!   4. The value-kind rule: a field holds one kind at a time, a getter
!      for the wrong kind answers with a zero-length array instead of
!      converting or guessing, and any setter replaces the values and
!      the kind together.
!   5. Every kind survives a round trip through a field - integer,
!      real, complex, logical and character.
!   6. Every kind survives a round trip through a functional, and in
!      particular the imaginary part does. A complex-step derivative
!      IS the imaginary part, and a real-only functional would throw
!      away the very number being computed.
!   7. A functional answers a yes-or-no question as true or false
!      rather than as a one or a zero.
!   8. Graph structure: the counts, and where each edge runs. A wall
!      face has no head at all, and still knows which cell it hangs
!      off.
!   9. The named sets. Interior and boundary are worked out from the
!      structure rather than stamped on by hand - a boundary edge is
!      one with no head, and a boundary vertex is one that touches
!      such an edge. A name nothing carries returns an empty set
!      rather than failing.
!  10. Walking the graph, with and without regard to direction. The
!      four direction queries are what replace the separate directed
!      graph type: the same graph answers them.
!  11. An uncut graph admits it carries no partition, and then answers
!      every partition question as the whole of itself - one part,
!      everything owned, nothing borrowed, identity maps - rather than
!      pretending or refusing.
!  12. Geometry rides on the graph and is fetched by name, which is
!      how a flux reaches a face normal without anyone threading it
!      down through every call.
!  13. Every reduction rule on numbers whose answer is obvious by
!      hand, including a measure turning a bare sum into an integral,
!      a predicate staying a predicate, and a complex sum keeping the
!      imaginary part a complex-step derivative lives in.
!  14. The case the four steps exist for: two parts averaging 2 and 7
!      combine to 4, not 4.5, whichever order they arrive in. If that
!      ever reads 4.5, a reduction has finished early on each part and
!      a parallel run has stopped agreeing with a serial one.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_graph_contract

  use iso_fortran_env       , only : dp => REAL64
  use abstract_graph_types  , only : GRAPH_SUPPORT_VERTEX, GRAPH_SUPPORT_EDGE
  use abstract_graph_types  , only : GRAPH_FIELD_INTEGER, GRAPH_FIELD_REAL
  use abstract_graph_types  , only : GRAPH_FIELD_COMPLEX, GRAPH_FIELD_LOGICAL
  use abstract_graph_types  , only : GRAPH_FIELD_CHARACTER
  use abstract_graph_types  , only : graph_data, graph_functional
  use abstract_graph_types  , only : graph_vertex_support, graph_edge_support
  use class_graph_support   , only : vertex_support, edge_support
  use class_graph           , only : stored_graph
  use class_graph_field     , only : vertex_field, edge_field
  use class_graph_functional, only : functional
  use class_graph_reduction , only : reduction
  use class_graph_reduction , only : REDUCE_SUM, REDUCE_AVERAGE, REDUCE_MINIMUM
  use class_graph_reduction , only : REDUCE_MAXIMUM, REDUCE_NORM, REDUCE_COUNT
  use class_graph_reduction , only : REDUCE_ALL, REDUCE_ANY

  implicit none

  integer :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "graph contract suite"
  write(*,'(1x,a)') "============================================="

  call check_supports(nfail)
  call check_ordering_law(nfail)
  call check_value_kind_rule(nfail)
  call check_field_round_trips(nfail)
  call check_functional_round_trips(nfail)
  call check_graph_structure(nfail)
  call check_graph_named_sets(nfail)
  call check_graph_walking(nfail)
  call check_graph_uncut(nfail)
  call check_graph_data(nfail)
  call check_reductions(nfail)
  call check_average_across_parts(nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all graph contract checks passed"
  else
     write(*,'(1x,a,i0,a)') "FAILED: ", nfail, " graph contract check(s)"
     error stop
  end if

contains

  !===================================================================!
  ! The teller prints one PASS/FAIL line per claim and keeps a running
  ! count of the broken ones.
  !===================================================================!

  subroutine report(ok, label, nfail)

    logical         , intent(in)    :: ok
    character(len=*), intent(in)    :: label
    integer         , intent(inout) :: nfail

    if (ok) then
       write(*,'(1x,a,a)') "PASS : ", label
    else
       write(*,'(1x,a,a)') "FAIL : ", label
       nfail = nfail + 1
    end if

  end subroutine report

  !===================================================================!
  ! A support is a set of ids. It gives back what it was given, it
  ! knows how many it holds, and it says which kind of id it carries.
  !===================================================================!

  subroutine check_supports(nfail)

    integer, intent(inout) :: nfail

    type(vertex_support)  :: vs, empty
    type(edge_support)    :: es
    integer, allocatable  :: got(:)

    vs = vertex_support([7, 3, 11])
    es = edge_support([2, 5])

    call vs % vertex_ids(got)
    call report(size(got) .eq. 3, "vertex support keeps its count", nfail)
    call report(all(got .eq. [7, 3, 11]), &
         & "vertex support keeps its ids, in order", nfail)
    call report(vs % size() .eq. 3, "vertex support reports its size", nfail)
    call report(vs % kind() .eq. GRAPH_SUPPORT_VERTEX, &
         & "a vertex support names itself a vertex support", nfail)

    call es % edge_ids(got)
    call report(all(got .eq. [2, 5]), "edge support keeps its ids", nfail)
    call report(es % kind() .eq. GRAPH_SUPPORT_EDGE, &
         & "an edge support names itself an edge support", nfail)

    ! An untouched support has nothing in it, and must still answer a
    ! caller who wants to loop over it without asking first.
    call report(empty % size() .eq. 0, "an empty support has size zero", nfail)
    call empty % vertex_ids(got)
    call report(allocated(got) .and. size(got) .eq. 0, &
         & "an empty support answers with a zero-length array", nfail)

  end subroutine check_supports

  !===================================================================!
  ! THE FIELD ORDERING LAW.
  !
  ! Three cells, two components each. The law says the values lie
  ! down like this:
  !
  !      support ids     7        7        3        3       11      11
  !      component       1        2        1        2        1       2
  !                   +-------+-------+-------+-------+-------+-------+
  !      values       | 71.0  | 72.0  | 31.0  | 32.0  | 111.0 | 112.0 |
  !                   +-------+-------+-------+-------+-------+-------+
  !      position         1       2       3       4       5       6
  !
  !      position = (entry - 1) * num_components + component
  !
  ! If a field ever stored these component-slowest instead, every
  ! other check in this library would still pass: a partition and an
  ! assembly carry the wrong layout out and bring it back unharmed.
  ! The damage would surface far away, when a linear solver or a file
  ! writer read the flat vector and believed the documented order.
  ! So the law is pinned here, on its own.
  !===================================================================!

  subroutine check_ordering_law(nfail)

    integer, intent(inout) :: nfail

    type(vertex_support)  :: on
    type(vertex_field)    :: f
    real(dp), allocatable :: got(:)
    integer               :: entry_position, component, position
    logical               :: ok

    on = vertex_support([7, 3, 11])
    f  = vertex_field('q', on, ncomp=2)

    call f % set_real_vector([71.0_dp, 72.0_dp, 31.0_dp, 32.0_dp, 111.0_dp, 112.0_dp])

    call report(f % num_entries() .eq. 3, &
         & "entries counts cells, not values", nfail)
    call report(f % num_components() .eq. 2, &
         & "components counts values per cell", nfail)

    call f % get_real_vector(got)
    call report(size(got) .eq. f % num_entries() * f % num_components(), &
         & "the flat vector is entries times components long", nfail)

    ! Walk the law itself. Entry 1 is id 7, entry 2 is id 3, entry 3
    ! is id 11; the tens digit records the cell and the units digit
    ! the component, so a scrambled layout cannot pass by accident.
    ok = .true.
    do entry_position = 1, 3
       do component = 1, 2
          position = (entry_position - 1) * f % num_components() + component
          select case (entry_position)
          case (1)
             ok = ok .and. abs(got(position) - (70.0_dp + component)) < 1.0d-13
          case (2)
             ok = ok .and. abs(got(position) - (30.0_dp + component)) < 1.0d-13
          case (3)
             ok = ok .and. abs(got(position) - (110.0_dp + component)) < 1.0d-13
          end select
       end do
    end do
    call report(ok, "values sit at (entry-1)*ncomp + component", nfail)

    ! Said the other way round: the two values of one cell are
    ! neighbours in the array. Component-slowest storage would put
    ! them three apart.
    call report(abs(got(1) - 71.0_dp) < 1.0d-13 .and. &
         &      abs(got(2) - 72.0_dp) < 1.0d-13, &
         & "one cell's components are adjacent, not strided", nfail)

  end subroutine check_ordering_law

  !===================================================================!
  ! A field holds one kind of value at a time. Asking for the wrong
  ! one gets a zero-length array - not a conversion, not a guess.
  ! These accessors are pure and have no way to complain, so the empty
  ! answer is the complaint.
  !===================================================================!

  subroutine check_value_kind_rule(nfail)

    integer, intent(inout) :: nfail

    type(vertex_support)  :: on
    type(vertex_field)    :: f
    real(dp), allocatable :: r(:)
    integer , allocatable :: i(:)

    on = vertex_support([1, 2])
    f  = vertex_field('q', on)

    call f % set_real_vector([1.0_dp, 2.0_dp])
    call report(f % value_kind() .eq. GRAPH_FIELD_REAL, &
         & "setting reals makes it a real field", nfail)

    call f % get_integer_vector(i)
    call report(size(i) .eq. 0, &
         & "the wrong getter answers empty, it does not convert", nfail)

    ! Setting a different kind replaces both the values and the kind.
    call f % set_integer_vector([4, 5])
    call report(f % value_kind() .eq. GRAPH_FIELD_INTEGER, &
         & "setting integers makes it an integer field", nfail)

    call f % get_real_vector(r)
    call report(size(r) .eq. 0, &
         & "the reals are gone once the kind has changed", nfail)

    call f % get_integer_vector(i)
    call report(size(i) .eq. 2 .and. all(i .eq. [4, 5]), &
         & "the integers are there", nfail)

  end subroutine check_value_kind_rule

  !===================================================================!
  ! Every kind a field may carry goes in and comes back out.
  ! Character fields are real data here - boundary group names,
  ! material names - not a curiosity.
  !===================================================================!

  subroutine check_field_round_trips(nfail)

    integer, intent(inout) :: nfail

    type(edge_support)                :: on
    type(edge_field)                  :: f
    integer         , allocatable     :: i(:)
    real(dp)        , allocatable     :: r(:)
    complex(dp)     , allocatable     :: c(:)
    logical         , allocatable     :: l(:)
    character(len=:), allocatable     :: s(:)

    on = edge_support([1, 2])
    f  = edge_field('F', on)

    call f % set_integer_vector([3, 4])
    call f % get_integer_vector(i)
    call report(all(i .eq. [3, 4]), "an edge field carries integers", nfail)

    call f % set_real_vector([0.5_dp, 1.5_dp])
    call f % get_real_vector(r)
    call report(all(abs(r - [0.5_dp, 1.5_dp]) < 1.0d-13), &
         & "an edge field carries reals", nfail)

    call f % set_complex_vector([(1.0_dp, 2.0_dp), (3.0_dp, 4.0_dp)])
    call f % get_complex_vector(c)
    call report(abs(aimag(c(2)) - 4.0_dp) < 1.0d-13, &
         & "an edge field carries the imaginary part too", nfail)

    call f % set_logical_vector([.true., .false.])
    call f % get_logical_vector(l)
    call report(l(1) .and. .not. l(2), "an edge field carries logicals", nfail)

    call f % set_character_vector(['wall', 'duct'])
    call f % get_character_vector(s)
    call report(s(1) .eq. 'wall' .and. s(2) .eq. 'duct', &
         & "an edge field carries boundary names", nfail)
    call report(f % value_kind() .eq. GRAPH_FIELD_CHARACTER, &
         & "and knows it is holding names, not numbers", nfail)

  end subroutine check_field_round_trips

  !===================================================================!
  ! A functional is one value. The two that matter most here are
  ! complex and logical, because the old contract could carry neither.
  !===================================================================!

  subroutine check_functional_round_trips(nfail)

    integer, intent(inout) :: nfail

    type(functional)              :: j
    integer                       :: i
    real(dp)                      :: r
    complex(dp)                   :: c
    logical                       :: l
    character(len=:), allocatable :: s

    j = functional('J')

    call j % set_integer_value(7)
    call j % get_integer_value(i)
    call report(i .eq. 7, "a functional carries an integer", nfail)

    call j % set_real_value(2.5_dp)
    call j % get_real_value(r)
    call report(abs(r - 2.5_dp) < 1.0d-13, "a functional carries a real", nfail)

    ! The complex-step case. The derivative is the imaginary part, and
    ! it is tiny on purpose: a functional that quietly rounded through
    ! a real would return exactly zero here.
    call j % set_complex_value((5.0_dp, 4.0d-20))
    call j % get_complex_value(c)
    call report(abs(real(c) - 5.0_dp) < 1.0d-13, &
         & "a functional carries the real part", nfail)
    call report(abs(aimag(c) - 4.0d-20) < 1.0d-33, &
         & "a functional carries the imaginary part - complex step lives", nfail)
    call report(j % value_kind() .eq. GRAPH_FIELD_COMPLEX, &
         & "and reports itself complex", nfail)

    ! A predicate answers true or false, not one or zero.
    call j % set_logical_value(.false.)
    call j % get_logical_value(l)
    call report(.not. l, "a functional answers a yes-or-no question", nfail)
    call report(j % value_kind() .eq. GRAPH_FIELD_LOGICAL, &
         & "and reports itself logical, not numeric", nfail)

    call j % set_character_value('converged')
    call j % get_character_value(s)
    call report(s .eq. 'converged', "a functional carries a word", nfail)

  end subroutine check_functional_round_trips

  !===================================================================!
  ! The graph every check below uses. A diamond with one boundary
  ! face hanging off the bottom cell:
  !
  !                        (1)
  !                       /   \            e1: 1 -> 2
  !                      v     v           e2: 1 -> 3
  !                    (2)     (3)         e3: 2 -> 4
  !                      \     /           e4: 3 -> 4
  !                       v   v            e5: 4 -> nowhere, 'wall'
  !                        (4)
  !                         |
  !                         o   e5
  !
  ! Edge 5 has no head. It is the wall: it hangs off cell 4 and stops,
  ! and no imaginary cell sits on the far side of it.
  !===================================================================!

  type(stored_graph) function diamond() result(g)

    g = stored_graph(4, tails=[1, 1, 2, 3, 4], heads=[2, 3, 4, 4, 0], &
         &           etags=['none', 'none', 'none', 'none', 'wall'])

  end function diamond

  !===================================================================!
  ! Counts, and where each edge goes.
  !===================================================================!

  subroutine check_graph_structure(nfail)

    integer, intent(inout) :: nfail

    type(stored_graph) :: g

    g = diamond()

    call report(g % num_vertices() .eq. 4, "the graph counts its vertices", nfail)
    call report(g % num_edges() .eq. 5, "the graph counts its edges", nfail)

    call report(g % edge_tail(3) .eq. 2 .and. g % edge_head(3) .eq. 4, &
         & "an interior edge runs from its tail to its head", nfail)

    call report(g % edge_has_head(1), "an interior edge has a head", nfail)
    call report(.not. g % edge_has_head(5), &
         & "a wall face has no head - nothing on the far side", nfail)
    call report(g % edge_tail(5) .eq. 4, &
         & "and it still knows which cell it hangs off", nfail)

  end subroutine check_graph_structure

  !===================================================================!
  ! The named sets. Interior and boundary are worked out from the
  ! structure, not stamped on by hand: a boundary edge is one with no
  ! head, and a boundary vertex is one that touches such an edge.
  !===================================================================!

  subroutine check_graph_named_sets(nfail)

    integer, intent(inout) :: nfail

    type(stored_graph)                       :: g
    class(graph_vertex_support), allocatable :: vs
    class(graph_edge_support)  , allocatable :: es
    integer, allocatable                     :: ids(:)

    g = diamond()

    call g % all_vertices(vs)
    call vs % vertex_ids(ids)
    call report(size(ids) .eq. 4 .and. all(ids .eq. [1, 2, 3, 4]), &
         & "all_vertices is every cell", nfail)

    call g % boundary_vertices(vs)
    call vs % vertex_ids(ids)
    call report(size(ids) .eq. 1 .and. ids(1) .eq. 4, &
         & "only the cell behind the wall is a boundary cell", nfail)

    call g % interior_vertices(vs)
    call vs % vertex_ids(ids)
    call report(size(ids) .eq. 3 .and. all(ids .eq. [1, 2, 3]), &
         & "the other three are interior", nfail)

    call g % all_edges(es)
    call es % edge_ids(ids)
    call report(size(ids) .eq. 5, "all_edges is every face", nfail)

    call g % boundary_edges(es)
    call es % edge_ids(ids)
    call report(size(ids) .eq. 1 .and. ids(1) .eq. 5, &
         & "the headless face is the boundary face", nfail)

    call g % interior_edges(es)
    call es % edge_ids(ids)
    call report(size(ids) .eq. 4, "the other four are interior faces", nfail)

    call g % tagged_edges('wall', es)
    call es % edge_ids(ids)
    call report(size(ids) .eq. 1 .and. ids(1) .eq. 5, &
         & "the wall answers to its name", nfail)

    call g % tagged_edges('inlet', es)
    call es % edge_ids(ids)
    call report(size(ids) .eq. 0, &
         & "a name nothing carries returns an empty set, not a failure", nfail)

    ! Nothing was tagged on the vertices, so every vertex query by
    ! name must come back empty rather than reaching into thin air.
    call g % tagged_vertices('heater', vs)
    call vs % vertex_ids(ids)
    call report(size(ids) .eq. 0, &
         & "an untagged graph answers no to every vertex name", nfail)

  end subroutine check_graph_named_sets

  !===================================================================!
  ! Walking the graph, with and without regard to direction.
  !===================================================================!

  subroutine check_graph_walking(nfail)

    integer, intent(inout) :: nfail

    type(stored_graph)   :: g
    integer, allocatable :: ids(:)

    g = diamond()

    call g % incident_edges(1, ids)
    call report(size(ids) .eq. 2 .and. all(ids .eq. [1, 2]), &
         & "two faces touch the top cell", nfail)

    call g % incident_edges(4, ids)
    call report(size(ids) .eq. 3 .and. all(ids .eq. [3, 4, 5]), &
         & "three touch the bottom cell, the wall among them", nfail)

    call g % adjacent_vertices(1, ids)
    call report(size(ids) .eq. 2 .and. all(ids .eq. [2, 3]), &
         & "the top cell has two neighbours", nfail)

    call g % adjacent_vertices(4, ids)
    call report(size(ids) .eq. 2 .and. all(ids .eq. [2, 3]), &
         & "the wall leads nowhere, so it adds no neighbour", nfail)

    ! Direction. The old library needed a second kind of graph for
    ! this; here it is four more questions the same graph answers.
    call g % outgoing_edges(1, ids)
    call report(size(ids) .eq. 2, "two faces lead out of the top cell", nfail)

    call g % incoming_edges(1, ids)
    call report(size(ids) .eq. 0, "and none lead into it", nfail)

    call g % incoming_edges(4, ids)
    call report(size(ids) .eq. 2 .and. all(ids .eq. [3, 4]), &
         & "two lead into the bottom cell", nfail)

    call g % outgoing_vertices(1, ids)
    call report(size(ids) .eq. 2 .and. all(ids .eq. [2, 3]), &
         & "following them out lands on the two middle cells", nfail)

    call g % incoming_vertices(4, ids)
    call report(size(ids) .eq. 2 .and. all(ids .eq. [2, 3]), &
         & "and they came from those same two", nfail)

    call g % outgoing_vertices(4, ids)
    call report(size(ids) .eq. 0, &
         & "the wall goes out but arrives nowhere", nfail)

  end subroutine check_graph_walking

  !===================================================================!
  ! A graph straight off a mesh file was never cut. It should say so,
  ! and then answer every partition question as the whole of itself
  ! rather than pretending or refusing.
  !===================================================================!

  subroutine check_graph_uncut(nfail)

    integer, intent(inout) :: nfail

    type(stored_graph)                       :: g
    class(graph_vertex_support), allocatable :: vs
    integer, allocatable                     :: ids(:)

    g = diamond()

    call report(.not. g % has_part_relation(), &
         & "an uncut graph admits it carries no partition", nfail)
    call report(g % num_parts() .eq. 1, "an uncut graph is one part", nfail)

    call g % owned_vertices(1, vs)
    call vs % vertex_ids(ids)
    call report(size(ids) .eq. 4, "and it owns every cell in it", nfail)

    call g % borrowed_vertices(1, vs)
    call vs % vertex_ids(ids)
    call report(size(ids) .eq. 0, "and borrows none", nfail)

    call report(g % vertex_owner_part(3) .eq. 1, &
         & "every cell belongs to the one part", nfail)
    call report(g % full_vertex_id(3) .eq. 3, &
         & "local numbering is whole-graph numbering", nfail)
    call report(g % part_vertex_id(3, 1) .eq. 3, &
         & "and the map reads the same backwards", nfail)

  end subroutine check_graph_uncut

  !===================================================================!
  ! Geometry rides on the graph and is fetched by name. This is how a
  ! flux reaches a face normal without anyone threading it down
  ! through every call.
  !===================================================================!

  subroutine check_graph_data(nfail)

    integer, intent(inout) :: nfail

    type(stored_graph)             :: g
    type(vertex_support)           :: cells
    type(edge_support)             :: faces
    type(vertex_field)             :: volume
    type(edge_field)               :: area
    class(graph_data), allocatable :: got
    real(dp), allocatable          :: r(:)

    cells = vertex_support([1, 2, 3, 4])
    faces = edge_support([1, 2, 3, 4, 5])

    volume = vertex_field('cell_volume', cells, unit_name='m3')
    call volume % set_real_vector([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp])

    area = edge_field('face_area', faces, unit_name='m2')
    call area % set_real_vector([0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp])

    g = stored_graph(4, tails=[1, 1, 2, 3, 4], heads=[2, 3, 4, 4, 0], &
         &           vdata=[volume], edata=[area])

    call report(g % has_data('cell_volume'), "the graph carries cell volumes", nfail)
    call report(g % has_data('face_area'), "and face areas", nfail)
    call report(.not. g % has_data('face_normal'), &
         & "and says no to what it was never given", nfail)

    call g % get_data('cell_volume', got)
    call report(allocated(got), "a named fetch hands something back", nfail)
    call report(got % name() .eq. 'cell_volume', "and it is the right one", nfail)
    call report(got % units() .eq. 'm3', "carrying its units", nfail)

    select type (got)
    class is (vertex_field)
       call got % get_real_vector(r)
       call report(size(r) .eq. 4 .and. abs(r(3) - 3.0_dp) < 1.0d-13, &
            & "with the values intact", nfail)
    class default
       call report(.false., "with the values intact", nfail)
    end select

    call g % get_data('face_area', got)
    select type (got)
    class is (edge_field)
       call got % get_real_vector(r)
       call report(abs(r(5) - 0.5_dp) < 1.0d-13, &
            & "an edge field comes back an edge field", nfail)
    class default
       call report(.false., "an edge field comes back an edge field", nfail)
    end select

  end subroutine check_graph_data

  !===================================================================!
  ! The reduction rules, each on numbers whose answer is obvious by
  ! hand. A measure turns a bare sum into an integral.
  !===================================================================!

  subroutine check_reductions(nfail)

    integer, intent(inout) :: nfail

    type(stored_graph)                   :: g
    type(vertex_support)                 :: on
    type(vertex_field)                   :: f, vol
    type(reduction)                      :: rule
    class(graph_functional), allocatable :: j
    real(dp)                             :: r
    complex(dp)                          :: c
    logical                              :: l
    integer                              :: i

    g  = diamond()
    on = vertex_support([1, 2, 3, 4])

    f = vertex_field('q', on)
    call f % set_real_vector([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp])

    rule = reduction(REDUCE_SUM)
    call rule % reduce(g, f, on, j)
    call j % get_real_value(r)
    call report(abs(r - 10.0_dp) < 1.0d-13, "sum adds the values up", nfail)

    ! With a measure the same rule becomes an integral. Each cell is
    ! weighted by its volume, so the answer stops depending on how
    ! finely the mesh was cut.
    vol = vertex_field('cell_volume', on)
    call vol % set_real_vector([2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp])
    call rule % reduce(g, f, on, j, measure=vol)
    call j % get_real_value(r)
    call report(abs(r - 20.0_dp) < 1.0d-13, &
         & "a measure turns the sum into an integral", nfail)

    rule = reduction(REDUCE_MINIMUM)
    call rule % reduce(g, f, on, j)
    call j % get_real_value(r)
    call report(abs(r - 1.0_dp) < 1.0d-13, "minimum finds the smallest", nfail)

    rule = reduction(REDUCE_MAXIMUM)
    call rule % reduce(g, f, on, j)
    call j % get_real_value(r)
    call report(abs(r - 4.0_dp) < 1.0d-13, "maximum finds the largest", nfail)

    ! Three-four-five, so the root is exact.
    call f % set_real_vector([3.0_dp, 4.0_dp])
    rule = reduction(REDUCE_NORM)
    call rule % reduce(g, f, on, j)
    call j % get_real_value(r)
    call report(abs(r - 5.0_dp) < 1.0d-13, "the two-norm takes its root once", nfail)

    call f % set_real_vector([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp])
    rule = reduction(REDUCE_COUNT)
    call rule % reduce(g, f, on, j)
    call j % get_integer_value(i)
    call report(i .eq. 4, "count counts, and answers as an integer", nfail)
    call report(j % value_kind() .eq. GRAPH_FIELD_INTEGER, &
         & "and says so", nfail)

    ! A predicate over a logical field. This is the shape a question
    ! like "is this graph acyclic" comes back in.
    call f % set_logical_vector([.true., .true., .false., .true.])
    rule = reduction(REDUCE_ALL)
    call rule % reduce(g, f, on, j)
    call j % get_logical_value(l)
    call report(.not. l, "all is false when one of them is", nfail)

    rule = reduction(REDUCE_ANY)
    call rule % reduce(g, f, on, j)
    call j % get_logical_value(l)
    call report(l, "any is true when one of them is", nfail)
    call report(j % value_kind() .eq. GRAPH_FIELD_LOGICAL, &
         & "and a predicate stays a predicate, not a one or a zero", nfail)

    ! The complex-step road. The imaginary parts are tiny on purpose:
    ! a reduction that rounded through a real would return zero.
    call f % set_complex_vector([(1.0_dp, 1.0d-20), (2.0_dp, 3.0d-20)])
    rule = reduction(REDUCE_SUM)
    call rule % reduce(g, f, on, j)
    call j % get_complex_value(c)
    call report(abs(real(c) - 3.0_dp) < 1.0d-13, &
         & "summing a complex field keeps the real part", nfail)
    call report(abs(aimag(c) - 4.0d-20) < 1.0d-33, &
         & "and the imaginary part - a complex-step objective survives", nfail)

  end subroutine check_reductions

  !===================================================================!
  ! THE CASE THE FOUR STEPS EXIST FOR.
  !
  ! Two parts, averaged separately, would give
  !
  !      part 1  [2 2 2]  ->  mean 2
  !      part 2  [5 9]    ->  mean 7          (2 + 7) / 2 = 4.5   WRONG
  !
  ! but the answer is twenty divided by five, which is 4. The running
  ! sum and the running count have to travel together and the division
  ! has to happen once, at the very end.
  !
  ! If this check ever reads 4.5, a reduction has finished early on
  ! each part - and a parallel run has quietly stopped agreeing with a
  ! serial one.
  !===================================================================!

  subroutine check_average_across_parts(nfail)

    integer, intent(inout) :: nfail

    type(stored_graph)                   :: g
    type(vertex_support)                 :: s1, s2
    type(vertex_field)                   :: f1, f2
    type(reduction)                      :: rule
    class(graph_functional), allocatable :: a, b, both, j
    real(dp)                             :: r

    g    = diamond()
    rule = reduction(REDUCE_AVERAGE)

    s1 = vertex_support([1, 2, 3])
    f1 = vertex_field('q', s1)
    call f1 % set_real_vector([2.0_dp, 2.0_dp, 2.0_dp])

    s2 = vertex_support([1, 2])
    f2 = vertex_field('q', s2)
    call f2 % set_real_vector([5.0_dp, 9.0_dp])

    ! Each part folds on its own, and neither one divides.
    call rule % identity(a)
    call rule % accumulate(g, f1, s1, a)

    call rule % identity(b)
    call rule % accumulate(g, f2, s2, b)

    call rule % combine(a, b, both)
    call rule % finalize(both, j)

    call j % get_real_value(r)
    call report(abs(r - 4.0_dp) < 1.0d-13, &
         & "two parts average to 4, not 4.5 - the division waits", nfail)

    ! Joining the parts the other way round must give the same answer,
    ! or a parallel run would depend on which image finished first.
    call rule % combine(b, a, both)
    call rule % finalize(both, j)
    call j % get_real_value(r)
    call report(abs(r - 4.0_dp) < 1.0d-13, &
         & "and the order the parts arrive in makes no difference", nfail)

    ! One part alone still has to come out right.
    call rule % finalize(a, j)
    call j % get_real_value(r)
    call report(abs(r - 2.0_dp) < 1.0d-13, &
         & "one part on its own averages to its own mean", nfail)

  end subroutine check_average_across_parts

end program test_graph_contract
