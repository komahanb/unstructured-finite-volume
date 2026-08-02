!=====================================================================!
! The graph-contract suite exercises the concrete classes that stand
! behind abstract_graph_types, with no mesh and no solver in sight.
!
! It checks:
!   1. Supports: a set of indices goes in and the same set comes back, the
!      size agrees, each kind names itself, and an empty support
!      returns a zero-length array rather than nothing at all.
!   2. Where a value sits: a field stores its values in the order the
!      support lists its cells, with the components of one cell next
!      to each other, so a value sits at
!      (entry - 1) * num_components + component. Nothing else in the
!      library can catch a mistake here - a round trip through
!      partition and assembly would carry a wrong layout out and back
!      unharmed - so it is checked directly.
!   3. The entry count is the support size, not the value count. A
!      three-cell field two components wide has three entries and six
!      values.
!   4. The value-kind rule: a field holds one kind at a time, a getter
!      for the wrong kind returns a zero-length array instead of
!      converting or inferring, and any setter replaces the values and
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
!      face has no head, and the cell it is attached to is still
!      recorded.
!   9. The named sets. Interior and boundary are worked out from the
!      structure rather than listed by hand - a boundary edge is
!      one with no head, and a boundary vertex is one that touches
!      such an edge. A tag present nowhere returns an empty set
!      rather than an error.
!  10. Walking the graph, with and without regard to direction. The
!      four direction queries act on the same graph; no separate
!      directed graph type exists.
!  11. An uncut graph reports no partition record, and every
!      partition query treats it as the whole of itself - one part,
!      everything owned, nothing borrowed, identity maps.
!  12. Geometry rides on the graph and is fetched by name, which is
!      how an edge operation reaches a face normal without any caller
!      threading it down through every call.
!  13. Every reduction rule on numbers whose answer is obvious by
!      hand, including a measure turning a bare sum into an integral,
!      a predicate staying a predicate, and a complex sum keeping the
!      imaginary part a complex-step derivative lives in.
!  14. The case the four steps exist for: two parts averaging 2 and 7
!      combine to 4, not 4.5, whichever order they arrive in. If that
!      ever reads 4.5, a reduction has finished early on each part and
!      a parallel run has stopped agreeing with a serial one.
!  15. The first law where it holds exactly: cut into one piece, the
!      round trip gives back the same cells, the same faces, the same
!      ends and the same values.
!  16. Every cell owned exactly once, under every cut rule and several
!      part counts. Own one twice and mass appears from nowhere; own
!      it never and it vanishes. Both only show up in parallel, near
!      a cut.
!  17. The pieces added together rebuild the whole field exactly. This
!      is the only form the first law can take when the contract hands
!      the assembler one part at a time.
!  18. The third law: a sum worked out on the whole graph agrees with
!      the same sum worked out piecewise and joined, provided only the
!      owned cells of each piece are folded in.
!  19. Coarsening and refinement: six cells glued in pairs make three
!      blocks, faces inside a block vanish, values average onto the
!      blocks or add when told to, and refinement holds them back
!      down to every child.
!  20. Edge contributions reduced through incidence exactly once. On a
!      closed ring with no walls the balance must sum to zero, however
!      lopsided the state and however many face terms there are.
!      Fold a face twice and the sum is wrong by that face; miss one
!      and it is wrong the other way. Neither crashes. Open the ring
!      and the sum stops cancelling, which is the negative control -
!      without it the check would also pass on a balance that returned
!      nothing but zeros.
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
  use class_graph_partitioner, only : partitioner, PARTITION_LINEAR
  use class_graph_partitioner, only : PARTITION_BREADTH_FIRST, PARTITION_ADOPTED
  use class_graph_assembler , only : assembler
  use abstract_graph_types  , only : graph, graph_vertex_field
  use class_graph_coarsener , only : coarsener, COARSEN_PAIRWISE, COARSEN_ADOPTED
  use class_graph_refiner   , only : refiner
  use class_graph_derivative, only : derivative
  use class_graph_balance   , only : balance
  use class_graph_walk      , only : walk, WALK_COLOURING, WALK_VISIT_ORDER
  use class_graph_walk      , only : WALK_COMPONENT, WALK_DEPTH

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
  call check_identity_law_one_part(nfail)
  call check_partition_covers_once(nfail)
  call check_assembly_rebuilds_data(nfail)
  call check_operation_consistency(nfail)
  call check_coarsen_and_refine(nfail)
  call check_balance_conserves(nfail)
  call check_walks(nfail)
  call check_doing_nothing_is_an_answer(nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all graph contract checks passed"
  else
     write(*,'(1x,a,i0,a)') "FAILED: ", nfail, " graph contract check(s)"
     error stop
  end if

contains

  !===================================================================!
  ! report prints one PASS/FAIL line per claim and counts the
  ! failures.
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
  ! A support is a set of indices. It returns what it was given,
  ! reports its size, and reports which kind of index it holds.
  !===================================================================!

  subroutine check_supports(nfail)

    integer, intent(inout) :: nfail

    type(vertex_support)  :: vs, empty
    type(edge_support)    :: es
    integer, allocatable  :: got(:)

    vs = vertex_support([7, 3, 11])
    es = edge_support([2, 5])

    call vs % vertex_indices(got)
    call report(size(got) .eq. 3, "vertex support keeps its count", nfail)
    call report(all(got .eq. [7, 3, 11]), &
         & "vertex support keeps its indices, in order", nfail)
    call report(vs % size() .eq. 3, "vertex support reports its size", nfail)
    call report(vs % kind() .eq. GRAPH_SUPPORT_VERTEX, &
         & "a vertex support reports the vertex kind", nfail)

    call es % edge_indices(got)
    call report(all(got .eq. [2, 5]), "edge support keeps its indices", nfail)
    call report(es % kind() .eq. GRAPH_SUPPORT_EDGE, &
         & "an edge support reports the edge kind", nfail)

    ! An untouched support has nothing in it, and must still return a
    ! zero-length array, so a loop over it requires no prior check.
    call report(empty % size() .eq. 0, "an empty support has size zero", nfail)
    call empty % vertex_indices(got)
    call report(allocated(got) .and. size(got) .eq. 0, &
         & "an empty support returns a zero-length array", nfail)

  end subroutine check_supports

  !===================================================================!
  ! THE FIELD ORDERING LAW.
  !
  ! Three cells, two components each. The rule: the values are
  ! stored in this order:
  !
  !      support indices     7        7        3        3       11      11
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

    ! Walk the law itself. Entry 1 is index 7, entry 2 is index 3, entry 3
    ! is index 11; the tens digit records the cell and the units digit
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
         & "the wrong getter returns empty, not a conversion", nfail)

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
    call report(all(i .eq. [3, 4]), "an edge field holds integers", nfail)

    call f % set_real_vector([0.5_dp, 1.5_dp])
    call f % get_real_vector(r)
    call report(all(abs(r - [0.5_dp, 1.5_dp]) < 1.0d-13), &
         & "an edge field holds reals", nfail)

    call f % set_complex_vector([(1.0_dp, 2.0_dp), (3.0_dp, 4.0_dp)])
    call f % get_complex_vector(c)
    call report(abs(aimag(c(2)) - 4.0_dp) < 1.0d-13, &
         & "an edge field holds the imaginary part too", nfail)

    call f % set_logical_vector([.true., .false.])
    call f % get_logical_vector(l)
    call report(l(1) .and. .not. l(2), "an edge field holds logicals", nfail)

    call f % set_character_vector(['wall', 'duct'])
    call f % get_character_vector(s)
    call report(s(1) .eq. 'wall' .and. s(2) .eq. 'duct', &
         & "an edge field holds boundary names", nfail)
    call report(f % value_kind() .eq. GRAPH_FIELD_CHARACTER, &
         & "and reports character kind, not numeric", nfail)

  end subroutine check_field_round_trips

  !===================================================================!
  ! A functional is one value, of any kind a field can hold. The two
  ! that matter most here are complex and logical: complex holds a
  ! complex-step derivative, and logical holds a true-or-false
  ! answer without encoding it as a number.
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
    call report(i .eq. 7, "a functional holds an integer", nfail)

    call j % set_real_value(2.5_dp)
    call j % get_real_value(r)
    call report(abs(r - 2.5_dp) < 1.0d-13, "a functional holds a real", nfail)

    ! The complex-step case. The derivative is the imaginary part, and
    ! it is tiny on purpose: a functional that quietly rounded through
    ! a real would return exactly zero here.
    call j % set_complex_value((5.0_dp, 4.0d-20))
    call j % get_complex_value(c)
    call report(abs(real(c) - 5.0_dp) < 1.0d-13, &
         & "a functional holds the real part", nfail)
    call report(abs(aimag(c) - 4.0d-20) < 1.0d-33, &
         & "a functional holds the imaginary part - complex step lives", nfail)
    call report(j % value_kind() .eq. GRAPH_FIELD_COMPLEX, &
         & "and reports itself complex", nfail)

    ! A predicate answers true or false, not one or zero.
    call j % set_logical_value(.false.)
    call j % get_logical_value(l)
    call report(.not. l, "a functional holds a true-or-false value", nfail)
    call report(j % value_kind() .eq. GRAPH_FIELD_LOGICAL, &
         & "and reports itself logical, not numeric", nfail)

    call j % set_character_value('converged')
    call j % get_character_value(s)
    call report(s .eq. 'converged', "a functional holds a word", nfail)

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
         & "and its one cell is still recorded", nfail)

  end subroutine check_graph_structure

  !===================================================================!
  ! The named sets. Interior and boundary are worked out from the
  ! structure, not listed by hand: a boundary edge is one with no
  ! head, and a boundary vertex is one that touches such an edge.
  !===================================================================!

  subroutine check_graph_named_sets(nfail)

    integer, intent(inout) :: nfail

    type(stored_graph)                       :: g
    class(graph_vertex_support), allocatable :: vs
    class(graph_edge_support)  , allocatable :: es
    integer, allocatable                     :: indices(:)

    g = diamond()

    call g % all_vertices(vs)
    call vs % vertex_indices(indices)
    call report(size(indices) .eq. 4 .and. all(indices .eq. [1, 2, 3, 4]), &
         & "all_vertices is every cell", nfail)

    call g % boundary_vertices(vs)
    call vs % vertex_indices(indices)
    call report(size(indices) .eq. 1 .and. indices(1) .eq. 4, &
         & "only the cell behind the wall is a boundary cell", nfail)

    call g % interior_vertices(vs)
    call vs % vertex_indices(indices)
    call report(size(indices) .eq. 3 .and. all(indices .eq. [1, 2, 3]), &
         & "the other three are interior", nfail)

    call g % all_edges(es)
    call es % edge_indices(indices)
    call report(size(indices) .eq. 5, "all_edges is every face", nfail)

    call g % boundary_edges(es)
    call es % edge_indices(indices)
    call report(size(indices) .eq. 1 .and. indices(1) .eq. 5, &
         & "the headless face is the boundary face", nfail)

    call g % interior_edges(es)
    call es % edge_indices(indices)
    call report(size(indices) .eq. 4, "the other four are interior faces", nfail)

    call g % tagged_edges('wall', es)
    call es % edge_indices(indices)
    call report(size(indices) .eq. 1 .and. indices(1) .eq. 5, &
         & "tagged_edges of wall returns the wall face", nfail)

    call g % tagged_edges('inlet', es)
    call es % edge_indices(indices)
    call report(size(indices) .eq. 0, &
         & "an unknown tag returns an empty set, not a failure", nfail)

    ! Nothing was tagged on the vertices, so every vertex query by
    ! name must come back empty rather than reaching into thin air.
    call g % tagged_vertices('heater', vs)
    call vs % vertex_indices(indices)
    call report(size(indices) .eq. 0, &
         & "an untagged graph returns an empty set for every vertex tag", nfail)

  end subroutine check_graph_named_sets

  !===================================================================!
  ! Walking the graph, with and without regard to direction.
  !===================================================================!

  subroutine check_graph_walking(nfail)

    integer, intent(inout) :: nfail

    type(stored_graph)   :: g
    integer, allocatable :: indices(:)

    g = diamond()

    call g % incident_edges(1, indices)
    call report(size(indices) .eq. 2 .and. all(indices .eq. [1, 2]), &
         & "two faces touch the top cell", nfail)

    call g % incident_edges(4, indices)
    call report(size(indices) .eq. 3 .and. all(indices .eq. [3, 4, 5]), &
         & "three touch the bottom cell, the wall among them", nfail)

    call g % adjacent_vertices(1, indices)
    call report(size(indices) .eq. 2 .and. all(indices .eq. [2, 3]), &
         & "the top cell has two neighbours", nfail)

    call g % adjacent_vertices(4, indices)
    call report(size(indices) .eq. 2 .and. all(indices .eq. [2, 3]), &
         & "the wall leads nowhere, so it adds no neighbour", nfail)

    ! Direction: four more queries on the same graph. No separate
    ! directed graph type exists.
    call g % outgoing_edges(1, indices)
    call report(size(indices) .eq. 2, "two faces lead out of the top cell", nfail)

    call g % incoming_edges(1, indices)
    call report(size(indices) .eq. 0, "and none lead into it", nfail)

    call g % incoming_edges(4, indices)
    call report(size(indices) .eq. 2 .and. all(indices .eq. [3, 4]), &
         & "two lead into the bottom cell", nfail)

    call g % outgoing_vertices(1, indices)
    call report(size(indices) .eq. 2 .and. all(indices .eq. [2, 3]), &
         & "following them out lands on the two middle cells", nfail)

    call g % incoming_vertices(4, indices)
    call report(size(indices) .eq. 2 .and. all(indices .eq. [2, 3]), &
         & "and they came from those same two", nfail)

    call g % outgoing_vertices(4, indices)
    call report(size(indices) .eq. 0, &
         & "the wall goes out but arrives nowhere", nfail)

  end subroutine check_graph_walking

  !===================================================================!
  ! A graph straight off a mesh file holds no partition record, and
  ! every partition query treats it as the whole of itself.
  !===================================================================!

  subroutine check_graph_uncut(nfail)

    integer, intent(inout) :: nfail

    type(stored_graph)                       :: g
    class(graph_vertex_support), allocatable :: vs
    integer, allocatable                     :: indices(:)

    g = diamond()

    call report(.not. g % has_part_relation(), &
         & "an uncut graph reports no partition record", nfail)
    call report(g % num_parts() .eq. 1, "an uncut graph is one part", nfail)

    call g % owned_vertices(1, vs)
    call vs % vertex_indices(indices)
    call report(size(indices) .eq. 4, "and it owns every cell in it", nfail)

    call g % borrowed_vertices(1, vs)
    call vs % vertex_indices(indices)
    call report(size(indices) .eq. 0, "and borrows none", nfail)

    call report(g % vertex_owner_part(3) .eq. 1, &
         & "every cell belongs to the one part", nfail)
    call report(g % full_vertex_index(3) .eq. 3, &
         & "its own numbering is the whole-graph numbering", nfail)
    call report(g % part_vertex_index(3, 1) .eq. 3, &
         & "and the map reads the same backwards", nfail)

  end subroutine check_graph_uncut

  !===================================================================!
  ! Geometry rides on the graph and is fetched by name. This is how a
  ! edge operation reaches a face normal without any caller threading
  ! it down through every call.
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

    call report(g % has_data('cell_volume'), "the graph holds cell volumes", nfail)
    call report(g % has_data('face_area'), "and face areas", nfail)
    call report(.not. g % has_data('face_normal'), &
         & "and has_data is false for absent names", nfail)

    call g % get_data('cell_volume', got)
    call report(allocated(got), "get_data returns an allocated object", nfail)
    call report(got % name() .eq. 'cell_volume', "and it is the right one", nfail)
    call report(got % units() .eq. 'm3', "holding its units", nfail)

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
         & "and reports integer kind", nfail)

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
    call rule % initialize(a)
    call rule % accumulate(g, f1, s1, a)

    call rule % initialize(b)
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

  !===================================================================!
  ! The chain the partition checks run on. Six cells in a row:
  !
  !      (1)--(2)--(3)--(4)--(5)--(6)
  !
  ! Simple enough that the right answer can be worked out by hand, and
  ! long enough that a cut has somewhere to fall.
  !===================================================================!

  type(stored_graph) function chain_of_six() result(g)

    g = stored_graph(6, tails=[1, 2, 3, 4, 5], heads=[2, 3, 4, 5, 6])

  end function chain_of_six

  !===================================================================!
  ! THE FIRST LAW, where it holds exactly.
  !
  !      assemble( partition( G ) )     ==  G
  !      assemble( partition( G, D ) )  ==  ( G, D )
  !
  ! Cut into one piece, that piece is the whole graph, and the round
  ! trip has to give back exactly what went in - same cells, same
  ! faces, same ends, same values.
  !
  ! This is the weakest case and the one worth checking first: if the
  ! maps are wrong here they are wrong everywhere, and nothing further
  ! along is worth reading.
  !===================================================================!

  subroutine check_identity_law_one_part(nfail)

    integer, intent(inout) :: nfail

    type(stored_graph)             :: g
    type(partitioner)              :: p
    type(assembler)                :: a
    class(graph), allocatable      :: part, back
    class(graph_data), allocatable :: pd, fd
    type(vertex_support)           :: on
    type(vertex_field)             :: d
    real(dp), allocatable          :: v(:)
    integer                        :: e
    logical                        :: same

    g = chain_of_six()
    p = partitioner(PARTITION_LINEAR, nparts=1, part=1)
    a = assembler()

    call report(p % defined_on_graph(g), "a partitioner accepts a real graph", nfail)

    call p % partition_graph(g, part)
    call report(part % num_vertices() .eq. 6, &
         & "cut into one, the piece has every cell", nfail)
    call report(part % num_edges() .eq. 5, "and every face", nfail)

    call report(a % defined_on_graph(part), &
         & "the piece holds its relation to the whole", nfail)

    call a % assemble_graph(part, back)
    call report(back % num_vertices() .eq. g % num_vertices(), &
         & "assemble(partition(G)) has G's cells", nfail)
    call report(back % num_edges() .eq. g % num_edges(), &
         & "and G's faces", nfail)

    same = .true.
    do e = 1, g % num_edges()
       same = same .and. back % edge_tail(e) .eq. g % edge_tail(e)
       same = same .and. back % edge_head(e) .eq. g % edge_head(e)
    end do
    call report(same, "and every face runs between the same two cells", nfail)

    ! Now the same round trip with values riding along.
    on = vertex_support([1, 2, 3, 4, 5, 6])
    d  = vertex_field('q', on)
    call d % set_real_vector([10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 50.0_dp, 60.0_dp])

    call p % partition_data(g, d, part, pd)
    call a % assemble_data(part, pd, g, fd)

    select type (fd)
    class is (vertex_field)
       call fd % get_real_vector(v)
       call report(size(v) .eq. 6, "the data returns at the right size", nfail)
       call report(all(abs(v - [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 50.0_dp, 60.0_dp]) < 1.0d-13), &
            & "assemble(partition(G,D)) gives back D unchanged", nfail)
    class default
       call report(.false., "assemble(partition(G,D)) gives back D unchanged", nfail)
    end select

  end subroutine check_identity_law_one_part

  !===================================================================!
  ! Cut into two, every cell must be owned by exactly one part.
  !
  !      (1)--(2)--(3) : (4)--(5)--(6)
  !                  \___/
  !             each side borrows one cell from the other,
  !             and neither side owns it twice
  !
  ! Owned once is what makes the sum in the next check right. Own a
  ! cell twice and mass appears from nowhere; own it never and it
  ! vanishes. Both only happen in parallel, and only near a cut.
  !===================================================================!

  subroutine check_partition_covers_once(nfail)

    integer, intent(inout) :: nfail

    type(stored_graph)                       :: g
    type(partitioner)                        :: p
    class(graph), allocatable                :: part
    class(graph_vertex_support), allocatable :: vs
    integer, allocatable                     :: indices(:)
    integer                                  :: k, l, f, times(6)

    g     = chain_of_six()
    times = 0

    do k = 1, 2
       p = partitioner(PARTITION_LINEAR, nparts=2, part=k)
       call p % partition_graph(g, part)

       call report(part % has_part_relation(), &
            & "has_part_relation is true on a cut piece", nfail)
       call report(part % num_parts() .eq. 2, "and how many pieces there are", nfail)

       ! Count how often each whole-graph cell is owned.
       call part % owned_vertices(k, vs)
       call vs % vertex_indices(indices)
       do l = 1, size(indices)
          f = part % full_vertex_index(indices(l))
          times(f) = times(f) + 1
       end do

       ! Each piece should also have borrowed the one cell across the
       ! cut, so a face term there has a value on both sides.
       call part % borrowed_vertices(k, vs)
       call vs % vertex_indices(indices)
       call report(size(indices) .eq. 1, &
            & "each piece borrows exactly one cell across the cut", nfail)
    end do

    call report(all(times .eq. 1), &
         & "every cell is owned exactly once, by exactly one piece", nfail)

    ! The same must hold however the cut was chosen. A rule that
    ! leaves a cell unclaimed, or claims one twice, breaks the sum in
    ! the next check - and it would only show up in parallel.
    call report(covers_once(g, PARTITION_BREADTH_FIRST, 2), &
         & "growing the pieces outward also covers every cell once", nfail)
    call report(covers_once(g, PARTITION_BREADTH_FIRST, 3), &
         & "and still does when cut three ways", nfail)
    call report(covers_once(g, PARTITION_LINEAR, 4), &
         & "and so does an uneven split, four ways into six cells", nfail)

  end subroutine check_partition_covers_once

  !===================================================================!
  ! Cut the graph every way the rule allows and count how often each
  ! cell comes out owned. The answer must be once, always.
  !===================================================================!

  logical function covers_once(g, rule, nparts) result(ok)

    type(stored_graph), intent(in) :: g
    integer           , intent(in) :: rule
    integer           , intent(in) :: nparts

    type(partitioner)                        :: p
    class(graph), allocatable                :: part
    class(graph_vertex_support), allocatable :: vs
    integer, allocatable                     :: indices(:)
    integer                                  :: times(g % num_vertices())
    integer                                  :: k, l, f

    times = 0

    do k = 1, nparts
       p = partitioner(rule, nparts=nparts, part=k)
       call p % partition_graph(g, part)
       call part % owned_vertices(k, vs)
       call vs % vertex_indices(indices)
       do l = 1, size(indices)
          f = part % full_vertex_index(indices(l))
          times(f) = times(f) + 1
       end do
    end do

    ok = all(times .eq. 1)

  end function covers_once

  !===================================================================!
  ! THE FIRST LAW AGAIN, in the only form one part can satisfy.
  !
  ! The contract hands the assembler a single piece, so no one call
  ! can rebuild a whole that was cut in two. What each call does is
  ! fill in its own share and leave the rest at zero. Add the shares
  ! and the whole comes back:
  !
  !      part 1  ->  [10 20 30  0  0  0]
  !      part 2  ->  [ 0  0  0 40 50 60]
  !                  ------------------
  !      sum         [10 20 30 40 50 60]   ==  D
  !
  ! That sum is only right because the owned sets do not overlap,
  ! which the previous check pins down.
  !===================================================================!

  subroutine check_assembly_rebuilds_data(nfail)

    integer, intent(inout) :: nfail

    type(stored_graph)             :: g
    type(partitioner)              :: p
    type(assembler)                :: a
    class(graph), allocatable      :: part
    class(graph_data), allocatable :: pd, fd
    type(vertex_support)           :: on
    type(vertex_field)             :: d
    real(dp), allocatable          :: v(:)
    real(dp)                       :: total(6)
    integer                        :: k

    g = chain_of_six()
    a = assembler()

    on = vertex_support([1, 2, 3, 4, 5, 6])
    d  = vertex_field('q', on)
    call d % set_real_vector([10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 50.0_dp, 60.0_dp])

    total = 0.0_dp

    do k = 1, 2
       p = partitioner(PARTITION_LINEAR, nparts=2, part=k)
       call p % partition_graph(g, part)
       call p % partition_data(g, d, part, pd)
       call a % assemble_data(part, pd, g, fd)

       select type (fd)
       class is (vertex_field)
          call fd % get_real_vector(v)
          total = total + v(1:6)
       end select
    end do

    call report(all(abs(total - [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 50.0_dp, 60.0_dp]) < 1.0d-13), &
         & "the pieces added together rebuild the whole field exactly", nfail)

  end subroutine check_assembly_rebuilds_data

  !===================================================================!
  ! THE THIRD LAW.
  !
  !      full_A( G, D )  ==  assemble( A( partition( G, D ) ) )
  !
  ! Working something out on the whole must give the same answer as
  ! working it out on each piece and putting the pieces together. Here
  ! A is a sum - the smallest operation that exercises the law - and
  ! the pieces are joined by the reduction's own combine.
  !
  ! Only the owned cells of each piece are folded in. Fold the
  ! borrowed ones too and the shared cell is counted twice - which is
  ! the whole reason a reduction has to know about ownership.
  !===================================================================!

  subroutine check_operation_consistency(nfail)

    integer, intent(inout) :: nfail

    type(stored_graph)                       :: g
    type(partitioner)                        :: p
    type(reduction)                          :: rule
    class(graph), allocatable                :: part
    class(graph_data), allocatable           :: pd
    class(graph_vertex_support), allocatable :: vs
    class(graph_functional), allocatable     :: whole, piece, running, joined
    type(vertex_support)                     :: on, owned_only
    type(vertex_field)                       :: d, owned_values
    real(dp), allocatable                    :: v(:), pick(:)
    integer, allocatable                     :: indices(:)
    real(dp)                                 :: a_whole, a_parts
    integer                                  :: k, l

    g    = chain_of_six()
    rule = reduction(REDUCE_SUM)

    on = vertex_support([1, 2, 3, 4, 5, 6])
    d  = vertex_field('q', on)
    call d % set_real_vector([10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 50.0_dp, 60.0_dp])

    ! A on the whole.
    call rule % reduce(g, d, on, whole)
    call whole % get_real_value(a_whole)
    call report(abs(a_whole - 210.0_dp) < 1.0d-13, &
         & "the sum over the whole graph is 210", nfail)

    ! A on each piece, then joined.
    call rule % initialize(running)

    do k = 1, 2
       p = partitioner(PARTITION_LINEAR, nparts=2, part=k)
       call p % partition_graph(g, part)
       call p % partition_data(g, d, part, pd)

       ! Fold in this piece's owned cells only.
       call part % owned_vertices(k, vs)
       call vs % vertex_indices(indices)

       select type (pd)
       class is (vertex_field)
          call pd % get_real_vector(v)
          allocate(pick(size(indices)))
          do l = 1, size(indices)
             pick(l) = v(indices(l))
          end do
       end select

       owned_only   = vertex_support(indices)
       owned_values = vertex_field('q', owned_only)
       call owned_values % set_real_vector(pick)
       deallocate(pick)

       call rule % initialize(piece)
       call rule % accumulate(g, owned_values, owned_only, piece)
       call rule % combine(running, piece, joined)
       call rule % initialize(running)
       call rule % combine(joined, running, running)
    end do

    call running % get_real_value(a_parts)
    call report(abs(a_parts - a_whole) < 1.0d-13, &
         & "working it out piecewise gives the same answer as on the whole", nfail)

  end subroutine check_operation_consistency

  !===================================================================!
  ! Coarsening and refinement: the transforms that change how much
  ! detail a graph holds, rather than how it is split up.
  !
  !      (1)--(2)--(3)--(4)--(5)--(6)      six cells
  !       \___/     \___/     \___/
  !         O    --    O   --    O         three blocks
  !
  ! Values averaged onto the blocks, then carried back down.
  !===================================================================!

  subroutine check_coarsen_and_refine(nfail)

    integer, intent(inout) :: nfail

    type(stored_graph)             :: g
    type(coarsener)                :: c
    type(refiner)                  :: r
    class(graph), allocatable      :: coarse, fine
    class(graph_data), allocatable :: cd, fd
    type(vertex_support)           :: on
    type(vertex_field)             :: d
    real(dp), allocatable          :: v(:)

    g = chain_of_six()
    c = coarsener(COARSEN_PAIRWISE)

    call report(c % defined_on_graph(g), "a coarsener accepts a graph with room to shrink", nfail)

    call c % coarsen_graph(g, coarse)
    call report(coarse % num_vertices() .eq. 3, &
         & "six cells glued in pairs make three blocks", nfail)
    call report(coarse % num_edges() .lt. g % num_edges(), &
         & "and some faces vanish inside the blocks", nfail)

    ! Averaging onto the blocks. Cells 1 and 2 hold 10 and 20, so
    ! their block holds 15.
    on = vertex_support([1, 2, 3, 4, 5, 6])
    d  = vertex_field('q', on)
    call d % set_real_vector([10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 50.0_dp, 60.0_dp])

    call c % coarsen_data(g, d, coarse, cd)
    select type (cd)
    class is (vertex_field)
       call cd % get_real_vector(v)
       call report(size(v) .eq. 3, "the coarse field has one value per block", nfail)
       call report(abs(v(1) - 15.0_dp) < 1.0d-13, &
            & "and each block holds the average of its cells", nfail)
    class default
       call report(.false., "and each block holds the average of its cells", nfail)
    end select

    ! Adding instead of averaging. A residual is a total, so its
    ! block holds 30, not 15. Choosing wrong here is quiet - a
    ! multigrid cycle just converges more slowly than it should.
    c = coarsener(COARSEN_PAIRWISE, average=.false.)
    call c % coarsen_data(g, d, coarse, cd)
    select type (cd)
    class is (vertex_field)
       call cd % get_real_vector(v)
       call report(abs(v(1) - 30.0_dp) < 1.0d-13, &
            & "told to add rather than average, a block holds the total", nfail)
    end select

    ! And back down. Every child starts from its parent's value.
    r = refiner(2)
    call r % refine_graph(coarse, fine)
    call report(fine % num_vertices() .eq. 6, &
         & "splitting three blocks in two gives six cells back", nfail)

    call c % coarsen_data(g, d, coarse, cd)
    call r % refine_data(coarse, cd, fine, fd)
    select type (fd)
    class is (vertex_field)
       call fd % get_real_vector(v)
       call report(size(v) .eq. 6, "the fine field has one value per cell", nfail)
       call report(abs(v(1) - v(2)) < 1.0d-13, &
            & "and both children start from the same parent value", nfail)
    class default
       call report(.false., "and both children start from the same parent value", nfail)
    end select

  end subroutine check_coarsen_and_refine

  !===================================================================!
  ! DOING NOTHING IS AN ANSWER.
  !
  ! A transform given work that is already done returns the same
  ! graph; it does not reject it. Rejecting would push a special case
  ! onto every caller - including the smallest level of a multigrid
  ! hierarchy, and a serial run of code written for parts.
  !===================================================================!

  subroutine check_doing_nothing_is_an_answer(nfail)

    integer, intent(inout) :: nfail

    type(stored_graph)             :: one_cell, g
    type(coarsener)                :: c
    type(assembler)                :: a
    type(partitioner)              :: p
    class(graph), allocatable      :: out, part, back
    integer, allocatable           :: no_edges(:)

    allocate(no_edges(0))
    one_cell = stored_graph(1, tails=no_edges, heads=no_edges)

    ! A single cell is as coarse as a graph gets.
    c = coarsener(COARSEN_PAIRWISE)
    call report(c % defined_on_graph(one_cell), &
         & "a coarsener accepts a graph already at one cell", nfail)

    call c % coarsen_graph(one_cell, out)
    call report(out % num_vertices() .eq. 1, &
         & "it returns the same single cell", nfail)

    ! A graph that was never cut is already assembled.
    g = chain_of_six()
    a = assembler()
    call report(a % defined_on_graph(g), &
         & "an assembler accepts a graph that was never cut", nfail)

    call a % assemble_graph(g, back)
    call report(back % num_vertices() .eq. 6 .and. back % num_edges() .eq. 5, &
         & "it returns the same graph", nfail)

    ! And cutting into one piece is the same story from the other end.
    p = partitioner(PARTITION_LINEAR, nparts=1, part=1)
    call p % partition_graph(g, part)
    call report(part % num_vertices() .eq. 6, &
         & "cutting into one piece returns the whole graph", nfail)

  end subroutine check_doing_nothing_is_an_answer

  !===================================================================!
  ! EDGE CONTRIBUTIONS REDUCED THROUGH INCIDENCE EXACTLY ONCE.
  !
  ! A ring of four cells with no walls anywhere:
  !
  !            (1) ---> (2)
  !             ^        |
  !             |        v
  !            (4) <--- (3)
  !
  ! Every face gives its number to the cell it leaves and takes the
  ! same number from the cell it enters. So however the numbers come
  ! out, they must cancel:
  !
  !            sum over all cells of the balance  ==  0
  !
  ! That is not an approximation and it does not depend on the
  ! physics. Walk one face twice and the sum is wrong by that face;
  ! miss one and it is wrong the other way. Neither crashes - both
  ! just make the answer quietly not the answer.
  !
  ! Checked with several different states, because a single state
  ! could cancel by luck.
  !===================================================================!

  subroutine check_balance_conserves(nfail)

    integer, intent(inout) :: nfail

    type(stored_graph)                       :: ring, open_chain
    type(derivative)                         :: edge_term
    type(balance)                            :: bal
    class(graph_vertex_field), allocatable   :: y
    type(vertex_support)                     :: on
    type(vertex_field)                       :: q
    real(dp), allocatable                    :: v(:)
    real(dp)                                 :: total
    logical                                  :: ok

    ! A closed ring: four cells, four faces, not a wall in sight.
    ring = stored_graph(4, tails=[1, 2, 3, 4], heads=[2, 3, 4, 1])

    edge_term = derivative(order=2, coefficient=1.0_dp)
    bal       = balance(edge_terms=[edge_term])

    on = vertex_support([1, 2, 3, 4])
    q  = vertex_field('q', on)

    ok = .true.

    call q % set_real_vector([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp])
    call bal % apply(ring, [q], y)
    call y % get_real_vector(v)
    total = sum(v)
    ok = ok .and. abs(total) < 1.0d-12
    call report(size(v) .eq. 4, "a balance answers one value per cell", nfail)

    ! A lopsided state, so nothing cancels by symmetry.
    call q % set_real_vector([0.0_dp, 7.0_dp, -3.0_dp, 11.5_dp])
    call bal % apply(ring, [q], y)
    call y % get_real_vector(v)
    ok = ok .and. abs(sum(v)) < 1.0d-12

    ! And one more, with the values in a different order again.
    call q % set_real_vector([100.0_dp, -0.25_dp, 4.0_dp, 0.0_dp])
    call bal % apply(ring, [q], y)
    call y % get_real_vector(v)
    ok = ok .and. abs(sum(v)) < 1.0d-12

    call report(ok, &
         & "on a closed ring the balance sums to zero - every face folded once", nfail)

    ! Two edge terms, not one. Each is folded through incidence in
    ! its turn, so the sum still has to cancel.
    bal = balance(edge_terms=[derivative(order=2, coefficient=1.0_dp), &
         &                    derivative(order=1, coefficient=0.5_dp)])
    call q % set_real_vector([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp])
    call bal % apply(ring, [q], y)
    call y % get_real_vector(v)
    call report(abs(sum(v)) < 1.0d-12, &
         & "two edge terms fold once each, and still cancel", nfail)

    ! Pass the same buffer twice. The second call replaces the result
    ! rather than adding to it, so the sum is still zero rather than
    ! double.
    call bal % apply(ring, [q], y)
    call y % get_real_vector(v)
    call report(abs(sum(v)) < 1.0d-12, &
         & "applying twice into one buffer overwrites, it does not accumulate", nfail)

    ! Open the ring and the balance stops cancelling, because now two
    ! faces are walls and their numbers leave the mesh. If the closed
    ! case passed only because the balance was returning zeros, this
    ! catches it.
    open_chain = stored_graph(4, tails=[1, 2, 3, 4], heads=[2, 3, 4, 0])
    call bal % apply(open_chain, [q], y)
    call y % get_real_vector(v)
    call report(abs(sum(v)) > 1.0d-12, &
         & "open the ring and the sum stops cancelling", nfail)

  end subroutine check_balance_conserves

  !===================================================================!
  ! The walks: colouring, visit order, components, depth.
  !
  ! Each reads the shape of the graph and returns a whole number
  ! per cell, so each is an operation over a graph rather than a
  ! procedure of one.
  !===================================================================!

  subroutine check_walks(nfail)

    integer, intent(inout) :: nfail

    type(stored_graph)                     :: g, split_mesh
    type(walk)                             :: w
    class(graph_vertex_field), allocatable :: f
    integer, allocatable                   :: c(:)
    integer                                :: e, t, h
    logical                                :: ok

    g = chain_of_six()

    ! COLOURING. The promise is not "few colours" - it is that no face
    ! has the same colour at both ends. That is what makes a colour
    ! safe to sweep in parallel, so that is what gets checked.
    w = walk(WALK_COLOURING)
    call w % apply(g, output=f)
    call f % get_integer_vector(c)

    call report(size(c) .eq. 6, "a colouring gives one colour per cell", nfail)
    call report(all(c >= 1), "and every cell gets one", nfail)

    ok = .true.
    do e = 1, g % num_edges()
       t = g % edge_tail(e)
       if (.not. g % edge_has_head(e)) cycle
       h = g % edge_head(e)
       ok = ok .and. c(t) /= c(h)
    end do
    call report(ok, "no face has the same colour at both ends", nfail)

    ! The same must hold on a ring, where a naive alternating colouring
    ! would fail at the join if the count is odd.
    split_mesh = stored_graph(5, tails=[1, 2, 3, 4, 5], heads=[2, 3, 4, 5, 1])
    call w % apply(split_mesh, output=f)
    call f % get_integer_vector(c)
    ok = .true.
    do e = 1, split_mesh % num_edges()
       t = split_mesh % edge_tail(e)
       h = split_mesh % edge_head(e)
       ok = ok .and. c(t) /= c(h)
    end do
    call report(ok, "and still not on an odd ring, where two colours cannot do", nfail)

    ! DEPTH. On a row of six starting at one, the distances are simply
    ! how far along each cell is.
    w = walk(WALK_DEPTH, seed=1)
    call w % apply(g, output=f)
    call f % get_integer_vector(c)
    call report(all(c .eq. [0, 1, 2, 3, 4, 5]), &
         & "depth counts the faces crossed to reach each cell", nfail)

    ! VISIT ORDER. Reached in order along the row.
    w = walk(WALK_VISIT_ORDER, seed=1)
    call w % apply(g, output=f)
    call f % get_integer_vector(c)
    call report(all(c .eq. [1, 2, 3, 4, 5, 6]), &
         & "visit order numbers the cells as the walk reaches them", nfail)

    ! COMPONENT. One connected row is one piece.
    w = walk(WALK_COMPONENT)
    call w % apply(g, output=f)
    call f % get_integer_vector(c)
    call report(all(c .eq. 1), "a connected mesh is one piece", nfail)

    ! A mesh in two halves says so, and the unreachable half is not
    ! given a distance it does not have.
    split_mesh = stored_graph(4, tails=[1, 3], heads=[2, 4])
    call w % apply(split_mesh, output=f)
    call f % get_integer_vector(c)
    call report(c(1) .eq. c(2) .and. c(3) .eq. c(4) .and. c(1) /= c(3), &
         & "a mesh in two halves comes back as two pieces", nfail)

    w = walk(WALK_DEPTH, seed=1)
    call w % apply(split_mesh, output=f)
    call f % get_integer_vector(c)
    call report(c(3) .eq. -1 .and. c(4) .eq. -1, &
         & "and what the seed cannot reach is marked unreachable, not zero", nfail)

  end subroutine check_walks

end program test_graph_contract
