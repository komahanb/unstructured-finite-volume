!=====================================================================!
! A deliberately nonlinear edge formula, defined here to prove a
! contract property: the edge leaf accepts ANY per-edge function of
! the two end values. Nothing in the contract assumes linearity.
!
!             q_t                 q_h
!            (i) ---------------> (j)
!
!            z_e = (q_t + q_h) * |q_h - q_t| / 2
!
! The formula is chosen to be nonlinear in both ends and to have no
! stored coefficients at all: everything state-dependent is computed
! inside apply, from the input, per edge - which is the shape of a
! wave-speed flux formula.
!=====================================================================!

module nonlinear_sample_support

  use iso_fortran_env    , only : dp => REAL64
  use operation_action, only : graph_operation
  use graph_directed_view, only : directed_graph
  use graph_field_calculus, only : graph_field
  use graph_directed_view, only : GRAPH_SIDE_VERTEX
  ! An operation names a domain and counts it. It asks no membership,
  ! so it imports no map: identity and count is the whole of what an
  ! operation is entitled to see.
  use graph_fractal      , only : set_graph => graph
  use class_graph_field  , only : field

  use relation_partition, only : partition_relation
  implicit none
  type(partition_relation) :: rel

  private
  public :: nonlinear_sample

  type, extends(graph_operation) :: nonlinear_sample
   contains
     procedure :: name   => nonlinear_sample_name
     procedure :: domain => nonlinear_sample_domain
     procedure :: apply  => nonlinear_sample_apply
  end type nonlinear_sample

contains

  pure function nonlinear_sample_name(this) result(name)
    class(nonlinear_sample), intent(in) :: this
    character(len=:), allocatable       :: name
    name = 'nonlinear sample'
  end function nonlinear_sample_name

  subroutine nonlinear_sample_domain(this, input_graph, domain, nentries)
    class(nonlinear_sample), intent(in)    :: this
    class(directed_graph), intent(in)               :: input_graph
    type(set_graph), intent(out) :: domain
    integer        , intent(out) :: nentries
    domain   = input_graph % all_edges()
    nentries = input_graph % num_edges()
  end subroutine nonlinear_sample_domain

  subroutine nonlinear_sample_apply(this, input_graph, input_data, output)

    class(nonlinear_sample), intent(in)                 :: this
    class(directed_graph), intent(in)                            :: input_graph
    class(graph_field), intent(in), optional             :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output

    type(field)      :: out
    type(set_graph)    :: on
    real(dp), allocatable :: q(:), z(:)
    real(dp)              :: qt, qh
    integer               :: ne, e, t, h

    ne = input_graph % num_edges()
    on = input_graph % edge_set()
    out = field(this % name(), on, ne)

    allocate(z(ne))
    z = 0.0_dp

    if (present(input_data)) then
       select type (state => input_data(1))
       class is (field)
          call state % get_real_vector(q)
          do e = 1, ne
             t  = input_graph % edge_tail(e)
             qt = q(t)
             if (input_graph % edge_has_head(e)) then
                h  = input_graph % edge_head(e)
                qh = q(h)
             else
                qh = 0.0_dp
             end if
             ! Nonlinear in both ends, computed on the spot.
             z(e) = 0.5_dp * (qt + qh) * abs(qh - qt)
          end do
       end select
    end if

    call out % set_real_vector(z)
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine nonlinear_sample_apply

end module nonlinear_sample_support

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
!      owned cells of each piece are counted in.
!  19. Coarsening and refinement: six cells glued in pairs make three
!      blocks, faces inside a block vanish, values average onto the
!      blocks or add when told to, and refinement holds them back
!      down to every child.
!  20. Edge contributions reduced through incidence exactly once. On a
!      closed ring with no walls the balance must sum to zero, however
!      lopsided the state and however many edge terms there are.
!      Count an edge twice and the sum is wrong by that edge; miss one
!      and it is wrong the other way. Neither crashes. Open the ring
!      and the sum stops cancelling, which is the negative control -
!      without it the check would also pass on a balance that returned
!      nothing but zeros.
!  21. A concrete graph mints its two carrier identities once, at
!      construction, and hands out the SAME identity at every asking.
!      Its two sides are two identities, and a graph built from
!      identical numbers owns its own. This is the law the whole
!      cutover rests on: a domain is WHICH set, and the graph is the
!      only thing that says which. Inherited from the retired carrier
!      suite, whose other laws now live in graph-set-view,
!      graph-inclusion and fractal-graph.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_graph_contract

  use iso_fortran_env       , only : dp => REAL64
  use graph_fractal        , only : set_graph => graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set        , only : set_map
  use map_label      , only : label_map
  use map_inclusion  , only : inclusion_map, declared_subobject
  use graph_field_calculus  , only : GRAPH_FIELD_INTEGER, GRAPH_FIELD_REAL
  use graph_field_calculus  , only : GRAPH_FIELD_COMPLEX, GRAPH_FIELD_LOGICAL
  use graph_field_calculus  , only : GRAPH_FIELD_CHARACTER
  use graph_directed_view   , only : GRAPH_SIDE_VERTEX
  use graph_field_calculus  , only : graph_functional
  use class_graph           , only : directed_stored_graph
  use class_graph_field     , only : field
  use class_graph_functional, only : functional
  use class_graph_reduction , only : reduction
  use class_graph_reduction , only : REDUCE_SUM, REDUCE_AVERAGE, REDUCE_MINIMUM
  use class_graph_reduction , only : REDUCE_MAXIMUM, REDUCE_NORM, REDUCE_COUNT
  use class_graph_reduction , only : REDUCE_ALL, REDUCE_ANY
  use class_graph_reduction , only : broadcast, BROADCAST_COPY, BROADCAST_SHARE
  use class_graph_partitioner, only : partitioner, PARTITION_LINEAR
  use class_graph_partitioner, only : PARTITION_BREADTH_FIRST, PARTITION_ADOPTED
  use class_graph_assembler , only : assembler
  use graph_directed_view   , only : directed_graph
  use graph_field_calculus  , only : graph_field
  use class_graph_coarsener , only : coarsener, COARSEN_PAIRWISE, COARSEN_ADOPTED
  use class_graph_refiner   , only : refiner
  use class_graph_differential_operator, only : differential_operator
  use class_graph_differential_operator, only : edge_differential_operator
  use class_graph_differential_operator, only : vertex_differential_operator
  use class_graph_differential_operator, only : gradient, interpolation
  use class_graph_differential_operator, only : divergence, laplacian
  use class_graph_differential_operator, only : stencil_of
  use class_graph_stencil, only : stencil_operator
  use class_graph_balance   , only : balance
  use class_graph_walk      , only : walk, WALK_COLOURING, WALK_VISIT_ORDER
  use class_graph_walk      , only : WALK_COMPONENT, WALK_DEPTH
  use nonlinear_sample_support, only : nonlinear_sample

  use relation_partition, only : partition_relation
  implicit none
  type(partition_relation) :: rel

  integer :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "graph contract suite"
  write(*,'(1x,a)') "============================================="

  call check_graph_mints_its_carriers(nfail)
  call check_supports(nfail)
  call check_ordering_law(nfail)
  call check_value_kind_rule(nfail)
  call check_field_round_trips(nfail)
  call check_functional_round_trips(nfail)
  call check_graph_structure(nfail)
  call check_graph_named_sets(nfail)
  call check_graph_walking(nfail)
  call check_graph_uncut(nfail)
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
  call check_differential_operators(nfail)
  call check_adjoints(nfail)
  call check_curl_on_border_graph(nfail)
  call check_components(nfail)
  call check_adjoint_is_the_transpose(nfail)
  call check_compiled_stencil_agrees(nfail)
  call check_nonlinear_edge_formulas(nfail)
  call check_inner_products(nfail)
  call check_broadcast(nfail)

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
  ! A concrete graph mints its own carrier identities and says how
  ! many members each holds, but it binds no representation - it owns
  ! no map. So the scope that holds the map describes them, once, at
  ! the point the graph enters that scope: both carriers are {1..n}.
  !
  ! Idempotent on purpose. A graph handed round several checks is
  ! still one identity, and an identity is described once.
  !===================================================================!

  subroutine describe(sets, g)

    type(set_map), intent(inout) :: sets
    class(directed_graph) , intent(in)    :: g

    if (.not. sets % describes(g % vertex_set())) &
         & call sets % bind(g % vertex_set(), &
         &      counted_set_representation(g % num_vertices()))
    if (.not. sets % describes(g % edge_set())) &
         & call sets % bind(g % edge_set(), &
         &      counted_set_representation(g % num_edges()))

  end subroutine describe

  !===================================================================!
  ! A concrete graph mints its two carrier identities at construction
  ! and hands out the SAME one at every asking. Nothing else in the
  ! library can say WHICH set a field or an operation is over, so if
  ! this ever answers a fresh identity per call, every question asked
  ! of a map afterwards is asked of the wrong set - and the map, being
  ! keyed on identity, would answer "nobody describes that" rather
  ! than answer wrongly. The map holds no opinion: the graph is the
  ! only source.
  !
  ! Inherited from the retired carrier suite. What it does NOT claim
  ! is as much of the ruling as what it does: the graph binds no
  ! representation and no label, so the count comes from the graph
  ! itself and the size only after a map has been told.
  !===================================================================!

  subroutine check_graph_mints_its_carriers(nfail)

    integer, intent(inout) :: nfail

    type(directed_stored_graph) :: g, h
    type(set_graph)    :: vs, es, vs2
    type(set_map)      :: sets

    g = directed_stored_graph(3, tails=[1, 2], heads=[2, 3])
    h = directed_stored_graph(3, tails=[1, 2], heads=[2, 3])

    vs = g % vertex_set()
    es = g % edge_set()

    call report(.not. vs % same_as(es), &
         & "a graph's two sides are two sets", nfail)

    vs2 = g % vertex_set()
    call report(vs % same_as(vs2), &
         & "asked twice, the graph answers the same set", nfail)

    call report(.not. vs % same_as(h % vertex_set()), &
         & "an identical twin graph still owns its own sets", nfail)

    ! The count is the graph's own answer, and a map agrees with it
    ! only once it has been told - the graph binds nothing.
    call describe(sets, g)
    call report(sets % size_of(vs) .eq. g % num_vertices() .and. &
         &      sets % size_of(es) .eq. g % num_edges(), &
         & "and once described, the map counts what the graph counts", nfail)

  end subroutine check_graph_mints_its_carriers

  !===================================================================!
  ! A support is a set of indices. It returns what it was given,
  ! reports its size, and reports which kind of index it holds.
  !===================================================================!

  subroutine check_supports(nfail)

    ! REWRITTEN at phase 5B: the support is a subobject S c--> A,
    ! itself a member set - the edgeless-graph fiction is retired.

    integer, intent(inout) :: nfail

    type(set_graph)    :: host
    type(set_graph)     :: vs, empty
    integer, allocatable :: got(:)
    type(set_map)     :: sets
    type(inclusion_map)     :: inclusions

    call host % declare()
    call sets % bind(host, counted_set_representation(12))
    call vs % declare()
    call sets       % bind(vs, listed_set_representation([7, 3, 11]))
    call inclusions % include_in(vs, host)

    call sets % members_of(vs, got)
    call report(size(got) .eq. 3, "a subset keeps its count", nfail)
    call report(all(got .eq. [7, 3, 11]), &
         & "and its members, in declaration order", nfail)
    call report(declared_subobject(vs, host, inclusions), &
         & "and stands embedded in its ambient", nfail)
    call report(sets % has_in(vs, 11) .and. .not. sets % has_in(vs, 5), &
         & "membership answers the chosen family alone", nfail)

    call empty % declare()
    call sets       % bind(empty, listed_set_representation([integer ::]))
    call inclusions % include_in(empty, host)
    call report(sets % size_of(empty) .eq. 0, &
         & "the empty subset is a domain", nfail)

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

    type(set_graph)  :: on
    type(field)    :: f
    real(dp), allocatable :: got(:)
    integer               :: entry_position, component, position
    logical               :: ok
    type(set_map)     :: sets

    call on % declare()
    call sets % bind(on, counted_set_representation(3))
    f  = field('q', on, sets % size_of(on), ncomp=2)

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

    type(set_graph)  :: on
    type(field)    :: f
    real(dp), allocatable :: r(:)
    integer , allocatable :: i(:)
    type(set_map)     :: sets

    call on % declare()
    call sets % bind(on, counted_set_representation(2))
    f  = field('q', on, sets % size_of(on))

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

    type(set_graph)                :: on
    type(field)                  :: f
    integer         , allocatable     :: i(:)
    real(dp)        , allocatable     :: r(:)
    complex(dp)     , allocatable     :: c(:)
    logical         , allocatable     :: l(:)
    character(len=:), allocatable     :: s(:)
    type(set_map)     :: sets

    call on % declare()
    call sets % bind(on, counted_set_representation(2))
    f  = field('F', on, sets % size_of(on))

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
    call value_integer(j, i)
    call report(i .eq. 7, "a functional holds an integer", nfail)

    call j % set_real_value(2.5_dp)
    call value_real(j, r)
    call report(abs(r - 2.5_dp) < 1.0d-13, "a functional holds a real", nfail)

    ! The complex-step case. The derivative is the imaginary part, and
    ! it is tiny on purpose: a functional that quietly rounded through
    ! a real would return exactly zero here.
    call j % set_complex_value((5.0_dp, 4.0d-20))
    call value_complex(j, c)
    call report(abs(real(c) - 5.0_dp) < 1.0d-13, &
         & "a functional holds the real part", nfail)
    call report(abs(aimag(c) - 4.0d-20) < 1.0d-33, &
         & "a functional holds the imaginary part - complex step lives", nfail)
    call report(j % value_kind() .eq. GRAPH_FIELD_COMPLEX, &
         & "and reports itself complex", nfail)

    ! A predicate answers true or false, not one or zero.
    call j % set_logical_value(.false.)
    call value_logical(j, l)
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

  type(directed_stored_graph) function diamond() result(g)

    g = directed_stored_graph(4, tails=[1, 1, 2, 3, 4], heads=[2, 3, 4, 4, 0], &
         &           etags=['none', 'none', 'none', 'none', 'wall'])

  end function diamond

  !===================================================================!
  ! Counts, and where each edge goes.
  !===================================================================!

  subroutine check_graph_structure(nfail)

    integer, intent(inout) :: nfail

    type(directed_stored_graph) :: g
    type(set_map)     :: sets

    g = diamond()
    call describe(sets, g)

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

    type(directed_stored_graph)                       :: g
    type(set_graph) :: vs
    type(set_graph) :: es
    integer, allocatable                     :: indices(:)
    type(set_map)     :: sets
    type(label_map)     :: labels
    type(inclusion_map)     :: inclusions

    g = diamond()
    call describe(sets, g)

    vs = g % all_vertices()
    call members_of(sets, vs, indices)
    call report(size(indices) .eq. 4 .and. all(indices .eq. [1, 2, 3, 4]), &
         & "all_vertices is every cell", nfail)

    call g % boundary_vertices(sets, labels, inclusions, vs)
    call members_of(sets, vs, indices)
    call report(size(indices) .eq. 1 .and. indices(1) .eq. 4, &
         & "only the cell behind the wall is a boundary cell", nfail)

    call g % interior_vertices(sets, labels, inclusions, vs)
    call members_of(sets, vs, indices)
    call report(size(indices) .eq. 3 .and. all(indices .eq. [1, 2, 3]), &
         & "the other three are interior", nfail)

    es = g % all_edges()
    call members_of(sets, es, indices)
    call report(size(indices) .eq. 5, "all_edges is every face", nfail)

    call g % boundary_edges(sets, labels, inclusions, es)
    call members_of(sets, es, indices)
    call report(size(indices) .eq. 1 .and. indices(1) .eq. 5, &
         & "the headless face is the boundary face", nfail)

    call g % interior_edges(sets, labels, inclusions, es)
    call members_of(sets, es, indices)
    call report(size(indices) .eq. 4, "the other four are interior faces", nfail)

    call g % tagged_edges('wall', sets, labels, inclusions, es)
    call members_of(sets, es, indices)
    call report(size(indices) .eq. 1 .and. indices(1) .eq. 5, &
         & "tagged_edges of wall returns the wall face", nfail)

    call g % tagged_edges('inlet', sets, labels, inclusions, es)
    call members_of(sets, es, indices)
    call report(size(indices) .eq. 0, &
         & "an unknown tag returns an empty set, not a failure", nfail)

    ! Nothing was tagged on the vertices, so every vertex query by
    ! name must come back empty rather than reaching into thin air.
    call g % tagged_vertices('heater', sets, labels, inclusions, vs)
    call members_of(sets, vs, indices)
    call report(size(indices) .eq. 0, &
         & "an untagged graph returns an empty set for every vertex tag", nfail)

  end subroutine check_graph_named_sets

  !===================================================================!
  ! Walking the graph, with and without regard to direction.
  !===================================================================!

  subroutine check_graph_walking(nfail)

    integer, intent(inout) :: nfail

    type(directed_stored_graph)   :: g
    integer, allocatable :: indices(:)
    type(set_map)     :: sets

    g = diamond()
    call describe(sets, g)

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

    type(directed_stored_graph)                       :: g
    type(set_graph) :: vs
    integer, allocatable                     :: indices(:)
    type(set_map)     :: sets
    type(label_map)     :: labels
    type(inclusion_map)     :: inclusions
    type(partition_relation) :: grel

    g = diamond()
    call describe(sets, g)

    ! The relation is a VALUE the graph carries, not a question it
    ! answers. A graph born whole stands in the identity relation, and
    ! that is what the next five checks read.
    grel = g % whole_relation()

    call report(grel % describes(g), &
         & "a graph's own relation describes it", nfail)
    call report(.not. grel % has_part_relation(), &
         & "an uncut graph reports no partition record", nfail)
    call report(grel % num_parts() .eq. 1, "an uncut graph is one part", nfail)

    call g % owned_vertices(1, sets, labels, inclusions, vs)
    call members_of(sets, vs, indices)
    call report(size(indices) .eq. 4, "and it owns every cell in it", nfail)

    call g % borrowed_vertices(1, sets, labels, inclusions, vs)
    call members_of(sets, vs, indices)
    call report(size(indices) .eq. 0, "and borrows none", nfail)

    call report(grel % vertex_owner_part(3) .eq. 1, &
         & "every cell belongs to the one part", nfail)
    call report(grel % global_vertex_index(3) .eq. 3, &
         & "its own numbering is the whole-graph numbering", nfail)
    call report(grel % part_vertex_index(3, 1) .eq. 3, &
         & "and the map reads the same backwards", nfail)

  end subroutine check_graph_uncut

  !===================================================================!
  ! Geometry rides on the graph and is fetched by name. This is how a
  ! edge operation reaches a face normal without any caller threading
  ! it down through every call.

  !===================================================================!
  ! The reduction rules, each on numbers whose answer is obvious by
  ! hand. A measure turns a bare sum into an integral.
  !===================================================================!

  subroutine check_reductions(nfail)

    integer, intent(inout) :: nfail

    type(directed_stored_graph)                   :: g
    type(set_graph)                 :: on
    type(field)                   :: f, vol
    type(reduction)                      :: rule
    class(graph_functional), allocatable :: j
    real(dp)                             :: r
    complex(dp)                          :: c
    logical                              :: l
    integer                              :: i
    type(set_map)     :: sets
    type(set_graph) :: pair_pair
    type(set_graph) :: pair_quad

    g = diamond()
    call describe(sets, g)
    call on % declare()
    call sets % bind(on, counted_set_representation(4))

    f = field('q', on, sets % size_of(on))
    call f % set_real_vector([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp])

    rule = reduction(REDUCE_SUM)
    call rule % reduce(f, j)
    call value_real(j, r)
    call report(abs(r - 10.0_dp) < 1.0d-13, "sum adds the values up", nfail)

    ! With a measure the same rule becomes an integral. Each cell is
    ! weighted by its volume, so the answer stops depending on how
    ! finely the mesh was cut.
    vol = field('cell_volume', on, sets % size_of(on))
    call vol % set_real_vector([2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp])
    call rule % reduce(f, j, measure=vol)
    call value_real(j, r)
    call report(abs(r - 20.0_dp) < 1.0d-13, &
         & "a measure turns the sum into an integral", nfail)

    rule = reduction(REDUCE_MINIMUM)
    call rule % reduce(f, j)
    call value_real(j, r)
    call report(abs(r - 1.0_dp) < 1.0d-13, "minimum finds the smallest", nfail)

    rule = reduction(REDUCE_MAXIMUM)
    call rule % reduce(f, j)
    call value_real(j, r)
    call report(abs(r - 4.0_dp) < 1.0d-13, "maximum finds the largest", nfail)

    ! Three-four-five, so the root is exact.
    ! The shape law now refuses reusing a four-cell field for two
    ! values: the norm check gets its own two-entry domain.
    call pair_pair % declare()
    call sets % bind(pair_pair, counted_set_representation(2))
    f = field('q', pair_pair, 2)
    call f % set_real_vector([3.0_dp, 4.0_dp])
    rule = reduction(REDUCE_NORM)
    call rule % reduce(f, j)
    call value_real(j, r)
    call report(abs(r - 5.0_dp) < 1.0d-13, "the two-norm takes its root once", nfail)

    call pair_quad % declare()
    call sets % bind(pair_quad, counted_set_representation(4))
    f = field('q', pair_quad, 4)
    call f % set_real_vector([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp])
    rule = reduction(REDUCE_COUNT)
    call rule % reduce(f, j)
    call value_integer(j, i)
    call report(i .eq. 4, "count counts, and answers as an integer", nfail)
    call report(j % value_kind() .eq. GRAPH_FIELD_INTEGER, &
         & "and reports integer kind", nfail)

    ! A predicate over a logical field. This is the shape a question
    ! like "is this graph acyclic" comes back in.
    call f % set_logical_vector([.true., .true., .false., .true.])
    rule = reduction(REDUCE_ALL)
    call rule % reduce(f, j)
    call value_logical(j, l)
    call report(.not. l, "all is false when one of them is", nfail)

    rule = reduction(REDUCE_ANY)
    call rule % reduce(f, j)
    call value_logical(j, l)
    call report(l, "any is true when one of them is", nfail)
    call report(j % value_kind() .eq. GRAPH_FIELD_LOGICAL, &
         & "and a predicate stays a predicate, not a one or a zero", nfail)

    ! The complex-step road. The imaginary parts are tiny on purpose:
    ! a reduction that rounded through a real would return zero.
    ! The two-entry domain is the one already signed above: an identity
    ! is assigned once, and {1,2} is the same set both times.
    f = field('q', pair_pair, 2)
    call f % set_complex_vector([(1.0_dp, 1.0d-20), (2.0_dp, 3.0d-20)])
    rule = reduction(REDUCE_SUM)
    call rule % reduce(f, j)
    call value_complex(j, c)
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

    type(directed_stored_graph)                   :: g
    type(set_graph)                 :: s1, s2
    type(field)                   :: f1, f2
    type(reduction)                      :: rule
    class(graph_functional), allocatable :: a, b, both, j
    real(dp)                             :: r
    type(set_map)     :: sets

    g = diamond()
    call describe(sets, g)
    rule = reduction(REDUCE_AVERAGE)

    ! Two part-sized domains: three cells here, two there.
    call s1 % declare()
    call sets % bind(s1, counted_set_representation(3))
    f1 = field('q', s1, sets % size_of(s1))
    call f1 % set_real_vector([2.0_dp, 2.0_dp, 2.0_dp])

    call s2 % declare()
    call sets % bind(s2, counted_set_representation(2))
    f2 = field('q', s2, sets % size_of(s2))
    call f2 % set_real_vector([5.0_dp, 9.0_dp])

    ! Each part accumulates on its own, and neither one divides.
    call rule % initialize(a)
    call rule % accumulate(f1, a)

    call rule % initialize(b)
    call rule % accumulate(f2, b)

    call rule % combine(a, b, both)
    call rule % finalize(both, j)

    call value_real(j, r)
    call report(abs(r - 4.0_dp) < 1.0d-13, &
         & "two parts average to 4, not 4.5 - the division waits", nfail)

    ! Joining the parts the other way round must give the same answer,
    ! or a parallel run would depend on which image finished first.
    call rule % combine(b, a, both)
    call rule % finalize(both, j)
    call value_real(j, r)
    call report(abs(r - 4.0_dp) < 1.0d-13, &
         & "and the order the parts arrive in makes no difference", nfail)

    ! One part alone still has to come out right.
    call rule % finalize(a, j)
    call value_real(j, r)
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

  type(directed_stored_graph) function chain_of_six() result(g)

    g = directed_stored_graph(6, tails=[1, 2, 3, 4, 5], heads=[2, 3, 4, 5, 6])

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

    type(directed_stored_graph)             :: g
    type(partitioner)              :: p
    type(assembler)                :: a
    class(directed_graph), allocatable      :: part, back
    class(graph_field), allocatable :: pd, fd
    type(set_graph)           :: on
    type(field)             :: d
    real(dp), allocatable          :: v(:)
    integer                        :: e
    logical                        :: same
    type(set_map)     :: sets
    type(label_map)     :: labels
    type(inclusion_map)     :: inclusions

    g = chain_of_six()
    call describe(sets, g)
    p = partitioner(PARTITION_LINEAR, nparts=1, part=1)
    a = assembler()

    call report(p % defined_on_graph(g), "a partitioner accepts a real graph", nfail)

    call p % partition_graph(g, part, rel)
    call describe(sets, part)
    call report(part % num_vertices() .eq. 6, &
         & "cut into one, the piece has every cell", nfail)
    call report(part % num_edges() .eq. 5, "and every face", nfail)

    call report(a % defined_on_relation(rel, part), &
         & "the relation the cut wrote is the one this part stands in", &
         & nfail)

    call a % assemble_graph(rel, part, back)
    call describe(sets, back)
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
    on = g % vertex_set()
    d  = field('q', on, sets % size_of(on))
    call d % set_real_vector([10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 50.0_dp, 60.0_dp])

    call p % partition_data(rel, g, d, part, sets, labels, inclusions, pd)
    call a % assemble_data(rel, part, pd, g, sets, labels, inclusions, fd)

    select type (fd)
    class is (field)
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

    type(directed_stored_graph)                       :: g
    type(partitioner)                        :: p
    class(directed_graph), allocatable                :: part
    type(set_graph) :: vs
    integer, allocatable                     :: indices(:)
    integer                                  :: k, l, f, times(6)
    type(set_map)     :: sets
    type(label_map)     :: labels
    type(inclusion_map)     :: inclusions

    g     = chain_of_six()
    call describe(sets, g)
    times = 0

    do k = 1, 2
       p = partitioner(PARTITION_LINEAR, nparts=2, part=k)
       call p % partition_graph(g, part, rel)
       call describe(sets, part)

       call report(rel % has_part_relation(), &
            & "has_part_relation is true on a cut piece", nfail)
       call report(rel % num_parts() .eq. 2, "and how many pieces there are", nfail)

       ! Count how often each whole-graph cell is owned.
       call part % owned_vertices(k, sets, labels, inclusions, vs)
       call members_of(sets, vs, indices)
       do l = 1, size(indices)
          f = rel % global_vertex_index(indices(l))
          times(f) = times(f) + 1
       end do

       ! Each piece should also have borrowed the one cell across the
       ! cut, so a face term there has a value on both sides.
       call part % borrowed_vertices(k, sets, labels, inclusions, vs)
       call members_of(sets, vs, indices)
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

    type(directed_stored_graph), intent(in) :: g
    integer           , intent(in) :: rule
    integer           , intent(in) :: nparts

    type(partitioner)                        :: p
    class(directed_graph), allocatable                :: part
    type(set_graph) :: vs
    integer, allocatable                     :: indices(:)
    integer                                  :: times(g % num_vertices())
    integer                                  :: k, l, f
    type(set_map)     :: sets
    type(label_map)     :: labels
    type(inclusion_map)     :: inclusions

    times = 0

    do k = 1, nparts
       p = partitioner(rule, nparts=nparts, part=k)
       call p % partition_graph(g, part, rel)
       call describe(sets, part)
       call part % owned_vertices(k, sets, labels, inclusions, vs)
       call members_of(sets, vs, indices)
       do l = 1, size(indices)
          f = rel % global_vertex_index(indices(l))
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

    type(directed_stored_graph)             :: g
    type(partitioner)              :: p
    type(assembler)                :: a
    class(directed_graph), allocatable      :: part
    class(graph_field), allocatable :: pd, fd
    type(set_graph)           :: on
    type(field)             :: d
    real(dp), allocatable          :: v(:)
    real(dp)                       :: total(6)
    integer                        :: k
    type(set_map)     :: sets
    type(label_map)     :: labels
    type(inclusion_map)     :: inclusions

    g = chain_of_six()
    call describe(sets, g)
    a = assembler()

    on = g % vertex_set()
    d  = field('q', on, sets % size_of(on))
    call d % set_real_vector([10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 50.0_dp, 60.0_dp])

    total = 0.0_dp

    do k = 1, 2
       p = partitioner(PARTITION_LINEAR, nparts=2, part=k)
       call p % partition_graph(g, part, rel)
       call describe(sets, part)
       call p % partition_data(rel, g, d, part, sets, labels, inclusions, pd)
       call a % assemble_data(rel, part, pd, g, sets, labels, inclusions, fd)

       select type (fd)
       class is (field)
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
  !      whole_A( G, D )  ==  assemble( A( partition( G, D ) ) )
  !
  ! Working something out on the whole must give the same answer as
  ! working it out on each piece and putting the pieces together. Here
  ! A is a sum - the smallest operation that exercises the law - and
  ! the pieces are joined by the reduction's own combine.
  !
  ! Only the owned cells of each piece are counted in. Count the
  ! borrowed ones too and the shared cell appears twice - which is
  ! the whole reason a reduction has to know about ownership.
  !===================================================================!

  subroutine check_operation_consistency(nfail)

    integer, intent(inout) :: nfail

    type(directed_stored_graph)                       :: g
    type(partitioner)                        :: p
    type(reduction)                          :: rule
    class(directed_graph), allocatable                :: part
    class(graph_field), allocatable           :: pd
    type(set_graph) :: vs
    class(graph_functional), allocatable     :: whole, piece, running, joined
    type(set_graph)                        :: on
    ! Two parts, two owned sets: an identity is assigned once, so the
    ! loop cannot re-sign one name. They are different sets anyway.
    type(set_graph)                         :: owned_only(2)
    type(field)                       :: d, owned_values
    real(dp), allocatable                    :: v(:), pick(:)
    integer, allocatable                     :: indices(:)
    real(dp)                                 :: a_whole, a_parts
    integer                                  :: k, l
    type(set_map)     :: sets
    type(label_map)     :: labels
    type(inclusion_map)     :: inclusions

    g    = chain_of_six()
    call describe(sets, g)
    rule = reduction(REDUCE_SUM)

    on = g % vertex_set()
    d  = field('q', on, sets % size_of(on))
    call d % set_real_vector([10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 50.0_dp, 60.0_dp])

    ! A on the whole.
    call rule % reduce(d, whole)
    call value_real(whole, a_whole)
    call report(abs(a_whole - 210.0_dp) < 1.0d-13, &
         & "the sum over the whole graph is 210", nfail)

    ! A on each piece, then joined.
    call rule % initialize(running)

    do k = 1, 2
       p = partitioner(PARTITION_LINEAR, nparts=2, part=k)
       call p % partition_graph(g, part, rel)
       call describe(sets, part)
       call p % partition_data(rel, g, d, part, sets, labels, inclusions, pd)

       ! Accumulate this piece's owned cells only.
       call part % owned_vertices(k, sets, labels, inclusions, vs)
       call members_of(sets, vs, indices)

       select type (pd)
       class is (field)
          call pd % get_real_vector(v)
          allocate(pick(size(indices)))
          do l = 1, size(indices)
             pick(l) = v(indices(l))
          end do
       end select

       call owned_only(k) % declare()
       call sets       % bind(owned_only(k), listed_set_representation(indices))
       call inclusions % include_in(owned_only(k), part % vertex_set())
       owned_values = field('q', owned_only(k), sets % size_of(owned_only(k)))
       call owned_values % set_real_vector(pick)
       deallocate(pick)

       call rule % initialize(piece)
       call rule % accumulate(owned_values, piece)
       call rule % combine(running, piece, joined)
       call rule % initialize(running)
       call rule % combine(joined, running, running)
    end do

    call value_real(running, a_parts)
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

    type(directed_stored_graph)             :: g
    type(coarsener)                :: c
    type(refiner)                  :: r
    class(directed_graph), allocatable      :: coarse, fine
    class(graph_field), allocatable :: cd, fd
    type(set_graph)           :: on
    type(field)             :: d
    real(dp), allocatable          :: v(:)
    type(set_map)     :: sets

    g = chain_of_six()
    call describe(sets, g)
    c = coarsener(COARSEN_PAIRWISE)

    call report(c % defined_on_graph(g), "a coarsener accepts a graph with room to shrink", nfail)

    call c % coarsen_graph(g, coarse)
    call describe(sets, coarse)
    call report(coarse % num_vertices() .eq. 3, &
         & "six cells glued in pairs make three blocks", nfail)
    call report(coarse % num_edges() .lt. g % num_edges(), &
         & "and some faces vanish inside the blocks", nfail)

    ! Averaging onto the blocks. Cells 1 and 2 hold 10 and 20, so
    ! their block holds 15.
    on = g % vertex_set()
    d  = field('q', on, sets % size_of(on))
    call d % set_real_vector([10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 50.0_dp, 60.0_dp])

    call c % coarsen_data(g, d, coarse, cd)
    select type (cd)
    class is (field)
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
    class is (field)
       call cd % get_real_vector(v)
       call report(abs(v(1) - 30.0_dp) < 1.0d-13, &
            & "told to add rather than average, a block holds the total", nfail)
    end select

    ! And back down. Every child starts from its parent's value.
    r = refiner(2)
    call r % refine_graph(coarse, fine)
    call describe(sets, fine)
    call report(fine % num_vertices() .eq. 6, &
         & "splitting three blocks in two gives six cells back", nfail)

    call c % coarsen_data(g, d, coarse, cd)
    call r % refine_data(coarse, cd, fine, fd)
    select type (fd)
    class is (field)
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

    type(directed_stored_graph)             :: one_cell, g
    type(coarsener)                :: c
    type(assembler)                :: a
    type(partitioner)              :: p
    class(directed_graph), allocatable      :: out, part, back
    integer, allocatable           :: no_edges(:)
    type(set_map)     :: sets
    type(partition_relation)       :: grel

    allocate(no_edges(0))
    one_cell = directed_stored_graph(1, tails=no_edges, heads=no_edges)
    call describe(sets, one_cell)

    ! A single cell is as coarse as a graph gets.
    c = coarsener(COARSEN_PAIRWISE)
    call report(c % defined_on_graph(one_cell), &
         & "a coarsener accepts a graph already at one cell", nfail)

    call c % coarsen_graph(one_cell, out)
    call describe(sets, out)
    call report(out % num_vertices() .eq. 1, &
         & "it returns the same single cell", nfail)

    ! A graph that was never cut is already assembled - and it says so
    ! through the IDENTITY RELATION it carries, which is handed to the
    ! assembler like any other. defined_on_relation is what makes the
    ! pairing checkable: another graph's relation is refused here.
    g = chain_of_six()
    call describe(sets, g)
    a    = assembler()
    grel = g % whole_relation()
    call report(a % defined_on_relation(grel, g), &
         & "an assembler accepts a graph that was never cut", nfail)

    call a % assemble_graph(grel, g, back)
    call describe(sets, back)
    call report(back % num_vertices() .eq. 6 .and. back % num_edges() .eq. 5, &
         & "it returns the same graph", nfail)

    ! And cutting into one piece is the same story from the other end.
    p = partitioner(PARTITION_LINEAR, nparts=1, part=1)
    call p % partition_graph(g, part, rel)
    call describe(sets, part)
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

    type(directed_stored_graph)                       :: ring, open_chain
    type(differential_operator)         :: edge_term
    type(balance)                            :: bal
    class(graph_field), allocatable   :: y
    type(set_graph)                     :: on
    type(field)                       :: q
    real(dp), allocatable                    :: v(:)
    real(dp)                                 :: total
    logical                                  :: ok
    type(set_map)     :: sets

    ! A closed ring: four cells, four faces, not a wall in sight.
    ring = directed_stored_graph(4, tails=[1, 2, 3, 4], heads=[2, 3, 4, 1])
    call describe(sets, ring)

    edge_term = gradient(coefficient=1.0_dp)
    bal       = balance(edge_terms=[edge_term])

    on = ring % vertex_set()
    q  = field('q', on, sets % size_of(on))

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
         & "on a closed ring the balance sums to zero - every edge once", nfail)

    ! Two edge terms, not one. Each is reduced through incidence in
    ! its turn, so the sum still has to cancel.
    bal = balance(edge_terms=[gradient(coefficient=1.0_dp), &
         &                    interpolation(coefficient=0.5_dp)])
    call q % set_real_vector([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp])
    call bal % apply(ring, [q], y)
    call y % get_real_vector(v)
    call report(abs(sum(v)) < 1.0d-12, &
         & "two edge terms reduce once each, and still cancel", nfail)

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
    open_chain = directed_stored_graph(4, tails=[1, 2, 3, 4], heads=[2, 3, 4, 0])
    call describe(sets, open_chain)
    ! Domains are identities now: the state must ride the graph it
    ! is applied against.
    q = field('q', open_chain % vertex_set(), open_chain % num_vertices())
    call q % set_real_vector([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp])
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

    type(directed_stored_graph)                     :: g, split_mesh
    type(walk)                             :: w
    class(graph_field), allocatable :: f
    integer, allocatable                   :: c(:)
    integer                                :: e, t, h
    logical                                :: ok
    type(set_map)     :: sets

    g = chain_of_six()
    call describe(sets, g)

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
    split_mesh = directed_stored_graph(5, tails=[1, 2, 3, 4, 5], heads=[2, 3, 4, 5, 1])
    call describe(sets, split_mesh)
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
    split_mesh = directed_stored_graph(4, tails=[1, 3], heads=[2, 4])
    call describe(sets, split_mesh)
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

  !===================================================================!
  ! The differential operators, checked against the numbers the
  ! theory guide proves.
  !
  !      (1)--(2)--(3)--(4)--(5)--(6)--(7)      spacing 1, measure 1
  !===================================================================!

  subroutine check_differential_operators(nfail)

    integer, intent(inout) :: nfail

    type(directed_stored_graph)                     :: g7, g3, ring
    type(differential_operator)     :: op
    class(graph_field), allocatable :: yf
    type(set_graph)                   :: on
    type(field)                     :: q
    type(set_graph)                     :: eon
    type(field)                       :: zf
    real(dp), allocatable                  :: y(:)
    real(dp)                               :: qv(7)
    integer                                :: v
    logical                                :: ok
    type(set_map)     :: sets

    g7 = directed_stored_graph(7, tails=[1,2,3,4,5,6], heads=[2,3,4,5,6,7])
    call describe(sets, g7)
    on = g7 % vertex_set()
    q  = field('q', on, sets % size_of(on))

    ! Order 0 is the term itself, coefficient and all.
    do v = 1, 7
       qv(v) = real(v, dp)
    end do
    call q % set_real_vector(qv)
    op = vertex_differential_operator(order=0, coefficient=3.0_dp)
    call op % apply(g7, [q], yf)
    call yf % get_real_vector(y)
    call report(all(abs(y - 3.0_dp * qv) < 1.0d-12), &
         & "order 0 returns c q at every vertex", nfail)

    ! Order 1 of a straight line is c times the slope, both signs.
    call q % set_real_vector(10.0_dp * qv)
    op = vertex_differential_operator(order=1, coefficient=1.0_dp)
    call op % apply(g7, [q], yf)
    call yf % get_real_vector(y)
    ok = .true.
    do v = 2, 6
       ok = ok .and. abs(y(v) - 10.0_dp) < 1.0d-12
    end do
    call report(ok, "order 1 of a straight line is the slope, inside", nfail)

    op = vertex_differential_operator(order=1, coefficient=-1.0_dp)
    call op % apply(g7, [q], yf)
    call yf % get_real_vector(y)
    ok = .true.
    do v = 2, 6
       ok = ok .and. abs(y(v) + 10.0_dp) < 1.0d-12
    end do
    call report(ok, "and with c negative, minus the slope - the other end", nfail)

    ! Order 2: zero on the line, two on the parabola. The guide's
    ! first claim, with its numbers.
    op = laplacian()
    call op % apply(g7, [q], yf)
    call yf % get_real_vector(y)
    ok = .true.
    do v = 2, 6
       ok = ok .and. abs(y(v)) < 1.0d-12
    end do
    call report(ok, "the second derivative of a straight line is zero", nfail)

    do v = 1, 7
       qv(v) = real(v, dp)**2
    end do
    call q % set_real_vector(qv)
    call op % apply(g7, [q], yf)
    call yf % get_real_vector(y)
    ok = .true.
    do v = 2, 6
       ok = ok .and. abs(y(v) - 2.0_dp) < 1.0d-12
    end do
    call report(ok, "and of a parabola, exactly two", nfail)

    ! The first vertex has one edge only, so its value is one-sided:
    ! the slope of the first edge, q(2) - q(1) = 3.
    call report(abs(y(1) - 3.0_dp) < 1.0d-12, &
         & "the chain's first vertex holds the one-sided value", nfail)

    ! Order 4: zero on the parabola, twenty-four on the fourth power,
    ! two rings inside. The guide's second claim.
    op = vertex_differential_operator(order=4)
    call op % apply(g7, [q], yf)
    call yf % get_real_vector(y)
    ok = .true.
    do v = 3, 5
       ok = ok .and. abs(y(v)) < 1.0d-12
    end do
    call report(ok, "the fourth derivative of a parabola is zero", nfail)

    do v = 1, 7
       qv(v) = real(v, dp)**4
    end do
    call q % set_real_vector(qv)
    call op % apply(g7, [q], yf)
    call yf % get_real_vector(y)
    ok = .true.
    do v = 3, 5
       ok = ok .and. abs(y(v) - 24.0_dp) < 1.0d-12
    end do
    call report(ok, "and of the fourth power, exactly twenty-four", nfail)

    ! Per-edge coefficients: three cells, weights 2 and 5, values
    ! 1, 2, 4. At the middle: 5*(4-2) - 2*(2-1) = 10 - 2 = 8.
    g3 = directed_stored_graph(3, tails=[1, 2], heads=[2, 3])
    call describe(sets, g3)
    q  = field('q', g3 % vertex_set(), g3 % num_vertices())
    call q % set_real_vector([1.0_dp, 2.0_dp, 4.0_dp])
    op = laplacian(coefficients=[2.0_dp, 5.0_dp])
    call op % apply(g3, [q], yf)
    call yf % get_real_vector(y)
    call report(abs(y(2) - 8.0_dp) < 1.0d-12, &
         & "per-edge coefficients: the middle vertex sums to eight", nfail)

    ! Divergence of a given edge field: a ring with samples 1,2,3,4.
    ! Out minus in at each vertex: -3, 1, 1, 1.
    ring = directed_stored_graph(4, tails=[1,2,3,4], heads=[2,3,4,1])
    call describe(sets, ring)
    eon = ring % edge_set()
    zf   = field('z', eon, sets % size_of(eon))
    call zf % set_real_vector([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp])
    op = divergence()
    call op % apply(ring, [zf], yf)
    call yf % get_real_vector(y)
    call report(abs(y(1) + 3.0_dp) < 1.0d-12 .and. &
         &      all(abs(y(2:4) - 1.0_dp) < 1.0d-12), &
         & "an edge field at order one is its incidence sum, hand-checked", nfail)

    ! An edge field at order two, with per-edge coefficients: the
    ! coefficient scales the given samples once, before the first
    ! incidence step, and never again in the deeper chain. By hand,
    ! with
    ! samples 1,2,3,4 and coefficients 2,1,1,1 around the ring:
    ! scaled samples 2,2,3,4; first reduction -2,0,1,1; then the order-1
    ! chain on that: -3, 2, 1, 0.
    op = vertex_differential_operator(order=2, &
         & coefficients=[2.0_dp, 1.0_dp, 1.0_dp, 1.0_dp])
    call op % apply(ring, [zf], yf)
    call yf % get_real_vector(y)
    call report(abs(y(1) + 3.0_dp) < 1.0d-12 .and. &
         &      abs(y(2) - 2.0_dp) < 1.0d-12 .and. &
         &      abs(y(3) - 1.0_dp) < 1.0d-12 .and. &
         &      abs(y(4)) < 1.0d-12, &
         & "a per-edge coefficient is spent once, not twice", nfail)

  end subroutine check_differential_operators

  !===================================================================!
  ! The adjoint is the reverse walk. For any q and p the pairing
  ! must balance:
  !
  !      sum of (A q) p  =  sum of q (A* p)
  !
  ! and for the order-2 operator with the same coefficients on both
  ! steps, A* is A itself.
  !===================================================================!

  subroutine check_adjoints(nfail)

    integer, intent(inout) :: nfail

    type(directed_stored_graph)                     :: g7
    type(differential_operator)     :: fwd, rev
    class(graph_field), allocatable :: yf
    type(set_graph)                   :: on
    type(field)                     :: qf, pf
    real(dp), allocatable                  :: aq(:), ap(:)
    real(dp)                               :: q(7), p(7), left, right
    integer                                :: v
    type(set_map)     :: sets

    g7 = directed_stored_graph(7, tails=[1,2,3,4,5,6], heads=[2,3,4,5,6,7])
    call describe(sets, g7)
    on = g7 % vertex_set()
    qf = field('q', on, sets % size_of(on))
    pf = field('p', on, sets % size_of(on))

    ! Two fields with nothing special about them.
    do v = 1, 7
       q(v) = real(v, dp)**2 - 3.0_dp * v
       p(v) = 2.0_dp * v + real(7 - v, dp)**2
    end do
    call qf % set_real_vector(q)
    call pf % set_real_vector(p)

    ! The order-1 operator against its reverse walk.
    fwd = vertex_differential_operator(order=1, coefficient=1.5_dp)
    rev = vertex_differential_operator(order=1, coefficient=1.5_dp, adjoint=.true.)

    call fwd % apply(g7, [qf], yf)
    call yf % get_real_vector(aq)
    call rev % apply(g7, [pf], yf)
    call yf % get_real_vector(ap)

    left  = sum(aq * p)
    right = sum(q * ap)
    call report(abs(left - right) < 1.0d-10, &
         & "order 1 pairs with its reverse walk exactly", nfail)

    ! The same identity read backwards: the reverse of the reverse
    ! acts as the original.
    call rev % apply(g7, [qf], yf)
    call yf % get_real_vector(aq)
    call fwd % apply(g7, [pf], yf)
    call yf % get_real_vector(ap)
    call report(abs(sum(aq * p) - sum(q * ap)) < 1.0d-10, &
         & "and the reverse of the reverse acts as the original", nfail)

    ! Order 2 with per-edge coefficients on its steps is its own
    ! adjoint: the two sign flips cancel in pairs.
    fwd = laplacian(coefficients=[2.0_dp, 5.0_dp, 1.0_dp, 4.0_dp, 3.0_dp, 6.0_dp])
    call fwd % apply(g7, [qf], yf)
    call yf % get_real_vector(aq)
    call fwd % apply(g7, [pf], yf)
    call yf % get_real_vector(ap)
    call report(abs(sum(aq * p) - sum(q * ap)) < 1.0d-9, &
         & "the weighted order-2 operator is its own adjoint", nfail)

  end subroutine check_adjoints

  !===================================================================!
  ! Curl, with nothing new. The square and its border graph:
  !
  !        (4) <--c-- (3)          a . . . +1 . .
  !         |          ^           b . . . +1 . . \
  !         d          b                            (f)
  !         v          |           c . . . +1 . . /
  !        (1) --a--> (2)          d . . . +1 . .
  !
  ! Walk the loop and add what the edges hold. A difference field
  ! cancels term against term; a circulating field does not.
  !===================================================================!

  subroutine check_curl_on_border_graph(nfail)

    integer, intent(inout) :: nfail

    type(directed_stored_graph)                     :: border
    type(differential_operator)     :: reduce_edges
    class(graph_field), allocatable :: yf
    type(set_graph)                     :: eon
    type(field)                       :: zf
    real(dp), allocatable                  :: y(:)
    real(dp)                               :: qsq(4)
    type(set_map)     :: sets

    ! The border graph: vertices 1..4 stand for the square's edges
    ! a..d, vertex 5 for the face; every border edge points at the
    ! face, and every square edge agrees with the walk, so every
    ! orientation coefficient is one.
    border = directed_stored_graph(5, tails=[1, 2, 3, 4], heads=[5, 5, 5, 5])
    call describe(sets, border)
    eon = border % edge_set()
    zf     = field('z', eon, sets % size_of(eon))
    reduce_edges = divergence()

    ! A difference field around the square: values q at the corners,
    ! each edge holding the difference of its ends.
    qsq = [3.0_dp, 7.0_dp, -2.0_dp, 5.0_dp]
    call zf % set_real_vector([qsq(2) - qsq(1), qsq(3) - qsq(2), &
         &                     qsq(4) - qsq(3), qsq(1) - qsq(4)])
    call reduce_edges % apply(border, [zf], yf)
    call yf % get_real_vector(y)
    call report(abs(y(5)) < 1.0d-12, &
         & "the curl of a difference field is exactly zero", nfail)

    ! A circulating field: the same value all the way around. The
    ! walk collects it four times, and the answer must not be zero -
    ! a check that cannot fail proves nothing.
    call zf % set_real_vector([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp])
    call reduce_edges % apply(border, [zf], yf)
    call yf % get_real_vector(y)
    call report(abs(y(5)) > 1.0d-12, &
         & "and a circulating field does not cancel", nfail)

  end subroutine check_curl_on_border_graph

  !===================================================================!
  ! Components. One field, two values per vertex - a straight line in
  ! the first slot, a parabola in the second, interleaved by the
  ! ordering rule:
  !
  !      vertex          1        1        2        2
  !      component       1        2        1        2
  !                   +--------+--------+--------+--------+--
  !      flat vector  |  line  |  v^2   |  line  |  v^2   |
  !                   +--------+--------+--------+--------+--
  !
  ! One application of the operator must answer for both at once:
  ! zero for the line, two for the parabola, each in its own slot.
  !===================================================================!

  subroutine check_components(nfail)

    integer, intent(inout) :: nfail

    type(directed_stored_graph)                     :: g7
    type(differential_operator)     :: op, fwd, rev
    class(graph_field), allocatable :: yf
    type(set_graph)                   :: on
    type(field)                     :: qf, pf
    real(dp), allocatable                  :: y(:), aq(:), ap(:)
    real(dp)                               :: q(14), p(14)
    integer                                :: v
    logical                                :: ok
    type(set_map)     :: sets

    g7 = directed_stored_graph(7, tails=[1,2,3,4,5,6], heads=[2,3,4,5,6,7])
    call describe(sets, g7)
    on = g7 % vertex_set()
    qf = field('q', on, sets % size_of(on), ncomp=2)

    do v = 1, 7
       q(2*v - 1) = 10.0_dp * v          ! component one: a line
       q(2*v)     = real(v, dp)**2       ! component two: a parabola
    end do
    call qf % set_real_vector(q)

    op = laplacian()
    call op % apply(g7, [qf], yf)

    call report(yf % num_components() .eq. 2, &
         & "the answer carries as many components as the question", nfail)

    call yf % get_real_vector(y)
    ok = .true.
    do v = 2, 6
       ok = ok .and. abs(y(2*v - 1)) < 1.0d-12          ! the line: zero
       ok = ok .and. abs(y(2*v) - 2.0_dp) < 1.0d-12     ! the parabola: two
    end do
    call report(ok, &
         & "one walk answers both components, each in its own slot", nfail)

    ! The adjoint pairing on the whole interleaved vector: with two
    ! components aboard, the sums must still balance.
    pf = field('p', on, sets % size_of(on), ncomp=2)
    do v = 1, 14
       p(v) = real(v, dp)**2 - 5.0_dp * v
    end do
    call pf % set_real_vector(p)

    fwd = vertex_differential_operator(order=1, coefficient=2.0_dp)
    rev = vertex_differential_operator(order=1, coefficient=2.0_dp, adjoint=.true.)

    call fwd % apply(g7, [qf], yf)
    call yf % get_real_vector(aq)
    call rev % apply(g7, [pf], yf)
    call yf % get_real_vector(ap)
    call report(abs(sum(aq * p) - sum(q * ap)) < 1.0d-9, &
         & "the reverse walk pairs exactly with two components aboard", nfail)

  end subroutine check_components

  !===================================================================!
  ! The compiled form agrees with the operator it compiles from.
  ! stencil_of returns the composed map as a stencil_operator -
  ! same triples, same constants - so applying either must give
  ! the same numbers, on a graph with a boundary edge so the
  ! affine part (the boundary value in the stencil's constants) is
  ! checked too, forward and adjoint.
  !===================================================================!

  subroutine check_compiled_stencil_agrees(nfail)

    integer, intent(inout) :: nfail

    type(directed_stored_graph) :: chain
    type(differential_operator) :: op, rev
    type(stencil_operator)      :: compiled
    type(field)                     :: qf
    class(graph_field), allocatable :: ya, yb
    real(dp), allocatable           :: a(:), b(:)
    type(set_graph)                 :: on
    integer :: v

    ! five vertices, the last edge headless: a boundary
    chain = directed_stored_graph(5, tails=[1,2,3,4,5], heads=[2,3,4,5,0])
    on = chain % vertex_set()

    qf = field('q', on, chain % num_vertices())
    call qf % set_real_vector([(real(v, dp)**2 + 3.0_dp * v, v = 1, 5)])

    op = laplacian(coefficients=[2.0_dp, 5.0_dp, 1.0_dp, 4.0_dp, 3.0_dp], &
         & boundary_value=7.0_dp)

    call op % apply(chain, [qf], ya)
    call ya % get_real_vector(a)

    compiled = stencil_of(op, chain)
    call compiled % apply(chain, [qf], yb)
    call yb % get_real_vector(b)

    call report(size(a) .eq. size(b) .and. maxval(abs(a - b)) .le. 0.0_dp, &
         & "the compiled stencil reproduces the operator, boundary included", &
         & nfail)

    rev = vertex_differential_operator(order=2, &
         & coefficients=[2.0_dp, 5.0_dp, 1.0_dp, 4.0_dp, 3.0_dp], &
         & adjoint=.true.)

    call rev % apply(chain, [qf], ya)
    call ya % get_real_vector(a)

    compiled = stencil_of(rev, chain)
    call compiled % apply(chain, [qf], yb)
    call yb % get_real_vector(b)

    call report(size(a) .eq. size(b) .and. maxval(abs(a - b)) .le. 0.0_dp, &
         & "the compiled adjoint is the transposed stencil", nfail)

  end subroutine check_compiled_stencil_agrees

  !===================================================================!
  ! The strongest statement the adjoint flag can make: the matrix.
  !
  ! Apply the operator to each unit field in turn and its columns
  ! appear; do the same with the flag raised; the second matrix must
  ! be the first one written sideways -
  !
  !      A*(i,j) = A(j,i)   for every entry, at every order,
  !
  ! on a graph that has a boundary edge, so the rows a boundary bends
  ! are checked too. A pairing identity can hold by luck on one pair
  ! of fields; twenty-five entries per order cannot.
  !===================================================================!

  subroutine check_adjoint_is_the_transpose(nfail)

    integer, intent(inout) :: nfail

    type(directed_stored_graph)                     :: g
    type(differential_operator)     :: fwd, rev
    class(graph_field), allocatable :: yf
    type(set_graph)                   :: on
    type(field)                     :: unit
    real(dp), allocatable                  :: col(:)
    real(dp)                               :: a(5,5), at(5,5), e(5)
    real(dp)                               :: cs(5)
    integer                                :: order, i, j
    logical                                :: ok
    type(set_map)     :: sets

    ! Five vertices, four interior edges, and one boundary edge
    ! hanging off the last vertex.
    g = directed_stored_graph(5, tails=[1,2,3,4,5], heads=[2,3,4,5,0])
    call describe(sets, g)
    on = g % vertex_set()
    unit = field('e', on, sets % size_of(on))

    ! Per-edge coefficients, so the transpose is tested with weights
    ! aboard, not just with ones.
    cs = [2.0_dp, 0.5_dp, 3.0_dp, 1.5_dp, 4.0_dp]

    ok = .true.

    do order = 1, 4

       fwd = vertex_differential_operator(order=order, coefficients=cs)
       rev = vertex_differential_operator(order=order, coefficients=cs, &
            &                             adjoint=.true.)

       do j = 1, 5
          e = 0.0_dp
          e(j) = 1.0_dp
          call unit % set_real_vector(e)

          call fwd % apply(g, [unit], yf)
          call yf % get_real_vector(col)
          a(:, j) = col

          call rev % apply(g, [unit], yf)
          call yf % get_real_vector(col)
          at(:, j) = col
       end do

       do i = 1, 5
          do j = 1, 5
             ok = ok .and. abs(at(i, j) - a(j, i)) < 1.0d-11
          end do
       end do

    end do

    call report(ok, &
         & "the flag is the matrix transpose, orders one to four", nfail)

  end subroutine check_adjoint_is_the_transpose

  !===================================================================!
  ! State-dependent formulas on edges. Two claims, both about the
  ! contract rather than about any formula:
  !
  ! First, the edge leaf accepts any per-edge function of the two end
  ! values - nothing assumes linearity - and its output reduces
  ! through incidence like any other edge field.
  !
  ! Second, conservation belongs to incidence, not to the formula.
  ! Whatever an edge holds, incidence gives it to one vertex and takes
  ! it from the other, so around a closed ring the total is zero for
  ! ANY formula. This is the property that lets a wave-speed formula
  ! be as nonlinear as it likes without breaking conservation.
  !
  ! Third, the lagged-coefficient path: coefficients computed FROM
  ! the state, loaded into a linear operator, applied - the loop a
  ! flow model runs each evaluation.
  !===================================================================!

  subroutine check_nonlinear_edge_formulas(nfail)

    integer, intent(inout) :: nfail

    type(directed_stored_graph)                     :: ring
    type(nonlinear_sample)                 :: formula
    type(differential_operator)     :: reduce_edges
    class(graph_field), allocatable   :: zf
    class(graph_field), allocatable :: yf
    type(set_graph)                   :: on
    type(field)                     :: qf
    type(set_graph)                     :: eon
    type(field)                       :: samples
    real(dp), allocatable                  :: z(:), y(:)
    real(dp)                               :: q(4), c(4)
    integer                                :: e, t, h
    type(set_map)     :: sets

    ring = directed_stored_graph(4, tails=[1,2,3,4], heads=[2,3,4,1])
    call describe(sets, ring)
    on = ring % vertex_set()
    qf   = field('q', on, sets % size_of(on))

    q = [3.0_dp, -1.0_dp, 4.0_dp, 1.5_dp]
    call qf % set_real_vector(q)

    ! The nonlinear formula runs through the ordinary leaf.
    call formula % apply(ring, [qf], zf)
    call zf % get_real_vector(z)
    call report(abs(z(1) - 0.5_dp * (q(1) + q(2)) * abs(q(2) - q(1))) < 1.0d-12, &
         & "an edge formula nonlinear in both ends runs on the leaf", nfail)

    ! Conservation belongs to incidence: the sum of ANY edge field
    ! over a closed ring sums to zero, this formula included. (The
    ! samples ride in a concrete field; gfortran 11 cannot build an
    ! array constructor from a polymorphic item.)
    eon = ring % edge_set()
    samples = field('z', eon, sets % size_of(eon))
    call samples % set_real_vector(z)
    reduce_edges = divergence()
    call reduce_edges % apply(ring, [samples], yf)
    call yf % get_real_vector(y)
    call report(abs(sum(y)) < 1.0d-12, &
         & "incidence conserves any nonlinear formula on a ring", nfail)

    ! The lagged-coefficient loop: compute per-edge numbers from the
    ! state - here the larger end magnitude, the shape of a wave
    ! speed - load them into the linear operator, apply.
    do e = 1, 4
       t = ring % edge_tail(e)
       h = ring % edge_head(e)
       c(e) = max(abs(q(t)), abs(q(h)))
    end do

    reduce_edges = vertex_differential_operator(order=2, coefficients=c)
    call reduce_edges % apply(ring, [qf], yf)
    call yf % get_real_vector(y)

    ! By hand at vertex 2: c2 (q3 - q2) - c1 (q2 - q1), with
    ! c1 = max(3,1) = 3 and c2 = max(1,4) = 4: 4*5 - 3*(-4) = 32.
    call report(abs(y(2) - 32.0_dp) < 1.0d-12, &
         & "coefficients computed from the state drive the operator", nfail)

    ! And the sum still vanishes: weighted incidence conserves too.
    call report(abs(sum(y)) < 1.0d-12, &
         & "and the state-dependent weights still conserve", nfail)

  end subroutine check_nonlinear_edge_formulas

  !===================================================================!
  ! THE INNER PRODUCT. It is not a new machine; it is the reduction
  ! with its measure seat occupied. Reducing u with measure v sums
  ! u times v entry by entry, which is the product <u, v>. These
  ! checks pin that reading, and use it to state integration by
  ! parts and the Laplacian's energy on a chain.
  !===================================================================!

  subroutine check_inner_products(nfail)

    integer, intent(inout) :: nfail

    type(directed_stored_graph)                     :: chain
    type(set_graph)                   :: on
    type(set_graph)                     :: eon
    type(field)                     :: uf, vf, wf
    type(field)                       :: zf, guf
    type(differential_operator)     :: opposite
    type(differential_operator)       :: slope
    type(reduction)                        :: total
    class(graph_functional), allocatable   :: answer
    class(graph_field), allocatable :: work_v
    class(graph_field), allocatable   :: work_e
    real(dp), allocatable                  :: values(:)
    real(dp)                               :: x, left, right
    type(set_map)     :: sets
    complex(dp)                            :: cx

    chain = directed_stored_graph(4, tails=[1, 2, 3], heads=[2, 3, 4])
    call describe(sets, chain)
    on = chain % vertex_set()
    eon = chain % edge_set()

    uf = field('u', on, sets % size_of(on))
    vf = field('v', on, sets % size_of(on))
    call uf % set_real_vector([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp])
    call vf % set_real_vector([2.0_dp, 1.0_dp, 0.0_dp, 3.0_dp])

    ! By hand: 1*2 + 2*1 + 3*0 + 4*3 = 16.
    total = reduction(REDUCE_SUM)
    call total % reduce(uf, answer, measure=vf)
    call value_real(answer, x)
    call report(abs(x - 16.0_dp) < 1.0d-12, &
         & "a sum with a measure is the inner product", nfail)

    ! Integration by parts. G reads the slope of u onto the edges,
    ! D gathers z back onto the vertices, and the two products
    ! cancel: the transpose of G is minus D.
    zf = field('z', eon, sets % size_of(eon))
    call zf % set_real_vector([2.0_dp, 0.0_dp, 1.0_dp])

    slope = gradient()
    call slope % apply(chain, [uf], work_e)
    call work_e % get_real_vector(values)
    guf = field('Gu', eon, sets % size_of(eon))
    call guf % set_real_vector(values)

    call total % reduce(zf, answer, measure=guf)
    call value_real(answer, left)

    opposite = divergence()
    call opposite % apply(chain, [zf], work_v)
    call work_v % get_real_vector(values)
    wf = field('Dz', on, sets % size_of(on))
    call wf % set_real_vector(values)

    call total % reduce(wf, answer, measure=uf)
    call value_real(answer, right)

    call report(abs(left + right) < 1.0d-12, &
         & "<z,Gu> + <Dz,u> = 0: integration by parts on the chain", nfail)

    ! The energy: u = [0,1,3,6] slopes to [1,2,3], and the product
    ! of the slope with itself is 1 + 4 + 9 = 14, the value of the
    ! quadratic form u^T L u.
    call uf % set_real_vector([0.0_dp, 1.0_dp, 3.0_dp, 6.0_dp])
    call slope % apply(chain, [uf], work_e)
    call work_e % get_real_vector(values)
    guf = field('Gu', eon, sets % size_of(eon))
    call guf % set_real_vector(values)
    call total % reduce(guf, answer, measure=guf)
    call value_real(answer, x)
    call report(abs(x - 14.0_dp) < 1.0d-12, &
         & "the product of the gradient with itself is u^T L u", nfail)

    ! The complex road: u + ih against a real measure v; the
    ! imaginary part of the answer is h times the sum of v.
    call uf % set_complex_vector([(1.0_dp, 0.001_dp), (2.0_dp, 0.001_dp), &
         &                        (3.0_dp, 0.001_dp), (4.0_dp, 0.001_dp)])
    call total % reduce(uf, answer, measure=vf)
    call value_complex(answer, cx)
    call report(abs(real(cx, dp) - 16.0_dp) < 1.0d-12 .and. &
         &      abs(aimag(cx) - 0.006_dp) < 1.0d-12, &
         & "a complex field carries its derivative through the product", nfail)

  end subroutine check_inner_products

  !===================================================================!
  ! THE BROADCAST. The reduction's transpose: one value fills a
  ! field. The round trips pin the two rule pairs, the pairing pins
  ! the transpose property, and a complex seed rides intact. No
  ! graph appears anywhere: the field's own support carries
  ! everything a broadcast needs.
  !===================================================================!

  subroutine check_broadcast(nfail)

    integer, intent(inout) :: nfail

    type(set_graph)                 :: on
    type(field)                   :: f, uf
    type(functional)                     :: seed
    type(broadcast)                      :: copy_rule, share_rule
    type(reduction)                      :: total
    class(graph_functional), allocatable :: answer
    real(dp)                             :: x
    type(set_map)     :: sets
    complex(dp)                          :: cx

    call on % declare()
    call sets % bind(on, counted_set_representation(4))
    f  = field('seed', on, sets % size_of(on))
    uf = field('u', on, sets % size_of(on))

    copy_rule  % rule = BROADCAST_COPY
    share_rule % rule = BROADCAST_SHARE
    call seed % set_real_value(6.0_dp)

    ! The average of four copies of 6 is 6.
    call copy_rule % broadcast(seed, f)
    total = reduction(REDUCE_AVERAGE)
    call total % reduce(f, answer)
    call value_real(answer, x)
    call report(abs(x - 6.0_dp) < 1.0d-12, &
         & "reduce(broadcast(J)) = J: the average undoes the copy", nfail)

    ! The sum of four shares of 6 is 6.
    call share_rule % broadcast(seed, f)
    total = reduction(REDUCE_SUM)
    call total % reduce(f, answer)
    call value_real(answer, x)
    call report(abs(x - 6.0_dp) < 1.0d-12, &
         & "reduce(broadcast(J)) = J: the sum undoes the share", nfail)

    ! The pairing <broadcast(J), u> = J * sum(u): 6 * 10 = 60,
    ! computed with the inner-product reading of the measure.
    call uf % set_real_vector([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp])
    call copy_rule % broadcast(seed, f)
    call total % reduce(uf, answer, measure=f)
    call value_real(answer, x)
    call report(abs(x - 60.0_dp) < 1.0d-12, &
         & "<broadcast(J), u> = J sum(u): the transpose pairing", nfail)

    ! A complex seed rides intact: four copies of 6 + 0.001i sum to
    ! 24 + 0.004i.
    call seed % set_complex_value((6.0_dp, 0.001_dp))
    call copy_rule % broadcast(seed, f)
    call total % reduce(f, answer)
    call value_complex(answer, cx)
    call report(abs(real(cx, dp) - 24.0_dp) < 1.0d-12 .and. &
         &      abs(aimag(cx) - 0.004_dp) < 1.0d-12, &
         & "a complex seed is copied and returns intact", nfail)

  end subroutine check_broadcast

  !===================================================================!
  ! Suite helpers. The member list of any graph is its global index
  ! column; one value comes out of any functional through the
  ! contract's vector adapters.
  !===================================================================!

  ! The map is a DUMMY, not a local: the identity handed in is
  ! described there and nowhere else.
  subroutine members_of(sets, g, indices)
    type(set_map)       , intent(in)  :: sets
    type(set_graph)     , intent(in)  :: g
    integer, allocatable, intent(out) :: indices(:)
    call sets % members_of(g, indices)
  end subroutine members_of

  subroutine value_integer(f, x)
    class(graph_functional), intent(in) :: f
    integer, intent(out) :: x
    integer, allocatable :: t(:)
    call f % get_integer_vector(t)
    x = 0
    if (size(t) >= 1) x = t(1)
  end subroutine value_integer

  subroutine value_real(f, x)
    class(graph_functional), intent(in) :: f
    real(dp), intent(out) :: x
    real(dp), allocatable :: t(:)
    call f % get_real_vector(t)
    x = 0.0_dp
    if (size(t) >= 1) x = t(1)
  end subroutine value_real

  subroutine value_complex(f, x)
    class(graph_functional), intent(in) :: f
    complex(dp), intent(out) :: x
    complex(dp), allocatable :: t(:)
    call f % get_complex_vector(t)
    x = (0.0_dp, 0.0_dp)
    if (size(t) >= 1) x = t(1)
  end subroutine value_complex

  subroutine value_logical(f, x)
    class(graph_functional), intent(in) :: f
    logical, intent(out) :: x
    logical, allocatable :: t(:)
    call f % get_logical_vector(t)
    x = .false.
    if (size(t) >= 1) x = t(1)
  end subroutine value_logical

end program test_graph_contract
