!=====================================================================!
! DRIVING A RULE OVER A DIGRAPH.
!
! A marcher traverses instants; a sweep traverses nodes; a smoother
! traverses colours; an adjoint traverses the same vertices in
! reverse. All four are ONE operation: visit the vertices of a
! digraph in an order its arcs admit, and apply a rule at each.
! Nothing in that operation is temporal. This module states it once.
!
!             THE TAXONOMY
!
!   VERTEX          What the rule is applied at. What a vertex
!                   denotes - an instant, a cell, a block - is
!                   decided by the caller and never by this module.
!
!   ARC             An ordered pair of vertices. The arc j -> k states
!                   that the rule at k READS the value at j.
!
!   SLOT            Which argument of k's rule the arc fills. Two
!                   arcs into the same vertex fill different slots:
!                   that is what distinguishes q_(n-1) from q_(n-2).
!
!   LABELLING       The map ell : A -> S storing the slot of every
!                   arc. Stored WITH the arcs, never beside them, so
!                   restricting or transposing the digraph transports
!                   the slots with the arcs and cannot separate them.
!
!   ORIENTATION     forward or reverse. A traversal in the forward
!                   orientation reaches a vertex only after every
!                   vertex it reads; reverse is the same order read
!                   backwards, which is the order an adjoint requires.
!
!   RULE            The operation stored at a vertex. Its arguments
!                   are the slots the labelling names.
!
!   DRIVER          The operation that applies the rule at every
!                   vertex in an admissible order. It is itself an
!                   operation, so a driver over a digraph of drivers
!                   is a driver, and nesting requires no new type.
!
!             THE TWO BRANCHES, AND WHAT THE DRIVER PAIRS
!
! A computation has two graphs over one skeleton, and they are the
! fractal graph's two branches: G = (B1, B2), symmetric in structure
! and independent in meaning.
!
!      B2   the OPERATION graph        what is computed
!
!             (r1) ---> (r2) ---> (r3) ---> (r4)
!               .        .         .         .
!               .        .         .         .        <- the PAIRING
!               v        v         v         v           one datum
!             [ d1 ]   [ d2 ]    [ d3 ]    [ d4 ]        per rule
!
!      B1   the DATA graph             what is computed ON
!
! The arcs are stored in B2: they state which rule reads which. The
! values are stored in B1: one datum per vertex. The PAIRING is the
! bijection between them, and the driver is the object that stores it.
!
! Neither branch is derivable from the other. B2 without B1 is a
! plan with nothing to compute on; B1 without B2 is an array with
! no operation over it. Stored apart and linked, each may be
! restricted, transposed or refined without altering the other -
! which is what lets one skeleton store a model and its correction,
! a stencil and its coefficients, a structure and what is learned on
! it.
!
!             WHY THE LINKAGE REDUCES MEMORY
!
! A datum must be stored from the step its rule writes it until the
! last rule that reads it has been applied - no longer. The driver
! stores both facts: the visiting order, and every arc that reads a
! vertex. So it can report, for each datum, the step after which
! nothing reads it again.
!
!             visiting order   1     2     3     4
!             d1 written       #-----+-----x
!                                          ^ last read by r3,
!                                            released after step 3
!             d2 written             #-----+-----x
!             d3 written                   #-----x
!
! A caller that calls `released_after` receives the vertices whose
! data may be deallocated at each step, so the peak storage is the
! widest live set rather than the whole trajectory. Nothing here
! deallocates anything: the driver reports the lifetime and the
! caller applies it.
!
! The lifetime is only as exact as the arcs. A pass that reads a
! datum without an arc recording that read is invisible here, and the
! result would deallocate what that pass still requires - so a
! reverse pass over the same vertices belongs in the graph as the
! TRANSPOSE, every forward arc with its ends exchanged, plus an arc
! from each datum to the transposed vertex that reads it again. A
! vertex of the transpose need store no rule for this purpose: its
! position in the order is what makes it a reader, and a step that
! computes nothing still releases whatever was last read at it.
!
!             WHAT ASSEMBLES, WHAT DRIVES
!
! Building the two branches is one task and evaluating them is
! another, and they belong to different layers:
!
!      an ASSEMBLER   states the vertices, the arcs and their slots,
!                     which rule is stored where, and which datum. It
!                     depends on the physics, the scheme and the mesh.
!
!      the DRIVER     receives that and evaluates it. It stores the
!                     order, the connection and the lifetimes, and
!                     nothing about what any of it denotes.
!
! So a caller assembles and then passes the pairing:
!
!      link = operations % pair(values)
!      call schedule % pair_with(link)
!      call schedule % evaluate(on)
!
! The separation is deliberate. An assembler may be rewritten - a
! different scheme, a different mesh, a coarser chain - without the
! driver changing, and the driver may be extended to evaluate two
! vertices at once without any assembler changing.
!
!             WHERE PARALLELISM IS DEFINED, AND ONLY THERE
!
! Two vertices with no arc between them may be driven at once. Only
! the driver can detect that, because only here are the order, the
! rules and the data all stored in one object. Every hardware
! decision - which device a rule runs on, where a datum resides, when
! it moves - is a decision about `evaluate` and about nothing above
! or below it.
!
!             THE DIGRAPH IS THE ONLY ORDER
!
! No routine here records a vertex's number, adds one to an
! index, or requests the "next" vertex. The order comes from
! `visiting_order`, which is the digraph's own topological order in
! the given orientation. A digraph with a cycle has no such order and
! is rejected, because a rule that reads its own result is not a rule
! this driver can drive.
!
!             WHAT IS NOT HERE
!
! The MEASURE along the coordinate the vertices discretise - the step
! between two instants, the volume of a cell - belongs to a grid, and
! a grid is an operation of its own. The driver reads a measure when
! a rule requires one and never computes it. That separation is why
! this module has no notion of time.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_driver

  use util_precision       , only : dp
  use operation_action     , only : operation, binding, moved_binding
  use view_directed        , only : directed_graph, forward, reverse
  use view_directed_stored , only : stored_directed_graph
  use view_read_write      , only : bipartite_digraph, FIRST_PART, SECOND_PART
  use field_calculus       , only : field
  use field_stored         , only : stored_field

  implicit none

  private

  public :: driver
  public :: rule_vertex, data_vertex, pairing
  public :: rule_graph, data_graph
  public :: driven_by, paired

  !===================================================================!
  ! ONE VERTEX OF B2. The rule computed there. A vertex whose rule is
  ! not allocated computes nothing, which is how a source vertex - one
  ! that is read but never written - is stated.
  !
  !        (r) ---> ...        the rule, and the arcs it is read by
  !===================================================================!

  type :: rule_vertex

     class(operation), allocatable :: rule

   contains

     procedure :: computes                ! whether a rule is stored here

  end type rule_vertex

  !===================================================================!
  ! ONE VERTEX OF B1. The datum stored there, and nothing about how it
  ! was computed. A datum not yet written is unallocated, and that is
  ! the only way to query whether a step has run.
  !
  !        [ d ]               the value, and nothing else
  !===================================================================!

  type :: data_vertex

     class(field), allocatable :: datum

   contains

     procedure :: written                 ! whether a value is stored here

  end type data_vertex

  !===================================================================!
  ! B2 WHOLE. The rules over every vertex, as one value. Connecting
  ! it to a data graph returns a pairing; which side is called is
  ! immaterial, because the branches are symmetric.
  !===================================================================!

  type :: rule_graph

     type(rule_vertex), allocatable :: at(:)

   contains

     procedure :: pair => rules_pair

  end type rule_graph

  !===================================================================!
  ! B1 WHOLE. The data over every vertex, as one value.
  !===================================================================!

  type :: data_graph

     type(data_vertex), allocatable :: at(:)

   contains

     procedure :: pair => data_pair

  end type data_graph

  !===================================================================!
  ! THE PAIRING. The bijection between B2's vertices and B1's, stored
  ! as one value so neither branch can be indexed without the other
  ! being referenced. Both branches store the same count of vertices
  ! by construction: that IS the bijection, and a constructor that is
  ! passed unequal counts stops the program.
  !
  !        B2   (r1)   (r2)   (r3)
  !               |      |      |          the linkage is the
  !               v      v      v          identity on positions
  !        B1   [ d1 ] [ d2 ] [ d3 ]
  !===================================================================!

  type :: pairing

     private

     type(rule_vertex), allocatable :: computed_by(:)
     type(data_vertex)     , allocatable :: stored_at(:)

   contains

     procedure :: paired_order            ! how many vertices are linked
     procedure :: rule_at                 ! the rule computing a vertex
     procedure :: datum_at                ! the value stored at a vertex
     procedure :: assign                   ! write a value at a vertex
     procedure :: release                 ! deallocate the value at a vertex
     procedure :: stored_data                ! the data branch, to retain
     procedure, private :: require_vertex

  end type pairing

  interface pairing
     module procedure paired
  end interface pairing

  !===================================================================!
  ! THE DRIVER. A rule, a digraph to drive it over, and the
  ! orientation to visit in. Applying it visits every vertex in an
  ! admissible order and applies the rule there.
  !===================================================================!

  type, extends(operation) :: driver

     private

     ! the two parts and the arcs crossing them: which datum every
     ! rule reads, and which it writes
     type(bipartite_digraph)       :: over
     class(operation), allocatable :: rule
     integer                       :: orientation = forward

     ! B1 and B2 paired: allocated once a driver is given them, and
     ! unallocated while a driver stores one rule for every vertex alike
     type(pairing), allocatable    :: stored_pairing

     ! The immutable incidence and orientation determine this schedule
     ! once. Release intervals contain each datum with a reader once;
     ! data without readers are outputs and remain available.
     integer, allocatable :: visiting(:), last_reader(:)
     integer, allocatable :: first_release(:), release_data(:)

   contains

     procedure :: name  => driver_name
     procedure :: apply => driver_apply

     procedure :: visits                  ! the order this driver visits in
     procedure :: driven_rule             ! the rule it stores
     procedure :: pair_with               ! give it B1 paired to B2
     procedure :: pairing_of              ! the pairing it stores
     procedure :: evaluate                ! the rules over the data, in order
     procedure :: last_reader_of          ! the step a datum is last read at
     procedure :: released_after          ! the data releasable after a step
     procedure, private :: neighbourhood

  end type driver

  interface driver
     module procedure driven_by
  end interface driver

contains

  !===================================================================!
  ! A RULE DRIVEN OVER A DIGRAPH, in an orientation. The result
  ! is an operation, so it composes with any other.
  !===================================================================!

  function driven_by(rule, over, orientation) result(this)

    class(operation)       , intent(in) :: rule
    type(bipartite_digraph), intent(in) :: over
    integer               , intent(in), optional :: orientation
    type(driver) :: this

    this % over = over
    allocate(this % rule, source=rule)
    this % orientation = forward
    if (present(orientation)) this % orientation = orientation
    if (this % orientation /= forward .and. this % orientation /= reverse) then
       error stop 'operation_driver: an orientation is forward or reverse'
    end if
    call scheduled(this)
    call this % declare_arguments(rule % num_arguments(), rule % contracts())

  end function driven_by

  pure function driver_name(this) result(name)
    class(driver), intent(in) :: this
    character(len=:), allocatable :: name
    if (allocated(this % rule)) then
       name = this % rule % name() // ' driven over a digraph'
    else
       name = 'driver'
    end if
  end function driver_name

  !===================================================================!
  ! THE ORDER THE RULES ARE VISITED IN: the topological order of the
  ! projection onto the first part. The arcs between the parts state
  ! which rule reads what another wrote; the projection turns that
  ! into an order over the rules alone, and no routine here states
  ! the order explicitly.
  !===================================================================!

  function visits(this) result(order)
    class(driver), intent(in) :: this
    integer, allocatable :: order(:)
    order = [integer ::]
    if (allocated(this % visiting)) order = this % visiting
  end function visits

  subroutine scheduled(this)
    class(driver), intent(inout) :: this
    type(stored_directed_graph) :: among_rules
    integer, allocatable :: readers(:), next(:), step_of_rule(:)
    integer :: k, d, i, n

    among_rules = this % over % projection(FIRST_PART)
    this % visiting = among_rules % loop(this % orientation)
    n = size(this % visiting)
    allocate(step_of_rule(n))
    do k = 1, n
       step_of_rule(this % visiting(k)) = k
    end do
    allocate(this % last_reader(this % over % order_of_part(SECOND_PART)), source=0)
    allocate(this % first_release(n + 1), source=0)
    do d = 1, size(this % last_reader)
       call this % neighbourhood(SECOND_PART, d, .false., readers)
       do i = 1, size(readers)
          this % last_reader(d) = max(this % last_reader(d), step_of_rule(readers(i)))
       end do
       k = this % last_reader(d)
       if (k > 0) this % first_release(k + 1) = this % first_release(k + 1) + 1
    end do
    this % first_release(1) = 1
    do k = 1, n
       this % first_release(k + 1) = this % first_release(k + 1) + this % first_release(k)
    end do
    allocate(this % release_data(this % first_release(n + 1) - 1))
    next = this % first_release(1:n)
    do d = 1, size(this % last_reader)
       k = this % last_reader(d)
       if (k == 0) cycle
       this % release_data(next(k)) = d
       next(k) = next(k) + 1
    end do
  end subroutine scheduled

  function driven_rule(this) result(rule)
    class(driver), intent(in) :: this
    class(operation), allocatable :: rule
    allocate(rule, source=this % rule)
  end function driven_rule

  !===================================================================!
  ! WHETHER A VERTEX STORES ANYTHING. A rule that is not allocated
  ! computes nothing; a datum that is not allocated has not been
  ! written. Both are valid states and neither is an error.
  !===================================================================!

  pure logical function computes(this)
    class(rule_vertex), intent(in) :: this
    computes = allocated(this % rule)
  end function computes

  pure logical function written(this)
    class(data_vertex), intent(in) :: this
    written = allocated(this % datum)
  end function written

  !===================================================================!
  ! PAIR THE TWO BRANCHES. A rule for every vertex of the first part
  ! and a datum for every vertex of the second. The counts need not
  ! agree: which datum a rule reads and which it writes is specified
  ! by the bipartite digraph, and given that specification the two
  ! parts may differ in size. An operation may write two data, a datum
  ! may be read by three operations, and a datum nothing writes is a
  ! source of the digraph.
  !
  ! The branches are symmetric, so either may be called, and there is
  ! nothing else either could be paired with:
  !
  !      link = operations % pair(values)
  !      link = values % pair(operations)
  !
  ! and the two return the same pairing.
  !===================================================================!

  function paired(computed_by, stored_at) result(this)

    type(rule_vertex), intent(in) :: computed_by(:)
    type(data_vertex)     , intent(in) :: stored_at(:)
    type(pairing) :: this

    this % computed_by = computed_by
    this % stored_at     = stored_at

  end function paired

  function rules_pair(this, values) result(link)
    class(rule_graph), intent(in) :: this
    type(data_graph)      , intent(in) :: values
    type(pairing) :: link
    link = paired(this % at, values % at)
  end function rules_pair

  function data_pair(this, operations) result(link)
    class(data_graph)     , intent(in) :: this
    type(rule_graph) , intent(in) :: operations
    type(pairing) :: link
    link = paired(operations % at, this % at)
  end function data_pair

  pure integer function paired_order(this, part)
    class(pairing), intent(in) :: this
    integer       , intent(in) :: part
    paired_order = 0
    if (part == SECOND_PART) then
       if (allocated(this % stored_at)) paired_order = size(this % stored_at)
    else
       if (allocated(this % computed_by)) paired_order = size(this % computed_by)
    end if
  end function paired_order

  !===================================================================!
  ! THE RULE AT A VERTEX and THE VALUE AT A VERTEX. Two reads of
  ! one position, which is what the pairing is for.
  !===================================================================!

  subroutine rule_at(this, vertex, rule)
    class(pairing)               , intent(in)  :: this
    integer                      , intent(in)  :: vertex
    class(operation), allocatable, intent(out) :: rule
    call this % require_vertex(FIRST_PART, vertex)
    if (this % computed_by(vertex) % computes()) &
         & allocate(rule, source=this % computed_by(vertex) % rule)
  end subroutine rule_at

  subroutine datum_at(this, vertex, datum)
    class(pairing)            , intent(in)  :: this
    integer                   , intent(in)  :: vertex
    class(field), allocatable , intent(out) :: datum
    call this % require_vertex(SECOND_PART, vertex)
    if (this % stored_at(vertex) % written()) &
         & allocate(datum, source=this % stored_at(vertex) % datum)
  end subroutine datum_at

  subroutine assign(this, vertex, datum)
    class(pairing), intent(inout) :: this
    integer       , intent(in)    :: vertex
    class(field)  , intent(in)    :: datum
    call this % require_vertex(SECOND_PART, vertex)
    if (allocated(this % stored_at(vertex) % datum)) deallocate(this % stored_at(vertex) % datum)
    allocate(this % stored_at(vertex) % datum, source=datum)
  end subroutine assign

  !===================================================================!
  ! DEALLOCATE THE VALUE AT A VERTEX. The driver reports when this is
  ! permitted; this routine does not check, because a caller may
  ! deallocate a datum it will not read for reasons the digraph does
  ! not record.
  !===================================================================!

  subroutine release(this, vertex)
    class(pairing), intent(inout) :: this
    integer       , intent(in)    :: vertex
    call this % require_vertex(SECOND_PART, vertex)
    if (allocated(this % stored_at(vertex) % datum)) deallocate(this % stored_at(vertex) % datum)
  end subroutine release

  !===================================================================!
  ! THE DATA BRANCH AS ONE VALUE, for a caller that requires the
  ! values after the rules have been applied. The pairing retains its
  ! own copy; this is the final data branch, returned whole.
  !===================================================================!

  function stored_data(this) result(values)
    class(pairing), intent(in) :: this
    type(data_graph) :: values
    if (allocated(this % stored_at)) values % at = this % stored_at
  end function stored_data

  subroutine require_vertex(this, part, vertex)
    class(pairing), intent(in) :: this
    integer       , intent(in) :: part, vertex
    if (vertex < 1 .or. vertex > this % paired_order(part)) then
       error stop 'operation_driver: a vertex is one the pairing stores'
    end if
  end subroutine require_vertex

  !===================================================================!
  ! GIVE A DRIVER ITS CONNECTED BRANCHES. Until connected, a driver
  ! stores one rule for every vertex alike; connected, each vertex
  ! stores its own rule and its own datum.
  !===================================================================!

  subroutine pair_with(this, connection)
    class(driver), intent(inout) :: this
    type(pairing), intent(in)    :: connection
    if (connection % paired_order(FIRST_PART) /= this % over % order_of_part(FIRST_PART) .or. &
        & connection % paired_order(SECOND_PART) /= this % over % order_of_part(SECOND_PART)) then
       error stop 'operation_driver: the pairing labels one vertex of each part'
    end if
    this % stored_pairing = connection
  end subroutine pair_with

  function pairing_of(this) result(stored)
    class(driver), intent(in) :: this
    type(pairing) :: stored
    if (.not. allocated(this % stored_pairing)) then
       error stop 'operation_driver: this driver is not paired with its data'
    end if
    stored = this % stored_pairing
  end function pairing_of

  !===================================================================!
  ! THE NEIGHBOURHOOD READ IN THIS DRIVER'S ORIENTATION. In the forward
  ! orientation a rule reads what entered it and writes what leaves;
  ! in the reverse the two exchange, because reversing an orientation
  ! exchanges previous with next. One traversal, read either way.
  !===================================================================!

  subroutine neighbourhood(this, part, vertex, entering, vertices)
    class(driver), intent(in)  :: this
    integer      , intent(in)  :: part, vertex
    logical      , intent(in)  :: entering   ! the in-neighbourhood, when forward
    integer, allocatable, intent(out) :: vertices(:)
    if ((this % orientation == forward) .eqv. entering) then
       call this % over % in_neighbourhood(part, vertex, vertices)
    else
       call this % over % out_neighbourhood(part, vertex, vertices)
    end if
  end subroutine neighbourhood

  !===================================================================!
  ! EVALUATE. Visit the vertices in an admissible order; at each, read
  ! the data its in-neighbours store, apply the rule stored there,
  ! and place the result. A datum whose last reader has been visited
  ! is reported as releasable by released_after, and is then released.
  !
  !      for each v in visiting order
  !          inputs  <- datum at in_neighbour(v, s), for each slot s
  !          value   <- rule at v, applied to those inputs
  !          place value at v
  !          release every vertex in released_after(this step)
  !
  ! THE ONE PLACE PARALLELISM BELONGS. Two vertices with no arc
  ! between them may be driven at once, and only this routine can
  ! detect that, because only here are the order, the rules and the
  ! data all available. Every hardware decision - where a datum is
  ! stored, which device a rule runs on - is a decision about this
  ! loop and about nothing above or below it.
  !===================================================================!

  subroutine evaluate(this, input_graph)

    class(driver)        , intent(inout) :: this
    class(directed_graph), intent(in)    :: input_graph

    class(operation), allocatable :: rule
    class(field)    , allocatable :: value, stored
    type(binding)   , allocatable :: inputs(:)
    integer, allocatable :: reads(:), writes(:)
    integer :: k, v, i, num_inputs

    if (.not. allocated(this % stored_pairing)) then
       error stop 'operation_driver: a driver is paired with its data before it evaluates'
    end if

    if (.not. allocated(this % visiting)) return
    do k = 1, size(this % visiting)

       v = this % visiting(k)
       call this % stored_pairing % rule_at(v, rule)

       ! A VERTEX STORING NO RULE STILL OCCUPIES A STEP. It computes
       ! nothing, but the arcs entering it are reads like any other,
       ! so a datum's last reader may be located there and the
       ! lifetimes below must be resolved at this step as well.
       if (allocated(rule)) then

          ! what the rule reads
          call this % neighbourhood(FIRST_PART, v, .true., reads)
          ! THE RULE'S ARGUMENT k IS THE DATUM AT THE k-TH VERTEX IT
          ! READS. The binding passes that identity into the rule, so
          ! a rule may receive a datum of its own type rather than
          ! a bare vector of values. A vertex nothing has written yet
          ! leaves its argument unbound.
          if (size(reads) > rule % num_arguments()) then
             error stop 'operation_driver: a rule declares an argument for every vertex it reads'
          end if
          allocate(inputs(size(reads)))
          num_inputs = 0
          do i = 1, size(reads)
             call this % stored_pairing % datum_at(reads(i), stored)
             if (.not. allocated(stored)) cycle
             num_inputs = num_inputs + 1
             inputs(num_inputs) = moved_binding(rule % argument(i), stored)
          end do

          if (num_inputs > 0) then
             call rule % apply(input_graph, inputs(1:num_inputs), value)
          else
             call rule % apply(input_graph, output=value)
          end if

          ! what the rule writes: the data on the other side of it
          if (allocated(value)) then
             call this % neighbourhood(FIRST_PART, v, .false., writes)
             do i = 1, size(writes)
                call this % stored_pairing % assign(writes(i), value)
             end do
          end if

          deallocate(inputs)
          deallocate(rule)

       end if

       ! the data nothing reads again
       do i = this % first_release(k), this % first_release(k + 1) - 1
          call this % stored_pairing % release(this % release_data(i))
       end do

    end do

  end subroutine evaluate

  !===================================================================!
  ! THE STEP A DATUM IS LAST READ AT. Zero records that no rule reads
  ! it. Such data remain available as outputs after evaluation.
  !
  !      order    1     2     3     4
  !      from v         .-----+-----'      last_reader_of(v) = 3
  !===================================================================!

  integer function last_reader_of(this, datum)
    class(driver), intent(in) :: this
    integer      , intent(in) :: datum
    if (datum < 1 .or. datum > this % over % order_of_part(SECOND_PART)) then
       error stop 'operation_driver: a datum belongs to the data part'
    end if
    last_reader_of = this % last_reader(datum)
  end function last_reader_of

  !===================================================================!
  ! THE DATA RELEASABLE AFTER ONE STEP: every vertex already visited
  ! whose last reader is this step or earlier. A caller that releases
  ! these at each step stores only the live set, never the trajectory.
  !
  ! This public query is cumulative and ordered by datum number.
  ! Evaluation uses the stored per-step intervals instead, so each
  ! datum is released exactly once during a traversal.
  !===================================================================!

  function released_after(this, step) result(vertices)
    class(driver), intent(in) :: this
    integer      , intent(in) :: step
    integer, allocatable :: vertices(:)
    integer :: d
    vertices = [integer ::]
    if (.not. allocated(this % last_reader)) return
    vertices = pack([(d, d = 1, size(this % last_reader))], &
         & this % last_reader > 0 .and. this % last_reader <= step)
  end function released_after

  !===================================================================!
  ! THE OPERATION INTERFACE IS NOT THE TRAVERSAL. A driver reads its
  ! inputs from the data it is paired with, not from bindings, so
  ! apply stops the program: the rules are driven by evaluate.
  !===================================================================!

  subroutine driver_apply(this, input_graph, inputs, output)

    class(driver)            , intent(in)    :: this
    class(directed_graph)    , intent(in)    :: input_graph
    type(binding)            , intent(in), optional :: inputs(:)
    class(field), allocatable, intent(inout) :: output

    associate (u1 => this, u2 => input_graph, u3 => present(inputs)); end associate
    if (allocated(output)) deallocate(output)

    error stop 'operation_driver: a driver is driven by evaluate, not applied'

  end subroutine driver_apply

end module operation_driver
