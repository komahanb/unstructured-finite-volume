!=====================================================================!
! DRIVING A RULE OVER A DIGRAPH.
!
! A marcher walks instants; a sweep walks nodes; a smoother walks
! colours; an adjoint walks the same vertices backwards. All four are
! ONE act: visit the vertices of a digraph in an order its arcs
! admit, and apply a rule at each. Nothing in that act is temporal.
! This module states it once.
!
!             THE TAXONOMY
!
!   VERTEX          What the rule is applied at. What a vertex
!                   denotes - an instant, a cell, a block - is the
!                   caller's business and never this module's.
!
!   ARC             An ordered pair of vertices. The arc j -> k says
!                   the rule at k READS the value at j.
!
!   SLOT            Which argument of k's rule the arc fills. Two
!                   arcs into the same vertex fill different slots:
!                   that is what distinguishes q_(n-1) from q_(n-2).
!
!   LABELLING       The map ell : A -> S carrying the slot of every
!                   arc. Held WITH the arcs, never beside them, so
!                   restricting or transposing the digraph carries
!                   the slots along and cannot leave them behind.
!
!   ORIENTATION     forward or reverse. A traversal in the forward
!                   orientation reaches a vertex only after every
!                   vertex it reads; reverse is the same order read
!                   backwards, which is what an adjoint wants.
!
!   RULE            The operation carried at a vertex. Its arguments
!                   are the slots the labelling names.
!
!   DRIVER          The operation that applies the rule at every
!                   vertex in an admissible order. It is itself an
!                   operation, so a driver over a digraph of drivers
!                   is a driver, and nothing new is needed to nest.
!
!             THE TWO BRANCHES, AND WHAT THE DRIVER STITCHES
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
! The arcs live in B2: they say which rule reads which. The values
! live in B1: one datum per vertex. The PAIRING is the bijection
! between them, and the driver is the thing that holds it.
!
! Neither branch is derivable from the other. B2 without B1 is a
! plan with nothing to compute on; B1 without B2 is an array with
! nothing to say about it. Held apart and linked, each may be
! restricted, transposed or refined without disturbing the other -
! which is what lets one skeleton carry a model and its correction,
! a stencil and its coefficients, a structure and what is learned on
! it.
!
!             WHY THE LINKAGE BUYS MEMORY
!
! A datum must live from the moment its rule writes it until the
! last rule that reads it has fired - no longer. The driver knows
! both facts: the visiting order, and every arc that reads a vertex.
! So it can say, for each datum, the step after which nothing will
! read it again.
!
!             visiting order   1     2     3     4
!             d1 written       #-----+-----x
!                                          ^ last read by r3,
!                                            released after step 3
!             d2 written             #-----+-----x
!             d3 written                   #-----x
!
! A caller that asks `released_after` gets the vertices whose data
! may be dropped at each step, so the high-water mark is the widest
! live set rather than the whole trajectory. Nothing here frees
! anything: the driver states the lifetime and the caller obeys it.
!
! The lifetime is only as honest as the arcs. A pass that reads a
! datum without an arc saying so is invisible here, and the answer
! would drop what that pass still wants - so a reverse pass over the
! same vertices belongs in the graph as the TRANSPOSE, every forward
! arc with its ends exchanged, plus an arc from each datum to the
! transposed vertex that reads it again. A vertex of the transpose
! need carry no rule to do its work here: standing in the order is
! what makes it a reader, and a step that computes nothing still
! retires whatever was last read at it.
!
!             WHO ASSEMBLES, WHO DRIVES
!
! Building the two branches is one job and running them is another,
! and they belong to different layers:
!
!      an ASSEMBLER   states the vertices, the arcs and their slots,
!                     which rule stands where, and which datum. It
!                     knows the physics, the scheme and the mesh.
!
!      the DRIVER     is handed that and runs it. It knows the order,
!                     the connection and the lifetimes, and nothing
!                     about what any of it means.
!
! So a caller assembles and then hands over:
!
!      link = operations % pair(values)
!      call runner % pair_with(link)
!      call runner % evaluate(on)
!
! The seam is deliberate. An assembler may be rewritten - a different
! scheme, a different mesh, a coarser chain - without the driver
! changing, and the driver may learn to run two vertices at once
! without any assembler knowing.
!
!             WHERE PARALLELISM LIVES, AND ONLY THERE
!
! Two vertices with no arc between them may be driven at once. Only
! the driver can see that, because only here are the order, the rules
! and the data all in one hand. Every hardware question - which
! device a rule runs on, where a datum resides, when it moves -
! is a question about `evaluate` and about nothing above or below it.
!
!             THE DIGRAPH IS THE ONLY ORDER
!
! No routine here writes a vertex's number down, adds one to an
! index, or asks which vertex is "next". The order comes from
! `visiting_order`, which is the digraph's own topological order in
! the given orientation. A digraph with a cycle has no such order and
! is refused, because a rule that reads its own result is not a rule
! this driver can drive.
!
!             WHAT IS NOT HERE
!
! The MEASURE along the coordinate the vertices discretise - the step
! between two instants, the volume of a cell - belongs to a grid, and
! a grid is an operation of its own. The driver reads a measure when
! a rule asks for one and never computes it. That separation is why
! this module has no notion of time.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_driver

  use util_precision       , only : dp
  use operation_action     , only : operation, argument, contract
  use view_directed        , only : directed_graph, forward, reverse
  use view_directed_stored , only : stored_directed_graph
  use view_read_write      , only : bipartite_digraph, FIRST_PART, SECOND_PART
  use field_calculus       , only : field
  use field_stored         , only : stored_field

  implicit none

  private

  public :: labelled_digraph, driver
  public :: rule_vertex, data_vertex, pairing
  public :: rule_graph, data_graph
  public :: reads, driven_by, paired

  !===================================================================!
  ! THE DIGRAPH WITH ITS SLOTS. A digraph and the labelling of its
  ! arcs, held as one value. The digraph answers what a digraph
  ! answers; the labelling answers which slot an arc fills.
  !===================================================================!

  type :: labelled_digraph

     private

     type(stored_directed_graph) :: arcs
     integer, allocatable        :: slot_of_arc(:)

   contains

     ! the digraph's own questions
     procedure :: order_of_digraph        ! how many vertices
     procedure :: size_of_digraph         ! how many arcs
     procedure :: tail_of                 ! the vertex an arc leaves
     procedure :: head_of                 ! the vertex an arc enters
     procedure :: in_arcs                 ! the arcs entering a vertex
     procedure :: in_degree               ! how many there are

     ! the labelling's question
     procedure :: slot_of                 ! which slot an arc fills
     procedure :: in_neighbour            ! the vertex filling a slot

     ! the order a traversal may visit in
     procedure :: visiting_order

     ! the digraph itself, for a caller that wants it bare
     procedure :: digraph

  end type labelled_digraph

  interface labelled_digraph
     module procedure reads
  end interface labelled_digraph

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

     procedure :: computes                ! whether a rule stands here

  end type rule_vertex

  !===================================================================!
  ! ONE VERTEX OF B1. The datum held there, and nothing about how it
  ! was computed. A datum not yet written is unallocated, and that is
  ! the only way to ask whether a step has run.
  !
  !        [ d ]               the value, and nothing else
  !===================================================================!

  type :: data_vertex

     class(field), allocatable :: datum

   contains

     procedure :: written                 ! whether a value stands here

  end type data_vertex

  !===================================================================!
  ! B2 WHOLE. The rules over every vertex, as one value. Connecting
  ! it to a data graph answers a linkage; which side is asked is
  ! immaterial, because the branches are symmetric.
  !===================================================================!

  type :: rule_graph

     type(rule_vertex), allocatable :: at(:)

   contains

     procedure :: pair => rules_pair
     procedure :: rule_graph_order

  end type rule_graph

  !===================================================================!
  ! B1 WHOLE. The data over every vertex, as one value.
  !===================================================================!

  type :: data_graph

     type(data_vertex), allocatable :: at(:)

   contains

     procedure :: pair => data_pair
     procedure :: data_graph_order

  end type data_graph

  !===================================================================!
  ! THE PAIRING. The bijection between B2's vertices and B1's, held
  ! as one value so neither branch can be indexed without the other
  ! being meant. Both branches carry the same count of vertices by
  ! construction: that IS the bijection, and a constructor that is
  ! handed unequal counts stops the program.
  !
  !        B2   (r1)   (r2)   (r3)
  !               |      |      |          the linkage is the
  !               v      v      v          identity on positions
  !        B1   [ d1 ] [ d2 ] [ d3 ]
  !===================================================================!

  type :: pairing

     private

     type(rule_vertex), allocatable :: computed_by(:)
     type(data_vertex)     , allocatable :: held_at(:)

   contains

     procedure :: paired_order            ! how many vertices are linked
     procedure :: rule_at                 ! the rule computing a vertex
     procedure :: datum_at                ! the value held at a vertex
     procedure :: place                   ! write a value at a vertex
     procedure :: release                 ! drop the value at a vertex
     procedure :: data_held                ! the data branch, to keep
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

     ! B1 and B2 stitched: present once a driver is given them, and
     ! absent while a driver carries one rule for every vertex alike
     type(pairing), allocatable    :: stitched

   contains

     procedure :: name  => driver_name
     procedure :: apply => driver_apply

     procedure :: visits                  ! the order this driver visits in
     procedure :: driven_rule             ! the rule it carries
     procedure :: pair_with               ! give it B1 paired to B2
     procedure :: pairing_of              ! the pairing it holds
     procedure :: evaluate                ! the rules over the data, in order
     procedure :: last_reader_of          ! the step a datum is last read at
     procedure :: released_after          ! the data droppable after a step

  end type driver

  interface driver
     module procedure driven_by
  end interface driver

contains

  !===================================================================!
  ! THE READS RELATION as a labelled digraph. Every arc is given as
  ! the triple it is: the vertex read, the slot it fills, the vertex
  ! reading. Lengths that disagree, or a slot below one, stop the
  ! program - a labelling that names no slot is not a labelling.
  !===================================================================!

  function reads(num_vertices, read_vertex, in_slot, by_vertex) result(this)

    integer, intent(in) :: num_vertices
    integer, intent(in) :: read_vertex(:)   ! the tail of each arc
    integer, intent(in) :: in_slot(:)       ! the slot that arc fills
    integer, intent(in) :: by_vertex(:)     ! the head of each arc
    type(labelled_digraph) :: this

    if (size(read_vertex) /= size(by_vertex) .or. size(in_slot) /= size(by_vertex)) then
       error stop 'operation_driver: one slot and one reader per vertex read'
    end if
    if (num_vertices < 1) then
       error stop 'operation_driver: a digraph carries a vertex at least'
    end if
    if (any(in_slot < 1)) then
       error stop 'operation_driver: a slot is one of the rule''s arguments'
    end if

    this % arcs = stored_directed_graph(num_vertices, tails=read_vertex, heads=by_vertex)
    this % slot_of_arc = in_slot

  end function reads

  !===================================================================!
  ! THE DIGRAPH'S OWN QUESTIONS. Order is the count of vertices and
  ! size the count of arcs, as graph theory names them.
  !===================================================================!

  pure integer function order_of_digraph(this)
    class(labelled_digraph), intent(in) :: this
    order_of_digraph = this % arcs % num_vertices()
  end function order_of_digraph

  pure integer function size_of_digraph(this)
    class(labelled_digraph), intent(in) :: this
    size_of_digraph = this % arcs % num_edges()
  end function size_of_digraph

  pure integer function tail_of(this, arc)
    class(labelled_digraph), intent(in) :: this
    integer                , intent(in) :: arc
    tail_of = this % arcs % edge_tail(arc)
  end function tail_of

  pure integer function head_of(this, arc)
    class(labelled_digraph), intent(in) :: this
    integer                , intent(in) :: arc
    head_of = this % arcs % edge_head(arc)
  end function head_of

  pure integer function slot_of(this, arc)
    class(labelled_digraph), intent(in) :: this
    integer                , intent(in) :: arc
    slot_of = this % slot_of_arc(arc)
  end function slot_of

  !===================================================================!
  ! The arcs entering a vertex, and how many. These are the arcs that
  ! fill the slots of the rule carried there.
  !===================================================================!

  subroutine in_arcs(this, vertex, arcs)
    class(labelled_digraph), intent(in)  :: this
    integer                , intent(in)  :: vertex
    integer, allocatable   , intent(out) :: arcs(:)
    integer, allocatable :: incident(:)
    integer :: e, kept
    call this % arcs % incident_edges(vertex, incident)
    allocate(arcs(size(incident)))
    kept = 0
    do e = 1, size(incident)
       if (this % arcs % edge_head(incident(e)) == vertex) then
          kept = kept + 1
          arcs(kept) = incident(e)
       end if
    end do
    arcs = arcs(1:kept)
  end subroutine in_arcs

  integer function in_degree(this, vertex)
    class(labelled_digraph), intent(in) :: this
    integer                , intent(in) :: vertex
    integer, allocatable :: arcs(:)
    call this % in_arcs(vertex, arcs)
    in_degree = size(arcs)
  end function in_degree

  !===================================================================!
  ! THE IN-NEIGHBOUR ALONG A LABEL. The digraph's in-neighbourhood of
  ! v is N-(v); this is the member of it reached by the arc labelled
  ! s, and zero when that arc is absent. The first instant of a march
  ! has an empty in-neighbourhood, which is a fact about the digraph
  ! and not a failure.
  !===================================================================!

  integer function in_neighbour(this, vertex, slot)
    class(labelled_digraph), intent(in) :: this
    integer                , intent(in) :: vertex, slot
    integer, allocatable :: arcs(:)
    integer :: e
    in_neighbour = 0
    call this % in_arcs(vertex, arcs)
    do e = 1, size(arcs)
       if (this % slot_of_arc(arcs(e)) == slot) then
          in_neighbour = this % arcs % edge_tail(arcs(e))
          return
       end if
    end do
  end function in_neighbour

  !===================================================================!
  ! AN ORDER THE ARCS ADMIT: the digraph's topological order in the
  ! given orientation. Forward reaches a vertex only after every
  ! vertex it reads. A cyclic digraph admits no such order.
  !===================================================================!

  function visiting_order(this, orientation) result(order)
    class(labelled_digraph), intent(in) :: this
    integer                , intent(in) :: orientation
    integer, allocatable :: order(:)
    order = this % arcs % loop(orientation)
  end function visiting_order

  function digraph(this) result(bare)
    class(labelled_digraph), intent(in) :: this
    type(stored_directed_graph) :: bare
    bare = this % arcs
  end function digraph

  !===================================================================!
  ! A RULE DRIVEN OVER A DIGRAPH, in an orientation. What comes back
  ! is an operation, so it composes with any other.
  !===================================================================!

  function driven_by(rule, over, orientation) result(this)

    class(operation)       , intent(in) :: rule
    type(bipartite_digraph), intent(in) :: over
    integer               , intent(in), optional :: orientation
    type(driver) :: this
    type(contract), allocatable :: contracts(:)
    type(argument) :: a
    integer :: k

    this % over = over
    allocate(this % rule, source=rule)
    this % orientation = forward
    if (present(orientation)) this % orientation = orientation
    if (this % orientation /= forward .and. this % orientation /= reverse) then
       error stop 'operation_driver: an orientation is forward or reverse'
    end if
    allocate(contracts(rule % num_arguments()))
    do k = 1, rule % num_arguments()
       a = rule % argument(k)
       contracts(k) = a % contract()
    end do
    call this % declare_arguments(rule % num_arguments(), contracts)

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
  ! projection onto the first part. The arcs between the parts say
  ! which rule reads what another wrote; the projection turns that
  ! into an order over the rules alone, and nothing here states it.
  !===================================================================!

  function visits(this) result(order)
    class(driver), intent(in) :: this
    integer, allocatable :: order(:)
    type(stored_directed_graph) :: among_rules
    among_rules = this % over % projection(FIRST_PART)
    order = among_rules % loop(this % orientation)
  end function visits

  function driven_rule(this) result(rule)
    class(driver), intent(in) :: this
    class(operation), allocatable :: rule
    allocate(rule, source=this % rule)
  end function driven_rule

  !===================================================================!
  ! WHETHER A VERTEX CARRIES ANYTHING. A rule that is not allocated
  ! computes nothing; a datum that is not allocated has not been
  ! written. Both are legitimate states and neither is an error.
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
  ! agree: which datum a rule reads and which it writes is the
  ! bipartite digraph's to say, and once it says it the two parts are
  ! free to differ in size. An operation may write two data, a datum
  ! may be read by three operations, and a datum nothing writes is a
  ! source of the digraph.
  !
  ! The branches are symmetric, so either may be asked, and there is
  ! nothing else either could be paired with:
  !
  !      link = operations % pair(values)
  !      link = values % pair(operations)
  !
  ! and the two answer the same pairing.
  !===================================================================!

  function paired(computed_by, held_at) result(this)

    type(rule_vertex), intent(in) :: computed_by(:)
    type(data_vertex)     , intent(in) :: held_at(:)
    type(pairing) :: this

    this % computed_by = computed_by
    this % held_at     = held_at

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

  pure integer function rule_graph_order(this)
    class(rule_graph), intent(in) :: this
    rule_graph_order = 0
    if (allocated(this % at)) rule_graph_order = size(this % at)
  end function rule_graph_order

  pure integer function data_graph_order(this)
    class(data_graph), intent(in) :: this
    data_graph_order = 0
    if (allocated(this % at)) data_graph_order = size(this % at)
  end function data_graph_order

  pure integer function paired_order(this, part)
    class(pairing), intent(in) :: this
    integer       , intent(in) :: part
    paired_order = 0
    if (part == SECOND_PART) then
       if (allocated(this % held_at)) paired_order = size(this % held_at)
    else
       if (allocated(this % computed_by)) paired_order = size(this % computed_by)
    end if
  end function paired_order

  !===================================================================!
  ! THE RULE AT A VERTEX and THE VALUE AT A VERTEX. Two readings of
  ! one position, which is what the linkage is for.
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
    if (this % held_at(vertex) % written()) &
         & allocate(datum, source=this % held_at(vertex) % datum)
  end subroutine datum_at

  subroutine place(this, vertex, datum)
    class(pairing), intent(inout) :: this
    integer       , intent(in)    :: vertex
    class(field)  , intent(in)    :: datum
    call this % require_vertex(SECOND_PART, vertex)
    if (allocated(this % held_at(vertex) % datum)) deallocate(this % held_at(vertex) % datum)
    allocate(this % held_at(vertex) % datum, source=datum)
  end subroutine place

  !===================================================================!
  ! DROP THE VALUE AT A VERTEX. The driver says when this is allowed;
  ! this routine does not check, because a caller may drop a datum it
  ! knows it will not read for reasons the digraph cannot see.
  !===================================================================!

  subroutine release(this, vertex)
    class(pairing), intent(inout) :: this
    integer       , intent(in)    :: vertex
    call this % require_vertex(SECOND_PART, vertex)
    if (allocated(this % held_at(vertex) % datum)) deallocate(this % held_at(vertex) % datum)
  end subroutine release

  !===================================================================!
  ! THE DATA BRANCH AS ONE VALUE, for a caller that wants the values
  ! after the rules have run. The pairing keeps its own copy; this is
  ! what came to rest, handed over whole.
  !===================================================================!

  function data_held(this) result(values)
    class(pairing), intent(in) :: this
    type(data_graph) :: values
    if (allocated(this % held_at)) values % at = this % held_at
  end function data_held

  subroutine require_vertex(this, part, vertex)
    class(pairing), intent(in) :: this
    integer       , intent(in) :: part, vertex
    if (vertex < 1 .or. vertex > this % paired_order(part)) then
       error stop 'operation_driver: a vertex is one the pairing carries'
    end if
  end subroutine require_vertex

  !===================================================================!
  ! GIVE A DRIVER ITS CONNECTED BRANCHES. Until connected, a driver
  ! carries one rule for every vertex alike; connected, each vertex
  ! carries its own rule and its own datum.
  !===================================================================!

  subroutine pair_with(this, connection)
    class(driver), intent(inout) :: this
    type(pairing), intent(in)    :: connection
    if (connection % paired_order(FIRST_PART) /= this % over % order_of_part(FIRST_PART) .or. &
        & connection % paired_order(SECOND_PART) /= this % over % order_of_part(SECOND_PART)) then
       error stop 'operation_driver: the pairing labels one vertex of each part'
    end if
    this % stitched = connection
  end subroutine pair_with

  function pairing_of(this) result(held)
    class(driver), intent(in) :: this
    type(pairing) :: held
    if (.not. allocated(this % stitched)) then
       error stop 'operation_driver: this driver is not paired with its data'
    end if
    held = this % stitched
  end function pairing_of

  !===================================================================!
  ! EVALUATE. Visit the vertices in an admissible order; at each, read
  ! the data its in-neighbours hold, apply the rule standing there,
  ! and place the result. A datum whose last reader has passed is
  ! reported as droppable, and the caller drops it.
  !
  !      for each v in visiting order
  !          inputs  <- datum at in_neighbour(v, s), for each slot s
  !          value   <- rule at v, applied to those inputs
  !          place value at v
  !          release every vertex in released_after(this step)
  !
  ! THE ONE PLACE PARALLELISM BELONGS. Two vertices with no arc
  ! between them may be driven at once, and only this routine knows
  ! that, because only here are the order, the rules and the data all
  ! in hand. Every hardware question - where a datum lives, which
  ! device a rule runs on - is a question about this loop and about
  ! nothing above or below it.
  !===================================================================!

  subroutine evaluate(this, input_graph)

    class(driver)        , intent(inout) :: this
    class(directed_graph), intent(in)    :: input_graph

    class(operation), allocatable :: rule
    class(field)    , allocatable :: value, held
    class(field), allocatable :: inputs(:)
    integer, allocatable :: order(:), reads(:), writes(:), droppable(:)
    integer :: k, v, i, filled

    if (.not. allocated(this % stitched)) then
       error stop 'operation_driver: a driver is paired with its data before it evaluates'
    end if

    order = this % visits()

    do k = 1, size(order)

       v = order(k)
       call this % stitched % rule_at(v, rule)

       ! A VERTEX CARRYING NO RULE STILL OCCUPIES A STEP. It computes
       ! nothing, but the arcs entering it are reads like any other,
       ! so a datum's last reader may stand there and the lifetimes
       ! below must be settled at this step all the same.
       if (allocated(rule)) then

          ! WHAT THE RULE READS. In the forward orientation a rule reads
          ! what entered it and writes what leaves; in the reverse the
          ! two exchange, because reversing an orientation exchanges
          ! previous with next. One traversal, read either way.
          if (this % orientation == forward) then
             call this % over % in_neighbourhood(FIRST_PART, v, reads)
          else
             call this % over % out_neighbourhood(FIRST_PART, v, reads)
          end if
          ! THE DATA A RULE IS HANDED SHARE ONE TYPE, which the first
          ! of them settles - so a rule may be given a datum of its
          ! own making rather than a bare vector of values. A vertex
          ! nothing has written yet is passed over.
          filled = 0
          do i = 1, size(reads)
             call this % stitched % datum_at(reads(i), held)
             if (.not. allocated(held)) cycle
             if (.not. allocated(inputs)) allocate(inputs(size(reads)), mold=held)
             filled = filled + 1
             call held % place_in(inputs(filled))
             deallocate(held)
          end do

          if (allocated(inputs)) then
             call rule % apply(input_graph, inputs(1:filled), value)
          else
             call rule % apply(input_graph, output=value)
          end if

          ! what the rule writes: the data on the other side of it
          if (allocated(value)) then
             if (this % orientation == forward) then
                call this % over % out_neighbourhood(FIRST_PART, v, writes)
             else
                call this % over % in_neighbourhood(FIRST_PART, v, writes)
             end if
             do i = 1, size(writes)
                call this % stitched % place(writes(i), value)
             end do
          end if

          if (allocated(inputs)) deallocate(inputs)
          deallocate(rule)

       end if

       ! what nothing will read again
       droppable = this % released_after(k)
       do i = 1, size(droppable)
          call this % stitched % release(droppable(i))
       end do

    end do

  end subroutine evaluate

  !===================================================================!
  ! THE STEP A DATUM IS LAST READ AT. Walk the visiting order; the
  ! answer is the latest step whose vertex has an arc from this one.
  ! Zero says nothing reads it, so it may be dropped as soon as it is
  ! written.
  !
  !      order    1     2     3     4
  !      from v         .-----+-----'      last_reader_of(v) = 3
  !===================================================================!

  integer function last_reader_of(this, datum)
    class(driver), intent(in) :: this
    integer      , intent(in) :: datum
    integer, allocatable :: order(:), readers(:)
    integer :: k, i
    last_reader_of = 0
    order = this % visits()
    if (this % orientation == forward) then
       call this % over % out_neighbourhood(SECOND_PART, datum, readers)
    else
       call this % over % in_neighbourhood(SECOND_PART, datum, readers)
    end if
    do k = 1, size(order)
       do i = 1, size(readers)
          if (readers(i) == order(k)) last_reader_of = k
       end do
    end do
  end function last_reader_of

  !===================================================================!
  ! THE DATA DROPPABLE AFTER ONE STEP: every vertex already visited
  ! whose last reader is this step or earlier. A caller that releases
  ! these as it goes holds only the live set, never the trajectory.
  !
  ! Every datum here is asked for its last reader, and each such
  ! answer orders the vertices afresh, so a traversal that calls this
  ! at every step costs the order once per step per datum. Against a
  ! graph whose vertices carry a linear solve apiece this does not
  ! show, and it is why a graph carrying its transpose - twice the
  ! vertices - measures the same. A graph of many cheap vertices
  ! wants the order held once instead.
  !===================================================================!

  function released_after(this, step) result(vertices)
    class(driver), intent(in) :: this
    integer      , intent(in) :: step
    integer, allocatable :: vertices(:), keep(:)
    integer :: d, kept, n
    n = this % over % order_of_part(SECOND_PART)
    allocate(keep(n))
    kept = 0
    do d = 1, n
       if (this % last_reader_of(d) > 0 .and. this % last_reader_of(d) <= step) then
          kept = kept + 1
          keep(kept) = d
       end if
    end do
    vertices = keep(1:kept)
  end function released_after

  !===================================================================!
  ! THE TRAVERSAL. Visit the vertices in an admissible order and
  ! apply the rule at each. The rule reads what its slots are filled
  ! by, which the digraph answers and this routine never computes.
  !===================================================================!

  subroutine driver_apply(this, input_graph, input_data, output)

    class(driver)            , intent(in)    :: this
    class(directed_graph)    , intent(in)    :: input_graph
    class(field)             , intent(in), optional :: input_data(:)
    class(field), allocatable, intent(inout) :: output

    integer, allocatable :: order(:)
    integer :: v, k

    if (.not. allocated(this % rule)) then
       error stop 'operation_driver: a driver carries the rule it drives'
    end if

    order = this % visits()

    do k = 1, size(order)
       v = order(k)
       ! the rule at v reads the vertices filling its slots; gathering
       ! those values is the caller's field layout and is not settled
       ! in this sketch
       call this % rule % apply(input_graph, input_data, output)
    end do

  end subroutine driver_apply

end module operation_driver
