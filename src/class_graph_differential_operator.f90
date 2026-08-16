!=====================================================================!
! Differential operators on a graph.
!
! Two operators live here, one per side, and everything both of them
! do is built from three elementary steps. The theory guide
! (doc/graph-differential-operators.pdf) states and proves what these
! comments draw.
!
!
!                     THE THREE ELEMENTARY STEPS
!
! Values live on vertices or on edges, and the steps move them:
!
!                 q_i                       q_j
!                (i) ---------------------> (j)
!                            z_e
!
!    S   the edge average step   vertices -> edges
!
!        z_e = (q_i + q_j) / 2,  or one end, chosen by the
!                                sign of the coefficient
!
!        The half is not arbitrary: the mean of a straight-line q
!        along the edge - its integral over the edge divided by the
!        edge's length - is exactly (q_i + q_j)/2. The one-sided
!        variant reads the same integral at one end.
!
!    G   the difference step     vertices -> edges
!
!        z_e = (q_j - q_i) / h_e
!
!    D   the incidence step      edges -> vertices
!
!             \    |    /
!              v   v   v                    out - in
!            ---> (v) --->        y_v  =  ------------
!              ^   ^   ^                      m_v
!             /    |    \
!
!        the values on edges leaving v, minus the values on
!        edges entering v, over the measure of v.
!        An edge with no head contributes to its tail alone.
!
!
!                    ANY ORDER, BY REPEATING THEM
!
! With q the vertex field and c the coefficient the operator carries:
!
!    vertex side                        edge side
!
!    order 0    c q                     order 0    S q
!    order 1    D S q                   order 1    G q
!    order 2    D G q                   order 2    G D S q
!    order 3    D G D S q               order 3    G D G q
!    order 4    D G D G q               order 4    G D G D S q
!
! The recurrence, stated once:
!
!    vertex(0) = c q          vertex(n) = D of edge(n - 1)
!    edge(0)   = S q          edge(n)   = G of vertex(n - 1)
!
! which closes by parity on the vertex side:
!
!    vertex(2k)     = (D G)^k
!    vertex(2k + 1) = (D G)^k D S
!
! The coefficient is applied once, at the innermost step - the first
! moment the values land on edges - so a per-edge coefficient makes
! order 2 the operator div(k grad q) with k varying edge to edge.
!
! Each step consults an edge's two ends, so order n reaches exactly
! n RINGS of neighbours. Ring r of a vertex is the set a walk of
! exactly r edges reaches and no shorter walk does - the level sets
! of the depth walk, which the engine already computes:
!
!             2   2   2
!           2   1   1   2         ring 1: the neighbours
!           2   1 (v) 1   2       ring 2: the neighbours' neighbours,
!           2   1   1   2                 minus what ring 1 took
!             2   2   2
!
! The operators are defined on any graph - every parameter is
! per-entity. The UNIFORM CHAIN is not an assumption; it is where
! exactness is claimed, because there the discrete formulas coincide
! with the calculus ones: a straight line has zero second derivative,
! x squared has second derivative two, x to the fourth has fourth
! derivative twenty-four. The test suite checks those numbers.
!
!
!                    TWO INCIDENCE DIRECTIONS
!
! The incidence step above counts out minus in, which gives the
! derivatives their textbook signs. The balance counts in minus out,
! because a balance measures what a vertex gains. One sign apart;
! each is stated where it is used.
!
!
!                  THE ADJOINT IS THE REVERSE WALK
!
!    the walk                       the walk, reversed
!
!    (o)-->(o)-->(o)                (o)<--(o)<--(o)
!
! Reversing every edge swaps tail with head, so out becomes in: the
! transposed incidence step is a difference, the transposed
! difference an incidence step, each with one sign flipped, and the
! transposed edge average returns each value - unsigned - to the end
! it was read from. An adjoint runs the transposed steps in reverse
! order.
!
! Parity decides the adjoint's character:
!
!    even orders    vertex(2k) = (D G)^k is its own adjoint whenever
!                   the same weights sit on its steps - the transpose
!                   of a power is the power of the transpose
!
!    odd orders     the adjoint reverses the sampled end: a one-sided
!                   operator's adjoint is one-sided the other way,
!                   the downstream walk of an upstream sample
!
! And no walk must run end to end. Transposes factor stage by stage,
!
!    (C B A) transposed = A* B* C*
!
! so any contiguous segment of a reverse walk is itself a valid
! adjoint - the sensitivity of the answer with respect to that
! segment's input. Checkpointed time marches and mid-chain
! sensitivities are segment walks; they belong to the pipeline,
! which walks stages by construction.
!
!
!                     CURL, WITH NOTHING NEW
!
! Around a loop of edges there is one more derivative: add the edge
! values as the loop is walked. Which edge borders which face is
! itself a graph - one vertex per edge, one per face, a connection
! where an edge borders a face, coefficient +1 or -1 by direction -
! and the incidence step on that graph is the curl. The suite walks a
! square: a difference field sums to zero around it, a circulating
! field does not.
!
!
!                        WHAT IS NOT HERE
!
! No physical names and no physical signs. This layer computes
! derivatives; models name them, each order carrying its own physics
! in its own class:
!
!    order 0    storage, reaction, mass
!    order 1    transport along a flow
!    order 2    diffusion, conduction, viscosity, pressure fields
!    order 3    dispersion - waves whose speed depends on length
!    order 4    bending - beams, plates, and the interfaces of
!               phase separation
!    order 6    pattern-forming films and crystals
!
! The signs those models require - a flux running down its gradient -
! are theirs to state, in their own classes.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_graph_differential_operator

  use iso_fortran_env    , only : dp => REAL64
  use graph_operation_view, only : graph_operation
  use graph_directed_view, only : directed_graph
  use graph_field_calculus, only : graph_field
  use fractal_graph      , only : set_graph => graph
  use graph_calculus     , only : GRAPH_SIDE_VERTEX, GRAPH_SIDE_EDGE
  use class_graph_field  , only : field

  implicit none

  private
  public :: differential_operator
  public :: edge_differential_operator
  public :: vertex_differential_operator
  public :: gradient, interpolation, divergence, laplacian

  !===================================================================!
  ! The shared parameters of both operators. Each is one number for
  ! the uniform case, with an optional per-entity array that wins
  ! when it is allocated:
  !
  !    coefficient   c     applied at the innermost step
  !                        (per edge; per vertex at order 0 on the
  !                        vertex side)
  !    spacing       h_e   the length of an edge, in G
  !    measure       m_v   the size of a vertex, in D
  !    boundary      b_e   the value standing in for the missing end
  !                        of an edge with no head, in S and G
  !===================================================================!

  !===================================================================!
  ! ONE OPERATOR, TWO LANDINGS. The derivative of order n, landing on
  ! the side the constructor chose:
  !
  !    landing = edge      order 0   S q    the average on each edge
  !                        order 1   G q    the slope along each edge
  !                        order n   G of the vertex result of n - 1
  !
  !    landing = vertex    order 0   c q    the term itself
  !                        order 1   D S q  first derivative
  !                        order 2   D G q  second derivative
  !                        order n   keep going
  !
  ! Handed an EDGE field on the vertex landing, the chain enters at
  ! the incidence step: order 1 is then the divergence of that field.
  ! With `adjoint` true the operator runs the reverse walk: the
  ! transposed steps, in reverse order.
  !
  ! INTERPRETED AND COMPILED. This operator and the stencil operator
  ! are the same mathematics in two execution styles. This one
  ! INTERPRETS: it reads the incidence at every apply, matrix-free,
  ! always fresh - the right default. The stencil is the COMPILED
  ! form, weights computed once and walked many times. Neither
  ! learns the other's business.
  !===================================================================!

  type, extends(graph_operation) :: differential_operator

     integer :: landing = GRAPH_SIDE_VERTEX
     integer :: order   = 2

     logical :: adjoint   = .false.
     logical :: one_sided = .false.

     real(dp)              :: coefficient    = 1.0_dp
     real(dp), allocatable :: coefficients(:)
     real(dp)              :: spacing        = 1.0_dp
     real(dp), allocatable :: spacings(:)
     real(dp)              :: measure        = 1.0_dp
     real(dp), allocatable :: measures(:)
     real(dp)              :: boundary_value = 0.0_dp
     real(dp), allocatable :: boundary_values(:)

     character(len=:), allocatable :: label

   contains

     procedure :: name   => operator_name
     procedure :: domain => operator_domain
     procedure :: apply  => operator_apply

  end type differential_operator

contains

  !===================================================================!
  ! Constructors. Order first; every parameter optional; each array
  ! wins over its scalar when given.
  !===================================================================!

  pure type(differential_operator) function edge_differential_operator &
       & (order, coefficient, coefficients, spacing, spacings, &
       &  measure, measures, boundary_value, boundary_values, one_sided, &
       &  label) result(this)

    integer         , intent(in)           :: order
    real(dp)        , intent(in), optional :: coefficient
    real(dp)        , intent(in), optional :: coefficients(:)
    real(dp)        , intent(in), optional :: spacing
    real(dp)        , intent(in), optional :: spacings(:)
    real(dp)        , intent(in), optional :: measure
    real(dp)        , intent(in), optional :: measures(:)
    real(dp)        , intent(in), optional :: boundary_value
    real(dp)        , intent(in), optional :: boundary_values(:)
    logical         , intent(in), optional :: one_sided
    character(len=*), intent(in), optional :: label

    this % landing = GRAPH_SIDE_EDGE
    this % order   = max(order, 0)

    if (present(one_sided)) this % one_sided = one_sided

    if (present(coefficient))     this % coefficient    = coefficient
    if (present(coefficients))    allocate(this % coefficients, source=coefficients)
    if (present(spacing))         this % spacing        = spacing
    if (present(spacings))        allocate(this % spacings, source=spacings)
    if (present(measure))         this % measure        = measure
    if (present(measures))        allocate(this % measures, source=measures)
    if (present(boundary_value))  this % boundary_value = boundary_value
    if (present(boundary_values)) allocate(this % boundary_values, source=boundary_values)
    if (present(label))           this % label          = label

  end function edge_differential_operator

  !===================================================================!
  ! The same constructor, for the vertex side, plus the adjoint flag:
  ! raised, the operator applies its transpose.
  !===================================================================!

  pure type(differential_operator) function vertex_differential_operator &
       & (order, coefficient, coefficients, spacing, spacings, &
       &  measure, measures, boundary_value, boundary_values, adjoint, label) result(this)

    integer         , intent(in)           :: order
    real(dp)        , intent(in), optional :: coefficient
    real(dp)        , intent(in), optional :: coefficients(:)
    real(dp)        , intent(in), optional :: spacing
    real(dp)        , intent(in), optional :: spacings(:)
    real(dp)        , intent(in), optional :: measure
    real(dp)        , intent(in), optional :: measures(:)
    real(dp)        , intent(in), optional :: boundary_value
    real(dp)        , intent(in), optional :: boundary_values(:)
    logical         , intent(in), optional :: adjoint
    character(len=*), intent(in), optional :: label

    this % landing = GRAPH_SIDE_VERTEX
    this % order   = max(order, 0)

    if (present(coefficient))     this % coefficient    = coefficient
    if (present(coefficients))    allocate(this % coefficients, source=coefficients)
    if (present(spacing))         this % spacing        = spacing
    if (present(spacings))        allocate(this % spacings, source=spacings)
    if (present(measure))         this % measure        = measure
    if (present(measures))        allocate(this % measures, source=measures)
    if (present(boundary_value))  this % boundary_value = boundary_value
    if (present(boundary_values)) allocate(this % boundary_values, source=boundary_values)
    if (present(adjoint))         this % adjoint        = adjoint
    if (present(label))           this % label          = label

  end function vertex_differential_operator

  !===================================================================!
  ! The named layer. Four operators every equation reaches for, as
  ! plain functions, so a call site reads as the calculus:
  !
  !    gradient(...)        the slope along each edge
  !    interpolation(...)   the average onto each edge
  !    divergence(...)      the incidence step on an edge field
  !    laplacian(...)       the second derivative at each vertex
  !===================================================================!

  pure type(differential_operator) function gradient(coefficient, coefficients, &
       & spacing, spacings, boundary_value, boundary_values) result(this)

    real(dp), intent(in), optional :: coefficient, coefficients(:)
    real(dp), intent(in), optional :: spacing, spacings(:)
    real(dp), intent(in), optional :: boundary_value, boundary_values(:)

    this = edge_differential_operator(order=1, coefficient=coefficient, &
         & coefficients=coefficients, spacing=spacing, spacings=spacings, &
         & boundary_value=boundary_value, boundary_values=boundary_values, &
         & label='gradient')

  end function gradient

  pure type(differential_operator) function interpolation(coefficient, coefficients, &
       & boundary_value, boundary_values) result(this)

    real(dp), intent(in), optional :: coefficient, coefficients(:)
    real(dp), intent(in), optional :: boundary_value, boundary_values(:)

    this = edge_differential_operator(order=0, coefficient=coefficient, &
         & coefficients=coefficients, boundary_value=boundary_value, &
         & boundary_values=boundary_values, label='interpolation')

  end function interpolation

  pure type(differential_operator) function divergence(coefficient, coefficients, &
       & measure, measures) result(this)

    real(dp), intent(in), optional :: coefficient, coefficients(:)
    real(dp), intent(in), optional :: measure, measures(:)

    this = vertex_differential_operator(order=1, coefficient=coefficient, &
         & coefficients=coefficients, measure=measure, measures=measures, &
         & label='divergence')

  end function divergence

  pure type(differential_operator) function laplacian(coefficient, coefficients, &
       & spacing, spacings, measure, measures, boundary_value, boundary_values) result(this)

    real(dp), intent(in), optional :: coefficient, coefficients(:)
    real(dp), intent(in), optional :: spacing, spacings(:)
    real(dp), intent(in), optional :: measure, measures(:)
    real(dp), intent(in), optional :: boundary_value, boundary_values(:)

    this = vertex_differential_operator(order=2, coefficient=coefficient, &
         & coefficients=coefficients, spacing=spacing, spacings=spacings, &
         & measure=measure, measures=measures, boundary_value=boundary_value, &
         & boundary_values=boundary_values, label='laplacian')

  end function laplacian

  !===================================================================!
  ! Names. A named operator answers with its name; any other answers
  ! with its order.
  !===================================================================!

  pure function operator_name(this) result(name)

    class(differential_operator), intent(in) :: this
    character(len=:), allocatable            :: name

    if (allocated(this % label)) then
       name = this % label
    else
       name = order_name(this % order)
    end if

  end function operator_name

  !===================================================================!
  ! The spelled-out order, for an operator with no name of its own.
  !===================================================================!

  pure function order_name(order) result(name)

    integer, intent(in)           :: order
    character(len=:), allocatable :: name

    character(len=12) :: digits

    write(digits, '(i0)') order
    name = 'derivative of order ' // trim(digits)

  end function order_name

  !===================================================================!
  ! Supports: every edge, every vertex. An operator aimed at a subset
  ! is a second instance handed that subset's graph.
  !===================================================================!

  subroutine operator_domain(this, input_graph, domain, nentries)

    class(differential_operator), intent(in) :: this
    class(directed_graph), intent(in)                 :: input_graph
    type(set_graph), intent(out) :: domain
    integer        , intent(out) :: nentries

    if (this % landing == GRAPH_SIDE_EDGE) then
       domain   = input_graph % all_edges()
       nentries = input_graph % num_edges()
    else
       domain   = input_graph % all_vertices()
       nentries = input_graph % num_vertices()
    end if

  end subroutine operator_domain

  !===================================================================!
  ! Parameter lookups. One line each, so the kernels read identically
  ! in the uniform and the varying case.
  !===================================================================!

  pure real(dp) function coefficient_at(uniform, varying, e)

    real(dp)             , intent(in) :: uniform
    real(dp), allocatable, intent(in) :: varying(:)
    integer              , intent(in) :: e

    if (allocated(varying)) then
       coefficient_at = varying(min(e, size(varying)))
    else
       coefficient_at = uniform
    end if

  end function coefficient_at

  !===================================================================!
  ! THE AVERAGE STEP, forward and reversed.
  !
  !    forward           q_i     q_j              reversed
  !                       \      /
  !                        v    v            z_e goes back, whole or
  !        z_e = (q_i + q_j) / 2             halved, to the end or
  !        or one end, by sign of c          ends it was read from
  !
  ! The coefficient rides along when this is the innermost step.
  ! An edge with no head reads the stored boundary value in place of
  ! the head; on the way back, that column belongs to no vertex and
  ! is dropped.
  !===================================================================!

  pure subroutine average_step(g, one_sided_by, with_c, op_c, op_cs, op_b, op_bs, q, z)

    class(directed_graph), intent(in)          :: g
    real(dp)    , intent(in)          :: one_sided_by   ! sign chooses the end; zero averages
    logical     , intent(in)          :: with_c
    real(dp)    , intent(in)          :: op_c
    real(dp), allocatable, intent(in) :: op_cs(:)
    real(dp)    , intent(in)          :: op_b
    real(dp), allocatable, intent(in) :: op_bs(:)
    real(dp)    , intent(in)          :: q(:)
    real(dp)    , intent(out)         :: z(:)

    real(dp) :: qt, qh, c, pick
    integer  :: e, t, h

    do e = 1, size(z)

       t  = g % edge_tail(e)
       qt = q(t)

       if (g % edge_has_head(e)) then
          h  = g % edge_head(e)
          qh = q(h)
       else
          qh = coefficient_at(op_b, op_bs, e)
       end if

       c = 1.0_dp
       if (with_c) c = coefficient_at(op_c, op_cs, e)

       pick = one_sided_by
       if (with_c) pick = sign(1.0_dp, c) * merge(1.0_dp, 0.0_dp, one_sided_by /= 0.0_dp)

       if (pick > 0.0_dp) then
          z(e) = c * qt                     ! the end the walk leaves
       else if (pick < 0.0_dp) then
          z(e) = c * qh                     ! the end the walk enters
       else
          z(e) = c * 0.5_dp * (qt + qh)     ! both ends, evenly
       end if

    end do

  end subroutine average_step

  !===================================================================!
  ! The transpose of the average step. Reading an edge from its two
  ! ends transposes to returning the edge value to those ends:
  !
  !    z_e ---> y_i and y_j     both ends, half each - or all of it
  !                             to the sampled end, unsigned, when
  !                             the step was one-sided
  !===================================================================!

  pure subroutine average_step_reversed(g, one_sided_by, with_c, op_c, op_cs, z, y)

    class(directed_graph), intent(in)          :: g
    real(dp)    , intent(in)          :: one_sided_by
    logical     , intent(in)          :: with_c
    real(dp)    , intent(in)          :: op_c
    real(dp), allocatable, intent(in) :: op_cs(:)
    real(dp)    , intent(in)          :: z(:)
    real(dp)    , intent(inout)       :: y(:)

    real(dp) :: c, pick
    integer  :: e, t, h

    y = 0.0_dp

    do e = 1, size(z)

       t = g % edge_tail(e)

       c = 1.0_dp
       if (with_c) c = coefficient_at(op_c, op_cs, e)

       pick = one_sided_by
       if (with_c) pick = sign(1.0_dp, c) * merge(1.0_dp, 0.0_dp, one_sided_by /= 0.0_dp)

       if (pick > 0.0_dp) then
          y(t) = y(t) + c * z(e)
       else if (pick < 0.0_dp) then
          if (g % edge_has_head(e)) then
             h = g % edge_head(e)
             y(h) = y(h) + c * z(e)
          end if
       else
          y(t) = y(t) + c * 0.5_dp * z(e)
          if (g % edge_has_head(e)) then
             h = g % edge_head(e)
             y(h) = y(h) + c * 0.5_dp * z(e)
          end if
       end if

    end do

  end subroutine average_step_reversed

  !===================================================================!
  ! THE DIFFERENCE STEP, forward and reversed.
  !
  !    forward                              reversed
  !
  !    (i) ------> (j)                  (i) <------ (j)
  !
  !         q_j - q_i                    -z_e / h to the tail,
  !    z_e = ---------                   +z_e / h to the head:
  !            h_e                       an incidence, one sign flipped
  !===================================================================!

  pure subroutine difference_step(g, with_c, op_c, op_cs, op_h, op_hs, op_b, op_bs, q, z)

    class(directed_graph), intent(in)          :: g
    logical     , intent(in)          :: with_c
    real(dp)    , intent(in)          :: op_c
    real(dp), allocatable, intent(in) :: op_cs(:)
    real(dp)    , intent(in)          :: op_h
    real(dp), allocatable, intent(in) :: op_hs(:)
    real(dp)    , intent(in)          :: op_b
    real(dp), allocatable, intent(in) :: op_bs(:)
    real(dp)    , intent(in)          :: q(:)
    real(dp)    , intent(out)         :: z(:)

    real(dp) :: qt, qh, c
    integer  :: e, t, h

    do e = 1, size(z)

       t  = g % edge_tail(e)
       qt = q(t)

       if (g % edge_has_head(e)) then
          h  = g % edge_head(e)
          qh = q(h)
       else
          qh = coefficient_at(op_b, op_bs, e)
       end if

       c = 1.0_dp
       if (with_c) c = coefficient_at(op_c, op_cs, e)

       z(e) = c * (qh - qt) / coefficient_at(op_h, op_hs, e)

    end do

  end subroutine difference_step

  !===================================================================!
  ! The transpose of the difference step. A difference read off two
  ! ends transposes to a signed return:
  !
  !    z_e / h_e   subtracted at the tail, added at the head
  !
  ! - an incidence, with the measure's seat taken by the spacing.
  !===================================================================!

  pure subroutine difference_step_reversed(g, with_c, op_c, op_cs, op_h, op_hs, z, y)

    class(directed_graph), intent(in)          :: g
    logical     , intent(in)          :: with_c
    real(dp)    , intent(in)          :: op_c
    real(dp), allocatable, intent(in) :: op_cs(:)
    real(dp)    , intent(in)          :: op_h
    real(dp), allocatable, intent(in) :: op_hs(:)
    real(dp)    , intent(in)          :: z(:)
    real(dp)    , intent(inout)       :: y(:)

    real(dp) :: c, w
    integer  :: e, t, h

    y = 0.0_dp

    do e = 1, size(z)

       c = 1.0_dp
       if (with_c) c = coefficient_at(op_c, op_cs, e)

       w = c * z(e) / coefficient_at(op_h, op_hs, e)

       t = g % edge_tail(e)
       y(t) = y(t) - w

       if (g % edge_has_head(e)) then
          h = g % edge_head(e)
          y(h) = y(h) + w
       end if

    end do

  end subroutine difference_step_reversed

  !===================================================================!
  ! THE INCIDENCE STEP, forward and reversed.
  !
  !    forward                              reversed
  !
  !         out - in                        q_tail   q_head
  !    y_v = --------                z_e =  ------ - ------
  !            m_v                          m_tail   m_head
  !
  ! Out minus in gives the derivatives their textbook signs; the
  ! balance counts in minus out because it measures what a vertex
  ! gains. An edge with no head contributes to its tail alone.
  !===================================================================!

  pure subroutine incidence_step(g, op_m, op_ms, z, y)

    class(directed_graph), intent(in)          :: g
    real(dp)    , intent(in)          :: op_m
    real(dp), allocatable, intent(in) :: op_ms(:)
    real(dp)    , intent(in)          :: z(:)
    real(dp)    , intent(out)         :: y(:)

    integer :: e, t, h, v

    y = 0.0_dp

    do e = 1, size(z)

       t = g % edge_tail(e)
       y(t) = y(t) + z(e)                   ! the edge leaves its tail

       if (g % edge_has_head(e)) then
          h = g % edge_head(e)
          y(h) = y(h) - z(e)                ! and enters its head
       end if

    end do

    do v = 1, size(y)
       y(v) = y(v) / coefficient_at(op_m, op_ms, v)
    end do

  end subroutine incidence_step

  !===================================================================!
  ! The transpose of the incidence step. Out minus in at each vertex
  ! transposes to, per edge: the tail's value over its measure minus
  ! the head's value over its measure - a difference, read the other
  ! way round.
  !===================================================================!

  pure subroutine incidence_step_reversed(g, op_m, op_ms, q, z)

    class(directed_graph), intent(in)          :: g
    real(dp)    , intent(in)          :: op_m
    real(dp), allocatable, intent(in) :: op_ms(:)
    real(dp)    , intent(in)          :: q(:)
    real(dp)    , intent(out)         :: z(:)

    integer :: e, t, h

    do e = 1, size(z)

       t = g % edge_tail(e)
       z(e) = q(t) / coefficient_at(op_m, op_ms, t)

       if (g % edge_has_head(e)) then
          h = g % edge_head(e)
          z(e) = z(e) - q(h) / coefficient_at(op_m, op_ms, h)
       end if

    end do

  end subroutine incidence_step_reversed

  !===================================================================!
  ! COMPONENTS. A field may carry several values per entry - the
  ! three parts of a velocity, the five of a conserved state. The
  ! ordering rule interleaves them:
  !
  !      entry           1        1        2        2
  !      component       1        2        1        2
  !                   +--------+--------+--------+--------+--
  !      flat vector  |  v(1)  |  v(2)  |  v(3)  |  v(4)  |
  !                   +--------+--------+--------+--------+--
  !
  ! The kernels stay scalar. The walk gathers one component into a
  ! contiguous work array, runs the chain, and scatters the result
  ! back into its slot - one pass of the graph per component. The
  ! coefficients, spacings, measures and boundary values are shared
  ! by all components; a component that needs its own gets its own
  ! operator instance.
  !===================================================================!

  pure subroutine gather_component(flat, ncomp, c, comp)

    real(dp), intent(in)  :: flat(:)
    integer , intent(in)  :: ncomp, c
    real(dp), intent(out) :: comp(:)

    integer :: i

    do i = 1, size(comp)
       comp(i) = flat((i - 1) * ncomp + c)
    end do

  end subroutine gather_component

  !===================================================================!
  ! The way back: write component c into its interleaved seats,
  ! flat((entry - 1) n + c) - the ordering law, inverted.
  !===================================================================!

  pure subroutine scatter_component(comp, ncomp, c, flat)

    real(dp), intent(in)    :: comp(:)
    integer , intent(in)    :: ncomp, c
    real(dp), intent(inout) :: flat(:)

    integer :: i

    do i = 1, size(comp)
       flat((i - 1) * ncomp + c) = comp(i)
    end do

  end subroutine scatter_component

  !===================================================================!
  ! THE EDGE OPERATOR. Order 0 is the average step; order 1 is the
  ! difference step; order n is the difference step applied to the
  ! vertex result of n - 1.
  !
  !      q on vertices
  !          |
  !          |  (vertex chain of order n - 1, when n > 1)
  !          v
  !      G or S  ------>  the derivative, sampled on edges
  !===================================================================!

  subroutine operator_apply(this, input_graph, input_data, output)

    class(differential_operator), intent(in)       :: this
    class(directed_graph), intent(in)                       :: input_graph
    class(graph_field), intent(in), optional       :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output

    type(field) :: out

    if (this % landing == GRAPH_SIDE_EDGE) then
       call apply_on_edges(this, input_graph, input_data, out)
    else
       call apply_on_vertices(this, input_graph, input_data, out)
    end if

    ! A supplied buffer is overwritten, never added to.
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine operator_apply

  !===================================================================!
  ! The edge landing, worked out.
  !===================================================================!

  subroutine apply_on_edges(this, input_graph, input_data, out)

    class(differential_operator), intent(in) :: this
    class(directed_graph), intent(in)                 :: input_graph
    class(graph_field), intent(in), optional :: input_data(:)
    type(field), intent(out)                 :: out

    real(dp), allocatable :: q(:), z(:), qc(:), zc(:), yc(:)
    integer               :: nv, ne, e, nc, c

    nv = input_graph % num_vertices()
    ne = input_graph % num_edges()

    call fetch_vertex_values(input_data, input_graph, nv, q, nc)

    out = field(this % name(), input_graph % edge_set(), input_graph % num_edges(), ncomp=max(nc, 1))

    allocate(z(ne * max(nc, 1)))
    z = 0.0_dp

    if (nc >= 1) then

       allocate(qc(nv), zc(ne))

       do c = 1, nc

          call gather_component(q, nc, c, qc)

          if (this % order <= 0) then
             call average_step(input_graph, &
                  & merge(1.0_dp, 0.0_dp, this % one_sided), .true., &
                  & this % coefficient, this % coefficients, &
                  & this % boundary_value, this % boundary_values, qc, zc)
          else if (this % order == 1) then
             call difference_step(input_graph, .true., &
                  & this % coefficient, this % coefficients, &
                  & this % spacing, this % spacings, &
                  & this % boundary_value, this % boundary_values, qc, zc)
          else
             ! The chain: the vertex result of order n - 1, then one
             ! difference step on top, coefficient already aboard.
             call run_vertex_chain(this % order - 1, input_graph, qc, &
                  & this % coefficient, this % coefficients, &
                  & this % spacing, this % spacings, &
                  & this % measure, this % measures, &
                  & this % boundary_value, this % boundary_values, yc)
             call difference_step(input_graph, .false., &
                  & this % coefficient, this % coefficients, &
                  & this % spacing, this % spacings, &
                  & this % boundary_value, this % boundary_values, yc, zc)
          end if

          call scatter_component(zc, nc, c, z)

       end do

    end if

    call out % set_real_vector(z)

  end subroutine apply_on_edges

  !===================================================================!
  ! THE VERTEX OPERATOR. The chain of steps, or - with `adjoint`
  ! true - the same chain transposed and walked backwards.
  !
  ! Handed an edge field, the chain enters at the incidence step: order 1 is
  ! then the divergence of the given field. Each component of the
  ! input makes the walk once.
  !===================================================================!

  subroutine apply_on_vertices(this, input_graph, input_data, out)

    class(differential_operator), intent(in) :: this
    class(directed_graph), intent(in)                 :: input_graph
    class(graph_field), intent(in), optional :: input_data(:)
    type(field), intent(out)                 :: out

    real(dp), allocatable :: q(:), z(:), y(:), qc(:), zc(:), yc(:), y2(:)
    real(dp), allocatable :: spent(:)   ! never allocated: the
                                        ! coefficient is applied once,
                                        ! before the incidence step, and
                                        ! empty array keeps the deeper
                                        ! chain from applying it again
    integer , allocatable :: indices(:)
    integer               :: nv, ne, v, e, nc, c

    nv = input_graph % num_vertices()
    ne = input_graph % num_edges()

    call fetch_vertex_values(input_data, input_graph, nv, q, nc)
    if (nc == 0) call fetch_edge_values(input_data, input_graph, ne, z, nc)

    out = field(this % name(), input_graph % vertex_set(), input_graph % num_vertices(), ncomp=max(nc, 1))

    allocate(y(nv * max(nc, 1)))
    y = 0.0_dp

    if (allocated(q) .and. size(q) > 0) then

       ! A vertex field: the whole chain per component, forward or
       ! reversed.
       allocate(qc(nv))

       do c = 1, nc

          call gather_component(q, nc, c, qc)

          if (this % order <= 0 .and. .not. this % adjoint) then
             allocate(yc(nv))
             do v = 1, nv
                yc(v) = coefficient_at(this % coefficient, this % coefficients, v) * qc(v)
             end do
          else if (this % adjoint) then
             call run_vertex_chain_reversed(this % order, input_graph, qc, &
                  & this % coefficient, this % coefficients, &
                  & this % spacing, this % spacings, &
                  & this % measure, this % measures, yc)
          else
             call run_vertex_chain(this % order, input_graph, qc, &
                  & this % coefficient, this % coefficients, &
                  & this % spacing, this % spacings, &
                  & this % measure, this % measures, &
                  & this % boundary_value, this % boundary_values, yc)
          end if

          call scatter_component(yc, nc, c, y)
          deallocate(yc)

       end do

    else if (allocated(z) .and. size(z) > 0) then

       ! An edge field: enter at the incidence step, per component.
       ! Order 1 is
       ! the divergence of the given samples; a higher order keeps
       ! walking.
       allocate(zc(ne), yc(nv))

       do c = 1, nc

          call gather_component(z, nc, c, zc)

          do e = 1, ne
             zc(e) = coefficient_at(this % coefficient, this % coefficients, e) * zc(e)
          end do
          call incidence_step(input_graph, this % measure, this % measures, zc, yc)

          if (this % order > 1) then
             call run_vertex_chain(this % order - 1, input_graph, yc, &
                  & 1.0_dp, spent, &
                  & this % spacing, this % spacings, &
                  & this % measure, this % measures, &
                  & this % boundary_value, this % boundary_values, y2)
             yc = y2
          end if

          call scatter_component(yc, nc, c, y)

       end do

    end if

    call out % set_real_vector(y)

  end subroutine apply_on_vertices

  !===================================================================!
  ! The forward chain, drawn once and run everywhere.
  !
  !    order 4:   q --G--> --D--> --G--> --D--> y
  !    order 3:   q --S--> --D--> --G--> --D--> y
  !    order 2:   q --G--> --D--> y
  !    order 1:   q --S--> --D--> y
  !    order 0:   y = c q
  !
  ! The innermost step - the first arrow - carries the coefficient;
  ! for odd orders it is the one-sided average, chosen by the sign of
  ! the coefficient, and for even orders the difference.
  !===================================================================!

  subroutine run_vertex_chain(order, g, q, c, cs, h, hs, m, ms, b, bs, y)

    integer     , intent(in)          :: order
    class(directed_graph), intent(in)          :: g
    real(dp)    , intent(in)          :: q(:)
    real(dp)    , intent(in)          :: c, h, m, b
    real(dp), allocatable, intent(in) :: cs(:), hs(:), ms(:), bs(:)
    real(dp), allocatable, intent(out):: y(:)

    real(dp), allocatable :: z(:)
    integer               :: nv, ne, k, v
    logical               :: innermost

    nv = g % num_vertices()
    ne = g % num_edges()

    allocate(y(nv), z(ne))

    if (order <= 0) then
       do v = 1, nv
          y(v) = coefficient_at(c, cs, v) * q(v)
       end do
       return
    end if

    y = q
    k = order
    innermost = .true.

    do while (k > 0)

       if (mod(k, 2) == 1) then
          ! The odd step: one-sided average, side by the sign of c.
          call average_step(g, 1.0_dp, innermost, c, cs, b, bs, y, z)
          k = k - 1
       else
          call difference_step(g, innermost, c, cs, h, hs, b, bs, y, z)
          k = k - 2
       end if

       call incidence_step(g, m, ms, z, y)
       innermost = .false.

    end do

  end subroutine run_vertex_chain

  !===================================================================!
  ! The reverse walk: the transposed steps, in reverse order.
  !
  !    forward, order 2:   q --G--> --D--> y
  !    adjoint, order 2:   q --D'--> --G'--> y
  !
  ! where D' is the transposed incidence step (a difference, sign
  ! flipped) and G' the transposed difference (an incidence step,
  ! sign flipped). The two sign
  ! flips cancel in pairs, which is why the order-2 operator with
  ! symmetric coefficients is its own adjoint - checked in the suite.
  !===================================================================!

  subroutine run_vertex_chain_reversed(order, g, q, c, cs, h, hs, m, ms, y)

    integer     , intent(in)          :: order
    class(directed_graph), intent(in)          :: g
    real(dp)    , intent(in)          :: q(:)
    real(dp)    , intent(in)          :: c, h, m
    real(dp), allocatable, intent(in) :: cs(:), hs(:), ms(:)
    real(dp), allocatable, intent(out):: y(:)

    real(dp), allocatable :: z(:), w(:)
    integer               :: nv, ne, pairs, i, v

    nv = g % num_vertices()
    ne = g % num_edges()

    allocate(y(nv), z(ne), w(nv))

    if (order <= 0) then
       do v = 1, nv
          y(v) = coefficient_at(c, cs, v) * q(v)
       end do
       return
    end if

    ! The forward chain, first step to last, is
    !
    !    even order:   G  D  G  D ... G  D      (coefficient on the
    !    odd order:    S  D  G  D ... G  D       first step)
    !
    ! so the reverse walk transposes each step and runs them last to
    ! first: the transposed incidence steps and differences make the
    ! pairs, and
    ! for an odd order the transposed average comes at the very end -
    ! carrying the coefficient, because its forward twin carried it.

    y = q

    if (mod(order, 2) == 0) then
       pairs = order / 2
    else
       pairs = (order - 1) / 2
    end if

    do i = 1, pairs
       call incidence_step_reversed(g, m, ms, y, z)
       call difference_step_reversed(g, &
            & mod(order, 2) == 0 .and. i == pairs, c, cs, h, hs, z, y)
    end do

    if (mod(order, 2) == 1) then
       call incidence_step_reversed(g, m, ms, y, z)
       call average_step_reversed(g, 1.0_dp, .true., c, cs, z, w)
       y = w
    end if

  end subroutine run_vertex_chain_reversed

  !===================================================================!
  ! Fetch the values once, before any loop, and report how many
  ! components ride in each entry. A wrong or missing field leaves a
  ! zero-length array and zero components, and the operator returns
  ! zeros rather than reading memory it was never given.
  !===================================================================!

  subroutine fetch_vertex_values(input_data, input_graph, nv, q, ncomp)

    class(graph_field), intent(in), optional :: input_data(:)
    class(directed_graph)     , intent(in)           :: input_graph
    integer          , intent(in)           :: nv
    real(dp), allocatable, intent(out)      :: q(:)
    integer          , intent(out)          :: ncomp

    type(set_graph) :: dom
    integer         :: n_dom

    ncomp = 0

    if (present(input_data)) then
       select type (state => input_data(1))
       class is (field)
          dom   = state % domain()
          n_dom = state % num_entries()
          ! Full coverage, by identity: this kernel indexes every
          ! vertex densely (routing is not admissibility).
          if (dom % same_as(input_graph % vertex_set())) then
             ncomp = max(state % num_components(), 1)
             call state % get_real_vector(q)
             if (size(q) == nv * ncomp) return
             ncomp = 0
          end if
       end select
    end if

    allocate(q(0))

  end subroutine fetch_vertex_values

  !===================================================================!
  ! The same fetch, for a field on edges.
  !===================================================================!

  subroutine fetch_edge_values(input_data, input_graph, ne, z, ncomp)

    class(graph_field), intent(in), optional :: input_data(:)
    class(directed_graph)     , intent(in)           :: input_graph
    integer          , intent(in)           :: ne
    real(dp), allocatable, intent(out)      :: z(:)
    integer          , intent(out)          :: ncomp

    type(set_graph) :: dom
    integer         :: n_dom

    ncomp = 0

    if (present(input_data)) then
       select type (state => input_data(1))
       class is (field)
          dom   = state % domain()
          n_dom = state % num_entries()
          if (dom % same_as(input_graph % edge_set())) then
             ncomp = max(state % num_components(), 1)
             call state % get_real_vector(z)
             if (size(z) == ne * ncomp) return
             ncomp = 0
          end if
       end select
    end if

    allocate(z(0))

  end subroutine fetch_edge_values


end module class_graph_differential_operator
