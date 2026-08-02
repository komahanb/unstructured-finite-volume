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
!    S   the average step        vertices -> edges
!
!        z_e = (q_i + q_j) / 2,  or one end, chosen by the
!                                sign of the coefficient
!
!    G   the difference step     vertices -> edges
!
!        z_e = (q_j - q_i) / h_e
!
!    D   the fold step           edges -> vertices
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
!    order 0        c q                the term itself
!    order 1        D S q             first derivative
!    order 2        D G q             second derivative
!    order 3        D G D S q         third
!    order 4        D G D G q         fourth
!    order n        one more step     reaches n rings of neighbours
!
! On the edge side, order 0 is S q, order 1 is G q (the slope along
! each edge), and order n applies G to the vertex result of n - 1.
! The coefficient is applied once, at the innermost step - the first
! moment the values land on edges - so a per-edge coefficient makes
! order 2 the operator div(k grad q) with k varying edge to edge.
!
! Exact where the formulas are exact: on a uniform chain a straight
! line has zero second derivative, x squared has second derivative
! two, and x to the fourth has fourth derivative twenty-four. The
! test suite checks those numbers.
!
!
!                       TWO FOLD DIRECTIONS
!
! The fold above is out minus in, which gives the derivatives their
! textbook signs. The balance folds in minus out, because a balance
! counts what a vertex gains. One sign apart; each is stated where
! it is used.
!
!
!                  THE ADJOINT IS THE REVERSE WALK
!
!    the walk                       the walk, reversed
!
!    (o)-->(o)-->(o)                (o)<--(o)<--(o)
!
! Reversing every edge swaps tail with head, so out becomes in:
! the transposed fold is a difference, the transposed difference is
! a fold, each with one sign flipped, and the transposed average
! folds - unsigned - onto the end it sampled. An adjoint runs the
! transposed steps in reverse order. Consequences, proved in the
! guide and checked in the suite: the order-2 operator is its own
! adjoint when the same per-edge coefficients sit on its steps, and
! the adjoint of the order-1 operator samples from the other end -
! the reversed-flow operator a sensitivity solve consumes.
!
!
!                     CURL, WITH NOTHING NEW
!
! Around a loop of edges there is one more derivative: add the edge
! values as the loop is walked. Which edge borders which face is
! itself a graph - one vertex per edge, one per face, a connection
! where an edge borders a face, coefficient +1 or -1 by direction -
! and the fold step on that graph is the curl. The suite walks a
! square: a difference field sums to zero around it, a circulating
! field does not.
!
!
!                        WHAT IS NOT HERE
!
! No physical names and no physical signs. A transport model uses
! order 1 with its velocity and order 2 with its conductivity, and
! writes its own sign conventions in its own class. This layer
! computes derivatives.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_graph_differential_operator

  use iso_fortran_env     , only : dp => REAL64
  use abstract_graph_types, only : graph_edge_field_operation
  use abstract_graph_types, only : graph_vertex_field_operation
  use abstract_graph_types, only : graph, graph_data
  use abstract_graph_types, only : graph_edge_field, graph_vertex_field
  use abstract_graph_types, only : graph_edge_support, graph_vertex_support
  use class_graph_support , only : vertex_support, edge_support
  use class_graph_field   , only : vertex_field, edge_field

  implicit none

  private
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
  ! The operator that returns an edge field: the derivative of order
  ! n, sampled on the edges.
  !
  !    order 0    S q     the average on each edge
  !    order 1    G q     the slope along each edge
  !    order n    G of the vertex result of n - 1
  !===================================================================!

  type, extends(graph_edge_field_operation) :: edge_differential_operator

     integer :: order = 1

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

     procedure :: name    => edge_differential_operator_name
     procedure :: support => edge_differential_operator_support
     procedure :: apply   => edge_differential_operator_apply

  end type edge_differential_operator

  !===================================================================!
  ! The operator that returns a vertex field: the derivative of order
  ! n at each vertex.
  !
  !    order 0    c q          the term itself
  !    order 1    D S q        first derivative
  !    order 2    D G q        second derivative
  !    order n    keep going
  !
  ! Handed an EDGE field instead of a vertex field, the chain enters
  ! at the fold: order 1 on an edge field is the divergence of that
  ! field.
  !
  ! With `adjoint` true the operator runs the reverse walk: the
  ! transposed steps, in reverse order.
  !===================================================================!

  type, extends(graph_vertex_field_operation) :: vertex_differential_operator

     integer :: order = 2

     logical :: adjoint = .false.

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

     procedure :: name    => vertex_differential_operator_name
     procedure :: support => vertex_differential_operator_support
     procedure :: apply   => vertex_differential_operator_apply

  end type vertex_differential_operator

  interface edge_differential_operator
     module procedure create_edge_operator
  end interface edge_differential_operator

  interface vertex_differential_operator
     module procedure create_vertex_operator
  end interface vertex_differential_operator

contains

  !===================================================================!
  ! Constructors. Order first; every parameter optional; each array
  ! wins over its scalar when given.
  !===================================================================!

  pure type(edge_differential_operator) function create_edge_operator &
       & (order, coefficient, coefficients, spacing, spacings, &
       &  measure, measures, boundary_value, boundary_values, label) result(this)

    integer         , intent(in)           :: order
    real(dp)        , intent(in), optional :: coefficient
    real(dp)        , intent(in), optional :: coefficients(:)
    real(dp)        , intent(in), optional :: spacing
    real(dp)        , intent(in), optional :: spacings(:)
    real(dp)        , intent(in), optional :: measure
    real(dp)        , intent(in), optional :: measures(:)
    real(dp)        , intent(in), optional :: boundary_value
    real(dp)        , intent(in), optional :: boundary_values(:)
    character(len=*), intent(in), optional :: label

    this % order = max(order, 0)

    if (present(coefficient))     this % coefficient    = coefficient
    if (present(coefficients))    allocate(this % coefficients, source=coefficients)
    if (present(spacing))         this % spacing        = spacing
    if (present(spacings))        allocate(this % spacings, source=spacings)
    if (present(measure))         this % measure        = measure
    if (present(measures))        allocate(this % measures, source=measures)
    if (present(boundary_value))  this % boundary_value = boundary_value
    if (present(boundary_values)) allocate(this % boundary_values, source=boundary_values)
    if (present(label))           this % label          = label

  end function create_edge_operator

  pure type(vertex_differential_operator) function create_vertex_operator &
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

    this % order = max(order, 0)

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

  end function create_vertex_operator

  !===================================================================!
  ! The named layer. Four operators every equation reaches for, as
  ! plain functions, so a call site reads as the calculus:
  !
  !    gradient(...)        the slope along each edge
  !    interpolation(...)   the average onto each edge
  !    divergence(...)      the fold of an edge field onto vertices
  !    laplacian(...)       the second derivative at each vertex
  !===================================================================!

  pure type(edge_differential_operator) function gradient(coefficient, coefficients, &
       & spacing, spacings, boundary_value, boundary_values) result(this)

    real(dp), intent(in), optional :: coefficient, coefficients(:)
    real(dp), intent(in), optional :: spacing, spacings(:)
    real(dp), intent(in), optional :: boundary_value, boundary_values(:)

    this = edge_differential_operator(order=1, coefficient=coefficient, &
         & coefficients=coefficients, spacing=spacing, spacings=spacings, &
         & boundary_value=boundary_value, boundary_values=boundary_values, &
         & label='gradient')

  end function gradient

  pure type(edge_differential_operator) function interpolation(coefficient, coefficients, &
       & boundary_value, boundary_values) result(this)

    real(dp), intent(in), optional :: coefficient, coefficients(:)
    real(dp), intent(in), optional :: boundary_value, boundary_values(:)

    this = edge_differential_operator(order=0, coefficient=coefficient, &
         & coefficients=coefficients, boundary_value=boundary_value, &
         & boundary_values=boundary_values, label='interpolation')

  end function interpolation

  pure type(vertex_differential_operator) function divergence(coefficient, coefficients, &
       & measure, measures) result(this)

    real(dp), intent(in), optional :: coefficient, coefficients(:)
    real(dp), intent(in), optional :: measure, measures(:)

    this = vertex_differential_operator(order=1, coefficient=coefficient, &
         & coefficients=coefficients, measure=measure, measures=measures, &
         & label='divergence')

  end function divergence

  pure type(vertex_differential_operator) function laplacian(coefficient, coefficients, &
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

  pure function edge_differential_operator_name(this) result(name)

    class(edge_differential_operator), intent(in) :: this
    character(len=:), allocatable                 :: name

    if (allocated(this % label)) then
       name = this % label
    else
       name = order_name(this % order)
    end if

  end function edge_differential_operator_name

  pure function vertex_differential_operator_name(this) result(name)

    class(vertex_differential_operator), intent(in) :: this
    character(len=:), allocatable                   :: name

    if (allocated(this % label)) then
       name = this % label
    else
       name = order_name(this % order)
    end if

  end function vertex_differential_operator_name

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

  subroutine edge_differential_operator_support(this, input_graph, support)

    class(edge_differential_operator), intent(in)       :: this
    class(graph), intent(in)                            :: input_graph
    class(graph_edge_support), allocatable, intent(out) :: support

    call input_graph % all_edges(support)

  end subroutine edge_differential_operator_support

  subroutine vertex_differential_operator_support(this, input_graph, support)

    class(vertex_differential_operator), intent(in)       :: this
    class(graph), intent(in)                              :: input_graph
    class(graph_vertex_support), allocatable, intent(out) :: support

    call input_graph % all_vertices(support)

  end subroutine vertex_differential_operator_support

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

    class(graph), intent(in)          :: g
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

  pure subroutine average_step_reversed(g, one_sided_by, with_c, op_c, op_cs, z, y)

    class(graph), intent(in)          :: g
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
  !            h_e                       a fold, one sign flipped
  !===================================================================!

  pure subroutine difference_step(g, with_c, op_c, op_cs, op_h, op_hs, op_b, op_bs, q, z)

    class(graph), intent(in)          :: g
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

  pure subroutine difference_step_reversed(g, with_c, op_c, op_cs, op_h, op_hs, z, y)

    class(graph), intent(in)          :: g
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
  ! THE FOLD STEP, forward and reversed.
  !
  !    forward                              reversed
  !
  !         out - in                        q_tail   q_head
  !    y_v = --------                z_e =  ------ - ------
  !            m_v                          m_tail   m_head
  !
  ! Out minus in gives the derivatives their textbook signs; the
  ! balance folds in minus out because it counts what a vertex
  ! gains. An edge with no head contributes to its tail alone.
  !===================================================================!

  pure subroutine fold_step(g, op_m, op_ms, z, y)

    class(graph), intent(in)          :: g
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

  end subroutine fold_step

  pure subroutine fold_step_reversed(g, op_m, op_ms, q, z)

    class(graph), intent(in)          :: g
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

  end subroutine fold_step_reversed

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

  subroutine edge_differential_operator_apply(this, input_graph, input_data, output)

    class(edge_differential_operator), intent(in)       :: this
    class(graph), intent(in)                            :: input_graph
    class(graph_data), intent(in), optional             :: input_data(:)
    class(graph_edge_field), allocatable, intent(inout) :: output

    type(edge_field)      :: out
    type(edge_support)    :: on
    real(dp), allocatable :: q(:), z(:), qc(:), zc(:), yc(:)
    integer , allocatable :: indices(:)
    integer               :: nv, ne, e, nc, c

    nv = input_graph % num_vertices()
    ne = input_graph % num_edges()

    call fetch_vertex_values(input_data, nv, q, nc)

    allocate(indices(ne))
    do e = 1, ne
       indices(e) = e
    end do
    on  = edge_support(indices)
    out = edge_field(this % name(), on, ncomp=max(nc, 1))

    allocate(z(ne * max(nc, 1)))
    z = 0.0_dp

    if (nc >= 1) then

       allocate(qc(nv), zc(ne))

       do c = 1, nc

          call gather_component(q, nc, c, qc)

          if (this % order <= 0) then
             call average_step(input_graph, 0.0_dp, .true., &
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

    ! A supplied buffer is overwritten, never added to.
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine edge_differential_operator_apply

  !===================================================================!
  ! THE VERTEX OPERATOR. The chain of steps, or - with `adjoint`
  ! true - the same chain transposed and walked backwards.
  !
  ! Handed an edge field, the chain enters at the fold: order 1 is
  ! then the divergence of the given field. Each component of the
  ! input makes the walk once.
  !===================================================================!

  subroutine vertex_differential_operator_apply(this, input_graph, input_data, output)

    class(vertex_differential_operator), intent(in)       :: this
    class(graph), intent(in)                              :: input_graph
    class(graph_data), intent(in), optional               :: input_data(:)
    class(graph_vertex_field), allocatable, intent(inout) :: output

    type(vertex_field)    :: out
    type(vertex_support)  :: on
    real(dp), allocatable :: q(:), z(:), y(:), qc(:), zc(:), yc(:), y2(:)
    real(dp), allocatable :: spent(:)   ! never allocated: the
                                        ! coefficient is applied once,
                                        ! before the fold, and this
                                        ! empty array keeps the deeper
                                        ! chain from applying it again
    integer , allocatable :: indices(:)
    integer               :: nv, ne, v, e, nc, c

    nv = input_graph % num_vertices()
    ne = input_graph % num_edges()

    call fetch_vertex_values(input_data, nv, q, nc)
    if (nc == 0) call fetch_edge_values(input_data, ne, z, nc)

    allocate(indices(nv))
    do v = 1, nv
       indices(v) = v
    end do
    on  = vertex_support(indices)
    out = vertex_field(this % name(), on, ncomp=max(nc, 1))

    allocate(y(nv * max(nc, 1)))
    y = 0.0_dp

    if (allocated(q) .and. size(q) > 0) then

       ! A vertex field: the full chain per component, forward or
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

       ! An edge field: enter at the fold, per component. Order 1 is
       ! the divergence of the given samples; a higher order keeps
       ! walking.
       allocate(zc(ne), yc(nv))

       do c = 1, nc

          call gather_component(z, nc, c, zc)

          do e = 1, ne
             zc(e) = coefficient_at(this % coefficient, this % coefficients, e) * zc(e)
          end do
          call fold_step(input_graph, this % measure, this % measures, zc, yc)

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

    ! A supplied buffer is overwritten, never added to.
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine vertex_differential_operator_apply

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
    class(graph), intent(in)          :: g
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

       call fold_step(g, m, ms, z, y)
       innermost = .false.

    end do

  end subroutine run_vertex_chain

  !===================================================================!
  ! The reverse walk: the transposed steps, in reverse order.
  !
  !    forward, order 2:   q --G--> --D--> y
  !    adjoint, order 2:   q --D'--> --G'--> y
  !
  ! where D' is the transposed fold (a difference, sign flipped) and
  ! G' the transposed difference (a fold, sign flipped). The two sign
  ! flips cancel in pairs, which is why the order-2 operator with
  ! symmetric coefficients is its own adjoint - checked in the suite.
  !===================================================================!

  subroutine run_vertex_chain_reversed(order, g, q, c, cs, h, hs, m, ms, y)

    integer     , intent(in)          :: order
    class(graph), intent(in)          :: g
    real(dp)    , intent(in)          :: q(:)
    real(dp)    , intent(in)          :: c, h, m
    real(dp), allocatable, intent(in) :: cs(:), hs(:), ms(:)
    real(dp), allocatable, intent(out):: y(:)

    real(dp), allocatable :: z(:), w(:)
    integer               :: nv, ne, k, v
    logical               :: outermost_pending

    nv = g % num_vertices()
    ne = g % num_edges()

    allocate(y(nv), z(ne), w(nv))

    if (order <= 0) then
       do v = 1, nv
          y(v) = coefficient_at(c, cs, v) * q(v)
       end do
       return
    end if

    y = q
    k = order

    ! The forward chain ends with a fold, so the reverse walk begins
    ! with the transposed fold, and the step that carried the
    ! coefficient - the forward chain's first - is met last.
    do while (k > 0)

       call fold_step_reversed(g, m, ms, y, z)

       outermost_pending = (k == 1 .or. k == 2)

       if (mod(k, 2) == 1) then
          call average_step_reversed(g, 1.0_dp, outermost_pending, c, cs, z, w)
          y = w
          k = k - 1
       else
          call difference_step_reversed(g, outermost_pending, c, cs, h, hs, z, y)
          k = k - 2
       end if

    end do

  end subroutine run_vertex_chain_reversed

  !===================================================================!
  ! Fetch the values once, before any loop. A wrong or missing field
  !===================================================================!
  ! Fetch the values once, before any loop, and report how many
  ! components ride in each entry. A wrong or missing field leaves a
  ! zero-length array and zero components, and the operator returns
  ! zeros rather than reading memory it was never given.
  !===================================================================!

  subroutine fetch_vertex_values(input_data, nv, q, ncomp)

    class(graph_data), intent(in), optional :: input_data(:)
    integer          , intent(in)           :: nv
    real(dp), allocatable, intent(out)      :: q(:)
    integer          , intent(out)          :: ncomp

    ncomp = 0

    if (present(input_data)) then
       select type (state => input_data(1))
       class is (vertex_field)
          ncomp = max(state % num_components(), 1)
          call state % get_real_vector(q)
          if (size(q) == nv * ncomp) return
          ncomp = 0
       end select
    end if

    allocate(q(0))

  end subroutine fetch_vertex_values

  subroutine fetch_edge_values(input_data, ne, z, ncomp)

    class(graph_data), intent(in), optional :: input_data(:)
    integer          , intent(in)           :: ne
    real(dp), allocatable, intent(out)      :: z(:)
    integer          , intent(out)          :: ncomp

    ncomp = 0

    if (present(input_data)) then
       select type (state => input_data(1))
       class is (edge_field)
          ncomp = max(state % num_components(), 1)
          call state % get_real_vector(z)
          if (size(z) == ne * ncomp) return
          ncomp = 0
       end select
    end if

    allocate(z(0))

  end subroutine fetch_edge_values


end module class_graph_differential_operator
