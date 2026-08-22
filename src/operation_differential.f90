!=====================================================================!
! Differential operators on a graph, compiled onto stencils.
!
! Values live on vertices or on edges; three elementary steps move
! them, and each step is an affine sparse map - a matrix plus a
! constant vector, the constant carrying what boundary values leave
! behind:
!
!    S   the average step        edges x vertices
!
!        z_e = (q_i + q_j)/2, or one end, chosen by the sign of
!        the coefficient when the step carries it. The half is the
!        mean of a straight-line q along the edge. An edge with no
!        head reads the stored boundary value in place of the
!        head, which lands in the constant.
!
!    G   the difference step     edges x vertices
!
!        z_e = c (q_j - q_i) / h_e
!
!    D   the incidence step      vertices x edges
!
!        y_v = (out - in) / m_v : the values on edges leaving v,
!        minus the values on edges entering v, over the measure of
!        v. An edge with no head contributes to its tail alone.
!
! The operator of order n is the composition of those maps by
! parity,
!
!    vertex(0) = C q            vertex(n) = D vertex-side(n-1)
!    edge(0)   = S q            edge(n)   = G vertex(n-1)
!
!    vertex(2k)   = (D G)^k
!    vertex(2k+1) = (D G)^k D S
!
! and a composition of affine maps is one affine map,
!
!    (A2, k2) o (A1, k1) = (A2 A1, A2 k1 + k2),
!
! computed once, by the stencil's composition: every step is a
! stencil and so is the chain. The coefficient is carried by the innermost
! step, so a per-edge coefficient makes order 2 the operator
! div(k grad q).
!
! THE ADJOINT IS THE TRANSPOSE. With `adjoint` true the operator
! applies the transpose of the composed matrix - rows and columns
! swapped, the constant dropped, because the adjoint acts on the
! linear part. No reversed step kernels exist: (C B A)^T =
! A^T B^T C^T is an identity of the composition, not code.
!
! THE STENCIL DOOR. The composed map is a stencil on either
! landing, its row domain the landing's set, and stencil_of returns
! it - the same triples, the same constant - so a minimizer can
! attach the compiled matrix directly. apply is the stencil's own
! apply, once per component.
!
! Handed an EDGE field on the vertex landing, the composition
! enters at the incidence step: order 1 is then the divergence of
! that field. Handed no field, or a field on the wrong domain, the
! operator returns zeros rather than reading memory it was never
! given.
!
! Each step consults an edge's two ends, so order n reaches
! exactly n rings of neighbours; the composed pattern states that
! reach explicitly. Exactness is claimed on the uniform chain,
! where the discrete formulas coincide with calculus, and the test
! suite checks those numbers. No physical names and no physical
! signs live here; models state their own.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_differential

  use iso_fortran_env    , only : dp => REAL64
  use operation_action, only : operation
  use view_directed, only : directed_graph
  use field_calculus, only : field
  use graph_fractal      , only : graph
  use view_directed     , only : SIDE_VERTEX, SIDE_EDGE
  use field_stored  , only : stored_field
  use operation_stencil, only : stencil, compose

  implicit none

  private
  public :: differential_operator
  public :: edge_derivative
  public :: vertex_derivative
  public :: gradient, interpolation, divergence, laplacian
  public :: stencil_of

  !===================================================================!
  ! The shared parameters. Each is one number for the uniform case,
  ! with an optional per-entity array that wins when allocated:
  !
  !    coefficient   c     applied at the innermost step (per edge;
  !                        per vertex at order 0 on the vertex side)
  !    spacing       h_e   the length of an edge, in G
  !    measure       m_v   the size of a vertex, in D
  !    boundary      b_e   the value standing in for the missing
  !                        end of an edge with no head, in S and G
  !===================================================================!

  type, extends(operation) :: differential_operator

     integer :: landing = SIDE_VERTEX
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

  type(differential_operator) function edge_derivative &
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

    this % landing = SIDE_EDGE
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

    ! one argument: the field differentiated
    call this % declare_arguments(1)

  end function edge_derivative

  !===================================================================!
  ! The same constructor, for the vertex landing, plus the adjoint
  ! flag: raised, the operator applies its transpose.
  !===================================================================!

  type(differential_operator) function vertex_derivative &
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

    this % landing = SIDE_VERTEX
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

    call this % declare_arguments(1)

  end function vertex_derivative

  !===================================================================!
  ! The named layer. Four operators every equation reaches for, as
  ! plain functions, so a call site reads as the calculus:
  !
  !    gradient(...)        the slope along each edge
  !    interpolation(...)   the average onto each edge
  !    divergence(...)      the incidence step on an edge field
  !    laplacian(...)       the second derivative at each vertex
  !===================================================================!

  type(differential_operator) function gradient(coefficient, coefficients, &
       & spacing, spacings, boundary_value, boundary_values) result(this)

    real(dp), intent(in), optional :: coefficient, coefficients(:)
    real(dp), intent(in), optional :: spacing, spacings(:)
    real(dp), intent(in), optional :: boundary_value, boundary_values(:)

    this = edge_derivative(order=1, coefficient=coefficient, &
         & coefficients=coefficients, spacing=spacing, spacings=spacings, &
         & boundary_value=boundary_value, boundary_values=boundary_values, &
         & label='gradient')

  end function gradient

  type(differential_operator) function interpolation(coefficient, coefficients, &
       & boundary_value, boundary_values) result(this)

    real(dp), intent(in), optional :: coefficient, coefficients(:)
    real(dp), intent(in), optional :: boundary_value, boundary_values(:)

    this = edge_derivative(order=0, coefficient=coefficient, &
         & coefficients=coefficients, boundary_value=boundary_value, &
         & boundary_values=boundary_values, label='interpolation')

  end function interpolation

  type(differential_operator) function divergence(coefficient, coefficients, &
       & measure, measures) result(this)

    real(dp), intent(in), optional :: coefficient, coefficients(:)
    real(dp), intent(in), optional :: measure, measures(:)

    this = vertex_derivative(order=1, coefficient=coefficient, &
         & coefficients=coefficients, measure=measure, measures=measures, &
         & label='divergence')

  end function divergence

  type(differential_operator) function laplacian(coefficient, coefficients, &
       & spacing, spacings, measure, measures, boundary_value, boundary_values) result(this)

    real(dp), intent(in), optional :: coefficient, coefficients(:)
    real(dp), intent(in), optional :: spacing, spacings(:)
    real(dp), intent(in), optional :: measure, measures(:)
    real(dp), intent(in), optional :: boundary_value, boundary_values(:)

    this = vertex_derivative(order=2, coefficient=coefficient, &
         & coefficients=coefficients, spacing=spacing, spacings=spacings, &
         & measure=measure, measures=measures, boundary_value=boundary_value, &
         & boundary_values=boundary_values, label='laplacian')

  end function laplacian

  !===================================================================!
  ! Names. A named operator answers with its name; any other with
  ! its order.
  !===================================================================!

  pure function operator_name(this) result(name)

    class(differential_operator), intent(in) :: this
    character(len=:), allocatable            :: name

    character(len=12) :: digits

    if (allocated(this % label)) then
       name = this % label
    else
       write(digits, '(i0)') this % order
       name = 'derivative of order ' // trim(digits)
    end if

  end function operator_name

  !===================================================================!
  ! Domains: every edge, or every vertex, by the landing. An
  ! operator aimed at a subset is a second instance handed that
  ! subset's graph.
  !===================================================================!

  subroutine operator_domain(this, input_graph, domain, num_entries)

    class(differential_operator), intent(in) :: this
    class(directed_graph), intent(in)                 :: input_graph
    type(graph), intent(out) :: domain
    integer        , intent(out) :: num_entries

    if (this % landing == SIDE_EDGE) then
       domain   = input_graph % edge_set()
       num_entries = input_graph % num_edges()
    else
       domain   = input_graph % vertex_set()
       num_entries = input_graph % num_vertices()
    end if

  end subroutine operator_domain

  !===================================================================!
  ! Parameter lookup: the array wins when allocated.
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
  ! THE THREE STEP MAPS, each one stencil: edges x vertices for the
  ! average and the difference, vertices x edges for the incidence.
  ! Each accumulates its triples and hands them to the stencil with
  ! both domains, so a composed chain lands on the set its outermost
  ! step names and reads the set its innermost step reads.
  !
  ! The average step S, edges x vertices. When it carries the
  ! coefficient and the one-sided choice is on, the sign of the
  ! per-edge coefficient picks the end: positive samples the tail,
  ! negative the head; otherwise both ends, half each. A missing
  ! head reads the boundary value into the constant.
  !===================================================================!

  function average_map(g, one_sided_by, with_c, op_c, op_cs, &
       & op_b, op_bs) result(a)

    class(directed_graph), intent(in) :: g
    real(dp)    , intent(in)          :: one_sided_by
    logical     , intent(in)          :: with_c
    real(dp)    , intent(in)          :: op_c
    real(dp), allocatable, intent(in) :: op_cs(:)
    real(dp)    , intent(in)          :: op_b
    real(dp), allocatable, intent(in) :: op_bs(:)
    type(stencil)                     :: a

    integer , allocatable :: rows(:), cols(:)
    real(dp), allocatable :: weights(:), constants(:)
    real(dp) :: c, pick, b
    integer  :: e, n, ne

    ne = g % num_edges()
    allocate(rows(2 * ne), cols(2 * ne), weights(2 * ne), constants(ne))
    constants = 0.0_dp

    n = 0
    do e = 1, ne

       c = 1.0_dp
       if (with_c) c = coefficient_at(op_c, op_cs, e)

       pick = one_sided_by
       if (with_c) pick = sign(1.0_dp, c) * merge(1.0_dp, 0.0_dp, one_sided_by /= 0.0_dp)

       b = coefficient_at(op_b, op_bs, e)

       if (pick > 0.0_dp) then
          ! the tail end
          call put(rows, cols, weights, n, e, g % edge_tail(e), c)
       else if (pick < 0.0_dp) then
          ! the head end, or the boundary value
          if (g % edge_has_head(e)) then
             call put(rows, cols, weights, n, e, g % edge_head(e), c)
          else
             constants(e) = c * b
          end if
       else
          ! both ends, evenly
          call put(rows, cols, weights, n, e, g % edge_tail(e), c * 0.5_dp)
          if (g % edge_has_head(e)) then
             call put(rows, cols, weights, n, e, g % edge_head(e), c * 0.5_dp)
          else
             constants(e) = c * 0.5_dp * b
          end if
       end if

    end do

    a = stencil(rows(1:n), cols(1:n), weights(1:n), constants, &
         & label='average', num_columns=g % num_vertices(), &
         & row_domain=g % edge_set(), column_domain=g % vertex_set())

  end function average_map

  !===================================================================!
  ! The difference step G, edges x vertices:
  !
  !    z_e = c (q_head - q_tail) / h_e,
  !
  ! the head replaced by the boundary value - into the constant -
  ! on an edge with no head.
  !===================================================================!

  function difference_map(g, with_c, op_c, op_cs, op_h, op_hs, &
       & op_b, op_bs) result(a)

    class(directed_graph), intent(in) :: g
    logical     , intent(in)          :: with_c
    real(dp)    , intent(in)          :: op_c
    real(dp), allocatable, intent(in) :: op_cs(:)
    real(dp)    , intent(in)          :: op_h
    real(dp), allocatable, intent(in) :: op_hs(:)
    real(dp)    , intent(in)          :: op_b
    real(dp), allocatable, intent(in) :: op_bs(:)
    type(stencil)                     :: a

    integer , allocatable :: rows(:), cols(:)
    real(dp), allocatable :: weights(:), constants(:)
    real(dp) :: c, w
    integer  :: e, n, ne

    ne = g % num_edges()
    allocate(rows(2 * ne), cols(2 * ne), weights(2 * ne), constants(ne))
    constants = 0.0_dp

    n = 0
    do e = 1, ne

       c = 1.0_dp
       if (with_c) c = coefficient_at(op_c, op_cs, e)
       w = c / coefficient_at(op_h, op_hs, e)

       call put(rows, cols, weights, n, e, g % edge_tail(e), -w)

       if (g % edge_has_head(e)) then
          call put(rows, cols, weights, n, e, g % edge_head(e), w)
       else
          constants(e) = w * coefficient_at(op_b, op_bs, e)
       end if

    end do

    a = stencil(rows(1:n), cols(1:n), weights(1:n), constants, &
         & label='difference', num_columns=g % num_vertices(), &
         & row_domain=g % edge_set(), column_domain=g % vertex_set())

  end function difference_map

  !===================================================================!
  ! The incidence step D, vertices x edges:
  !
  !    y_v = (out - in) / m_v.
  !
  ! Out minus in gives the derivatives their textbook signs; an
  ! edge with no head contributes to its tail alone.
  !===================================================================!

  function incidence_map(g, op_m, op_ms) result(a)

    class(directed_graph), intent(in) :: g
    real(dp)    , intent(in)          :: op_m
    real(dp), allocatable, intent(in) :: op_ms(:)
    type(stencil)                     :: a

    integer , allocatable :: rows(:), cols(:)
    real(dp), allocatable :: weights(:), constants(:)
    integer :: e, t, h, n, ne, nv

    ne = g % num_edges()
    nv = g % num_vertices()
    allocate(rows(2 * ne), cols(2 * ne), weights(2 * ne), constants(nv))
    constants = 0.0_dp

    n = 0
    do e = 1, ne

       t = g % edge_tail(e)
       call put(rows, cols, weights, n, t, e, 1.0_dp / coefficient_at(op_m, op_ms, t))

       if (g % edge_has_head(e)) then
          h = g % edge_head(e)
          call put(rows, cols, weights, n, h, e, -1.0_dp / coefficient_at(op_m, op_ms, h))
       end if

    end do

    a = stencil(rows(1:n), cols(1:n), weights(1:n), constants, &
         & label='incidence', num_columns=ne, &
         & row_domain=g % vertex_set(), column_domain=g % edge_set())

  end function incidence_map

  !===================================================================!
  ! The diagonal map on one set, n x n: y_i = c_i q_i. Order 0 on
  ! the vertex landing, and the per-edge coefficient of an
  ! edge-field entry.
  !===================================================================!

  function diagonal_map(on, n, op_c, op_cs) result(a)

    type(graph) , intent(in)          :: on
    integer     , intent(in)          :: n
    real(dp)    , intent(in)          :: op_c
    real(dp), allocatable, intent(in) :: op_cs(:)
    type(stencil)                     :: a

    integer , allocatable :: rows(:)
    real(dp), allocatable :: weights(:), constants(:)
    integer :: i

    allocate(rows(n), weights(n), constants(n))
    constants = 0.0_dp

    do i = 1, n
       rows(i)    = i
       weights(i) = coefficient_at(op_c, op_cs, i)
    end do

    a = stencil(rows, rows, weights, constants, label='diagonal', &
         & row_domain=on, column_domain=on)

  end function diagonal_map

  !===================================================================!
  ! Triple bookkeeping: append one entry.
  !===================================================================!

  pure subroutine put(rows, cols, weights, n, r, c, w)

    integer , intent(inout) :: rows(:), cols(:)
    real(dp), intent(inout) :: weights(:)
    integer , intent(inout) :: n
    integer , intent(in)    :: r, c
    real(dp), intent(in)    :: w

    n = n + 1
    rows(n)    = r
    cols(n)    = c
    weights(n) = w

  end subroutine put

  !===================================================================!
  ! The parity chain, stated as the law reads: vertex(0) = C,
  ! vertex(n) = D edge(n-1); edge(0) = S, edge(1) = G, edge(n) =
  ! G vertex(n-1) with G bare. The innermost step carries the
  ! coefficient: the average one-sided by the coefficient's sign
  ! (one_sided_by = 1) wherever a vertex chain reaches it, and by
  ! the operator's own flag only at compiled_map's direct edge
  ! order 0. The composition associates as D (G (D (... S_c))),
  ! innermost first. The incidence map is formed at each vertex
  ! level - once for every order the tree uses.
  !===================================================================!

  recursive function vertex_chain(order, g, c, cs, h, hs, m, ms, b, bs) &
       & result(a)

    integer     , intent(in)          :: order
    class(directed_graph), intent(in) :: g
    real(dp)    , intent(in)          :: c, h, m, b
    real(dp), allocatable, intent(in) :: cs(:), hs(:), ms(:), bs(:)
    type(stencil)                     :: a

    if (order <= 0) then
       a = diagonal_map(g % vertex_set(), g % num_vertices(), c, cs)
    else
       a = compose(incidence_map(g, m, ms), &
            & edge_chain(order - 1, g, 1.0_dp, c, cs, h, hs, m, ms, b, bs))
    end if

  end function vertex_chain

  recursive function edge_chain(order, g, one_sided_by, c, cs, h, hs, &
       & m, ms, b, bs) result(a)

    integer     , intent(in)          :: order
    class(directed_graph), intent(in) :: g
    real(dp)    , intent(in)          :: one_sided_by, c, h, m, b
    real(dp), allocatable, intent(in) :: cs(:), hs(:), ms(:), bs(:)
    type(stencil)                     :: a

    if (order <= 0) then
       a = average_map(g, one_sided_by, .true., c, cs, b, bs)
    else if (order == 1) then
       a = difference_map(g, .true., c, cs, h, hs, b, bs)
    else
       a = compose(difference_map(g, .false., c, cs, h, hs, b, bs), &
            & vertex_chain(order - 1, g, c, cs, h, hs, m, ms, b, bs))
    end if

  end function edge_chain

  !===================================================================!
  ! The whole operator as one affine map, by landing and entry:
  !
  !    edge landing            the edge chain of the order: S, G,
  !                            or G (bare) after the vertex chain
  !                            of n - 1
  !    vertex landing          the vertex chain, transposed when
  !                            adjoint
  !    vertex landing, entry   diag(c) then D, then the bare vertex
  !    on an edge field        chain of n - 1: order 1 is the
  !                            divergence of the given field
  !===================================================================!

  function compiled_map(this, g, enters_on_edges) result(a)

    class(differential_operator), intent(in) :: this
    class(directed_graph)       , intent(in) :: g
    logical                     , intent(in) :: enters_on_edges
    type(stencil)                            :: a

    real(dp), allocatable :: spent(:)   ! never allocated: the
                                        ! coefficient is applied once

    if (this % landing == SIDE_EDGE) then

       ! the operator's one-sided flag reaches only its own order 0;
       ! every average inside a chain is one-sided by its sign
       a = edge_chain(this % order, g, merge(1.0_dp, 0.0_dp, this % one_sided), &
            & this % coefficient, this % coefficients, &
            & this % spacing, this % spacings, &
            & this % measure, this % measures, &
            & this % boundary_value, this % boundary_values)

    else if (enters_on_edges) then

       a = compose(incidence_map(g, this % measure, this % measures), &
            & diagonal_map(g % edge_set(), g % num_edges(), &
            & this % coefficient, this % coefficients))
       if (this % order > 1) then
          a = compose(vertex_chain(this % order - 1, g, &
               & 1.0_dp, spent, &
               & this % spacing, this % spacings, &
               & this % measure, this % measures, &
               & this % boundary_value, this % boundary_values), a)
       end if

    else

       a = vertex_chain(this % order, g, &
            & this % coefficient, this % coefficients, &
            & this % spacing, this % spacings, &
            & this % measure, this % measures, &
            & this % boundary_value, this % boundary_values)
       if (this % adjoint) a = a % transpose()

    end if

  end function compiled_map

  !===================================================================!
  ! THE STENCIL DOOR: the compiled operator as a stencil, on either
  ! landing, named after the operator. on_edge_field compiles the
  ! vertex landing's entry on an edge field - apply's choice when
  ! handed one - so the bare incidence step is reachable.
  !===================================================================!

  function stencil_of(operator, input_graph, on_edge_field) result(compiled)

    type(differential_operator), intent(in) :: operator
    class(directed_graph)      , intent(in) :: input_graph
    logical, intent(in), optional           :: on_edge_field
    type(stencil)                  :: compiled

    logical :: enters_on_edges

    enters_on_edges = .false.
    if (present(on_edge_field)) enters_on_edges = on_edge_field &
         & .and. operator % landing == SIDE_VERTEX

    compiled = compiled_map(operator, input_graph, enters_on_edges)
    compiled % label = operator % name()

  end function stencil_of

  !===================================================================!
  ! COMPONENTS. A field may carry several values per entry,
  ! interleaved entry-fastest:
  !
  !      flat((entry - 1) * num_components + component)
  !
  ! The map is compiled once; each component is gathered, swept,
  ! and scattered back. The parameters are shared by all
  ! components; a component that needs its own gets its own
  ! operator instance.
  !===================================================================!

  pure subroutine gather_component(flat, num_components, c, comp)

    real(dp), intent(in)  :: flat(:)
    integer , intent(in)  :: num_components, c
    real(dp), intent(out) :: comp(:)

    integer :: i

    do i = 1, size(comp)
       comp(i) = flat((i - 1) * num_components + c)
    end do

  end subroutine gather_component

  pure subroutine scatter_component(comp, num_components, c, flat)

    real(dp), intent(in)    :: comp(:)
    integer , intent(in)    :: num_components, c
    real(dp), intent(inout) :: flat(:)

    integer :: i

    do i = 1, size(comp)
       flat((i - 1) * num_components + c) = comp(i)
    end do

  end subroutine scatter_component

  !===================================================================!
  ! Apply: fetch the input, compile the map once, apply it per
  ! component. A vertex field enters the chain at its innermost
  ! step; an edge field on the vertex landing enters at the
  ! incidence step; no field, or a field on the wrong domain,
  ! returns zeros rather than reading memory it was never given.
  !===================================================================!

  subroutine operator_apply(this, input_graph, input_data, output)

    class(differential_operator), intent(in)       :: this
    class(directed_graph), intent(in)                       :: input_graph
    class(field), intent(in), optional       :: input_data(:)
    class(field), allocatable, intent(inout) :: output

    type(stencil)             :: a
    type(stored_field)        :: out, component
    class(field), allocatable :: swept
    real(dp), allocatable :: q(:), y(:), qc(:), yc(:)
    integer :: nv, ne, nout, nc, c
    logical :: enters_on_edges

    nv = input_graph % num_vertices()
    ne = input_graph % num_edges()

    ! the input: vertex values first; on the vertex landing an edge
    ! field is also lawful and enters at the incidence step
    enters_on_edges = .false.
    call fetch_values(input_data, input_graph, .false., nv, q, nc)
    if (nc == 0 .and. this % landing == SIDE_VERTEX) then
       call fetch_values(input_data, input_graph, .true., ne, q, nc)
       enters_on_edges = nc > 0
    end if

    if (this % landing == SIDE_EDGE) then
       nout = ne
       out  = stored_field(this % name(), input_graph % edge_set(), ne, &
            & num_components=max(nc, 1))
    else
       nout = nv
       out  = stored_field(this % name(), input_graph % vertex_set(), nv, &
            & num_components=max(nc, 1))
    end if

    allocate(y(nout * max(nc, 1)))
    y = 0.0_dp

    if (nc >= 1) then

       a = compiled_map(this, input_graph, enters_on_edges)

       allocate(qc(a % num_columns), yc(a % num_rows))

       do c = 1, nc
          call gather_component(q, nc, c, qc)
          component = stored_field('component', a % column_domain, a % num_columns)
          call component % set_real_vector(qc)
          call a % apply(input_graph, [component], swept)
          call swept % real_vector(yc)
          call scatter_component(yc, nc, c, y)
       end do

    end if

    call out % set_real_vector(y)

    ! a supplied buffer is overwritten, never added to
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine operator_apply

  !===================================================================!
  ! Fetch the input values once and report how many components are carried
  ! in each entry. The field must cover the named side's whole set,
  ! by identity, because the sweep indexes it densely; anything
  ! else leaves a zero-length array and zero components.
  !===================================================================!

  subroutine fetch_values(input_data, input_graph, on_edges, n, q, num_components)

    class(field), intent(in), optional :: input_data(:)
    class(directed_graph)     , intent(in)           :: input_graph
    logical          , intent(in)           :: on_edges
    integer          , intent(in)           :: n
    real(dp), allocatable, intent(out)      :: q(:)
    integer          , intent(out)          :: num_components

    type(graph) :: dom, expected

    num_components = 0

    if (present(input_data)) then
       select type (state => input_data(1))
       class is (stored_field)
          dom = state % domain()
          if (on_edges) then
             expected = input_graph % edge_set()
          else
             expected = input_graph % vertex_set()
          end if
          if (dom % same_as(expected)) then
             num_components = max(state % num_components(), 1)
             call state % real_vector(q)
             if (size(q) == n * num_components) return
             num_components = 0
          end if
       end select
    end if

    allocate(q(0))

  end subroutine fetch_values

end module operation_differential
