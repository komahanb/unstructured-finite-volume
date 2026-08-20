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
! computed here by sparse triple composition with duplicate (row,
! column) entries combined. The coefficient rides the innermost
! step, so a per-edge coefficient makes order 2 the operator
! div(k grad q).
!
! THE ADJOINT IS THE TRANSPOSE. With `adjoint` true the operator
! applies the transpose of the composed matrix - rows and columns
! swapped, the constant dropped, because the adjoint acts on the
! linear part. No reversed step kernels exist: (C B A)^T =
! A^T B^T C^T is an identity of the composition, not code.
!
! THE STENCIL DOOR. The composed map on the vertex landing is
! square, and stencil_of returns it as a stencil_operator - the
! same triples, the same constant - so a minimizer can attach the
! compiled matrix directly. The edge landing is a rectangular
! relation (edges x vertices) and is refused there, because a
! stencil's input and output share one vertex set. apply walks the
! composed triples the same way the stencil walks its edges.
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

module class_graph_differential_operator

  use iso_fortran_env    , only : dp => REAL64
  use operation_action, only : graph_operation
  use graph_directed_view, only : directed_graph
  use graph_field_calculus, only : graph_field
  use graph_fractal      , only : set_graph => graph
  use graph_directed_view     , only : GRAPH_SIDE_VERTEX, GRAPH_SIDE_EDGE
  use class_graph_field  , only : field
  use class_graph_stencil, only : stencil_operator, combine_triples

  implicit none

  private
  public :: differential_operator
  public :: edge_differential_operator
  public :: vertex_differential_operator
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

  !===================================================================!
  ! One affine sparse map, y = A q + k: the triples of A and the
  ! constant k, with the two extents. Private: maps are how this
  ! module computes, not what it promises.
  !===================================================================!

  type :: affine_map

     integer :: nrows = 0
     integer :: ncols = 0

     integer , allocatable :: rows(:)
     integer , allocatable :: cols(:)
     real(dp), allocatable :: weights(:)
     real(dp), allocatable :: constants(:)

  end type affine_map

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
  ! The same constructor, for the vertex landing, plus the adjoint
  ! flag: raised, the operator applies its transpose.
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
  ! THE THREE STEP MAPS, as affine triples.
  !
  ! The average step S, edges x vertices. When it carries the
  ! coefficient and the one-sided choice is on, the sign of the
  ! per-edge coefficient picks the end: positive samples the tail,
  ! negative the head; otherwise both ends, half each. A missing
  ! head reads the boundary value into the constant.
  !===================================================================!

  pure function average_map(g, one_sided_by, with_c, op_c, op_cs, &
       & op_b, op_bs) result(a)

    class(directed_graph), intent(in) :: g
    real(dp)    , intent(in)          :: one_sided_by
    logical     , intent(in)          :: with_c
    real(dp)    , intent(in)          :: op_c
    real(dp), allocatable, intent(in) :: op_cs(:)
    real(dp)    , intent(in)          :: op_b
    real(dp), allocatable, intent(in) :: op_bs(:)
    type(affine_map)                  :: a

    real(dp) :: c, pick, b
    integer  :: e, n

    a % nrows = g % num_edges()
    a % ncols = g % num_vertices()

    allocate(a % rows(2 * a % nrows), a % cols(2 * a % nrows))
    allocate(a % weights(2 * a % nrows))
    allocate(a % constants(a % nrows))
    a % constants = 0.0_dp

    n = 0
    do e = 1, a % nrows

       c = 1.0_dp
       if (with_c) c = coefficient_at(op_c, op_cs, e)

       pick = one_sided_by
       if (with_c) pick = sign(1.0_dp, c) * merge(1.0_dp, 0.0_dp, one_sided_by /= 0.0_dp)

       b = coefficient_at(op_b, op_bs, e)

       if (pick > 0.0_dp) then
          ! the end the walk leaves
          call put(a, n, e, g % edge_tail(e), c)
       else if (pick < 0.0_dp) then
          ! the end the walk enters, or the boundary value
          if (g % edge_has_head(e)) then
             call put(a, n, e, g % edge_head(e), c)
          else
             a % constants(e) = c * b
          end if
       else
          ! both ends, evenly
          call put(a, n, e, g % edge_tail(e), c * 0.5_dp)
          if (g % edge_has_head(e)) then
             call put(a, n, e, g % edge_head(e), c * 0.5_dp)
          else
             a % constants(e) = c * 0.5_dp * b
          end if
       end if

    end do

    call shrink(a, n)

  end function average_map

  !===================================================================!
  ! The difference step G, edges x vertices:
  !
  !    z_e = c (q_head - q_tail) / h_e,
  !
  ! the head replaced by the boundary value - into the constant -
  ! on an edge with no head.
  !===================================================================!

  pure function difference_map(g, with_c, op_c, op_cs, op_h, op_hs, &
       & op_b, op_bs) result(a)

    class(directed_graph), intent(in) :: g
    logical     , intent(in)          :: with_c
    real(dp)    , intent(in)          :: op_c
    real(dp), allocatable, intent(in) :: op_cs(:)
    real(dp)    , intent(in)          :: op_h
    real(dp), allocatable, intent(in) :: op_hs(:)
    real(dp)    , intent(in)          :: op_b
    real(dp), allocatable, intent(in) :: op_bs(:)
    type(affine_map)                  :: a

    real(dp) :: c, w
    integer  :: e, n

    a % nrows = g % num_edges()
    a % ncols = g % num_vertices()

    allocate(a % rows(2 * a % nrows), a % cols(2 * a % nrows))
    allocate(a % weights(2 * a % nrows))
    allocate(a % constants(a % nrows))
    a % constants = 0.0_dp

    n = 0
    do e = 1, a % nrows

       c = 1.0_dp
       if (with_c) c = coefficient_at(op_c, op_cs, e)
       w = c / coefficient_at(op_h, op_hs, e)

       call put(a, n, e, g % edge_tail(e), -w)

       if (g % edge_has_head(e)) then
          call put(a, n, e, g % edge_head(e), w)
       else
          a % constants(e) = w * coefficient_at(op_b, op_bs, e)
       end if

    end do

    call shrink(a, n)

  end function difference_map

  !===================================================================!
  ! The incidence step D, vertices x edges:
  !
  !    y_v = (out - in) / m_v.
  !
  ! Out minus in gives the derivatives their textbook signs; an
  ! edge with no head contributes to its tail alone.
  !===================================================================!

  pure function incidence_map(g, op_m, op_ms) result(a)

    class(directed_graph), intent(in) :: g
    real(dp)    , intent(in)          :: op_m
    real(dp), allocatable, intent(in) :: op_ms(:)
    type(affine_map)                  :: a

    integer :: e, t, h, n

    a % nrows = g % num_vertices()
    a % ncols = g % num_edges()

    allocate(a % rows(2 * a % ncols), a % cols(2 * a % ncols))
    allocate(a % weights(2 * a % ncols))
    allocate(a % constants(a % nrows))
    a % constants = 0.0_dp

    n = 0
    do e = 1, a % ncols

       t = g % edge_tail(e)
       call put(a, n, t, e, 1.0_dp / coefficient_at(op_m, op_ms, t))

       if (g % edge_has_head(e)) then
          h = g % edge_head(e)
          call put(a, n, h, e, -1.0_dp / coefficient_at(op_m, op_ms, h))
       end if

    end do

    call shrink(a, n)

  end function incidence_map

  !===================================================================!
  ! The diagonal map, n x n: y_i = c_i q_i. Order 0 on the vertex
  ! landing, and the per-edge coefficient of an edge-field entry.
  !===================================================================!

  pure function diagonal_map(n, op_c, op_cs) result(a)

    integer     , intent(in)          :: n
    real(dp)    , intent(in)          :: op_c
    real(dp), allocatable, intent(in) :: op_cs(:)
    type(affine_map)                  :: a

    integer :: i

    a % nrows = n
    a % ncols = n

    allocate(a % rows(n), a % cols(n), a % weights(n), a % constants(n))
    a % constants = 0.0_dp

    do i = 1, n
       a % rows(i)    = i
       a % cols(i)    = i
       a % weights(i) = coefficient_at(op_c, op_cs, i)
    end do

  end function diagonal_map

  !===================================================================!
  ! Triple bookkeeping: append one entry; trim to the count.
  !===================================================================!

  pure subroutine put(a, n, r, c, w)

    type(affine_map), intent(inout) :: a
    integer         , intent(inout) :: n
    integer         , intent(in)    :: r, c
    real(dp)        , intent(in)    :: w

    n = n + 1
    a % rows(n)    = r
    a % cols(n)    = c
    a % weights(n) = w

  end subroutine put

  pure subroutine shrink(a, n)

    type(affine_map), intent(inout) :: a
    integer         , intent(in)    :: n

    a % rows    = a % rows(1:n)
    a % cols    = a % cols(1:n)
    a % weights = a % weights(1:n)

  end subroutine shrink

  !===================================================================!
  ! Composition of affine maps:
  !
  !    (A2, k2) o (A1, k1) = (A2 A1, A2 k1 + k2),
  !
  ! by sparse triple product - the inner extent of A2 must equal
  ! A1's row count, checked because a mismatch means the chain was
  ! assembled wrong. Duplicate (row, column) entries are combined,
  ! so the result is a matrix, one entry per pair.
  !===================================================================!

  pure function compose(a2, a1) result(a)

    type(affine_map), intent(in) :: a2, a1
    type(affine_map)             :: a

    integer , allocatable :: ptr(:), r(:), c(:)
    real(dp), allocatable :: w(:)
    integer :: k2, j, n, row1

    if (a2 % ncols /= a1 % nrows) then
       error stop 'differential_operator: composed maps agree on the inner extent'
    end if

    a % nrows = a2 % nrows
    a % ncols = a1 % ncols

    ! group A1's entries by row, counting-sort style
    allocate(ptr(a1 % nrows + 1))
    ptr = 0
    do j = 1, size(a1 % rows)
       ptr(a1 % rows(j) + 1) = ptr(a1 % rows(j) + 1) + 1
    end do
    ptr(1) = 1
    do j = 1, a1 % nrows
       ptr(j + 1) = ptr(j + 1) + ptr(j)
    end do

    block
      integer , allocatable :: order(:), cursor(:)
      allocate(order(size(a1 % rows)), cursor(a1 % nrows))
      cursor = ptr(1:a1 % nrows)
      do j = 1, size(a1 % rows)
         order(cursor(a1 % rows(j))) = j
         cursor(a1 % rows(j)) = cursor(a1 % rows(j)) + 1
      end do

      ! emit one product entry per (A2 entry, matching A1 entry)
      n = 0
      do k2 = 1, size(a2 % rows)
         row1 = a2 % cols(k2)
         n = n + ptr(row1 + 1) - ptr(row1)
      end do

      allocate(r(n), c(n), w(n))
      n = 0
      do k2 = 1, size(a2 % rows)
         row1 = a2 % cols(k2)
         do j = ptr(row1), ptr(row1 + 1) - 1
            n = n + 1
            r(n) = a2 % rows(k2)
            c(n) = a1 % cols(order(j))
            w(n) = a2 % weights(k2) * a1 % weights(order(j))
         end do
      end do
    end block

    call combine_triples(a % nrows, a % ncols, r, c, w, &
         & a % rows, a % cols, a % weights)

    ! the constant travels through the outer map
    allocate(a % constants(a % nrows))
    a % constants = a2 % constants
    do k2 = 1, size(a2 % rows)
       a % constants(a2 % rows(k2)) = a % constants(a2 % rows(k2)) &
            & + a2 % weights(k2) * a1 % constants(a2 % cols(k2))
    end do

  end function compose

  !===================================================================!
  ! The transpose: rows and columns swapped, the constant dropped -
  ! the adjoint acts on the linear part.
  !===================================================================!

  pure function transpose_of(a) result(t)

    type(affine_map), intent(in) :: a
    type(affine_map)             :: t

    t % nrows = a % ncols
    t % ncols = a % nrows

    t % rows    = a % cols
    t % cols    = a % rows
    t % weights = a % weights

    allocate(t % constants(t % nrows))
    t % constants = 0.0_dp

  end function transpose_of

  !===================================================================!
  ! The vertex chain of one order, composed:
  !
  !    order 0:   C                     (the diagonal)
  !    order 2k:  (D G)^k               coefficient on the first G
  !    order 2k+1: (D G)^k D S          coefficient on S, one-sided
  !                                     by its sign
  !
  ! The innermost step carries the coefficient; every later step
  ! runs bare.
  !===================================================================!

  pure function vertex_chain_map(order, g, c, cs, h, hs, m, ms, b, bs) result(a)

    integer     , intent(in)          :: order
    class(directed_graph), intent(in) :: g
    real(dp)    , intent(in)          :: c, h, m, b
    real(dp), allocatable, intent(in) :: cs(:), hs(:), ms(:), bs(:)
    type(affine_map)                  :: a

    type(affine_map) :: d
    integer          :: k
    logical          :: innermost

    if (order <= 0) then
       a = diagonal_map(g % num_vertices(), c, cs)
       return
    end if

    d = incidence_map(g, m, ms)

    k = order
    innermost = .true.

    do while (k > 0)

       if (mod(k, 2) == 1) then
          ! the odd step: the one-sided average, side by the sign
          ! of the coefficient it carries
          if (innermost) then
             a = compose(d, average_map(g, 1.0_dp, .true., c, cs, b, bs))
          else
             a = compose(d, compose(average_map(g, 1.0_dp, .true., &
                  & c, cs, b, bs), a))
          end if
          k = k - 1
       else
          if (innermost) then
             a = compose(d, difference_map(g, .true., c, cs, h, hs, b, bs))
          else
             a = compose(d, compose(difference_map(g, .false., c, cs, &
                  & h, hs, b, bs), a))
          end if
          k = k - 2
       end if

       innermost = .false.

    end do

  end function vertex_chain_map

  !===================================================================!
  ! The whole operator as one affine map, by landing and entry:
  !
  !    edge landing            order 0: S; order 1: G; order n:
  !                            G (bare) after the vertex chain of
  !                            n - 1
  !    vertex landing          the vertex chain, transposed when
  !                            adjoint
  !    vertex landing, entry   diag(c) then D, then the bare vertex
  !    on an edge field        chain of n - 1: order 1 is the
  !                            divergence of the given field
  !===================================================================!

  pure function compiled_map(this, g, enters_on_edges) result(a)

    class(differential_operator), intent(in) :: this
    class(directed_graph)       , intent(in) :: g
    logical                     , intent(in) :: enters_on_edges
    type(affine_map)                         :: a

    real(dp), allocatable :: spent(:)   ! never allocated: the
                                        ! coefficient is applied once

    if (this % landing == GRAPH_SIDE_EDGE) then

       if (this % order <= 0) then
          a = average_map(g, merge(1.0_dp, 0.0_dp, this % one_sided), &
               & .true., this % coefficient, this % coefficients, &
               & this % boundary_value, this % boundary_values)
       else if (this % order == 1) then
          a = difference_map(g, .true., this % coefficient, &
               & this % coefficients, this % spacing, this % spacings, &
               & this % boundary_value, this % boundary_values)
       else
          a = compose(difference_map(g, .false., this % coefficient, &
               & this % coefficients, this % spacing, this % spacings, &
               & this % boundary_value, this % boundary_values), &
               & vertex_chain_map(this % order - 1, g, &
               & this % coefficient, this % coefficients, &
               & this % spacing, this % spacings, &
               & this % measure, this % measures, &
               & this % boundary_value, this % boundary_values))
       end if

    else if (enters_on_edges) then

       a = compose(incidence_map(g, this % measure, this % measures), &
            & diagonal_map(g % num_edges(), this % coefficient, &
            & this % coefficients))
       if (this % order > 1) then
          a = compose(vertex_chain_map(this % order - 1, g, &
               & 1.0_dp, spent, &
               & this % spacing, this % spacings, &
               & this % measure, this % measures, &
               & this % boundary_value, this % boundary_values), a)
       end if

    else

       a = vertex_chain_map(this % order, g, &
            & this % coefficient, this % coefficients, &
            & this % spacing, this % spacings, &
            & this % measure, this % measures, &
            & this % boundary_value, this % boundary_values)
       if (this % adjoint) a = transpose_of(a)

    end if

  end function compiled_map

  !===================================================================!
  ! THE STENCIL DOOR: the compiled operator as a stencil_operator.
  ! Only the vertex landing compiles to one, because a stencil's
  ! input and output share one vertex set; the edge landing is a
  ! rectangular relation and stops the program here.
  !===================================================================!

  impure function stencil_of(operator, input_graph) result(compiled)

    type(differential_operator), intent(in) :: operator
    class(directed_graph)      , intent(in) :: input_graph
    type(stencil_operator)                  :: compiled

    type(affine_map) :: a

    if (operator % landing /= GRAPH_SIDE_VERTEX) then
       error stop 'differential_operator: a stencil is square - only the &
            &vertex landing compiles to one'
    end if

    a = compiled_map(operator, input_graph, enters_on_edges=.false.)

    compiled = stencil_operator(a % rows, a % cols, a % weights, &
         & a % constants, label=operator % name())

  end function stencil_of

  !===================================================================!
  ! The affine sweep, the same walk the stencil's apply performs:
  ! y = k, then every triple carries its weight times the column's
  ! value onto its row.
  !===================================================================!

  pure subroutine sweep(a, q, y)

    type(affine_map), intent(in)  :: a
    real(dp)        , intent(in)  :: q(:)
    real(dp)        , intent(out) :: y(:)

    integer :: j

    y = a % constants
    do j = 1, size(a % rows)
       y(a % rows(j)) = y(a % rows(j)) + a % weights(j) * q(a % cols(j))
    end do

  end subroutine sweep

  !===================================================================!
  ! COMPONENTS. A field may carry several values per entry,
  ! interleaved entry-fastest:
  !
  !      flat((entry - 1) * ncomp + component)
  !
  ! The map is compiled once; each component is gathered, swept,
  ! and scattered back. The parameters are shared by all
  ! components; a component that needs its own gets its own
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
  ! Apply: fetch the input, compile the map once, sweep it per
  ! component. A vertex field enters the chain at its innermost
  ! step; an edge field on the vertex landing enters at the
  ! incidence step; no field, or a field on the wrong domain,
  ! returns zeros rather than reading memory it was never given.
  !===================================================================!

  subroutine operator_apply(this, input_graph, input_data, output)

    class(differential_operator), intent(in)       :: this
    class(directed_graph), intent(in)                       :: input_graph
    class(graph_field), intent(in), optional       :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output

    type(affine_map) :: a
    type(field)      :: out
    real(dp), allocatable :: q(:), y(:), qc(:), yc(:)
    integer :: nv, ne, nout, nc, c
    logical :: enters_on_edges

    nv = input_graph % num_vertices()
    ne = input_graph % num_edges()

    ! the input: vertex values first; on the vertex landing an edge
    ! field is also lawful and enters at the incidence step
    enters_on_edges = .false.
    call fetch_values(input_data, input_graph, .false., nv, q, nc)
    if (nc == 0 .and. this % landing == GRAPH_SIDE_VERTEX) then
       call fetch_values(input_data, input_graph, .true., ne, q, nc)
       enters_on_edges = nc > 0
    end if

    if (this % landing == GRAPH_SIDE_EDGE) then
       nout = ne
       out  = field(this % name(), input_graph % edge_set(), ne, &
            & ncomp=max(nc, 1))
    else
       nout = nv
       out  = field(this % name(), input_graph % vertex_set(), nv, &
            & ncomp=max(nc, 1))
    end if

    allocate(y(nout * max(nc, 1)))
    y = 0.0_dp

    if (nc >= 1) then

       a = compiled_map(this, input_graph, enters_on_edges)

       allocate(qc(a % ncols), yc(a % nrows))

       do c = 1, nc
          call gather_component(q, nc, c, qc)
          call sweep(a, qc, yc)
          call scatter_component(yc, nc, c, y)
       end do

    end if

    call out % set_real_vector(y)

    ! a supplied buffer is overwritten, never added to
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine operator_apply

  !===================================================================!
  ! Fetch the input values once and report how many components ride
  ! in each entry. The field must cover the named side's whole set,
  ! by identity, because the sweep indexes it densely; anything
  ! else leaves a zero-length array and zero components.
  !===================================================================!

  subroutine fetch_values(input_data, input_graph, on_edges, n, q, ncomp)

    class(graph_field), intent(in), optional :: input_data(:)
    class(directed_graph)     , intent(in)           :: input_graph
    logical          , intent(in)           :: on_edges
    integer          , intent(in)           :: n
    real(dp), allocatable, intent(out)      :: q(:)
    integer          , intent(out)          :: ncomp

    type(set_graph) :: dom, expected

    ncomp = 0

    if (present(input_data)) then
       select type (state => input_data(1))
       class is (field)
          dom = state % domain()
          if (on_edges) then
             expected = input_graph % edge_set()
          else
             expected = input_graph % vertex_set()
          end if
          if (dom % same_as(expected)) then
             ncomp = max(state % num_components(), 1)
             call state % get_real_vector(q)
             if (size(q) == n * ncomp) return
             ncomp = 0
          end if
       end select
    end if

    allocate(q(0))

  end subroutine fetch_values

end module class_graph_differential_operator
