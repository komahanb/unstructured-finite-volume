!=====================================================================!
! The time-marching family: an edge function on a scheme's coupling
! that also declares how the scheme marches.
!
! One type stores every family as data. A multistep family (Adams,
! BDF) is an order p and the functional its coefficients are read
! from; a Runge-Kutta family (DIRK) is a tableau; Newmark is the
! coefficient pair beta and gamma. The constructors below fix the
! data, and every query is one procedure
! that reads it. The coefficients produced are dimensionless: the
! weight an edge finally stores is the coefficient times a power of
! the step, and that power is fixed by the two degrees the edge
! joins, which operation_weight multiplies in.
!
!=====================================================================!
!
!                    THE TWO LAGRANGE FUNCTIONALS
!
! A multistep family interpolates through the instants behind k at
! the nodes
!
!      u_j  =  -(t_k - t_(k-j)) / dt_k ,     u_0 = 0 ,
!
! and every coefficient is either the slope at zero of a Lagrange
! basis function through those nodes (a difference) or its integral
! over the last step (a quadrature). On a uniform grid every u_j is
! -j and the tabulated coefficients result; they are the uniform
! value of the formula, not a separate case.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_family

  use util_precision  , only : dp
  use operation_edge_function , only : edge_function
  use view_directed_connectivity, only : connectivity_graph
  use util_derivative_terms   , only : derivative_terms, value, &
       & operator(+), operator(-), operator(*), operator(/)

  implicit none

  private
  public :: family
  public :: adams_family, bdf_family, dirk_family
  public :: newmark_family
  public :: implicit_midpoint, crouzeix_two_stage, crouzeix_three_stage, &
       & hairer_wanner_five_stage
  public :: slope_at_zero, integral_over_step
  public :: adams, bdf, dirk, chain

  ! the families by order: adams(p) and bdf(p) the multistep families
  ! of order p, dirk(p) the tabulated tableau of order p
  interface adams
     module procedure adams_family
  end interface adams

  interface bdf
     module procedure bdf_family
  end interface bdf

  interface dirk
     module procedure dirk_of_order
  end interface dirk

  !===================================================================!
  ! A CHAIN OF FAMILIES over the instants: scheme(k) approximates the
  ! derivative from instant from(k) to the instant before from(k+1),
  ! the last to the end. from(1) is one and the list rises.
  !===================================================================!

  type :: chain

     type(family), allocatable :: scheme(:)
     integer     , allocatable :: from(:)

   contains

     procedure :: num_blocks

  end type chain

  interface chain
     module procedure create_chain
  end interface chain

  ! The coupling geometries a family reads its edges by.

  integer, parameter :: FAMILY_ADAMS   = 1
  integer, parameter :: FAMILY_BDF     = 2
  integer, parameter :: FAMILY_DIRK    = 3
  integer, parameter :: FAMILY_NEWMARK = 4

  !===================================================================!
  ! A family is its label, its geometry and its data: the order of a
  ! multistep family, the tableau (a, b) of a DIRK family, or beta and
  ! gamma of a Newmark family. Every unstaged family has one stage
  ! with weight one, so b = [1]. The order is read by the multistep
  ! geometries only; a one-step family stores one, its history depth.
  !===================================================================!

  type, extends(edge_function) :: family

     integer , private :: geometry = FAMILY_BDF
     integer , private :: order    = 1
     real(dp), private, allocatable :: a(:,:)
     real(dp), private, allocatable :: b(:)
     real(dp), private :: beta  = 0.0_dp
     real(dp), private :: gamma = 0.0_dp

   contains

     procedure :: history_depth    => family_history_depth
     procedure :: num_stages       => family_num_stages
     procedure :: stage_weight     => family_stage_weight
     procedure :: step_quadrature  => family_step_quadrature
     procedure :: primary_degree   => family_primary_degree
     procedure :: row_pattern      => family_row_pattern
     procedure :: block_connectivity => family_block_connectivity
     procedure :: stage_connectivity => family_stage_connectivity
     procedure :: edge_coefficient => family_edge_coefficient

  end type family

contains

  !===================================================================!
  ! The constructors. An order below one, a non-square tableau, or a
  ! tableau entry above the diagonal stops the program.
  !===================================================================!

  function create(label, geometry, order, a, b, beta, gamma) result(this)

    character(len=*), intent(in) :: label
    integer         , intent(in) :: geometry, order
    real(dp)        , intent(in) :: a(:,:), b(:)
    real(dp), optional, intent(in) :: beta, gamma
    type(family) :: this

    this % geometry = geometry
    this % order    = order
    this % a        = a
    this % b        = b
    if (present(beta))  this % beta  = beta
    if (present(gamma)) this % gamma = gamma
    call this % declare_edge_arguments(label)

  end function create

  function adams_family(order) result(this)

    integer, intent(in) :: order
    type(family) :: this

    character(len=250) :: message

    if (order < 1) then
       write(message,'(a,i0)') 'operation_family: the order must be positive; order = ', order
       error stop trim(message)
    end if
    this = create('adams-moulton', FAMILY_ADAMS, order, &
         & reshape([1.0_dp], [1, 1]), [1.0_dp])

  end function adams_family

  function bdf_family(order) result(this)

    integer, intent(in) :: order
    type(family) :: this

    character(len=250) :: message

    if (order < 1) then
       write(message,'(a,i0)') 'operation_family: the order must be positive; order = ', order
       error stop trim(message)
    end if
    this = create('bdf', FAMILY_BDF, order, reshape([1.0_dp], [1, 1]), [1.0_dp])

  end function bdf_family

  function dirk_family(a, b) result(this)

    real(dp), intent(in) :: a(:,:)
    real(dp), intent(in) :: b(:)
    type(family) :: this

    integer :: i, j
    character(len=250) :: message

    if (size(a, 1) /= size(a, 2) .or. size(b) /= size(a, 1)) then
       write(message,'(a,i0,a,i0,a,i0)') 'operation_family: the tableau must be square with one &
            &weight per stage; shape(a) = [', size(a, 1), ',', size(a, 2), '], size(b) = ', size(b)
       error stop trim(message)
    end if
    do i = 1, size(a, 1)
       do j = i + 1, size(a, 2)
          if (a(i, j) /= 0.0_dp) then
             write(message,'(a,i0,a,i0,a,es12.5)') 'operation_family: a diagonally implicit &
                  &tableau must have no entry above its diagonal; a(', i, ',', j, ') = ', a(i, j)
             error stop trim(message)
          end if
       end do
    end do
    this = create('dirk', FAMILY_DIRK, 1, a, b)

  end function dirk_family

  !===================================================================!
  ! Newmark is the pair (beta, gamma) and nothing else: its rows
  !
  !    q_k  = q_(k-1) + h q'_(k-1) + h^2 (1/2 - beta) q''_(k-1)
  !                                + h^2 beta q''_k
  !    q'_k = q'_(k-1) + h (1 - gamma) q''_(k-1) + h gamma q''_k
  !
  ! reproduce t^0, t^1, t^2 for every pair and are second order in h
  ! when gamma = 1/2, first order otherwise; no integer states more.
  ! The pair (0, 0) is the explicit Taylor step of the jet.
  !===================================================================!

  function newmark_family(beta, gamma) result(this)

    real(dp), intent(in) :: beta, gamma
    type(family) :: this

    this = create('newmark', FAMILY_NEWMARK, 1, reshape([1.0_dp], [1, 1]), &
         & [1.0_dp], beta=beta, gamma=gamma)

  end function newmark_family

  !===================================================================!
  ! The queries, each one read from the data.
  !===================================================================!

  pure subroutine require_newmark_equation(this, equation_degree)

    class(family), intent(in) :: this
    integer      , intent(in) :: equation_degree

    if (this % geometry == FAMILY_NEWMARK .and. equation_degree /= 2) then
       error stop 'operation_family: Newmark requires a second-order (degree-two) state; a &
            &different equation degree was passed'
    end if

  end subroutine require_newmark_equation

  pure integer function family_history_depth(this, equation_degree)

    class(family), intent(in) :: this
    integer      , intent(in) :: equation_degree

    call require_newmark_equation(this, equation_degree)
    select case (this % geometry)
    case (FAMILY_ADAMS)
       family_history_depth = max(this % order - 1, 1)
    case (FAMILY_BDF)
       family_history_depth = this % order
    case default
       family_history_depth = 1
    end select

  end function family_history_depth

  pure integer function family_primary_degree(this, equation_degree)

    class(family), intent(in) :: this
    integer      , intent(in) :: equation_degree

    call require_newmark_equation(this, equation_degree)
    family_primary_degree = merge(0, equation_degree, this % geometry == FAMILY_BDF)

  end function family_primary_degree

  pure integer function family_num_stages(this)

    class(family), intent(in) :: this

    family_num_stages = size(this % b)

  end function family_num_stages

  pure real(dp) function family_stage_weight(this, i)

    class(family), intent(in) :: this
    integer      , intent(in) :: i

    if (i < 1 .or. i > size(this % b)) then
       error stop 'operation_family: this stage is not one of the tableau'
    end if
    family_stage_weight = this % b(i)

  end function family_stage_weight

  !===================================================================!
  ! The interpolatory quadrature over the step [t_(k-1), t_k], on the
  ! instants the family's rows read, so that it is of the order of the
  ! scheme: min(order, k) for a multistep family; for Newmark the two
  ! instants k - 1 and k its rows read (offsets 1 and 0), which is the
  ! trapezoidal rule, weights (1/2, 1/2), exact on linear integrands,
  ! global error -(h^2/12) [f'(T) - f'(0)] + O(h^4), the second order
  ! of the pair with gamma = 1/2. A rule on a third instant would read
  ! history the family does not hold. At k = 1 the single node has
  ! the measure dt_1 = 0. A staged family integrates over its stages
  ! with the tableau weights b, not over instants, and stops the
  ! program here. An instant outside the block stops the program.
  !
  ! The complete rule keeps the full node count at every step after
  ! the first: at a step k with fewer than n instants behind it the
  ! nodes are the instants 1 .. n, those ahead of the step included,
  ! and the weights are the integrals over [t_(k-1), t_k] of the
  ! Lagrange basis through them; the rule is then of order n at every
  ! step. It is the quadrature of an enriched functional (the
  ! functional-error estimator), never of the family's own F_h. The
  ! instant each weight multiplies is returned in instant, ordered k,
  ! k - 1, .., 1, then k + 1, .., n.
  !===================================================================!

  pure subroutine family_step_quadrature(this, dt, k, weight, complete, instant)

    class(family)         , intent(in) :: this
    type(derivative_terms), intent(in) :: dt(:)
    integer               , intent(in) :: k
    type(derivative_terms), allocatable, intent(out) :: weight(:)
    logical               , intent(in) , optional :: complete
    integer  , allocatable, intent(out), optional :: instant(:)

    type(derivative_terms), allocatable :: u(:)
    integer :: num_nodes, n, j
    logical :: ahead

    if (k < 1 .or. k > size(dt)) then
       error stop 'operation_family: this quadrature was evaluated at an instant outside the block'
    end if
    select case (this % geometry)
    case (FAMILY_ADAMS, FAMILY_BDF)
       n = this % order
    case (FAMILY_NEWMARK)
       n = 2
    case default
       error stop 'operation_family: this staged family cannot be integrated over instants; it &
            &integrates over its stages by the tableau weights'
    end select
    ahead = .false.
    if (present(complete)) ahead = complete .and. k > 1 .and. k < n
    num_nodes = min(n, k)
    if (ahead) num_nodes = n
    if (num_nodes > size(dt)) then
       error stop 'operation_family: this complete rule reads instants beyond the block'
    end if
    if (ahead) then
       u = nodes_ahead(dt, k, n)
    else
       u = nodes(dt, k, num_nodes)
    end if
    allocate(weight(num_nodes))
    do j = 1, num_nodes
       weight(j) = integral_over_step(u, j - 1)
    end do
    if (present(instant)) instant = [(k - j + 1, j = 1, min(n, k)), (j, j = k + 1, num_nodes)]

  end subroutine family_step_quadrature

  !===================================================================!
  ! The row pattern: the offsets an edge reads and the degree at each,
  ! for the row that head_degree a degree. Adams reads the value at
  ! offset one and the derivative at the p previous instants; BDF
  ! reads the degree below at offsets 0..p; Newmark reads q and q'
  ! from the previous instant and q'' from the previous and current
  ! instants; DIRK has no derived rows.
  !===================================================================!

  pure subroutine family_row_pattern(this, head_degree, equation_degree, &
       & offset, tail_degree)

    class(family), intent(in) :: this
    integer      , intent(in) :: head_degree, equation_degree
    integer, allocatable, intent(out) :: offset(:), tail_degree(:)

    integer :: i

    call require_newmark_equation(this, equation_degree)
    select case (this % geometry)
    case (FAMILY_ADAMS)
       if (head_degree < 0 .or. head_degree >= equation_degree) then
          allocate(offset(0), tail_degree(0))
          return
       end if
       offset        = [1, (i, i = 0, this % order - 1)]
       tail_degree = [head_degree, (head_degree + 1, i = 0, this % order - 1)]
    case (FAMILY_BDF)
       if (head_degree < 1 .or. head_degree > equation_degree) then
          allocate(offset(0), tail_degree(0))
          return
       end if
       offset = [(i, i = 0, this % order)]
       allocate(tail_degree(this % order + 1), source=head_degree - 1)
    case (FAMILY_NEWMARK)
       select case (head_degree)
       case (0)
          offset      = [1, 1, 1, 0]
          tail_degree = [0, 1, 2, 2]
       case (1)
          offset      = [1, 1, 0]
          tail_degree = [1, 2, 2]
       case default
          allocate(offset(0), tail_degree(0))
       end select
    case default
       allocate(offset(0), tail_degree(0))
    end select

  end subroutine family_row_pattern

  !===================================================================!
  ! The connectivity of the family over a block of n instants at nd
  ! degrees: the row pattern placed at every instant it fits behind,
  ! instant outer, degree inner, the primary degree omitted. The edge
  ! order is the order the coupling's relation is formed in.
  !===================================================================!

  function family_block_connectivity(this, nd, n) result(connectivity)

    class(family), intent(in) :: this
    integer      , intent(in) :: nd, n
    type(connectivity_graph) :: connectivity

    integer, allocatable :: tails(:), heads(:), tail_degree(:), head_degree(:)
    integer, allocatable :: offset(:), degrees_of(:)
    integer :: primary, kk, d, e, counted, pass

    primary = this % primary_degree(nd - 1)
    do pass = 1, 2
       counted = 0
       do kk = 1, n
          do d = 0, nd - 1
             if (d == primary) cycle
             call this % row_pattern(d, nd - 1, offset, degrees_of)
             if (size(offset) == 0) cycle
             if (kk - maxval(offset) < 1) cycle
             do e = 1, size(offset)
                counted = counted + 1
                if (pass == 2) then
                   tails(counted)         = kk - offset(e)
                   heads(counted)         = kk
                   tail_degree(counted) = degrees_of(e)
                   head_degree(counted)    = d
                end if
             end do
          end do
       end do
       if (pass == 1) allocate(tails(counted), heads(counted), &
            & tail_degree(counted), head_degree(counted))
    end do

    connectivity = connectivity_graph(n, tails, heads, tail_degree, head_degree)

  end function family_block_connectivity

  !===================================================================!
  ! The connectivity of the tableau over one step at nd degrees, on
  ! the vertices 1 (the instant behind), 2..s+1 (the stages) and s+2
  ! (the instant ahead). Below the top degree, stage i reads its own
  ! degree at the instant behind and the degree above at stages 1..i,
  ! and the instant ahead reads its own degree at the instant behind
  ! and the degree above at every stage:
  !
  !    Q_i = q_(k-1) + h sum_(j<=i) a_ij Q'_j ,  q_k = q_(k-1) + h sum_j b_j Q'_j .
  !
  ! The top degree has no row of the family: at every stage and at
  ! the instant ahead it is the law's own row. Degree outer, stage
  ! inner; within a row the instant behind precedes the stages.
  !===================================================================!

  function family_stage_connectivity(this, nd) result(connectivity)

    class(family), intent(in) :: this
    integer      , intent(in) :: nd
    type(connectivity_graph) :: connectivity

    integer, allocatable :: tails(:), heads(:), tail_degree(:), head_degree(:)
    integer :: s, d, i, j, at

    s  = size(this % b)
    at = (nd - 1) * (s * (s + 1) / 2 + 2 * s + 1)
    allocate(tails(at), heads(at), tail_degree(at), head_degree(at))
    at = 0
    do d = 0, nd - 2
       do i = 1, s
          call behind(1 + i)
          do j = 1, i
             at = at + 1
             tails(at) = 1 + j
             heads(at) = 1 + i
             tail_degree(at) = d + 1
             head_degree(at) = d
          end do
       end do
       call behind(2 + s)
       do j = 1, s
          at = at + 1
          tails(at) = 1 + j
          heads(at) = 2 + s
          tail_degree(at) = d + 1
          head_degree(at) = d
       end do
    end do

    connectivity = connectivity_graph(s + 2, tails, heads, tail_degree, head_degree)

  contains

    subroutine behind(head)
      integer, intent(in) :: head
      at = at + 1
      tails(at) = 1
      heads(at) = head
      tail_degree(at) = d
      head_degree(at) = d
    end subroutine behind

  end function family_stage_connectivity

  !===================================================================!
  ! The dimensionless coefficient on one edge. An edge from an earlier
  ! instant, a source degree the geometry does not read, or an offset
  ! beyond the order stops the program.
  !===================================================================!

  pure function family_edge_coefficient(this, dt, tail, head, &
       & tail_degree, head_degree) result(c)

    class(family)         , intent(in) :: this
    type(derivative_terms), intent(in) :: dt(:)
    integer               , intent(in) :: tail, head, tail_degree, head_degree
    type(derivative_terms) :: c

    integer :: i, j, s

    select case (this % geometry)

    case (FAMILY_ADAMS)
       i = head - tail
       if (i < 0) error stop 'operation_family: an Adams edge runs from an instant later than its head'
       if (head_degree < 0) then
          error stop 'operation_family: an Adams constraint names a negative head_degree'
       end if
       if (tail_degree == head_degree) then
          if (i /= 1) error stop 'operation_family: the same degree must be read at offset one, &
               &but was read at a different offset'
          c = derivative_terms(1.0_dp, dt(head))
       else if (tail_degree == head_degree + 1) then
          if (i >= this % order) error stop 'operation_family: the quadrature must span at most &
               &p instants, but this offset reaches beyond p'
          c = integral_over_step(nodes(dt, head, this % order), i)
       else
          error stop 'operation_family: this source names neither the constraint''s degree nor &
               &one above it'
       end if

    case (FAMILY_BDF)
       j = head - tail
       if (j < 0) error stop 'operation_family: a BDF edge runs from an instant later than its head'
       if (head_degree < 1) error stop 'operation_family: a derived BDF row''s head_degree must &
            &name a derivative, but head_degree is not positive'
       if (tail_degree /= head_degree - 1) then
          error stop 'operation_family: this source does not name the degree below the one determined'
       end if
       if (j > this % order) error stop 'operation_family: this BDF row reaches beyond p instants'
       c = slope_at_zero(nodes(dt, head, this % order + 1), j)

    case (FAMILY_NEWMARK)
       if (head_degree < 0 .or. head_degree > 1) then
          error stop 'operation_family: a Newmark row must determine q or qdot; this head_degree &
               &names neither'
       end if
       if (tail /= head .and. tail /= head - 1) then
          error stop 'operation_family: a Newmark row must read the current or previous instant; &
               &this tail names neither'
       end if
       if (tail == head) then
          if (tail_degree /= 2) then
             error stop 'operation_family: Newmark must read the current acceleration here, but &
                  &tail_degree names a different derivative'
          end if
          if (head_degree == 0) then
             c = derivative_terms(this % beta, dt(head))
          else
             c = derivative_terms(this % gamma, dt(head))
          end if
       else
          select case (head_degree)
          case (0)
             select case (tail_degree)
             case (0, 1)
                c = derivative_terms(1.0_dp, dt(head))
             case (2)
                c = derivative_terms(0.5_dp - this % beta, dt(head))
             case default
                error stop 'operation_family: a Newmark value row must read q, qdot or qddot; &
                     &this tail_degree names none of them'
             end select
          case (1)
             select case (tail_degree)
             case (1)
                c = derivative_terms(1.0_dp, dt(head))
             case (2)
                c = derivative_terms(1.0_dp - this % gamma, dt(head))
             case default
                error stop 'operation_family: a Newmark velocity row must read qdot or qddot; &
                     &this tail_degree names neither'
             end select
          end select
       end if

    case default
       s = size(this % b)
       if (head == 1 .or. tail == 2 + s) then
          error stop 'operation_family: this edge runs from the initial instant or a stage into &
               &a later vertex, which a staged family does not permit'
       end if
       if (tail == 1) then
          if (tail_degree /= head_degree) then
             error stop 'operation_family: the instant behind must be read at the constraint''s &
                  &degree; tail_degree and head_degree differ here'
          end if
          c = derivative_terms(1.0_dp, dt(head))
          return
       end if
       if (tail_degree /= head_degree + 1) then
          error stop 'operation_family: a stage must be read at the degree above the &
               &constraint''s; tail_degree does not name that degree here'
       end if
       j = tail - 1
       if (head == 2 + s) then
          c = derivative_terms(this % b(j), dt(head))
          return
       end if
       i = head - 1
       if (j > i) error stop 'operation_family: a stage must read stages at or before it; this &
            &one reads a later stage'
       c = derivative_terms(this % a(i, j), dt(head))

    end select

  end function family_edge_coefficient

  !===================================================================!
  ! The n Lagrange nodes ending at instant k, scaled by dt(k) and
  ! negated so that the step just taken is the interval [-1, 0]. A
  ! step read that is not positive stops the program.
  !===================================================================!

  pure function nodes(dt, k, n) result(u)

    type(derivative_terms), intent(in) :: dt(:)
    integer               , intent(in) :: k, n
    type(derivative_terms) :: u(0:n-1)

    integer :: j

    do j = k - n + 2, k
       if (value(dt(j)) <= 0.0_dp) then
          error stop 'operation_family: one of the time steps behind instant k is not positive'
       end if
    end do
    u(0) = derivative_terms(0.0_dp, dt(k))
    do j = 1, n - 1
       u(j) = u(j - 1) - dt(k - j + 1) / dt(k)
    end do

  end function nodes

  !===================================================================!
  ! The n Lagrange nodes at the instants 1 .. n of a step k below n,
  ! in the order of the complete rule (k, k - 1, .., 1, k + 1, .., n),
  ! scaled by dt(k) so that the step just taken is [-1, 0]: the
  ! instants behind by the backward recurrence of nodes, those ahead
  ! by the steps following k. A step read that is not positive stops
  ! the program.
  !===================================================================!

  pure function nodes_ahead(dt, k, n) result(u)

    type(derivative_terms), intent(in) :: dt(:)
    integer               , intent(in) :: k, n
    type(derivative_terms) :: u(0:n-1)

    integer :: j

    do j = 2, n
       if (value(dt(j)) <= 0.0_dp) then
          error stop 'operation_family: one of the time steps among the first n is not positive'
       end if
    end do
    u(0:k-1) = nodes(dt, k, k)
    u(k) = derivative_terms(0.0_dp, dt(k)) + dt(k + 1) / dt(k)
    do j = k + 2, n
       u(j - 1) = u(j - 2) + dt(j) / dt(k)
    end do

  end function nodes_ahead

  !===================================================================!
  ! The slope at zero of the j-th Lagrange basis function through the
  ! nodes u.
  !===================================================================!

  pure function slope_at_zero(u, j) result(s)

    type(derivative_terms), intent(in) :: u(0:)
    integer               , intent(in) :: j
    type(derivative_terms) :: s

    type(derivative_terms) :: term
    integer :: m, i, n

    n = ubound(u, 1)
    s = derivative_terms(0.0_dp, u(j))
    do m = 0, n
       if (m == j) cycle
       term = derivative_terms(1.0_dp, u(j)) / (u(j) - u(m))
       do i = 0, n
          if (i == j .or. i == m) cycle
          term = term * (-u(i)) / (u(j) - u(i))
       end do
       s = s + term
    end do

  end function slope_at_zero

  !===================================================================!
  ! The monomial coefficients of the j-th Lagrange basis function
  ! through the nodes u, by repeated multiplication of the linear
  ! factors.
  !===================================================================!

  pure function basis_polynomial(u, j) result(c)

    type(derivative_terms), intent(in) :: u(0:)
    integer               , intent(in) :: j
    type(derivative_terms) :: c(0:ubound(u, 1))

    type(derivative_terms) :: shifted(0:ubound(u, 1))
    integer :: m, l, n

    n = ubound(u, 1)
    do l = 0, n
       c(l) = derivative_terms(0.0_dp, u(j))
    end do
    c(0) = derivative_terms(1.0_dp, u(j))
    do m = 0, n
       if (m == j) cycle
       shifted(0) = derivative_terms(0.0_dp, u(j))
       do l = 1, n
          shifted(l) = c(l - 1)
       end do
       do l = 0, n
          c(l) = (shifted(l) - u(m) * c(l)) / (u(j) - u(m))
       end do
    end do

  end function basis_polynomial

  !===================================================================!
  ! The integral over [-1, 0] of the j-th Lagrange basis function
  ! through the nodes u.
  !===================================================================!

  pure function integral_over_step(u, j) result(w)

    type(derivative_terms), intent(in) :: u(0:)
    integer               , intent(in) :: j
    type(derivative_terms) :: w

    type(derivative_terms) :: c(0:ubound(u, 1))
    integer :: n

    c = basis_polynomial(u, j)
    w = derivative_terms(0.0_dp, u(j))
    do n = 0, ubound(u, 1)
       w = w + (real((-1)**n, dp) / real(n + 1, dp)) * c(n)
    end do

  end function integral_over_step

  !===================================================================!
  ! The tabulated DIRK families.
  !===================================================================!

  function implicit_midpoint() result(this)

    type(family) :: this

    this = dirk_family(reshape([0.5_dp], [1, 1]), [1.0_dp])

  end function implicit_midpoint

  function crouzeix_two_stage() result(this)

    type(family) :: this
    real(dp) :: g

    g = (3.0_dp + sqrt(3.0_dp)) / 6.0_dp
    this = dirk_family(reshape([g, 1.0_dp - 2.0_dp * g, 0.0_dp, g], [2, 2]), &
         & [0.5_dp, 0.5_dp])

  end function crouzeix_two_stage

  function crouzeix_three_stage() result(this)

    type(family) :: this
    real(dp) :: g, w, pi

    pi = acos(-1.0_dp)
    g  = cos(pi / 18.0_dp) / sqrt(3.0_dp) + 0.5_dp
    w  = 1.0_dp / (6.0_dp * (1.0_dp - 2.0_dp * g)**2)
    this = dirk_family(reshape( &
         & [g,             0.5_dp - g, 2.0_dp * g,      &
         &  0.0_dp,        g,          1.0_dp - 4.0_dp * g, &
         &  0.0_dp,        0.0_dp,     g], [3, 3]),    &
         & [w, 1.0_dp - 2.0_dp * w, w])

  end function crouzeix_three_stage

  !===================================================================!
  ! The tabulated DIRK family of an order: the implicit midpoint rule
  ! at two, Crouzeix's two stages at three, Crouzeix's three stages at
  ! four. Another order stops the program.
  !===================================================================!

  function dirk_of_order(order) result(this)

    integer, intent(in) :: order
    type(family) :: this

    character(len=250) :: message

    select case (order)
    case (2)
       this = implicit_midpoint()
    case (3)
       this = crouzeix_two_stage()
    case (4)
       this = crouzeix_three_stage()
    case default
       write(message,'(a,i0)') 'operation_family: a DIRK tableau is tabulated at orders 2, 3 and 4; &
            &order = ', order
       error stop trim(message)
    end select

  end function dirk_of_order

  function create_chain(schemes, from) result(this)

    type(family), intent(in) :: schemes(:)
    integer     , intent(in) :: from(:)
    type(chain) :: this

    integer :: k
    character(len=250) :: message

    if (size(schemes) /= size(from)) then
       write(message,'(a,i0,a,i0)') 'operation_family: a chain requires one first instant per &
            &scheme; size(schemes) = ', size(schemes), ', size(from) = ', size(from)
       error stop trim(message)
    end if
    if (size(from) < 1 .or. from(1) /= 1) then
       error stop 'operation_family: a chain begins at the first instant'
    end if
    do k = 2, size(from)
       if (from(k) <= from(k - 1)) then
          write(message,'(a,i0,a,i0,a,i0)') 'operation_family: the first instants of a chain rise; &
               &from(', k - 1, ') = ', from(k - 1), ', from(k) = ', from(k)
          error stop trim(message)
       end if
    end do
    this % scheme = schemes
    this % from   = from

  end function create_chain

  pure integer function num_blocks(this)
    class(chain), intent(in) :: this
    num_blocks = size(this % scheme)
  end function num_blocks

  function hairer_wanner_five_stage() result(this)

    type(family) :: this
    real(dp) :: a(5, 5)

    a = 0.0_dp
    a(1, 1:1) = [1.0_dp / 4.0_dp]
    a(2, 1:2) = [1.0_dp / 2.0_dp, 1.0_dp / 4.0_dp]
    a(3, 1:3) = [17.0_dp / 50.0_dp, -1.0_dp / 25.0_dp, 1.0_dp / 4.0_dp]
    a(4, 1:4) = [371.0_dp / 1360.0_dp, -137.0_dp / 2720.0_dp, 15.0_dp / 544.0_dp, 1.0_dp / 4.0_dp]
    a(5, 1:5) = [25.0_dp / 24.0_dp, -49.0_dp / 48.0_dp, 125.0_dp / 16.0_dp, -85.0_dp / 12.0_dp, &
         &       1.0_dp / 4.0_dp]
    this = dirk_family(a, a(5, :))

  end function hairer_wanner_five_stage

end module operation_family
