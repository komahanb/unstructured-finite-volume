!=====================================================================!
! The stencil operator: a matrix as the tower says it.
!
! A sparse matrix is a graph with a number on every edge. This type
! is that sentence made literal: it HOLDS a stored graph - one
! directed edge per dependency, column to row - a field of weights
! on that graph's edges, and a field of constants on its vertices,
! the affine part boundary values leave behind. Its apply walks the
! edges once:
!
!      y(head) += weight * q(tail)
!
! Because the pattern IS a graph, the structure questions come free:
! the sparsity answers adjacency, the colouring walk runs on it for
! probing and sweeps, and a coarsener applied to it is the Galerkin
! road to a coarse operator. This is the spatial concretion of the
! discretization operator, and the family contract holds: the
! pattern is exposed by law. Nothing here knows where the weights
! came from: a scheme fills them from geometry, a multigrid from a
! product, a test by hand.
!
! INTERPRETED AND COMPILED. The calculus's differential operator and
! this type are the same mathematics in two execution styles. The
! differential operator INTERPRETS: it reads the host's incidence at
! every apply, matrix-free, always fresh - the right default. The
! stencil is the COMPILED form: weights computed once and walked
! many times - coarse levels, preconditioners, assembled exactness.
! Neither learns the other's business.
!
! A stencil is also compiled from any operation by evaluation on the
! standard basis (zero state -> constant, basis vector minus constant
! -> column), and transposed from another stencil (edges reversed,
! constants dropped: the affine part of a map has no transpose).
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_stencil

  use iso_fortran_env    , only : dp => REAL64
  use view_directed, only : directed_graph
  use field_calculus, only : field
  use operation_action, only : operation
  use operation_discretization     , only : discretization
  use relation_binary, only : group_by_key
  use field_stored  , only : stored_field
  use view_directed_stored        , only : stored_directed_graph
  use graph_fractal      , only : graph

  implicit none

  private
  public :: stencil
  public :: combine_triples

  type, extends(discretization) :: stencil

     type(stored_directed_graph) :: pattern

     type(stored_field) :: weights
     type(stored_field) :: constants

     character(len=:), allocatable :: label

   contains

     procedure :: name         => stencil_name
     procedure :: domain       => stencil_domain
     procedure :: apply        => stencil_apply
     procedure :: dependencies => stencil_dependencies
     procedure :: transpose     => stencil_transpose
     procedure :: max_degree     => stencil_max_degree
     procedure :: partial_action => stencil_partial_action

  end type stencil

  interface stencil
     module procedure create
     module procedure create_dense
     module procedure create_compiled
  end interface stencil

contains

  !===================================================================!
  ! Build from the triples and the constants. The triples become the
  ! dependency graph - one edge per (row, column) pair, tail at the
  ! column, head at the row - and the numbers become its fields.
  !===================================================================!

  type(stencil) function create(rows, columns, weights, &
       & constant, label) result(this)

    integer , intent(in) :: rows(:)
    integer , intent(in) :: columns(:)
    real(dp), intent(in) :: weights(:)
    real(dp), intent(in) :: constant(:)
    character(len=*), intent(in), optional :: label

    integer :: nv

    nv = size(constant)

    this % pattern = stored_directed_graph(nv, tails=columns, heads=rows)

    this % weights   = stored_field('stencil weights', this % pattern % edge_set(), this % pattern % num_edges())
    call this % weights % set_real_vector(weights)
    this % constants = stored_field('stencil constants', this % pattern % vertex_set(), this % pattern % num_vertices())
    call this % constants % set_real_vector(constant)

    if (present(label)) then
       this % label = label
    else
       this % label = 'stencil'
    end if

  end function create

  !===================================================================!
  ! Build from a dense matrix: each entry becomes one weighted edge,
  ! column to row, constants zero. The matrix must be square,
  ! because a stencil's input and output share one vertex set;
  ! a rectangular array stops the program.
  !===================================================================!

  type(stencil) function create_dense(a, label) result(this)

    real(dp)        , intent(in)           :: a(:,:)
    character(len=*), intent(in), optional :: label

    integer , allocatable :: rows(:), columns(:)
    real(dp), allocatable :: weights(:), constant(:)
    integer :: n, i, j, e

    n = size(a, 1)
    if (size(a, 2) /= n) then
       error stop 'stencil: a dense matrix is square'
    end if

    allocate(rows(n * n), columns(n * n), weights(n * n), constant(n))
    constant = 0.0_dp

    e = 0
    do j = 1, n
       do i = 1, n
          e          = e + 1
          rows(e)    = i
          columns(e) = j
          weights(e) = a(i, j)
       end do
    end do

    this = create(rows, columns, weights, constant, label)

  end function create_dense

  !===================================================================!
  ! Build from an operation by evaluation on the standard basis: the
  ! operation applied to the zero state is the constant, and applied
  ! to each basis vector minus that constant is one column. width is
  ! the number of values a state carries; it must be a positive
  ! whole multiple of the operation's domain size, and every apply
  ! must return exactly width values; both are checked and stop the
  ! program, because a mismatched column cannot be placed in the
  ! matrix. The label defaults to the operation's name.
  !===================================================================!

  type(stencil) function create_compiled(action, on, width, label) &
       & result(this)

    class(operation)     , intent(in) :: action
    class(directed_graph), intent(in) :: on
    integer              , intent(in) :: width
    character(len=*), intent(in), optional :: label

    type(stored_field)        :: state
    class(field), allocatable :: output
    type(graph) :: dom
    real(dp), allocatable :: a(:,:), e(:), y(:), constant(:)
    integer :: n_dom, num_components, j

    call action % domain(on, dom, n_dom)

    if (n_dom <= 0) then
       error stop 'stencil: the operation''s domain is nonempty'
    end if
    if (width <= 0 .or. mod(width, n_dom) /= 0) then
       error stop 'stencil: the width carries a whole number per member'
    end if

    num_components = width / n_dom

    allocate(a(width, width), e(width))

    ! j = 0 is the zero state, whose value is the constant
    do j = 0, width
       e = 0.0_dp
       if (j > 0) e(j) = 1.0_dp
       state = stored_field('basis', dom, n_dom, num_components=num_components)
       call state % set_real_vector(e)
       call action % apply(on, [state], output)
       call output % real_vector(y)
       if (size(y) /= width) then
          error stop 'stencil: the operation result matches the width'
       end if
       if (j == 0) then
          constant = y
       else
          a(:, j) = y - constant
       end if
    end do

    if (present(label)) then
       this = create_dense(a, label)
    else
       this = create_dense(a, action % name())
    end if
    call this % constants % set_real_vector(constant)

  end function create_compiled

  pure function stencil_name(this) result(name)

    class(stencil), intent(in) :: this
    character(len=:), allocatable :: name

    name = this % label

  end function stencil_name

  subroutine stencil_domain(this, input_graph, domain, num_entries)

    class(stencil), intent(in)    :: this
    class(directed_graph), intent(in)               :: input_graph
    type(graph), intent(out) :: domain
    integer        , intent(out) :: num_entries

    associate (u1 => this); end associate

    domain   = input_graph % all_vertices()
    num_entries = input_graph % num_vertices()

  end subroutine stencil_domain

  !===================================================================!
  ! y = constants + the dependency edges, walked once: each edge
  ! carries its weight times the tail's value onto its head.
  !===================================================================!

  subroutine stencil_apply(this, input_graph, input_data, output)

    class(stencil), intent(in)            :: this
    class(directed_graph), intent(in)                       :: input_graph
    class(field), intent(in), optional       :: input_data(:)
    class(field), allocatable, intent(inout) :: output

    type(stored_field)   :: out
    real(dp), allocatable :: q(:), y(:)

    call this % constants % real_vector(y)

    if (present(input_data)) then
       call input_data(1) % real_vector(q)
       call accumulate_edges(this, q, y)
    end if

    out = stored_field(this % label, input_graph % vertex_set(), input_graph % num_vertices())
    call out % set_real_vector(y)

    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine stencil_apply

  !===================================================================!
  ! The one edge walk, shared by apply and the tangent: each edge
  ! carries its weight times the tail's value onto the head.
  !===================================================================!

  subroutine accumulate_edges(this, q, y)

    class(stencil), intent(in)    :: this
    real(dp)      , intent(in)    :: q(:)
    real(dp)      , intent(inout) :: y(:)

    real(dp), allocatable :: w(:)
    integer :: e

    call this % weights % real_vector(w)
    do e = 1, this % pattern % num_edges()
       y(this % pattern % edge_head(e)) = &
            & y(this % pattern % edge_head(e)) &
            & + w(e) * q(this % pattern % edge_tail(e))
    end do

  end subroutine accumulate_edges

  !===================================================================!
  ! A stencil is linear, so it is its own tangent: the first partial
  ! action in its one input slot is the edge walk on the direction,
  ! without the constants. An order past one or a slot other than
  ! one stops the program, because a linear map in one slot has no
  ! other partial.
  !===================================================================!

  pure function stencil_max_degree(this) result(degree)

    class(stencil), intent(in) :: this
    integer :: degree

    associate (u1 => this); end associate

    degree = 1

  end function stencil_max_degree

  subroutine stencil_partial_action(this, input_graph, input_data, slots, &
       & directions, output)

    class(stencil), intent(in)               :: this
    class(directed_graph), intent(in)        :: input_graph
    class(field), intent(in)                 :: input_data(:)
    integer, intent(in)                      :: slots(:)
    class(field), intent(in)                 :: directions(:)
    class(field), allocatable, intent(inout) :: output

    type(stored_field)   :: out
    real(dp), allocatable :: v(:), y(:)

    associate (u1 => input_data); end associate

    if (size(slots) /= 1) then
       error stop 'stencil: the requested order is within max_degree'
    end if
    if (slots(1) /= 1) then
       error stop 'stencil: the partial action is taken in the one input slot'
    end if

    call directions(1) % real_vector(v)
    allocate(y(this % pattern % num_vertices()))
    y = 0.0_dp
    call accumulate_edges(this, v, y)

    out = stored_field(this % label, input_graph % vertex_set(), input_graph % num_vertices())
    call out % set_real_vector(y)

    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine stencil_partial_action

  !===================================================================!
  ! The contract's answer: the pattern IS a graph, handed out whole.
  !===================================================================!

  subroutine stencil_dependencies(this, pattern)

    class(stencil), intent(in)    :: this
    class(directed_graph), allocatable, intent(out) :: pattern

    allocate(pattern, source=this % pattern)

  end subroutine stencil_dependencies

  !===================================================================!
  ! The transpose: every edge reversed, so the weight that carried
  ! the tail's value onto the head now carries the head's onto the
  ! tail. The constants are dropped, because the affine part of a
  ! map has no transpose.
  !===================================================================!

  type(stencil) function stencil_transpose(this) result(transposed)

    class(stencil), intent(in) :: this

    integer , allocatable :: tails(:), heads(:)
    real(dp), allocatable :: weights(:), zeros(:)
    integer :: e, ne, nv

    ne = this % pattern % num_edges()
    nv = this % pattern % num_vertices()

    allocate(tails(ne), heads(ne), zeros(nv))
    do e = 1, ne
       tails(e) = this % pattern % edge_tail(e)
       heads(e) = this % pattern % edge_head(e)
    end do
    zeros = 0.0_dp

    call this % weights % real_vector(weights)

    transposed = create(rows=tails, columns=heads, weights=weights, &
         & constant=zeros, label='transpose of ' // this % label)

  end function stencil_transpose

  !===================================================================!
  ! Combine duplicate (row, column) entries of a weighted triple
  ! list: a matrix has one entry per pair, so equal pairs sum. Two
  ! stable groupings (by column, then by row) bring equal pairs
  ! adjacent; one pass merges them. Used wherever triples are
  ! produced with repeats - a sparse product's emitted terms, a
  ! Galerkin coarsening's aggregated edges - before they become a
  ! stencil.
  !===================================================================!

  pure subroutine combine_triples(nrows, ncols, r, c, w, rows, cols, weights)

    integer , intent(in) :: nrows, ncols
    integer , intent(in) :: r(:), c(:)
    real(dp), intent(in) :: w(:)
    integer , allocatable, intent(out) :: rows(:)
    integer , allocatable, intent(out) :: cols(:)
    real(dp), allocatable, intent(out) :: weights(:)

    integer, allocatable :: identity(:), ptr(:), by_c(:), by_rc(:)
    integer :: j, n, m

    n = size(r)

    allocate(identity(n))
    identity = [(j, j = 1, n)]
    call group_by_key(ncols, c, identity, ptr, by_c)

    block
      integer, allocatable :: rkey(:), order(:)
      allocate(rkey(n))
      do j = 1, n
         rkey(j) = r(by_c(j))
      end do
      call group_by_key(nrows, rkey, identity, ptr, order)
      allocate(by_rc(n))
      do j = 1, n
         by_rc(j) = by_c(order(j))
      end do
    end block

    allocate(rows(n), cols(n), weights(n))
    m = 0
    do j = 1, n
       if (m > 0) then
          if (rows(m) == r(by_rc(j)) .and. cols(m) == c(by_rc(j))) then
             weights(m) = weights(m) + w(by_rc(j))
             cycle
          end if
       end if
       m = m + 1
       rows(m)    = r(by_rc(j))
       cols(m)    = c(by_rc(j))
       weights(m) = w(by_rc(j))
    end do

    rows    = rows(1:m)
    cols    = cols(1:m)
    weights = weights(1:m)

  end subroutine combine_triples

end module operation_stencil
