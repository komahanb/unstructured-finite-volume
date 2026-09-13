!=====================================================================!
! The stencil operator: a matrix as the tower represents it.
!
! A sparse matrix is a graph with a number on every edge. This type
! is that definition, stored: it CONTAINS a stored graph - one
! directed edge per dependency, column to row - a field of weights
! on that graph's edges, and a field of constants on its vertices,
! the affine part contributed by boundary values. Its apply
! traverses the edges once:
!
!      y(head) += weight * q(tail)
!
! Because the pattern IS a graph, the structure queries require no
! additional operations: the sparsity gives adjacency, the colouring
! traversal runs on it for probing and sweeps, and a coarsener
! applied to it is the Galerkin construction of a coarse operator.
! This is the spatial concretion of the discretization, and
! the family contract is satisfied: the pattern is exposed by contract.
! Nothing here records where the weights came from: a scheme
! computes them from geometry, a multigrid from a product, a test
! from explicit values.
!
! INTERPRETED AND COMPILED. The calculus's differential operator and
! this type are the same mathematics in two execution styles. The
! differential operator INTERPRETS: it reads the host's incidence at
! every apply, matrix-free, always recomputed - the default. The
! stencil is the COMPILED form: weights computed once and traversed
! many times - coarse levels, preconditioners, assembled exactness.
! Neither depends on the other's implementation.
!
! A stencil is also compiled from any operation by evaluation on the
! standard basis (zero state -> constant, basis vector minus constant
! -> column), and transposed from another stencil (edges reversed,
! constants removed: the affine part of a map has no transpose).
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_stencil

  use util_precision  , only : dp
  use view_directed, only : directed_graph
  use field_calculus, only : field, FIELD_REAL
  use operation_action, only : operation, variation, contract
  use operation_action, only : binding, bound_real_vector
  use operation_action, only : emit, dense_of_triples
  use relation_binary, only : group_by_key
  use field_stored  , only : stored_field, typed_field_domain
  use view_directed_stored        , only : stored_directed_graph
  use graph_fractal      , only : graph

  implicit none

  private
  public :: stencil
  public :: combine_triples
  public :: triple_list
  public :: compile_matrix_from_action

  type, extends(operation) :: stencil

     type(stored_directed_graph) :: pattern

     type(stored_field) :: weights
     type(stored_field) :: constants

   contains

     procedure :: apply        => stencil_apply
     procedure :: entries       => stencil_entries
     procedure :: transpose     => stencil_transpose
     procedure :: reverse       => stencil_reverse
     procedure :: restricted    => stencil_restricted
     procedure :: partial_action => stencil_partial_action

  end type stencil

  !===================================================================!
  ! A list of (row, column, weight) triples being gathered, the
  ! capacity doubling when it is exhausted. Every assembly on the tower
  ! appends to one of these and reads the filled entries out; none
  ! maintains its own counter.
  !===================================================================!

  type :: triple_list

     integer :: filled = 0
     integer , allocatable :: rows(:), columns(:)
     real(dp), allocatable :: weights(:)

   contains

     procedure :: assign
     procedure :: entries

  end type triple_list

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

    character(len=:), allocatable :: stencil_label
    integer :: nv

    nv = size(constant)

    this % pattern = stored_directed_graph(nv, tails=columns, heads=rows)

    this % weights   = stored_field('stencil weights', this % pattern % edge_set(), this % pattern % num_edges())
    call this % weights % set_real_vector(weights)
    this % constants = stored_field('stencil constants', this % pattern % vertex_set(), this % pattern % num_vertices())
    call this % constants % set_real_vector(constant)

    stencil_label = 'stencil'
    if (present(label)) stencil_label = label

    ! one argument: the state the matrix multiplies; linear, so its
    ! one exact partial action is of degree one
    call this % declare_arguments(1, [contract(FIELD_REAL, 1)], label=stencil_label, max_degree=1)

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
    character(len=250) :: message

    n = size(a, 1)
    if (size(a, 2) /= n) then
       write(message,'(a,i0,a,i0)') 'stencil: a dense matrix must be square; size(a,1) = ', n, &
            & ', size(a,2) = ', size(a, 2)
       error stop trim(message)
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
  ! the number of values a state stores; it must be a positive
  ! whole multiple of the operation's domain size, and every apply
  ! must return exactly width values; both are checked and stop the
  ! program, because a mismatched column cannot be placed in the
  ! matrix. The label defaults to the operation's name.
  !===================================================================!

  type(stencil) function create_compiled(action, context, width, label) &
       & result(this)

    class(operation)     , intent(in) :: action
    class(directed_graph), intent(in) :: context
    integer              , intent(in) :: width
    character(len=*), intent(in), optional :: label

    type(graph) :: dom
    real(dp), allocatable :: a(:,:), constant(:)
    integer :: n_dom, num_components
    character(len=250) :: message

    call action % domain(context, dom, n_dom)

    if (n_dom <= 0) then
       write(message,'(a,i0)') 'stencil: the operation''s domain must be nonempty; n_dom = ', n_dom
       error stop trim(message)
    end if
    if (width <= 0 .or. mod(width, n_dom) /= 0) then
       write(message,'(a,i0,a,i0)') 'stencil: the width must be a positive whole multiple of &
            &n_dom; width = ', width, ', n_dom = ', n_dom
       error stop trim(message)
    end if

    num_components = width / n_dom

    call compile_matrix_from_action(action, context, dom, n_dom, width, &
         & num_components, a, constant)

    if (present(label)) then
       this = create_dense(a, label)
    else
       this = create_dense(a, action % name())
    end if
    call this % constants % set_real_vector(constant)

  end function create_compiled

  !===================================================================!
  ! The columns of a linear operation, one per basis vector: column j
  ! is F(e_j, h) - F(0, h), the affine part F(0, h) read at j = 0 and
  ! retained as the constant, h the inputs fixed if any. The compiled
  ! stencil and the direct solver both read their matrix here.
  !===================================================================!

  subroutine compile_matrix_from_action(action, context, dom, n_dom, width, &
       & num_components, a, constant, stored)

    class(operation)     , intent(in) :: action
    class(directed_graph), intent(in) :: context
    type(graph)          , intent(in) :: dom
    integer              , intent(in) :: n_dom, width, num_components

    real(dp), allocatable, intent(out) :: a(:,:), constant(:)
    type(stored_field)   , intent(in), optional :: stored(:)

    type(stored_field)        :: state
    class(field), allocatable :: output
    real(dp), allocatable :: e(:), y(:), weights(:)
    integer , allocatable :: rows(:), columns(:)
    integer :: j
    character(len=250) :: message

    ! A STENCIL IS ITS OWN MATRIX. Evaluating it column by column costs
    ! one application per column to recover numbers the pattern
    ! already stores, so a stencil returns its matrix from its edges
    ! and only an operation of another kind is evaluated.
    select type (action)
    type is (stencil)
       if (action % pattern % num_vertices() == width) then
          call action % constants % real_vector(constant)
          call action % entries(rows, columns, weights)
          if (size(constant) == width .and. size(weights) == size(rows)) then
             call dense_of_triples(width, rows, columns, weights, a)
             return
          end if
       end if
    end select

    allocate(a(width, width), e(width))

    do j = 0, width
       e = 0.0_dp
       if (j > 0) e(j) = 1.0_dp
       state = stored_field('basis', dom, n_dom, num_components=num_components)
       call state % set_real_vector(e)
       if (present(stored)) then
          call action % apply(context, action % bind([state, stored]), output)
       else
          call action % apply(context, action % bind([state]), output)
       end if
       call output % real_vector(y)
       if (size(y) /= width) then
          write(message,'(a,i0,a,i0)') 'stencil: the operation result must match the width; &
               &size(y) = ', size(y), ', width = ', width
          error stop trim(message)
       end if
       if (j == 0) then
          constant = y
       else
          a(:, j) = y - constant
       end if
    end do

  end subroutine compile_matrix_from_action

  !===================================================================!
  ! THE STENCIL RESTRICTED to a subset of its vertices, everything
  ! outside the subset fixed at the values given: the rows retained
  ! are those of the subset, a column inside it remains a dependency,
  ! and a column outside it is added to the row's constant as its
  ! weight times the fixed value. The result is the same linear map,
  ! restricted to the subset, with the outside as an affine part.
  ! A member outside the vertices, or values of the wrong extent,
  ! stops the program.
  !===================================================================!

  function stencil_restricted(this, retained, values) result(sub)

    class(stencil), intent(in) :: this
    integer       , intent(in) :: retained(:)
    real(dp)      , intent(in) :: values(:)
    type(stencil) :: sub

    integer , allocatable :: sub_of(:), rows(:), columns(:), heads(:), tails(:)
    real(dp), allocatable :: weights(:), w(:), constant(:), stored(:)
    type(triple_list) :: triples
    integer :: n, m, e, row
    character(len=250) :: message

    n = this % pattern % num_vertices()
    m = size(retained)

    if (size(values) /= n) then
       write(message,'(a,i0,a,i0)') 'stencil: one fixed value is required per vertex; &
            &size(values) = ', size(values), ', num_vertices = ', n
       error stop trim(message)
    end if
    if (any(retained < 1) .or. any(retained > n)) then
       write(message,'(a,i0,a,i0,a,i0)') 'stencil: every retained member must be one of the &
            &vertices 1..', n, '; retained ranges from ', minval(retained), ' to ', maxval(retained)
       error stop trim(message)
    end if

    allocate(sub_of(n), source=0)
    do e = 1, m
       sub_of(retained(e)) = e
    end do

    call this % entries(heads, tails, w)
    call this % constants % real_vector(stored)

    allocate(constant(m))
    constant = stored(retained)

    do e = 1, size(heads)
       row = sub_of(heads(e))
       if (row == 0) cycle
       if (sub_of(tails(e)) > 0) then
          call triples % assign(row, sub_of(tails(e)), w(e))
       else
          constant(row) = constant(row) + w(e) * values(tails(e))
       end if
    end do

    call triples % entries(rows, columns, weights)
    sub = stencil(rows, columns, weights, constant, label=this % name() // ' restricted')

  end function stencil_restricted

  !===================================================================!
  ! The (row, column, weight) triples of the matrix, one per edge in
  ! the pattern's edge order: the row is the head, the column the
  ! tail. Weights not yet stored give an empty weight list.
  !===================================================================!

  pure subroutine stencil_entries(this, rows, columns, weights)

    class(stencil)       , intent(in)  :: this
    integer , allocatable, intent(out) :: rows(:), columns(:)
    real(dp), allocatable, intent(out) :: weights(:)

    integer :: e

    rows    = [(this % pattern % edge_head(e), e = 1, this % pattern % num_edges())]
    columns = [(this % pattern % edge_tail(e), e = 1, this % pattern % num_edges())]
    call this % weights % real_vector(weights)

  end subroutine stencil_entries

  !===================================================================!
  ! y = constants + the dependency edges, traversed once: each edge
  ! adds its weight times the tail's value to its head.
  !===================================================================!

  subroutine stencil_apply(this, input_graph, inputs, output)

    class(stencil), intent(in)            :: this
    class(directed_graph), intent(in)                       :: input_graph
    type(binding), intent(in), optional       :: inputs(:)
    class(field), allocatable, intent(inout) :: output

    type(stored_field)   :: out
    type(typed_field_domain)   :: image
    real(dp), allocatable :: q(:), y(:)

    call this % constants % real_vector(y)

    if (present(inputs)) then
       call bound_real_vector(inputs, this % argument(1), q)
       call accumulate_edges(this, q, y)
    end if

    image = typed_field_domain(input_graph % vertex_set(), input_graph % num_vertices())
    out   = image % real_field(this % name(), y)

    call emit(out, output)

  end subroutine stencil_apply

  !===================================================================!
  ! The edge traversal, shared by apply and the tangent: each edge
  ! adds its weight times the tail's value to the head.
  !===================================================================!

  subroutine accumulate_edges(this, q, y)

    class(stencil), intent(in)    :: this
    real(dp)      , intent(in)    :: q(:)
    real(dp)      , intent(inout) :: y(:)

    real(dp), pointer :: w(:)

    ! the weights are read in place: an earlier apply copied the
    ! whole edge vector before reading it
    w => this % weights % real_values()
    if (.not. associated(w)) return

    ! The graph selects its orientation once and supplies its grouped
    ! incidence to one read. Each row is accumulated locally and written
    ! once, with no copy of the topology and no per-edge procedure call.
    call this % pattern % read_incoming(accumulate_rows)

  contains

    subroutine accumulate_rows(offsets, indices, sources)

      integer, intent(in) :: offsets(:), indices(:), sources(:)
      real(dp) :: acc
      integer :: e, p, v

      do v = 1, size(offsets) - 1
         acc = 0.0_dp
         do p = offsets(v), offsets(v + 1) - 1
            e = indices(p)
            acc = acc + w(e) * q(sources(e))
         end do
         y(v) = y(v) + acc
      end do

    end subroutine accumulate_rows

  end subroutine accumulate_edges

  !===================================================================!
  ! A stencil is linear, so it is its own tangent: the first partial
  ! action in its one argument is the edge traversal on the direction,
  ! without the constants. An order past one or an argument other
  ! than one stops the program, because a linear map in one argument
  ! has no other partial.
  !===================================================================!

  subroutine stencil_partial_action(this, input_graph, inputs, &
       & variations, output)

    class(stencil), intent(in)               :: this
    class(directed_graph), intent(in)        :: input_graph
    type(binding), intent(in)                 :: inputs(:)
    type(variation), intent(in)              :: variations(:)
    class(field), allocatable, intent(inout) :: output

    type(stored_field)   :: out
    type(typed_field_domain)   :: image
    real(dp), allocatable :: v(:), y(:)
    character(len=250) :: message

    associate (u1 => inputs); end associate

    call this % require_owned(variations)

    if (size(variations) /= 1) then
       write(message,'(a,i0)') 'stencil: a stencil''s max_degree is 1, so exactly one variation &
            &is required; size(variations) = ', size(variations)
       error stop trim(message)
    end if
    if (.not. variations(1) % argument_is(this % argument(1))) then
       error stop 'stencil: variations(1) does not vary the stencil''s one argument'
    end if

    call variations(1) % direction(v)
    allocate(y(this % pattern % num_vertices()))
    y = 0.0_dp
    call accumulate_edges(this, v, y)

    image = typed_field_domain(input_graph % vertex_set(), input_graph % num_vertices())
    out   = image % real_field(this % name(), y)

    call emit(out, output)

  end subroutine stencil_partial_action

  !===================================================================!
  ! The transpose: the same pattern read in the reverse direction, so
  ! the weight that multiplied the tail's value into the head now
  ! multiplies the head's into the tail - no edge is rebuilt and no
  ! weight is moved, and the transpose of the transpose is this stencil
  ! exactly. The constants are removed, because the affine part of a
  ! map has no transpose.
  !===================================================================!

  type(stencil) function stencil_transpose(this) result(transposed)

    class(stencil), intent(in) :: this

    real(dp), allocatable :: zeros(:)

    transposed % pattern   = this % pattern % transpose()
    transposed % weights   = this % weights

    allocate(zeros(this % pattern % num_vertices()))
    zeros = 0.0_dp
    transposed % constants = stored_field('stencil constants', &
         & transposed % pattern % vertex_set(), transposed % pattern % num_vertices())
    call transposed % constants % set_real_vector(zeros)

    ! the transpose is an operation of one argument as the original is;
    ! attached as a matrix-vector product, it is called with that argument
    call transposed % declare_arguments(1, [contract(FIELD_REAL, 1)], &
         & label='transpose of ' // this % name(), max_degree=1)

  end function stencil_transpose

  !===================================================================!
  ! The transpose taken in place: the pattern's orientation reversed
  ! at constant cost, the constants removed, no weight and no edge
  ! copied. The result equals the owning transpose of the former value.
  !===================================================================!

  subroutine stencil_reverse(this)

    class(stencil), intent(inout) :: this

    real(dp), allocatable :: zeros(:)

    call this % pattern % reverse()

    allocate(zeros(this % pattern % num_vertices()))
    zeros = 0.0_dp
    this % constants = stored_field('stencil constants', &
         & this % pattern % vertex_set(), this % pattern % num_vertices())
    call this % constants % set_real_vector(zeros)

    call this % declare_arguments(1, [contract(FIELD_REAL, 1)], &
         & label='transpose of ' // this % name(), max_degree=1)

  end subroutine stencil_reverse

  !===================================================================!
  ! Combine duplicate (row, column) entries of a weighted triple
  ! list: a matrix has one entry per pair, so equal pairs sum. Two
  ! stable groupings (by column, then by row) place equal pairs
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

  !===================================================================!
  ! One triple appended, the capacity doubling when it is exhausted.
  !===================================================================!

  pure subroutine assign(this, row, column, weight)

    class(triple_list), intent(inout) :: this
    integer           , intent(in)    :: row, column
    real(dp)          , intent(in)    :: weight

    integer , allocatable :: wider(:)
    real(dp), allocatable :: heavier(:)
    integer :: capacity

    if (.not. allocated(this % rows)) then
       allocate(this % rows(16), this % columns(16), this % weights(16))
       this % filled = 0
    end if

    if (this % filled == size(this % rows)) then
       capacity = 2 * size(this % rows)
       allocate(wider(capacity))
       wider(1:this % filled) = this % rows
       call move_alloc(wider, this % rows)
       allocate(wider(capacity))
       wider(1:this % filled) = this % columns
       call move_alloc(wider, this % columns)
       allocate(heavier(capacity))
       heavier(1:this % filled) = this % weights
       call move_alloc(heavier, this % weights)
    end if

    this % filled = this % filled + 1
    this % rows(this % filled)    = row
    this % columns(this % filled) = column
    this % weights(this % filled) = weight

  end subroutine assign

  !===================================================================!
  ! The triples placed so far, exactly, in the order placed.
  !===================================================================!

  pure subroutine entries(this, rows, columns, weights)

    class(triple_list)   , intent(in)  :: this
    integer , allocatable, intent(out) :: rows(:), columns(:)
    real(dp), allocatable, intent(out) :: weights(:)

    if (allocated(this % rows)) then
       rows    = this % rows(1:this % filled)
       columns = this % columns(1:this % filled)
       weights = this % weights(1:this % filled)
    else
       allocate(rows(0), columns(0), weights(0))
    end if

  end subroutine entries

end module operation_stencil
