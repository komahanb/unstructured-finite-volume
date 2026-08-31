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

  use util_precision  , only : dp
  use view_directed, only : directed_graph
  use field_calculus, only : field, FIELD_REAL
  use operation_action, only : operation, variation, contract
  use operation_binding, only : binding, bound_inputs, bound_real_vector
  use operation_action, only : emit
  use operation_discretization     , only : discretization
  use relation_binary, only : group_by_key
  use field_stored  , only : stored_field
  use view_directed_stored        , only : stored_directed_graph
  use graph_fractal      , only : graph

  implicit none

  private
  public :: stencil
  public :: combine_triples
  public :: triple_list
  public :: compile_matrix_from_action

  type, extends(discretization) :: stencil

     type(stored_directed_graph) :: pattern

     type(stored_field) :: weights
     type(stored_field) :: constants

     character(len=:), allocatable :: label

   contains

     procedure :: name         => stencil_name
     procedure :: apply        => stencil_apply
     procedure :: transpose     => stencil_transpose
     procedure :: restricted    => stencil_restricted
     procedure :: max_degree     => stencil_max_degree
     procedure :: partial_action => stencil_partial_action

  end type stencil

  !===================================================================!
  ! A list of (row, column, weight) triples being gathered, the room
  ! doubling when it runs out. Every assembly on the tower appends to
  ! one of these and reads the filled entries out; none keeps its own
  ! counter.
  !===================================================================!

  type :: triple_list

     integer :: filled = 0
     integer , allocatable :: rows(:), columns(:)
     real(dp), allocatable :: weights(:)

   contains

     procedure :: place
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

    ! one argument: the state the matrix multiplies
    call this % declare_arguments(1, [contract(FIELD_REAL, 1)])

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

    type(graph) :: dom
    real(dp), allocatable :: a(:,:), constant(:)
    integer :: n_dom, num_components

    call action % domain(on, dom, n_dom)

    if (n_dom <= 0) then
       error stop 'stencil: the operation''s domain is nonempty'
    end if
    if (width <= 0 .or. mod(width, n_dom) /= 0) then
       error stop 'stencil: the width carries a whole number per member'
    end if

    num_components = width / n_dom

    call compile_matrix_from_action(action, on, dom, n_dom, width, &
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
  ! kept as the constant, h the inputs held fixed if any. The compiled
  ! stencil and the direct solver both read their matrix here.
  !===================================================================!

  subroutine compile_matrix_from_action(action, on, dom, n_dom, width, &
       & num_components, a, constant, held)

    class(operation)     , intent(in) :: action
    class(directed_graph), intent(in) :: on
    type(graph)          , intent(in) :: dom
    integer              , intent(in) :: n_dom, width, num_components

    real(dp), allocatable, intent(out) :: a(:,:), constant(:)
    type(stored_field)   , intent(in), optional :: held(:)

    type(stored_field)        :: state
    class(field), allocatable :: output
    real(dp), allocatable :: e(:), y(:)
    real(dp), pointer     :: wgt(:)
    integer :: j, k, ne, row, column

    ! A STENCIL IS ITS OWN MATRIX. Probing it column by column costs
    ! one application per column to recover numbers the pattern
    ! already carries, so a stencil answers from its edges and only
    ! an operation of another kind is probed.
    select type (action)
    type is (stencil)
       if (action % pattern % num_vertices() == width) then
          call action % constants % real_vector(constant)
          wgt => action % weights % real_values()
          if (size(constant) == width .and. associated(wgt)) then
             allocate(a(width, width))
             a  = 0.0_dp
             ne = action % pattern % num_edges()
             do k = 1, ne
                row    = action % pattern % edge_head(k)
                column = action % pattern % edge_tail(k)
                a(row, column) = a(row, column) + wgt(k)
             end do
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
       if (present(held)) then
          call action % apply(on, [state, held], output)
       else
          call action % apply(on, [state], output)
       end if
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

  end subroutine compile_matrix_from_action

  !===================================================================!
  ! THE STENCIL RESTRICTED to a subset of its vertices, everything
  ! outside the subset held at the values given: the rows kept are
  ! those of the subset, a column inside it stays a dependency, and a
  ! column outside it is taken into the row's constant as its weight
  ! times the held value. What comes back is the same linear map,
  ! seen from inside the subset, with the outside as an affine part.
  ! A member outside the vertices, or values of the wrong extent,
  ! stops the program.
  !===================================================================!

  function stencil_restricted(this, kept, values) result(sub)

    class(stencil), intent(in) :: this
    integer       , intent(in) :: kept(:)
    real(dp)      , intent(in) :: values(:)
    type(stencil) :: sub

    integer , allocatable :: sub_of(:), rows(:), columns(:)
    real(dp), allocatable :: weights(:), w(:), constant(:), held(:)
    type(triple_list) :: triples
    integer :: n, m, e, row, column

    n = this % pattern % num_vertices()
    m = size(kept)

    if (size(values) /= n) then
       error stop 'stencil: one held value per vertex'
    end if
    if (any(kept < 1) .or. any(kept > n)) then
       error stop 'stencil: a kept member is one of the vertices'
    end if

    allocate(sub_of(n), source=0)
    do e = 1, m
       sub_of(kept(e)) = e
    end do

    call this % weights   % real_vector(w)
    call this % constants % real_vector(held)

    allocate(constant(m))
    constant = held(kept)

    do e = 1, this % pattern % num_edges()
       row    = sub_of(this % pattern % edge_head(e))
       if (row == 0) cycle
       column = this % pattern % edge_tail(e)
       if (sub_of(column) > 0) then
          call triples % place(row, sub_of(column), w(e))
       else
          constant(row) = constant(row) + w(e) * values(column)
       end if
    end do

    call triples % entries(rows, columns, weights)
    sub = stencil(rows, columns, weights, constant, label=this % label // ' restricted')

  end function stencil_restricted

  pure function stencil_name(this) result(name)

    class(stencil), intent(in) :: this
    character(len=:), allocatable :: name

    name = this % label

  end function stencil_name
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
    type(binding), allocatable :: bound(:)
    real(dp), allocatable :: q(:), y(:)

    call this % constants % real_vector(y)

    if (present(input_data)) then
       call bound_inputs(this, input_data, bound)
       call bound_real_vector(bound, this % argument(1), q)
       call accumulate_edges(this, q, y)
    end if

    out = stored_field(this % label, input_graph % vertex_set(), input_graph % num_vertices())
    call out % set_real_vector(y)

    call emit(out, output)

  end subroutine stencil_apply

  !===================================================================!
  ! The one edge walk, shared by apply and the tangent: each edge
  ! carries its weight times the tail's value onto the head.
  !===================================================================!

  subroutine accumulate_edges(this, q, y)

    class(stencil), intent(in)    :: this
    real(dp)      , intent(in)    :: q(:)
    real(dp)      , intent(inout) :: y(:)

    real(dp), pointer :: w(:)
    real(dp) :: acc
    integer :: e, p, v, nv

    ! the weights are read where they are held: one apply copied the
    ! whole edge vector before reading it
    w => this % weights % real_values()
    if (.not. associated(w)) return

    ! a row at a time, through the lists the pattern already groups by
    ! endpoint: the row's sum lives in a register and each row is
    ! written once, where an edge at a time wrote to a scattered
    ! subscript. which list holds the rows is the reversal's question,
    ! and it is asked once
    associate (g => this % pattern)
      nv = g % num_vertices()
      if (g % reversed) then
         do v = 1, nv
            acc = 0.0_dp
            do p = g % xout(v), g % xout(v + 1) - 1
               e   = g % eout(p)
               acc = acc + w(e) * q(g % head(e))
            end do
            y(v) = y(v) + acc
         end do
      else
         do v = 1, nv
            acc = 0.0_dp
            do p = g % xin(v), g % xin(v + 1) - 1
               e   = g % ein(p)
               acc = acc + w(e) * q(g % tail(e))
            end do
            y(v) = y(v) + acc
         end do
      end if
    end associate

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

  subroutine stencil_partial_action(this, input_graph, input_data, &
       & variations, output)

    class(stencil), intent(in)               :: this
    class(directed_graph), intent(in)        :: input_graph
    class(field), intent(in)                 :: input_data(:)
    type(variation), intent(in)              :: variations(:)
    class(field), allocatable, intent(inout) :: output

    type(stored_field)   :: out
    real(dp), allocatable :: v(:), y(:)

    associate (u1 => input_data); end associate

    call this % require_owned(variations)

    if (size(variations) /= 1) then
       error stop 'stencil: the requested order is within max_degree'
    end if
    if (.not. variations(1) % argument_is(this % argument(1))) then
       error stop 'stencil: the partial action is taken in the one argument'
    end if

    call variations(1) % direction(v)
    allocate(y(this % pattern % num_vertices()))
    y = 0.0_dp
    call accumulate_edges(this, v, y)

    out = stored_field(this % label, input_graph % vertex_set(), input_graph % num_vertices())
    call out % set_real_vector(y)

    call emit(out, output)

  end subroutine stencil_partial_action

  !===================================================================!
  ! The transpose: the same pattern read the other way, so the weight
  ! that carried the tail's value onto the head now carries the head's
  ! onto the tail - no edge is rebuilt and no weight is moved, and the
  ! transpose of the transpose is this stencil exactly. The constants
  ! are dropped, because the affine part of a map has no transpose.
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

    transposed % label = 'transpose of ' // this % label

    ! the transpose is an operation of one argument as the original is;
    ! attached as a matvec, it is asked for that argument
    call transposed % declare_arguments(1, [contract(FIELD_REAL, 1)])

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

  !===================================================================!
  ! One triple appended, the room doubling when it runs out.
  !===================================================================!

  pure subroutine place(this, row, column, weight)

    class(triple_list), intent(inout) :: this
    integer           , intent(in)    :: row, column
    real(dp)          , intent(in)    :: weight

    integer , allocatable :: wider(:)
    real(dp), allocatable :: heavier(:)
    integer :: room

    if (.not. allocated(this % rows)) then
       allocate(this % rows(16), this % columns(16), this % weights(16))
       this % filled = 0
    end if

    if (this % filled == size(this % rows)) then
       room = 2 * size(this % rows)
       allocate(wider(room))
       wider(1:this % filled) = this % rows
       call move_alloc(wider, this % rows)
       allocate(wider(room))
       wider(1:this % filled) = this % columns
       call move_alloc(wider, this % columns)
       allocate(heavier(room))
       heavier(1:this % filled) = this % weights
       call move_alloc(heavier, this % weights)
    end if

    this % filled = this % filled + 1
    this % rows(this % filled)    = row
    this % columns(this % filled) = column
    this % weights(this % filled) = weight

  end subroutine place

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
