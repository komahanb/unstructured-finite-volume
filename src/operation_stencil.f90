!=====================================================================!
! The stencil operator: an affine sparse map, y = A q + k, as the
! tower says it.
!
! A sparse matrix is a weighted relation between two finite sets:
! one entry per dependent (row, column) pair, its weight, and a
! constant per row - the affine part boundary values leave behind.
! The type holds the triples, the constants, and the two sets as
! declared domains when its maker knows them, so a map between
! different sets - edges x vertices - is one stencil exactly as a
! map within one set is. Built from bare arrays it has no domains
! and lands on the host graph's vertex set.
!
! THE TWO LAWS, WRITTEN ONCE.
!
!    composition   (A2, k2) o (A1, k1) = (A2 A1, A2 k1 + k2)
!    transpose     (A, k)^T = (A^T, 0)
!
! Every map the tower compiles - a parity chain, a scheme's tangent,
! a Galerkin coarse operator - is a composition of stencils, and
! every adjoint is the transpose of one. (C B A)^T = A^T B^T C^T is
! an identity of the composition, so no reversed kernel exists: the
! transpose is taken once, here, after the chain is composed.
!
! THE PATTERN IS DERIVED. The dependency graph - one edge per entry,
! column to row - is what a colouring and a Galerkin coarsening
! read; dependencies builds it when asked and nothing stores it, so
! a stencil composed inside a chain never pays for neighbour lists
! it does not read - the law the directed graph applies to its own
! relations.
!
! The differential operator composes its chain fresh at every apply;
! the stencil is the compiled form, entries computed once and
! applied many times. A stencil is also compiled from any operation
! by evaluation on the standard basis: the zero state gives the
! constant, each basis vector minus the constant gives one column.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_stencil

  use iso_fortran_env    , only : dp => REAL64
  use view_directed, only : directed_graph
  use field_calculus, only : field
  use operation_action, only : operation, variation
  use operation_discretization     , only : discretization
  use relation_binary, only : group_by_key
  use field_stored  , only : stored_field
  use view_directed_stored        , only : stored_directed_graph
  use graph_fractal      , only : graph
  use token_identity, only : token

  implicit none

  private
  public :: stencil
  public :: compose
  public :: combine_triples

  type, extends(discretization) :: stencil

     integer :: num_rows    = 0
     integer :: num_columns = 0

     integer , allocatable :: rows(:)
     integer , allocatable :: columns(:)
     real(dp), allocatable :: weights(:)
     real(dp), allocatable :: constants(:)

     type(graph) :: row_domain
     type(graph) :: column_domain
     integer :: row_entries    = 0
     integer :: column_entries = 0

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
  ! Build from the triples and the constants. The row count is the
  ! length of the constants; the column count defaults to it, so a
  ! square map names one extent. A side's domain comes with its
  ! member count, the rows or columns divided among the members as
  ! components; a domain without an assigned identity is the same
  ! as none, and the result lands on the host. Triple arrays of
  ! unequal length or an entry outside the extents stop the program.
  !===================================================================!

  function create(rows, columns, weights, constant, label, num_columns, &
       & row_domain, column_domain, row_entries, column_entries) result(this)

    integer , intent(in) :: rows(:)
    integer , intent(in) :: columns(:)
    real(dp), intent(in) :: weights(:)
    real(dp), intent(in) :: constant(:)
    character(len=*), intent(in), optional :: label
    integer         , intent(in), optional :: num_columns
    type(graph)     , intent(in), optional :: row_domain
    type(graph)     , intent(in), optional :: column_domain
    integer         , intent(in), optional :: row_entries
    integer         , intent(in), optional :: column_entries
    type(stencil) :: this

    if (size(rows) /= size(columns) .or. size(rows) /= size(weights)) then
       error stop 'stencil: one row, one column and one weight per entry'
    end if

    this % num_rows    = size(constant)
    this % num_columns = this % num_rows
    if (present(num_columns)) this % num_columns = num_columns

    if (size(rows) > 0) then
       if (minval(rows) < 1 .or. maxval(rows) > this % num_rows .or. &
            & minval(columns) < 1 .or. maxval(columns) > this % num_columns) then
          error stop 'stencil: every entry lies within the extents'
       end if
    end if

    this % rows      = rows
    this % columns   = columns
    this % weights   = weights
    this % constants = constant

    if (present(row_domain)) then
       this % row_domain  = row_domain
       this % row_entries = this % num_rows
       if (present(row_entries)) this % row_entries = row_entries
    end if
    if (present(column_domain)) then
       this % column_domain  = column_domain
       this % column_entries = this % num_columns
       if (present(column_entries)) this % column_entries = column_entries
    end if

    if (present(label)) then
       this % label = label
    else
       this % label = 'stencil'
    end if

    ! one argument: the state the matrix multiplies
    call this % declare_arguments(1)

  end function create

  !===================================================================!
  ! Build from a dense matrix: each entry becomes one weighted
  ! triple, constants zero. A rectangular array is lawful; its
  ! column count is the second extent.
  !===================================================================!

  function create_dense(a, label) result(this)

    real(dp)        , intent(in)           :: a(:,:)
    character(len=*), intent(in), optional :: label
    type(stencil) :: this

    integer , allocatable :: rows(:), columns(:)
    real(dp), allocatable :: weights(:), constant(:)
    integer :: n, m, i, j, e

    n = size(a, 1)
    m = size(a, 2)

    allocate(rows(n * m), columns(n * m), weights(n * m), constant(n))
    constant = 0.0_dp

    e = 0
    do j = 1, m
       do i = 1, n
          e          = e + 1
          rows(e)    = i
          columns(e) = j
          weights(e) = a(i, j)
       end do
    end do

    this = create(rows, columns, weights, constant, label, num_columns=m)

  end function create_dense

  !===================================================================!
  ! Build from an operation by evaluation on the standard basis: the
  ! operation applied to the zero state is the constant, and applied
  ! to each basis vector minus that constant is one column. width is
  ! the number of values the input carries, on the operation's own
  ! domain unless a column domain and its member count are given -
  ! the rectangular case, a tangent in an argument of another
  ! domain. The rows are the operation's domain, as many values per
  ! member as the zero state reports. A width that does not divide
  ! among the column members, or a result whose length changes
  ! between applications, stops the program, because that column
  ! cannot be placed in the matrix. The label defaults to the
  ! operation's name.
  !===================================================================!

  function create_compiled(action, on, width, label, column_domain, &
       & column_entries) result(this)

    class(operation)     , intent(in) :: action
    class(directed_graph), intent(in) :: on
    integer              , intent(in) :: width
    character(len=*), intent(in), optional :: label
    type(graph)     , intent(in), optional :: column_domain
    integer         , intent(in), optional :: column_entries
    type(stencil) :: this

    type(stored_field)        :: state
    class(field), allocatable :: output
    type(graph) :: dom, along
    real(dp), allocatable :: a(:,:), e(:), y(:), constant(:)
    integer :: n_dom, n_along, num_components, j

    call action % domain(on, dom, n_dom)

    if (n_dom <= 0) then
       error stop 'stencil: the operation''s domain is nonempty'
    end if

    along   = dom
    n_along = n_dom
    if (present(column_domain)) then
       along   = column_domain
       n_along = width
       if (present(column_entries)) n_along = column_entries
    end if
    if (width <= 0 .or. n_along <= 0 .or. mod(width, n_along) /= 0) then
       error stop 'stencil: the width carries a whole number per member'
    end if

    num_components = width / n_along
    allocate(e(width))

    ! j = 0 is the zero state, whose value is the constant and whose
    ! length is the row count
    do j = 0, width
       e = 0.0_dp
       if (j > 0) e(j) = 1.0_dp
       state = stored_field('basis', along, n_along, num_components=num_components)
       call state % set_real_vector(e)
       call action % apply(on, [state], output)
       call output % real_vector(y)
       if (j == 0) then
          constant = y
          if (mod(size(y), n_dom) /= 0) then
             error stop 'stencil: the operation result carries a whole number per member'
          end if
          allocate(a(size(y), width))
       else
          if (size(y) /= size(constant)) then
             error stop 'stencil: the operation result keeps its length'
          end if
          a(:, j) = y - constant
       end if
    end do

    if (present(label)) then
       this = create_dense(a, label)
    else
       this = create_dense(a, action % name())
    end if
    this % constants = constant

    this % row_domain     = dom
    this % row_entries    = n_dom
    this % column_domain  = along
    this % column_entries = n_along

  end function create_compiled

  pure function stencil_name(this) result(name)

    class(stencil), intent(in) :: this
    character(len=:), allocatable :: name

    name = this % label

  end function stencil_name

  !===================================================================!
  ! Where a result lives: the declared row domain, one value per
  ! row; without one, the host's vertex set, the rows divided among
  ! its vertices as components. Rows that do not divide among the
  ! host's vertices stop the program, because no field could hold
  ! them.
  !===================================================================!

  subroutine result_extent(this, input_graph, domain, num_entries, &
       & num_components)

    class(stencil), intent(in)        :: this
    class(directed_graph), intent(in) :: input_graph
    type(graph), intent(out)          :: domain
    integer    , intent(out)          :: num_entries
    integer    , intent(out)          :: num_components

    type(token) :: key

    key = this % row_domain % id()
    if (key % declared()) then
       domain      = this % row_domain
       num_entries = this % row_entries
    else
       domain      = input_graph % vertex_set()
       num_entries = input_graph % num_vertices()
    end if

    if (num_entries <= 0 .or. mod(this % num_rows, num_entries) /= 0) then
       error stop 'stencil: the rows carry a whole number of values per member &
            &of the row domain'
    end if
    num_components = this % num_rows / num_entries

  end subroutine result_extent

  subroutine stencil_domain(this, input_graph, domain, num_entries)

    class(stencil), intent(in)    :: this
    class(directed_graph), intent(in)               :: input_graph
    type(graph), intent(out) :: domain
    integer        , intent(out) :: num_entries

    integer :: num_components

    call result_extent(this, input_graph, domain, num_entries, num_components)

  end subroutine stencil_domain

  !===================================================================!
  ! The result field on the row extent, holding y.
  !===================================================================!

  subroutine result_field(this, input_graph, y, output)

    class(stencil), intent(in)               :: this
    class(directed_graph), intent(in)        :: input_graph
    real(dp), intent(in)                     :: y(:)
    class(field), allocatable, intent(inout) :: output

    type(stored_field) :: out
    type(graph) :: domain
    integer :: num_entries, num_components

    call result_extent(this, input_graph, domain, num_entries, num_components)

    out = stored_field(this % label, domain, num_entries, &
         & num_components=num_components)
    call out % set_real_vector(y)

    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine result_field

  !===================================================================!
  ! Check that an input lives on the column domain, by identity,
  ! whenever the stencil declares one; an input of the right length
  ! on another set is not the same input. Violation stops the
  ! program. Without a declared column domain only the length is
  ! checked, in accumulate.
  !===================================================================!

  subroutine require_column_domain(this, domain)

    class(stencil), intent(in) :: this
    type(graph)   , intent(in) :: domain

    type(token) :: key

    key = this % column_domain % id()
    if (.not. key % declared()) return

    if (.not. domain % same_as(this % column_domain)) then
       error stop 'stencil: the input lives on the column domain'
    end if

  end subroutine require_column_domain

  !===================================================================!
  ! y = constants + A q, the entries traversed once: each carries
  ! its weight times the column's value onto its row.
  !===================================================================!

  subroutine stencil_apply(this, input_graph, input_data, output)

    class(stencil), intent(in)            :: this
    class(directed_graph), intent(in)                       :: input_graph
    class(field), intent(in), optional       :: input_data(:)
    class(field), allocatable, intent(inout) :: output

    real(dp), allocatable :: q(:), y(:)

    y = this % constants

    if (present(input_data)) then
       call require_column_domain(this, input_data(1) % domain())
       call input_data(1) % real_vector(q)
       call accumulate(this, q, y)
    end if

    call result_field(this, input_graph, y, output)

  end subroutine stencil_apply

  !===================================================================!
  ! The one entry traversal, shared by apply and the tangent. An
  ! input that does not fill the columns stops the program, because
  ! an entry would read past it.
  !===================================================================!

  pure subroutine accumulate(this, q, y)

    class(stencil), intent(in)    :: this
    real(dp)      , intent(in)    :: q(:)
    real(dp)      , intent(inout) :: y(:)

    integer :: e

    if (size(q) /= this % num_columns) then
       error stop 'stencil: the input fills the columns'
    end if

    do e = 1, size(this % rows)
       y(this % rows(e)) = y(this % rows(e)) &
            & + this % weights(e) * q(this % columns(e))
    end do

  end subroutine accumulate

  !===================================================================!
  ! A stencil is linear, so it is its own tangent: the first partial
  ! action in its one input slot is the entry traversal on the
  ! direction, without the constants. An order past one or a slot
  ! other than one stops the program, because a linear map in one
  ! slot has no other partial.
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

    real(dp), allocatable :: v(:), y(:)

    associate (u1 => input_data); end associate

    call this % require_owned(variations)

    if (size(variations) /= 1) then
       error stop 'stencil: the requested order is within max_degree'
    end if
    if (.not. variations(1) % argument_is(this % argument(1))) then
       error stop 'stencil: the partial action is taken in the one argument'
    end if

    call require_column_domain(this, variations(1) % domain())
    call variations(1) % direction(v)
    allocate(y(this % num_rows))
    y = 0.0_dp
    call accumulate(this, v, y)

    call result_field(this, input_graph, y, output)

  end subroutine stencil_partial_action

  !===================================================================!
  ! The dependency graph, derived when asked: one vertex per row,
  ! one directed edge per entry, column to row. Only a square map
  ! has one - its rows and columns are one set; a rectangular map
  ! stops the program here.
  !===================================================================!

  subroutine stencil_dependencies(this, pattern)

    class(stencil), intent(in)    :: this
    class(directed_graph), allocatable, intent(out) :: pattern

    if (this % num_rows /= this % num_columns) then
       error stop 'stencil: a dependency graph is square'
    end if

    allocate(pattern, source=stored_directed_graph(this % num_rows, &
         & tails=this % columns, heads=this % rows))

  end subroutine stencil_dependencies

  !===================================================================!
  ! The transpose: rows and columns exchanged, the domains with
  ! them, the constants dropped - the affine part of a map has no
  ! transpose.
  !===================================================================!

  function stencil_transpose(this) result(transposed)

    class(stencil), intent(in) :: this
    type(stencil) :: transposed

    real(dp), allocatable :: zeros(:)

    allocate(zeros(this % num_columns))
    zeros = 0.0_dp

    transposed = create(rows=this % columns, columns=this % rows, &
         & weights=this % weights, constant=zeros, &
         & label='transpose of ' // this % label, &
         & num_columns=this % num_rows)

    transposed % row_domain     = this % column_domain
    transposed % row_entries    = this % column_entries
    transposed % column_domain  = this % row_domain
    transposed % column_entries = this % row_entries

  end function stencil_transpose

  !===================================================================!
  ! The composition outer o inner, by sparse triple product: one
  ! product entry per (entry of outer, entry of inner in the row it
  ! reads), equal (row, column) pairs combined, the inner constant
  ! carried through the outer map. The column count of outer must
  ! equal the row count of inner; a mismatch means a chain was
  ! assembled wrong and stops the program. The result reads inner's
  ! columns and lands on outer's rows.
  !===================================================================!

  function compose(outer, inner) result(composed)

    type(stencil), intent(in) :: outer
    type(stencil), intent(in) :: inner
    type(stencil) :: composed

    integer , allocatable :: ptr(:), order(:), identity(:), r(:), c(:)
    integer , allocatable :: rows(:), columns(:)
    real(dp), allocatable :: w(:), weights(:), constant(:)
    integer :: k2, j, n, row1

    if (outer % num_columns /= inner % num_rows) then
       error stop 'stencil: composed maps agree on the inner extent'
    end if

    ! group inner's entries by row with the one counting sort: order
    ! lists inner's entry indices row by row, in arrival order, and
    ! ptr(row)..ptr(row+1)-1 is each row's range
    allocate(identity(size(inner % rows)))
    identity = [(j, j = 1, size(inner % rows))]
    call group_by_key(inner % num_rows, inner % rows, identity, ptr, order)

    n = 0
    do k2 = 1, size(outer % rows)
       row1 = outer % columns(k2)
       n = n + ptr(row1 + 1) - ptr(row1)
    end do

    allocate(r(n), c(n), w(n))
    n = 0
    do k2 = 1, size(outer % rows)
       row1 = outer % columns(k2)
       do j = ptr(row1), ptr(row1 + 1) - 1
          n = n + 1
          r(n) = outer % rows(k2)
          c(n) = inner % columns(order(j))
          w(n) = outer % weights(k2) * inner % weights(order(j))
       end do
    end do

    call combine_triples(outer % num_rows, inner % num_columns, r, c, w, &
         & rows, columns, weights)

    ! the inner constant travels through this map
    constant = outer % constants
    do k2 = 1, size(outer % rows)
       constant(outer % rows(k2)) = constant(outer % rows(k2)) &
            & + outer % weights(k2) * inner % constants(outer % columns(k2))
    end do

    composed = create(rows, columns, weights, constant, &
         & label=outer % label // ' o ' // inner % label, &
         & num_columns=inner % num_columns)

    composed % row_domain     = outer % row_domain
    composed % row_entries    = outer % row_entries
    composed % column_domain  = inner % column_domain
    composed % column_entries = inner % column_entries

  end function compose

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
