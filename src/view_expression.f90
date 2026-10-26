!=====================================================================!
! Module that writes the printed forms of an expression.
!
! The vertex table of an expression is read once and written five
! ways: as a tree from the root down, as a diagram with operands
! centred under their operator, as a formula in infix notation, as
! nested brackets, and as a digraph in the DOT language. For a table
!
!      1  q1        2  x        3  sin(2)        4  *(1, 3)
!
! the five read:
!
!      tree             diagram          formula      brackets
!         4  *             *             q1*sin(x)    [* [q1] [sin [x]]]
!         1  ├── q1     ┌──┴──┐
!         3  └── sin    q1   sin         digraph
!         2      └── x        │          v3 -> v2;  v4 -> v1;  v4 -> v3;
!                             x
!
! Leaves. A leaf is written by the application's own component names,
! a constant by its value and an operator by its symbol:
!
!      a field's value, field f       q<f>
!      its derivatives                one t per order along the instants,
!                                     one x, y or z per spatial order: q1tx
!      the design                     nu
!      a multiplier j                 lambda<j>
!      a coordinate                   t, x, y, z, then c<n>
!      a time map along a coordinate  x_map
!      a constant                     2.0, 0.333333
!      an operator                    +  -  *  /  **2  sin  dx
!
! A view given the names of the unknowns and of the design coordinates
! writes those in place of q<f> and nu. The vertex index printed is
! the one root() reports.
!
! The composition of views. Every vertex is written when its view is
! formed, so that views of expressions named differently compose into
! one: the sum and the product of two views, and a vertex of text over
! no view (a leaf), over one (an operator written before its operand,
! "integral over t of"), or over two (a separator written after the
! first operand, a comma). The vertices of the operands come first,
! the new root last:
!
!      first       1 .. n1
!      second      1 .. n2      renumbered n1+1 .. n1+n2
!      composed    1 .. n1 | n1+1 .. n1+n2 | n1+n2+1
!                  first's   second's        the new root
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module view_expression

  use util_precision      , only : dp
  use util_string         , only : string
  use operation_expression, only : expression, expression_vertex, FIRST_COORDINATE
  use operation_expression, only : VERTEX_LEAF, VERTEX_CONSTANT, VERTEX_SUM, VERTEX_DIFFERENCE
  use operation_expression, only : VERTEX_PRODUCT, VERTEX_QUOTIENT, VERTEX_INTEGER_POWER
  use operation_expression, only : VERTEX_REAL_POWER, VERTEX_FUNCTION, VERTEX_DERIVATIVE
  use operation_expression, only : ARGUMENT_DESIGN, ARGUMENT_MULTIPLIER, ARGUMENT_COORDINATE, ARGUMENT_TIME_MAP
  use operation_expression, only : SINE, COSINE, EXPONENTIAL, LOGARITHM, SQUARE_ROOT, UNIT_STEP

  implicit none

  private

  public :: expression_view, operator(+), operator(*)

  ! A vertex of a view that no expression has, written by its text;
  ! the kinds of an expression's vertices are positive
  integer, parameter :: VERTEX_TEXT = 0

  !-------------------------------------------------------------------!
  ! Type that models a view: the vertex table of an expression and
  ! the text of each vertex, written once when the view is formed
  !-------------------------------------------------------------------!

  type :: expression_view

     type(expression_vertex), allocatable, private :: v(:)
     type(string)           , allocatable, private :: label(:)

   contains

     ! The five printed forms
     procedure :: tree
     procedure :: diagram
     procedure :: formula
     procedure :: brackets
     procedure :: digraph

  end type expression_view

  interface expression_view
     module procedure create
     module procedure text_view
  end interface expression_view

  interface operator(+)
     module procedure view_sum
  end interface operator(+)

  interface operator(*)
     module procedure view_product
  end interface operator(*)

contains

  !===================================================================!
  ! Creates the view of an expression, its leaves named by the names
  ! of the unknowns and of the design coordinates where given
  !===================================================================!

  function create(e, unknowns, designs) result(this)

    type(expression)     , intent(in)           :: e
    character(len=*)     , intent(in), optional :: unknowns(:), designs(:)
    type(expression_view)                       :: this

    integer :: i

    this % v = e % vertices()
    allocate(this % label(size(this % v)))
    do i = 1, size(this % v)
       this % label(i) = string(vertex_label(this % v(i), unknowns, designs))
    end do

  end function create

  !===================================================================!
  ! Creates a vertex of text over no view, one or two
  !===================================================================!

  function text_view(text, first, second) result(this)

    character(len=*)     , intent(in)           :: text
    type(expression_view), intent(in), optional :: first, second
    type(expression_view)                       :: this

    this = composition(VERTEX_TEXT, text, first, second)

  end function text_view

  !===================================================================!
  ! Returns the sum of two views
  !===================================================================!

  function view_sum(a, b) result(this)

    type(expression_view), intent(in) :: a, b
    type(expression_view)             :: this

    this = composition(VERTEX_SUM, '+', a, b)

  end function view_sum

  !===================================================================!
  ! Returns the product of two views
  !===================================================================!

  function view_product(a, b) result(this)

    type(expression_view), intent(in) :: a, b
    type(expression_view)             :: this

    this = composition(VERTEX_PRODUCT, '*', a, b)

  end function view_product

  !===================================================================!
  ! Returns the vertices of first, those of second numbered after
  ! them, and a root of the kind and text given over the two roots
  !===================================================================!

  function composition(kind, text, first, second) result(this)

    integer              , intent(in)           :: kind
    character(len=*)     , intent(in)           :: text
    type(expression_view), intent(in), optional :: first, second
    type(expression_view)                       :: this

    type(expression_vertex)              :: root
    type(expression_vertex), allocatable :: w(:)
    integer                              :: n

    root % kind = kind
    allocate(this % v(0), this % label(0))
    if (present(first)) then
       this % v     = first % v
       this % label = first % label
       root % first = size(this % v)
    end if
    if (present(second)) then
       n = size(this % v)
       w = second % v
       where (w % first  > 0) w % first  = w % first  + n
       where (w % second > 0) w % second = w % second + n
       this % v      = [this % v, w]
       this % label  = [this % label, second % label]
       root % second = size(this % v)
    end if
    this % v     = [this % v, root]
    this % label = [this % label, string(text)]

  end function composition

  !===================================================================!
  ! Returns the text of one vertex of an expression
  !===================================================================!

  function vertex_label(u, unknowns, designs) result(text)

    type(expression_vertex), intent(in)           :: u
    character(len=*)       , intent(in), optional :: unknowns(:), designs(:)
    character(len=:)       , allocatable          :: text

    select case (u % kind)
    case (VERTEX_LEAF)
       select case (u % position)
       case (ARGUMENT_DESIGN)
          text = leaf_name(designs, u % field, 'nu')
       case (ARGUMENT_MULTIPLIER)
          text = 'lambda' // written(u % field)
       case (ARGUMENT_COORDINATE)
          text = axis(u % field)
       case (ARGUMENT_TIME_MAP)
          text = axis(u % field) // '_map'
       case default
          text = leaf_name(unknowns, u % field, 'q' // written(u % field)) &
               & // repeat(axis(FIRST_COORDINATE), u % mixed_order) // repeat(axis(u % along), u % order)
       end select
    case (VERTEX_CONSTANT)
       text = number(u % coefficient)
    case (VERTEX_SUM)
       text = '+'
    case (VERTEX_DIFFERENCE)
       text = '-'
    case (VERTEX_PRODUCT)
       text = '*'
    case (VERTEX_QUOTIENT)
       text = '/'
    case (VERTEX_INTEGER_POWER)
       text = '**' // written(nint(u % coefficient))
    case (VERTEX_REAL_POWER)
       text = '**' // number(u % coefficient)
    case (VERTEX_FUNCTION)
       select case (u % order)
       case (SINE);        text = 'sin'
       case (COSINE);      text = 'cos'
       case (EXPONENTIAL); text = 'exp'
       case (LOGARITHM);   text = 'log'
       case (SQUARE_ROOT); text = 'sqrt'
       case (UNIT_STEP);   text = 'heaviside'
       case default
          error stop 'view_expression: this vertex names a function that is not one of those defined'
       end select
    case (VERTEX_DERIVATIVE)
       text = 'd' // repeat(axis(u % along), u % order)
    case default
       error stop 'view_expression: this vertex names a kind that is not one of those defined'
    end select

  end function vertex_label

  !===================================================================!
  ! Returns the name of index i among the names given; the default
  ! where no names are given or i is beyond them
  !===================================================================!

  pure function leaf_name(names, i, default) result(text)

    character(len=*), intent(in), optional :: names(:)
    integer         , intent(in)           :: i
    character(len=*), intent(in)           :: default
    character(len=:), allocatable          :: text

    text = default
    if (.not. present(names)) return
    if (i <= size(names)) text = trim(names(i))

  end function leaf_name

  !===================================================================!
  ! Returns the letter of a coordinate: t for the first, x, y, z for
  ! the next three, c<n> beyond
  !===================================================================!

  pure function axis(coordinate) result(letter)

    integer         , intent(in)  :: coordinate
    character(len=:), allocatable :: letter

    select case (coordinate - FIRST_COORDINATE)
    case (0);      letter = 't'
    case (1);      letter = 'x'
    case (2);      letter = 'y'
    case (3);      letter = 'z'
    case default;  letter = 'c' // written(coordinate)
    end select

  end function axis

  !===================================================================!
  ! Returns an integer written without padding
  !===================================================================!

  pure function written(n) result(text)

    integer         , intent(in)  :: n
    character(len=:), allocatable :: text

    character(len=12) :: buffer

    write(buffer, '(i0)') n
    text = trim(buffer)

  end function written

  !===================================================================!
  ! Returns a real written as an integral value with one decimal, 2.0,
  ! or any other to six digits, 0.333333
  !===================================================================!

  pure function number(x) result(text)

    real(dp)        , intent(in)  :: x
    character(len=:), allocatable :: text

    character(len=32) :: buffer

    if (abs(x) < 1.0e9_dp .and. x == aint(x)) then
       write(buffer, '(i0,a)') nint(x), '.0'
    else
       write(buffer, '(g0.6)') x
    end if
    text = trim(adjustl(buffer))

  end function number

  !===================================================================!
  ! Returns the tree from the root down, one vertex per line, each
  ! operand indented below its operator, the vertex index at the left
  !===================================================================!

  function tree(this) result(text)

    class(expression_view), intent(in)  :: this
    character(len=:)      , allocatable :: text

    text = ''
    call tree_lines(this, size(this % v), '', '', '', text)
    text = text(1:len(text) - 1)

  end function tree

  !===================================================================!
  ! Appends the line of vertex i and those of its operands below it.
  ! The branch is the connector of this line; the continuation is the
  ! prefix of its operands' lines, a bar where a later sibling
  ! follows and spaces where none does.
  !===================================================================!

  recursive subroutine tree_lines(this, i, prefix, branch, continuation, text)

    class(expression_view)       , intent(in)    :: this
    integer                      , intent(in)    :: i
    character(len=*)             , intent(in)    :: prefix, branch, continuation
    character(len=:), allocatable, intent(inout) :: text

    character(len=6) :: column

    write(column, '(i4,2x)') i
    text = text // column // prefix // branch // this % label(i) % str // new_line('a')
    if (this % v(i) % second > 0) then
       call tree_lines(this, this % v(i) % first,  continuation, '├── ', continuation // '│   ', text)
       call tree_lines(this, this % v(i) % second, continuation, '└── ', continuation // '    ', text)
    else if (this % v(i) % first > 0) then
       call tree_lines(this, this % v(i) % first,  continuation, '└── ', continuation // '    ', text)
    end if

  end subroutine tree_lines

  !===================================================================!
  ! Returns the diagram from the root down, operands centred under
  ! their operator. The drawing is built as cells of one visual column
  ! each, so that the box characters, three bytes wide, align with the
  ! letters.
  !===================================================================!

  function diagram(this) result(text)

    class(expression_view), intent(in)  :: this
    character(len=:)      , allocatable :: text

    character(len=3), allocatable :: cells(:,:)
    integer                       :: centre, r, c, n

    call drawn(this, size(this % v), cells, centre)
    text = ''
    do r = 1, size(cells, 1)
       do c = 1, size(cells, 2)
          n = len_trim(cells(r, c))
          if (n == 0) then
             text = text // ' '
          else
             text = text // cells(r, c)(1:n)
          end if
       end do
       if (r < size(cells, 1)) text = text // new_line('a')
    end do

  end function diagram

  !===================================================================!
  ! Returns the cells of the drawing of vertex i and the column of its
  ! centre. A leaf is its text. One operand is drawn below its
  ! operator, joined by a bar. Two are drawn side by side, three
  ! columns apart, and the operator is centred over the midpoint of
  ! their centres:
  !
  !            op            the operator's text, centred on m
  !         ┌──┴──┐          m = (cl + cr) / 2
  !       left   right       gap of three columns between the blocks
  !
  ! Where the text is wider than what it is centred over, every column
  ! is shifted right by the offset that keeps the first in the drawing.
  !===================================================================!

  recursive subroutine drawn(this, i, cells, centre)

    class(expression_view)       , intent(in)  :: this
    integer                      , intent(in)  :: i
    character(len=3), allocatable, intent(out) :: cells(:,:)
    integer                      , intent(out) :: centre

    character(len=3), allocatable :: left(:,:), right(:,:), below(:,:)
    character(len=:), allocatable :: text
    integer, parameter            :: gap = 3
    integer                       :: cl, cr, wl, wr, hl, hr, w, h, m, k, offset, start

    text = this % label(i) % str

    if (this % v(i) % first == 0) then
       allocate(cells(1, len(text)))
       do k = 1, len(text)
          cells(1, k) = text(k:k)
       end do
       centre = (len(text) + 1) / 2
       return
    end if

    if (this % v(i) % second == 0) then
       call drawn(this, this % v(i) % first, below, cl)
       start  = cl - (len(text) - 1) / 2
       offset = max(0, 1 - start)
       w = max(size(below, 2), start + len(text) - 1) + offset
       h = size(below, 1) + 2
       allocate(cells(h, w))
       cells = ' '
       do k = 1, len(text)
          cells(1, start + offset + k - 1) = text(k:k)
       end do
       cells(2, cl + offset) = '│'
       cells(3:, 1 + offset:size(below, 2) + offset) = below
       centre = cl + offset
       return
    end if

    call drawn(this, this % v(i) % first,  left,  cl)
    call drawn(this, this % v(i) % second, right, cr)
    wl = size(left, 2)
    wr = size(right, 2)
    hl = size(left, 1)
    hr = size(right, 1)
    cr = cr + wl + gap
    m  = (cl + cr) / 2
    start  = m - (len(text) - 1) / 2
    offset = max(0, 1 - start)
    w = max(wl + gap + wr, start + len(text) - 1) + offset
    h = max(hl, hr) + 2
    allocate(cells(h, w))
    cells = ' '
    do k = 1, len(text)
       cells(1, start + offset + k - 1) = text(k:k)
    end do
    do k = cl + 1, cr - 1
       cells(2, k + offset) = '─'
    end do
    cells(2, cl + offset) = '┌'
    cells(2, m + offset)  = '┴'
    cells(2, cr + offset) = '┐'
    cells(3:2 + hl, 1 + offset:wl + offset) = left
    cells(3:2 + hr, wl + gap + 1 + offset:wl + gap + wr + offset) = right
    centre = m + offset

  end subroutine drawn

  !===================================================================!
  ! Returns the formula in infix notation, parenthesized only where
  ! the tree requires it: a sum or difference under a product, quotient
  ! or power; a product or quotient under a power; the right operand of
  ! a difference or quotient at its own level. A vertex of text with an
  ! operand binds below every operator, so that it is parenthesized
  ! under each; its operand is parenthesized as that of a product, a
  ! separator's operands never.
  !
  ! Each vertex has a level, and each operand is written with the
  ! level it must reach; one below it is parenthesized:
  !
  !      level   vertex                  operand must reach
  !        5     leaf, constant          -
  !        4     function, derivative    0, already in parentheses
  !        3     power                   4
  !        2     product, quotient       2; the right of a quotient 3
  !        1     sum, difference         1; the right of a difference 2
  !        0     text over operands      2 for one, 0 for a separator's
  !
  !      (q1 + q2)*x        q1 - (q2 - x)        (q1*x)**2
  !===================================================================!

  function formula(this) result(text)

    class(expression_view), intent(in)  :: this
    character(len=:)      , allocatable :: text

    text = infix(this, size(this % v), 0)

  end function formula

  !===================================================================!
  ! Returns vertex i in infix notation, parenthesized where its level
  ! is below above. The operands are written into locals before the
  ! concatenation: gfortran overwrites the result of a nested call to a
  ! recursive function of deferred length within one expression.
  !===================================================================!

  recursive function infix(this, i, above) result(text)

    class(expression_view), intent(in)  :: this
    integer               , intent(in)  :: i, above
    character(len=:)      , allocatable :: text

    character(len=:), allocatable :: first, second, operator
    integer                       :: level

    operator = this % label(i) % str
    select case (this % v(i) % kind)
    case (VERTEX_LEAF, VERTEX_CONSTANT)
       level = 5
       text  = operator
    case (VERTEX_TEXT)
       level = 0
       text  = operator
       if (this % v(i) % second > 0) then
          first  = infix(this, this % v(i) % first, 0)
          second = infix(this, this % v(i) % second, 0)
          text   = first // operator // ' ' // second
       else if (this % v(i) % first > 0) then
          first = infix(this, this % v(i) % first, 2)
          text  = operator // ' ' // first
       else
          level = 5
       end if
    case (VERTEX_FUNCTION, VERTEX_DERIVATIVE)
       level = 4
       first = infix(this, this % v(i) % first, 0)
       text  = operator // '(' // first // ')'
    case (VERTEX_INTEGER_POWER, VERTEX_REAL_POWER)
       level = 3
       first = infix(this, this % v(i) % first, 4)
       text  = first // operator
    case (VERTEX_PRODUCT, VERTEX_QUOTIENT)
       level  = 2
       first  = infix(this, this % v(i) % first, 2)
       second = infix(this, this % v(i) % second, merge(3, 2, this % v(i) % kind == VERTEX_QUOTIENT))
       text   = first // operator // second
    case default
       level  = 1
       first  = infix(this, this % v(i) % first, 1)
       second = infix(this, this % v(i) % second, merge(2, 1, this % v(i) % kind == VERTEX_DIFFERENCE))
       text   = first // ' ' // operator // ' ' // second
    end select
    if (level < above) text = '(' // text // ')'

  end function infix

  !===================================================================!
  ! Returns the brackets: each vertex as [label operand operand], the
  ! input syntax of a forest figure
  !===================================================================!

  function brackets(this) result(text)

    class(expression_view), intent(in)  :: this
    character(len=:)      , allocatable :: text

    text = bracketed(this, size(this % v))

  end function brackets

  !===================================================================!
  ! Returns vertex i and its operands in brackets
  !===================================================================!

  recursive function bracketed(this, i) result(text)

    class(expression_view), intent(in)  :: this
    integer               , intent(in)  :: i
    character(len=:)      , allocatable :: text

    character(len=:), allocatable :: operand

    text = '[' // this % label(i) % str
    if (this % v(i) % first > 0) then
       operand = bracketed(this, this % v(i) % first)
       text = text // ' ' // operand
    end if
    if (this % v(i) % second > 0) then
       operand = bracketed(this, this % v(i) % second)
       text = text // ' ' // operand
    end if
    text = text // ']'

  end function bracketed

  !===================================================================!
  ! Returns the digraph in the DOT language: one node per vertex,
  ! named by its index and labelled by its text, then one arc from each
  ! operator to each operand
  !
  !      digraph expression {
  !        v1 [label="q1"];
  !        ...
  !        v4 -> v1;
  !      }
  !===================================================================!

  function digraph(this) result(text)

    class(expression_view), intent(in)  :: this
    character(len=:)      , allocatable :: text

    integer :: i

    text = 'digraph expression {' // new_line('a')
    do i = 1, size(this % v)
       text = text // '  v' // written(i) // ' [label="' // this % label(i) % str // '"];' // new_line('a')
    end do
    do i = 1, size(this % v)
       if (this % v(i) % first > 0) then
          text = text // '  v' // written(i) // ' -> v' // written(this % v(i) % first) // ';' // new_line('a')
       end if
       if (this % v(i) % second > 0) then
          text = text // '  v' // written(i) // ' -> v' // written(this % v(i) % second) // ';' // new_line('a')
       end if
    end do
    text = text // '}'

  end function digraph

end module view_expression
