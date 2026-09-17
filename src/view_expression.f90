!=====================================================================!
! THE PRINTED FORMS OF AN EXPRESSION. The vertex table of an
! expression is read once and written as a tree from the root down,
! as a diagram with operands centred under their operator, as a
! formula in infix notation, as nested brackets, and as a digraph in
! the DOT language.
!
! A leaf is written by the application's own component names: q<f>
! for a field's value, one t per order along the instants, one x, y
! or z per order along a spatial coordinate, nu for the design,
! lambda<j> for a multiplier; a constant by its value; an operator by
! its symbol. The vertex index printed is the one root() reports.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module view_expression

  use util_precision      , only : dp
  use operation_expression, only : expression, expression_vertex, FIRST_COORDINATE
  use operation_expression, only : VERTEX_LEAF, VERTEX_CONSTANT, VERTEX_SUM, VERTEX_DIFFERENCE
  use operation_expression, only : VERTEX_PRODUCT, VERTEX_QUOTIENT, VERTEX_INTEGER_POWER
  use operation_expression, only : VERTEX_REAL_POWER, VERTEX_FUNCTION
  use operation_expression, only : ARGUMENT_DESIGN, ARGUMENT_MULTIPLIER, ARGUMENT_COORDINATE
  use operation_expression, only : SINE, COSINE, EXPONENTIAL, LOGARITHM, SQUARE_ROOT

  implicit none

  private
  public :: expression_view

  type :: expression_view

     type(expression_vertex), allocatable, private :: v(:)

   contains

     procedure :: tree
     procedure :: diagram
     procedure :: formula
     procedure :: brackets
     procedure :: digraph
     procedure, private :: label

  end type expression_view

  interface expression_view
     module procedure create
  end interface expression_view

contains

  function create(e) result(this)

    type(expression), intent(in) :: e
    type(expression_view) :: this

    this % v = e % vertices()

  end function create

  !===================================================================!
  ! The text of one vertex.
  !===================================================================!

  function label(this, i) result(text)

    class(expression_view), intent(in) :: this
    integer               , intent(in) :: i
    character(len=:), allocatable :: text

    associate (u => this % v(i))
      select case (u % kind)
      case (VERTEX_LEAF)
         select case (u % position)
         case (ARGUMENT_DESIGN)
            text = 'nu'
         case (ARGUMENT_MULTIPLIER)
            text = 'lambda' // written(u % field)
         case (ARGUMENT_COORDINATE)
            text = axis(u % field)
         case default
            text = 'q' // written(u % field) // repeat(axis(u % along), u % order)
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
         case default
            error stop 'view_expression: this vertex names a function that is not one of those defined'
         end select
      case default
         error stop 'view_expression: this vertex names a kind that is not one of those defined'
      end select
    end associate

  end function label

  ! the letter of a coordinate: t for the first, x, y, z for the next
  ! three, c<n> beyond
  pure function axis(coordinate) result(letter)

    integer, intent(in) :: coordinate
    character(len=:), allocatable :: letter

    select case (coordinate - FIRST_COORDINATE)
    case (0);      letter = 't'
    case (1);      letter = 'x'
    case (2);      letter = 'y'
    case (3);      letter = 'z'
    case default;  letter = 'c' // written(coordinate)
    end select

  end function axis

  pure function written(n) result(text)

    integer, intent(in) :: n
    character(len=:), allocatable :: text

    character(len=12) :: buffer

    write(buffer, '(i0)') n
    text = trim(buffer)

  end function written

  ! an integral value with one decimal, any other to six digits
  pure function number(x) result(text)

    real(dp), intent(in) :: x
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
  ! THE TREE from the root down, one vertex per line, each operand
  ! indented below its operator, the vertex index at the left.
  !===================================================================!

  function tree(this) result(text)

    class(expression_view), intent(in) :: this
    character(len=:), allocatable :: text

    text = ''
    call tree_lines(this, size(this % v), '', '', '', text)
    text = text(1:len(text) - 1)

  end function tree

  recursive subroutine tree_lines(this, i, prefix, branch, continuation, text)

    class(expression_view), intent(in)    :: this
    integer               , intent(in)    :: i
    character(len=*)      , intent(in)    :: prefix, branch, continuation
    character(len=:), allocatable, intent(inout) :: text

    character(len=6) :: column

    write(column, '(i4,2x)') i
    text = text // column // prefix // branch // this % label(i) // new_line('a')
    if (this % v(i) % second > 0) then
       call tree_lines(this, this % v(i) % first,  continuation, '├── ', continuation // '│   ', text)
       call tree_lines(this, this % v(i) % second, continuation, '└── ', continuation // '    ', text)
    else if (this % v(i) % first > 0) then
       call tree_lines(this, this % v(i) % first,  continuation, '└── ', continuation // '    ', text)
    end if

  end subroutine tree_lines

  !===================================================================!
  ! THE DIAGRAM from the root down, operands centred under their
  ! operator. The drawing is built as cells of one visual column each,
  ! so that the box characters, three bytes wide, align with the
  ! letters.
  !===================================================================!

  function diagram(this) result(text)

    class(expression_view), intent(in) :: this
    character(len=:), allocatable :: text

    character(len=3), allocatable :: cells(:,:)
    integer :: centre, r, c, n

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

  recursive subroutine drawn(this, i, cells, centre)

    class(expression_view), intent(in)  :: this
    integer               , intent(in)  :: i
    character(len=3), allocatable, intent(out) :: cells(:,:)
    integer               , intent(out) :: centre

    character(len=3), allocatable :: left(:,:), right(:,:), below(:,:)
    character(len=:), allocatable :: text
    integer, parameter :: gap = 3
    integer :: cl, cr, wl, wr, hl, hr, w, h, m, k, offset, start

    text = this % label(i)

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
  ! THE FORMULA in infix notation, parenthesized only where the tree
  ! requires it: a sum or difference under a product, quotient or
  ! power; a product or quotient under a power; the right operand of
  ! a difference or quotient at its own level.
  !===================================================================!

  function formula(this) result(text)

    class(expression_view), intent(in) :: this
    character(len=:), allocatable :: text

    text = infix(this, size(this % v), 0)

  end function formula

  ! the operands are written into locals before the concatenation:
  ! gfortran overwrites the result of a nested call to a recursive
  ! function of deferred length within one expression
  recursive function infix(this, i, above) result(text)

    class(expression_view), intent(in) :: this
    integer               , intent(in) :: i, above
    character(len=:), allocatable :: text

    character(len=:), allocatable :: first, second, operator
    integer :: level

    operator = this % label(i)
    select case (this % v(i) % kind)
    case (VERTEX_LEAF, VERTEX_CONSTANT)
       level = 5
       text  = operator
    case (VERTEX_FUNCTION)
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
  ! THE BRACKETS: each vertex as [label operand operand], the input
  ! syntax of a forest figure.
  !===================================================================!

  function brackets(this) result(text)

    class(expression_view), intent(in) :: this
    character(len=:), allocatable :: text

    text = bracketed(this, size(this % v))

  end function brackets

  recursive function bracketed(this, i) result(text)

    class(expression_view), intent(in) :: this
    integer               , intent(in) :: i
    character(len=:), allocatable :: text

    character(len=:), allocatable :: operand

    text = '[' // this % label(i)
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
  ! THE DIGRAPH in the DOT language: one node per vertex, named by its
  ! index and labelled by its text, one arc from each operator to
  ! each operand.
  !===================================================================!

  function digraph(this) result(text)

    class(expression_view), intent(in) :: this
    character(len=:), allocatable :: text

    integer :: i

    text = 'digraph expression {' // new_line('a')
    do i = 1, size(this % v)
       text = text // '  v' // written(i) // ' [label="' // this % label(i) // '"];' // new_line('a')
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
