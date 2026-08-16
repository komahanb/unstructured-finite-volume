!=====================================================================!
! GRAPH VIEWS
!
! The kernel carries structure and identity. Attributes are bound
! here, in a partial map keyed by identity:
!
!     attribute_map : identity -> (number, symbol, index)
!
! VIEWS. Four maps over one graph type, introducing no graph subtype:
!
!     relation_view    G -> {(i,j)}             tuples, by index
!     residual         (G, Q) -> R(Q)           operator
!     compile_csr      G -> (rowptr, colidx)    compiled representation
!     evaluate         G -> number              expression
!
! ENCODING. A graph has exactly two branches, so a finite sequence is
! encoded by nesting:
!
!     sequence node    branch(1) = element      branch(2) = tail | NULL
!     pair node        branch(1) = left member  branch(2) = right member
!
! NULL terminates a sequence, and NULL denotes an absent member. Only
! NULL denotes absence: an UNKNOWN member is refused, because the
! kernel keeps NULL and UNKNOWN distinct and a view must not merge
! them.
!
! csr is a representation, not a subtype of graph.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module graph_views

  use graph_identity, only : token
  use fractal_graph , only : graph, GRAPH_NULL, GRAPH_UNKNOWN, GRAPH_KNOWN

  implicit none

  private
  public :: dp, attribute_map, csr
  public :: evaluate, num_atoms, relation_view, residual, compile_csr

  integer, parameter :: dp = kind(1.0d0)

  real(dp)        , parameter :: NO_NUMBER = -huge(1.0_dp)
  character(len=1), parameter :: NO_SYMBOL = ' '
  integer         , parameter :: NO_INDEX  = 0

  !===================================================================!
  ! One row per graph in the domain of the map. Absence of a column
  ! is absence from the domain: the map refuses rather than answers.
  !===================================================================!

  type :: attribute_map

     type(token)     , allocatable, private :: key(:)
     real(dp)        , allocatable, private :: number(:)
     character(len=1), allocatable, private :: symbol(:)
     integer         , allocatable, private :: index(:)

   contains

     procedure :: bind
     procedure :: number_of
     procedure :: symbol_of
     procedure :: index_of
     procedure, private :: row

  end type attribute_map

  !===================================================================!
  ! Compressed sparse row adjacency. No branch, no identity, no
  ! extension of graph.
  !===================================================================!

  type :: csr

     integer, allocatable :: rowptr(:)
     integer, allocatable :: colidx(:)

  end type csr

contains

  !===================================================================!
  ! Extend the map at one identity. Columns are independent.
  !===================================================================!

  subroutine bind(this, g, number, symbol, index)

    class(attribute_map), intent(inout)        :: this
    type(graph)         , intent(in)           :: g
    real(dp)            , intent(in), optional :: number
    character(len=1)    , intent(in), optional :: symbol
    integer             , intent(in), optional :: index

    integer :: k

    if (.not. allocated(this % key)) then
       allocate(this % key(0), this % number(0))
       allocate(this % symbol(0), this % index(0))
    end if

    k = this % row(g)
    if (k .eq. 0) then
       this % key    = [this % key   , g % id() ]
       this % number = [this % number, NO_NUMBER]
       this % symbol = [this % symbol, NO_SYMBOL]
       this % index  = [this % index , NO_INDEX ]
       k = size(this % key)
    end if

    if (present(number)) this % number(k) = number
    if (present(symbol)) this % symbol(k) = symbol
    if (present(index )) this % index(k)  = index

  end subroutine bind

  pure integer function row(this, g)

    class(attribute_map), intent(in) :: this
    type(graph)         , intent(in) :: g

    integer :: k

    row = 0
    if (.not. allocated(this % key)) return
    do k = 1, size(this % key)
       if (this % key(k) % matches(g % id())) row = k
    end do

  end function row

  real(dp) function number_of(this, g)

    class(attribute_map), intent(in) :: this
    type(graph)         , intent(in) :: g

    integer :: k

    k = this % row(g)
    if (k .eq. 0) error stop 'graph_views: no number bound to this identity'
    if (this % number(k) .eq. NO_NUMBER) then
       error stop 'graph_views: no number bound to this identity'
    end if
    number_of = this % number(k)

  end function number_of

  character(len=1) function symbol_of(this, g)

    class(attribute_map), intent(in) :: this
    type(graph)         , intent(in) :: g

    integer :: k

    k = this % row(g)
    if (k .eq. 0) error stop 'graph_views: no symbol bound to this identity'
    if (this % symbol(k) .eq. NO_SYMBOL) then
       error stop 'graph_views: no symbol bound to this identity'
    end if
    symbol_of = this % symbol(k)

  end function symbol_of

  integer function index_of(this, g)

    class(attribute_map), intent(in) :: this
    type(graph)         , intent(in) :: g

    integer :: k

    k = this % row(g)
    if (k .eq. 0) error stop 'graph_views: no index bound to this identity'
    if (this % index(k) .eq. NO_INDEX) then
       error stop 'graph_views: no index bound to this identity'
    end if
    index_of = this % index(k)

  end function index_of

  !===================================================================!
  ! EXPRESSION VIEW.
  !
  !     (NULL, NULL)     -> number_of(G)
  !     (KNOWN, KNOWN)   -> symbol_of(G) applied to both branches
  !     otherwise          refused
  !===================================================================!

  recursive real(dp) function evaluate(g, m) result(v)

    type(graph)        , intent(in) :: g
    type(attribute_map), intent(in) :: m

    character(len=1) :: op

    if (g % branch(1) % status() .eq. GRAPH_NULL .and. &
         & g % branch(2) % status() .eq. GRAPH_NULL) then
       v = m % number_of(g)
       return
    end if

    if (g % branch(1) % status() .ne. GRAPH_KNOWN .or. &
         & g % branch(2) % status() .ne. GRAPH_KNOWN) then
       error stop 'graph_views: UNKNOWN branch has no value'
    end if

    op = m % symbol_of(g)
    select case (op)
    case ('+')
       v = evaluate(g % branch(1) % known(), m) + evaluate(g % branch(2) % known(), m)
    case ('*')
       v = evaluate(g % branch(1) % known(), m) * evaluate(g % branch(2) % known(), m)
    case default
       error stop 'graph_views: no operator defined for this symbol'
    end select

  end function evaluate

  !===================================================================!
  ! Number of (NULL, NULL) graphs reachable through KNOWN branches.
  ! Defined on acyclic graphs.
  !===================================================================!

  recursive integer function num_atoms(g) result(n)

    type(graph), intent(in) :: g

    integer :: s

    if (g % branch(1) % status() .eq. GRAPH_NULL .and. &
         & g % branch(2) % status() .eq. GRAPH_NULL) then
       n = 1
       return
    end if

    n = 0
    do s = 1, 2
       if (g % branch(s) % status() .eq. GRAPH_KNOWN) then
          n = n + num_atoms(g % branch(s) % known())
       end if
    end do

  end function num_atoms

  !===================================================================!
  ! RELATION VIEW. The sequence of pairs, reported by index.
  !
  !     NULL member    -> 0, absence
  !     KNOWN member   -> index_of(member)
  !     UNKNOWN member -> refused, since only NULL denotes absence
  !===================================================================!

  subroutine relation_view(sequence, m, left, right)

    type(graph)        , target, intent(in)  :: sequence
    type(attribute_map)        , intent(in)  :: m
    integer, allocatable       , intent(out) :: left(:)
    integer, allocatable       , intent(out) :: right(:)

    type(graph), pointer :: node, pair

    allocate(left(0), right(0))

    node => sequence
    do
       if (node % branch(1) % status() .ne. GRAPH_KNOWN) then
          error stop 'graph_views: sequence node requires a KNOWN element'
       end if
       pair => node % branch(1) % known()
       left  = [left , member(pair, 1)]
       right = [right, member(pair, 2)]
       if (node % branch(2) % status() .eq. GRAPH_NULL) exit
       if (node % branch(2) % status() .ne. GRAPH_KNOWN) then
          error stop 'graph_views: sequence tail must be KNOWN or NULL'
       end if
       node => node % branch(2) % known()
    end do

  contains

    integer function member(p, s)
      type(graph), intent(in) :: p
      integer    , intent(in) :: s
      select case (p % branch(s) % status())
      case (GRAPH_NULL)
         member = 0
      case (GRAPH_KNOWN)
         member = m % index_of(p % branch(s) % known())
      case default
         error stop 'graph_views: UNKNOWN member is not NULL; only NULL denotes absence'
      end select
    end function member

  end subroutine relation_view

  !===================================================================!
  ! OPERATOR VIEW over the same sequence.
  !
  !     R_i(Q) = sum over pairs (i,j) of ( Q_j - Q_i )
  !
  ! A pair with a NULL member contributes nothing. Q constant implies
  ! R(Q) = 0, and sum(R) = 0 for every Q.
  !===================================================================!

  subroutine residual(sequence, m, q, res)

    type(graph)        , target, intent(in)  :: sequence
    type(attribute_map)        , intent(in)  :: m
    real(dp)                   , intent(in)  :: q(:)
    real(dp), allocatable      , intent(out) :: res(:)

    integer, allocatable :: left(:), right(:)
    integer              :: f, i, j

    call relation_view(sequence, m, left, right)

    allocate(res(size(q)))
    res = 0.0_dp

    do f = 1, size(left)
       i = left(f)
       j = right(f)
       if (i .eq. 0 .or. j .eq. 0) cycle
       res(i) = res(i) + (q(j) - q(i))
       res(j) = res(j) + (q(i) - q(j))
    end do

  end subroutine residual

  !===================================================================!
  ! COMPILED REPRESENTATION over the same sequence. The dependency
  ! runs graph_views -> fractal_graph and not back: the kernel does not
  ! name csr, and csr does not name graph.
  !===================================================================!

  type(csr) function compile_csr(sequence, m, n) result(this)

    type(graph)        , target, intent(in) :: sequence
    type(attribute_map)        , intent(in) :: m
    integer                    , intent(in) :: n

    integer, allocatable :: left(:), right(:), degree(:), fill(:)
    integer              :: f, i, j

    call relation_view(sequence, m, left, right)

    allocate(degree(n)); degree = 0
    do f = 1, size(left)
       i = left(f); j = right(f)
       if (i .eq. 0 .or. j .eq. 0) cycle
       degree(i) = degree(i) + 1
       degree(j) = degree(j) + 1
    end do

    allocate(this % rowptr(n + 1))
    this % rowptr(1) = 1
    do i = 1, n
       this % rowptr(i + 1) = this % rowptr(i) + degree(i)
    end do

    allocate(this % colidx(this % rowptr(n + 1) - 1))
    allocate(fill(n)); fill = this % rowptr(1:n)
    do f = 1, size(left)
       i = left(f); j = right(f)
       if (i .eq. 0 .or. j .eq. 0) cycle
       this % colidx(fill(i)) = j; fill(i) = fill(i) + 1
       this % colidx(fill(j)) = i; fill(j) = fill(j) + 1
    end do

  end function compile_csr

end module graph_views
