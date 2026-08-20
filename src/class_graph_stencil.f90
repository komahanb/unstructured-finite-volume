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
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_graph_stencil

  use iso_fortran_env    , only : dp => REAL64
  use graph_directed_view, only : directed_graph
  use graph_field_calculus, only : graph_field
  use graph_discretization     , only : discretization_operator
  use relation_binary, only : group_by_key
  use class_graph_field  , only : field
  use class_graph        , only : directed_stored_graph
  use graph_fractal      , only : set_graph => graph

  implicit none

  private
  public :: stencil_operator
  public :: combine_triples

  type, extends(discretization_operator) :: stencil_operator

     type(directed_stored_graph) :: pattern

     type(field) :: weights
     type(field) :: constants

     character(len=:), allocatable :: label

   contains

     procedure :: name         => stencil_name
     procedure :: domain       => stencil_domain
     procedure :: apply        => stencil_apply
     procedure :: dependencies => stencil_dependencies

  end type stencil_operator

  interface stencil_operator
     module procedure create
     module procedure create_dense
  end interface stencil_operator

contains

  !===================================================================!
  ! Build from the triples and the constants. The triples become the
  ! dependency graph - one edge per (row, column) pair, tail at the
  ! column, head at the row - and the numbers become its fields.
  !===================================================================!

  type(stencil_operator) function create(rows, columns, weights, &
       & constant, label) result(this)

    integer , intent(in) :: rows(:)
    integer , intent(in) :: columns(:)
    real(dp), intent(in) :: weights(:)
    real(dp), intent(in) :: constant(:)
    character(len=*), intent(in), optional :: label

    integer :: nv

    nv = size(constant)

    this % pattern = directed_stored_graph(nv, tails=columns, heads=rows)

    this % weights   = field('stencil weights', this % pattern % edge_set(), this % pattern % num_edges())
    call this % weights % set_real_vector(weights)
    this % constants = field('stencil constants', this % pattern % vertex_set(), this % pattern % num_vertices())
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

  type(stencil_operator) function create_dense(a, label) result(this)

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

  pure function stencil_name(this) result(name)

    class(stencil_operator), intent(in) :: this
    character(len=:), allocatable :: name

    name = this % label

  end function stencil_name

  subroutine stencil_domain(this, input_graph, domain, nentries)

    class(stencil_operator), intent(in)    :: this
    class(directed_graph), intent(in)               :: input_graph
    type(set_graph), intent(out) :: domain
    integer        , intent(out) :: nentries

    associate (u1 => this); end associate

    domain   = input_graph % all_vertices()
    nentries = input_graph % num_vertices()

  end subroutine stencil_domain

  !===================================================================!
  ! y = constants + the dependency edges, walked once: each edge
  ! carries its weight times the tail's value onto its head.
  !===================================================================!

  subroutine stencil_apply(this, input_graph, input_data, output)

    class(stencil_operator), intent(in)            :: this
    class(directed_graph), intent(in)                       :: input_graph
    class(graph_field), intent(in), optional       :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output

    type(field)   :: out
    real(dp), allocatable :: q(:), y(:), w(:)
    integer :: nv, e

    nv = input_graph % num_vertices()

    call this % constants % get_real_vector(y)

    if (present(input_data)) then
       call input_data(1) % get_real_vector(q)
       call this % weights % get_real_vector(w)
       do e = 1, this % pattern % num_edges()
          y(this % pattern % edge_head(e)) = &
               & y(this % pattern % edge_head(e)) &
               & + w(e) * q(this % pattern % edge_tail(e))
       end do
    end if

    out = field(this % label, input_graph % vertex_set(), input_graph % num_vertices())
    call out % set_real_vector(y)

    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine stencil_apply

  !===================================================================!
  ! The contract's answer: the pattern IS a graph, handed out whole.
  !===================================================================!

  subroutine stencil_dependencies(this, pattern)

    class(stencil_operator), intent(in)    :: this
    class(directed_graph), allocatable, intent(out) :: pattern

    allocate(pattern, source=this % pattern)

  end subroutine stencil_dependencies

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

end module class_graph_stencil
