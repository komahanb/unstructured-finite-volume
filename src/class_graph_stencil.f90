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
  use graph_ordinary_view, only : graph
  use graph_field_calculus, only : graph_field
  use graph_calculus     , only : discretization_operator
  use class_graph_field  , only : field
  use class_graph        , only : stored_graph
  use fractal_graph      , only : set_graph => graph

  implicit none

  private
  public :: stencil_operator

  type, extends(discretization_operator) :: stencil_operator

     type(stored_graph) :: pattern

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

    this % pattern = stored_graph(nv, tails=columns, heads=rows)

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

  pure function stencil_name(this) result(name)

    class(stencil_operator), intent(in) :: this
    character(len=:), allocatable :: name

    name = this % label

  end function stencil_name

  subroutine stencil_domain(this, input_graph, domain, nentries)

    class(stencil_operator), intent(in)    :: this
    class(graph), intent(in)               :: input_graph
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
    class(graph), intent(in)                       :: input_graph
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
    class(graph), allocatable, intent(out) :: pattern

    allocate(pattern, source=this % pattern)

  end subroutine stencil_dependencies

end module class_graph_stencil
