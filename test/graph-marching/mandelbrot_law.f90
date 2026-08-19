!=====================================================================!
! The mandelbrot law, on the tower: the escape-time fractal recast
! as a reacting source, ported from the old world's most elegant
! demo.
!
! The iteration z -> z^2 + c is forward euler with step one on
!
!      dz/dt = z^2 + c - z
!
! (write the step z + h*(z^2 + c - z), set h = 1: the z cancels and
! the map remains - an identity, not an approximation). Split
! z = (u, v) and c = (c1, c2): a two-variable reacting source whose
! coupling mixes the parts both ways. The marcher's convention
! moves the state AGAINST the statement, so the law answers minus
! the velocity,
!
!      S_u = u - (u*u - v*v) - c1
!      S_v = v - (2*u*v)     - c2
!
! and the marched map comes out z -> z^2 + c on the nose. Each cell
! carries its own c: the cell IS a point of the complex plane.
!=====================================================================!

module mandelbrot_law_fixture

  use iso_fortran_env    , only : dp => REAL64
  use graph_operation_view, only : graph_operation
  use graph_directed_view, only : directed_graph
  use graph_field_calculus, only : graph_field
  use graph_directed_view     , only : GRAPH_SIDE_VERTEX
  use fractal_graph         , only : set_graph => graph
  use class_graph_field  , only : field

  implicit none

  private
  public :: mandelbrot_law

  type, extends(graph_operation) :: mandelbrot_law

     real(dp), allocatable :: creal(:)
     real(dp), allocatable :: cimag(:)

   contains

     procedure :: name   => law_name
     procedure :: domain => law_domain
     procedure :: apply  => law_apply

  end type mandelbrot_law

contains

  pure function law_name(this) result(name)
    class(mandelbrot_law), intent(in) :: this
    character(len=:), allocatable :: name
    associate (u1 => this); end associate
    name = 'mandelbrot law'
  end function law_name

  subroutine law_domain(this, input_graph, domain, nentries)
    class(mandelbrot_law), intent(in)      :: this
    class(directed_graph), intent(in)               :: input_graph
    type(set_graph), intent(out) :: domain
    integer        , intent(out) :: nentries
    associate (u1 => this); end associate
    domain   = input_graph % all_vertices()
    nentries = input_graph % num_vertices()
  end subroutine law_domain

  subroutine law_apply(this, input_graph, input_data, output)

    class(mandelbrot_law), intent(in)              :: this
    class(directed_graph), intent(in)                       :: input_graph
    class(graph_field), intent(in), optional       :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output

    type(set_graph) :: cells
    type(field)   :: out
    real(dp), allocatable :: q(:), s(:)
    real(dp) :: u, v
    integer :: nv, k

    nv = input_graph % num_vertices()
    allocate(s(2 * nv))
    s = 0.0_dp

    if (present(input_data)) then
       call input_data(1) % get_real_vector(q)
       do k = 1, nv
          u = q(2 * k - 1)
          v = q(2 * k)
          s(2 * k - 1) = u - (u * u - v * v) - this % creal(k)
          s(2 * k)     = v - (2.0_dp * u * v) - this % cimag(k)
       end do
    end if

    cells = input_graph % vertex_set()
    out = field('velocity', cells, input_graph % num_vertices(), ncomp=2)
    call out % set_real_vector(s)

    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine law_apply

end module mandelbrot_law_fixture
