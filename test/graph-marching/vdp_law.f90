!=====================================================================!
! Van der Pol, three ways, all from one graph: the law, its tangent
! as an AUGMENTED law marched forward on the same chain, and its
! transposed linearization walked in reverse. The paper's witness,
! ported: agreement between the roads is by construction, not by
! adjustment.
!
!      du/dt = v                     J = | 0            1          |
!      dv/dt = mu (1-u^2) v - u         | -2 mu u v - 1  mu (1-u^2) |
!
! Every law answers MINUS its velocity, per the house convention.
!=====================================================================!

module vdp_fixture

  use iso_fortran_env    , only : dp => REAL64
  use operation_action, only : graph_operation
  use graph_directed_view, only : directed_graph
  use graph_field_calculus, only : graph_field
  use graph_directed_view     , only : GRAPH_SIDE_VERTEX
  use fractal_graph         , only : set_graph => graph
  use class_graph_field  , only : field

  implicit none

  private
  public :: vdp_law, vdp_tangent_law, vdp_adjoint_law

  type, extends(graph_operation) :: vdp_law
     real(dp) :: mu = 1.0_dp
   contains
     procedure :: name   => law_name
     procedure :: domain => law_domain
     procedure :: apply  => law_apply
  end type vdp_law

  ! The augmented statement: (q, dq) marched together; the tangent
  ! rides the same forward walk as the state.
  type, extends(graph_operation) :: vdp_tangent_law
     real(dp) :: mu = 1.0_dp
   contains
     procedure :: name   => tangent_name
     procedure :: domain => law_domain2
     procedure :: apply  => tangent_apply
  end type vdp_tangent_law

  ! The transposed linearization at a stored state, for the reverse
  ! walk; the test moves `at` along the stored trajectory.
  type, extends(graph_operation) :: vdp_adjoint_law
     real(dp) :: mu = 1.0_dp
     real(dp) :: at(2) = 0.0_dp
   contains
     procedure :: name   => adjoint_name
     procedure :: domain => law_domain3
     procedure :: apply  => adjoint_apply
  end type vdp_adjoint_law

contains

  pure function law_name(this) result(name)
    class(vdp_law), intent(in) :: this
    character(len=:), allocatable :: name
    associate (u1 => this); end associate
    name = 'van der pol'
  end function law_name

  pure function tangent_name(this) result(name)
    class(vdp_tangent_law), intent(in) :: this
    character(len=:), allocatable :: name
    associate (u1 => this); end associate
    name = 'van der pol, augmented'
  end function tangent_name

  pure function adjoint_name(this) result(name)
    class(vdp_adjoint_law), intent(in) :: this
    character(len=:), allocatable :: name
    associate (u1 => this); end associate
    name = 'van der pol, transposed'
  end function adjoint_name

  subroutine law_domain(this, input_graph, domain, nentries)
    class(vdp_law), intent(in) :: this
    class(directed_graph), intent(in) :: input_graph
    type(set_graph), intent(out) :: domain
    integer        , intent(out) :: nentries
    associate (u1 => this); end associate
    domain   = input_graph % all_vertices()
    nentries = input_graph % num_vertices()
  end subroutine law_domain

  subroutine law_domain2(this, input_graph, domain, nentries)
    class(vdp_tangent_law), intent(in) :: this
    class(directed_graph), intent(in) :: input_graph
    type(set_graph), intent(out) :: domain
    integer        , intent(out) :: nentries
    associate (u1 => this); end associate
    domain   = input_graph % all_vertices()
    nentries = input_graph % num_vertices()
  end subroutine law_domain2

  subroutine law_domain3(this, input_graph, domain, nentries)
    class(vdp_adjoint_law), intent(in) :: this
    class(directed_graph), intent(in) :: input_graph
    type(set_graph), intent(out) :: domain
    integer        , intent(out) :: nentries
    associate (u1 => this); end associate
    domain   = input_graph % all_vertices()
    nentries = input_graph % num_vertices()
  end subroutine law_domain3

  subroutine law_apply(this, input_graph, input_data, output)

    class(vdp_law), intent(in)                     :: this
    class(directed_graph), intent(in)                       :: input_graph
    class(graph_field), intent(in), optional       :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output

    real(dp), allocatable :: q(:)
    real(dp) :: s(2), u, v

    s = 0.0_dp
    if (present(input_data)) then
       call input_data(1) % get_real_vector(q)
       u = q(1)
       v = q(2)
       s(1) = -v
       s(2) = -(this % mu * (1.0_dp - u * u) * v - u)
    end if

    call pack_answer(input_graph, s, output)

  end subroutine law_apply

  subroutine tangent_apply(this, input_graph, input_data, output)

    class(vdp_tangent_law), intent(in)             :: this
    class(directed_graph), intent(in)                       :: input_graph
    class(graph_field), intent(in), optional       :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output

    real(dp), allocatable :: q(:)
    real(dp) :: s(4), u, v, du, dv

    s = 0.0_dp
    if (present(input_data)) then
       call input_data(1) % get_real_vector(q)
       u  = q(1)
       v  = q(2)
       du = q(3)
       dv = q(4)
       s(1) = -v
       s(2) = -(this % mu * (1.0_dp - u * u) * v - u)
       s(3) = -dv
       s(4) = -((-2.0_dp * this % mu * u * v - 1.0_dp) * du &
            &   + this % mu * (1.0_dp - u * u) * dv)
    end if

    call pack_answer(input_graph, s, output)

  end subroutine tangent_apply

  subroutine adjoint_apply(this, input_graph, input_data, output)

    class(vdp_adjoint_law), intent(in)             :: this
    class(directed_graph), intent(in)                       :: input_graph
    class(graph_field), intent(in), optional       :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output

    real(dp), allocatable :: lam(:)
    real(dp) :: s(2), u, v

    s = 0.0_dp
    if (present(input_data)) then
       call input_data(1) % get_real_vector(lam)
       u = this % at(1)
       v = this % at(2)
       ! Minus J^T lambda: the transpose of the velocity jacobian.
       s(1) = -((-2.0_dp * this % mu * u * v - 1.0_dp) * lam(2))
       s(2) = -(lam(1) + this % mu * (1.0_dp - u * u) * lam(2))
    end if

    call pack_answer(input_graph, s, output)

  end subroutine adjoint_apply

  subroutine pack_answer(input_graph, s, output)

    class(directed_graph), intent(in) :: input_graph
    real(dp), intent(in) :: s(:)
    class(graph_field), allocatable, intent(inout) :: output

    type(set_graph) :: cells
    type(field)   :: out
    integer :: v

    cells = input_graph % vertex_set()
    out = field('velocity', cells, input_graph % num_vertices(), ncomp=size(s))
    call out % set_real_vector(s)

    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine pack_answer

end module vdp_fixture
