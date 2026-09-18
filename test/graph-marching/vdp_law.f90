!=====================================================================!
! Van der Pol, two ways: the law, stating its Jacobian once as
! partial_action, and its tangent as an AUGMENTED law marched
! forward on the same chain as an independent check. The reverse
! traversal derives the transposed linearization from the law
! itself, so agreement between the two passes is by construction.
!
!      du/dt = v                     J = | 0            1          |
!      dv/dt = mu (1-u^2) v - u         | -2 mu u v - 1  mu (1-u^2) |
!
! Every law returns MINUS its velocity, per the convention.
!=====================================================================!

module vdp_fixture

  use iso_fortran_env    , only : dp => REAL64
  use operation_action, only : operation, variation, binding, bound_real_vector, emit_real
  use view_directed, only : directed_graph
  use field_calculus, only : field
  use graph_fractal         , only : graph

  implicit none

  private
  public :: vdp_law, vdp_tangent_law

  type, extends(operation) :: vdp_law
     real(dp) :: mu = 1.0_dp
   contains
     procedure :: name           => law_name
     procedure :: domain         => law_domain
     procedure :: apply          => law_apply
     procedure :: max_degree     => law_max_degree
     procedure :: partial_action => law_partial_action
  end type vdp_law

  interface vdp_law
     module procedure create_law
  end interface vdp_law

  ! The augmented statement: (q, dq) marched together; the tangent
  ! is computed by the same forward traversal as the state.
  type, extends(operation) :: vdp_tangent_law
     real(dp) :: mu = 1.0_dp
   contains
     procedure :: name   => tangent_name
     procedure :: domain => law_domain2
     procedure :: apply  => tangent_apply
  end type vdp_tangent_law

  interface vdp_tangent_law
     module procedure create_tangent_law
  end interface vdp_tangent_law

contains

  ! The constructors declare the one argument, the state.
  function create_law(mu) result(this)
    real(dp), intent(in), optional :: mu
    type(vdp_law) :: this
    if (present(mu)) this % mu = mu
    call this % declare_arguments(1)
  end function create_law

  function create_tangent_law(mu) result(this)
    real(dp), intent(in), optional :: mu
    type(vdp_tangent_law) :: this
    if (present(mu)) this % mu = mu
    call this % declare_arguments(1)
  end function create_tangent_law

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

  subroutine law_domain(this, input_graph, domain, num_entries)
    class(vdp_law), intent(in) :: this
    class(directed_graph), intent(in) :: input_graph
    type(graph), intent(out) :: domain
    integer        , intent(out) :: num_entries
    associate (u1 => this); end associate
    domain   = input_graph % vertex_set()
    num_entries = input_graph % num_vertices()
  end subroutine law_domain

  subroutine law_domain2(this, input_graph, domain, num_entries)
    class(vdp_tangent_law), intent(in) :: this
    class(directed_graph), intent(in) :: input_graph
    type(graph), intent(out) :: domain
    integer        , intent(out) :: num_entries
    associate (u1 => this); end associate
    domain   = input_graph % vertex_set()
    num_entries = input_graph % num_vertices()
  end subroutine law_domain2

  subroutine law_apply(this, input_graph, inputs, output)

    class(vdp_law), intent(in)                     :: this
    class(directed_graph), intent(in)                       :: input_graph
    type(binding), intent(in), optional       :: inputs(:)
    class(field), allocatable, intent(inout) :: output

    real(dp), allocatable :: q(:)
    real(dp) :: s(2), u, v

    s = 0.0_dp
    if (present(inputs)) then
       call bound_real_vector(inputs, this % argument(1), q)
       u = q(1)
       v = q(2)
       s(1) = -v
       s(2) = -(this % mu * (1.0_dp - u * u) * v - u)
    end if

    call pack_velocity(input_graph, s, output)

  end subroutine law_apply

  subroutine tangent_apply(this, input_graph, inputs, output)

    class(vdp_tangent_law), intent(in)             :: this
    class(directed_graph), intent(in)                       :: input_graph
    type(binding), intent(in), optional       :: inputs(:)
    class(field), allocatable, intent(inout) :: output

    real(dp), allocatable :: q(:)
    real(dp) :: s(4), u, v, du, dv

    s = 0.0_dp
    if (present(inputs)) then
       call bound_real_vector(inputs, this % argument(1), q)
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

    call pack_velocity(input_graph, s, output)

  end subroutine tangent_apply

  pure function law_max_degree(this) result(degree)
    class(vdp_law), intent(in) :: this
    integer :: degree
    associate (u1 => this); end associate
    degree = 1
  end function law_max_degree

  ! Minus J dq at the state, the one argument: the velocity jacobian,
  ! written once; jacobian_of reads it and the reverse traversal
  ! transposes it.
  subroutine law_partial_action(this, input_graph, inputs, &
       & variations, output)

    class(vdp_law), intent(in)               :: this
    class(directed_graph), intent(in)        :: input_graph
    type(binding), intent(in)                :: inputs(:)
    type(variation), intent(in)              :: variations(:)
    class(field), allocatable, intent(inout) :: output

    real(dp), allocatable :: q(:), d(:)
    real(dp) :: s(2), u, v

    call this % require_owned(variations)
    if (size(variations) /= 1) error stop 'van der pol: the law is exact to first order'
    if (.not. variations(1) % argument_is(this % argument(1))) then
       error stop 'van der pol: the law takes one argument'
    end if

    call bound_real_vector(inputs, this % argument(1), q)
    call variations(1) % direction(d)
    u = q(1)
    v = q(2)
    s(1) = -d(2)
    s(2) = -((-2.0_dp * this % mu * u * v - 1.0_dp) * d(1) &
         &   + this % mu * (1.0_dp - u * u) * d(2))

    call pack_velocity(input_graph, s, output)

  end subroutine law_partial_action

  subroutine pack_velocity(input_graph, s, output)

    class(directed_graph), intent(in) :: input_graph
    real(dp), intent(in) :: s(:)
    class(field), allocatable, intent(inout) :: output

    call emit_real('velocity', input_graph % vertex_set(), input_graph % num_vertices(), &
         & s, output, num_components=size(s))

  end subroutine pack_velocity

end module vdp_fixture
