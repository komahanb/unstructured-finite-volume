!=====================================================================!
! THE TEMPORAL STEP FIXTURE - earned at LEVEL 6.
!
! One step of a multistep discretisation of qdot = -S(q), as the
! residual
!
!      R(q_n, q_(n-1), ..., q_(n-r)) = sum_(j=0..r) c_j q_(n-j)
!                                    + h [ theta S(q_n)
!                                        + (1 - theta) S(q_(n-1)) ]
!
! The coefficients c_j are READ FROM PRODUCTION: operation_family's
! BDF family evaluates the slope at zero of the Lagrange basis
! through the r + 1 instants behind n, on the uniform step h. This
! module writes no coefficient of its own. theta = 1 is backward
! euler (r = 1) and bdf-k (r = k); theta = 0 is forward euler.
!
! operation_step, which stated this residual in src, was deleted in
! b9944c3 as unreachable; the family that supplied its coefficients
! is what remains, and the residual is restated here so the levels
! above can check that a temporal discretisation preserves the
! DOMAIN of the action it discretises.
!
! Argument 1 is the state q_n; argument 1 + j is the history q_(n-j)
! for j = 1..r. The step stores no state: its domain is its action's
! domain, whatever host it is applied on.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module temporal_step_fixture

  use iso_fortran_env      , only : dp => REAL64
  use operation_action     , only : operation, binding, contract
  use operation_action     , only : bound_value, emit_real
  use view_directed        , only : directed_graph
  use field_calculus       , only : field, FIELD_REAL
  use graph_fractal        , only : graph
  use operation_family     , only : family, bdf_family
  use util_derivative_terms, only : derivative_terms, value

  implicit none

  private
  public :: temporal_step
  public :: backward_euler, forward_euler, bdf

  !===================================================================!
  ! The step: the action it discretises, the uniform step h, theta,
  ! the reach r and the r + 1 coefficients c(0:r).
  !===================================================================!

  type, extends(operation) :: temporal_step

     class(operation), allocatable :: action
     real(dp) :: h     = 0.0_dp
     real(dp) :: theta = 1.0_dp
     integer  :: reach = 0
     real(dp), allocatable :: c(:)

   contains

     procedure :: name   => step_name
     procedure :: domain => step_domain
     procedure :: apply  => step_apply

  end type temporal_step

contains

  !===================================================================!
  ! The three constructors: one family, theta and the order.
  !===================================================================!

  function backward_euler(action, h) result(this)

    class(operation), intent(in) :: action
    real(dp)        , intent(in) :: h
    type(temporal_step) :: this

    this = bdf(1, action, h)

  end function backward_euler

  function forward_euler(action, h) result(this)

    class(operation), intent(in) :: action
    real(dp)        , intent(in) :: h
    type(temporal_step) :: this

    this = bdf(1, action, h)
    this % theta = 0.0_dp

  end function forward_euler

  !===================================================================!
  ! bdf-k on the uniform step h: the reach and every coefficient
  ! read from the production family. A step that is not positive
  ! stops the program.
  !===================================================================!

  function bdf(k, action, h) result(this)

    integer         , intent(in) :: k
    class(operation), intent(in) :: action
    real(dp)        , intent(in) :: h
    type(temporal_step) :: this

    type(family) :: scheme
    type(derivative_terms), allocatable :: dt(:)
    integer :: j, head

    if (h <= 0.0_dp) error stop 'temporal_step: the step is positive'

    scheme = bdf_family(k)

    allocate(this % action, source=action)
    this % h     = h
    this % theta = 1.0_dp
    this % reach = scheme % history_depth(1)

    ! the r + 1 instants behind n on the uniform step, the row
    ! determining the derivative (degree 1) from the values (degree 0)
    head = this % reach + 1
    allocate(dt(head))
    do j = 1, head
       dt(j) = derivative_terms(h, 0)
    end do
    allocate(this % c(0:this % reach))
    do j = 0, this % reach
       this % c(j) = value(scheme % edge_coefficient(dt, head - j, head, 0, 1))
    end do

    call this % declare_arguments(1 + this % reach, &
         & [(contract(FIELD_REAL, 1), j = 0, this % reach)])

  end function bdf

  pure function step_name(this) result(name)

    class(temporal_step), intent(in) :: this
    character(len=:), allocatable    :: name

    associate (u1 => this); end associate

    name = 'temporal step'

  end function step_name

  !===================================================================!
  ! The residual is defined where the action is: on the action's
  ! own domain, never on the host's vertices.
  !===================================================================!

  subroutine step_domain(this, input_graph, domain, num_entries)

    class(temporal_step) , intent(in)  :: this
    class(directed_graph), intent(in)  :: input_graph
    type(graph)          , intent(out) :: domain
    integer              , intent(out) :: num_entries

    call this % action % domain(input_graph, domain, num_entries)

  end subroutine step_domain

  !===================================================================!
  ! The position of an argument's binding; an unbound argument stops
  ! the program.
  !===================================================================!

  integer function position_of(this, inputs, k) result(at)

    class(temporal_step), intent(in) :: this
    type(binding)       , intent(in) :: inputs(:)
    integer             , intent(in) :: k

    integer :: i

    at = 0
    do i = 1, size(inputs)
       if (inputs(i) % argument_is(this % argument(k))) at = i
    end do
    if (at == 0) error stop 'temporal_step: every instant the scheme reaches is bound'

  end function position_of

  !===================================================================!
  ! The values bound to argument k, required on the action's domain.
  !===================================================================!

  subroutine values_at(this, inputs, k, expected, q)

    class(temporal_step), intent(in)  :: this
    type(binding)       , intent(in)  :: inputs(:)
    integer             , intent(in)  :: k
    type(graph)         , intent(in)  :: expected
    real(dp), allocatable, intent(out) :: q(:)

    class(field), allocatable :: given

    call bound_value(inputs, this % argument(k), given)
    if (.not. given % defined_on(expected)) then
       error stop 'temporal_step: every state is defined on the action''s own domain'
    end if
    call given % real_vector(q)

  end subroutine values_at

  !===================================================================!
  ! S at the state bound to argument k, required on the action's
  ! domain.
  !===================================================================!

  subroutine action_at(this, input_graph, inputs, k, expected, s)

    class(temporal_step) , intent(in)  :: this
    class(directed_graph), intent(in)  :: input_graph
    type(binding)        , intent(in)  :: inputs(:)
    integer              , intent(in)  :: k
    type(graph)          , intent(in)  :: expected
    real(dp), allocatable, intent(out) :: s(:)

    class(field), allocatable :: velocity
    integer :: at

    at = position_of(this, inputs, k)
    call this % action % apply(input_graph, this % action % bind(inputs(at:at)), velocity)
    if (.not. velocity % defined_on(expected)) then
       error stop 'temporal_step: the action result is defined on its stated domain'
    end if
    call velocity % real_vector(s)

  end subroutine action_at

  !===================================================================!
  ! The residual. Without inputs it is zero on the action's domain.
  !===================================================================!

  subroutine step_apply(this, input_graph, inputs, output)

    class(temporal_step)     , intent(in)    :: this
    class(directed_graph)    , intent(in)    :: input_graph
    type(binding)            , intent(in), optional :: inputs(:)
    class(field), allocatable, intent(inout) :: output

    type(graph) :: expected
    real(dp), allocatable :: q(:), s(:), y(:)
    integer :: n_expected, j

    call this % action % domain(input_graph, expected, n_expected)

    if (present(inputs)) then

       call values_at(this, inputs, 1, expected, q)
       y = this % c(0) * q
       do j = 1, this % reach
          call values_at(this, inputs, 1 + j, expected, q)
          y = y + this % c(j) * q
       end do

       if (this % theta /= 0.0_dp) then
          call action_at(this, input_graph, inputs, 1, expected, s)
          y = y + this % h * this % theta * s
       end if
       if (this % theta /= 1.0_dp) then
          call action_at(this, input_graph, inputs, 2, expected, s)
          y = y + this % h * (1.0_dp - this % theta) * s
       end if

    else
       allocate(y(n_expected))
       y = 0.0_dp
    end if

    call emit_real('step residual', expected, n_expected, y, output, num_components=1)

  end subroutine step_apply

end module temporal_step_fixture
