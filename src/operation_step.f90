!=====================================================================!
! The time discretization: the residual of one step of the theta
! family over an action S, with the history of up to two instants,
!
!      R(q_n, xi, q_(n-1), q_(n-2)) = a0 q_n + a1 q_(n-1) + a2 q_(n-2)
!                                   + h [ theta S(q_n, xi)
!                                       + (1 - theta) S(q_(n-1), xi) ]
!
! theta = 1 is backward euler (reach 1) and BDF2 (reach 2); theta = 0
! is forward euler; theta = 1/2 is Crank-Nicolson. Constructors:
!
!      backward_euler(action, h)          reach 1, [1, -1]
!      bdf(k, action, h)                  reach k, uniform
!      bdf_variable(k, action, steps)     reach k, coefficients
!                                         from the steps taken
!
! The action returns minus the velocity, so the residual is zero at
! the discrete solution.
!
! ARGUMENTS. The scheme's argument space is the action's followed by
! its history: arguments 1..m stand for the action's arguments 1..m
! (the state, then any auxiliaries such as a parameter xi); argument
! m+k is history(k), the state k instants back, k <= reach. The
! history may arrive as inputs m+1.. or, for now, from the stored
! qold/qolder the marcher assigns; the inputs win when given.
!
! DERIVATIVES. Every partial of R is derived from the action's own
! partial actions by the formula above: the state block is
! a0 I + h theta D_q S(q_n), the first history block is
! a1 I + h (1 - theta) D_q S(q_(n-1)), the second is a2 I, and an
! auxiliary block draws from both states. A variation on a scheme
! argument is restated on the action's argument before the action is
! asked; the action never sees a scheme argument.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_step

  use iso_fortran_env    , only : dp => REAL64
  use operation_action, only : operation, argument, variation
  use operation_discretization     , only : discretization
  use view_directed, only : directed_graph
  use field_calculus, only : field
  use graph_fractal      , only : graph
  use field_stored  , only : stored_field
  use view_directed_stored        , only : stored_directed_graph

  implicit none

  private
  public :: scheme
  public :: backward_euler, bdf, bdf_variable

  type, extends(discretization) :: scheme

     class(operation), allocatable :: action

     real(dp), allocatable :: qold(:)
     real(dp), allocatable :: qolder(:)

     real(dp) :: a0 = 1.0_dp
     real(dp) :: a1 = -1.0_dp
     real(dp) :: a2 = 0.0_dp
     real(dp) :: hs = 0.0_dp
     real(dp) :: theta = 1.0_dp

     integer :: reach = 1

   contains

     procedure :: name         => step_name
     procedure :: domain       => step_domain
     procedure :: apply        => step_apply
     procedure :: dependencies => step_dependencies
     procedure :: set_bdf

     procedure :: max_degree     => step_max_degree
     procedure :: partial_action => step_partial_action

     procedure :: state
     procedure :: auxiliary
     procedure :: history
     procedure :: action_argument
     procedure :: from_action

  end type scheme

contains

  !===================================================================!
  ! Backward euler: one instant back.
  !===================================================================!

  function backward_euler(action, h) result(this)

    class(operation), intent(in) :: action
    real(dp), intent(in)               :: h
    type(scheme)                :: this

    allocate(this % action, source=action)
    call this % set_bdf(1, [h])

  end function backward_euler

  !===================================================================!
  ! BDF of order k at a uniform step. Only orders 1 and 2 have
  ! coefficient tables; any other order stops the program.
  !===================================================================!

  function bdf(k, action, h) result(this)

    integer, intent(in)                :: k
    class(operation), intent(in) :: action
    real(dp), intent(in)               :: h
    type(scheme)                :: this

    allocate(this % action, source=action)

    select case (k)
    case (1)
       call this % set_bdf(1, [h])
    case (2)
       call this % set_bdf(2, [h, h])
    case default
       error stop 'bdf: only orders one and two carry tables so far'
    end select

  end function bdf

  !===================================================================!
  ! BDF of order k with the steps actually taken: steps(1) is the
  ! current step h_n, steps(2) the one before it. The order-2
  ! coefficients are the derivative of the interpolating quadratic
  ! at the newest instant, exact on quadratics for any spacing.
  !===================================================================!

  function bdf_variable(k, action, steps) result(this)

    integer, intent(in)                :: k
    class(operation), intent(in) :: action
    real(dp), intent(in)               :: steps(:)
    type(scheme)                :: this

    allocate(this % action, source=action)
    call this % set_bdf(k, steps)

  end function bdf_variable

  !===================================================================!
  ! Reconfigure the operator to order k with the given steps,
  ! scaled so the residual reads a0 q + a1 qold + a2 qolder
  ! + h_n S(q). With h0 = h_n and h1 = h_{n-1} the order-2 row is
  !
  !      a0 = (2 h0 + h1)/(h0 + h1)
  !      a1 = -(h0 + h1)/h1
  !      a2 = h0^2 / (h1 (h0 + h1))
  !
  ! assigned as the exact constants (3/2, -2, 1/2) when h0 = h1 so
  ! no rounding enters the uniform case. Checks, each stopping the
  ! program: only orders 1 and 2 have tables; size(steps) must
  ! equal k; every step must be positive. The argument space is
  ! declared here, once the action is known: the action's arguments
  ! then k history arguments; a later call changes the count only.
  !===================================================================!

  subroutine set_bdf(this, k, steps)

    class(scheme), intent(inout) :: this
    integer             , intent(in)    :: k
    real(dp)            , intent(in)    :: steps(:)

    real(dp) :: h0, h1
    integer  :: j

    if (k /= 1 .and. k /= 2) then
       error stop 'step: only orders one and two carry tables so far'
    end if

    if (size(steps) /= k) then
       error stop 'step: the step count matches the order'
    end if

    do j = 1, k
       if (steps(j) <= 0.0_dp) then
          error stop 'step: a time step is positive'
       end if
    end do

    this % hs    = steps(1)
    this % reach = k

    if (allocated(this % action)) then
       if (this % action % num_arguments() < 1) then
          error stop 'step: the action declares its state argument'
       end if
       call this % declare_arguments(this % action % num_arguments() + k)
    end if

    if (k == 1) then
       this % a0 = 1.0_dp
       this % a1 = -1.0_dp
       this % a2 = 0.0_dp
       return
    end if

    h0 = steps(1)
    h1 = steps(2)

    if (h0 == h1) then
       this % a0 = 1.5_dp
       this % a1 = -2.0_dp
       this % a2 = 0.5_dp
    else
       this % a0 = (2.0_dp * h0 + h1) / (h0 + h1)
       this % a1 = -(h0 + h1) / h1
       this % a2 = (h0 * h0) / (h1 * (h0 + h1))
    end if

  end subroutine set_bdf

  !===================================================================!
  ! The arguments, by name. state is the action's state; auxiliary(j)
  ! the action's argument 1+j; history(k) the state k instants back.
  ! A history index past the reach, or an auxiliary the action does
  ! not declare, stops the program.
  !===================================================================!

  function state(this) result(a)

    class(scheme), intent(in) :: this
    type(argument) :: a

    a = this % argument(1)

  end function state

  function auxiliary(this, j) result(a)

    class(scheme), intent(in) :: this
    integer      , intent(in) :: j
    type(argument) :: a

    if (j < 1 .or. 1 + j > this % action % num_arguments()) then
       error stop 'step: the auxiliary argument is one the action declares'
    end if

    a = this % argument(1 + j)

  end function auxiliary

  function history(this, k) result(a)

    class(scheme), intent(in) :: this
    integer      , intent(in) :: k
    type(argument) :: a

    if (k < 1 .or. k > this % reach) then
       error stop 'step: the history argument is within the reach'
    end if

    a = this % argument(this % action % num_arguments() + k)

  end function history

  !===================================================================!
  ! The action's argument a scheme argument stands for: argument k
  ! of the scheme is argument k of the action for k <= m, and
  ! history(1) stands for the action's state, evaluated one instant
  ! back. history(2) stands for no action argument and an unnamed
  ! argument is returned. An argument of another operation stops the
  ! program.
  !===================================================================!

  function action_argument(this, a) result(b)

    class(scheme)  , intent(in) :: this
    type(argument) , intent(in) :: a
    type(argument) :: b

    type(argument) :: unnamed
    integer :: k, m

    if (.not. this % owns(a)) then
       error stop 'step: the argument is one of the scheme''s'
    end if

    m = this % action % num_arguments()

    do k = 1, m
       if (a % matches(this % argument(k))) then
          b = this % action % argument(k)
          return
       end if
    end do

    if (a % matches(this % history(1))) then
       b = this % action % argument(1)
       return
    end if

    ! history(2): no action argument stands behind it
    b = unnamed

  end function action_argument

  !===================================================================!
  ! The scheme argument standing for an action argument: the same
  ! position. An argument the action does not own stops the program.
  !===================================================================!

  function from_action(this, a) result(b)

    class(scheme)  , intent(in) :: this
    type(argument) , intent(in) :: a
    type(argument) :: b

    integer :: k

    do k = 1, this % action % num_arguments()
       if (a % matches(this % action % argument(k))) then
          b = this % argument(k)
          return
       end if
    end do

    error stop 'step: the argument is one of the action''s'

  end function from_action

  !===================================================================!
  ! The dependency pattern on the time axis. The residual at the
  ! newest instant reads every instant its coefficients reach,
  !
  !      R_n = a0 q_n + a1 q_(n-1) + a2 q_(n-2) + h S(q_n),
  !
  ! so the pattern is a fan-in of reach + 1 vertices onto the
  ! last, including the self-edge for the implicit unknown:
  !
  !      backward euler        bdf-2
  !          1 --> 2               1 --> 3
  !          2 --> 2               2 --> 3
  !                                3 --> 3
  !
  ! This records which instants the residual reads, not which
  ! instant follows which. A stencil's dependencies are
  ! the same pattern on the dependent-variable axis.
  !===================================================================!

  subroutine step_dependencies(this, pattern)

    class(scheme), intent(in)       :: this
    class(directed_graph), allocatable, intent(out) :: pattern

    integer :: n, newest

    newest = this % reach + 1

    allocate(pattern, source=stored_directed_graph(newest, &
         & tails=[(n, n = 1, newest)], &
         & heads=[(newest, n = 1, newest)]))

  end subroutine step_dependencies

  !===================================================================!
  ! The operation interface.
  !===================================================================!

  pure function step_name(this) result(name)

    class(scheme), intent(in) :: this
    character(len=:), allocatable :: name

    associate (u1 => this); end associate

    name = 'step'

  end function step_name

  !===================================================================!
  ! The step's domain is the domain of the action it discretizes:
  ! the residual is a statement about the same unknown S is a
  ! statement about, so the question is delegated to the action.
  !===================================================================!

  subroutine step_domain(this, input_graph, domain, num_entries)

    class(scheme), intent(in)       :: this
    class(directed_graph), intent(in)               :: input_graph
    type(graph), intent(out) :: domain
    integer        , intent(out) :: num_entries

    call this % action % domain(input_graph, domain, num_entries)

  end subroutine step_domain

  !===================================================================!
  ! The history state k instants back: input m+k when the caller
  ! supplied it (checked to live on the action's domain), else the
  ! stored qold/qolder. Neither given stops the program.
  !===================================================================!

  subroutine history_values(this, k, input_data, expected, values)

    class(scheme), intent(in)          :: this
    integer      , intent(in)          :: k
    class(field) , intent(in)          :: input_data(:)
    type(graph)  , intent(in)          :: expected
    real(dp), allocatable, intent(out) :: values(:)

    type(graph) :: given
    integer     :: at

    at = this % action % num_arguments() + k

    if (size(input_data) >= at) then
       given = input_data(at) % domain()
       if (.not. given % same_as(expected)) then
          error stop 'step: a history state lives on the action''s own domain'
       end if
       call input_data(at) % real_vector(values)
       return
    end if

    if (k == 1 .and. allocated(this % qold)) then
       values = this % qold
    else if (k == 2 .and. allocated(this % qolder)) then
       values = this % qolder
    else
       error stop 'step: the history state is given'
    end if

  end subroutine history_values

  !===================================================================!
  ! The action's input list at a state: the state as input 1 on the
  ! action's domain, then the auxiliaries the caller supplied in
  ! inputs 2..m, copied as stored fields. Fewer auxiliaries than the
  ! action declares are passed as fewer, as a direct call would.
  !===================================================================!

  subroutine inputs_at(this, input_data, expected, n_expected, &
       & num_components, values, inputs)

    class(scheme), intent(in) :: this
    class(field) , intent(in) :: input_data(:)
    type(graph)  , intent(in) :: expected
    integer      , intent(in) :: n_expected, num_components
    real(dp)     , intent(in) :: values(:)
    type(stored_field), allocatable, intent(out) :: inputs(:)

    real(dp), allocatable :: v(:)
    integer :: j, m

    m = min(this % action % num_arguments(), size(input_data))

    allocate(inputs(max(m, 1)))

    inputs(1) = stored_field('state', expected, n_expected, num_components=num_components)
    call inputs(1) % set_real_vector(values)

    do j = 2, m
       inputs(j) = stored_field(input_data(j) % name(), input_data(j) % domain(), &
            & input_data(j) % num_entries(), num_components=input_data(j) % num_components())
       call input_data(j) % real_vector(v)
       call inputs(j) % set_real_vector(v)
    end do

  end subroutine inputs_at

  !===================================================================!
  ! Evaluate the residual on the action's domain. Two identity
  ! checks stop the program: the state must live on the action's
  ! domain, and the action's result must land on that same domain.
  ! Equal length alone would let a field from a different domain
  ! through and produce plausible but wrong numbers. The action is
  ! evaluated at the newest state, and - when theta < 1 - at the
  ! state one instant back as well.
  !===================================================================!

  subroutine step_apply(this, input_graph, input_data, output)

    class(scheme), intent(in)               :: this
    class(directed_graph), intent(in)                       :: input_graph
    class(field), intent(in), optional       :: input_data(:)
    class(field), allocatable, intent(inout) :: output

    type(stored_field)   :: out
    class(field), allocatable :: velocity
    type(stored_field), allocatable :: inputs(:)
    type(graph) :: expected, given
    integer         :: n_expected, m
    real(dp), allocatable :: q(:), s(:), y(:), qold(:), qolder(:)
    integer :: num_components

    call this % action % domain(input_graph, expected, n_expected)

    if (present(input_data)) then

       given   = input_data(1) % domain()
       if (.not. given % same_as(expected)) then
          error stop 'step: the state must live on the action''s own domain'
       end if
       num_components = input_data(1) % num_components()

       call input_data(1) % real_vector(q)
       call history_values(this, 1, input_data, expected, qold)

       y = this % a0 * q + this % a1 * qold

       if (abs(this % a2) > 0.0_dp) then
          call history_values(this, 2, input_data, expected, qolder)
          y = y + this % a2 * qolder
       end if

       m = min(this % action % num_arguments(), size(input_data))

       if (this % theta /= 0.0_dp) then
          call this % action % apply(input_graph, input_data(1:m), velocity)
          given = velocity % domain()
          if (.not. given % same_as(expected)) then
             error stop 'step: the action result lives on its stated domain'
          end if
          call velocity % real_vector(s)
          y = y + this % hs * this % theta * s
       end if

       if (this % theta /= 1.0_dp) then
          call inputs_at(this, input_data, expected, n_expected, num_components, qold, inputs)
          call this % action % apply(input_graph, inputs, velocity)
          given = velocity % domain()
          if (.not. given % same_as(expected)) then
             error stop 'step: the action result lives on its stated domain'
          end if
          call velocity % real_vector(s)
          y = y + this % hs * (1.0_dp - this % theta) * s
       end if

    else
       num_components = 1
       allocate(y(n_expected))
       y = 0.0_dp
    end if

    out = stored_field('step residual', expected, n_expected, num_components=num_components)
    call out % set_real_vector(y)

    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine step_apply

  !===================================================================!
  ! The calculus reaches as deep as the action's. A scheme without
  ! an action stops the program, because there is nothing to
  ! differentiate.
  !===================================================================!

  pure function step_max_degree(this) result(degree)

    class(scheme), intent(in) :: this
    integer :: degree

    if (.not. allocated(this % action)) then
       error stop 'step: the action is attached'
    end if

    degree = this % action % max_degree()

  end function step_max_degree

  !===================================================================!
  ! One mixed partial of the residual, derived from the action's.
  ! With V the variation list, k its length, and nq, nh1, nh2 the
  ! counts of factors on the state, history(1) and history(2):
  !
  !      nh2 > 0              a2 v when k = 1, else zero
  !      nq > 0 and nh1 > 0   zero
  !      nq > 0               h theta D^k S(q_n)[V]       (+ a0 v, k = 1)
  !      nh1 > 0              h (1-theta) D^k S(q_(n-1))[V] (+ a1 v, k = 1)
  !      auxiliaries only     h theta D^k S(q_n)[V] + h (1-theta) D^k S(q_(n-1))[V]
  !
  ! A term whose coefficient is zero is not evaluated. Every factor
  ! is restated on the action's argument before the action is asked;
  ! history(1) stands for the action's state at q_(n-1). Checks,
  ! each stopping the program: the variations are the scheme's; the
  ! action's partial lives on its stated domain; a first-order
  ! direction on the state or history lives on that domain and
  ! matches its width, because a0 v or a1 v is added entry by entry.
  !===================================================================!

  subroutine step_partial_action(this, input_graph, input_data, &
       & variations, output)

    class(scheme), intent(in)                :: this
    class(directed_graph), intent(in)        :: input_graph
    class(field), intent(in)                 :: input_data(:)
    type(variation), intent(in)              :: variations(:)
    class(field), allocatable, intent(inout) :: output

    type(stored_field)        :: out
    type(variation), allocatable :: restated(:)
    type(graph) :: expected
    integer     :: n_expected, num_components, k, j
    integer     :: nq, nh1, nh2
    real(dp), allocatable :: y(:), v(:), partial(:), qold(:)
    type(argument) :: q_n, h1, h2

    call this % require_owned(variations)

    if (.not. allocated(this % action)) then
       error stop 'step: the action is attached'
    end if
    if (size(input_data) < 1) then
       error stop 'step: the state is given'
    end if

    call this % action % domain(input_graph, expected, n_expected)
    num_components = input_data(1) % num_components()

    k = size(variations)

    q_n = this % state()
    h1  = this % history(1)
    nq  = 0
    nh1 = 0
    nh2 = 0
    do j = 1, k
       if (variations(j) % argument_is(q_n)) nq = nq + 1
       if (variations(j) % argument_is(h1))  nh1 = nh1 + 1
    end do
    if (this % reach >= 2) then
       h2 = this % history(2)
       do j = 1, k
          if (variations(j) % argument_is(h2)) nh2 = nh2 + 1
       end do
    end if

    allocate(y(n_expected * num_components))
    y = 0.0_dp

    ! every factor on the action's own argument; history(1) becomes
    ! the action's state, evaluated one instant back
    allocate(restated(k))
    if (nh2 == 0) then
       do j = 1, k
          restated(j) = variations(j) % with_argument( &
               & this % action_argument(variations(j) % argument()))
       end do
    end if

    if (nh2 > 0) then

       if (k == 1) then
          call first_order_direction(variations(1), expected, size(y), v)
          y = this % a2 * v
       end if

    else if (nq > 0 .and. nh1 > 0) then

       continue

    else if (nq > 0) then

       if (this % theta /= 0.0_dp) then
          call action_partial_at_newest(this, input_graph, input_data, expected, &
               & restated, partial)
          y = this % hs * this % theta * partial
       end if
       if (k == 1) then
          call first_order_direction(variations(1), expected, size(y), v)
          y = y + this % a0 * v
       end if

    else if (nh1 > 0) then

       if (this % theta /= 1.0_dp) then
          call history_values(this, 1, input_data, expected, qold)
          call action_partial_at(this, input_graph, input_data, expected, &
               & n_expected, num_components, qold, restated, partial)
          y = this % hs * (1.0_dp - this % theta) * partial
       end if
       if (k == 1) then
          call first_order_direction(variations(1), expected, size(y), v)
          y = y + this % a1 * v
       end if

    else

       if (this % theta /= 0.0_dp) then
          call action_partial_at_newest(this, input_graph, input_data, expected, &
               & restated, partial)
          y = y + this % hs * this % theta * partial
       end if
       if (this % theta /= 1.0_dp) then
          call history_values(this, 1, input_data, expected, qold)
          call action_partial_at(this, input_graph, input_data, expected, &
               & n_expected, num_components, qold, restated, partial)
          y = y + this % hs * (1.0_dp - this % theta) * partial
       end if

    end if

    out = stored_field('step tangent', expected, n_expected, &
         & num_components=num_components)
    call out % set_real_vector(y)

    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine step_partial_action

  !===================================================================!
  ! The action's mixed partial at the newest state, over the
  ! caller's first m inputs.
  !===================================================================!

  subroutine action_partial_at_newest(this, input_graph, input_data, expected, &
       & restated, partial)

    class(scheme), intent(in)          :: this
    class(directed_graph), intent(in)  :: input_graph
    class(field) , intent(in)          :: input_data(:)
    type(graph)  , intent(in)          :: expected
    type(variation), intent(in)        :: restated(:)
    real(dp), allocatable, intent(out) :: partial(:)

    class(field), allocatable :: answer
    type(graph) :: given
    integer :: m

    m = min(this % action % num_arguments(), size(input_data))

    call this % action % partial_action(input_graph, input_data(1:m), restated, answer)
    given = answer % domain()
    if (.not. given % same_as(expected)) then
       error stop 'step: the action partial lives on its stated domain'
    end if
    call answer % real_vector(partial)

  end subroutine action_partial_at_newest

  !===================================================================!
  ! The action's mixed partial at another state: the given values
  ! as the state, the caller's auxiliaries beside it.
  !===================================================================!

  subroutine action_partial_at(this, input_graph, input_data, expected, &
       & n_expected, num_components, values, restated, partial)

    class(scheme), intent(in)          :: this
    class(directed_graph), intent(in)  :: input_graph
    class(field) , intent(in)          :: input_data(:)
    type(graph)  , intent(in)          :: expected
    integer      , intent(in)          :: n_expected, num_components
    real(dp)     , intent(in)          :: values(:)
    type(variation), intent(in)        :: restated(:)
    real(dp), allocatable, intent(out) :: partial(:)

    type(stored_field), allocatable :: inputs(:)
    class(field), allocatable :: answer
    type(graph) :: given

    call inputs_at(this, input_data, expected, n_expected, num_components, values, inputs)

    call this % action % partial_action(input_graph, inputs, restated, answer)
    given = answer % domain()
    if (.not. given % same_as(expected)) then
       error stop 'step: the action partial lives on its stated domain'
    end if
    call answer % real_vector(partial)

  end subroutine action_partial_at

  !===================================================================!
  ! The direction of a first-order factor on the state or a history
  ! state: on the action's domain, of the residual's width.
  !===================================================================!

  subroutine first_order_direction(factor, expected, width, v)

    type(variation), intent(in)        :: factor
    type(graph)    , intent(in)        :: expected
    integer        , intent(in)        :: width
    real(dp), allocatable, intent(out) :: v(:)

    type(graph) :: given

    given = factor % domain()
    if (.not. given % same_as(expected)) then
       error stop 'step: the direction lives on the action''s own domain'
    end if
    call factor % direction(v)
    if (size(v) /= width) then
       error stop 'step: the direction matches the action result'
    end if

  end subroutine first_order_direction

end module operation_step
