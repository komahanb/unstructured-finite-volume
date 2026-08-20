!=====================================================================!
! The step operator: one step of a time recurrence as an
! operation,
!
!      a0 q  +  a1 qold  +  a2 qolder  +  hs S(q),
!
! the temporal counterpart of the spatial stencil. Constructors:
!
!      backward_euler(action, h)          reach 1, [1, -1]
!      bdf(k, action, h)                  reach k, uniform
!                                         coefficients
!      bdf_variable(k, action, steps)     reach k, coefficients
!                                         from the steps actually
!                                         taken - steps(1) = h_n,
!                                         steps(2) = h_{n-1}
!
! Every coefficient table lives in set_bdf, which reconfigures a
! standing step operator between edges; nothing else defines a BDF
! coefficient. At equal steps set_bdf assigns the uniform
! constants directly, with no arithmetic. A diagonal RK stage and
! an Adams corrector are this same backward step with a scaled hs
! against an externally assembled base (hs = h gamma, hs = h
! beta0), so they need no constructors here. The action returns
! minus the velocity.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_graph_step

  use iso_fortran_env    , only : dp => REAL64
  use operation_action, only : graph_operation
  use view_directed, only : directed_graph
  use field_calculus, only : graph_field
  use graph_fractal      , only : set_graph => graph
  use graph_discretization     , only : discretization_operator
  use field_stored  , only : field
  use view_directed_stored        , only : directed_stored_graph

  implicit none

  private
  public :: step_operator
  public :: backward_euler, bdf, bdf_variable

  type, extends(discretization_operator) :: step_operator

     class(graph_operation), allocatable :: action

     real(dp), allocatable :: qold(:)
     real(dp), allocatable :: qolder(:)

     real(dp) :: a0 = 1.0_dp
     real(dp) :: a1 = -1.0_dp
     real(dp) :: a2 = 0.0_dp
     real(dp) :: hs = 0.0_dp

     integer :: reach = 1

   contains

     procedure :: name         => step_name
     procedure :: domain       => step_domain
     procedure :: apply        => step_apply
     procedure :: dependencies => step_dependencies
     procedure :: set_bdf

  end type step_operator

contains

  !===================================================================!
  ! Backward euler: one instant back.
  !===================================================================!

  function backward_euler(action, h) result(this)

    class(graph_operation), intent(in) :: action
    real(dp), intent(in)               :: h
    type(step_operator)                :: this

    allocate(this % action, source=action)
    call this % set_bdf(1, [h])

  end function backward_euler

  !===================================================================!
  ! BDF of order k at a uniform step. Only orders 1 and 2 have
  ! coefficient tables; any other order stops the program.
  !===================================================================!

  function bdf(k, action, h) result(this)

    integer, intent(in)                :: k
    class(graph_operation), intent(in) :: action
    real(dp), intent(in)               :: h
    type(step_operator)                :: this

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
    class(graph_operation), intent(in) :: action
    real(dp), intent(in)               :: steps(:)
    type(step_operator)                :: this

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
  ! equal k; every step must be positive.
  !===================================================================!

  subroutine set_bdf(this, k, steps)

    class(step_operator), intent(inout) :: this
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
  ! instant follows which. A stencil_operator's dependencies are
  ! the same pattern on the dependent-variable axis.
  !===================================================================!

  subroutine step_dependencies(this, pattern)

    class(step_operator), intent(in)       :: this
    class(directed_graph), allocatable, intent(out) :: pattern

    integer :: n, newest

    newest = this % reach + 1

    allocate(pattern, source=directed_stored_graph(newest, &
         & tails=[(n, n = 1, newest)], &
         & heads=[(newest, n = 1, newest)]))

  end subroutine step_dependencies

  !===================================================================!
  ! The operation interface.
  !===================================================================!

  pure function step_name(this) result(name)

    class(step_operator), intent(in) :: this
    character(len=:), allocatable :: name

    associate (u1 => this); end associate

    name = 'step'

  end function step_name

  !===================================================================!
  ! The step's domain is the domain of the action it discretizes:
  ! the residual is a statement about the same unknown S is a
  ! statement about, so the question is delegated to the action.
  !===================================================================!

  subroutine step_domain(this, input_graph, domain, nentries)

    class(step_operator), intent(in)       :: this
    class(directed_graph), intent(in)               :: input_graph
    type(set_graph), intent(out) :: domain
    integer        , intent(out) :: nentries

    call this % action % domain(input_graph, domain, nentries)

  end subroutine step_domain

  !===================================================================!
  ! Evaluate the residual on the action's domain. Two identity
  ! checks stop the program: the state must live on the action's
  ! domain, and the action's result must land on that same domain.
  ! Equal length alone would let a field from a different domain
  ! through and produce plausible but wrong numbers.
  !===================================================================!

  subroutine step_apply(this, input_graph, input_data, output)

    class(step_operator), intent(in)               :: this
    class(directed_graph), intent(in)                       :: input_graph
    class(graph_field), intent(in), optional       :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output

    type(field)   :: out
    class(graph_field), allocatable :: velocity
    type(set_graph) :: expected, given
    integer         :: n_expected, n_given
    real(dp), allocatable :: q(:), s(:), y(:)
    integer :: ncomp

    call this % action % domain(input_graph, expected, n_expected)

    if (present(input_data)) then

       given   = input_data(1) % domain()
       n_given = input_data(1) % num_entries()
       if (.not. given % same_as(expected)) then
          error stop 'step: the state must live on the action''s own domain'
       end if
       ncomp = input_data(1) % num_components()

       call input_data(1) % get_real_vector(q)

       call this % action % apply(input_graph, input_data, velocity)
       given   = velocity % domain()
       n_given = velocity % num_entries()
       if (.not. given % same_as(expected)) then
          error stop 'step: the action result lives on its stated domain'
       end if
       call velocity % get_real_vector(s)

       y = this % a0 * q + this % a1 * this % qold + this % hs * s
       if (allocated(this % qolder) .and. abs(this % a2) > 0.0_dp) then
          y = y + this % a2 * this % qolder
       end if

    else
       ncomp = 1
       allocate(y(n_expected))
       y = 0.0_dp
    end if

    out = field('step residual', expected, n_expected, ncomp=ncomp)
    call out % set_real_vector(y)

    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine step_apply

end module class_graph_step
