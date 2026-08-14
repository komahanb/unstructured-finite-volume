!=====================================================================!
! The step operator: time's discretization stencil.
!
! LEVEL 1 OF THE STRATIFICATION. One step of a time recurrence, as
! an operation - the temporal concretion of the discretization
! operator, sibling of the spatial stencil:
!
!      a0 q  +  a1 qold  +  a2 qolder  +  h S(q)
!
! The scheme is a MOTIF stamped along the chain: backward euler
! reaches one instant back, bdf-2 reaches two, and the dependency
! pattern this type answers by contract is exactly that motif - a
! little chain of reach + 1 instants. The triangularity a causal
! solver exploits is made here, by the maker, not assumed there, by
! the solver.
!
! The shelf of schemes is this module's constructor roster - the
! module is the library and the compiler is its index:
!
!      backward_euler(action, h)          reach 1, [1, -1]
!      bdf(k, action, h)                  reach k, the k-step table
!
! One name per family, the order riding as an argument, never as a
! type. The statement it holds returns MINUS the velocity, matching
! the house convention that a balance measures what a cell has left
! over.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_graph_step

  use iso_fortran_env    , only : dp => REAL64
  use graph_grammar      , only : graph, graph_field, graph_operation
  use graph_carrier      , only : member_set
  use graph_calculus     , only : discretization_operator
  use class_graph_field  , only : field
  use class_graph        , only : stored_graph

  implicit none

  private
  public :: step_operator
  public :: backward_euler, bdf

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

  end type step_operator

contains

  !===================================================================!
  ! The shelf. Backward euler: one instant back.
  !===================================================================!

  function backward_euler(action, h) result(this)

    class(graph_operation), intent(in) :: action
    real(dp), intent(in)               :: h
    type(step_operator)                :: this

    allocate(this % action, source=action)
    this % a0    = 1.0_dp
    this % a1    = -1.0_dp
    this % a2    = 0.0_dp
    this % hs    = h
    this % reach = 1

  end function backward_euler

  !===================================================================!
  ! The bdf family: k instants back, one name, k absorbed. Order one
  ! is backward euler by another road; order two is the second-order
  ! table; higher orders join here when their tables are ruled.
  !===================================================================!

  function bdf(k, action, h) result(this)

    integer, intent(in)                :: k
    class(graph_operation), intent(in) :: action
    real(dp), intent(in)               :: h
    type(step_operator)                :: this

    allocate(this % action, source=action)
    this % hs    = h
    this % reach = k

    select case (k)
    case (1)
       this % a0 = 1.0_dp
       this % a1 = -1.0_dp
       this % a2 = 0.0_dp
    case (2)
       this % a0 = 1.5_dp
       this % a1 = -2.0_dp
       this % a2 = 0.5_dp
    case default
       error stop 'bdf: only orders one and two carry tables so far'
    end select

  end function bdf

  !===================================================================!
  ! The contract's answer: THE STENCIL ON THE INDEPENDENT AXIS.
  !
  ! A step's residual at the newest instant reads every instant its
  ! table reaches back to,
  !
  !      R_n = a0 q_n + a1 q_(n-1) + a2 q_(n-2) + h S(q_n)
  !
  ! so the dependency is a FAN-IN onto the instant being solved for -
  ! reach + 1 instants, every one of them an arrow into the last:
  !
  !      backward euler        bdf-2
  !          1 --> 2               1 --> 3
  !          2 --> 2               2 --> 3
  !                                3 --> 3
  !
  ! A STENCIL IS NOT A CHRONOLOGY. The succession
  !
  !      1 -> 2 -> 3
  !
  ! is a true and useful relation - it is what the time integration
  ! tower's instants and control chain describe - but it says which
  ! instant FOLLOWS which, not which instants the residual READS. This
  ! contract owes the second. The self-arrow is the implicit part, and
  ! it is what makes the newest instant an unknown rather than data.
  !
  ! The axis here is the INDEPENDENT variable; a stencil_operator's
  ! answer to the same verb is the stencil on the DEPENDENT variable.
  ! One contract, one meaning - the stencil on whichever axis the
  ! concrete type represents.
  !===================================================================!

  subroutine step_dependencies(this, pattern)

    class(step_operator), intent(in)       :: this
    class(graph), allocatable, intent(out) :: pattern

    integer :: n, newest

    newest = this % reach + 1

    allocate(pattern, source=stored_graph(newest, &
         & tails=[(n, n = 1, newest)], &
         & heads=[(newest, n = 1, newest)]))

  end subroutine step_dependencies

  !===================================================================!
  ! The step's three answers as an operation.
  !===================================================================!

  pure function step_name(this) result(name)

    class(step_operator), intent(in) :: this
    character(len=:), allocatable :: name

    associate (u1 => this); end associate

    name = 'step'

  end function step_name

  !===================================================================!
  ! THE DOMAIN OF A STEP IS THE DOMAIN OF THE ACTION IT DISCRETIZES.
  !
  ! A step is an operation BUILT FROM another operation, and its
  ! residual
  !
  !      a0 q + a1 qold + a2 qolder + h S(q)
  !
  ! is a statement about the same unknown S is a statement about.
  ! So the question is delegated: the action answers, and the host
  ! is the conduit that carries it there.
  !
  ! For every action that reads its domain off the graph - which is
  ! all of them on the ordinary-graph road - this returns exactly
  ! what asking the graph returned, so no marching caller sees a
  ! change. For an action that carries its own domain, it returns
  ! that domain, which asking the graph never could.
  !===================================================================!

  subroutine step_domain(this, input_graph, domain)

    class(step_operator), intent(in)       :: this
    class(graph), intent(in)               :: input_graph
    class(member_set), allocatable, intent(out) :: domain

    call this % action % domain(input_graph, domain)

  end subroutine step_domain

  !===================================================================!
  ! The step's residual, on the action's own domain.
  !
  ! Three things travel with the DATA and the ACTION rather than
  ! with the host: the domain, the component width, and the carrier
  ! the answer lands on. The state's width is read from the state -
  ! a coordinate carrying several numbers is measured whole - and
  ! never inferred by dividing by a vertex count.
  !
  ! Both domain checks are identity questions, and both are refusals
  ! the caller wants early: a state on a foreign carrier, or an
  ! action that answers somewhere other than where it said it
  ! would, are wrong in ways that produce plausible numbers if left
  ! alone.
  !===================================================================!

  subroutine step_apply(this, input_graph, input_data, output)

    class(step_operator), intent(in)               :: this
    class(graph), intent(in)                       :: input_graph
    class(graph_field), intent(in), optional       :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output

    type(field)   :: out
    class(graph_field), allocatable :: velocity
    class(member_set), allocatable  :: expected, given
    real(dp), allocatable :: q(:), s(:), y(:)
    integer :: ncomp

    call this % action % domain(input_graph, expected)

    if (present(input_data)) then

       call input_data(1) % domain(given)
       if (.not. given % same_as(expected)) then
          error stop 'step: the state must live on the action''s own domain'
       end if
       ncomp = input_data(1) % num_components()

       call input_data(1) % get_real_vector(q)

       call this % action % apply(input_graph, input_data, velocity)
       call velocity % domain(given)
       if (.not. given % same_as(expected)) then
          error stop 'step: the action must answer on its stated domain'
       end if
       call velocity % get_real_vector(s)

       y = this % a0 * q + this % a1 * this % qold + this % hs * s
       if (allocated(this % qolder) .and. abs(this % a2) > 0.0_dp) then
          y = y + this % a2 * this % qolder
       end if

    else
       ncomp = 1
       allocate(y(expected % size()))
       y = 0.0_dp
    end if

    out = field('step residual', expected, ncomp=ncomp)
    call out % set_real_vector(y)

    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine step_apply

end module class_graph_step
