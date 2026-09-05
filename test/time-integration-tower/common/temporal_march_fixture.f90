!=====================================================================!
! THE TEMPORAL MARCH FIXTURE - earned at LEVEL 8.
!
! The march is a graph the production driver evaluates (b9944c3):
! a bipartite digraph whose first part is the steps and whose
! second part is the states,
!
!      [ q0 ] --> (s1) --> [ q1 ] --> (s2) --> [ q2 ] --> ...
!                            \___________________^        bdf-2
!
! where an arc from a state into a step is a read and an arc from a
! step into a state is a write. The driver takes the topological
! order of the projection onto the steps, applies the rule stored
! at each, and places the result at the state it writes. No routine
! here records an instant's index.
!
! THE RULE AT A STEP reads the history and returns the next state:
! the temporal step's residual (Level 6) driven to zero by newton
! over gmres (Level 7) for an implicit scheme, and one evaluation of
! the residual at zero for an explicit one, since R(q) = q + R(0)
! when theta = 0.
!
! Three rules, as the deleted marcher (operation_marching, b9944c3)
! stated them:
!
!      MARCH_FORWARD    theta 0, reach 1     q_n = q_(n-1) - h S(q_(n-1))
!      MARCH_BACKWARD   theta 1, reach 1     q_n - q_(n-1) + h S(q_n) = 0
!      MARCH_BDF2       theta 1, reach 2     (3 q_n - 4 q_(n-1) + q_(n-2))/2
!                                                        + h S(q_n) = 0
!
! and bdf-2's first step is a backward-euler step, because one
! history state is all that exists there.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module temporal_march_fixture

  use iso_fortran_env      , only : dp => REAL64
  use operation_action     , only : operation, binding, contract
  use operation_action     , only : bound_value, emit_real
  use view_directed        , only : directed_graph, forward
  use field_calculus       , only : field, FIELD_REAL
  use field_stored         , only : stored_field
  use graph_fractal        , only : graph
  use operation_newton     , only : newton
  use operation_gmres      , only : gmres
  use operation_driver     , only : driver, rule_graph, data_graph, pairing
  use view_read_write      , only : bipartite_digraph, FIRST_PART, SECOND_PART
  use temporal_step_fixture, only : temporal_step, backward_euler, &
       &                            forward_euler, bdf

  implicit none

  private
  public :: advance
  public :: march_incidence, marched
  public :: MARCH_FORWARD, MARCH_BACKWARD, MARCH_BDF2

  integer, parameter :: MARCH_FORWARD  = 1
  integer, parameter :: MARCH_BACKWARD = 2
  integer, parameter :: MARCH_BDF2     = 3

  !===================================================================!
  ! The rule at one step: the temporal step it solves and the
  ! tolerances of the implicit solve. Its arguments are the history,
  ! one per instant the step reaches; its result is the next state.
  !===================================================================!

  type, extends(operation) :: advance

     type(temporal_step) :: step
     real(dp) :: newton_tolerance = 1.0e-13_dp
     real(dp) :: linear_tolerance = 1.0e-14_dp

   contains

     procedure :: name   => advance_name
     procedure :: domain => advance_domain
     procedure :: apply  => advance_apply

  end type advance

  interface advance
     module procedure create_advance
  end interface advance

contains

  function create_advance(step, newton_tolerance, linear_tolerance) result(this)

    type(temporal_step), intent(in) :: step
    real(dp)           , intent(in) :: newton_tolerance, linear_tolerance
    type(advance) :: this

    integer :: j

    this % step             = step
    this % newton_tolerance = newton_tolerance
    this % linear_tolerance = linear_tolerance

    call this % declare_arguments(step % reach, &
         & [(contract(FIELD_REAL, 1), j = 1, step % reach)])

  end function create_advance

  pure function advance_name(this) result(name)

    class(advance), intent(in)    :: this
    character(len=:), allocatable :: name

    associate (u1 => this); end associate

    name = 'advance one step'

  end function advance_name

  subroutine advance_domain(this, input_graph, domain, num_entries)

    class(advance)       , intent(in)  :: this
    class(directed_graph), intent(in)  :: input_graph
    type(graph)          , intent(out) :: domain
    integer              , intent(out) :: num_entries

    call this % step % domain(input_graph, domain, num_entries)

  end subroutine advance_domain

  !===================================================================!
  ! The next state from the history. A missing history state, or one
  ! outside the action's domain, stops the program.
  !===================================================================!

  subroutine advance_apply(this, input_graph, inputs, output)

    class(advance)           , intent(in)    :: this
    class(directed_graph)    , intent(in)    :: input_graph
    type(binding)            , intent(in), optional :: inputs(:)
    class(field), allocatable, intent(inout) :: output

    type(stored_field), allocatable :: history(:)
    class(field)      , allocatable :: given, residual
    type(stored_field) :: zero
    type(newton)       :: governor
    type(graph)        :: expected
    real(dp), allocatable :: v(:), x(:), zeros(:)
    real(dp) :: achieved
    integer  :: n_expected, k

    call this % step % domain(input_graph, expected, n_expected)

    if (.not. present(inputs)) error stop 'advance: a step reads its history'
    if (size(inputs) /= this % step % reach) then
       error stop 'advance: one history state per instant the scheme reaches'
    end if

    allocate(history(this % step % reach))
    do k = 1, this % step % reach
       call bound_value(inputs, this % argument(k), given)
       if (.not. given % defined_on(expected)) then
          error stop 'advance: every history state is defined on the action''s own domain'
       end if
       call given % real_vector(v)
       history(k) = stored_field('history', expected, n_expected, num_components=1)
       call history(k) % set_real_vector(v)
    end do

    allocate(zeros(n_expected))
    zeros = 0.0_dp

    if (this % step % theta == 0.0_dp) then

       ! R(q) = q + R(0): the next state is -R(0)
       zero = stored_field('zero', expected, n_expected, num_components=1)
       call zero % set_real_vector(zeros)
       call this % step % apply(input_graph, this % step % bind([zero, history]), residual)
       call residual % real_vector(x)
       x = -x

    else

       ! newton over gmres, the residual driven to zero from the
       ! previous state
       governor % tolerance = this % newton_tolerance
       allocate(governor % inner, source=gmres())
       governor % inner % tolerance = this % linear_tolerance
       call governor % state(this % step, input_graph, expected, n_expected, &
            & num_components=1, stored_inputs=history)
       call history(1) % real_vector(x)
       call governor % solve(zeros, x, achieved)

    end if

    call emit_real('state', expected, n_expected, x, output, num_components=1)

  end subroutine advance_apply

  !===================================================================!
  ! The march digraph over nsteps steps at the given reach: step k
  ! reads the states k, k-1, ..., k-reach+1 that exist, in that order,
  ! and writes state k+1. State 1 is read and never written: the
  ! initial state, a source of the digraph.
  !===================================================================!

  function march_incidence(nsteps, reach) result(incidence)

    integer, intent(in) :: nsteps, reach
    type(bipartite_digraph) :: incidence

    integer, allocatable :: from_part(:), from_vertex(:), to_part(:), to_vertex(:)
    integer :: k, j, n

    from_part   = [integer ::]
    from_vertex = [integer ::]
    to_part     = [integer ::]
    to_vertex   = [integer ::]

    do k = 1, nsteps
       do j = 1, reach
          n = k - j + 1
          if (n < 1) cycle
          from_part   = [from_part  , SECOND_PART]
          from_vertex = [from_vertex, n]
          to_part     = [to_part    , FIRST_PART]
          to_vertex   = [to_vertex  , k]
       end do
       from_part   = [from_part  , FIRST_PART]
       from_vertex = [from_vertex, k]
       to_part     = [to_part    , SECOND_PART]
       to_vertex   = [to_vertex  , k + 1]
    end do

    incidence = bipartite_digraph(nsteps, nsteps + 1, from_part, from_vertex, &
         & to_part, to_vertex)

  end function march_incidence

  !===================================================================!
  ! The state after nsteps steps of the rule from the initial state,
  ! by the production driver over the march digraph. The host is
  ! passed to every rule and read by none of them; the tolerances
  ! are the implicit solve's, newton's and its inner gmres's.
  !===================================================================!

  function marched(rule, action, h, nsteps, host, initial, &
       & newton_tolerance, linear_tolerance) result(final)

    integer              , intent(in) :: rule
    class(operation)     , intent(in) :: action
    real(dp)             , intent(in) :: h
    integer              , intent(in) :: nsteps
    class(directed_graph), intent(in) :: host
    type(stored_field)   , intent(in) :: initial
    real(dp)             , intent(in) :: newton_tolerance, linear_tolerance
    class(field), allocatable :: final

    type(bipartite_digraph) :: incidence
    type(rule_graph) :: rules
    type(data_graph) :: values
    type(driver)     :: schedule
    type(pairing)    :: pairs
    type(advance)    :: first, later
    integer :: k

    select case (rule)
    case (MARCH_FORWARD)
       first = advance(forward_euler(action, h), newton_tolerance, linear_tolerance)
       later = first
    case (MARCH_BACKWARD)
       first = advance(backward_euler(action, h), newton_tolerance, linear_tolerance)
       later = first
    case (MARCH_BDF2)
       first = advance(backward_euler(action, h), newton_tolerance, linear_tolerance)
       later = advance(bdf(2, action, h), newton_tolerance, linear_tolerance)
    case default
       error stop 'temporal_march: the rule is forward, backward or bdf-2'
    end select

    incidence = march_incidence(nsteps, later % step % reach)

    allocate(rules % at(nsteps), values % at(nsteps + 1))
    allocate(rules % at(1) % rule, source=first)
    do k = 2, nsteps
       allocate(rules % at(k) % rule, source=later)
    end do
    allocate(values % at(1) % datum, source=initial)

    schedule = driver(rules % at(1) % rule, incidence, forward)
    call schedule % pair_with(rules % pair(values))
    call schedule % evaluate(host)

    pairs = schedule % pairing_of()
    call pairs % datum_at(nsteps + 1, final)
    if (.not. allocated(final)) error stop 'temporal_march: the last step writes its state'

  end function marched

end module temporal_march_fixture
