!=====================================================================!
! Time integration over a chain graph: instants are vertices,
! steps are edges, and a step size is a number stored per edge.
! march takes an optional per-edge step array; absent, every edge
! uses the marcher's one step.
!
! Rules:
!
!      MARCH_FORWARD    q <- q - h * action(q)          explicit
!      MARCH_BACKWARD   q - qold + h * action(q) = 0    implicit,
!                                                       one solve
!                                                       per edge
!      MARCH_BDF2       (3q - 4qold + qolder)/2
!                              + h * action(q) = 0      implicit,
!                                                       started by
!                                                       one
!                                                       backward
!                                                       step
!
! The implicit rules use the held minimizer (inner) on one
! scheme per edge; every BDF coefficient comes from
! scheme's set_bdf, and configure_edge is the one place a scheme
! is set up for an edge. The action returns minus the velocity.
!
! march_adjoint runs the chain in reverse over the recorded
! trajectory; the adjoint of either rule is derived from the primal,
! never supplied. Under MARCH_FORWARD the tangent of the action at
! each recorded state is compiled to a stencil and transposed. Under
! the implicit rules it is backward substitution: the tangent of the
! step equation, tangent_of(scheme), is compiled, transposed, and
! solved by the dense direct minimizer; the couplings to earlier
! instants are the scheme's constant coefficients a1 and a2.
!
! march_directional computes forward directional derivatives of
! any order along the recorded trajectory. march_adaptive chooses
! the steps at run time under a step_policy.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_marching

  use iso_fortran_env    , only : dp => REAL64
  use operation_action, only : operation, argument
  use view_directed, only : directed_graph
  use field_calculus, only : field
  use graph_fractal      , only : graph
  use field_stored  , only : stored_field
  use view_directed_stored        , only : stored_directed_graph
  use operation_step   , only : scheme, bdf
  use operation_minimization , only : minimizer
  use operation_linearization, only : linearization, tangent_of
  use operation_chain_rule, only : chain_rule, argument_path, &
       & path_derivative
  use operation_step_policy, only : step_policy
  use operation_stencil, only : stencil
  use operation_dense_direct, only : dense_direct

  implicit none

  private
  public :: marcher
  public :: MARCH_FORWARD, MARCH_BACKWARD, MARCH_BDF2

  integer, parameter :: MARCH_FORWARD  = 1
  integer, parameter :: MARCH_BACKWARD = 2
  integer, parameter :: MARCH_BDF2     = 3

  type :: marcher

     integer  :: rule = MARCH_FORWARD
     real(dp) :: step = 1.0_dp
     real(dp) :: singular_tolerance = 1.0e-14_dp

     class(minimizer), allocatable :: inner

   contains

     procedure :: instants
     procedure :: march
     procedure :: march_adjoint
     procedure :: march_directional
     procedure :: march_adaptive

  end type marcher

contains

  !===================================================================!
  ! Build the chain graph: one vertex per instant, one edge per
  ! step.
  !===================================================================!

  subroutine instants(this, nsteps, chain)

    class(marcher), intent(in) :: this
    integer, intent(in) :: nsteps
    type(stored_directed_graph), intent(out) :: chain

    integer :: n

    associate (u1 => this); end associate

    chain = stored_directed_graph(nsteps + 1, &
         & tails=[(n, n = 1, nsteps)], heads=[(n + 1, n = 1, nsteps)])

  end subroutine instants

  !===================================================================!
  ! Integrate nsteps edges forward under the marcher's rule,
  ! optionally recording the state at every instant.
  !===================================================================!

  subroutine march(this, action, on, q, nsteps, steps, trajectory)

    class(marcher), intent(inout)      :: this
    class(operation), intent(in) :: action
    class(directed_graph), intent(in)           :: on
    real(dp), intent(inout)            :: q(:)
    integer, intent(in)                :: nsteps
    real(dp), intent(in), optional     :: steps(:)
    real(dp), allocatable, intent(out), optional :: trajectory(:,:)

    type(stored_directed_graph) :: chain
    type(scheme) :: statement
    type(graph) :: state_domain
    integer         :: n_state_domain
    real(dp), allocatable :: s(:), qold(:), qolder(:), zeros(:)
    real(dp) :: achieved, h_edge, h_previous
    integer :: e, num_components

    call require_valid_steps(steps, nsteps)

    call this % instants(nsteps, chain)

    if (present(trajectory)) then
       allocate(trajectory(size(q), nsteps + 1))
       trajectory(:, 1) = q
    end if

    if (this % rule == MARCH_FORWARD) then

       do e = 1, chain % num_edges()
          call read_statement(action, on, q, s)
          q = q - edge_step(this, steps, e) * s
          if (present(trajectory)) trajectory(:, e + 1) = q
       end do
       return

    end if

    ! Implicit rules: one minimizer solve per edge. BDF2 needs two
    ! history states, so its first edge is a backward-euler step:
    ! qolder is unallocated until the second edge, and configure_edge
    ! reads an absent qolder as order 1.
    statement = bdf(1, action, this % step)

    h_previous = this % step

    ! the unknown of each implicit solve is the state, so the
    ! solve is attached on the action's domain, read once here
    call read_state_domain(action, on, q, state_domain, n_state_domain, num_components)

    allocate(zeros(size(q)))
    zeros = 0.0_dp

    qold = q

    do e = 1, chain % num_edges()

       h_edge = edge_step(this, steps, e)

       call configure_edge(this, statement, h_edge, qold, h_previous, qolder)

       ! num_components is the state's component count, so a multi-component
       ! entry is solved whole
       call this % inner % attach(statement, on, state_domain, &
            & n_state_domain, num_components = num_components)
       call this % inner % solve(zeros, q, achieved)

       qolder     = qold
       qold       = q
       h_previous = h_edge

       if (present(trajectory)) trajectory(:, e + 1) = q

    end do

  end subroutine march

  !===================================================================!
  ! Check an optional per-edge step array: if present it must hold
  ! exactly one positive value per edge. A wrong length or a
  ! nonpositive step stops the program, because a silent default
  ! would integrate a different chain than the caller specified.
  !===================================================================!

  pure subroutine require_valid_steps(steps, nsteps)

    real(dp), intent(in), optional :: steps(:)
    integer , intent(in)           :: nsteps

    integer :: e

    if (.not. present(steps)) return

    if (size(steps) /= nsteps) then
       error stop 'marcher: one step per edge'
    end if

    do e = 1, size(steps)
       if (steps(e) <= 0.0_dp) then
          error stop 'marcher: a time step is positive'
       end if
    end do

  end subroutine require_valid_steps

  pure function edge_step(this, steps, e) result(h)

    class(marcher), intent(in)     :: this
    real(dp), intent(in), optional :: steps(:)
    integer , intent(in)           :: e
    real(dp) :: h

    if (present(steps)) then
       h = steps(e)
    else
       h = this % step
    end if

  end function edge_step

  !===================================================================!
  ! The scheme for one edge: order 2 with the previous step when the
  ! rule is MARCH_BDF2 and a previous state and step are given,
  ! order 1 otherwise; qold is the state the edge starts from. Every
  ! coefficient comes from set_bdf.
  !===================================================================!

  subroutine configure_edge(this, statement, h, qold, h_previous, qolder)

    class(marcher), intent(in)     :: this
    type(scheme), intent(inout)    :: statement
    real(dp), intent(in)           :: h
    real(dp), intent(in)           :: qold(:)
    real(dp), intent(in), optional :: h_previous
    real(dp), intent(in), optional :: qolder(:)

    if (this % rule == MARCH_BDF2 .and. present(qolder) .and. &
         & present(h_previous)) then
       call statement % set_bdf(2, [h, h_previous])
       statement % qolder = qolder
    else
       call statement % set_bdf(1, [h])
    end if

    statement % qold = qold

  end subroutine configure_edge

  !===================================================================!
  ! The scheme for a recorded edge: the step of edge e with the
  ! state at instant e as qold and, past the first edge, the
  ! previous step and the state at instant e - 1 as qolder.
  !===================================================================!

  subroutine recorded_edge(this, statement, steps, e, trajectory)

    class(marcher), intent(in)     :: this
    type(scheme), intent(inout)    :: statement
    real(dp), intent(in), optional :: steps(:)
    integer , intent(in)           :: e
    real(dp), intent(in)           :: trajectory(:,:)

    if (e > 1) then
       call configure_edge(this, statement, edge_step(this, steps, e), &
            & trajectory(:, e), edge_step(this, steps, e - 1), &
            & trajectory(:, e - 1))
    else
       call configure_edge(this, statement, edge_step(this, steps, e), &
            & trajectory(:, e))
    end if

  end subroutine recorded_edge

  !===================================================================!
  ! Reverse traversal over the recorded trajectory. The adjoint of
  ! either rule is derived from the primal: at each edge the tangent
  ! of the statement the forward step took is compiled to a stencil
  ! and transposed. Under MARCH_FORWARD the reverse step is
  !
  !      lambda <- lambda - h_e S'(q_e)^T lambda,
  !
  ! edge by edge in reverse order; under the implicit rules it is
  ! backward substitution. seeds(:, k), if present, is added when
  ! instant k is reached. The trajectory must hold one state per
  ! instant and seeds one entry per instant; both are checked and
  ! stop the program, because a misaligned array would pair states
  ! with the wrong instants.
  !===================================================================!

  subroutine march_adjoint(this, action, on, lambda, nsteps, trajectory, &
       & steps, seeds)

    class(marcher), intent(inout)      :: this
    class(operation), intent(in)       :: action
    class(directed_graph), intent(in)  :: on
    real(dp), intent(inout)            :: lambda(:)
    integer, intent(in)                :: nsteps
    real(dp), intent(in)               :: trajectory(:,:)
    real(dp), intent(in), optional     :: steps(:)
    real(dp), intent(in), optional     :: seeds(:,:)

    type(stored_directed_graph) :: chain
    type(linearization) :: tangent
    type(stencil) :: compiled, adjoint
    real(dp), allocatable :: s(:)
    integer :: e, n

    call require_valid_steps(steps, nsteps)

    call this % instants(nsteps, chain)

    n = size(lambda)

    if (size(trajectory, 1) /= n .or. &
         & size(trajectory, 2) /= chain % num_vertices()) then
       error stop 'march_adjoint: the trajectory carries one state per instant'
    end if

    if (present(seeds)) then
       if (size(seeds, 1) /= n .or. &
            & size(seeds, 2) /= chain % num_vertices()) then
          error stop 'march_adjoint: the seeds carry one entry per instant'
       end if
    end if

    if (this % rule /= MARCH_FORWARD) then
       call substitute_backward(this, action, on, lambda, chain, steps, &
            & trajectory, seeds)
       return
    end if

    ! the explicit step read the state it started from, so its
    ! tangent is frozen there
    tangent = tangent_of(action)

    do e = chain % num_edges(), 1, -1

       if (e < chain % num_edges() .and. present(seeds)) then
          lambda = lambda + seeds(:, e + 1)
       end if

       call tangent % freeze(trajectory(:, e))
       compiled = stencil(tangent, on, n)
       adjoint  = compiled % transpose()
       call read_statement(adjoint, adjoint % pattern, lambda, s)
       lambda = lambda - edge_step(this, steps, e) * s

    end do

    if (present(seeds)) lambda = lambda + seeds(:, 1)

  end subroutine march_adjoint

  !===================================================================!
  ! Backward substitution over the chain: at each edge, in reverse
  ! order, solve
  !
  !      J_e^T lambda_e = seed_e,      J_e = tangent_of(scheme),
  !
  ! with the tangent of the step equation taken at the recorded
  ! state, compiled to a stencil, transposed, and solved by the
  ! dense direct minimizer. The couplings to earlier instants are
  ! the scheme's constant coefficients a1 and a2, so the seed
  ! carried to edge e-1 is -a1 lambda_e plus the -a2 term from the
  ! edge after it. On return, lambda holds the sensitivity at the
  ! first instant.
  !===================================================================!

  subroutine substitute_backward(this, action, on, lambda, chain, steps, &
       & trajectory, seeds)

    class(marcher), intent(inout)      :: this
    class(operation), intent(in) :: action
    class(directed_graph), intent(in)  :: on
    real(dp), intent(inout)            :: lambda(:)
    type(stored_directed_graph), intent(in) :: chain
    real(dp), intent(in), optional     :: steps(:)
    real(dp), intent(in)               :: trajectory(:,:)
    real(dp), intent(in), optional     :: seeds(:,:)

    type(scheme)       :: statement
    type(dense_direct) :: direct
    type(linearization) :: tangent
    type(stencil) :: compiled, adjoint
    real(dp), allocatable :: seed(:), lambda_e(:), carry_one(:), carry_two(:)
    real(dp) :: achieved
    integer :: e, n

    n = size(lambda)

    statement = bdf(1, action, this % step)
    direct % singular_tolerance = this % singular_tolerance

    allocate(carry_one(n), carry_two(n), lambda_e(n))
    carry_one = 0.0_dp
    carry_two = 0.0_dp

    seed = lambda

    do e = chain % num_edges(), 1, -1

       call recorded_edge(this, statement, steps, e, trajectory)

       ! seeds(:, k), if given, is added when instant k is reached
       if (e < chain % num_edges()) then
          seed = carry_one
          if (present(seeds)) seed = seed + seeds(:, e + 1)
       end if

       ! the tangent copies the statement, so it is taken after the
       ! edge is configured, and frozen before it is compiled
       tangent = tangent_of(statement)
       call tangent % freeze(trajectory(:, e + 1))
       compiled = stencil(tangent, on, n)
       adjoint  = compiled % transpose()

       call direct % attach(adjoint, adjoint % pattern, &
            & adjoint % pattern % vertex_set(), &
            & adjoint % pattern % num_vertices())
       lambda_e = 0.0_dp
       call direct % solve(seed, lambda_e, achieved)

       ! carry the couplings to the two earlier instants
       carry_one = carry_two - statement % a1 * lambda_e
       carry_two = -statement % a2 * lambda_e

    end do

    lambda = carry_one
    if (present(seeds)) lambda = lambda + seeds(:, 1)

  end subroutine substitute_backward

  !===================================================================!
  ! Read the action's domain and check the state fits it: the
  ! domain must be nonempty and size(q) must be a whole multiple
  ! of its entry count; both violations stop the program, because
  ! num_components is derived from that division.
  !===================================================================!

  subroutine read_state_domain(action, on, q, state_domain, n_state_domain, &
       & num_components)

    class(operation), intent(in)  :: action
    class(directed_graph)          , intent(in)  :: on
    real(dp)              , intent(in)  :: q(:)
    type(graph)       , intent(out) :: state_domain
    integer               , intent(out) :: n_state_domain
    integer               , intent(out) :: num_components

    integer :: n

    call action % domain(on, state_domain, n_state_domain)

    n = n_state_domain
    if (n <= 0) then
       error stop 'marcher: the action''s state domain is empty'
    end if
    if (mod(size(q), n) /= 0) then
       error stop 'marcher: the state must carry a whole number per member &
            &of the action''s domain'
    end if

    num_components = size(q) / n

  end subroutine read_state_domain

  !===================================================================!
  ! Evaluate the action at the state q. The result must live on
  ! the same domain as the state, not merely have the same length:
  ! a different domain of equal size would otherwise pass, and
  ! q - h s would mix fields from two domains. Both mismatches
  ! stop the program.
  !===================================================================!

  subroutine read_statement(action, on, q, s)

    class(operation), intent(in) :: action
    class(directed_graph), intent(in)           :: on
    real(dp), intent(in)               :: q(:)
    real(dp), allocatable, intent(out) :: s(:)

    type(stored_field)   :: state
    class(field), allocatable :: output
    type(graph) :: state_domain, given
    integer         :: n_state_domain, n_given
    integer :: num_components

    call read_state_domain(action, on, q, state_domain, n_state_domain, num_components)

    state = stored_field('state', state_domain, n_state_domain, num_components=num_components)
    call state % set_real_vector(q)

    call action % apply(on, [state], output)

    given   = output % domain()
    n_given = output % num_entries()
    if (.not. given % same_as(state_domain)) then
       error stop 'marcher: the action result lives on the state''s domain'
    end if

    call output % real_vector(s)
    if (size(s) /= size(q)) then
       error stop 'marcher: the action result matches the state''s width'
    end if

  end subroutine read_statement

  !===================================================================!
  ! Forward directional derivatives of any order along the
  ! recorded trajectory. Under an implicit rule the derivative of
  ! the step equation at edge e splits into
  !
  !      J_e q_e^(s) = -( a1 qprev^(s) + a2 qolder^(s)
  !                       + degree-s composition over the scheme
  !                         with the order-s state derivative zero ),
  !      J_e = tangent_of(scheme) at the solved state,
  !
  ! because the couplings to earlier instants are the scheme's
  ! constant coefficients while everything else is assembled by
  ! chain_rule over the scheme itself. One tangent per edge is
  ! attached to the dense direct minimizer and solved once per
  ! order; only the right-hand side changes with s. Under
  ! MARCH_FORWARD the same recursion over the action needs no solve.
  !
  ! Parameter paths arrive as argument_path values on input slots
  ! 2 and higher, each covered by a supplied parameter field. A
  ! path on slot 1 stops the program, because the state's path is
  ! what this routine computes; a path on an uncovered slot also
  ! stops the program. The initial state is held fixed, so its
  ! derivatives are zero.
  !===================================================================!

  subroutine march_directional(this, action, on, nsteps, trajectory, order, &
       & sensitivities, steps, parameters, paths)

    class(marcher), intent(inout)               :: this
    class(operation), intent(in)                :: action
    class(directed_graph), intent(in)           :: on
    integer, intent(in)                         :: nsteps
    real(dp), intent(in)                        :: trajectory(:,:)
    integer, intent(in)                         :: order
    real(dp), allocatable, intent(out)          :: sensitivities(:,:,:)
    real(dp), intent(in), optional              :: steps(:)
    type(stored_field), intent(in), optional           :: parameters(:)
    type(argument_path), intent(in), optional   :: paths(:)

    type(stored_directed_graph) :: chain
    type(chain_rule)            :: composer
    type(scheme)                :: statement
    type(dense_direct)          :: direct
    type(linearization) :: tangent
    type(argument_path), allocatable :: assembled(:)
    type(stored_field), allocatable         :: inputs(:)
    class(field), allocatable  :: total_field
    type(stored_field)     :: state
    type(graph) :: state_domain
    real(dp), allocatable :: total(:), b(:), q_s(:)
    real(dp) :: achieved, h_edge
    integer :: n_state_domain, num_components
    integer :: e, s_order, k, n, at, npaths

    if (order < 1) then
       error stop 'march_directional: the order is positive'
    end if

    ! degree s always requests the order-s partial (the all-ones
    ! partition), so a shallower calculus is refused before any work
    if (action % max_degree() < order) then
       error stop 'march_directional: the action''s max_degree covers the &
            &requested order'
    end if

    call require_valid_steps(steps, nsteps)
    call this % instants(nsteps, chain)

    n = size(trajectory, 1)
    if (size(trajectory, 2) /= chain % num_vertices()) then
       error stop 'march_directional: the trajectory carries one state per instant'
    end if

    call read_state_domain(action, on, trajectory(:, 1), state_domain, &
         & n_state_domain, num_components)

    ! check the parameter paths: each must name an auxiliary argument
    ! of the action that a supplied parameter field covers -
    ! parameters(j) is the field for the action's argument 1 + j
    npaths = 0
    if (present(paths)) then
       npaths = size(paths)
       do k = 1, npaths
          if (paths(k) % wrt % matches(action % argument(1))) then
             error stop 'march_directional: the state path is computed, &
                  &not supplied'
          end if
          if (.not. present(parameters)) then
             error stop 'march_directional: a parameter path names a supplied argument'
          end if
          if (.not. covered_by_parameters(action, paths(k) % wrt, size(parameters))) then
             error stop 'march_directional: a parameter path names a supplied argument'
          end if
       end do
    end if

    allocate(sensitivities(n, order, nsteps + 1))
    sensitivities = 0.0_dp

    if (this % rule /= MARCH_FORWARD) then
       statement = bdf(1, action, this % step)
       direct % singular_tolerance = this % singular_tolerance
       allocate(q_s(n))
    end if

    ! for each edge, assemble the composition degree by degree and
    ! advance each order's derivative
    do e = 1, chain % num_edges()

       h_edge = edge_step(this, steps, e)

       ! the composition is evaluated at the solved instant under
       ! an implicit rule and at the previous instant under the
       ! explicit rule, matching the state each rule differentiates
       if (this % rule == MARCH_FORWARD) then
          at = e
       else
          at = e + 1
          call recorded_edge(this, statement, steps, e, trajectory)
          ! the tangent copies the statement, so it is taken after
          ! the edge is configured, and frozen before it is attached
          tangent = tangent_of(statement)
          call tangent % freeze(trajectory(:, e + 1))
          call direct % attach(tangent, on, state_domain, n_state_domain, &
               & num_components = num_components)
       end if

       state = stored_field('state', state_domain, n_state_domain, num_components=num_components)
       call state % set_real_vector(trajectory(:, at))

       if (allocated(inputs)) deallocate(inputs)
       if (present(parameters)) then
          allocate(inputs(1 + size(parameters)))
          inputs(1) = state
          do k = 1, size(parameters)
             inputs(1 + k) = parameters(k)
          end do
       else
          allocate(inputs(1))
          inputs(1) = state
       end if

       do s_order = 1, order

          ! the state path names the state argument of the statement
          ! the composition runs over; under an implicit rule the
          ! parameter paths are restated in the scheme's space
          if (this % rule == MARCH_FORWARD) then
             call build_paths(this, sensitivities, state_domain, &
                  & n_state_domain, num_components, at, s_order, npaths, paths, &
                  & action % argument(1), assembled)
          else
             call build_paths(this, sensitivities, state_domain, &
                  & n_state_domain, num_components, at, s_order, npaths, paths, &
                  & statement % state(), assembled, statement)
          end if

          if (this % rule == MARCH_FORWARD) then

             call composer % assemble(action, on, inputs, s_order, assembled, &
                  & total_field)
             call total_field % real_vector(total)

             sensitivities(:, s_order, e + 1) = &
                  & sensitivities(:, s_order, e) - h_edge * total

          else

             ! over the scheme, with the order-s state derivative
             ! zero, the total is hs times the action's composition
             call composer % assemble(statement, on, inputs, s_order, &
                  & assembled, total_field)
             call total_field % real_vector(total)

             b = statement % a1 * sensitivities(:, s_order, e) + total
             if (e > 1) then
                b = b + statement % a2 * sensitivities(:, s_order, e - 1)
             end if

             q_s = 0.0_dp
             call direct % solve(-b, q_s, achieved)
             sensitivities(:, s_order, e + 1) = q_s

          end if

       end do

    end do

  end subroutine march_directional

  !===================================================================!
  ! Whether a parameter path's argument is one of the action's
  ! auxiliaries that a supplied parameter field covers: parameters(j)
  ! stands for the action's argument 1 + j.
  !===================================================================!

  logical function covered_by_parameters(action, wrt, num_parameters) result(covered)

    class(operation), intent(in) :: action
    type(argument)  , intent(in) :: wrt
    integer         , intent(in) :: num_parameters

    integer :: j

    covered = .false.
    do j = 1, min(num_parameters, action % num_arguments() - 1)
       if (wrt % matches(action % argument(1 + j))) covered = .true.
    end do

  end function covered_by_parameters

  !===================================================================!
  ! Build the argument paths for one edge and one order: the state
  ! path, on the given state argument, holds the solved derivatives
  ! below the current order; under an implicit rule its order-s
  ! entry is set to zero, which makes the assembled total the
  ! right-hand side for the unknown q^(s). The caller's parameter
  ! paths are appended, restated in the scheme's argument space when
  ! the composition runs over a scheme.
  !===================================================================!

  subroutine build_paths(this, sensitivities, state_domain, n_state_domain, &
       & num_components, at, s_order, npaths, parameter_paths, state_argument, &
       & assembled, statement)

    class(marcher), intent(in) :: this
    real(dp)      , intent(in) :: sensitivities(:,:,:)
    type(graph), intent(in) :: state_domain
    integer       , intent(in) :: n_state_domain, num_components, at, s_order, npaths
    type(argument_path), intent(in), optional :: parameter_paths(:)
    type(argument), intent(in) :: state_argument
    type(argument_path), allocatable, intent(out) :: assembled(:)
    type(scheme), intent(in), optional :: statement

    type(stored_field) :: derivative_field
    real(dp), allocatable :: derivative_values(:)
    integer :: k

    allocate(assembled(1 + npaths))

    assembled(1) % wrt = state_argument
    allocate(assembled(1) % derivative(s_order))

    do k = 1, s_order

       if (k == s_order .and. this % rule /= MARCH_FORWARD) then
          ! the unknown order's entry is zero while the total is
          ! assembled
          allocate(derivative_values(size(sensitivities, 1)))
          derivative_values = 0.0_dp
       else
          derivative_values = sensitivities(:, k, at)
       end if

       derivative_field = stored_field('state path', state_domain, n_state_domain, &
            & num_components=num_components)
       call derivative_field % set_real_vector(derivative_values)
       assembled(1) % derivative(k) % occupied  = .true.
       assembled(1) % derivative(k) % direction = derivative_field
       deallocate(derivative_values)

    end do

    do k = 1, npaths
       assembled(1 + k) = parameter_paths(k)
       if (present(statement)) then
          assembled(1 + k) % wrt = statement % from_action(parameter_paths(k) % wrt)
       end if
    end do

  end subroutine build_paths

  !===================================================================!
  ! Adaptive integration: the policy proposes a step, the marcher
  ! computes a trial state under its rule without modifying q,
  ! measures the distance between the trial and the extrapolating
  ! predictor through the accepted history, and the policy judges.
  ! An accepted edge commits the trial and records its step; a
  ! rejected edge changes nothing. When max_attempts trials in a
  ! row are rejected the routine returns with completed = .false.
  ! and q at the last accepted state.
  !
  ! Checks, each stopping the program: duration and max_attempts
  ! must be positive, and the policy must propose positive steps.
  ! What is returned is steps_taken; a caller wanting the
  ! trajectory or the adjoint marches again with those steps.
  !===================================================================!

  subroutine march_adaptive(this, action, on, q, duration, policy, &
       & max_attempts, steps_taken, completed)

    class(marcher), intent(inout)      :: this
    class(operation), intent(in) :: action
    class(directed_graph), intent(in)  :: on
    real(dp), intent(inout)            :: q(:)
    real(dp), intent(in)               :: duration
    class(step_policy), intent(inout)  :: policy
    integer, intent(in)                :: max_attempts
    real(dp), allocatable, intent(out) :: steps_taken(:)
    logical, intent(out)               :: completed

    type(scheme) :: statement
    type(graph)     :: state_domain
    integer             :: n_state_domain, num_components
    real(dp), allocatable :: trial(:), predictor(:), qprev(:), s(:), zeros(:)
    real(dp), allocatable :: grown(:)
    real(dp) :: t, h, h_previous, estimate, achieved
    integer  :: attempt, taken
    logical  :: accepted, have_previous

    if (duration <= 0.0_dp) then
       error stop 'march_adaptive: the duration is positive'
    end if
    if (max_attempts < 1) then
       error stop 'march_adaptive: the attempt budget is positive'
    end if

    call read_state_domain(action, on, q, state_domain, n_state_domain, num_components)

    if (this % rule /= MARCH_FORWARD) then
       statement = bdf(1, action, this % step)
       allocate(zeros(size(q)))
       zeros = 0.0_dp
    end if

    allocate(steps_taken(0))
    t             = 0.0_dp
    taken         = 0
    have_previous = .false.
    h_previous    = 0.0_dp
    completed     = .false.

    do while (duration - t > 1.0e-12_dp * duration)

       call policy % propose(h)

       accepted = .false.
       attempt  = 0

       do while (.not. accepted .and. attempt < max_attempts)

          attempt = attempt + 1

          if (h <= 0.0_dp) then
             error stop 'march_adaptive: the policy proposes a positive step'
          end if
          h = min(h, duration - t)

          ! compute the trial state; q is not modified
          if (this % rule == MARCH_FORWARD) then

             call read_statement(action, on, q, s)
             trial = q - h * s

          else

             ! qprev is unallocated until the first accepted edge, so
             ! configure_edge reads it as absent and takes order 1
             call configure_edge(this, statement, h, q, h_previous, qprev)

             trial = q
             call this % inner % attach(statement, on, state_domain, &
                  & n_state_domain, num_components = num_components)
             call this % inner % solve(zeros, trial, achieved)

          end if

          ! error estimate: distance from the extrapolating
          ! predictor - constant with one accepted state behind,
          ! linear with two. The accept decision is the policy's.
          if (have_previous) then
             predictor = q + (h / h_previous) * (q - qprev)
          else
             predictor = q
          end if
          estimate = norm2(trial - predictor)

          call policy % judge(estimate, h, attempt, accepted)

          if (.not. accepted .and. attempt < max_attempts) then
             call policy % retry(estimate, h)
          end if

       end do

       if (.not. accepted) return

       qprev         = q
       q             = trial
       t             = t + h
       h_previous    = h
       have_previous = .true.

       allocate(grown(taken + 1))
       grown(1:taken) = steps_taken
       grown(taken + 1) = h
       call move_alloc(grown, steps_taken)
       taken = taken + 1

    end do

    completed = .true.

  end subroutine march_adaptive

end module operation_marching
