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
! step_operator per edge; every BDF coefficient comes from
! step_operator's set_bdf. The action returns minus the velocity.
!
! march_adjoint runs the chain in reverse. Under MARCH_FORWARD it
! applies the caller's transposed statement edge by edge. Under
! the implicit rules it is backward substitution: at each edge the
! step Jacobian a0 I + h S'(q_e) is assembled at the recorded
! state via tangent_of and dense_matrix_of, its transpose is
! solved by the dense direct minimizer, and the couplings to
! earlier instants are the constant coefficients a1 and a2. This
! requires the action and the forward trajectory as arguments.
!
! march_directional computes forward directional derivatives of
! any order along the recorded trajectory. march_adaptive chooses
! the steps at run time under a step_policy.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_graph_marcher

  use iso_fortran_env    , only : dp => REAL64
  use graph_operation_view, only : graph_operation
  use graph_directed_view, only : directed_graph
  use graph_field_calculus, only : graph_field
  use fractal_graph      , only : set_graph => graph
  use class_graph_field  , only : field
  use class_graph        , only : directed_stored_graph
  use class_graph_step   , only : step_operator, bdf
  use graph_minimization , only : minimizer
  use graph_discretization     , only : linearization_operator, &
       & differentiable_operation
  use class_graph_linearization, only : tangent_of
  use class_graph_chain_rule, only : chain_rule, argument_path, &
       & path_derivative
  use class_graph_step_policy, only : step_policy
  use class_graph_dense_direct, only : solve_dense_matrix_with_dense_direct, &
       & dense_matrix_of

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
    type(directed_stored_graph), intent(out) :: chain

    integer :: n

    associate (u1 => this); end associate

    chain = directed_stored_graph(nsteps + 1, &
         & tails=[(n, n = 1, nsteps)], heads=[(n + 1, n = 1, nsteps)])

  end subroutine instants

  !===================================================================!
  ! Integrate nsteps edges forward under the marcher's rule,
  ! optionally recording the state at every instant.
  !===================================================================!

  subroutine march(this, action, on, q, nsteps, steps, trajectory)

    class(marcher), intent(inout)      :: this
    class(graph_operation), intent(in) :: action
    class(directed_graph), intent(in)           :: on
    real(dp), intent(inout)            :: q(:)
    integer, intent(in)                :: nsteps
    real(dp), intent(in), optional     :: steps(:)
    real(dp), allocatable, intent(out), optional :: trajectory(:,:)

    type(directed_stored_graph) :: chain
    type(step_operator) :: statement
    type(set_graph) :: state_domain
    integer         :: n_state_domain
    real(dp), allocatable :: s(:), qold(:), qolder(:), zeros(:)
    real(dp) :: answered, h_edge, h_previous
    integer :: e, ncomp

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
    ! history states, so its first edge is a backward-euler step.
    ! All coefficients come from set_bdf; none are written here.
    statement = bdf(1, action, this % step)

    h_previous = this % step

    ! the unknown of each implicit solve is the state, so the
    ! solve is attached on the action's domain, read once here
    call read_state_domain(action, on, q, state_domain, n_state_domain, ncomp)

    allocate(zeros(size(q)))
    zeros = 0.0_dp

    qold = q

    do e = 1, chain % num_edges()

       h_edge = edge_step(this, steps, e)

       if (this % rule == MARCH_BDF2 .and. e > 1) then
          call statement % set_bdf(2, [h_edge, h_previous])
          statement % qolder = qolder
       else
          call statement % set_bdf(1, [h_edge])
       end if

       statement % qold = qold

       ! ncomp is the state's component count, so a multi-component
       ! entry is solved whole
       call this % inner % attach(statement, on, state_domain, &
            & n_state_domain, ncomp = ncomp)
       call this % inner % solve(zeros, q, answered)

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
  ! Reverse traversal. Under MARCH_FORWARD the caller's transposed
  ! statement is applied edge by edge in reverse order. Under the
  ! implicit rules the reverse traversal must linearize the solves
  ! actually taken, so it needs the action and the forward
  ! trajectory; calling without them stops the program, because an
  ! explicit reverse pass over an implicit forward march would
  ! return sensitivities of a different scheme.
  !===================================================================!

  subroutine march_adjoint(this, transposed, on, lambda, nsteps, steps, &
       & action, trajectory, seeds)

    class(marcher), intent(inout)      :: this
    class(graph_operation), intent(in) :: transposed
    class(directed_graph), intent(in)           :: on
    real(dp), intent(inout)            :: lambda(:)
    integer, intent(in)                :: nsteps
    real(dp), intent(in), optional     :: steps(:)
    class(graph_operation), intent(in), optional :: action
    real(dp), intent(in), optional     :: trajectory(:,:)
    real(dp), intent(in), optional     :: seeds(:,:)

    type(directed_stored_graph) :: chain
    real(dp), allocatable :: s(:)
    integer :: e

    call require_valid_steps(steps, nsteps)

    call this % instants(nsteps, chain)

    if (this % rule == MARCH_FORWARD) then

       do e = chain % num_edges(), 1, -1
          call read_statement(transposed, on, lambda, s)
          lambda = lambda - edge_step(this, steps, e) * s
       end do
       return

    end if

    if (.not. (present(action) .and. present(trajectory))) then
       error stop 'march_adjoint: the implicit reverse traversal needs the &
            &action and its forward trajectory'
    end if

    call substitute_backward(this, action, on, lambda, chain, steps, &
         & trajectory, seeds)

  end subroutine march_adjoint

  !===================================================================!
  ! Backward substitution over the chain: at each edge, in reverse
  ! order, solve
  !
  !      (a0 I + h_e S'(q_e))^T lambda_e = seed_e,
  !
  ! with the Jacobian assembled at the recorded state through
  ! tangent_of and dense_matrix_of and the transpose solved by the
  ! dense direct minimizer. The couplings to earlier instants are
  ! the constant coefficients a1 and a2, so the seed carried to
  ! edge e-1 is -a1 lambda_e plus the -a2 term from the edge after
  ! it. On return, lambda holds the sensitivity at the first
  ! instant.
  !
  ! The trajectory must hold one state per instant, and seeds, if
  ! present, one entry per instant; both are checked and stop the
  ! program, because a misaligned array would pair states with the
  ! wrong instants.
  !===================================================================!

  subroutine substitute_backward(this, action, on, lambda, chain, steps, &
       & trajectory, seeds)

    class(marcher), intent(inout)      :: this
    class(graph_operation), intent(in) :: action
    class(directed_graph), intent(in)  :: on
    real(dp), intent(inout)            :: lambda(:)
    type(directed_stored_graph), intent(in) :: chain
    real(dp), intent(in), optional     :: steps(:)
    real(dp), intent(in)               :: trajectory(:,:)
    real(dp), intent(in), optional     :: seeds(:,:)

    class(linearization_operator), allocatable :: tangent
    type(step_operator) :: table
    real(dp), allocatable :: jac(:,:), jstep(:,:)
    real(dp), allocatable :: seed(:), lambda_e(:), carry_one(:), carry_two(:)
    real(dp) :: answered, h_edge, h_previous
    integer :: e, n, j

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

    tangent = tangent_of(action)

    allocate(carry_one(n), carry_two(n))
    carry_one = 0.0_dp
    carry_two = 0.0_dp

    seed = lambda

    do e = chain % num_edges(), 1, -1

       h_edge = edge_step(this, steps, e)

       ! same coefficients as the forward march
       if (this % rule == MARCH_BDF2 .and. e > 1) then
          h_previous = edge_step(this, steps, e - 1)
          call table % set_bdf(2, [h_edge, h_previous])
       else
          call table % set_bdf(1, [h_edge])
       end if

       ! seeds(:, k), if given, is added when instant k is reached
       if (e < chain % num_edges()) then
          seed = carry_one
          if (present(seeds)) seed = seed + seeds(:, e + 1)
       end if

       call tangent % freeze(trajectory(:, e + 1))
       call dense_matrix_of(tangent, on, n, jac)

       jstep = table % hs * jac
       do j = 1, n
          jstep(j, j) = jstep(j, j) + table % a0
       end do

       call solve_dense_matrix_with_dense_direct(transpose(jstep), seed, &
            & this % singular_tolerance, lambda_e, answered)

       ! carry the couplings to the two earlier instants
       carry_one = carry_two - table % a1 * lambda_e
       carry_two = -table % a2 * lambda_e

    end do

    lambda = carry_one
    if (present(seeds)) lambda = lambda + seeds(:, 1)

  end subroutine substitute_backward

  !===================================================================!
  ! Read the action's domain and check the state fits it: the
  ! domain must be nonempty and size(q) must be a whole multiple
  ! of its entry count; both violations stop the program, because
  ! ncomp is derived from that division.
  !===================================================================!

  subroutine read_state_domain(action, on, q, state_domain, n_state_domain, &
       & ncomp)

    class(graph_operation), intent(in)  :: action
    class(directed_graph)          , intent(in)  :: on
    real(dp)              , intent(in)  :: q(:)
    type(set_graph)       , intent(out) :: state_domain
    integer               , intent(out) :: n_state_domain
    integer               , intent(out) :: ncomp

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

    ncomp = size(q) / n

  end subroutine read_state_domain

  !===================================================================!
  ! Evaluate the action at the state q. The result must live on
  ! the same domain as the state, not merely have the same length:
  ! a different domain of equal size would otherwise pass, and
  ! q - h s would mix fields from two domains. Both mismatches
  ! stop the program.
  !===================================================================!

  subroutine read_statement(action, on, q, s)

    class(graph_operation), intent(in) :: action
    class(directed_graph), intent(in)           :: on
    real(dp), intent(in)               :: q(:)
    real(dp), allocatable, intent(out) :: s(:)

    type(field)   :: state
    class(graph_field), allocatable :: output
    type(set_graph) :: state_domain, given
    integer         :: n_state_domain, n_given
    integer :: ncomp

    call read_state_domain(action, on, q, state_domain, n_state_domain, ncomp)

    state = field('state', state_domain, n_state_domain, ncomp=ncomp)
    call state % set_real_vector(q)

    call action % apply(on, [state], output)

    given   = output % domain()
    n_given = output % num_entries()
    if (.not. given % same_as(state_domain)) then
       error stop 'marcher: the action result lives on the state''s domain'
    end if

    call output % get_real_vector(s)
    if (size(s) /= size(q)) then
       error stop 'marcher: the action result matches the state''s width'
    end if

  end subroutine read_statement

  !===================================================================!
  ! Forward directional derivatives of any order along the
  ! recorded trajectory. At each edge the derivative of the step
  ! equation splits into
  !
  !      J_e q_e^(s) = -( a1 qprev^(s) + a2 qolder^(s)
  !                       + h_e * (degree-s composition with the
  !                         order-s state derivative set to zero) ),
  !      J_e = a0 I + h_e S'(q_e),
  !
  ! because the couplings to earlier instants are the constant
  ! coefficients a1 and a2 while everything inside S is assembled
  ! by chain_rule. One Jacobian per edge serves every order; only
  ! the right-hand side changes with s. Under MARCH_FORWARD the
  ! same recursion needs no solve.
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
    class(differentiable_operation), intent(in) :: action
    class(directed_graph), intent(in)           :: on
    integer, intent(in)                         :: nsteps
    real(dp), intent(in)                        :: trajectory(:,:)
    integer, intent(in)                         :: order
    real(dp), allocatable, intent(out)          :: sensitivities(:,:,:)
    real(dp), intent(in), optional              :: steps(:)
    type(field), intent(in), optional           :: parameters(:)
    type(argument_path), intent(in), optional   :: paths(:)

    type(directed_stored_graph) :: chain
    type(chain_rule)            :: composer
    type(step_operator)         :: table
    class(linearization_operator), allocatable :: tangent
    type(argument_path), allocatable :: assembled(:)
    type(field), allocatable         :: inputs(:)
    class(graph_field), allocatable  :: total_field
    type(field)     :: state
    type(set_graph) :: state_domain
    real(dp), allocatable :: jac(:,:), jstep(:,:), total(:), b(:), q_s(:)
    real(dp) :: answered, h_edge, h_previous
    integer :: n_state_domain, ncomp
    integer :: e, s_order, k, j, n, at, npaths

    if (order < 1) then
       error stop 'march_directional: the order is positive'
    end if

    call require_valid_steps(steps, nsteps)
    call this % instants(nsteps, chain)

    n = size(trajectory, 1)
    if (size(trajectory, 2) /= chain % num_vertices()) then
       error stop 'march_directional: the trajectory carries one state per instant'
    end if

    call read_state_domain(action, on, trajectory(:, 1), state_domain, &
         & n_state_domain, ncomp)

    ! check the parameter paths: each must name a slot >= 2 that a
    ! supplied parameter field covers
    npaths = 0
    if (present(paths)) then
       npaths = size(paths)
       do k = 1, npaths
          if (paths(k) % slot <= 1) then
             error stop 'march_directional: the state path is computed, &
                  &not supplied'
          end if
          if (.not. present(parameters)) then
             error stop 'march_directional: a parameter path names a supplied slot'
          end if
          if (paths(k) % slot > 1 + size(parameters)) then
             error stop 'march_directional: a parameter path names a supplied slot'
          end if
       end do
    end if

    allocate(sensitivities(n, order, nsteps + 1))
    sensitivities = 0.0_dp

    tangent = tangent_of(action)

    ! for each edge, assemble the composition degree by degree and
    ! advance each order's derivative
    do e = 1, chain % num_edges()

       h_edge = edge_step(this, steps, e)

       if (this % rule == MARCH_BDF2 .and. e > 1) then
          h_previous = edge_step(this, steps, e - 1)
          call table % set_bdf(2, [h_edge, h_previous])
       else
          call table % set_bdf(1, [h_edge])
       end if

       ! the composition is evaluated at the solved instant under
       ! an implicit rule and at the previous instant under the
       ! explicit rule, matching the state each rule differentiates
       if (this % rule == MARCH_FORWARD) then
          at = e
       else
          at = e + 1
       end if

       state = field('state', state_domain, n_state_domain, ncomp=ncomp)
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

       if (this % rule /= MARCH_FORWARD) then
          call tangent % freeze(trajectory(:, e + 1))
          call dense_matrix_of(tangent, on, n, jac)
          jstep = h_edge * jac
          do j = 1, n
             jstep(j, j) = jstep(j, j) + table % a0
          end do
       end if

       do s_order = 1, order

          call build_paths(this, sensitivities, state_domain, &
               & n_state_domain, ncomp, at, s_order, npaths, paths, assembled)

          call composer % assemble(action, on, inputs, s_order, assembled, &
               & total_field)
          call total_field % get_real_vector(total)

          if (this % rule == MARCH_FORWARD) then

             sensitivities(:, s_order, e + 1) = &
                  & sensitivities(:, s_order, e) - h_edge * total

          else

             b = table % a1 * sensitivities(:, s_order, e) + h_edge * total
             if (this % rule == MARCH_BDF2 .and. e > 1) then
                b = b + table % a2 * sensitivities(:, s_order, e - 1)
             end if

             call solve_dense_matrix_with_dense_direct(jstep, -b, &
                  & this % singular_tolerance, q_s, answered)
             sensitivities(:, s_order, e + 1) = q_s

          end if

       end do

    end do

  end subroutine march_directional

  !===================================================================!
  ! Build the argument paths for one edge and one order: the state
  ! path holds the solved derivatives below the current order;
  ! under an implicit rule its order-s entry is set to zero, which
  ! makes the assembled total the right-hand side for the unknown
  ! q^(s). The caller's parameter paths are appended unchanged.
  !===================================================================!

  subroutine build_paths(this, sensitivities, state_domain, n_state_domain, &
       & ncomp, at, s_order, npaths, parameter_paths, assembled)

    class(marcher), intent(in) :: this
    real(dp)      , intent(in) :: sensitivities(:,:,:)
    type(set_graph), intent(in) :: state_domain
    integer       , intent(in) :: n_state_domain, ncomp, at, s_order, npaths
    type(argument_path), intent(in), optional :: parameter_paths(:)
    type(argument_path), allocatable, intent(out) :: assembled(:)

    type(field) :: derivative_field
    real(dp), allocatable :: derivative_values(:)
    integer :: k

    allocate(assembled(1 + npaths))

    assembled(1) % slot = 1
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

       derivative_field = field('state path', state_domain, n_state_domain, &
            & ncomp=ncomp)
       call derivative_field % set_real_vector(derivative_values)
       assembled(1) % derivative(k) % occupied  = .true.
       assembled(1) % derivative(k) % direction = derivative_field
       deallocate(derivative_values)

    end do

    do k = 1, npaths
       assembled(1 + k) = parameter_paths(k)
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
    class(graph_operation), intent(in) :: action
    class(directed_graph), intent(in)  :: on
    real(dp), intent(inout)            :: q(:)
    real(dp), intent(in)               :: duration
    class(step_policy), intent(inout)  :: policy
    integer, intent(in)                :: max_attempts
    real(dp), allocatable, intent(out) :: steps_taken(:)
    logical, intent(out)               :: completed

    type(step_operator) :: statement
    type(set_graph)     :: state_domain
    integer             :: n_state_domain, ncomp
    real(dp), allocatable :: trial(:), predictor(:), qprev(:), s(:), zeros(:)
    real(dp), allocatable :: grown(:)
    real(dp) :: t, h, h_previous, estimate, answered
    integer  :: attempt, taken
    logical  :: accepted, have_previous

    if (duration <= 0.0_dp) then
       error stop 'march_adaptive: the duration is positive'
    end if
    if (max_attempts < 1) then
       error stop 'march_adaptive: the attempt budget is positive'
    end if

    call read_state_domain(action, on, q, state_domain, n_state_domain, ncomp)

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

             if (this % rule == MARCH_BDF2 .and. have_previous) then
                call statement % set_bdf(2, [h, h_previous])
                statement % qolder = qprev
             else
                call statement % set_bdf(1, [h])
             end if
             statement % qold = q

             trial = q
             call this % inner % attach(statement, on, state_domain, &
                  & n_state_domain, ncomp = ncomp)
             call this % inner % solve(zeros, trial, answered)

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

end module class_graph_marcher
