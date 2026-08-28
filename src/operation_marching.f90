!=====================================================================!
! Time integration over the reads graph: instants are vertices, and
! an edge j -> k says the statement determining instant k reads
! instant j, its slot - which history argument - carried on the
! edge; the first instant is read and determines nothing. The steps
! a statement takes are differences of the instants' times, and the
! graph is the only thing that says which instant comes before
! which: a march is the loop over the graph in its forward
! orientation, the adjoint the loop in reverse, and no routine here
! writes an instant's number down. march takes an optional per-step
! array of step sizes; absent, every step is the marcher's one step.
!
! Every rule is one scheme - the theta residual of operation_step -
! and every verb is one traversal of the chain with that scheme:
!
!      MARCH_FORWARD    theta 0, reach 1     q_n = q_(n-1) - h S(q_(n-1))
!      MARCH_BACKWARD   theta 1, reach 1     q_n - q_(n-1) + h S(q_n) = 0
!      MARCH_BDF2       theta 1, reach 2     (3 q_n - 4 q_(n-1) + q_(n-2))/2
!                                              + h S(q_n) = 0, the first
!                                              edge backward euler
!
! The scheme's inputs for an edge are its argument tuple: the state,
! every auxiliary argument of the action as a supplied parameter
! field, then the history the statement reaches. A march over an
! action with auxiliaries is refused unless all of them are supplied,
! because the history is placed after them. With theta = 0 the state
! block of the residual is a0 I, so the step is the residual at the
! zero state divided by -a0 and no minimizer runs; otherwise the held
! minimizer (inner) solves the step with the rest of the tuple held.
!
! march_adjoint loops over the graph in reverse over the recorded
! trajectory and applies the dual of every local differential: at
! each instant the state block of the tangent is compiled, transposed
! and solved, and the dual of each history block - also compiled and
! transposed - is subtracted from the seed of the instant it reads,
! wherever that instant lies.
! Nothing here names a1 or a2. march_directional computes forward
! directional derivatives of any order: the chain rule over the
! scheme with the state, history and parameter paths gives the
! right-hand side, and the state block is solved. march_adaptive
! chooses the steps at run time under a step_policy.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_marching

  use util_precision  , only : dp
  use operation_action, only : operation, argument
  use view_directed, only : directed_graph, forward, reverse
  use field_calculus, only : field
  use graph_fractal      , only : graph
  use field_stored  , only : stored_field
  use view_directed_stored        , only : stored_directed_graph
  use operation_step   , only : scheme, bdf
  use operation_minimization , only : minimizer
  use operation_linearization, only : linearization, tangent_of
  use operation_chain_rule, only : chain_rule, argument_path
  use operation_solve_connectivity, only : assembly_connectivity
  use operation_solve_connectivity, only : tangent_residual_connectivity
  use operation_solve_connectivity, only : tangent_jacobian_connectivity
  use operation_solve_connectivity, only : adjoint_residual_connectivity
  use operation_solve_connectivity, only : adjoint_jacobian_connectivity
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

  subroutine instants(this, nsteps, steps, chain, slot, time)

    class(marcher), intent(in) :: this
    integer, intent(in) :: nsteps
    real(dp), intent(in), optional :: steps(:)
    type(stored_directed_graph), intent(out) :: chain
    integer , allocatable, intent(out) :: slot(:)
    real(dp), allocatable, intent(out) :: time(:)

    integer, allocatable :: tails(:), heads(:)
    integer :: k, e

    ! every instant past the first reads the one before it; under
    ! the two-step rule, the one before that as well, once there is one
    allocate(tails(2 * nsteps), heads(2 * nsteps), slot(2 * nsteps))
    e = 0
    do k = 2, nsteps + 1
       e = e + 1
       tails(e) = k - 1
       heads(e) = k
       slot(e)  = 1
       if (this % rule == MARCH_BDF2 .and. k >= 3) then
          e = e + 1
          tails(e) = k - 2
          heads(e) = k
          slot(e)  = 2
       end if
    end do
    chain = stored_directed_graph(nsteps + 1, tails=tails(1:e), heads=heads(1:e))
    slot  = slot(1:e)

    allocate(time(nsteps + 1))
    time(1) = 0.0_dp
    do k = 1, nsteps
       time(k + 1) = time(k) + edge_step(this, steps, k)
    end do

  end subroutine instants

  !===================================================================!
  ! What the statement determining an instant reads, by slot: the
  ! instant at each history slot, zero where the slot is not read.
  ! An instant reading nothing is given, not determined.
  !===================================================================!

  subroutine reads_of(chain, slot, k, read)

    type(stored_directed_graph), intent(in)  :: chain
    integer                    , intent(in)  :: slot(:), k
    integer                    , intent(out) :: read(2)

    integer, allocatable :: edges(:)
    integer :: e

    read = 0
    call chain % incoming_edges(k, edges)
    do e = 1, size(edges)
       read(slot(edges(e))) = chain % edge_tail(edges(e))
    end do

  end subroutine reads_of

  !===================================================================!
  ! The statement configured for the instant it determines: its step
  ! is the time from the instant it reads first, and the step before
  ! that, where a second read exists, the time between its two reads.
  !===================================================================!

  subroutine configured(this, statement, time, k, read)

    class(marcher), intent(in)  :: this
    type(scheme), intent(inout) :: statement
    real(dp)    , intent(in)    :: time(:)
    integer     , intent(in)    :: k, read(2)

    if (read(2) > 0) then
       call configure_edge(this, statement, time(k) - time(read(1)), &
            & time(read(1)) - time(read(2)), .true.)
    else
       call configure_edge(this, statement, time(k) - time(read(1)), &
            & time(k) - time(read(1)), .false.)
    end if

  end subroutine configured

  !===================================================================!
  ! Integrate nsteps edges forward under the marcher's rule,
  ! optionally recording the state at every instant.
  !===================================================================!

  subroutine march(this, action, on, q, nsteps, steps, trajectory, parameters)

    class(marcher), intent(inout)      :: this
    class(operation), intent(in) :: action
    class(directed_graph), intent(in)           :: on
    real(dp), intent(inout)            :: q(:)
    integer, intent(in)                :: nsteps
    real(dp), intent(in), optional     :: steps(:)
    real(dp), allocatable, intent(out), optional :: trajectory(:,:)
    type(stored_field), intent(in), optional :: parameters(:)

    type(stored_directed_graph) :: chain
    type(scheme) :: statement
    type(graph) :: state_domain
    type(stored_field), allocatable :: inputs(:)
    integer         :: n_state_domain
    real(dp), allocatable :: states(:,:), zeros(:), time(:)
    integer , allocatable :: order(:), slot(:)
    real(dp) :: achieved
    integer :: i, k, read(2), num_components

    call require_valid_steps(steps, nsteps)
    call require_parameters(action, parameters)
    call this % instants(nsteps, steps, chain, slot, time)

    statement = prepared_statement(this, action)

    ! the unknown of each step is the state, so the solve is attached
    ! on the action's domain, read once here
    call read_state_domain(action, on, q, state_domain, n_state_domain, num_components)

    allocate(zeros(size(q)))
    zeros = 0.0_dp
    allocate(states(size(q), chain % num_vertices()))

    ! the loop over the reads graph: an instant reading nothing is the
    ! one given; every other is determined from the instants it reads
    order = chain % loop(forward)
    do i = 1, size(order)
       k = order(i)
       call reads_of(chain, slot, k, read)
       if (read(1) == 0) then
          states(:, k) = q
          cycle
       end if
       call configured(this, statement, time, k, read)
       if (read(2) > 0) then
          call edge_inputs(statement, q, states(:, read(1)), state_domain, n_state_domain, &
               & num_components, inputs, qolder=states(:, read(2)), parameters=parameters)
       else
          call edge_inputs(statement, q, states(:, read(1)), state_domain, n_state_domain, &
               & num_components, inputs, parameters=parameters)
       end if
       call advance(this, statement, on, state_domain, n_state_domain, &
            & num_components, inputs, zeros, q, achieved)
       states(:, k) = q
    end do

    if (present(trajectory)) trajectory = states

  end subroutine march

  !===================================================================!
  ! One step of the statement from its input tuple: with theta = 0
  ! the state block is a0 I, so the new state is the residual at the
  ! zero state over -a0; otherwise the held minimizer solves the step
  ! with the rest of the tuple held fixed.
  !===================================================================!

  subroutine advance(this, statement, on, state_domain, n_state_domain, &
       & num_components, inputs, zeros, q, achieved)

    class(marcher), intent(inout)     :: this
    type(scheme), intent(in)          :: statement
    class(directed_graph), intent(in) :: on
    type(graph), intent(in)           :: state_domain
    integer, intent(in)               :: n_state_domain, num_components
    type(stored_field), intent(inout) :: inputs(:)
    real(dp), intent(in)              :: zeros(:)
    real(dp), intent(inout)           :: q(:)
    real(dp), intent(out)             :: achieved

    class(field), allocatable :: residual
    real(dp), allocatable :: c(:)

    if (statement % theta == 0.0_dp) then
       call inputs(1) % set_real_vector(zeros)
       call statement % apply(on, inputs, residual)
       call residual % real_vector(c)
       q = -c / statement % a0
       achieved = 0.0_dp
       return
    end if

    ! the step statement is affine in the state, A q + c = 0 with c
    ! read at the zero state; a linear solver answers the linear part
    ! alone, so it is asked for A q = -c
    call this % inner % attach(statement, on, state_domain, &
         & n_state_domain, num_components = num_components, &
         & held_inputs = inputs(2:))
    call this % inner % constant(c)
    call this % inner % solve(-c, q, achieved)

  end subroutine advance

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

  !===================================================================!
  ! Every auxiliary argument of the action must be supplied as a
  ! parameter field, one per argument in the action's order, because
  ! the scheme places the history after them; a missing or surplus
  ! parameter stops the program.
  !===================================================================!

  subroutine require_parameters(action, parameters)

    class(operation), intent(in) :: action
    type(stored_field), intent(in), optional :: parameters(:)

    integer :: np

    np = 0
    if (present(parameters)) np = size(parameters)

    if (np /= action % num_arguments() - 1) then
       error stop 'marcher: every auxiliary argument of the action is supplied as a parameter'
    end if

  end subroutine require_parameters

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
  ! The scheme of the marcher's rule over the action: order 1 at the
  ! marcher's step, theta 0 for the explicit rule and 1 otherwise.
  !===================================================================!

  function prepared_statement(this, action) result(statement)

    class(marcher), intent(in)   :: this
    class(operation), intent(in) :: action
    type(scheme) :: statement

    statement = bdf(1, action, this % step)

    if (this % rule == MARCH_FORWARD) then
       statement % theta = 0.0_dp
    else
       statement % theta = 1.0_dp
    end if

  end function prepared_statement

  !===================================================================!
  ! The scheme's coefficients for one edge: order 2 with the
  ! previous step when the rule is MARCH_BDF2 and an older state
  ! exists, order 1 otherwise. Every coefficient comes from set_bdf.
  !===================================================================!

  subroutine configure_edge(this, statement, h, h_previous, have_older)

    class(marcher), intent(in)  :: this
    type(scheme), intent(inout) :: statement
    real(dp), intent(in)        :: h
    real(dp), intent(in)        :: h_previous
    logical, intent(in)         :: have_older

    if (this % rule == MARCH_BDF2 .and. have_older) then
       call statement % set_bdf(2, [h, h_previous])
    else
       call statement % set_bdf(1, [h])
    end if

  end subroutine configure_edge


  !===================================================================!
  ! The statement's input tuple for one edge, in the scheme's argument
  ! order: the state, the parameter fields - one per auxiliary
  ! argument of the action - and the history states the statement
  ! reaches. This tuple is what the minimizer holds during the solve
  ! and what a tangent is frozen at, so residual and tangent see one
  ! function. A parameter count other than the action's auxiliaries,
  ! or a reach-2 statement without its second history state, stops
  ! the program.
  !===================================================================!

  subroutine edge_inputs(statement, state_values, qold, state_domain, &
       & n_state_domain, num_components, inputs, qolder, parameters)

    type(scheme), intent(in) :: statement
    real(dp)    , intent(in) :: state_values(:)
    real(dp)    , intent(in) :: qold(:)
    type(graph) , intent(in) :: state_domain
    integer     , intent(in) :: n_state_domain, num_components
    type(stored_field), allocatable, intent(out) :: inputs(:)
    real(dp)    , intent(in), optional :: qolder(:)
    type(stored_field), intent(in), optional :: parameters(:)

    integer :: m, np, j

    m  = statement % action % num_arguments()
    np = 0
    if (present(parameters)) np = size(parameters)
    if (np /= m - 1) then
       error stop 'marcher: every auxiliary argument of the action is supplied as a parameter'
    end if

    allocate(inputs(m + statement % reach))

    inputs(1) = stored_field('state', state_domain, n_state_domain, num_components=num_components)
    call inputs(1) % set_real_vector(state_values)

    do j = 1, m - 1
       inputs(1 + j) = parameters(j)
    end do

    inputs(m + 1) = stored_field('history 1', state_domain, n_state_domain, num_components=num_components)
    call inputs(m + 1) % set_real_vector(qold)

    if (statement % reach >= 2) then
       if (.not. present(qolder)) then
          error stop 'marcher: a reach-2 statement is given two history states'
       end if
       inputs(m + 2) = stored_field('history 2', state_domain, n_state_domain, num_components=num_components)
       call inputs(m + 2) % set_real_vector(qolder)
    end if

  end subroutine edge_inputs

  !===================================================================!
  ! The input tuple at a recorded instant: its state, and the states
  ! at the instants it reads.
  !===================================================================!

  subroutine recorded_inputs(statement, k, read, trajectory, state_domain, &
       & n_state_domain, num_components, inputs, parameters)

    type(scheme), intent(in) :: statement
    integer     , intent(in) :: k, read(2)
    real(dp)    , intent(in) :: trajectory(:,:)
    type(graph) , intent(in) :: state_domain
    integer     , intent(in) :: n_state_domain, num_components
    type(stored_field), allocatable, intent(out) :: inputs(:)
    type(stored_field), intent(in), optional :: parameters(:)

    if (read(2) > 0) then
       call edge_inputs(statement, trajectory(:, k), trajectory(:, read(1)), &
            & state_domain, n_state_domain, num_components, inputs, &
            & qolder=trajectory(:, read(2)), parameters=parameters)
    else
       call edge_inputs(statement, trajectory(:, k), trajectory(:, read(1)), &
            & state_domain, n_state_domain, num_components, inputs, &
            & parameters=parameters)
    end if

  end subroutine recorded_inputs

  !===================================================================!
  ! The loop over the reads graph in reverse over the recorded
  ! trajectory, applying the dual of every local differential. At
  ! each determined instant k, with A the state block of the
  ! residual's tangent at its input tuple,
  !
  !      A^T lambda_k = seed_k
  !      seed of the instant read at slot s  -=  (D_history(s) R_k)^T lambda_k
  !
  ! for every slot the statement reads; each block is compiled to a
  ! stencil and transposed. With theta = 0 the state block is a0 I and
  ! no solve runs. The seed of the instant the loop begins at is
  ! lambda as given; seeds(:, k), if present, is what every other
  ! instant's seed begins as. On return, lambda holds the sensitivity
  ! at the given instant, the one reading nothing. The trajectory must
  ! hold one state per instant and seeds one entry per instant; both
  ! are checked and stop the program, because a misaligned array
  ! would pair states with the wrong instants.
  !===================================================================!

  subroutine march_adjoint(this, action, on, lambda, nsteps, trajectory, &
       & steps, seeds, parameters)

    class(marcher), intent(inout)      :: this
    class(operation), intent(in)       :: action
    class(directed_graph), intent(in)  :: on
    real(dp), intent(inout)            :: lambda(:)
    integer, intent(in)                :: nsteps
    real(dp), intent(in)               :: trajectory(:,:)
    real(dp), intent(in), optional     :: steps(:)
    real(dp), intent(in), optional     :: seeds(:,:)
    type(stored_field), intent(in), optional :: parameters(:)

    type(stored_directed_graph) :: chain
    type(scheme)       :: statement
    type(assembly_connectivity) :: residual_connectivity, jacobian_connectivity
    type(stored_field), allocatable :: inputs(:)
    type(graph) :: state_domain
    real(dp), allocatable :: seed(:,:), lambda_k(:), g(:), time(:)
    integer , allocatable :: order(:), slot(:)
    real(dp) :: achieved
    integer :: i, k, s, n, history, read(2), n_state_domain, num_components

    call require_valid_steps(steps, nsteps)
    call require_parameters(action, parameters)
    call this % instants(nsteps, steps, chain, slot, time)

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

    call read_state_domain(action, on, lambda, state_domain, n_state_domain, num_components)
    statement = prepared_statement(this, action)

    allocate(seed(n, chain % num_vertices()), lambda_k(n))
    seed = 0.0_dp
    if (present(seeds)) seed = seeds

    order = chain % loop(reverse)
    seed(:, order(1)) = lambda

    do i = 1, size(order)
       k = order(i)
       call reads_of(chain, slot, k, read)

       ! the given instant: its sensitivity is what was seeded into it
       if (read(1) == 0) then
          lambda_k = seed(:, k)
          cycle
       end if

       call configured(this, statement, time, k, read)
       call recorded_inputs(statement, k, read, trajectory, state_domain, &
            & n_state_domain, num_components, inputs, parameters)
       residual_connectivity = adjoint_residual_connectivity(statement % reach)
       jacobian_connectivity = adjoint_jacobian_connectivity(size(inputs))

       ! The adjoint Jacobian connectivity names the state block as
       ! the differentiated slot and marks it as transposed. The
       ! residual connectivity names the history slots whose residual
       ! partials push seed backwards.
       call solve_adjoint_state_block(jacobian_connectivity, statement, &
            & inputs, on, n, seed(:, k), lambda_k, achieved)

       ! the dual of each history block, subtracted from the seed of
       ! the instant that slot reads. The slots come from the supplied
       ! connectivity, not from an inline reach loop.
       do s = 1, residual_connectivity % num_history()
          history = residual_connectivity % history(s)
          if (read(history) == 0) cycle
          call transposed_block(statement, residual_connectivity, history, &
               & inputs, on, n, lambda_k, g)
          seed(:, read(history)) = seed(:, read(history)) - g
       end do
    end do

    lambda = lambda_k

  end subroutine march_adjoint

  subroutine solve_adjoint_state_block(connectivity, statement, inputs, on, &
       & n, seed, lambda, achieved)

    type(assembly_connectivity), intent(in) :: connectivity
    type(scheme)              , intent(in) :: statement
    type(stored_field)        , intent(in) :: inputs(:)
    class(directed_graph)     , intent(in) :: on
    integer                   , intent(in) :: n
    real(dp)                  , intent(in) :: seed(:)
    real(dp), allocatable     , intent(out) :: lambda(:)
    real(dp)                  , intent(out) :: achieved

    type(dense_direct) :: direct
    type(linearization) :: tangent
    type(stencil) :: compiled, adjoint

    allocate(lambda(size(seed)))

    if (.not. connectivity % transposed()) then
       error stop 'march_adjoint: the state block connectivity is transposed'
    end if

    if (statement % theta == 0.0_dp) then
       lambda = seed / statement % a0
       achieved = 0.0_dp
       return
    end if

    tangent = tangent_of(statement, &
         & statement % argument(connectivity % differentiated()))
    call tangent % freeze(inputs)
    compiled = stencil(tangent, on, n)
    adjoint  = compiled % transpose()
    call direct % attach(adjoint, adjoint % pattern, &
         & adjoint % pattern % vertex_set(), &
         & adjoint % pattern % num_vertices())
    lambda = 0.0_dp
    call direct % solve(seed, lambda, achieved)

  end subroutine solve_adjoint_state_block

  !===================================================================!
  ! (D_a R)^T lambda for one argument a of the statement at its input
  ! tuple: the tangent in a, compiled to a stencil, transposed and
  ! applied. The block is square - a history state lives on the
  ! state's domain.
  !===================================================================!

  subroutine transposed_block(statement, connectivity, history, inputs, on, &
       & n, lambda, g)

    type(scheme), intent(in)          :: statement
    type(assembly_connectivity), intent(in) :: connectivity
    integer, intent(in)               :: history
    type(stored_field), intent(in)    :: inputs(:)
    class(directed_graph), intent(in) :: on
    integer, intent(in)               :: n
    real(dp), intent(in)              :: lambda(:)
    real(dp), allocatable, intent(out) :: g(:)

    type(linearization) :: tangent
    type(stencil) :: compiled, adjoint
    type(stored_field) :: lambda_field
    class(field), allocatable :: answer

    if (.not. connectivity % transposed()) then
       error stop 'march_adjoint: history block connectivity is transposed'
    end if

    tangent = tangent_of(statement, statement % history(history), &
         & at_inputs=inputs)
    compiled = stencil(tangent, on, n)
    adjoint  = compiled % transpose()

    lambda_field = stored_field('lambda', inputs(1) % domain(), inputs(1) % num_entries(), &
         & num_components=inputs(1) % num_components())
    call lambda_field % set_real_vector(lambda)
    call adjoint % apply(adjoint % pattern, [lambda_field], answer)
    call answer % real_vector(g)

  end subroutine transposed_block

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
  ! Forward directional derivatives of any order along the recorded
  ! trajectory. At edge e the derivative of the step equation is
  !
  !      A q_(e+1)^(s) = -( chain rule total over the scheme with the
  !                         order-s state derivative zero, the history
  !                         paths known, the parameter paths supplied )
  !
  ! with A the state block of the residual's tangent at the edge's
  ! input tuple; with theta = 0 A is a0 I and no solve runs. One
  ! tangent per edge is attached to the dense direct minimizer and
  ! solved once per order; only the right-hand side changes with s.
  !
  ! Parameter paths name auxiliary arguments of the action, each
  ! covered by a supplied parameter field. A path on the state stops
  ! the program, because the state's path is what this routine
  ! computes; a path on an uncovered argument also stops the program.
  ! The initial state is held fixed, so its derivatives are zero.
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
    type(assembly_connectivity) :: residual_connectivity, jacobian_connectivity
    type(argument_path), allocatable :: assembled(:)
    type(stored_field), allocatable  :: inputs(:)
    class(field), allocatable  :: total_field
    type(graph) :: state_domain
    real(dp), allocatable :: total(:), q_s(:), time(:)
    integer , allocatable :: walk(:), slot(:)
    real(dp) :: achieved
    integer :: n_state_domain, num_components
    integer :: i, s_order, k, n, npaths, read(2)

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
    call require_parameters(action, parameters)
    call this % instants(nsteps, steps, chain, slot, time)

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
          if (.not. covered_by_parameters(action, paths(k) % wrt, &
               & action % num_arguments() - 1)) then
             error stop 'march_directional: a parameter path names a supplied argument'
          end if
       end do
    end if

    allocate(sensitivities(n, order, nsteps + 1))
    sensitivities = 0.0_dp

    statement = prepared_statement(this, action)
    allocate(q_s(n))

    ! for each determined instant, in the loop's order, assemble the
    ! composition degree by degree and advance each order's derivative
    walk = chain % loop(forward)
    do i = 1, size(walk)
       k = walk(i)
       call reads_of(chain, slot, k, read)
       if (read(1) == 0) cycle

       call configured(this, statement, time, k, read)
       call recorded_inputs(statement, k, read, trajectory, state_domain, &
            & n_state_domain, num_components, inputs, parameters)
       residual_connectivity = &
            & tangent_residual_connectivity(order, statement % reach, npaths)
       jacobian_connectivity = tangent_jacobian_connectivity(size(inputs))

       ! The tangent Jacobian connectivity supplies the diagonal state
       ! block; the residual connectivity supplies the chain-rule path
       ! tuple for the right hand side below.
       call attach_tangent_state_block(jacobian_connectivity, statement, &
            & inputs, on, state_domain, n_state_domain, num_components, direct)

       do s_order = 1, order
          call build_paths(statement, residual_connectivity, sensitivities, &
               & state_domain, n_state_domain, num_components, k, read, &
               & s_order, npaths, paths, assembled)
          call composer % assemble(statement, on, inputs, s_order, &
               & assembled, total_field, &
               & connectivity=residual_connectivity % partial(s_order))
          call total_field % real_vector(total)
          if (statement % theta == 0.0_dp) then
             q_s = -total / statement % a0
          else
             q_s = 0.0_dp
             call direct % solve(-total, q_s, achieved)
          end if
          sensitivities(:, s_order, k) = q_s
       end do
    end do

  end subroutine march_directional

  subroutine attach_tangent_state_block(connectivity, statement, inputs, on, &
       & state_domain, n_state_domain, num_components, direct)

    type(assembly_connectivity), intent(in) :: connectivity
    type(scheme)              , intent(in) :: statement
    type(stored_field)        , intent(in) :: inputs(:)
    class(directed_graph)     , intent(in) :: on
    type(graph)               , intent(in) :: state_domain
    integer                   , intent(in) :: n_state_domain, num_components
    type(dense_direct)        , intent(inout) :: direct

    type(linearization) :: tangent

    if (connectivity % differentiated() /= connectivity % state()) then
       error stop 'march_directional: tangent connectivity differentiates the state'
    end if

    if (statement % theta == 0.0_dp) return

    tangent = tangent_of(statement, &
         & statement % argument(connectivity % differentiated()))
    call tangent % freeze(inputs)
    call direct % attach(tangent, on, state_domain, n_state_domain, &
         & num_components = num_components)

  end subroutine attach_tangent_state_block

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
  ! The argument paths of the scheme at one instant and one order:
  ! the supplied connectivity names the state, history, and parameter
  ! path slots. The state path holds the solved derivatives below the
  ! current order and zero at the order, which makes the assembled
  ! total the right-hand side for the unknown q^(s); each history path
  ! holds every derivative of the state at the instant its slot reads;
  ! the caller's parameter paths are restated in the scheme's argument
  ! space.
  !===================================================================!

  subroutine build_paths(statement, connectivity, sensitivities, &
       & state_domain, n_state_domain, num_components, at, read, s_order, &
       & npaths, parameter_paths, assembled)

    type(scheme)  , intent(in) :: statement
    type(assembly_connectivity), intent(in) :: connectivity
    real(dp)      , intent(in) :: sensitivities(:,:,:)
    type(graph), intent(in) :: state_domain
    integer       , intent(in) :: n_state_domain, num_components, at, read(2), s_order, npaths
    type(argument_path), intent(in), optional :: parameter_paths(:)
    type(argument_path), allocatable, intent(out) :: assembled(:)

    real(dp), allocatable :: zero(:)
    integer :: k, j, slot, history

    if (connectivity % num_inputs() /= 1 + statement % reach + npaths) then
       error stop 'march_directional: tangent connectivity matches the path tuple'
    end if
    if (connectivity % num_history() /= statement % reach) then
       error stop 'march_directional: tangent connectivity matches the history reach'
    end if
    if (connectivity % num_parameters() /= npaths) then
       error stop 'march_directional: tangent connectivity matches parameter paths'
    end if

    allocate(assembled(connectivity % num_inputs()))
    allocate(zero(size(sensitivities, 1)))
    zero = 0.0_dp

    ! the state path: the unknown order's entry is zero while the
    ! total is assembled
    slot = connectivity % state()
    assembled(slot) % wrt = statement % state()
    allocate(assembled(slot) % derivative(s_order))
    do k = 1, s_order - 1
       call occupy(assembled(slot), k, sensitivities(:, k, at))
    end do
    call occupy(assembled(slot), s_order, zero)

    ! the history paths: every derivative of the earlier instants
    do j = 1, connectivity % num_history()
       history = connectivity % history(j)
       slot = connectivity % input(1 + j)
       assembled(slot) % wrt = statement % history(history)
       allocate(assembled(slot) % derivative(s_order))
       do k = 1, s_order
          call occupy(assembled(slot), k, sensitivities(:, k, read(history)))
       end do
    end do

    do k = 1, npaths
       slot = connectivity % parameter(k)
       assembled(slot) = parameter_paths(k)
       assembled(slot) % wrt = &
            & statement % from_action(parameter_paths(k) % wrt)
    end do

  contains

    subroutine occupy(path, k, values)
      type(argument_path), intent(inout) :: path
      integer, intent(in) :: k
      real(dp), intent(in) :: values(:)
      path % derivative(k) % occupied  = .true.
      path % derivative(k) % direction = stored_field('path', state_domain, &
           & n_state_domain, num_components=num_components)
      call path % derivative(k) % direction % set_real_vector(values)
    end subroutine occupy

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
       & max_attempts, steps_taken, completed, parameters)

    class(marcher), intent(inout)      :: this
    class(operation), intent(in) :: action
    class(directed_graph), intent(in)  :: on
    real(dp), intent(inout)            :: q(:)
    real(dp), intent(in)               :: duration
    class(step_policy), intent(inout)  :: policy
    integer, intent(in)                :: max_attempts
    real(dp), allocatable, intent(out) :: steps_taken(:)
    logical, intent(out)               :: completed
    type(stored_field), intent(in), optional :: parameters(:)

    type(scheme) :: statement
    type(graph)     :: state_domain
    type(stored_field), allocatable :: inputs(:)
    integer             :: n_state_domain, num_components
    real(dp), allocatable :: trial(:), predictor(:), qprev(:), zeros(:)
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

    call require_parameters(action, parameters)
    call read_state_domain(action, on, q, state_domain, n_state_domain, num_components)

    statement = prepared_statement(this, action)
    allocate(zeros(size(q)))
    zeros = 0.0_dp

    allocate(steps_taken(0))
    t             = 0.0_dp
    taken         = 0
    have_previous = .false.
    h_previous    = 0.0_dp
    completed     = .false.

    do while (duration - t > spacing(duration))

       call policy % propose(h)

       accepted = .false.
       attempt  = 0

       do while (.not. accepted .and. attempt < max_attempts)

          attempt = attempt + 1

          if (h <= 0.0_dp) then
             error stop 'march_adaptive: the policy proposes a positive step'
          end if
          h = min(h, duration - t)

          ! the trial state from the accepted state; q is not modified.
          ! qprev is unallocated until the first accepted edge, so the
          ! statement takes order 1 until then
          call configure_edge(this, statement, h, h_previous, allocated(qprev))
          call edge_inputs(statement, q, q, state_domain, n_state_domain, &
               & num_components, inputs, qolder=qprev, parameters=parameters)

          trial = q
          call advance(this, statement, on, state_domain, n_state_domain, &
               & num_components, inputs, zeros, trial, achieved)

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
