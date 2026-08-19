!=====================================================================!
! GTI TIME ADAPTIVE CONTROLLER
!
! The policy layer's first inhabitant: given a solved time graph,
! one step forward is DECIDED here and only executed below -
!
!      policy proposes (step, order)
!        -> the controller builds the candidate and tries
!           reversible growth
!        -> the controller measures the evidence
!        -> the policy judges: accept keeps the grown graph,
!           reject discards it whole and proposes a retry
!
! - the loop the roadmap promised, arriving exactly as promised: a
! client of try_candidate, adding not one line to any existing
! driver. The mechanism/policy split is the whole architecture of
! this seat. The CONTROLLER owns the lifecycle and the
! measurement; the POLICY - an abstract verb-triple the caller
! extends - owns every number: the first proposal, the acceptance
! threshold, the retry step and order. A different adaptivity is a
! different policy, never a different controller.
!
! The candidate is built from the graph's own tail: a warm-started
! vertex at t_last + h, and a BDF relation minted by the
! variable-step builders - order 1 from [h], order 2 from
! [h, h_prev] with h_prev read off the graph's times, so
! nonuniform history is priced exactly, not assumed uniform.
!
! The evidence is measured, not judged, here: the error estimate
! is the distance between the accepted-candidate solve and the
! extrapolating predictor through the graph's tail - constant
! through one history vertex, linear through two. A non-converged
! attempt carries no estimate and is never kept; the growth seat
! has already rolled it back.
!
! An advance whose every attempt is rejected or unconverged is a
! LAWFUL answer, never an error stop: the budget is spent, the
! record says so, and the graph stands exactly as it was.
!
! The controller carries nothing: no graph, no forms, no solver
! state, no policy, no map - and it never estimates a step size
! of its own.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_time_adaptive_controllers

  use iso_fortran_env      , only : dp => REAL64
  use gti_value_buffers    , only : gti_value_buffer
  use gti_state_bundles    , only : GTI_STATE_Q
  use gti_design_bundles   , only : gti_design_bundle
  use gti_form_interface   , only : gti_differentiable_form
  use gti_time_graphs      , only : gti_time_graph
  use gti_time_motif_builders, only : gti_time_motif_builder
  use gti_time_forward_drivers, only : gti_time_forward_options
  use gti_time_adaptive_growth_drivers, only : gti_time_growth_candidate, &
       & gti_time_growth_step_result, gti_time_adaptive_growth_driver

  implicit none

  private
  public :: gti_time_step_proposal
  public :: gti_time_attempt_evidence
  public :: gti_time_step_policy
  public :: gti_time_halving_step_policy
  public :: gti_time_adaptive_controller_options
  public :: gti_time_adaptive_attempt
  public :: gti_time_adaptive_advance_result
  public :: gti_time_adaptive_controller

  !===================================================================!
  ! One proposed step: how far, and at which BDF order.
  !===================================================================!

  type :: gti_time_step_proposal

     real(dp) :: step  = 0.0_dp
     integer  :: order = 1

  end type gti_time_step_proposal

  !===================================================================!
  ! What one attempt showed: whether the local solve converged,
  ! and - when it did - how far the solution landed from the
  ! extrapolating predictor. Measured by the controller, judged by
  ! the policy.
  !===================================================================!

  type :: gti_time_attempt_evidence

     logical  :: converged      = .false.
     logical  :: has_estimate   = .false.
     real(dp) :: error_estimate = huge(1.0_dp)
     real(dp) :: step           = 0.0_dp
     integer  :: order          = 0
     integer  :: attempt        = 0

  end type gti_time_attempt_evidence

  !===================================================================!
  ! The abstract policy: three deferred verbs that own every
  ! number - the first proposal, the acceptance verdict, and the
  ! retry proposal. The controller never decides; a policy never
  ! executes.
  !===================================================================!

  type, abstract :: gti_time_step_policy

   contains

     procedure(policy_first_proposal), deferred :: first_proposal
     procedure(policy_judge)         , deferred :: judge
     procedure(policy_retry_proposal), deferred :: retry_proposal

  end type gti_time_step_policy

  abstract interface

     subroutine policy_first_proposal(this, history_count, proposal)
       import :: gti_time_step_policy, gti_time_step_proposal
       class(gti_time_step_policy)  , intent(inout) :: this
       integer                      , intent(in)    :: history_count
       type(gti_time_step_proposal) , intent(out)   :: proposal
     end subroutine policy_first_proposal

     subroutine policy_judge(this, evidence, accept)
       import :: gti_time_step_policy, gti_time_attempt_evidence
       class(gti_time_step_policy)    , intent(inout) :: this
       type(gti_time_attempt_evidence), intent(in)    :: evidence
       logical                        , intent(out)   :: accept
     end subroutine policy_judge

     subroutine policy_retry_proposal(this, evidence, proposal)
       import :: gti_time_step_policy, gti_time_step_proposal, &
            & gti_time_attempt_evidence
       class(gti_time_step_policy)    , intent(inout) :: this
       type(gti_time_attempt_evidence), intent(in)    :: evidence
       type(gti_time_step_proposal)   , intent(out)   :: proposal
     end subroutine policy_retry_proposal

  end interface

  !===================================================================!
  ! The reference policy: start at the caller's step, accept when
  ! the estimate is within the caller's tolerance, halve on every
  ! rejection, and prefer the caller's order once two history
  ! vertices exist. A policy example - not the policy.
  !===================================================================!

  type, extends(gti_time_step_policy) :: gti_time_halving_step_policy

     real(dp) :: initial_step    = 1.0_dp
     real(dp) :: error_tolerance = 1.0e-2_dp
     integer  :: preferred_order = 1

   contains

     procedure :: first_proposal => halving_first_proposal
     procedure :: judge          => halving_judge
     procedure :: retry_proposal => halving_retry_proposal

  end type gti_time_halving_step_policy

  !===================================================================!
  ! The one knob of the controller: how many attempts one advance
  ! may spend before failure is the lawful answer.
  !===================================================================!

  type :: gti_time_adaptive_controller_options

     integer :: max_attempts = 8

   contains

     procedure :: validate => options_validate

  end type gti_time_adaptive_controller_options

  !===================================================================!
  ! What one advance reports: every attempt's proposal and
  ! evidence, which one was kept, where it appended, and the
  ! growth step that carried it.
  !===================================================================!

  type :: gti_time_adaptive_attempt

     type(gti_time_step_proposal)    :: proposal
     type(gti_time_attempt_evidence) :: evidence
     logical :: kept = .false.

  end type gti_time_adaptive_attempt

  type :: gti_time_adaptive_advance_result

     logical :: accepted = .false.
     integer :: attempts = 0
     integer :: appended_vertex = 0
     integer :: appended_relation = 0

     type(gti_time_adaptive_attempt), allocatable :: attempt(:)

     type(gti_time_growth_step_result) :: growth

  end type gti_time_adaptive_advance_result

  !===================================================================!
  ! The stateless controller verb. The types keep their public
  ! singular names; Fortran denies a type its host module's name,
  ! so the module speaks in the plural.
  !===================================================================!

  type :: gti_time_adaptive_controller

   contains

     procedure :: advance

  end type gti_time_adaptive_controller

contains

  pure subroutine options_validate(this)

    class(gti_time_adaptive_controller_options), intent(in) :: this

    if (this % max_attempts < 1) then
       error stop 'gti_time_adaptive_controller: attempt budget is positive'
    end if

  end subroutine options_validate

  !===================================================================!
  ! One adaptive advance: propose, try reversible growth, measure,
  ! judge, keep or discard and retry - until acceptance or a spent
  ! budget. A spent budget is a lawful answer; the graph stands
  ! exactly as it was.
  !===================================================================!

  subroutine advance(this, residual_form, graph, design, policy, &
       & forward_options, options, result)

    class(gti_time_adaptive_controller)        , intent(in)    :: this
    class(gti_differentiable_form)             , intent(in)    :: residual_form
    type(gti_time_graph)                       , intent(inout) :: graph
    type(gti_design_bundle)                    , intent(in)    :: design
    class(gti_time_step_policy)                , intent(inout) :: policy
    type(gti_time_forward_options)             , intent(in)    :: forward_options
    type(gti_time_adaptive_controller_options) , intent(in)    :: options
    type(gti_time_adaptive_advance_result)     , intent(inout) :: result

    type(gti_time_adaptive_growth_driver) :: grower
    type(gti_time_growth_candidate)       :: candidate
    type(gti_time_growth_step_result)     :: step
    type(gti_time_step_proposal)          :: proposal
    type(gti_time_attempt_evidence)       :: evidence

    integer :: attempt_number, history_count
    logical :: accept

    call options % validate()

    if (graph % num_vertices() < 1) then
       error stop 'gti_time_adaptive_controller: the graph carries an initial vertex'
    end if

    result % accepted = .false.
    result % attempts = 0
    result % appended_vertex = 0
    result % appended_relation = 0
    if (allocated(result % attempt)) deallocate(result % attempt)
    allocate(result % attempt(0))
    result % growth = gti_time_growth_step_result()

    history_count = graph % num_vertices()

    call policy % first_proposal(history_count, proposal)

    do attempt_number = 1, options % max_attempts

       call require_lawful_proposal(proposal, history_count)

       call build_candidate(graph, proposal, candidate)

       !--------------------------------------------------------------!
       ! Try, committed: a converged candidate joins the graph, and
       ! the verdict comes AFTER - a rejection discards the
       ! committed step whole, through the same transactional seat.
       !--------------------------------------------------------------!

       call grower % try_candidate(residual_form, graph, candidate, design, &
            & forward_options, .true., step)

       evidence = gti_time_attempt_evidence()
       evidence % converged = step % converged
       evidence % step      = proposal % step
       evidence % order     = proposal % order
       evidence % attempt   = attempt_number

       if (step % converged) then
          call measure_estimate(graph, history_count, proposal, step, evidence)
          call policy % judge(evidence, accept)
       else
          ! an unconverged candidate is never kept, and the growth
          ! seat has already rolled it back
          accept = .false.
       end if

       call record_attempt(result, proposal, evidence, &
            & step % converged .and. accept)
       result % attempts = attempt_number

       if (step % converged .and. accept) then
          result % accepted          = .true.
          result % appended_vertex   = step % appended_vertex
          result % appended_relation = step % appended_relation
          result % growth            = step
          return
       end if

       if (step % converged) then
          call grower % discard_last_candidate(graph, step % appended_vertex, &
               & step % appended_relation)
       end if

       if (attempt_number < options % max_attempts) then
          call policy % retry_proposal(evidence, proposal)
       end if

    end do

  end subroutine advance

  !===================================================================!
  ! The proposal gate: a policy that proposes a nonpositive step,
  ! an unsupported order, or an order the history cannot carry has
  ! contradicted itself, and the contradiction dies here - not in
  ! a motif builder's message far from the deciding seat.
  !===================================================================!

  pure subroutine require_lawful_proposal(proposal, history_count)

    type(gti_time_step_proposal), intent(in) :: proposal
    integer                     , intent(in) :: history_count

    if (proposal % step <= 0.0_dp) then
       error stop 'gti_time_adaptive_controller: policy proposes a positive step'
    end if

    if (proposal % order /= 1 .and. proposal % order /= 2) then
       error stop 'gti_time_adaptive_controller: policy proposes a supported order'
    end if

    if (proposal % order == 2 .and. history_count < 2) then
       error stop 'gti_time_adaptive_controller: order two needs two history vertices'
    end if

  end subroutine require_lawful_proposal

  !===================================================================!
  ! Build one candidate from the graph's tail: a warm-started
  ! vertex at t_last + h, and a BDF relation whose variable-step
  ! rows price the graph's own spacing - order 1 from [h], order 2
  ! from [h, h_prev].
  !===================================================================!

  subroutine build_candidate(graph, proposal, candidate)

    type(gti_time_graph)           , intent(in)  :: graph
    type(gti_time_step_proposal)   , intent(in)  :: proposal
    type(gti_time_growth_candidate), intent(out) :: candidate

    type(gti_time_motif_builder) :: builder
    real(dp) :: t_last, h_previous
    integer  :: last

    last   = graph % num_vertices()
    t_last = graph % vertex(last) % sample % time

    candidate % vertex = graph % vertex(last)
    candidate % vertex % sample % time = t_last + proposal % step
    candidate % vertex % has_solution  = .false.

    select case (proposal % order)
    case (1)
       candidate % relation % sample_vertex  = [last, -1]
       candidate % relation % unknown_sample = 2
       call builder % bdf_variable(1, [proposal % step], &
            & candidate % relation % motif)
    case default
       h_previous = t_last - graph % vertex(last - 1) % sample % time
       candidate % relation % sample_vertex  = [last - 1, last, -1]
       candidate % relation % unknown_sample = 3
       call builder % bdf_variable(2, [proposal % step, h_previous], &
            & candidate % relation % motif)
    end select

    candidate % relation % evaluation_time = t_last + proposal % step

  end subroutine build_candidate

  !===================================================================!
  ! Measure one converged attempt: the estimate is the distance
  ! between the solved candidate and the extrapolating predictor
  ! through the graph's tail - constant through one history
  ! vertex, linear through two. Measurement only; the judgment is
  ! the policy's.
  !===================================================================!

  subroutine measure_estimate(graph, history_count, proposal, step, evidence)

    type(gti_time_graph)             , intent(in)    :: graph
    integer                          , intent(in)    :: history_count
    type(gti_time_step_proposal)     , intent(in)    :: proposal
    type(gti_time_growth_step_result), intent(in)    :: step
    type(gti_time_attempt_evidence)  , intent(inout) :: evidence

    real(dp), allocatable :: q_new(:), q_last(:), q_previous(:), q_predicted(:)
    real(dp) :: h_previous

    ! a defined descriptor quiets a gfortran maybe-uninitialized
    ! false positive; both branches replace it wholesale
    allocate(q_predicted(0))

    call step % forward_step % q % get_real(q_new)
    call read_vertex_q(graph, history_count, q_last)

    if (history_count >= 2) then
       call read_vertex_q(graph, history_count - 1, q_previous)
       h_previous = graph % vertex(history_count) % sample % time - &
            &       graph % vertex(history_count - 1) % sample % time
       q_predicted = q_last + (proposal % step / h_previous) * &
            & (q_last - q_previous)
    else
       q_predicted = q_last
    end if

    evidence % error_estimate = norm2(q_new - q_predicted)
    evidence % has_estimate   = .true.

  end subroutine measure_estimate

  subroutine read_vertex_q(graph, vertex_index, values)

    type(gti_time_graph) , intent(in)  :: graph
    integer              , intent(in)  :: vertex_index
    real(dp), allocatable, intent(out) :: values(:)

    call graph % vertex(vertex_index) % sample % state % &
         & component(1 + GTI_STATE_Q) % value % get_real_vector(values)

  end subroutine read_vertex_q

  !===================================================================!
  ! Land one attempt on the record.
  !===================================================================!

  subroutine record_attempt(result, proposal, evidence, kept)

    type(gti_time_adaptive_advance_result), intent(inout) :: result
    type(gti_time_step_proposal)          , intent(in)    :: proposal
    type(gti_time_attempt_evidence)       , intent(in)    :: evidence
    logical                               , intent(in)    :: kept

    type(gti_time_adaptive_attempt), allocatable :: grown(:)
    integer :: n

    n = size(result % attempt)
    allocate(grown(n + 1))
    grown(1:n) = result % attempt
    grown(n + 1) % proposal = proposal
    grown(n + 1) % evidence = evidence
    grown(n + 1) % kept     = kept
    call move_alloc(grown, result % attempt)

  end subroutine record_attempt

  !===================================================================!
  ! The reference policy's three verbs.
  !===================================================================!

  subroutine halving_first_proposal(this, history_count, proposal)

    class(gti_time_halving_step_policy), intent(inout) :: this
    integer                            , intent(in)    :: history_count
    type(gti_time_step_proposal)       , intent(out)   :: proposal

    if (this % initial_step <= 0.0_dp) then
       error stop 'gti_time_halving_step_policy: initial step is positive'
    end if

    if (this % error_tolerance <= 0.0_dp) then
       error stop 'gti_time_halving_step_policy: error tolerance is positive'
    end if

    if (this % preferred_order /= 1 .and. this % preferred_order /= 2) then
       error stop 'gti_time_halving_step_policy: preferred order is supported'
    end if

    proposal % step = this % initial_step

    if (history_count >= 2) then
       proposal % order = this % preferred_order
    else
       proposal % order = 1
    end if

  end subroutine halving_first_proposal

  subroutine halving_judge(this, evidence, accept)

    class(gti_time_halving_step_policy), intent(inout) :: this
    type(gti_time_attempt_evidence)    , intent(in)    :: evidence
    logical                            , intent(out)   :: accept

    accept = evidence % converged .and. evidence % has_estimate .and. &
         & evidence % error_estimate <= this % error_tolerance

  end subroutine halving_judge

  subroutine halving_retry_proposal(this, evidence, proposal)

    class(gti_time_halving_step_policy), intent(inout) :: this
    type(gti_time_attempt_evidence)    , intent(in)    :: evidence
    type(gti_time_step_proposal)       , intent(out)   :: proposal

    associate(unread => this)
    end associate

    proposal % step  = 0.5_dp * evidence % step
    proposal % order = evidence % order

  end subroutine halving_retry_proposal

end module gti_time_adaptive_controllers
