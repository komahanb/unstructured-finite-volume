!=====================================================================!
! GTI TIME-LOCAL NEWTON DRIVER
!
! The first local implicit solve: for one unknown sample of a
! local time motif, drive
!
!      residual(q_trial) = 0
!
! by Newton iteration, entirely through the seats already built -
! the unknown problem answers the residual, the motif rows and the
! form's partial actions answer the exact Jacobian action, and a
! private dense elimination answers the step. Nothing here is a
! scheme: the same three verbs serve BDF, DIRK, ABM, or any motif
! a builder mints.
!
! The Jacobian law. All samples but the unknown are fixed, so the
! residual depends on the trial alone, and for each motif row
!
!      component_m = sum_i w_m(i) q_i
!
! the unknown contributes  w_m(unknown) dq  to component m. Hence
! the exact first-order action is
!
!      dR/dq_unknown [dq]  =  sum_m  D_m R [ w_m(unknown) dq ],
!
! one order-1 partial action per active row - the local chain rule
! for the unknown sample only, never a finite difference, never a
! scheme-specific formula.
!
! The dense system J step = -r is built column by column from that
! action on basis directions and solved by Gaussian elimination
! with partial pivoting, a private helper for this small local
! driver only - no linear-solver abstraction is introduced, and no
! existing solver is called.
!
! Failure to converge is an answer, not an error: the result
! carries converged/.false., the iterations spent, and the norms
! where the iteration stopped. Unlawful options, valueless trials,
! misshapen directions, non-vector residuals, and singular pivots
! are refused loudly.
!
! The driver carries nothing: no form, no motif, no samples, no
! graph, no scheme table, no map. Still local: no time graph, no
! chaining, no adaptivity, no tangent, no adjoint.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_time_local_newton_drivers

  use iso_fortran_env      , only : dp => REAL64
  use gti_value_buffers    , only : gti_value_buffer
  use gti_design_bundles   , only : gti_design_bundle
  use gti_form_interface   , only : gti_differentiable_form, &
       & gti_partial_request, gti_direction_bundle, GTI_ARG_STATE
  use gti_form_evaluators  , only : gti_form_evaluator
  use gti_time_local_schemes, only : gti_time_sample, gti_time_motif, &
       & gti_time_local_residual_evaluator, gti_evaluation_point
  use gti_time_local_unknown_problems, only : gti_time_local_unknown_problem

  implicit none

  private
  public :: gti_time_local_newton_driver
  public :: gti_time_local_newton_options
  public :: gti_time_local_newton_result

  !===================================================================!
  ! The knobs of one local solve, each with a lawful range the
  ! validation refuses to leave.
  !===================================================================!

  type :: gti_time_local_newton_options

     integer  :: max_iterations = 20
     real(dp) :: residual_tolerance = 1.0e-10_dp
     real(dp) :: step_tolerance     = 1.0e-12_dp
     real(dp) :: singular_tolerance = 1.0e-14_dp

   contains

     procedure :: validate => options_validate

  end type gti_time_local_newton_options

  !===================================================================!
  ! What one local solve reports: whether it converged, how much
  ! it spent, where it stopped, and the trial and residual there.
  !===================================================================!

  type :: gti_time_local_newton_result

     logical  :: converged = .false.
     integer  :: iterations = 0
     real(dp) :: residual_norm = huge(1.0_dp)
     real(dp) :: step_norm     = huge(1.0_dp)
     type(gti_value_buffer) :: q
     type(gti_value_buffer) :: residual

  end type gti_time_local_newton_result

  !===================================================================!
  ! The stateless verb-set. The types keep their public singular
  ! names; Fortran denies a type its host module's name, so the
  ! module speaks in the plural.
  !===================================================================!

  type :: gti_time_local_newton_driver

   contains

     procedure :: residual        => driver_residual
     procedure :: jacobian_action => driver_jacobian_action
     procedure :: solve           => driver_solve

  end type gti_time_local_newton_driver

contains

  !===================================================================!
  ! The options laws, refused in order.
  !===================================================================!

  pure subroutine options_validate(this)

    class(gti_time_local_newton_options), intent(in) :: this

    if (this % max_iterations < 1) then
       error stop 'gti_time_local_newton_driver: max_iterations is positive'
    end if

    if (this % residual_tolerance < 0.0_dp) then
       error stop 'gti_time_local_newton_driver: residual tolerance is nonnegative'
    end if

    if (this % step_tolerance < 0.0_dp) then
       error stop 'gti_time_local_newton_driver: step tolerance is nonnegative'
    end if

    if (this % singular_tolerance <= 0.0_dp) then
       error stop 'gti_time_local_newton_driver: singular tolerance is positive'
    end if

  end subroutine options_validate

  !===================================================================!
  ! residual(q_trial), held to Newton's needs: the answer must be
  ! a flat real vector - one component per entry - or it cannot
  ! seed a square system.
  !===================================================================!

  subroutine driver_residual(this, form, motif, samples, unknown_index, &
       & q_trial, design, time, output)

    class(gti_time_local_newton_driver), intent(in)    :: this
    class(gti_differentiable_form)     , intent(in)    :: form
    type(gti_time_motif)               , intent(in)    :: motif
    type(gti_time_sample)              , intent(in)    :: samples(:)
    integer                            , intent(in)    :: unknown_index
    type(gti_value_buffer)             , intent(in)    :: q_trial
    type(gti_design_bundle)            , intent(in)    :: design
    real(dp)                           , intent(in)    :: time
    type(gti_value_buffer)             , intent(inout) :: output

    type(gti_time_local_unknown_problem) :: problem

    call problem % value(form, motif, samples, unknown_index, q_trial, &
         & design, time, output)

    call require_vector(output)

  end subroutine driver_residual

  !===================================================================!
  ! The exact first-order action dR/dq_unknown [dq]: one order-1
  ! partial action per motif row with a live weight at the unknown,
  ! each along  w_m(unknown) dq  in that row's state component,
  ! summed. The local chain rule for the unknown sample only -
  ! never a finite difference.
  !===================================================================!

  subroutine driver_jacobian_action(this, form, motif, samples, unknown_index, &
       & q_trial, dq_trial, design, time, output)

    class(gti_time_local_newton_driver), intent(in)    :: this
    class(gti_differentiable_form)     , intent(in)    :: form
    type(gti_time_motif)               , intent(in)    :: motif
    type(gti_time_sample)              , intent(in)    :: samples(:)
    integer                            , intent(in)    :: unknown_index
    type(gti_value_buffer)             , intent(in)    :: q_trial
    type(gti_value_buffer)             , intent(in)    :: dq_trial
    type(gti_design_bundle)            , intent(in)    :: design
    real(dp)                           , intent(in)    :: time
    type(gti_value_buffer)             , intent(inout) :: output

    type(gti_time_local_unknown_problem)    :: problem
    type(gti_time_local_residual_evaluator) :: point_builder
    type(gti_form_evaluator)                :: evaluator
    type(gti_time_sample), allocatable      :: work_samples(:)
    type(gti_evaluation_point)              :: point
    type(gti_partial_request)               :: request
    type(gti_direction_bundle)              :: direction
    type(gti_value_buffer)                  :: term

    real(dp), allocatable :: direction_values(:)
    real(dp), allocatable :: sum_values(:), term_values(:)
    integer , allocatable :: signature(:)
    integer :: m, ntrial, ndirection, term_ncomp
    logical :: started
    real(dp) :: weight

    ntrial = 0
    if (allocated(q_trial % rvals)) ntrial = size(q_trial % rvals)
    if (ntrial == 0) then
       error stop 'gti_time_local_newton_driver: trial q has values'
    end if

    ndirection = 0
    if (allocated(dq_trial % rvals)) ndirection = size(dq_trial % rvals)
    if (dq_trial % nentries /= q_trial % nentries .or. &
         & dq_trial % ncomp /= q_trial % ncomp .or. ndirection /= ntrial) then
       error stop 'gti_time_local_newton_driver: direction q shape matches trial q'
    end if

    ! a defined descriptor quiets a gfortran maybe-uninitialized
    ! false positive; the injection replaces it wholesale
    allocate(work_samples(0))

    call problem % inject_trial_q(samples, unknown_index, q_trial, work_samples)
    call point_builder % build_point(motif, work_samples, design, time, point)

    call dq_trial % get_real(direction_values)

    started    = .false.
    term_ncomp = 1

    do m = 1, motif % size()

       weight = motif % rule(m) % weight(unknown_index)
       if (weight == 0.0_dp) cycle

       direction % argument_kind   = GTI_ARG_STATE
       direction % state_component = motif % rule(m) % state_component
       call direction % values % set_real(weight * direction_values, &
            & ncomp=dq_trial % ncomp)

       request = gti_partial_request(order=1, &
            & argument_kind  =[GTI_ARG_STATE], &
            & state_component=[motif % rule(m) % state_component])

       call evaluator % partial_action(form, point, request, &
            & [direction], term)

       call term % get_real(term_values)
       term_ncomp = term % ncomp

       if (started) then
          sum_values = sum_values + term_values
       else
          sum_values = term_values
          started    = .true.
       end if

    end do

    if (.not. started) then
       signature = form % output_signature()
       if (size(signature) /= 2) then
          error stop 'gti_time_local_newton_driver: residual is a vector'
       end if
       sum_values = spread(0.0_dp, dim=1, ncopies=signature(1) * signature(2))
       term_ncomp = signature(2)
    end if

    call output % set_real(sum_values, ncomp=term_ncomp)

    call require_vector(output)

  end subroutine driver_jacobian_action

  !===================================================================!
  ! Newton: residual, dense Jacobian by columns of the exact
  ! action, eliminate, step, repeat. Convergence is a small
  ! residual; a stalled step without a small residual ends the
  ! iteration unconverged - an answer, never an error stop.
  !===================================================================!

  subroutine driver_solve(this, form, motif, samples, unknown_index, &
       & q_initial, design, time, options, result)

    class(gti_time_local_newton_driver), intent(in)    :: this
    class(gti_differentiable_form)     , intent(in)    :: form
    type(gti_time_motif)               , intent(in)    :: motif
    type(gti_time_sample)              , intent(in)    :: samples(:)
    integer                            , intent(in)    :: unknown_index
    type(gti_value_buffer)             , intent(in)    :: q_initial
    type(gti_design_bundle)            , intent(in)    :: design
    real(dp)                           , intent(in)    :: time
    type(gti_time_local_newton_options), intent(in)    :: options
    type(gti_time_local_newton_result) , intent(inout) :: result

    type(gti_value_buffer) :: q_current, r_current, dq_basis, column
    real(dp), allocatable  :: q_values(:), r_values(:), column_values(:)
    real(dp), allocatable  :: jacobian(:,:), step(:), basis(:)
    integer :: n, j

    call options % validate()

    call q_initial % get_real(q_values)
    if (size(q_values) == 0) then
       error stop 'gti_time_local_newton_driver: trial q has values'
    end if
    n = size(q_values)

    result % converged     = .false.
    result % iterations    = 0
    result % residual_norm = huge(1.0_dp)
    result % step_norm     = huge(1.0_dp)

    q_current = q_initial

    call this % residual(form, motif, samples, unknown_index, q_current, &
         & design, time, r_current)
    call r_current % get_real(r_values)
    if (size(r_values) /= n) then
       error stop 'gti_time_local_newton_driver: residual size matches unknown size'
    end if
    result % residual_norm = norm2(r_values)

    if (result % residual_norm <= options % residual_tolerance) then
       result % converged = .true.
    end if

    allocate(jacobian(n, n), basis(n))

    do while (.not. result % converged .and. &
         &    result % iterations < options % max_iterations)

       !-------------------------------------------------------------!
       ! The dense Jacobian, one exact action per basis column.
       !-------------------------------------------------------------!

       do j = 1, n
          basis    = 0.0_dp
          basis(j) = 1.0_dp
          call dq_basis % set_real(basis, ncomp=q_initial % ncomp)
          call this % jacobian_action(form, motif, samples, unknown_index, &
               & q_current, dq_basis, design, time, column)
          call column % get_real(column_values)
          jacobian(:, j) = column_values
       end do

       call dense_solve(jacobian, -r_values, options % singular_tolerance, step)

       q_values = q_values + step
       call q_current % set_real(q_values, ncomp=q_initial % ncomp)

       result % iterations = result % iterations + 1
       result % step_norm  = norm2(step)

       call this % residual(form, motif, samples, unknown_index, q_current, &
            & design, time, r_current)
       call r_current % get_real(r_values)
       result % residual_norm = norm2(r_values)

       if (result % residual_norm <= options % residual_tolerance) then
          result % converged = .true.
       else if (result % step_norm <= options % step_tolerance) then
          exit
       end if

    end do

    result % q        = q_current
    result % residual = r_current

  end subroutine driver_solve

  !===================================================================!
  ! Newton's shape law: the residual is a flat real vector, one
  ! component per entry.
  !===================================================================!

  pure subroutine require_vector(output)

    type(gti_value_buffer), intent(in) :: output

    if (output % ncomp /= 1 .or. .not. allocated(output % rvals)) then
       error stop 'gti_time_local_newton_driver: residual is a vector'
    end if

  end subroutine require_vector

  !===================================================================!
  ! Solve A x = b by Gaussian elimination with partial pivoting,
  ! on private copies. A pivot below the singular tolerance is
  ! refused - a flat Jacobian cannot point at a root. This helper
  ! serves the local driver alone; it is not a solver abstraction.
  !===================================================================!

  pure subroutine dense_solve(a, b, singular_tolerance, x)

    real(dp)             , intent(in)  :: a(:,:), b(:)
    real(dp)             , intent(in)  :: singular_tolerance
    real(dp), allocatable, intent(out) :: x(:)

    real(dp), allocatable :: m(:,:), r(:), row(:)
    real(dp) :: pivot_value, factor
    integer :: n, k, p, i

    n = size(b)
    m = a
    r = b
    allocate(row(n))

    do k = 1, n

       p = k - 1 + maxloc(abs(m(k:n, k)), dim=1)

       if (abs(m(p, k)) <= singular_tolerance) then
          error stop 'gti_time_local_newton_driver: dense Jacobian pivot is nonsingular'
       end if

       if (p /= k) then
          row       = m(k, :)
          m(k, :)   = m(p, :)
          m(p, :)   = row
          pivot_value = r(k)
          r(k)        = r(p)
          r(p)        = pivot_value
       end if

       do i = k + 1, n
          factor    = m(i, k) / m(k, k)
          m(i, k:n) = m(i, k:n) - factor * m(k, k:n)
          r(i)      = r(i) - factor * r(k)
       end do

    end do

    allocate(x(n))
    do k = n, 1, -1
       x(k) = (r(k) - dot_product(m(k, k+1:n), x(k+1:n))) / m(k, k)
    end do

  end subroutine dense_solve

end module gti_time_local_newton_drivers
