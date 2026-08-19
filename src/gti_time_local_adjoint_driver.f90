!=====================================================================!
! GTI TIME-LOCAL ADJOINT DRIVER
!
! The first local adjoint solve. At a converged local state,
!
!      R(q*, xi, t) = 0,      F(q*, xi, t) scalar,
!
! the total and constraint variations
!
!      dF = F_q [dq] + F_xi [eta]
!      J_q dq + R_xi [eta] = 0
!
! eliminate dq through one adjoint vector: define lambda by
!
!      J_q^T lambda = F_q^T,
!
! and the design-direction gradient action follows without ever
! solving for dq:
!
!      dF/dxi [eta]  =  F_xi [eta]  -  lambda^T R_xi [eta].
!
! J_q is the very Jacobian Newton and the tangent used - the
! motif-row action reused column by column from the local Newton
! driver, then TRANSPOSED for the adjoint solve. The functional's
! q-gradient is assembled the same way the Jacobian is: one
! order-1 state partial action per live motif row per basis
! direction, each answer a scalar because F is a scalar - never a
! finite difference, never a stored tensor.
!
! The dense transposed system is solved through the existing
! graph minimization tower: transpose(jacobian) rides a stencil -
! a matrix is a graph with weights on its edges - and dense_direct
! eliminates it. A singular pivot is refused loudly: one linear
! solve has no lawful way to limp.
!
! The driver carries nothing: no form, no motif, no samples, no
! graph, no scheme table, no map. Still local: no time graph, no
! chaining, no adaptivity, no higher-order chain-rule assembly.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_time_local_adjoint_drivers

  use iso_fortran_env      , only : dp => REAL64
  use gti_value_buffers    , only : gti_value_buffer
  use gti_design_bundles   , only : gti_design_bundle
  use gti_form_interface   , only : gti_differentiable_form, &
       & gti_partial_request, gti_direction_bundle, &
       & GTI_ARG_STATE, GTI_ARG_DESIGN
  use gti_form_evaluators  , only : gti_form_evaluator
  use gti_time_local_schemes, only : gti_time_sample, gti_time_motif, &
       & gti_time_local_residual_evaluator, gti_evaluation_point
  use gti_time_local_unknown_problems, only : gti_time_local_unknown_problem
  use gti_time_local_newton_drivers  , only : gti_time_local_newton_driver
  ! the dense Jacobian rides class_graph_stencil as a weighted-edge
  ! graph inside this adapter, and dense_direct minimizes it there
  use class_graph_dense_direct, only : solve_dense_matrix_with_dense_direct

  implicit none

  private
  public :: gti_time_local_adjoint_driver
  public :: gti_time_local_adjoint_result

  !===================================================================!
  ! What one local adjoint solve reports: that it solved, how
  ! exactly the transposed system was satisfied, the adjoint
  ! vector, and the functional q-gradient it answered.
  !===================================================================!

  type :: gti_time_local_adjoint_result

     logical  :: solved = .false.
     real(dp) :: linear_residual_norm = huge(1.0_dp)
     type(gti_value_buffer) :: lambda
     type(gti_value_buffer) :: functional_q_gradient

  end type gti_time_local_adjoint_result

  !===================================================================!
  ! The stateless verb-set. The types keep their public singular
  ! names; Fortran denies a type its host module's name, so the
  ! module speaks in the plural.
  !===================================================================!

  type :: gti_time_local_adjoint_driver

   contains

     procedure :: functional_q_gradient  => driver_functional_q_gradient
     procedure :: residual_design_action => driver_residual_design_action
     procedure :: functional_design_action => driver_functional_design_action
     procedure :: solve_adjoint          => driver_solve_adjoint
     procedure :: design_gradient_action => driver_design_gradient_action

  end type gti_time_local_adjoint_driver

contains

  !===================================================================!
  ! The dense vector g with  g^T dq = dF/dq_unknown [dq]: for each
  ! basis direction of the unknown, sum one order-1 state partial
  ! action of the functional per live motif row - the functional's
  ! side of the very chain rule the Jacobian action uses. Every
  ! answer must be a scalar, for F is a scalar; the assembled
  ! gradient must wear the unknown's shape.
  !===================================================================!

  subroutine driver_functional_q_gradient(this, functional, motif, samples, &
       & unknown_index, q_star, design, time, gradient)

    class(gti_time_local_adjoint_driver), intent(in)    :: this
    class(gti_differentiable_form)      , intent(in)    :: functional
    type(gti_time_motif)                , intent(in)    :: motif
    type(gti_time_sample)               , intent(in)    :: samples(:)
    integer                             , intent(in)    :: unknown_index
    type(gti_value_buffer)              , intent(in)    :: q_star
    type(gti_design_bundle)             , intent(in)    :: design
    real(dp)                            , intent(in)    :: time
    type(gti_value_buffer)              , intent(inout) :: gradient

    type(gti_time_local_unknown_problem)    :: problem
    type(gti_time_local_residual_evaluator) :: point_builder
    type(gti_form_evaluator)                :: evaluator
    type(gti_time_sample), allocatable      :: work_samples(:)
    type(gti_evaluation_point)              :: point
    type(gti_partial_request)               :: request
    type(gti_direction_bundle)              :: direction
    type(gti_value_buffer)                  :: term

    real(dp), allocatable :: basis(:), term_values(:), gradient_values(:)
    real(dp) :: entry_value, weight
    integer :: n, j, m

    n = 0
    if (allocated(q_star % rvals)) n = size(q_star % rvals)
    if (n == 0) then
       error stop 'gti_time_local_adjoint_driver: q star has values'
    end if

    ! a defined descriptor quiets a gfortran maybe-uninitialized
    ! false positive; the injection replaces it wholesale
    allocate(work_samples(0))

    call problem % inject_trial_q(samples, unknown_index, q_star, work_samples)
    call point_builder % build_point(motif, work_samples, design, time, point)

    allocate(basis(n), gradient_values(n))

    do j = 1, n

       basis    = 0.0_dp
       basis(j) = 1.0_dp

       entry_value = 0.0_dp

       do m = 1, motif % size()

          weight = motif % rule(m) % weight(unknown_index)
          if (weight == 0.0_dp) cycle

          direction % argument_kind   = GTI_ARG_STATE
          direction % state_component = motif % rule(m) % state_component
          call direction % values % set_real(weight * basis, ncomp=q_star % ncomp)

          request = gti_partial_request(order=1, &
               & argument_kind  =[GTI_ARG_STATE], &
               & state_component=[motif % rule(m) % state_component])

          call evaluator % partial_action(functional, point, request, &
               & [direction], term)

          if (term % nentries /= 1 .or. term % ncomp /= 1) then
             error stop 'gti_time_local_adjoint_driver: functional action is scalar'
          end if

          call term % get_real(term_values)
          entry_value = entry_value + term_values(1)

       end do

       gradient_values(j) = entry_value

    end do

    call gradient % set_real(gradient_values)

    if (gradient % nentries /= q_star % nentries .or. &
         & gradient % ncomp /= q_star % ncomp) then
       error stop 'gti_time_local_adjoint_driver: functional q gradient size matches unknown size'
    end if

  end subroutine driver_functional_q_gradient

  !===================================================================!
  ! R_xi [eta] at the point built from q*: one order-1 design
  ! partial action of the residual form, unnegated - the raw
  ! constraint sensitivity the gradient contraction consumes.
  !===================================================================!

  subroutine driver_residual_design_action(this, residual_form, motif, samples, &
       & unknown_index, q_star, design, design_direction, time, action)

    class(gti_time_local_adjoint_driver), intent(in)    :: this
    class(gti_differentiable_form)      , intent(in)    :: residual_form
    type(gti_time_motif)                , intent(in)    :: motif
    type(gti_time_sample)               , intent(in)    :: samples(:)
    integer                             , intent(in)    :: unknown_index
    type(gti_value_buffer)              , intent(in)    :: q_star
    type(gti_design_bundle)             , intent(in)    :: design
    type(gti_value_buffer)              , intent(in)    :: design_direction
    real(dp)                            , intent(in)    :: time
    type(gti_value_buffer)              , intent(inout) :: action

    call design_partial_action(residual_form, motif, samples, unknown_index, &
         & q_star, design, design_direction, time, action)

    if (action % ncomp /= 1 .or. .not. allocated(action % rvals)) then
       error stop 'gti_time_local_adjoint_driver: residual design action is a vector'
    end if

  end subroutine driver_residual_design_action

  !===================================================================!
  ! F_xi [eta] at the point built from q*: the same design partial
  ! action, asked of the functional - and a scalar, for F is one.
  !===================================================================!

  subroutine driver_functional_design_action(this, functional, motif, samples, &
       & unknown_index, q_star, design, design_direction, time, action)

    class(gti_time_local_adjoint_driver), intent(in)    :: this
    class(gti_differentiable_form)      , intent(in)    :: functional
    type(gti_time_motif)                , intent(in)    :: motif
    type(gti_time_sample)               , intent(in)    :: samples(:)
    integer                             , intent(in)    :: unknown_index
    type(gti_value_buffer)              , intent(in)    :: q_star
    type(gti_design_bundle)             , intent(in)    :: design
    type(gti_value_buffer)              , intent(in)    :: design_direction
    real(dp)                            , intent(in)    :: time
    type(gti_value_buffer)              , intent(inout) :: action

    call design_partial_action(functional, motif, samples, unknown_index, &
         & q_star, design, design_direction, time, action)

    if (action % nentries /= 1 .or. action % ncomp /= 1) then
       error stop 'gti_time_local_adjoint_driver: functional action is scalar'
    end if

  end subroutine driver_functional_design_action

  !===================================================================!
  ! Solve  J_q^T lambda = F_q^T  for one unknown sample: the
  ! functional's q-gradient as the seed, the dense J_q from
  ! Newton's exact action on basis directions - the same Jacobian
  ! the primal and tangent solves used - transposed, eliminated,
  ! and the linear residual held up for inspection. One linear
  ! solve, never an iteration.
  !===================================================================!

  subroutine driver_solve_adjoint(this, residual_form, functional, motif, &
       & samples, unknown_index, q_star, design, time, singular_tolerance, result)

    class(gti_time_local_adjoint_driver), intent(in)    :: this
    class(gti_differentiable_form)      , intent(in)    :: residual_form
    class(gti_differentiable_form)      , intent(in)    :: functional
    type(gti_time_motif)                , intent(in)    :: motif
    type(gti_time_sample)               , intent(in)    :: samples(:)
    integer                             , intent(in)    :: unknown_index
    type(gti_value_buffer)              , intent(in)    :: q_star
    type(gti_design_bundle)             , intent(in)    :: design
    real(dp)                            , intent(in)    :: time
    real(dp)                            , intent(in)    :: singular_tolerance
    type(gti_time_local_adjoint_result) , intent(inout) :: result

    type(gti_time_local_newton_driver) :: newton
    type(gti_value_buffer)             :: seed, dq_basis, column
    real(dp) :: achieved

    real(dp), allocatable :: seed_values(:), column_values(:)
    real(dp), allocatable :: jacobian(:,:), transposed(:,:), lambda_values(:), basis(:)
    integer :: n, j

    if (singular_tolerance <= 0.0_dp) then
       error stop 'gti_time_local_adjoint_driver: singular tolerance is positive'
    end if

    n = 0
    if (allocated(q_star % rvals)) n = size(q_star % rvals)
    if (n == 0) then
       error stop 'gti_time_local_adjoint_driver: q star has values'
    end if

    result % solved = .false.
    result % linear_residual_norm = huge(1.0_dp)

    call this % functional_q_gradient(functional, motif, samples, unknown_index, &
         & q_star, design, time, seed)
    call seed % get_real(seed_values)

    !----------------------------------------------------------------!
    ! The dense J_q, one exact Newton action per basis column,
    ! each held to the unknown's size.
    !----------------------------------------------------------------!

    allocate(jacobian(n, n), basis(n))

    do j = 1, n
       basis    = 0.0_dp
       basis(j) = 1.0_dp
       call dq_basis % set_real(basis, ncomp=q_star % ncomp)
       call newton % jacobian_action(residual_form, motif, samples, unknown_index, &
            & q_star, dq_basis, design, time, column)
       call column % get_real(column_values)
       if (size(column_values) /= n) then
          error stop 'gti_time_local_adjoint_driver: residual size matches unknown size'
       end if
       jacobian(:, j) = column_values
    end do

    transposed = transpose(jacobian)

    call solve_dense_matrix_with_dense_direct(transposed, seed_values, &
         & singular_tolerance, lambda_values, achieved)

    result % linear_residual_norm = &
         & norm2(matmul(transposed, lambda_values) - seed_values)
    result % solved = .true.

    call result % lambda % set_real(lambda_values, ncomp=q_star % ncomp)
    result % functional_q_gradient = seed

  end subroutine driver_solve_adjoint

  !===================================================================!
  ! The scalar design-gradient action,
  !
  !      dF/dxi [eta] = F_xi [eta] - lambda^T R_xi [eta],
  !
  ! answered as a [1, 1] buffer. The two explicit actions are
  ! taken first, so a misshapen residual sensitivity is refused
  ! with its own law before any linear algebra runs; the adjoint
  ! solve then supplies lambda, and one contraction finishes it.
  !===================================================================!

  subroutine driver_design_gradient_action(this, residual_form, functional, &
       & motif, samples, unknown_index, q_star, design, design_direction, &
       & time, singular_tolerance, gradient_action)

    class(gti_time_local_adjoint_driver), intent(in)    :: this
    class(gti_differentiable_form)      , intent(in)    :: residual_form
    class(gti_differentiable_form)      , intent(in)    :: functional
    type(gti_time_motif)                , intent(in)    :: motif
    type(gti_time_sample)               , intent(in)    :: samples(:)
    integer                             , intent(in)    :: unknown_index
    type(gti_value_buffer)              , intent(in)    :: q_star
    type(gti_design_bundle)             , intent(in)    :: design
    type(gti_value_buffer)              , intent(in)    :: design_direction
    real(dp)                            , intent(in)    :: time
    real(dp)                            , intent(in)    :: singular_tolerance
    type(gti_value_buffer)              , intent(inout) :: gradient_action

    type(gti_time_local_adjoint_result) :: result
    type(gti_value_buffer)              :: r_action, f_action

    real(dp), allocatable :: r_values(:), f_values(:), lambda_values(:)
    integer :: n

    n = 0
    if (allocated(q_star % rvals)) n = size(q_star % rvals)
    if (n == 0) then
       error stop 'gti_time_local_adjoint_driver: q star has values'
    end if

    call this % residual_design_action(residual_form, motif, samples, &
         & unknown_index, q_star, design, design_direction, time, r_action)
    call r_action % get_real(r_values)
    if (size(r_values) /= n) then
       error stop 'gti_time_local_adjoint_driver: residual design action size matches unknown size'
    end if

    call this % functional_design_action(functional, motif, samples, &
         & unknown_index, q_star, design, design_direction, time, f_action)
    call f_action % get_real(f_values)

    call this % solve_adjoint(residual_form, functional, motif, samples, &
         & unknown_index, q_star, design, time, singular_tolerance, result)
    call result % lambda % get_real(lambda_values)

    call gradient_action % set_real([f_values(1) - dot_product(lambda_values, r_values)])

  end subroutine driver_design_gradient_action

  !===================================================================!
  ! The shared design partial: inject q*, build the point, ask one
  ! order-1 GTI_ARG_DESIGN action of the given form.
  !===================================================================!

  subroutine design_partial_action(form, motif, samples, unknown_index, &
       & q_star, design, design_direction, time, action)

    class(gti_differentiable_form), intent(in)    :: form
    type(gti_time_motif)          , intent(in)    :: motif
    type(gti_time_sample)         , intent(in)    :: samples(:)
    integer                       , intent(in)    :: unknown_index
    type(gti_value_buffer)        , intent(in)    :: q_star
    type(gti_design_bundle)       , intent(in)    :: design
    type(gti_value_buffer)        , intent(in)    :: design_direction
    real(dp)                      , intent(in)    :: time
    type(gti_value_buffer)        , intent(inout) :: action

    type(gti_time_local_unknown_problem)    :: problem
    type(gti_time_local_residual_evaluator) :: point_builder
    type(gti_form_evaluator)                :: evaluator
    type(gti_time_sample), allocatable      :: work_samples(:)
    type(gti_evaluation_point)              :: point
    type(gti_partial_request)               :: request
    type(gti_direction_bundle)              :: direction

    integer :: nstar, ndirection

    nstar = 0
    if (allocated(q_star % rvals)) nstar = size(q_star % rvals)
    if (nstar == 0) then
       error stop 'gti_time_local_adjoint_driver: q star has values'
    end if

    ndirection = 0
    if (allocated(design_direction % rvals)) ndirection = size(design_direction % rvals)
    if (ndirection == 0) then
       error stop 'gti_time_local_adjoint_driver: design direction has values'
    end if

    ! a defined descriptor quiets a gfortran maybe-uninitialized
    ! false positive; the injection replaces it wholesale
    allocate(work_samples(0))

    call problem % inject_trial_q(samples, unknown_index, q_star, work_samples)
    call point_builder % build_point(motif, work_samples, design, time, point)

    direction % argument_kind = GTI_ARG_DESIGN
    direction % values        = design_direction

    request = gti_partial_request(order=1, argument_kind=[GTI_ARG_DESIGN])

    call evaluator % partial_action(form, point, request, [direction], action)

  end subroutine design_partial_action

end module gti_time_local_adjoint_drivers
