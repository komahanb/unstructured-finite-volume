!=====================================================================!
! GTI TIME-LOCAL TANGENT DRIVER
!
! The first local tangent solve. A converged local implicit state
! satisfies
!
!      R(q*, xi, t) = 0,
!
! and differentiating along a design direction eta gives
!
!      R_q [q_eta]  +  R_xi [eta]  =  0,
!
! so the tangent of the unknown sample solves the LINEAR system
!
!      J_q q_eta  =  - R_xi [eta],
!
! with J_q the very Jacobian Newton used - the motif-row action
!
!      J_q [dq] = sum_m D_m R [ w_m(unknown) dq ],
!
! reused verbatim from the local Newton driver, never rederived,
! never finite-differenced. The right-hand side is the explicit
! design part alone: one order-1 design partial action at the
! point built from q*, negated - no q-tangent effect belongs in
! it, and no adjoint, no higher-order chain-rule assembly, and no
! time graph appear anywhere here.
!
! The dense system is solved by a private pivoted elimination for
! this small local driver only - a deliberate duplicate of the
! Newton driver's helper; a later cleanup may factor a shared
! local dense solver, and neither is a solver abstraction. A
! singular pivot is refused loudly: this is one linear solve, not
! an iteration, and it has no lawful way to limp.
!
! The result carries the tangent, the right-hand side it answered,
! and the norm of the linear residual J_q q_eta - rhs, so a caller
! can hold the solve to account without repeating it.
!
! The driver carries nothing: no form, no motif, no samples, no
! graph, no scheme table, no map.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_time_local_tangent_drivers

  use iso_fortran_env      , only : dp => REAL64
  use gti_value_buffers    , only : gti_value_buffer
  use gti_design_bundles   , only : gti_design_bundle
  use gti_form_interface   , only : gti_differentiable_form, &
       & gti_partial_request, gti_direction_bundle, GTI_ARG_DESIGN
  use gti_form_evaluators  , only : gti_form_evaluator
  use gti_time_local_schemes, only : gti_time_sample, gti_time_motif, &
       & gti_time_local_residual_evaluator, gti_evaluation_point
  use gti_time_local_unknown_problems, only : gti_time_local_unknown_problem
  use gti_time_local_newton_drivers  , only : gti_time_local_newton_driver

  implicit none

  private
  public :: gti_time_local_tangent_driver
  public :: gti_time_local_tangent_result

  !===================================================================!
  ! What one local tangent solve reports: that it solved, how
  ! exactly the linear system was satisfied, and the tangent and
  ! right-hand side themselves.
  !===================================================================!

  type :: gti_time_local_tangent_result

     logical  :: solved = .false.
     real(dp) :: linear_residual_norm = huge(1.0_dp)
     type(gti_value_buffer) :: q_tangent
     type(gti_value_buffer) :: rhs

  end type gti_time_local_tangent_result

  !===================================================================!
  ! The stateless verb-pair. The types keep their public singular
  ! names; Fortran denies a type its host module's name, so the
  ! module speaks in the plural.
  !===================================================================!

  type :: gti_time_local_tangent_driver

   contains

     procedure :: design_rhs
     procedure :: solve_design_tangent

  end type gti_time_local_tangent_driver

contains

  !===================================================================!
  ! The explicit design right-hand side,
  !
  !      rhs = - R_xi [eta],
  !
  ! one order-1 design partial action at the point built from q*,
  ! negated. Only the explicit design part: the q-tangent effect
  ! is the left-hand side's business.
  !===================================================================!

  subroutine design_rhs(this, form, motif, samples, unknown_index, q_star, &
       & design, design_direction, time, rhs)

    class(gti_time_local_tangent_driver), intent(in)    :: this
    class(gti_differentiable_form)      , intent(in)    :: form
    type(gti_time_motif)                , intent(in)    :: motif
    type(gti_time_sample)               , intent(in)    :: samples(:)
    integer                             , intent(in)    :: unknown_index
    type(gti_value_buffer)              , intent(in)    :: q_star
    type(gti_design_bundle)             , intent(in)    :: design
    type(gti_value_buffer)              , intent(in)    :: design_direction
    real(dp)                            , intent(in)    :: time
    type(gti_value_buffer)              , intent(inout) :: rhs

    type(gti_time_local_unknown_problem)    :: problem
    type(gti_time_local_residual_evaluator) :: point_builder
    type(gti_form_evaluator)                :: evaluator
    type(gti_time_sample), allocatable      :: work_samples(:)
    type(gti_evaluation_point)              :: point
    type(gti_partial_request)               :: request
    type(gti_direction_bundle)              :: direction
    type(gti_value_buffer)                  :: term

    real(dp), allocatable :: term_values(:)
    integer :: nstar, ndirection

    nstar = 0
    if (allocated(q_star % rvals)) nstar = size(q_star % rvals)
    if (nstar == 0) then
       error stop 'gti_time_local_tangent_driver: q star has values'
    end if

    ndirection = 0
    if (allocated(design_direction % rvals)) ndirection = size(design_direction % rvals)
    if (ndirection == 0) then
       error stop 'gti_time_local_tangent_driver: design direction has values'
    end if

    ! a defined descriptor quiets a gfortran maybe-uninitialized
    ! false positive; the injection replaces it wholesale
    allocate(work_samples(0))

    call problem % inject_trial_q(samples, unknown_index, q_star, work_samples)
    call point_builder % build_point(motif, work_samples, design, time, point)

    direction % argument_kind = GTI_ARG_DESIGN
    direction % values        = design_direction

    request = gti_partial_request(order=1, argument_kind=[GTI_ARG_DESIGN])

    call evaluator % partial_action(form, point, request, [direction], term)

    call term % get_real(term_values)
    call rhs % set_real(-term_values, ncomp=term % ncomp)

    if (rhs % ncomp /= 1 .or. .not. allocated(rhs % rvals)) then
       error stop 'gti_time_local_tangent_driver: tangent rhs is a vector'
    end if

  end subroutine design_rhs

  !===================================================================!
  ! Solve  J_q q_eta = - R_xi [eta]  for one unknown sample: the
  ! explicit design right-hand side, the dense J_q from Newton's
  ! exact action on basis directions, one pivoted elimination, and
  ! the linear residual held up for inspection. One linear solve,
  ! never an iteration.
  !===================================================================!

  subroutine solve_design_tangent(this, form, motif, samples, unknown_index, &
       & q_star, design, design_direction, time, singular_tolerance, result)

    class(gti_time_local_tangent_driver), intent(in)    :: this
    class(gti_differentiable_form)      , intent(in)    :: form
    type(gti_time_motif)                , intent(in)    :: motif
    type(gti_time_sample)               , intent(in)    :: samples(:)
    integer                             , intent(in)    :: unknown_index
    type(gti_value_buffer)              , intent(in)    :: q_star
    type(gti_design_bundle)             , intent(in)    :: design
    type(gti_value_buffer)              , intent(in)    :: design_direction
    real(dp)                            , intent(in)    :: time
    real(dp)                            , intent(in)    :: singular_tolerance
    type(gti_time_local_tangent_result) , intent(inout) :: result

    type(gti_time_local_newton_driver) :: newton
    type(gti_value_buffer)             :: rhs_buffer, dq_basis, column

    real(dp), allocatable :: rhs_values(:), column_values(:)
    real(dp), allocatable :: jacobian(:,:), tangent(:), basis(:)
    integer :: n, j

    if (singular_tolerance <= 0.0_dp) then
       error stop 'gti_time_local_tangent_driver: singular tolerance is positive'
    end if

    n = 0
    if (allocated(q_star % rvals)) n = size(q_star % rvals)
    if (n == 0) then
       error stop 'gti_time_local_tangent_driver: q star has values'
    end if

    result % solved = .false.
    result % linear_residual_norm = huge(1.0_dp)

    call this % design_rhs(form, motif, samples, unknown_index, q_star, &
         & design, design_direction, time, rhs_buffer)

    call rhs_buffer % get_real(rhs_values)
    if (size(rhs_values) /= n) then
       error stop 'gti_time_local_tangent_driver: rhs size matches unknown size'
    end if

    !----------------------------------------------------------------!
    ! The dense J_q, one exact Newton action per basis column -
    ! the same Jacobian the primal solve used, reused, never
    ! rederived.
    !----------------------------------------------------------------!

    allocate(jacobian(n, n), basis(n))

    do j = 1, n
       basis    = 0.0_dp
       basis(j) = 1.0_dp
       call dq_basis % set_real(basis, ncomp=q_star % ncomp)
       call newton % jacobian_action(form, motif, samples, unknown_index, &
            & q_star, dq_basis, design, time, column)
       call column % get_real(column_values)
       jacobian(:, j) = column_values
    end do

    call dense_solve(jacobian, rhs_values, singular_tolerance, tangent)

    result % linear_residual_norm = norm2(matmul(jacobian, tangent) - rhs_values)
    result % solved = .true.

    call result % q_tangent % set_real(tangent, ncomp=q_star % ncomp)
    result % rhs = rhs_buffer

  end subroutine solve_design_tangent

  !===================================================================!
  ! Solve A x = b by Gaussian elimination with partial pivoting,
  ! on private copies. A pivot below the singular tolerance is
  ! refused. A deliberate local duplicate of the Newton driver's
  ! helper - serving this driver alone, not an abstraction.
  !===================================================================!

  pure subroutine dense_solve(a, b, singular_tolerance, x)

    real(dp)             , intent(in)  :: a(:,:), b(:)
    real(dp)             , intent(in)  :: singular_tolerance
    real(dp), allocatable, intent(out) :: x(:)

    real(dp), allocatable :: m(:,:), r(:), row(:)
    real(dp) :: swap_value, factor
    integer :: n, k, p, i

    n = size(b)
    m = a
    r = b
    allocate(row(n))

    do k = 1, n

       p = k - 1 + maxloc(abs(m(k:n, k)), dim=1)

       if (abs(m(p, k)) <= singular_tolerance) then
          error stop 'gti_time_local_tangent_driver: dense Jacobian pivot is nonsingular'
       end if

       if (p /= k) then
          row        = m(k, :)
          m(k, :)    = m(p, :)
          m(p, :)    = row
          swap_value = r(k)
          r(k)       = r(p)
          r(p)       = swap_value
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

end module gti_time_local_tangent_drivers
