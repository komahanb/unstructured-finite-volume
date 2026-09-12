program gti_contract
  use util_precision, only : dp
  use operation_minimization, only : relative, absolute, by_count, SOLVE_NOT_STARTED
  use operation_family, only : implicit_midpoint, crouzeix_three_stage, bdf_family
  use operation_grid, only : uniform_grid, fixed_grid
  use operation_expression, only : expression
  use gti_physics, only : van_der_pol, van_der_pol_energy
  use gti_expansion, only : expansion, family_container
  use gti_chain, only : chain_block, march_chain, chain_expansion, chain_versions, &
       & chain_derivative, instant_components, grid_stationary_partition, &
       & functional_error, functional_error_estimate
  use gti_march, only : march_context, imbalance, consistent_state, frozen_inputs, solve_linear
  use gti_sweeps, only : reverse_pass
  use gti_adaptive, only : adaptive_partition
  use field_stored, only : stored_field
  implicit none
  type(march_context) :: context
  type(family_container) :: schemes(1), enriched(1)
  type(functional_error_estimate), allocatable :: estimates(:)
  type(chain_block), allocatable :: chain(:)
  type(expansion), allocatable, target :: tower
  type(expression) :: physics, functional(1)
  type(imbalance) :: final_imbalance
  type(stored_field), allocatable :: inputs(:)
  real(dp), allocatable :: initial(:), dt(:), t(:), values(:), table(:,:), rhs(:), direction(:), weights(:)
  integer, allocatable :: versions(:)
  ! the Newton stopping rule of every solve here: relative to the initial residual
  real(dp), parameter :: solver_tolerance = 1.0e-12_dp
  ! implicit midpoint on q'' + q = 0 from (1, 0): one step is a rotation, |y_k| = 1
  ! exactly, the stage value has |Y_k|^2 = 1 / (1 + h_k^2/4), and the one-stage
  ! quadrature gives E_h = sum_k h_k / (2 (1 + h_k^2/4)), one function of each step:
  ! on the uniform grid dE_h/dw_k is the same for every k and the projected grid
  ! gradient is zero, while at T = 2, N = 8, h = 1/4 the value is E_h = 64/65 and
  ! E_h - 1 = -1/65 = -1.5384615384615385e-2, the exact value E = T/2 = 1 missed by 1.5e-2
  integer , parameter :: stationary_steps = 8
  real(dp), parameter :: stationary_duration = 2.0_dp
  real(dp), parameter :: stationary_energy = 64.0_dp / 65.0_dp
  real(dp) :: achieved, errors(3), slope, design_value, stationarity_defect, energy, floor
  integer :: k, n, rejects
  character(len=32) :: mode
  call get_command_argument(1, mode)
  call context % set_time_coupling('coupled')
  call context % set_stopping(solver_tolerance, relative, by_count, 100)
  physics = van_der_pol(2)
  functional = [van_der_pol_energy(2)]
  design_value = 0.8_dp
  if (trim(mode) == 'accuracy' .or. trim(mode) == 'grid_stationary') design_value = 0.0_dp
  initial = consistent_state(physics, 3, [1.0_dp, 0.0_dp], design_value, context=context)
  allocate(schemes(1) % scheme, source=implicit_midpoint())
  select case (trim(mode))
  case ('accuracy')
     do k = 1, 3
        n = 10 * 2**(k - 1) + 1
        call march_chain(schemes, [n], physics, 3, uniform_grid(1.0_dp), &
             & design_value, initial, chain, tower, dt, t, achieved, final_imbalance=final_imbalance, context=context)
        if (.not. final_imbalance % converged) error stop 'accuracy: primal solve must converge'
        values = instant_components(chain, size(dt))
        errors(k) = norm2(values(1:2) - [cos(1.0_dp), -sin(1.0_dp)])
     end do
     do k = 1, 2
        slope = log(errors(k) / errors(k + 1)) / log(2.0_dp)
        if (slope < 1.8_dp .or. slope > 2.2_dp) error stop 'accuracy: midpoint has second order'
     end do
     print *, 'PASS: midpoint second order against the analytic oscillator'
  case ('grid_stationary')
     ! the grid-stationary criterion accepts the uniform seed at the solver tolerance,
     ! far above roundoff, while the functional misses its exact value by 1/65: the
     ! stationarity defect is not an error bound
     weights = grid_stationary_partition(implicit_midpoint(), physics, functional(1), 3, &
          & stationary_duration, [1.0_dp, 0.0_dp], design_value, solver_tolerance, .true., rejects, &
          & stationarity_defect=stationarity_defect, context=context)
     if (rejects /= 0) error stop 'grid_stationary: the uniform seed is accepted at the first attempt'
     if (size(weights) /= stationary_steps) error stop 'grid_stationary: the seed has eight steps'
     if (any(abs(weights - stationary_duration / real(stationary_steps, dp)) > &
          & epsilon(1.0_dp) * stationary_duration)) error stop 'grid_stationary: the accepted grid is uniform'
     if (stationarity_defect > solver_tolerance) error stop 'grid_stationary: the defect is within the tolerance'
     ! one block of N + 1 instants: the first has no step, each weight is one step
     call march_chain(schemes, [stationary_steps + 1], physics, 3, fixed_grid(weights), &
          & design_value, initial, chain, tower, dt, t, achieved, final_imbalance=final_imbalance, context=context)
     if (size(dt) /= stationary_steps + 1 .or. abs(t(size(t)) - stationary_duration) > &
          & epsilon(1.0_dp) * stationary_duration) error stop 'grid_stationary: the march covers T in N steps'
     if (.not. final_imbalance % converged) error stop 'grid_stationary: primal solve must converge'
     call chain_expansion(chain, tower, functional, 3, 0, table, context=context)
     energy = table(0, 1)
     ! the physics at nu = 0 is linear in the state, so Newton converges in one iteration
     ! and the solve's own relative tolerance is the floor on every value it returns
     floor  = solver_tolerance * stationary_energy
     if (abs(energy - stationary_energy) > floor) error stop 'grid_stationary: E_h = 64/65 on the uniform grid'
     values = instant_components(chain, size(dt))
     if (abs(values(1) ** 2 + values(2) ** 2 - 1.0_dp) > solver_tolerance) &
          & error stop 'grid_stationary: the midpoint step is a rotation, |y_N| = 1'
     if (abs(energy - stationary_duration / 2.0_dp) < 1.0e-2_dp) &
          & error stop 'grid_stationary: the accepted grid misses E = T/2 by 1/65'
     write(*,'(a,es9.2,a,es9.2)') ' PASS: grid stationary yet inaccurate: stationarity defect ', &
          & stationarity_defect, ' accepted at ', solver_tolerance
     write(*,'(a,es23.16,a,es9.2,a,es9.2)') '       E_h = ', energy, ' with E_h - 64/65 = ', &
          & energy - stationary_energy, ' and E_h - T/2 = ', energy - stationary_duration / 2.0_dp
     ! the functional-error estimator on the accepted grid, enriched by bdf3 on the
     ! instants: its quadrature part is F+(P Q_h) - F_h = 1 - 64/65 exactly (the instant
     ! jets have (q^2 + q'^2)/2 = 1/2 and the bdf3 weights sum to each step), so the
     ! estimate reports the 1/65 the stationarity check misses, with the sign of E - E_h
     allocate(enriched(1) % scheme, source=bdf_family(3))
     call functional_error(chain, tower, enriched, physics, 3, fixed_grid(weights), design_value, &
          & functional, estimates, context=context)
     if (abs(estimates(1) % quadrature_part - (1.0_dp - stationary_energy)) > floor) &
          & error stop 'grid_stationary: the quadrature part is 1 - 64/65'
     if (abs(estimates(1) % estimate) <= 1.0e-2_dp) &
          & error stop 'grid_stationary: the estimator reports the error the stationarity check misses'
     if (estimates(1) % estimate * (stationary_duration / 2.0_dp - energy) <= 0.0_dp) &
          & error stop 'grid_stationary: the estimate has the sign of E - E_h'
     write(*,'(a,es12.4,a,es12.4,a,f8.4)') '       functional error estimate ', estimates(1) % estimate, &
          & ' against E - E_h = ', stationary_duration / 2.0_dp - energy, ' effectivity ', &
          & estimates(1) % estimate / (stationary_duration / 2.0_dp - energy)
  case ('adaptive_failure')
     call context % set_stopping(1.0e-14_dp, relative, by_count, 1)
     dt = adaptive_partition(crouzeix_three_stage(), 4, physics, 3, 1.0_dp, &
          & [1.0_dp, 0.0_dp], design_value, 1.0e-6_dp, .false., context=context)
  case ('minimum_step')
     dt = adaptive_partition(implicit_midpoint(), 2, physics, 3, 1.0_dp, &
          & [1.0_dp, 0.0_dp], design_value, 1.0e-14_dp, .false., minimum_step=0.125_dp, context=context)
  case ('status', 'forward', 'reverse')
     if (trim(mode) == 'status') then
        final_imbalance = imbalance()
        if (final_imbalance % converged .or. final_imbalance % outcome % converged()) &
             & error stop 'status: an uncomputed imbalance is unconverged'
        if (final_imbalance % outcome % reason /= SOLVE_NOT_STARTED) &
             & error stop 'status: an uncomputed imbalance is unstarted'
        print *, 'PASS: an uncomputed imbalance is unstarted and unconverged'
     end if
     call context % set_stopping(1.0e-14_dp, relative, by_count, 1)
     call march_chain(schemes, [11], physics, 3, uniform_grid(2.0_dp), &
          & design_value, initial, chain, tower, dt, t, achieved, final_imbalance=final_imbalance, context=context)
     if (final_imbalance % converged) error stop 'status: under-iterated primal was accepted'
     if (final_imbalance % outcome % iterations /= 1) error stop 'status: count completed iterations'
     if (abs(achieved - final_imbalance % norm) > 1.0e-14_dp) error stop 'status: report block residual'
     if (trim(mode) == 'status') then
        print *, 'PASS: unsuccessful primal result survives driver rule copies'
     else if (trim(mode) == 'forward') then
        call chain_expansion(chain, tower, functional, 3, 1, table, context=context)
     else
        call chain_versions(chain, tower, functional, 3, versions, context=context)
        call chain_derivative(chain, tower, versions, functional, 3, 1, reverse_pass, table, context=context)
     end if
  case ('linear_forward', 'linear_reverse')
     call march_chain(schemes, [11], physics, 3, uniform_grid(1.0_dp), &
          & design_value, initial, chain, tower, dt, t, achieved, final_imbalance=final_imbalance, context=context)
     if (.not. final_imbalance % converged) error stop 'linear: primal solve must converge'
     call frozen_inputs(chain(1) % rows, chain(1) % state, design_value, inputs)
     rhs = [(real(k, dp), k = 1, size(chain(1) % state))]
     call context % set_linear_solver('iterative')
     call context % set_storage('sparse')
     call context % set_linear_limits(1, 1, 1)
     call context % set_stopping(1.0e-14_dp, absolute, by_count, 1)
     call solve_linear(chain(1) % rows, inputs, rhs, trim(mode) == 'linear_reverse', 1, direction, &
         & context=context)
  case default
     error stop 'gti_contract: an existing case is selected'
  end select
end program gti_contract
