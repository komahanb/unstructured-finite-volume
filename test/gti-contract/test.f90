program gti_contract
  use util_precision, only : dp
  use operation_minimization, only : relative, absolute, by_count, SOLVE_NOT_STARTED
  use operation_family, only : implicit_midpoint, crouzeix_three_stage
  use operation_grid, only : uniform_grid
  use operation_expression, only : expression
  use gti_physics, only : van_der_pol, van_der_pol_energy
  use gti_expansion, only : expansion, family_container
  use gti_chain, only : chain_block, march_chain, chain_expansion, chain_versions, &
       & chain_derivative, instant_components
  use gti_march, only : imbalance, consistent_state, set_stopping, set_time_coupling, &
       & frozen_inputs, solve_linear
  use gti_sweeps, only : reverse_pass, set_linear_solver, set_storage, set_linear_limits
  use gti_adaptive, only : adaptive_partition
  use view_directed_stored, only : stored_directed_graph
  use field_stored, only : stored_field
  implicit none
  type(family_container) :: schemes(1)
  type(chain_block), allocatable :: chain(:)
  type(expansion), allocatable, target :: tower
  type(expression) :: physics, functional(1)
  type(imbalance) :: final_imbalance
  type(stored_directed_graph) :: domain
  type(stored_field), allocatable :: inputs(:)
  real(dp), allocatable :: initial(:), dt(:), t(:), values(:), table(:,:), rhs(:), direction(:)
  integer, allocatable :: versions(:)
  real(dp) :: achieved, errors(3), slope, design_value
  integer :: k, n
  character(len=32) :: mode
  call get_command_argument(1, mode)
  call set_time_coupling('coupled')
  call set_stopping(1.0e-12_dp, relative, by_count, 100)
  physics = van_der_pol(2)
  functional = [van_der_pol_energy(2)]
  design_value = 0.8_dp
  if (trim(mode) == 'accuracy') design_value = 0.0_dp
  initial = consistent_state(physics, 3, [1.0_dp, 0.0_dp], design_value)
  allocate(schemes(1) % scheme, source=implicit_midpoint())
  select case (trim(mode))
  case ('accuracy')
     do k = 1, 3
        n = 10 * 2**(k - 1) + 1
        call march_chain(schemes, [n], physics, 3, uniform_grid(1.0_dp), &
             & design_value, initial, chain, tower, dt, t, achieved, final_imbalance=final_imbalance)
        if (.not. final_imbalance % converged) error stop 'accuracy: primal solve must converge'
        values = instant_components(chain, size(dt))
        errors(k) = norm2(values(1:2) - [cos(1.0_dp), -sin(1.0_dp)])
     end do
     do k = 1, 2
        slope = log(errors(k) / errors(k + 1)) / log(2.0_dp)
        if (slope < 1.8_dp .or. slope > 2.2_dp) error stop 'accuracy: midpoint has second order'
     end do
     print *, 'PASS: midpoint second order against the analytic oscillator'
  case ('adaptive_failure')
     call set_stopping(1.0e-14_dp, relative, by_count, 1)
     dt = adaptive_partition(crouzeix_three_stage(), 4, physics, 3, 1.0_dp, &
          & [1.0_dp, 0.0_dp], design_value, 1.0e-6_dp, .false.)
  case ('minimum_step')
     dt = adaptive_partition(implicit_midpoint(), 2, physics, 3, 1.0_dp, &
          & [1.0_dp, 0.0_dp], design_value, 1.0e-14_dp, .false., minimum_step=0.125_dp)
  case ('status', 'forward', 'reverse')
     if (trim(mode) == 'status') then
        final_imbalance = imbalance()
        if (final_imbalance % converged .or. final_imbalance % outcome % converged()) &
             & error stop 'status: an uncomputed imbalance is unconverged'
        if (final_imbalance % outcome % reason /= SOLVE_NOT_STARTED) &
             & error stop 'status: an uncomputed imbalance is unstarted'
        print *, 'PASS: an uncomputed imbalance is unstarted and unconverged'
     end if
     call set_stopping(1.0e-14_dp, relative, by_count, 1)
     call march_chain(schemes, [11], physics, 3, uniform_grid(2.0_dp), &
          & design_value, initial, chain, tower, dt, t, achieved, final_imbalance=final_imbalance)
     if (final_imbalance % converged) error stop 'status: under-iterated primal was accepted'
     if (final_imbalance % outcome % iterations /= 1) error stop 'status: count completed iterations'
     if (abs(achieved - final_imbalance % norm) > 1.0e-14_dp) error stop 'status: report block residual'
     if (trim(mode) == 'status') then
        print *, 'PASS: unsuccessful primal result survives driver rule copies'
     else if (trim(mode) == 'forward') then
        call chain_expansion(chain, tower, functional, 3, 1, table)
     else
        call chain_versions(chain, tower, functional, 3, versions)
        call chain_derivative(chain, tower, versions, functional, 3, 1, reverse_pass, table)
     end if
  case ('linear_forward', 'linear_reverse')
     call march_chain(schemes, [11], physics, 3, uniform_grid(1.0_dp), &
          & design_value, initial, chain, tower, dt, t, achieved, final_imbalance=final_imbalance)
     if (.not. final_imbalance % converged) error stop 'linear: primal solve must converge'
     call frozen_inputs(chain(1) % state, design_value, chain(1) % rows % num_points(), domain, inputs)
     rhs = [(real(k, dp), k = 1, size(chain(1) % state))]
     call set_linear_solver('iterative')
     call set_storage('sparse')
     call set_linear_limits(1, 1, 1)
     call set_stopping(1.0e-14_dp, absolute, by_count, 1)
     call solve_linear(chain(1) % rows, domain, inputs, rhs, trim(mode) == 'linear_reverse', 1, direction)
  case default
     error stop 'gti_contract: an existing case is selected'
  end select
end program gti_contract
