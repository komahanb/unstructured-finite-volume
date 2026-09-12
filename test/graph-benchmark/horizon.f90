!=====================================================================!
! One complete GTI consumer over a growing horizon: van der Pol under
! the Crouzeix three-stage DIRK, one instant per block, so the number
! of blocks equals the number of instants and the dependency schedule
! grows with the horizon. Initialization (dependency incidence and the
! driver schedule), advance until complete, and taking the results are
! measured separately; the Taylor mode streams the functionals'
! derivatives to max_order and reports the live and total tower and
! state storage. The oracle is the same horizon marched as one block
! and expanded afterwards: the demo's agreement law tau^(2/3).
!
! usage: horizon n mode max_order tau
!=====================================================================!

program horizon_scaling

  use util_precision        , only : dp
  use operation_minimization, only : relative, by_rate
  use operation_family      , only : crouzeix_three_stage
  use operation_grid        , only : uniform_grid
  use operation_expression  , only : expression
  use gti_physics           , only : van_der_pol, van_der_pol_energy, van_der_pol_dissipation
  use gti_expansion         , only : expansion, family_container
  use gti_march             , only : march_context, imbalance, consistent_state
  use gti_chain             , only : chain_execution, chain_block, chain_expansion
  use benchmark_measurement , only : phase, begin_phase, end_phase, record, peak_rss_kilobytes, &
       &                             argument_integer, argument_string, argument_real

  implicit none

  integer, parameter :: state_degree = 2, degrees = state_degree + 1
  real(dp), parameter :: duration = 8.0_dp, design = 0.8_dp
  character(len=:), allocatable :: mode, tokens
  character(len=96) :: written
  type(march_context) :: context
  type(chain_execution) :: execution
  type(family_container), allocatable :: schemes(:)
  type(expression) :: functionals(2)
  type(chain_block), allocatable :: chain(:), whole_chain(:)
  type(expansion), allocatable, target :: tower, whole_tower
  type(imbalance) :: final_imbalance
  real(dp), allocatable :: q0(:), dt(:), t(:), streamed(:,:), expanded(:,:), final_state(:)
  real(dp) :: tau, agreement, achieved, departure
  integer :: n, max_order, b, advances, num_failures, tower_storage(2), state_storage(2)
  type(phase) :: measured

  n         = argument_integer(1, 40)
  mode      = argument_string(2, 'primal')
  max_order = argument_integer(3, 4)
  tau       = argument_real(4, 1.0e-12_dp)
  if (n < 2) error stop 'horizon: at least two instants'
  write(written, '(a,i0,a,a,a,i0)') 'suite=horizon n=', n, ' mode=', mode, ' max_order=', max_order
  tokens = trim(written)
  num_failures = 0
  tower_storage = 0
  state_storage = 0

  call context % set_stopping(tau, relative, by_rate, 100)
  agreement = tau ** (2.0_dp / 3.0_dp)
  functionals(1) = van_der_pol_energy(state_degree)
  functionals(2) = van_der_pol_dissipation(state_degree)
  q0 = consistent_state(van_der_pol(state_degree), degrees, [1.0_dp, 0.0_dp], design, context=context)

  allocate(schemes(n - 1))
  do b = 1, n - 1
     allocate(schemes(b) % scheme, source=crouzeix_three_stage())
  end do

  call begin_phase(measured)
  select case (mode)
  case ('primal')
     call execution % initialize(schemes, [2, (1, b = 2, n - 1)], van_der_pol(state_degree), degrees, &
          & uniform_grid(duration), design, q0, context=context)
  case ('taylor')
     call execution % initialize(schemes, [2, (1, b = 2, n - 1)], van_der_pol(state_degree), degrees, &
          & uniform_grid(duration), design, q0, functionals=functionals, derivative_order=max_order, &
          & context=context)
  case default
     error stop 'horizon: the mode is primal or taylor'
  end select
  call end_phase(measured)
  call record(tokens, 'initialize', measured)

  advances = 0
  call begin_phase(measured)
  do while (.not. execution % complete())
     call execution % advance()
     advances = advances + 1
  end do
  call end_phase(measured)
  call record(tokens, 'advance', measured)

  call begin_phase(measured)
  if (mode == 'taylor') then
     call execution % take_results(chain, tower, dt, t, achieved, final_imbalance=final_imbalance, &
          & f=streamed, tower_storage=tower_storage, state_storage=state_storage)
  else
     call execution % take_results(chain, tower, dt, t, achieved, final_imbalance=final_imbalance)
  end if
  call end_phase(measured)
  call record(tokens, 'results', measured)
  call verified(final_imbalance % converged, 'pipelined_march_converged')
  call verified(advances == n - 1, 'one_advance_per_block')

  ! the oracle: the same horizon as one block
  deallocate(schemes)
  allocate(schemes(1))
  allocate(schemes(1) % scheme, source=crouzeix_three_stage())
  call begin_phase(measured)
  call execution % initialize(schemes, [n], van_der_pol(state_degree), degrees, uniform_grid(duration), &
       & design, q0, context=context)
  do while (.not. execution % complete())
     call execution % advance()
  end do
  call execution % take_results(whole_chain, whole_tower, dt, t, achieved, final_imbalance=final_imbalance)
  call end_phase(measured)
  call record(tokens, 'whole_horizon', measured)
  call verified(final_imbalance % converged, 'whole_march_converged')

  if (mode == 'taylor') then
     call begin_phase(measured)
     call chain_expansion(whole_chain, whole_tower, functionals, degrees, max_order, expanded, context=context)
     call end_phase(measured)
     call record(tokens, 'whole_expansion', measured)
     departure = maxval(abs(streamed - expanded) / max(1.0_dp, abs(expanded)))
     call verified(all(shape(streamed) == shape(expanded)) .and. departure <= agreement, &
          & 'streamed_derivatives_agree_with_whole_horizon')
  else
     final_state = chain(size(chain)) % state
     departure = maxval(abs(final_state - whole_chain(1) % state(size(whole_chain(1) % state) - size(final_state) + 1:)) &
          & / max(1.0_dp, abs(final_state)))
     call verified(departure <= agreement, 'final_state_agrees_with_whole_horizon')
  end if

  write(*, '(a,1x,a,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,es12.4e3,1x,a,es12.4e3,1x,a,i0)') &
       & 'summary', tokens, 'blocks=', n - 1, 'advances=', advances, 'tower_live=', tower_storage(1), &
       & 'tower_total=', tower_storage(2), 'state_live=', state_storage(1), 'state_total=', state_storage(2), &
       & 'departure=', departure, 'agreement=', agreement, 'peak_rss_kilobytes=', peak_rss_kilobytes()
  if (num_failures > 0) error stop 'horizon: a verification failed'

contains

  subroutine verified(condition, name)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: name
    if (.not. condition) num_failures = num_failures + 1
    write(*, '(a,1x,a,1x,a,a,1x,a,l1)') 'verification', tokens, 'name=', name, 'satisfied=', condition
  end subroutine verified

end program horizon_scaling
