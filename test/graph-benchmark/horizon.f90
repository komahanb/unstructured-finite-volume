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
! The reverse mode is a streamed reverse execution of the derivatives
! to max_order within the storage limit admitting restart_count stored
! restart states (retained: no limit), over one of three chains: crouzeix, the
! one-instant blocks above; mixed, cycles of bdf 3 / adams 3 /
! crouzeix / newmark blocks (4, 3, 1, 1 instants) over a startup of
! four fine steps; designed, bdf 2 blocks of three instants over the
! same startup on a designed grid with the step weights as a second
! design. The march
! and the reverse phase are recorded separately with the tally's Newton
! and linear solves of each; the accounts, the schedule's quantities,
! the peak resident size before the references, and the table to 17
! digits are the summary. The references: the post-hoc reverse pass
! over a retained primal, equal bitwise, and the forward expansion
! along the physics' parameter at tau^(2/3).
!
! usage: horizon n mode max_order tau [restart_count] [chain]
!=====================================================================!

program horizon_scaling

  use util_precision        , only : dp
  use operation_minimization, only : relative, by_rate
  use operation_family      , only : crouzeix_three_stage, bdf_family, adams_family, newmark_family
  use operation_grid        , only : uniform_grid, designed_grid
  use operation_expression  , only : expression
  use util_tally            , only : tally_open, tally_close, tally_amount, tally_num_levels, newton_solves, linear_solves
  use gti_configuration     , only : hierarchy_levels
  use gti_physics           , only : van_der_pol, van_der_pol_energy, van_der_pol_dissipation
  use gti_expansion         , only : expansion, family_container
  use gti_march             , only : march_context, imbalance, consistent_state
  use gti_sweeps            , only : reverse_pass
  use gti_chain             , only : chain_execution, chain_block, chain_expansion, chain_versions, chain_derivative, &
       &                             derivative_storage
  use benchmark_measurement , only : phase, begin_phase, end_phase, record, peak_rss_kilobytes, &
       &                             argument_integer, argument_string, argument_real

  implicit none

  integer, parameter :: state_degree = 2, degrees = state_degree + 1
  real(dp), parameter :: duration = 8.0_dp, design = 0.8_dp
  character(len=:), allocatable :: mode, tokens, restart_count, chain_kind
  character(len=160) :: written
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
  restart_count    = argument_string(5, 'retained')
  chain_kind = argument_string(6, 'crouzeix')
  if (n < 2) error stop 'horizon: at least two instants'
  if (mode == 'reverse') then
     write(written, '(a,i0,a,a,a,i0,a,a,a,a)') 'suite=reverse n=', n, ' mode=', mode, ' max_order=', max_order, &
          & ' restart_count=', restart_count, ' chain=', chain_kind
  else
     write(written, '(a,i0,a,a,a,i0)') 'suite=horizon n=', n, ' mode=', mode, ' max_order=', max_order
  end if
  tokens = trim(written)
  num_failures = 0
  tower_storage = 0
  state_storage = 0

  call context % set_stopping(tau, relative, by_rate, 100)
  agreement = tau ** (2.0_dp / 3.0_dp)
  functionals(1) = van_der_pol_energy(state_degree)
  functionals(2) = van_der_pol_dissipation(state_degree)
  q0 = consistent_state(van_der_pol(state_degree), degrees, [1.0_dp, 0.0_dp], design, context=context)

  if (mode == 'reverse') then
     call reverse_horizon()
     if (num_failures > 0) error stop 'horizon: a verification failed'
     stop
  end if

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
     error stop 'horizon: the mode is primal, taylor or reverse'
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

  !-------------------------------------------------------------------!
  ! The chain of one kind over about n instants: the families, the
  ! instants each block adds, the startup refinement, the grid and its
  ! design, the designs the reverse runs over.
  !-------------------------------------------------------------------!
  subroutine chain_of(kind, added, startup, weights, designs, instants)
    character(len=*), intent(in) :: kind
    integer, allocatable, intent(out) :: added(:)
    integer, intent(out) :: startup, designs, instants
    real(dp), allocatable, intent(out) :: weights(:)
    integer :: k, cycles, blocks
    startup = 0
    designs = 1
    weights = [real(dp) ::]
    select case (kind)
    case ('crouzeix')
       added = [2, (1, b = 2, n - 1)]
       allocate(schemes(n - 1))
       do b = 1, n - 1
          allocate(schemes(b) % scheme, source=crouzeix_three_stage())
       end do
       instants = n
    case ('mixed')
       cycles = max(1, (n - 1) / 9)
       added = [5, 3, 1, 1, (4, 3, 1, 1, k = 2, cycles)]
       allocate(schemes(4 * cycles))
       do b = 1, 4 * cycles, 4
          allocate(schemes(b) % scheme, source=bdf_family(3))
          allocate(schemes(b + 1) % scheme, source=adams_family(3))
          allocate(schemes(b + 2) % scheme, source=crouzeix_three_stage())
          allocate(schemes(b + 3) % scheme, source=newmark_family(0.25_dp, 0.5_dp))
       end do
       startup = 4
       instants = sum(added)
    case ('designed')
       blocks = max(1, n / 3)
       added = [(3, b = 1, blocks)]
       allocate(schemes(blocks))
       do b = 1, blocks
          allocate(schemes(b) % scheme, source=bdf_family(2))
       end do
       instants = 3 * blocks
       weights = [(1.0_dp + 0.05_dp * real(mod(k, 3), dp), k = 1, instants - 1)]
       designs = 2
       startup = 4
    case default
       error stop 'horizon: the chain is crouzeix, mixed or designed'
    end select
  end subroutine chain_of

  subroutine started(added, startup, weights, designs)
    integer, intent(in) :: added(:), startup, designs
    real(dp), intent(in) :: weights(:)
    if (size(weights) > 0) then
       call execution % initialize(schemes, added, van_der_pol(state_degree), degrees, designed_grid(duration), &
            & design, q0, grid_design=weights, startup=startup, functionals=functionals, derivative_order=max_order, &
            & context=context, pass_kind=reverse_pass, designs=designs)
    else if (startup > 0) then
       call execution % initialize(schemes, added, van_der_pol(state_degree), degrees, uniform_grid(duration), &
            & design, q0, startup=startup, functionals=functionals, derivative_order=max_order, &
            & context=context, pass_kind=reverse_pass, designs=designs)
    else
       call execution % initialize(schemes, added, van_der_pol(state_degree), degrees, uniform_grid(duration), &
            & design, q0, functionals=functionals, derivative_order=max_order, &
            & context=context, pass_kind=reverse_pass, designs=designs)
    end if
  end subroutine started

  subroutine retained_primal(added, startup, weights)
    integer, intent(in) :: added(:), startup
    real(dp), intent(in) :: weights(:)
    if (size(weights) > 0) then
       call execution % initialize(schemes, added, van_der_pol(state_degree), degrees, designed_grid(duration), &
            & design, q0, grid_design=weights, startup=startup, context=context)
    else if (startup > 0) then
       call execution % initialize(schemes, added, van_der_pol(state_degree), degrees, uniform_grid(duration), &
            & design, q0, startup=startup, context=context)
    else
       call execution % initialize(schemes, added, van_der_pol(state_degree), degrees, uniform_grid(duration), &
            & design, q0, context=context)
    end if
  end subroutine retained_primal

  real(dp) function counted(event)
    integer, intent(in) :: event
    integer :: level, order
    counted = 0.0_dp
    do level = 1, tally_num_levels()
       do order = 0, max_order + 1
          counted = counted + tally_amount(level, order, event)
       end do
    end do
  end function counted

  subroutine reverse_horizon()
    type(derivative_storage) :: quantities, accounts
    integer, allocatable :: added(:), versions(:)
    real(dp), allocatable :: weights(:), table(:,:), reference(:,:)
    real(dp) :: before_newton, before_linear, march_newton, march_linear, reverse_newton, reverse_linear
    integer :: startup, designs, instants, c, rss, i
    integer(8) :: limit
    call chain_of(chain_kind, added, startup, weights, designs, instants)
    call tally_open(max_order + 1, hierarchy_levels)
    ! the schedule's quantities, read from an initialization without a limit
    call started(added, startup, weights, designs)
    quantities = execution % streamed_storage()
    limit = huge(1)
    if (restart_count /= 'retained') then
       read(restart_count, *) c
       limit = quantities % restart + quantities % leaf + quantities % window + quantities % terms &
            & + int(c, 8) * quantities % restart
       if (limit > huge(1)) error stop 'horizon: the limit is an entry count'
       call context % set_reverse_limit(int(limit))
    end if
    call begin_phase(measured)
    call started(added, startup, weights, designs)
    call end_phase(measured)
    call record(tokens, 'initialize', measured)
    advances = 0
    before_newton = counted(newton_solves)
    before_linear = counted(linear_solves)
    call begin_phase(measured)
    do while (.not. execution % complete())
       call execution % advance()
       advances = advances + 1
    end do
    call end_phase(measured)
    call record(tokens, 'advance', measured)
    march_newton = counted(newton_solves) - before_newton
    march_linear = counted(linear_solves) - before_linear
    before_newton = counted(newton_solves)
    before_linear = counted(linear_solves)
    call begin_phase(measured)
    call execution % derivative(functionals, max_order, reverse_pass, table, designs=designs, storage=accounts)
    call end_phase(measured)
    call record(tokens, 'reverse', measured)
    reverse_newton = counted(newton_solves) - before_newton
    reverse_linear = counted(linear_solves) - before_linear
    ! the peak resident size of the bounded execution, before the references
    rss = int(peak_rss_kilobytes())
    call verified(accounts % peak <= accounts % limit, 'peak_within_limit')
    if (restart_count == 'retained' .or. limit >= accounts % retained) then
       call verified(accounts % checkpoints == 0 .and. accounts % evaluations == 0, 'retention_recomputes_nothing')
    else
       call verified(accounts % checkpoints == c, 'the_limit_admits_the_restart_count')
    end if
    ! the reference: the post-hoc reverse pass over a retained primal
    call begin_phase(measured)
    call retained_primal(added, startup, weights)
    do while (.not. execution % complete())
       call execution % advance()
    end do
    call execution % take_results(chain, tower, dt, t, achieved, final_imbalance=final_imbalance)
    call chain_versions(chain, tower, functionals, degrees, versions, context=context)
    call chain_derivative(chain, tower, versions, functionals, degrees, max_order, reverse_pass, reference, &
         & designs=designs, context=context)
    call end_phase(measured)
    call record(tokens, 'reference', measured)
    call verified(final_imbalance % converged, 'retained_march_converged')
    call verified(all(shape(table) == shape(reference)) .and. all(table == reference), &
         & 'bounded_reverse_equals_retained_post_hoc_bitwise')
    call chain_expansion(chain, tower, functionals, degrees, max_order, expanded, context=context)
    departure = maxval(abs(table(:, 1) - expanded(max_order, :)) / max(1.0_dp, abs(expanded(max_order, :))))
    call verified(departure <= agreement, 'reverse_derivatives_agree_with_forward_expansion')
    call tally_close()
    write(*, '(a,1x,a,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0)', &
         & advance='no') 'summary', tokens, 'instants=', instants, 'blocks=', size(chain), 'advances=', advances, &
         & 'checkpoints=', accounts % checkpoints, 'evaluations=', accounts % evaluations, 'limit=', accounts % limit, &
         & 'restart=', accounts % restart, 'leaf=', accounts % leaf, 'window=', accounts % window, &
         & 'terms=', accounts % terms, 'retained=', accounts % retained, 'peak=', accounts % peak
    write(*, '(1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0)', advance='no') &
         & 'state_high=', accounts % state % high, 'tower_high=', accounts % tower % high, &
         & 'costate_high=', accounts % costate % high, 'checkpoint_high=', accounts % checkpoint % high, &
         & 'scalars_high=', accounts % scalars % high, 'march_newton_solves=', nint(march_newton), &
         & 'march_linear_solves=', nint(march_linear), 'reverse_newton_solves=', nint(reverse_newton), &
         & 'reverse_linear_solves=', nint(reverse_linear)
    write(*, '(1x,a,es12.4e3,1x,a,es12.4e3,1x,a,i0)', advance='no') 'departure=', departure, 'agreement=', agreement, &
         & 'peak_rss_kilobytes=', rss
    do i = 1, size(functionals)
       write(*, '(1x,a,i0,a,es25.17e3)', advance='no') 'table_', i, '=', table(i, 1)
    end do
    write(*, '(a)') ''
  end subroutine reverse_horizon

end program horizon_scaling
