program execution_independence

  use iso_fortran_env, only : int64
  use util_precision, only : dp
  use operation_minimization, only : relative, by_count
  use operation_family, only : bdf_family, adams_family, crouzeix_three_stage, newmark_family
  use operation_grid, only : grid, uniform_grid, designed_grid
  use util_tally, only : tally, linear_solves, newton_solves
  use gti_configuration, only : hierarchy_levels
  use operation_expression, only : expression
  use gti_physics, only : van_der_pol, van_der_pol_energy
  use gti_expansion, only : expansion, family_container
  use gti_march, only : march_context, imbalance, consistent_state
  use gti_sweeps, only : forward_pass, reverse_pass
  use gti_chain, only : chain_execution, chain_block, derivative_storage, multiset_count

  implicit none

  type :: execution_result
     real(dp), allocatable :: states(:), steps(:), times(:)
     real(dp), allocatable :: forward_derivative(:,:), reverse_derivative(:,:), functional(:,:)
     integer :: tower_storage(2) = 0, state_storage(2) = 0
  end type execution_result

  ! A containing object: its copy copies the execution it contains.
  type :: execution_holder
     type(chain_execution) :: execution
     integer :: label = 0
  end type execution_holder

  type(march_context) :: context_a, context_b
  type(chain_execution) :: execution_a, execution_b
  type(chain_execution), allocatable :: twin
  type(execution_result) :: reference_a, reference_b, interleaved_a, interleaved_b, repeated
  type(execution_result) :: stream_a, stream_b, interleaved_stream_a, interleaved_stream_b
  integer :: advances_a, advances_b, order
  character(len=32) :: mode

  call configure(context_a, 1)
  call configure(context_b, 2)
  call get_command_argument(1, mode)
  select case (trim(mode))
  case ('source_twin')
     call initialize_case(execution_a, 1, context_a)
     allocate(twin, source=execution_a)
     deallocate(twin)
     call execution_a % advance()
  case ('advance_uninitialized')
     call execution_a % advance()
  case ('results_incomplete')
     call initialize_case(execution_a, 1, context_a)
     call extract(execution_a, repeated)
  case ('derivative_incomplete')
     call initialize_case(execution_a, 1, context_a)
     call differentiate_order(execution_a, 2, 1, forward_pass, repeated)
  case ('derivative_streamed')
     call initialize_case(execution_a, 1, context_a, streamed=.true.)
     call finish(execution_a)
     call differentiate_order(execution_a, 2, 1, forward_pass, repeated)
  case ('reverse_limit_insufficient')
     call refused_limit()
  case ('reverse_streamed_forward_pass')
     call initialize_case(execution_a, 1, context_a, reverse_order=2)
     call finish(execution_a)
     call differentiate_order(execution_a, 2, 2, forward_pass, repeated)
  case ('reverse_streamed_results')
     call initialize_case(execution_a, 1, context_a, reverse_order=2)
     call finish(execution_a)
     call extract(execution_a, repeated, streamed=.true.)
  end select
  if (len_trim(mode) > 0) error stop 'execution_independence: invalid execution transition was accepted'
  call isolated(1, context_a, reference_a)
  call isolated(2, context_b, reference_b)

  call initialize_case(execution_a, 1, context_a)
  call initialize_case(execution_b, 2, context_b)
  call check(.not. execution_a % complete() .and. .not. execution_b % complete(), &
       & 'initialization leaves the primal rules unevaluated')
  ! Altering the input configurations cannot change either owned execution.
  call context_a % set_stopping(0.2_dp, relative, by_count, 1)
  call context_b % set_stopping(0.3_dp, relative, by_count, 1)
  advances_a = 0
  advances_b = 0
  do while (.not. execution_a % complete() .or. .not. execution_b % complete())
     if (.not. execution_a % complete()) then
        call execution_a % advance()
        advances_a = advances_a + 1
     end if
     if (.not. execution_b % complete()) then
        call execution_b % advance()
        advances_b = advances_b + 1
     end if
     call check(advances_a + advances_b <= 10, 'each advance evaluates one scheduled block')
  end do
  call check(advances_a == 4 .and. advances_b == 2, &
       & 'the startup block participates in the same incremental schedule')

  ! Alternate derivative work as well as primal work, including both orientations.
  do order = 1, 2
     call differentiate_order(execution_a, 2, order, forward_pass, interleaved_a)
     call differentiate_order(execution_b, 3, order, forward_pass, interleaved_b)
     call differentiate_order(execution_a, 2, order, reverse_pass, interleaved_a)
     call differentiate_order(execution_b, 3, order, reverse_pass, interleaved_b)
  end do
  call extract(execution_a, interleaved_a)
  call extract(execution_b, interleaved_b)
  call compare(interleaved_a, reference_a)
  call compare(interleaved_b, reference_b)
  print *, 'PASS: different physics, dimensions, schemes and solver settings interleave independently'
  print *, 'PASS: first and second forward and reverse derivatives survive interleaved execution'

  ! Reinitialize a partially evaluated execution with another problem.
  call configure(context_a, 1)
  call configure(context_b, 2)
  call initialize_case(execution_a, 2, context_b)
  call execution_a % advance()
  call initialize_case(execution_a, 1, context_a)
  call finish(execution_a)
  call differentiate(execution_a, 2, repeated)
  call extract(execution_a, repeated)
  call compare(repeated, reference_a)
  print *, 'PASS: reinitialization removes an unfinished trajectory and its solver state'

  call initialize_case(execution_a, 1, context_a, streamed=.true.)
  call finish(execution_a)
  call extract(execution_a, stream_a, streamed=.true.)
  call initialize_case(execution_b, 2, context_b, streamed=.true.)
  call finish(execution_b)
  call extract(execution_b, stream_b, streamed=.true.)

  call initialize_case(execution_a, 1, context_a, streamed=.true.)
  call initialize_case(execution_b, 2, context_b, streamed=.true.)
  do while (.not. execution_a % complete() .or. .not. execution_b % complete())
     if (.not. execution_b % complete()) call execution_b % advance()
     if (.not. execution_a % complete()) call execution_a % advance()
  end do
  call extract(execution_a, interleaved_stream_a, streamed=.true.)
  call extract(execution_b, interleaved_stream_b, streamed=.true.)
  call compare_stream(interleaved_stream_a, stream_a, reference_a)
  call compare_stream(interleaved_stream_b, stream_b, reference_b)
  print *, 'PASS: streamed Taylor coefficients and live storage are independent across executions'

  call initialize_case(execution_a, 2, context_b, streamed=.true.)
  call execution_a % advance()
  call initialize_case(execution_a, 1, context_a, streamed=.true.)
  call finish(execution_a)
  call extract(execution_a, repeated, streamed=.true.)
  call compare_stream(repeated, stream_a, reference_a)
  print *, 'PASS: reinitialization removes an unfinished Taylor expansion'

  call check_copies()
  call check_bounded_reverse()

contains

  !-------------------------------------------------------------------!
  ! Bounded reverse storage. A reverse derivative under a storage limit
  ! recomputes the tangent towers from restart states: the tables of
  ! every order, the lower orders read from the same products and the
  ! entries equal the retained pass bitwise; the peak of the accounts
  ! is within the limit; the forward block evaluations equal the count
  ! of the schedule's recurrence for the chosen restart states; the
  ! extra linear solves are the recomputed towers alone. Retention: a
  ! limit at the retained peak stores everything and evaluates each
  ! block once. A streamed reverse execution stores the restart states
  ! in its march, states included, and recomputes the primal with the
  ! towers: the same tables bitwise, the recomputations counted as
  ! Newton solves.
  !-------------------------------------------------------------------!
  subroutine check_bounded_reverse()
    integer :: case_index
    call check(recurrence_evaluations(2, 0) == 3 .and. recurrence_evaluations(4, 0) == 10 .and. &
         & recurrence_evaluations(2, 1) == 2 .and. recurrence_evaluations(4, 1) == 6 .and. &
         & recurrence_evaluations(16, 15) == 16 .and. recurrence_evaluations(16, 3) <= recurrence_evaluations(16, 2), &
         & 'the recurrence reproduces the counts computed by hand')
    do case_index = 1, 6
       call bounded_case(case_index)
    end do
    print *, 'PASS: retained and recomputed reverse derivatives agree bitwise over limits, orders, families and designs'
    print *, 'PASS: a streamed reverse execution recomputes states and towers from the restart states of its march'
  end subroutine check_bounded_reverse

  subroutine bounded_case(case_index)
    integer, intent(in) :: case_index
    type(march_context) :: context, limited
    type(chain_execution) :: execution, other, copy
    type(expression) :: functional(1)
    real(dp), allocatable :: table(:,:), by_order(:,:,:), entries(:,:,:)
    real(dp), allocatable :: whole_table(:,:), whole_by_order(:,:,:), whole_entries(:,:,:)
    type(derivative_storage) :: whole, bounded, streamed
    real(dp) :: before, solves_retained, solves_bounded
    integer :: order, c, blocks, top, nd, designs, m, expected, exercised
    integer(int64) :: limit
    character(len=8) :: label
    write(label, '(a,i0)') 'case ', case_index
    call configure(context, solver_of(case_index))
    call initialize_case(execution, case_index, context)
    call finish(execution)
    functional = [van_der_pol_energy(degree_of(case_index))]
    designs = merge(2, 1, case_index == 4)
    nd = designs
    do order = 1, 3
       top = order - 1
       before = counted_solves(execution % account(), linear_solves)
       call execution % derivative(functional, order, reverse_pass, whole_table, by_order=whole_by_order, &
            & entries=whole_entries, designs=designs, storage=whole)
       solves_retained = counted_solves(execution % account(), linear_solves) - before
       blocks = whole % evaluations
       call check(whole % checkpoints == 0 .and. whole % peak <= whole % retained .and. whole % limit == huge(1), &
            & label // ': without a limit every tower is retained and the peak is within the retained account')
       call check(order == 1 .eqv. whole % restart == 0, label // ': order one stores no tower')
       m = 0
       do top = 1, order - 1
          m = m + multiset_count(nd, top)
       end do
       top = order - 1
       if (order > 1) then
          ! the bound at every count of stored restart states from none up
          ! to every position; a limit reaching the retained peak stores everything
          exercised = 0
          do c = 0, blocks - 1
             limit = whole % restart + whole % leaf + whole % window + whole % terms + c * whole % restart
             expected = c
             if (limit >= whole % retained) expected = -1
             if (expected >= 0) exercised = exercised + 1
             call configure(limited, solver_of(case_index))
             call limited % set_reverse_limit(int(limit))
             call initialize_case(other, case_index, limited)
             call finish(other)
             before = counted_solves(other % account(), linear_solves)
             call other % derivative(functional, order, reverse_pass, table, by_order=by_order, entries=entries, &
                  & designs=designs, storage=bounded)
             solves_bounded = counted_solves(other % account(), linear_solves) - before
             call check(bounded % checkpoints == max(expected, 0), &
                  & label // ': the limit admits the count of restart states it was set for')
             call check(bounded % peak <= limit .and. bounded % checkpoint % high <= c * bounded % restart, &
                  & label // ': the peak of the accounts is within the limit')
             if (expected >= 0) then
                call check(bounded % evaluations == recurrence_evaluations(blocks, c), &
                     & label // ': the forward evaluations are the count of the recurrence')
             else
                call check(bounded % evaluations == blocks, label // ': retention evaluates each block once')
             end if
             ! the coupled direct solver counts one linear solve per right
             ! side, so the extra solves are the recomputed towers exactly;
             ! a partitioned coupling solves each member of a tower
             if (solver_of(case_index) == 3) then
                call check(nint(solves_bounded - solves_retained) == (bounded % evaluations - blocks) * m, &
                     & label // ': the extra linear solves are the recomputed towers')
             else
                call check(nint(solves_bounded - solves_retained) >= (bounded % evaluations - blocks) * m .and. &
                     & ((nint(solves_bounded - solves_retained) > 0) .eqv. (bounded % evaluations > blocks)), &
                     & label // ': the extra linear solves are those of the recomputed towers')
             end if
             call check(all(shape(table) == shape(whole_table)) .and. all(table == whole_table), &
                  & label // ': the recomputed table equals the retained table bitwise')
             call check(all(by_order == whole_by_order) .and. all(entries == whole_entries), &
                  & label // ': the lower orders and the entries equal the retained ones bitwise')
          end do
          call check(case_index /= 5 .or. exercised >= 10, label // ': the long chain exercises ten checkpoint counts')
          ! a limit at the retained peak recomputes nothing; repeated on one execution and on a copy
          call configure(limited, solver_of(case_index))
          call limited % set_reverse_limit(int(whole % retained))
          call initialize_case(other, case_index, limited)
          call finish(other)
          call other % derivative(functional, order, reverse_pass, table, designs=designs, storage=bounded)
          call check(bounded % checkpoints == 0 .and. bounded % evaluations == blocks .and. all(table == whole_table), &
               & label // ': a limit at the retained peak stores every tower and evaluates each block once')
          call limited % set_reverse_limit(int(whole % minimum))
          call initialize_case(other, case_index, limited)
          call finish(other)
          copy = other
          call other % derivative(functional, order, reverse_pass, table, designs=designs, storage=bounded)
          call check(all(table == whole_table) .and. bounded % checkpoints == 0 .and. bounded % peak <= whole % minimum, &
               & label // ': the smallest admissible limit is met')
          call other % derivative(functional, order, reverse_pass, table, designs=designs, storage=bounded)
          call check(all(table == whole_table) .and. bounded % peak <= whole % minimum, &
               & label // ': a repeated bounded derivative on one execution equals the retained table')
          call copy % derivative(functional, order, reverse_pass, table, designs=designs, storage=bounded)
          call check(all(table == whole_table), label // ': a bounded derivative on a copy equals the retained table')
       end if
       call streamed_case(case_index, order, m, blocks, functional, designs, whole_table, whole_by_order, &
            & whole_entries, streamed, label)
    end do
  end subroutine bounded_case

  !-------------------------------------------------------------------!
  ! The streamed reverse execution of one case and order: the march
  ! stores the restart states its schedule names and retains the final
  ! leaf; the derivative recomputes states and towers. Tables, lower
  ! orders and entries equal the post-hoc retained pass bitwise at
  ! every checkpoint count; the peak is within the limit; the reverse
  ! phase's evaluations are the recurrence's count less the march's
  ! own sweep; every recomputed block is a Newton solve (exactly the
  ! march's count per block on the uniform chain); a second derivative
  ! finds no restart state and recomputes from the initial state; a
  ! copy taken before the derivative and a reinitialized execution
  ! reach the same table.
  !-------------------------------------------------------------------!
  subroutine streamed_case(case_index, order, m, blocks, functional, designs, whole_table, whole_by_order, &
       & whole_entries, streamed, label)
    integer, intent(in) :: case_index, order, m, blocks, designs
    type(expression), intent(in) :: functional(:)
    real(dp), intent(in) :: whole_table(:,:), whole_by_order(:,:,:), whole_entries(:,:,:)
    type(derivative_storage), intent(out) :: streamed
    character(len=*), intent(in) :: label
    type(march_context) :: limited
    type(chain_execution) :: other, copy
    type(derivative_storage) :: bounded
    real(dp), allocatable :: table(:,:), by_order(:,:,:), entries(:,:,:)
    real(dp) :: before_newton, before_linear, march_newton, march_linear, extra_newton, extra_linear
    real(dp) :: newton_retained, linear_retained
    integer :: c, expected
    integer(int64) :: limit
    logical :: uniform
    uniform = case_index == 5
    ! retention: the march keeps every state and tower, the derivative recomputes nothing
    call configure(limited, solver_of(case_index))
    call initialize_case(other, case_index, limited, reverse_order=order)
    before_newton = counted_solves(other % account(), newton_solves)
    before_linear = counted_solves(other % account(), linear_solves)
    call finish(other)
    march_newton = counted_solves(other % account(), newton_solves) - before_newton
    march_linear = counted_solves(other % account(), linear_solves) - before_linear
    ! the solves of the reverse phase without recomputation: the
    ! costate solves, which a partitioned coupling counts as Newton
    ! solves of one iteration; the extra solves of a bounded run are
    ! the deltas against these
    before_newton = counted_solves(other % account(), newton_solves)
    before_linear = counted_solves(other % account(), linear_solves)
    call other % derivative(functional, order, reverse_pass, table, by_order=by_order, entries=entries, &
         & designs=designs, storage=streamed)
    newton_retained = counted_solves(other % account(), newton_solves) - before_newton
    linear_retained = counted_solves(other % account(), linear_solves) - before_linear
    call check(streamed % checkpoints == 0 .and. streamed % evaluations == 0 .and. streamed % peak <= streamed % retained, &
         & label // ': a streamed reverse execution without a limit retains its march')
    call check(streamed % restart > 0, label // ': the restart state of a streamed execution stores states at every order')
    call check(all(shape(table) == shape(whole_table)) .and. all(table == whole_table) .and. &
         & all(by_order == whole_by_order) .and. all(entries == whole_entries), &
         & label // ': the streamed reverse tables equal the post-hoc retained tables bitwise')
    do c = 0, blocks - 1
       limit = streamed % restart + streamed % leaf + streamed % window + streamed % terms + c * streamed % restart
       expected = c
       if (limit >= streamed % retained) expected = -1
       call configure(limited, solver_of(case_index))
       call limited % set_reverse_limit(int(limit))
       call initialize_case(other, case_index, limited, reverse_order=order)
       before_newton = counted_solves(other % account(), newton_solves)
       before_linear = counted_solves(other % account(), linear_solves)
       call finish(other)
       call check(nint(counted_solves(other % account(), newton_solves) - before_newton) == nint(march_newton) .and. &
            & nint(counted_solves(other % account(), linear_solves) - before_linear) == nint(march_linear), &
            & label // ': the march of a bounded streamed execution solves what the retained march solves')
       if (c == 0) copy = other
       before_newton = counted_solves(other % account(), newton_solves)
       before_linear = counted_solves(other % account(), linear_solves)
       call other % derivative(functional, order, reverse_pass, table, by_order=by_order, entries=entries, &
            & designs=designs, storage=bounded)
       extra_newton = counted_solves(other % account(), newton_solves) - before_newton - newton_retained
       extra_linear = counted_solves(other % account(), linear_solves) - before_linear - linear_retained
       call check(bounded % checkpoints == max(expected, 0), &
            & label // ': the streamed limit admits the count of restart states it was set for')
       call check(bounded % peak <= limit .and. bounded % checkpoint % high <= c * bounded % restart, &
            & label // ': the streamed peak of the accounts is within the limit')
       if (expected >= 0) then
          call check(bounded % evaluations == recurrence_evaluations(blocks, c) - blocks, &
               & label // ': the reverse phase evaluates the recurrence''s count less the march')
       else
          call check(bounded % evaluations == 0, label // ': streamed retention recomputes nothing')
       end if
       call check((nint(extra_newton) > 0) .eqv. (bounded % evaluations > 0), &
            & label // ': a recomputed block is a Newton solve')
       if (uniform) then
          call check(mod(nint(march_newton), blocks) == 0 .and. &
               & nint(extra_newton) == bounded % evaluations * (nint(march_newton) / blocks), &
               & label // ': the extra Newton solves are the recomputed blocks')
       end if
       call check(nint(extra_linear) >= bounded % evaluations * m .and. &
            & nint(extra_linear) <= bounded % evaluations * nint(march_linear), &
            & label // ': the extra linear solves are the recomputed towers and Newton iterations')
       call check(all(shape(table) == shape(whole_table)) .and. all(table == whole_table), &
            & label // ': the streamed recomputed table equals the retained table bitwise')
       call check(all(by_order == whole_by_order) .and. all(entries == whole_entries), &
            & label // ': the streamed lower orders and entries equal the retained ones bitwise')
       if (c == 0) then
          ! a second derivative finds no restart state stored: recomputed from the initial state
          call other % derivative(functional, order, reverse_pass, table, designs=designs, storage=bounded)
          call check(all(table == whole_table) .and. bounded % peak <= limit, &
               & label // ': a repeated streamed derivative recomputes from the initial state')
          call copy % derivative(functional, order, reverse_pass, table, designs=designs, storage=bounded)
          call check(all(table == whole_table) .and. bounded % peak <= limit, &
               & label // ': a copy of a streamed execution reaches the same table')
          call initialize_case(other, case_index, limited, reverse_order=order)
          call finish(other)
          call other % derivative(functional, order, reverse_pass, table, designs=designs, storage=bounded)
          call check(all(table == whole_table) .and. bounded % peak <= limit, &
               & label // ': a reinitialized streamed execution reaches the same table')
       end if
    end do
  end subroutine streamed_case

  ! The schedule's recurrence, written independently: t(0) = 0, t(1) = 1,
  ! t(n, 0) = n (n + 1) / 2, t(n, c) = min over m of m + t(n - m, c - 1)
  ! + t(m - 1, c), a restart state stored after block m that block's forward
  ! data included. By hand: t(2, 0) = 3, t(4, 0) = 10, t(2, 1) = 2,
  ! t(4, 1) = 6 (advance 1..3, store, evaluate 4, reverse 4 and 3,
  ! evaluate 1, store, evaluate 2, reverse 2 and 1).
  recursive integer function recurrence_evaluations(n, c) result(t)
    integer, intent(in) :: n, c
    integer :: m, at
    if (n <= 1) then
       t = n
    else if (c == 0) then
       t = n * (n + 1) / 2
    else
       t = huge(1)
       do m = 1, n - 1
          at = m + recurrence_evaluations(n - m, c - 1) + recurrence_evaluations(m - 1, c)
          if (at < t) t = at
       end do
    end if
  end function recurrence_evaluations

  ! the solves an execution counted in its own account
  real(dp) function counted_solves(account, event)
    type(tally), intent(in) :: account
    integer, intent(in) :: event
    integer :: level, order
    counted_solves = 0.0_dp
    do level = 1, account % num_levels()
       do order = 0, 3
          counted_solves = counted_solves + account % amount(level, order, event)
       end do
    end do
  end function counted_solves

  ! the solver settings of a case: the long chain counts one linear
  ! solve per right side, so the extra solves are the recomputed towers
  integer function solver_of(case_index)
    integer, intent(in) :: case_index
    solver_of = 1
    if (case_index == 2) solver_of = 2
    if (case_index == 5) solver_of = 3
  end function solver_of

  integer function degree_of(case_index)
    integer, intent(in) :: case_index
    degree_of = merge(3, 2, case_index == 2)
  end function degree_of

  ! the refusal: a limit one entry below the smallest admissible
  subroutine refused_limit()
    type(march_context) :: context, limited
    type(chain_execution) :: execution
    type(derivative_storage) :: whole
    real(dp), allocatable :: table(:,:)
    call configure(context, 1)
    call initialize_case(execution, 1, context)
    call finish(execution)
    call execution % derivative([van_der_pol_energy(2)], 2, reverse_pass, table, storage=whole)
    call configure(limited, 1)
    call limited % set_reverse_limit(int(whole % minimum) - 1)
    call initialize_case(execution, 1, limited)
    call finish(execution)
    call execution % derivative([van_der_pol_energy(2)], 2, reverse_pass, table, storage=whole)
  end subroutine refused_limit

  !-------------------------------------------------------------------!
  ! Copy semantics. A copy of an execution is an independent execution
  ! with equal state: it reaches the same results by the same
  ! arithmetic, so states, derivative tables and Taylor coefficients of
  ! a copy and its source are equal exactly, and both equal the
  ! isolated reference within the declared tolerance.
  !-------------------------------------------------------------------!
  subroutine check_copies()
    type(chain_execution) :: source, copy, elements(2)
    type(chain_execution), allocatable :: first_destroyed, survivor
    type(execution_holder) :: holder, holder_copy
    type(execution_result) :: from_source, from_copy, from_elements(2), from_holder
    type(execution_result) :: streamed_source, streamed_copy

    ! before any advance
    call initialize_case(source, 1, context_a)
    copy = source
    call check(.not. copy % complete(), 'a copy of an unevaluated execution is unevaluated')
    call finish(source)
    call check(.not. copy % complete(), 'advancing the source does not advance the copy')
    call finish(copy)
    call differentiate(source, 2, from_source)
    call differentiate(copy, 2, from_copy)
    call extract(source, from_source)
    call extract(copy, from_copy)
    call compare_exact(from_copy, from_source)
    call compare(from_copy, reference_a)
    print *, 'PASS: a copy taken before the march reaches the same states and derivatives'

    ! during the march, then alternating progress and reinitialization
    call initialize_case(source, 1, context_a)
    call source % advance()
    copy = source
    call source % advance()
    call check(.not. copy % complete(), 'the copy keeps its own progress')
    call initialize_case(source, 2, context_b)
    call copy % advance()
    call source % advance()
    call finish(copy)
    call finish(source)
    call differentiate(copy, 2, from_copy)
    call extract(copy, from_copy)
    call compare(from_copy, reference_a)
    call differentiate(source, 3, from_source)
    call extract(source, from_source)
    call compare(from_source, reference_b)
    print *, 'PASS: a copy taken during the march continues independently of the reinitialized source'

    ! after the march, with derivatives already computed on the source
    call initialize_case(source, 1, context_a)
    call finish(source)
    call differentiate_order(source, 2, 2, forward_pass, from_source)
    copy = source
    call differentiate(copy, 2, from_copy)
    call differentiate(source, 2, from_source)
    call check(all(from_copy % forward_derivative == from_source % forward_derivative) .and. &
         & all(from_copy % reverse_derivative == from_source % reverse_derivative), &
         & 'forward and reverse derivatives of a copy taken after the march equal the source exactly')
    call extract(copy, from_copy)
    call extract(source, from_source)
    call compare_exact(from_copy, from_source)
    print *, 'PASS: a copy taken after the march has equal forward and reverse derivatives'

    ! the streamed Taylor march, copied during the march
    call initialize_case(source, 1, context_a, streamed=.true.)
    call source % advance()
    copy = source
    call finish(source)
    call finish(copy)
    call extract(source, streamed_source, streamed=.true.)
    call extract(copy, streamed_copy, streamed=.true.)
    call check(all(streamed_copy % functional == streamed_source % functional) .and. &
         & all(streamed_copy % tower_storage == streamed_source % tower_storage) .and. &
         & all(streamed_copy % state_storage == streamed_source % state_storage), &
         & 'a streamed copy owns its Taylor coefficients and storage accounting')
    call compare_stream(streamed_copy, stream_a, reference_a)
    print *, 'PASS: a copy taken during a streamed Taylor march reaches the same coefficients'

    ! source destroyed first; destination destroyed first
    allocate(first_destroyed, survivor)
    call initialize_case(first_destroyed, 1, context_a)
    call first_destroyed % advance()
    survivor = first_destroyed
    deallocate(first_destroyed)
    call finish(survivor)
    call differentiate(survivor, 2, from_copy)
    call extract(survivor, from_copy)
    call compare(from_copy, reference_a)
    allocate(first_destroyed)
    call initialize_case(survivor, 1, context_a)
    call survivor % advance()
    first_destroyed = survivor
    deallocate(first_destroyed)
    call finish(survivor)
    call differentiate(survivor, 2, from_source)
    call extract(survivor, from_source)
    call compare_exact(from_source, from_copy)
    deallocate(survivor)
    print *, 'PASS: the survivor of either destruction order completes the march'

    ! a containing object, array elements and a function result
    call initialize_case(source, 1, context_a)
    call source % advance()
    holder % execution = source
    holder % label = 1
    holder_copy = holder
    elements(1) = source
    elements(2) = elements(1)
    copy = execution_copy(source)
    call finish(source)
    call finish(holder_copy % execution)
    call finish(elements(1))
    call finish(elements(2))
    call finish(copy)
    call check(.not. holder % execution % complete(), 'the contained source keeps its own progress')
    call finish(holder % execution)
    call differentiate(source, 2, from_source)
    call differentiate(holder_copy % execution, 2, from_holder)
    call differentiate(elements(1), 2, from_elements(1))
    call differentiate(elements(2), 2, from_elements(2))
    call differentiate(copy, 2, from_copy)
    call extract(source, from_source)
    call extract(holder_copy % execution, from_holder)
    call extract(elements(1), from_elements(1))
    call extract(elements(2), from_elements(2))
    call extract(copy, from_copy)
    call compare_exact(from_holder, from_source)
    call compare_exact(from_elements(1), from_source)
    call compare_exact(from_elements(2), from_source)
    call compare_exact(from_copy, from_source)
    call compare(from_source, reference_a)
    print *, 'PASS: copies through a container, array elements and a function result are complete executions'

    ! a source= twin is not an owner: it computes the same results while
    ! both live, and taking its results releases the one binding both
    ! carry; run.sh checks that the source is refused after that
    call initialize_case(source, 1, context_a)
    call source % advance()
    allocate(twin, source=source)
    call finish(twin)
    call finish(source)
    call differentiate(twin, 2, from_copy)
    call differentiate(source, 2, from_source)
    call check(all(from_copy % forward_derivative == from_source % forward_derivative) .and. &
         & all(from_copy % reverse_derivative == from_source % reverse_derivative), &
         & 'a source= twin has exactly the derivatives of its source while both live')
    call extract(twin, from_copy)
    call compare(from_copy, reference_a)
    deallocate(twin)
    print *, 'PASS: a source= twin computes the same results while both live; run.sh checks its refusal after'
  end subroutine check_copies

  type(chain_execution) function execution_copy(source)
    type(chain_execution), intent(in) :: source
    execution_copy = source
  end function execution_copy

  subroutine compare_exact(actual, expected)
    type(execution_result), intent(in) :: actual, expected
    call check(size(actual % states) == size(expected % states), 'copied trajectory extents agree')
    call check(all(actual % states == expected % states), 'a copy reaches exactly the states of its source')
    call check(all(actual % steps == expected % steps) .and. all(actual % times == expected % times), &
         & 'a copy has exactly the grid of its source')
    call check(all(shape(actual % forward_derivative) == shape(expected % forward_derivative)) .and. &
         & all(actual % forward_derivative == expected % forward_derivative), &
         & 'a copy has exactly the forward derivatives of its source')
    call check(all(shape(actual % reverse_derivative) == shape(expected % reverse_derivative)) .and. &
         & all(actual % reverse_derivative == expected % reverse_derivative), &
         & 'a copy has exactly the reverse derivatives of its source')
  end subroutine compare_exact

  subroutine configure(context, case_index)
    type(march_context), intent(out) :: context
    integer, intent(in) :: case_index
    call context % account % open(3, hierarchy_levels)
    call context % set_stopping(1.0e-11_dp, relative, by_count, 100)
    if (case_index == 1) then
       call context % set_linear_solver('direct')
       call context % set_storage('dense')
       call context % set_time_coupling('sequential')
    else if (case_index == 3) then
       ! one linear solve per right side: the direct solver over the coupled block
       call context % set_linear_solver('direct')
       call context % set_storage('dense')
       call context % set_time_coupling('coupled')
    else
       call context % set_linear_solver('iterative')
       call context % set_storage('sparse')
       call context % set_preconditioner('gauss_seidel')
       call context % set_linear_limits(40, 2, 200)
       call context % set_time_coupling('coupled')
    end if
  end subroutine configure

  subroutine initialize_case(execution, case_index, context, streamed, reverse_order)
    type(chain_execution), intent(inout) :: execution
    integer, intent(in) :: case_index
    type(march_context), intent(inout) :: context
    logical, intent(in), optional :: streamed
    ! a streamed reverse execution of this order
    integer, intent(in), optional :: reverse_order
    type(family_container), allocatable :: schemes(:)
    type(expression) :: physics, functional(1)
    integer, allocatable :: added(:)
    real(dp), allocatable :: initial(:), lower(:)
    real(dp) :: duration, design
    integer :: degree, order, b
    if (case_index >= 3) then
       call initialize_mixed(execution, case_index, context, reverse_order)
       return
    end if
    if (case_index == 1) then
       degree = 2
       order = 2
       duration = 1.0_dp
       design = 0.4_dp
       lower = [1.0_dp, 0.0_dp]
       added = [5, 4, 3]
    else
       degree = 3
       order = 1
       duration = 0.75_dp
       design = 0.9_dp
       lower = [0.5_dp, -0.25_dp, 0.125_dp]
       added = [4, 3]
    end if
    allocate(schemes(size(added)))
    do b = 1, size(schemes)
       allocate(schemes(b) % scheme, source=bdf_family(order))
    end do
    physics = van_der_pol(degree)
    functional = [van_der_pol_energy(degree)]
    initial = consistent_state(physics, degree + 1, lower, design, context=context)
    call started(execution, schemes, added, physics, degree + 1, uniform_grid(duration), design, initial, &
         & functional, context, startup=4, streamed=streamed, reverse_order=reverse_order)
  end subroutine initialize_case

  ! one initialization: primal, streamed forward Taylor to order two,
  ! or streamed reverse to the order given
  subroutine started(execution, schemes, added, physics, degrees, steps, design, initial, functional, context, &
       & startup, grid_design, streamed, reverse_order, designs)
    type(chain_execution), intent(inout) :: execution
    type(family_container), intent(in) :: schemes(:)
    integer, intent(in) :: added(:), degrees
    type(expression), intent(in) :: physics, functional(:)
    class(grid), intent(in) :: steps
    real(dp), intent(in) :: design, initial(:)
    type(march_context), intent(in) :: context
    integer, intent(in), optional :: startup, reverse_order, designs
    real(dp), intent(in), optional :: grid_design(:)
    logical, intent(in), optional :: streamed
    logical :: stream
    stream = .false.
    if (present(streamed)) stream = streamed
    if (present(reverse_order)) then
       call execution % initialize(schemes, added, physics, degrees, steps, design, initial, grid_design=grid_design, &
            & startup=startup, functionals=functional, derivative_order=reverse_order, context=context, &
            & pass_kind=reverse_pass, designs=designs)
    else if (stream) then
       call execution % initialize(schemes, added, physics, degrees, steps, design, initial, grid_design=grid_design, &
            & startup=startup, functionals=functional, derivative_order=2, context=context)
    else
       call execution % initialize(schemes, added, physics, degrees, steps, design, initial, grid_design=grid_design, &
            & startup=startup, context=context)
    end if
  end subroutine started

  ! case 3: a mixed chain, bdf 3 / adams 3 / crouzeix / newmark over a
  ! startup of four fine steps; case 4: bdf 2 blocks on a designed grid
  ! with a second design (the step weights); case 5: sixteen
  ! one-instant crouzeix blocks; case 6: the mixed chain three times
  ! over, thirteen blocks with the startup, long enough for its
  ! restart states (two blocks where bdf 3 reaches back) to be stored
  ! and recomputed
  subroutine initialize_mixed(execution, case_index, context, reverse_order)
    type(chain_execution), intent(inout) :: execution
    integer, intent(in) :: case_index
    type(march_context), intent(inout) :: context
    integer, intent(in), optional :: reverse_order
    type(family_container), allocatable :: schemes(:)
    type(expression) :: physics, functional(1)
    integer, allocatable :: added(:)
    real(dp), allocatable :: initial(:), weights(:)
    real(dp) :: duration, design
    integer :: b, k
    duration = 1.0_dp
    design = 0.4_dp
    physics = van_der_pol(2)
    functional = [van_der_pol_energy(2)]
    initial = consistent_state(physics, 3, [1.0_dp, 0.0_dp], design, context=context)
    if (case_index == 5) then
       ! sixteen one-instant crouzeix blocks: the reach is one block
       added = [2, (1, b = 2, 16)]
       allocate(schemes(16))
       do b = 1, 16
          allocate(schemes(b) % scheme, source=crouzeix_three_stage())
       end do
       call started(execution, schemes, added, physics, 3, uniform_grid(duration), design, initial, functional, &
            & context, reverse_order=reverse_order)
    else if (case_index == 3 .or. case_index == 6) then
       if (case_index == 3) then
          added = [5, 4, 3, 3]
       else
          added = [5, 3, 1, 1, 4, 3, 1, 1, 4, 3, 1, 1]
       end if
       allocate(schemes(size(added)))
       do b = 1, size(added), 4
          allocate(schemes(b) % scheme, source=bdf_family(3))
          allocate(schemes(b + 1) % scheme, source=adams_family(3))
          allocate(schemes(b + 2) % scheme, source=crouzeix_three_stage())
          allocate(schemes(b + 3) % scheme, source=newmark_family(0.25_dp, 0.5_dp))
       end do
       call started(execution, schemes, added, physics, 3, uniform_grid(duration), design, initial, functional, &
            & context, startup=4, reverse_order=reverse_order)
    else
       added = [5, 4, 3]
       allocate(schemes(3))
       do b = 1, 3
          allocate(schemes(b) % scheme, source=bdf_family(2))
       end do
       weights = [(1.0_dp + 0.05_dp * real(mod(k, 3), dp), k = 1, sum(added) - 1)]
       call started(execution, schemes, added, physics, 3, designed_grid(duration), design, initial, functional, &
            & context, startup=4, grid_design=weights, reverse_order=reverse_order, designs=2)
    end if
  end subroutine initialize_mixed

  subroutine finish(execution)
    type(chain_execution), intent(inout) :: execution
    integer :: advances
    advances = 0
    do while (.not. execution % complete())
       call execution % advance()
       advances = advances + 1
       call check(advances <= 32, 'a finite chain reaches completion')
    end do
  end subroutine finish

  subroutine differentiate(execution, degree, result)
    type(chain_execution), intent(inout) :: execution
    integer, intent(in) :: degree
    type(execution_result), intent(inout) :: result
    integer :: order
    do order = 1, 2
       call differentiate_order(execution, degree, order, forward_pass, result)
       call differentiate_order(execution, degree, order, reverse_pass, result)
    end do
  end subroutine differentiate

  subroutine differentiate_order(execution, degree, order, orientation, result)
    type(chain_execution), intent(inout) :: execution
    integer, intent(in) :: degree, order, orientation
    type(execution_result), intent(inout) :: result
    real(dp), allocatable :: table(:,:)
    call execution % derivative([van_der_pol_energy(degree)], order, orientation, table)
    call check(all(shape(table) == [1, 1]), 'one functional has one derivative in the single design')
    if (orientation == forward_pass) then
       if (.not. allocated(result % forward_derivative)) allocate(result % forward_derivative(2, 1))
       result % forward_derivative(order, 1) = table(1, 1)
    else
       if (.not. allocated(result % reverse_derivative)) allocate(result % reverse_derivative(2, 1))
       result % reverse_derivative(order, 1) = table(1, 1)
    end if
  end subroutine differentiate_order

  subroutine isolated(case_index, context, result)
    integer, intent(in) :: case_index
    type(march_context), intent(inout) :: context
    type(execution_result), intent(inout) :: result
    type(chain_execution) :: execution
    call initialize_case(execution, case_index, context)
    call finish(execution)
    call differentiate(execution, case_index + 1, result)
    call extract(execution, result)
  end subroutine isolated

  subroutine extract(execution, result, streamed)
    type(chain_execution), intent(inout) :: execution
    type(execution_result), intent(inout) :: result
    logical, intent(in), optional :: streamed
    type(chain_block), allocatable :: chain(:)
    type(expansion), allocatable, target :: tower
    type(imbalance) :: outcome
    real(dp) :: achieved
    integer :: b
    logical :: stream
    stream = .false.
    if (present(streamed)) stream = streamed
    if (stream) then
       call execution % take_results(chain, tower, result % steps, result % times, achieved, &
            & final_imbalance=outcome, f=result % functional, tower_storage=result % tower_storage, &
            & state_storage=result % state_storage)
       call check(all(result % tower_storage > 0) .and. all(result % state_storage > 0), &
            & 'streamed storage accounting describes nonempty states and Taylor coefficients')
       do b = 1, size(chain)
          call check(.not. allocated(chain(b) % state), &
               & 'every streamed primal state is retired after its final functional or history reader')
       end do
    else
       call execution % take_results(chain, tower, result % steps, result % times, achieved, final_imbalance=outcome)
       result % states = [real(dp) ::]
       do b = 1, size(chain)
          result % states = [result % states, chain(b) % state]
       end do
    end if
    call check(outcome % converged, 'the execution reports a converged primal result')
  end subroutine extract

  subroutine compare(actual, expected)
    type(execution_result), intent(in) :: actual, expected
    call check(size(actual % states) == size(expected % states), 'trajectory extents agree')
    call check(maxval(abs(actual % states - expected % states)) < 2.0e-10_dp, 'all trajectory states agree')
    call check(all(actual % steps == expected % steps), 'time increments agree')
    call check(all(actual % times == expected % times), 'physical instants agree')
    call compare_matrix(actual % forward_derivative, expected % forward_derivative, 2.0e-10_dp)
    call compare_matrix(actual % reverse_derivative, expected % reverse_derivative, 2.0e-10_dp)
    call compare_matrix(actual % forward_derivative, actual % reverse_derivative, 2.0e-7_dp)
  end subroutine compare

  subroutine compare_stream(actual, expected, retained)
    type(execution_result), intent(in) :: actual, expected, retained
    integer :: first
    call compare_matrix(actual % functional, expected % functional, 2.0e-10_dp)
    call check(all(actual % tower_storage == expected % tower_storage), 'Taylor storage is independent')
    call check(all(actual % state_storage == expected % state_storage), 'primal storage is independent')
    call check(actual % tower_storage(1) <= actual % tower_storage(2), 'live Taylor storage cannot exceed total storage')
    call check(actual % state_storage(1) <= actual % state_storage(2), 'live primal storage cannot exceed total storage')
    first = lbound(actual % functional, 1)
    call compare_matrix(actual % functional(first + 1:, :), retained % forward_derivative, 2.0e-7_dp)
  end subroutine compare_stream

  subroutine compare_matrix(actual, expected, tolerance)
    real(dp), intent(in) :: actual(:,:), expected(:,:), tolerance
    call check(all(shape(actual) == shape(expected)), 'derivative table extents agree')
    call check(maxval(abs(actual - expected)) <= tolerance * max(1.0_dp, maxval(abs(expected))), &
         & 'the same derivative table is obtained independently')
  end subroutine compare_matrix

  subroutine check(satisfied, description)
    logical, intent(in) :: satisfied
    character(len=*), intent(in) :: description
    if (.not. satisfied) error stop 'execution_independence: ' // description
  end subroutine check

end program execution_independence
