!=====================================================================!
! Independent concurrent executions.
!
! A set of heterogeneous executions E_1..E_n (state degrees 2 and 3;
! BDF, Adams, DIRK and Newmark families; direct, GMRES with Gauss-Seidel,
! GMRES with one multigrid cycle and numerical-elimination solvers;
! designs and random grids; forward
! and reverse derivatives of orders 1 and 2; one streamed Taylor
! execution; one execution configured not to converge) is run serially
! into the reference S, then over `!$omp parallel do schedule(dynamic,1)`
! at the requested thread count, five times. Every per-execution
! quantity - states, grid, derivative tables, Taylor coefficients,
! storage pairs, convergence, the failed outcome, and the counted
! accounting events at every level and order - equals S at tolerance
! zero: each execution's arithmetic runs in the program order of the
! serial build. A shared template execution is copied by every
! iteration and is still advanceable after the join. Nothing is
! written inside the parallel region but each execution's own file;
! every line of standard output is written after the join in index
! order. Without -fopenmp the directives are comments and the same
! program is the serial fallback.
!
! usage: run <output directory> [threads] [measure]
!=====================================================================!
program concurrent_executions

  !$ use omp_lib, only : omp_set_num_threads, omp_get_max_threads
  use util_precision, only : dp
  use operation_minimization, only : relative, by_count, solve_result
  use operation_family, only : bdf_family, adams_family, newmark_family, &
       & crouzeix_three_stage, implicit_midpoint
  use operation_grid, only : uniform_grid, random_grid
  use operation_expression, only : expression
  use gti_physics, only : van_der_pol, van_der_pol_energy
  use gti_expansion, only : expansion, family_container
  use gti_march, only : march_context, imbalance, consistent_state
  use gti_sweeps, only : forward_pass, reverse_pass
  use gti_chain, only : chain_execution, chain_block
  use gti_configuration, only : hierarchy_levels
  use util_tally, only : tally, tally_num_events, elapsed_time

  implicit none

  integer, parameter :: num_executions = 10, num_repetitions = 5, highest_order = 2
  integer, parameter :: streamed_case = 7, failing_case = 8

  type :: execution_record
     logical :: converged = .false.
     logical :: derivative_converged = .true.
     integer :: failure_reason = 0
     real(dp), allocatable :: states(:), steps(:), times(:)
     real(dp) :: forward(highest_order) = 0.0_dp, reverse(highest_order) = 0.0_dp
     real(dp), allocatable :: functional(:,:)
     integer :: tower_storage(2) = 0, state_storage(2) = 0
     real(dp), allocatable :: counts(:,:,:)
  end type execution_record

  type(execution_record) :: reference(num_executions), concurrent(num_executions)
  type(chain_execution) :: template(num_executions)
  type(march_context) :: context
  character(len=256) :: out_dir, argument
  character(len=64) :: label
  integer :: i, repetition, threads_requested, threads_in_use
  logical :: measure

  call get_command_argument(1, out_dir)
  if (len_trim(out_dir) == 0) error stop 'concurrent_executions: an output directory is given'
  threads_requested = 1
  call get_command_argument(2, argument)
  if (len_trim(argument) > 0) read(argument, *) threads_requested
  call get_command_argument(3, argument)
  measure = trim(argument) == 'measure'

  threads_in_use = 1
  !$ call omp_set_num_threads(threads_requested)
  !$ threads_in_use = omp_get_max_threads()
  write(*,'(a,i0,a,i0)') ' threads requested ', threads_requested, ' in use ', threads_in_use

  call execute_command_line('mkdir -p ' // trim(out_dir) // '/serial ' // trim(out_dir) // '/measure')
  do repetition = 1, num_repetitions
     write(label, '(a,i0,a,i0)') '/threads_', threads_in_use, '_repetition_', repetition
     call execute_command_line('mkdir -p ' // trim(out_dir) // trim(label))
  end do

  if (measure) then
     call measured()
     stop
  end if

  ! the serial reference S
  do i = 1, num_executions
     call run_case(i, reference(i), trim(out_dir) // '/serial')
  end do
  do i = 1, num_executions
     call summarised(i, reference(i))
  end do
  call check(.not. reference(failing_case) % converged, 'the failing execution reports a non-converged primal')
  call check(.not. reference(failing_case) % derivative_converged, &
       & 'a derivative of the non-converged primal is a non-converged result')
  do i = 1, num_executions
     if (i /= failing_case) call check(reference(i) % converged, 'every other execution converges')
  end do
  print *, 'PASS: the serial reference set converges except the execution configured not to'

  ! the same set concurrently, five times
  do repetition = 1, num_repetitions
     write(label, '(a,i0,a,i0)') '/threads_', threads_in_use, '_repetition_', repetition
     !$omp parallel do schedule(dynamic,1) default(shared) private(i)
     do i = 1, num_executions
        call run_case(i, concurrent(i), trim(out_dir) // trim(label))
     end do
     !$omp end parallel do
     do i = 1, num_executions
        call compare_exact(concurrent(i), reference(i), i)
     end do
  end do
  write(*,'(a,i0,a)') ' PASS: ', num_repetitions, &
       & ' concurrent runs of the set equal the serial reference at tolerance zero'
  print *, 'PASS: the failing execution leaves the other executions intact and the process alive'
  print *, 'PASS: per-execution accounting counts equal their serial counts'

  ! a shared template: every iteration copies one execution
  do i = 1, num_executions
     call initialize_case(template(i), i, context)
  end do
  !$omp parallel do schedule(dynamic,1) default(shared) private(i)
  do i = 1, num_executions
     call run_copy(template(i), i, concurrent(i))
  end do
  !$omp end parallel do
  do i = 1, num_executions
     call compare_exact(concurrent(i), reference(i), i)
     call check(.not. template(i) % complete(), 'the template is unevaluated after its copies complete')
  end do
  ! the template's cells are still owned: it marches after the join
  call run_copy(template(1), 1, concurrent(1))
  call finish(template(1))
  call check(template(1) % complete(), 'the shared template marches after the join')
  print *, 'PASS: copies of a shared template execution reach the reference; the template survives the join'

  write(*,'(a,es24.16)') ' sum of forward first derivatives over the set, formed after the join in index order: ', &
       & total_after_join()

contains

  !-------------------------------------------------------------------!
  ! One execution from initialization to its results and account.
  !-------------------------------------------------------------------!
  subroutine run_case(case_index, record, directory)
    integer, intent(in) :: case_index
    type(execution_record), intent(out) :: record
    character(len=*), intent(in) :: directory
    type(chain_execution) :: execution
    type(march_context) :: context
    call initialize_case(execution, case_index, context)
    call run_copy(execution, case_index, record)
    call written(record, case_index, directory)
  end subroutine run_case

  !-------------------------------------------------------------------!
  ! An initialized execution (or a copy of one) marched, differentiated
  ! and read; the copy's ownership ends with this call.
  !-------------------------------------------------------------------!
  subroutine run_copy(source, case_index, record)
    type(chain_execution), intent(in) :: source
    integer, intent(in) :: case_index
    type(execution_record), intent(out) :: record
    type(chain_execution) :: execution
    execution = source
    call finish(execution)
    if (case_index /= streamed_case) call differentiate(execution, case_index, record)
    call extract(execution, case_index, record)
  end subroutine run_copy

  subroutine initialize_case(execution, case_index, context)
    type(chain_execution), intent(inout) :: execution
    integer, intent(in) :: case_index
    type(march_context), intent(out) :: context
    type(family_container), allocatable :: schemes(:)
    type(expression) :: physics, functional(1)
    integer, allocatable :: added(:)
    real(dp), allocatable :: initial(:), lower(:)
    real(dp) :: duration, design
    integer :: degree, startup, seed, b
    call context % set_stopping(1.0e-11_dp, relative, by_count, 100)
    call context % account % open(highest_order, hierarchy_levels)
    startup = 0
    seed = 0
    select case (case_index)
    case (1, failing_case)
       degree = 2; duration = 1.0_dp; design = 0.4_dp; lower = [1.0_dp, 0.0_dp]; added = [5, 4, 3]; startup = 4
       allocate(schemes(size(added)))
       do b = 1, size(added); allocate(schemes(b) % scheme, source=bdf_family(2)); end do
       call context % set_linear_solver('direct'); call context % set_storage('dense')
       call context % set_time_coupling('sequential')
       if (case_index == failing_case) call context % set_stopping(1.0e-14_dp, relative, by_count, 1)
    case (2)
       degree = 3; duration = 0.75_dp; design = 0.9_dp; lower = [0.5_dp, -0.25_dp, 0.125_dp]; added = [4, 3]
       allocate(schemes(size(added)))
       do b = 1, size(added); allocate(schemes(b) % scheme, source=bdf_family(1)); end do
       call context % set_linear_solver('iterative'); call context % set_storage('sparse')
       call context % set_preconditioner('gauss_seidel'); call context % set_linear_limits(40, 2, 200)
       call context % set_time_coupling('coupled')
    case (3)
       degree = 2; duration = 1.0_dp; design = 0.1_dp; lower = [0.8_dp, 0.2_dp]; added = [6, 4]; startup = 3; seed = 11
       allocate(schemes(size(added)))
       do b = 1, size(added); allocate(schemes(b) % scheme, source=adams_family(2)); end do
       call context % set_linear_solver('direct'); call context % set_rows('states')
       call context % set_elimination('numerical')
    case (4)
       degree = 2; duration = 0.8_dp; design = 0.3_dp; lower = [1.0_dp, 0.5_dp]; added = [4, 4]
       allocate(schemes(size(added)))
       do b = 1, size(added); allocate(schemes(b) % scheme, source=crouzeix_three_stage()); end do
       call context % set_linear_solver('iterative'); call context % set_storage('sparse')
       call context % set_preconditioner('gauss_seidel'); call context % set_linear_limits(30, 2, 200)
    case (5)
       degree = 2; duration = 1.0_dp; design = 0.2_dp; lower = [1.0_dp, 0.0_dp]; added = [5, 5]
       allocate(schemes(size(added)))
       do b = 1, size(added); allocate(schemes(b) % scheme, source=newmark_family(0.25_dp, 0.5_dp)); end do
       call context % set_linear_solver('direct')
    case (6)
       degree = 3; duration = 0.5_dp; design = 0.6_dp; lower = [0.5_dp, 0.0_dp, 0.25_dp]; added = [4, 3, 3]
       startup = 2; seed = 7
       allocate(schemes(size(added)))
       do b = 1, size(added); allocate(schemes(b) % scheme, source=bdf_family(2)); end do
       call context % set_linear_solver('direct'); call context % set_elimination('numerical')
    case (streamed_case)
       degree = 2; duration = 1.0_dp; design = 0.4_dp; lower = [1.0_dp, 0.0_dp]; added = [5, 4, 3]; startup = 4
       allocate(schemes(size(added)))
       do b = 1, size(added); allocate(schemes(b) % scheme, source=bdf_family(2)); end do
       call context % set_linear_solver('direct'); call context % set_time_coupling('sequential')
    case (9)
       degree = 2; duration = 2.0_dp; design = 0.0_dp; lower = [1.0_dp, 0.0_dp]; added = [8]
       allocate(schemes(1)); allocate(schemes(1) % scheme, source=implicit_midpoint())
       call context % set_linear_solver('iterative'); call context % set_storage('sparse')
       call context % set_linear_limits(20, 1, 100); call context % set_time_coupling('sequential')
    case (10)
       degree = 2; duration = 0.6_dp; design = 0.5_dp; lower = [0.9_dp, 0.1_dp]; added = [5, 4]; startup = 2
       allocate(schemes(size(added)))
       do b = 1, size(added); allocate(schemes(b) % scheme, source=bdf_family(2)); end do
       call context % set_linear_solver('iterative'); call context % set_storage('sparse')
       call context % set_preconditioner('multigrid'); call context % set_linear_limits(30, 2, 200)
    case default
       error stop 'concurrent_executions: a case of the set is selected'
    end select
    physics = van_der_pol(degree)
    functional = [van_der_pol_energy(degree)]
    initial = consistent_state(physics, degree + 1, lower, design, context=context)
    if (case_index == streamed_case) then
       call execution % initialize(schemes, added, physics, degree + 1, uniform_grid(duration), design, initial, &
            & startup=startup, functionals=functional, derivative_order=highest_order, context=context)
    else if (seed > 0) then
       call execution % initialize(schemes, added, physics, degree + 1, random_grid(duration, seed), design, initial, &
            & startup=startup, context=context)
    else
       call execution % initialize(schemes, added, physics, degree + 1, uniform_grid(duration), design, initial, &
            & startup=startup, context=context)
    end if
  end subroutine initialize_case

  subroutine finish(execution)
    type(chain_execution), intent(inout) :: execution
    integer :: advances
    advances = 0
    do while (.not. execution % complete())
       call execution % advance()
       advances = advances + 1
       if (advances > 20) error stop 'concurrent_executions: a finite chain reaches completion'
    end do
  end subroutine finish

  !-------------------------------------------------------------------!
  ! Forward and reverse derivatives of the energy functional, the
  ! orders and passes selected by the case; a non-converged outcome
  ! is recorded and does not stop the process.
  !-------------------------------------------------------------------!
  subroutine differentiate(execution, case_index, record)
    type(chain_execution), intent(inout) :: execution
    integer, intent(in) :: case_index
    type(execution_record), intent(inout) :: record
    type(expression) :: functional(1)
    real(dp), allocatable :: table(:,:)
    type(solve_result) :: outcome
    integer :: order, top
    logical :: reverse
    top = highest_order
    reverse = .true.
    select case (case_index)
    case (3, failing_case); top = 1; reverse = .false.
    case (4); top = 1
    case (9); top = 2; reverse = .false.
    end select
    functional = [van_der_pol_energy(physics_degree(case_index))]
    do order = 1, top
       call execution % derivative(functional, order, forward_pass, table, outcome=outcome)
       call noted(outcome, record)
       record % forward(order) = table(1, 1)
       if (reverse) then
          call execution % derivative(functional, order, reverse_pass, table, outcome=outcome)
          call noted(outcome, record)
          record % reverse(order) = table(1, 1)
       end if
    end do
  end subroutine differentiate

  subroutine noted(outcome, record)
    type(solve_result), intent(in) :: outcome
    type(execution_record), intent(inout) :: record
    if (outcome % converged()) return
    record % derivative_converged = .false.
    if (record % failure_reason == 0) record % failure_reason = outcome % reason
  end subroutine noted

  pure integer function physics_degree(case_index)
    integer, intent(in) :: case_index
    select case (case_index)
    case (2, 6); physics_degree = 3
    case default; physics_degree = 2
    end select
  end function physics_degree

  subroutine extract(execution, case_index, record)
    type(chain_execution), intent(inout) :: execution
    integer, intent(in) :: case_index
    type(execution_record), intent(inout) :: record
    type(chain_block), allocatable :: chain(:)
    type(expansion), allocatable, target :: tower
    type(imbalance) :: final_imbalance
    type(tally) :: account
    real(dp) :: achieved
    integer :: b, level, order, event
    account = execution % account()
    if (case_index == streamed_case) then
       call execution % take_results(chain, tower, record % steps, record % times, achieved, &
            & final_imbalance=final_imbalance, f=record % functional, tower_storage=record % tower_storage, &
            & state_storage=record % state_storage)
       record % states = [real(dp) ::]
    else
       call execution % take_results(chain, tower, record % steps, record % times, achieved, &
            & final_imbalance=final_imbalance)
       record % states = [real(dp) ::]
       do b = 1, size(chain)
          record % states = [record % states, chain(b) % state]
       end do
       allocate(record % functional(0, 0))
    end if
    record % converged = final_imbalance % converged
    allocate(record % counts(account % num_levels(), 0:highest_order, tally_num_events()), source=0.0_dp)
    do event = 1, tally_num_events()
       if (event == elapsed_time) cycle
       do order = 0, highest_order
          do level = 1, account % num_levels()
             record % counts(level, order, event) = account % amount(level, order, event)
          end do
       end do
    end do
  end subroutine extract

  !-------------------------------------------------------------------!
  ! Every quantity of one execution equals the reference exactly.
  !-------------------------------------------------------------------!
  subroutine compare_exact(actual, expected, case_index)
    type(execution_record), intent(in) :: actual, expected
    integer, intent(in) :: case_index
    character(len=48) :: which
    write(which, '(a,i0)') 'execution ', case_index
    call check(actual % converged .eqv. expected % converged, trim(which) // ': convergence equals the reference')
    call check(actual % derivative_converged .eqv. expected % derivative_converged .and. &
         & actual % failure_reason == expected % failure_reason, trim(which) // ': derivative outcome equals the reference')
    call check(size(actual % states) == size(expected % states), trim(which) // ': trajectory extents agree')
    call check(all(actual % states == expected % states), trim(which) // ': states equal the reference exactly')
    call check(all(actual % steps == expected % steps) .and. all(actual % times == expected % times), &
         & trim(which) // ': the grid equals the reference exactly')
    call check(all(actual % forward == expected % forward) .and. all(actual % reverse == expected % reverse), &
         & trim(which) // ': derivative tables equal the reference exactly')
    call check(all(shape(actual % functional) == shape(expected % functional)), trim(which) // ': Taylor extents agree')
    call check(all(actual % functional == expected % functional), trim(which) // ': Taylor coefficients equal exactly')
    call check(all(actual % tower_storage == expected % tower_storage) .and. &
         & all(actual % state_storage == expected % state_storage), trim(which) // ': storage pairs agree')
    call check(all(actual % counts == expected % counts), trim(which) // ': counted accounting events agree')
  end subroutine compare_exact

  !-------------------------------------------------------------------!
  ! Each execution's own file, named by its index: written inside the
  ! parallel region through its own unit.
  !-------------------------------------------------------------------!
  subroutine written(record, case_index, directory)
    type(execution_record), intent(in) :: record
    integer, intent(in) :: case_index
    character(len=*), intent(in) :: directory
    character(len=320) :: path
    integer :: unit, i, order
    write(path, '(a,a,i0,a)') trim(directory), '/execution_', case_index, '.txt'
    open(newunit=unit, file=trim(path), status='replace', action='write')
    write(unit, '(a,l1,a,l1,a,i0)') 'converged ', record % converged, ' derivative_converged ', &
         & record % derivative_converged, ' failure_reason ', record % failure_reason
    do i = 1, size(record % states)
       write(unit, '(a,i0,a,es24.16)') 'state ', i, ' ', record % states(i)
    end do
    do i = 1, size(record % steps)
       write(unit, '(a,i0,a,es24.16,a,es24.16)') 'grid ', i, ' ', record % steps(i), ' ', record % times(i)
    end do
    do order = 1, highest_order
       write(unit, '(a,i0,a,es24.16,a,es24.16)') 'derivative ', order, ' forward ', record % forward(order), &
            & ' reverse ', record % reverse(order)
    end do
    do order = 0, size(record % functional, 1) - 1
       write(unit, '(a,i0,a,es24.16)') 'taylor ', order, ' ', record % functional(order, 1)
    end do
    write(unit, '(a,4(1x,i0))') 'storage', record % tower_storage, record % state_storage
    write(unit, '(a,*(1x,i0))') 'counts', nint(record % counts)
    close(unit)
  end subroutine written

  subroutine summarised(case_index, record)
    integer, intent(in) :: case_index
    type(execution_record), intent(in) :: record
    write(*,'(a,i0,a,l1,a,l1,a,i0,a,i0,a,i0)') '   execution ', case_index, ': converged ', record % converged, &
         & ' derivatives converged ', record % derivative_converged, ' states ', size(record % states), &
         & ' newton solves ', nint(sum(record % counts(:, :, 5))), ' linear solves ', nint(sum(record % counts(:, :, 6)))
  end subroutine summarised

  real(dp) function total_after_join() result(total)
    integer :: k
    total = 0.0_dp
    do k = 1, num_executions
       total = total + reference(k) % forward(1)
    end do
  end function total_after_join

  !-------------------------------------------------------------------!
  ! Throughput and peak memory: the set run num_repetitions times at
  ! the requested thread count; every line is written after each join.
  !-------------------------------------------------------------------!
  subroutine measured()
    type(execution_record) :: records(num_executions)
    integer(8) :: started, finished, rate
    real(dp) :: seconds
    integer :: k, r
    do r = 1, num_repetitions
       call system_clock(started, rate)
       !$omp parallel do schedule(dynamic,1) default(shared) private(k)
       do k = 1, num_executions
          call run_case(k, records(k), trim(out_dir) // '/measure')
       end do
       !$omp end parallel do
       call system_clock(finished)
       seconds = real(finished - started, dp) / real(rate, dp)
       write(*,'(a,i0,a,i0,a,f10.4,a,f10.3,a,i0)') ' MEASURE threads=', threads_in_use, ' repetition=', r, &
            & ' seconds=', seconds, ' executions_per_second=', real(num_executions, dp) / seconds, &
            & ' vmhwm_kb=', peak_resident_kb()
    end do
  end subroutine measured

  integer function peak_resident_kb() result(kb)
    character(len=128) :: line
    integer :: unit, status
    kb = 0
    open(newunit=unit, file='/proc/self/status', action='read', status='old', iostat=status)
    if (status /= 0) return
    do
       read(unit, '(a)', iostat=status) line
       if (status /= 0) exit
       if (line(1:6) == 'VmHWM:') then
          read(line(7:), *) kb
          exit
       end if
    end do
    close(unit)
  end function peak_resident_kb

  subroutine check(satisfied, description)
    logical, intent(in) :: satisfied
    character(len=*), intent(in) :: description
    if (.not. satisfied) then
       print *, 'FAIL: ', description
       error stop 'concurrent_executions: a law failed'
    end if
  end subroutine check

end program concurrent_executions
