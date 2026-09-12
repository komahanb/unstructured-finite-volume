program execution_independence

  use util_precision, only : dp
  use operation_minimization, only : relative, by_count
  use operation_family, only : bdf_family
  use operation_grid, only : uniform_grid
  use operation_expression, only : expression
  use gti_physics, only : van_der_pol, van_der_pol_energy
  use gti_expansion, only : expansion, family_container
  use gti_march, only : march_context, imbalance, consistent_state
  use gti_sweeps, only : forward_pass, reverse_pass
  use gti_chain, only : chain_execution, chain_block

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

contains

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
    call context % set_stopping(1.0e-11_dp, relative, by_count, 100)
    if (case_index == 1) then
       call context % set_linear_solver('direct')
       call context % set_storage('dense')
       call context % set_time_coupling('sequential')
    else
       call context % set_linear_solver('iterative')
       call context % set_storage('sparse')
       call context % set_preconditioner('gauss_seidel')
       call context % set_linear_limits(40, 2, 200)
       call context % set_time_coupling('coupled')
    end if
  end subroutine configure

  subroutine initialize_case(execution, case_index, context, streamed)
    type(chain_execution), intent(inout) :: execution
    integer, intent(in) :: case_index
    type(march_context), intent(inout) :: context
    logical, intent(in), optional :: streamed
    type(family_container), allocatable :: schemes(:)
    type(expression) :: physics, functional(1)
    integer, allocatable :: added(:)
    real(dp), allocatable :: initial(:), lower(:)
    real(dp) :: duration, design
    integer :: degree, order, b
    logical :: stream
    stream = .false.
    if (present(streamed)) stream = streamed
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
    if (stream) then
       call execution % initialize(schemes, added, physics, degree + 1, uniform_grid(duration), design, initial, &
            & startup=4, functionals=functional, derivative_order=2, context=context)
    else
       call execution % initialize(schemes, added, physics, degree + 1, uniform_grid(duration), design, initial, &
            & startup=4, context=context)
    end if
  end subroutine initialize_case

  subroutine finish(execution)
    type(chain_execution), intent(inout) :: execution
    integer :: advances
    advances = 0
    do while (.not. execution % complete())
       call execution % advance()
       advances = advances + 1
       call check(advances <= 10, 'a finite chain reaches completion')
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
