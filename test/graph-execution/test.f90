program graph_execution
  use util_precision, only : dp
  use operation_action, only : operation
  use operation_driver, only : driver, pairing, rule_graph, data_graph
  use operation_temporal_minimization, only : temporal_minimizer
  use operation_minimization, only : solve_result, SOLVE_NOT_STARTED, SOLVE_EVALUATED
  use view_read_write, only : bipartite_digraph, FIRST_PART, SECOND_PART
  use view_directed_stored, only : stored_directed_graph
  use field_calculus, only : field
  use field_stored, only : stored_field
  use linear_rules, only : linear_rule
  implicit none

  type(stored_directed_graph) :: domain
  integer :: num_failures = 0

  domain = stored_directed_graph(1, [integer ::], [integer ::])
  call check_interleaving()
  call check_lifetimes()
  call check_absent_rule()
  call check_argument_identity()
  call check_output_absence()
  call check_replay()
  call check_temporal()
  call check_empty()
  if (num_failures /= 0) error stop 'graph execution checks failed'
  print *, 'All graph execution checks passed.'

contains

  subroutine assert_all(satisfied, description)
    logical, intent(in) :: satisfied
    character(len=*), intent(in) :: description
    if (satisfied) then
       print *, 'PASS : ', description
    else
       num_failures = num_failures + 1
       print *, 'FAIL : ', description
    end if
  end subroutine assert_all

  function sequence_incidence(order) result(incidence)
    integer, intent(in) :: order(:)
    type(bipartite_digraph) :: incidence
    integer, allocatable :: from_part(:), from_vertex(:), to_part(:), to_vertex(:)
    integer :: i, n
    n = size(order)
    allocate(from_part(2*n), from_vertex(2*n), to_part(2*n), to_vertex(2*n))
    do i = 1, n
       from_part(2*i-1:2*i) = [SECOND_PART, FIRST_PART]
       from_vertex(2*i-1:2*i) = [i, order(i)]
       to_part(2*i-1:2*i) = [FIRST_PART, SECOND_PART]
       to_vertex(2*i-1:2*i) = [order(i), i+1]
    end do
    incidence = bipartite_digraph(n, n+1, from_part, from_vertex, to_part, to_vertex)
  end function sequence_incidence

  function sequence_pairing(num_rules, initial_value) result(connection)
    integer, intent(in) :: num_rules
    real(dp), intent(in) :: initial_value
    type(pairing) :: connection
    type(rule_graph) :: rules
    type(data_graph) :: values
    type(stored_field) :: initial
    integer :: i
    allocate(rules % at(num_rules), values % at(num_rules+1))
    do i = 1, num_rules
       allocate(rules % at(i) % rule, source=linear_rule([2.0_dp], 1.0_dp))
    end do
    initial = stored_field('initial value', domain % vertex_set(), 1)
    call initial % set_real_vector([initial_value])
    allocate(values % at(1) % datum, source=initial)
    connection = rules % pair(values)
  end function sequence_pairing

  real(dp) function datum_value(connection, vertex) result(value)
    type(pairing), intent(in) :: connection
    integer, intent(in) :: vertex
    class(field), allocatable :: datum
    real(dp), allocatable :: values(:)
    call connection % datum_at(vertex, datum)
    if (.not. allocated(datum)) error stop 'the expected datum was not written'
    call datum % real_vector(values)
    value = values(1)
  end function datum_value

  logical function written(connection, vertex)
    type(pairing), intent(in) :: connection
    integer, intent(in) :: vertex
    class(field), allocatable :: datum
    call connection % datum_at(vertex, datum)
    written = allocated(datum)
  end function written

  subroutine check_interleaving()
    type(driver) :: first, second, copied
    type(bipartite_digraph) :: incidence
    type(pairing) :: connection
    class(operation), allocatable :: rule
    integer, target :: first_evaluations, second_evaluations
    integer :: executed

    first_evaluations = 0
    second_evaluations = 0
    incidence = sequence_incidence([3,1,2])
    first = driver(linear_rule([2.0_dp], 1.0_dp), incidence)
    second = first
    call first % pair_with(sequence_pairing(3, 1.0_dp))
    call second % pair_with(sequence_pairing(3, 10.0_dp))
    call assert_all(first % next_rule() == 3 .and. first % num_completed() == 0, &
         & 'execution starts at the dependency order, not the first vertex identity')

    call first % set_rule(3, linear_rule([2.0_dp], 1.0_dp, first_evaluations))
    call first % advance(domain, executed)
    call assert_all(executed == 3 .and. first % num_completed() == 1 .and. first % next_rule() == 1, &
         & 'one advance executes exactly one scheduled vertex')
    call assert_all(second % num_completed() == 0 .and. second % next_rule() == 3, &
         & 'advancing one execution does not advance another')
    call first % clear_rule(3)
    connection = first % pairing_of()
    call connection % rule_at(3, rule)
    call assert_all(.not. allocated(rule), 'completed transient rules can be removed without rebuilding the pairing')
    copied = first

    call second % set_rule(3, linear_rule([2.0_dp], 1.0_dp, second_evaluations))
    call second % advance(domain)
    call second % clear_rule(3)
    call assert_all(datum_value(first % pairing_of(), 2) == 3.0_dp .and. &
         & datum_value(second % pairing_of(), 2) == 21.0_dp, 'interleaved executions own independent values')

    call first % set_rule(1, linear_rule([2.0_dp], 3.0_dp, first_evaluations))
    call second % set_rule(1, linear_rule([3.0_dp], 0.0_dp, second_evaluations))
    call first % advance(domain)
    call second % advance(domain)
    call first % clear_rule(1)
    call second % clear_rule(1)
    call assert_all(copied % num_completed() == 1 .and. datum_value(copied % pairing_of(), 2) == 3.0_dp, &
         & 'copying an execution preserves its position and independent data')

    call first % set_rule(2, linear_rule([0.5_dp], 0.0_dp, first_evaluations))
    call second % set_rule(2, linear_rule([1.0_dp], 1.0_dp, second_evaluations))
    call first % advance(domain)
    call second % advance(domain)
    call first % clear_rule(2)
    call second % clear_rule(2)
    call assert_all(first % complete() .and. second % complete() .and. first % next_rule() == 0, &
         & 'all scheduled positions exhaust both executions independently')
    call assert_all(datum_value(first % pairing_of(), 4) == 4.5_dp .and. &
         & datum_value(second % pairing_of(), 4) == 64.0_dp, 'single-rule replacement changes only the intended execution')
    call first % advance(domain, executed)
    call assert_all(executed == 0 .and. first_evaluations == 3 .and. second_evaluations == 3, &
         & 'advancing a completed execution performs no additional operation')

    call copied % set_rule(1, linear_rule([1.0_dp], 1.0_dp))
    call copied % advance(domain)
    call copied % advance(domain)
    call assert_all(datum_value(copied % pairing_of(), 4) == 9.0_dp .and. &
         & datum_value(first % pairing_of(), 4) == 4.5_dp, 'a copied partial execution can continue independently')
  end subroutine check_interleaving

  subroutine check_lifetimes()
    type(bipartite_digraph) :: incidence
    type(driver) :: schedule
    type(rule_graph) :: rules
    type(data_graph) :: values
    type(stored_field) :: initial
    integer, allocatable :: released(:)
    integer :: executed, num_released, step

    ! Two forward rules followed by two reverse rules. The first source
    ! remains a declared input of the final reverse rule.
    incidence = bipartite_digraph(4, 6, &
         & [2,1,2,1,2,2,2,1,2,2,1], [1,1,2,2,4,2,3,3,5,1,4], &
         & [1,2,1,2,1,1,1,2,1,1,2], [1,2,2,3,3,3,3,5,4,4,6])
    allocate(rules % at(4), values % at(6))
    allocate(rules % at(1) % rule, source=linear_rule([2.0_dp], 0.0_dp))
    allocate(rules % at(2) % rule, source=linear_rule([2.0_dp], 0.0_dp))
    allocate(rules % at(3) % rule, source=linear_rule([1.0_dp,10.0_dp,0.0_dp], 0.0_dp))
    allocate(rules % at(4) % rule, source=linear_rule([1.0_dp,100.0_dp], 0.0_dp))
    initial = stored_field('source', domain % vertex_set(), 1)
    call initial % set_real_vector([3.0_dp])
    allocate(values % at(1) % datum, source=initial)
    call initial % set_real_vector([1.0_dp])
    allocate(values % at(4) % datum, source=initial)
    schedule = driver(linear_rule([1.0_dp], 0.0_dp), incidence)
    call schedule % pair_with(rules % pair(values))
    call assert_all(schedule % last_reader_of(1) == 4, 'the initial source remains live through its reverse dependency')
    do step = 1, 3
       call schedule % advance(domain, executed)
       call assert_all(executed == step .and. written(schedule % pairing_of(), 1), &
            & 'earlier advances preserve a source with a later reader')
    end do
    released = schedule % released_at(3)
    call assert_all(all(released == [2,3,4]), 'one release interval contains only that position''s final reads')
    call assert_all(.not. written(schedule % pairing_of(), 2) .and. &
         & .not. written(schedule % pairing_of(), 3) .and. .not. written(schedule % pairing_of(), 4), &
         & 'forward and terminal data are released after their final reverse use')
    call schedule % advance(domain)
    call assert_all(.not. written(schedule % pairing_of(), 1) .and. &
         & datum_value(schedule % pairing_of(), 6) == 361.0_dp, 'the final reader consumes the initial source before its release')
    num_released = 0
    do step = 1, schedule % num_completed()
       released = schedule % released_at(step)
       num_released = num_released + size(released)
    end do
    call assert_all(num_released == 5, 'each consumed datum occurs in exactly one release interval')
  end subroutine check_lifetimes

  subroutine check_absent_rule()
    type(driver) :: schedule
    type(pairing) :: connection
    integer :: executed
    schedule = driver(linear_rule([1.0_dp], 0.0_dp), sequence_incidence([1,2,3]))
    connection = sequence_pairing(3, 2.0_dp)
    call connection % clear_rule(3)
    call schedule % pair_with(connection)
    call schedule % advance(domain)
    call schedule % advance(domain)
    call assert_all(written(schedule % pairing_of(), 3), 'a datum remains live until an absent rule''s declared read')
    call schedule % advance(domain, executed)
    call assert_all(executed == 3 .and. schedule % complete() .and. &
         & .not. written(schedule % pairing_of(), 3) .and. .not. written(schedule % pairing_of(), 4), &
         & 'an absent rule occupies its scheduled position and resolves its final reads without writing')
  end subroutine check_absent_rule

  subroutine check_argument_identity()
    type(driver) :: schedule
    type(rule_graph) :: rules
    type(data_graph) :: values
    type(stored_field) :: initial
    type(bipartite_digraph) :: incidence
    incidence = bipartite_digraph(1, 3, [2,2,1], [1,2,1], [1,1,2], [1,1,3])
    allocate(rules % at(1), values % at(3))
    allocate(rules % at(1) % rule, source=linear_rule([10.0_dp,1.0_dp], 0.0_dp))
    initial = stored_field('second input', domain % vertex_set(), 1)
    call initial % set_real_vector([4.0_dp])
    allocate(values % at(2) % datum, source=initial)
    schedule = driver(linear_rule([10.0_dp,1.0_dp], 0.0_dp), incidence)
    call schedule % pair_with(rules % pair(values))
    call schedule % advance(domain)
    call assert_all(datum_value(schedule % pairing_of(), 3) == 4.0_dp, &
         & 'an unwritten input remains unbound without changing other argument identities')
  end subroutine check_argument_identity

  subroutine check_output_absence()
    type(driver) :: schedule
    type(rule_graph) :: rules
    type(data_graph) :: values
    type(bipartite_digraph) :: incidence
    incidence = bipartite_digraph(2, 2, [1,1], [1,2], [2,2], [1,2])
    allocate(rules % at(2), values % at(2))
    allocate(rules % at(1) % rule, source=linear_rule([real(dp) ::], 5.0_dp))
    allocate(rules % at(2) % rule, source=linear_rule([real(dp) ::], 0.0_dp, defined=.false.))
    schedule = driver(linear_rule([real(dp) ::], 0.0_dp), incidence)
    call schedule % pair_with(rules % pair(values))
    call schedule % evaluate(domain)
    call assert_all(datum_value(schedule % pairing_of(), 1) == 5.0_dp .and. &
         & .not. written(schedule % pairing_of(), 2), 'a rule without output does not reuse a previous rule''s value')
  end subroutine check_output_absence

  subroutine check_replay()
    type(driver) :: schedule
    type(rule_graph) :: rules
    type(data_graph) :: values
    type(pairing) :: connection
    integer, target :: evaluations
    evaluations = 0
    allocate(rules % at(1), values % at(1))
    allocate(rules % at(1) % rule, source=linear_rule([real(dp) ::], 7.0_dp, evaluations))
    connection = rules % pair(values)
    schedule = driver(linear_rule([real(dp) ::], 0.0_dp), &
         & bipartite_digraph(1, 1, [1], [1], [2], [1]))
    call schedule % pair_with(connection)
    call schedule % advance(domain)
    call schedule % evaluate(domain)
    call schedule % evaluate(domain)
    call assert_all(evaluations == 3 .and. schedule % num_completed() == 1, &
         & 'evaluate replays the complete schedule using the incremental kernel')
    call schedule % pair_with(connection)
    call assert_all(schedule % num_completed() == 0 .and. .not. schedule % complete(), &
         & 'newly paired data begin a new execution')
  end subroutine check_replay

  subroutine check_temporal()
    type(temporal_minimizer) :: first, second
    type(driver) :: schedule
    type(solve_result) :: outcome
    real(dp) :: achieved, rhs(0), x(0)
    integer :: executed
    schedule = driver(linear_rule([2.0_dp], 1.0_dp), sequence_incidence([2,1]))
    call first % state(schedule, domain, domain % vertex_set(), 0)
    call second % state(schedule, domain, domain % vertex_set(), 0)
    call first % pair_with(sequence_pairing(2, 2.0_dp))
    call second % pair_with(sequence_pairing(2, 7.0_dp))
    call first % set_rule(2, linear_rule([2.0_dp], 1.0_dp))
    call first % advance(executed)
    call first % clear_rule(executed)
    outcome = first % result()
    call assert_all(executed == 2 .and. first % next_rule() == 1 .and. second % next_rule() == 2 .and. &
         & outcome % reason == SOLVE_NOT_STARTED, 'temporal executions delegate independent progress without premature completion')
    call second % solve(rhs, x, achieved)
    call first % advance()
    outcome = first % result()
    call assert_all(first % complete() .and. first % num_completed() == 2 .and. outcome % iterations == 2 .and. &
         & outcome % reason == SOLVE_EVALUATED, 'temporal completion records the number of executed schedule positions')
    call assert_all(datum_value(first % pairing_of(), 3) == 11.0_dp .and. &
         & datum_value(second % pairing_of(), 3) == 31.0_dp .and. achieved == 0.0_dp, &
         & 'incremental and complete temporal solves retain independent results')
  end subroutine check_temporal

  subroutine check_empty()
    type(driver) :: schedule
    type(rule_graph) :: rules
    type(data_graph) :: values
    integer :: executed
    integer, allocatable :: released(:)
    allocate(rules % at(0), values % at(0))
    schedule = driver(linear_rule([real(dp) ::], 0.0_dp), &
         & bipartite_digraph(0, 0, [integer ::], [integer ::], [integer ::], [integer ::]))
    call schedule % pair_with(rules % pair(values))
    call schedule % advance(domain, executed)
    call schedule % evaluate(domain)
    released = schedule % released_at(0)
    call assert_all(schedule % complete() .and. schedule % next_rule() == 0 .and. &
         & schedule % num_completed() == 0 .and. executed == 0 .and. size(released) == 0, &
         & 'an empty execution is complete without an operation or release')
  end subroutine check_empty

end program graph_execution
