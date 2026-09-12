!=====================================================================!
! Schedule construction against execution over one dependency graph.
!
! Rules 1..n write data 2..n+1; datum 1 is the initial value. Shapes:
!
!      chain     rule i reads datum i
!      history   rule i reads data max(1, i-r+1)..i          (reach r)
!      tree      rule i reads datum i/2 + 1                  (branching)
!
! Every rule is the linear combination 1 + (1/r) sum of what it reads,
! so the values stay bounded by n and an oracle recomputes them from
! the same arrays in the same operation order. Measured separately:
! incidence construction, driver construction (the schedule), pairing
! with data, and execution by advance until complete.
!
! usage: schedule shape n reach
!=====================================================================!

program schedule_scaling

  use util_precision       , only : dp
  use operation_driver     , only : driver, pairing, rule_graph, data_graph
  use view_read_write      , only : bipartite_digraph, FIRST_PART, SECOND_PART
  use view_directed_stored , only : stored_directed_graph
  use field_calculus       , only : field
  use field_stored         , only : stored_field
  use linear_rules         , only : linear_rule
  use benchmark_measurement, only : phase, begin_phase, end_phase, record, peak_rss_kilobytes, &
       &                            argument_integer, argument_string

  implicit none

  character(len=:), allocatable :: shape, tokens
  character(len=96) :: written
  integer :: n, reach, width, i, k, d, a, num_arcs, num_data, executed, steps, num_failures
  integer, allocatable :: first_read(:), reads(:), from_part(:), from_vertex(:), to_part(:), to_vertex(:)
  integer, allocatable :: order(:), position(:), released(:), release_count(:), num_readers(:)
  real(dp), allocatable :: coefficients(:), oracle(:), values(:)
  type(bipartite_digraph) :: incidence
  type(driver) :: execution
  type(stored_directed_graph) :: domain
  type(rule_graph) :: rules
  type(data_graph) :: data
  type(pairing) :: connection
  type(stored_field) :: initial
  class(field), allocatable :: datum
  type(phase) :: measured
  integer, target :: evaluations
  logical :: satisfied

  shape = argument_string(1, 'chain')
  n     = argument_integer(2, 1000)
  reach = argument_integer(3, 1)
  if (n < 1 .or. reach < 1) error stop 'schedule: n and reach are positive'
  num_data = n + 1
  write(written, '(a,a,a,i0,a,i0)') 'suite=schedule shape=', shape, ' n=', n, ' reach=', reach
  tokens = trim(written)

  ! the reads of every rule, ascending datum
  allocate(first_read(n + 1))
  first_read(1) = 1
  do i = 1, n
     first_read(i + 1) = first_read(i) + num_reads(i)
  end do
  allocate(reads(first_read(n + 1) - 1))
  do i = 1, n
     select case (shape)
     case ('chain')
        reads(first_read(i)) = i
     case ('history')
        do k = 1, num_reads(i)
           reads(first_read(i) + k - 1) = i - num_reads(i) + k
        end do
     case ('tree')
        reads(first_read(i)) = i / 2 + 1
     case default
        error stop 'schedule: the shape is chain, history or tree'
     end select
  end do
  width = 1
  if (shape == 'history') width = reach
  allocate(coefficients(width))
  coefficients = 1.0_dp / real(width, dp)

  num_arcs = size(reads) + n
  allocate(from_part(num_arcs), from_vertex(num_arcs), to_part(num_arcs), to_vertex(num_arcs))
  a = 0
  do i = 1, n
     do k = first_read(i), first_read(i + 1) - 1
        a = a + 1
        from_part(a) = SECOND_PART
        from_vertex(a) = reads(k)
        to_part(a) = FIRST_PART
        to_vertex(a) = i
     end do
     a = a + 1
     from_part(a) = FIRST_PART
     from_vertex(a) = i
     to_part(a) = SECOND_PART
     to_vertex(a) = i + 1
  end do

  call begin_phase(measured)
  incidence = bipartite_digraph(n, num_data, from_part, from_vertex, to_part, to_vertex)
  call end_phase(measured)
  call record(tokens, 'incidence', measured)

  evaluations = 0
  call begin_phase(measured)
  execution = driver(linear_rule(coefficients, 1.0_dp, evaluations), incidence)
  call end_phase(measured)
  call record(tokens, 'schedule', measured)

  domain = stored_directed_graph(1, [integer ::], [integer ::])
  initial = stored_field('initial value', domain % vertex_set(), 1)
  call initial % set_real_vector([1.0_dp])
  call begin_phase(measured)
  allocate(rules % at(n), data % at(num_data))
  do i = 1, n
     allocate(rules % at(i) % rule, source=linear_rule(coefficients, 1.0_dp, evaluations))
  end do
  allocate(data % at(1) % datum, source=initial)
  call execution % pair_with(rules % pair(data))
  call end_phase(measured)
  call record(tokens, 'pairing', measured)

  steps = 0
  call begin_phase(measured)
  do while (.not. execution % complete())
     call execution % advance(domain, executed)
     steps = steps + 1
  end do
  call end_phase(measured)
  call record(tokens, 'execution', measured)

  ! ---- verification against independent constructions ----
  num_failures = 0

  ! the visiting order is a permutation with every writer before its readers
  order = execution % visits()
  allocate(position(n), source=0)
  satisfied = size(order) == n
  if (satisfied) then
     do k = 1, n
        if (order(k) < 1 .or. order(k) > n) then
           satisfied = .false.
           exit
        end if
        position(order(k)) = k
     end do
     satisfied = satisfied .and. all(position > 0)
  end if
  if (satisfied) then
     do i = 1, n
        do k = first_read(i), first_read(i + 1) - 1
           d = reads(k)
           if (d > 1) satisfied = satisfied .and. position(d - 1) < position(i)
        end do
     end do
  end if
  call verified(satisfied, 'topological_order')

  ! every datum with a reader is released exactly once, at its last reader
  allocate(num_readers(num_data), source=0)
  allocate(release_count(num_data), source=0)
  do k = 1, size(reads)
     num_readers(reads(k)) = num_readers(reads(k)) + 1
  end do
  satisfied = .true.
  do k = 1, n
     released = execution % released_at(k)
     do i = 1, size(released)
        d = released(i)
        release_count(d) = release_count(d) + 1
        satisfied = satisfied .and. execution % last_reader_of(d) == k
     end do
  end do
  satisfied = satisfied .and. all((release_count == 1) .eqv. (num_readers > 0))
  do i = 1, n
     do k = first_read(i), first_read(i + 1) - 1
        satisfied = satisfied .and. execution % last_reader_of(reads(k)) >= position(i)
     end do
  end do
  call verified(satisfied, 'release_intervals')

  ! the values, recomputed in the rule's operation order
  allocate(oracle(num_data))
  oracle(1) = 1.0_dp
  do i = 1, n
     oracle(i + 1) = 1.0_dp
     do k = first_read(i), first_read(i + 1) - 1
        oracle(i + 1) = oracle(i + 1) + coefficients(k - first_read(i) + 1) * oracle(reads(k))
     end do
  end do
  connection = execution % pairing_of()
  satisfied = .true.
  a = 0
  do d = 1, num_data
     call connection % datum_at(d, datum)
     if (.not. allocated(datum)) then
        satisfied = satisfied .and. num_readers(d) > 0
        cycle
     end if
     a = a + 1
     call datum % real_vector(values)
     satisfied = satisfied .and. size(values) == 1
     if (size(values) == 1) satisfied = satisfied .and. values(1) == oracle(d)
     satisfied = satisfied .and. num_readers(d) == 0
  end do
  call verified(satisfied .and. a > 0, 'values_equal_oracle')
  call verified(evaluations == n .and. steps == n, 'one_evaluation_per_rule')

  write(*, '(a,1x,a,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0)') 'summary', tokens, &
       & 'arcs=', num_arcs, 'reads=', size(reads), 'steps=', steps, 'evaluations=', evaluations, &
       & 'outputs=', a, 'peak_rss_kilobytes=', peak_rss_kilobytes()
  if (num_failures > 0) error stop 'schedule: a verification failed'

contains

  integer function num_reads(rule)
    integer, intent(in) :: rule
    num_reads = 1
    if (shape == 'history') num_reads = min(rule, reach)
  end function num_reads

  subroutine verified(condition, name)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: name
    if (.not. condition) num_failures = num_failures + 1
    write(*, '(a,1x,a,1x,a,a,1x,a,l1)') 'verification', tokens, 'name=', name, 'satisfied=', condition
  end subroutine verified

end program schedule_scaling
