program execution_refusal
  use util_precision, only : dp
  use operation_driver, only : driver, rule_graph, data_graph
  use view_read_write, only : bipartite_digraph
  use view_directed_stored, only : stored_directed_graph
  use linear_rules, only : linear_rule, placed_rule
  implicit none
  type(driver) :: schedule
  type(stored_directed_graph) :: domain
  type(rule_graph) :: rules
  type(data_graph) :: values
  type(placed_rule) :: rule
  character(len=32) :: case_name

  call get_command_argument(1, case_name)
  domain = stored_directed_graph(1, [integer ::], [integer ::])
  schedule = driver(linear_rule([real(dp) ::], 0.0_dp), &
       & bipartite_digraph(1, 1, [1], [1], [2], [1]))
  allocate(rules % at(1), values % at(1))
  select case(trim(case_name))
  case('unpaired_advance')
     call schedule % advance(domain)
  case('unpaired_advance_with')
     call schedule % advance_with(domain, rule)
  case('unpaired_set')
     call schedule % set_rule(1, linear_rule([real(dp) ::], 0.0_dp))
  case('unpaired_clear')
     call schedule % clear_rule(1)
  case('outside_set')
     call schedule % pair_with(rules % pair(values))
     call schedule % set_rule(2, linear_rule([real(dp) ::], 0.0_dp))
  case('outside_clear')
     call schedule % pair_with(rules % pair(values))
     call schedule % clear_rule(0)
  case('outside_position')
     call schedule % pair_with(rules % pair(values), position=2)
  case('outside_live')
     call schedule % pair_with(rules % pair(values))
     print *, schedule % live_after(-1)
  case default
     error stop 'unknown refusal case'
  end select
  print *, 'Invalid execution was admitted.'
end program execution_refusal
