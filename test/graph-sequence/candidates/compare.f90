! The same three consumers written against each candidate. What is
! counted is not the module but the CALL SITE: how much the consumer
! must know about emptiness.
program compare
  use graph_fractal, only : graph, GRAPH_NULL, null_branch, known_branch
  use branch_form  , only : b_size => size_of
  use graph_form   , only : g_size => size_of, g_empty => is_empty
  implicit none
  type(graph), target :: holder, c1, c2, c3
  integer :: n

  call holder % declare(); call c1 % declare()
  call c2 % declare(); call c3 % declare()

  ! empty
  holder % branch(1) = null_branch()
  n = b_size(holder % branch(1))                       ! branch form: one call
  print *, 'branch form, empty  =', n
  if (g_empty(holder % branch(1))) then                ! graph form: guard
     n = 0
  else
     n = g_size(holder % branch(1) % known())
  end if
  print *, 'graph  form, empty  =', n

  ! length 3
  c3 % branch(1) = known_branch(c3); c3 % branch(2) = null_branch()
  c2 % branch(1) = known_branch(c2); c2 % branch(2) = known_branch(c3)
  c1 % branch(1) = known_branch(c1); c1 % branch(2) = known_branch(c2)
  holder % branch(1) = known_branch(c1)

  n = b_size(holder % branch(1))
  print *, 'branch form, len 3  =', n
  if (g_empty(holder % branch(1))) then
     n = 0
  else
     n = g_size(holder % branch(1) % known())
  end if
  print *, 'graph  form, len 3  =', n
end program compare
