! INVARIANT PROBE. External assignment of the branch status.
! EXPECTED: rejected - 'status_' is a PRIVATE component.
! Admitting it would allow status = KNOWN with known disassociated.
program private_status_write
  use fractal_graph
  implicit none
  type(graph), target :: a
  call a % declare()
  a % branch(1) % status_ = GRAPH_KNOWN
end program private_status_write
