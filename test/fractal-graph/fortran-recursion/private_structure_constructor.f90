! INVARIANT PROBE. The structure constructor, which would set status
! and reference independently. EXPECTED: rejected - 'status_' is a
! PRIVATE component. The three branch constructors are therefore the
! only introductions of a branch value.
program private_structure_constructor
  use fractal_graph
  implicit none
  type(graph)        :: a
  type(graph_branch) :: x
  call a % declare()
  x = graph_branch(GRAPH_KNOWN, null())
  a % branch(1) = x
end program private_structure_constructor
