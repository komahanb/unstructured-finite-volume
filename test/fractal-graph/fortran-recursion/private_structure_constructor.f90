! INVARIANT PROBE. The structure constructor, which would set status
! and reference independently. EXPECTED: rejected - 'status_' is a
! PRIVATE component. The three branch constructors are therefore the
! only introductions of a branch value.
program private_structure_constructor
  use graph_fractal
  implicit none
  type(graph)        :: a
  type(branch) :: x
  call a % declare()
  x = branch(GRAPH_KNOWN, null())
  a % branch(1) = x
end program private_structure_constructor
