! INVARIANT PROBE. External rebinding of the branch reference.
! EXPECTED: rejected - 'known_' is a PRIVATE component.
! Admitting it would allow KNOWN with known disassociated.
program private_reference_write
  use graph_fractal
  implicit none
  type(graph), target :: a, b
  call a % declare(); call b % declare()
  a % branch(1) = known_branch(b)
  nullify(a % branch(1) % known_)
end program private_reference_write
