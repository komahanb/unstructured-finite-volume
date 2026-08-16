! NAVIGATION PROBE. The admitted public forms, against the shipped
! kernel. EXPECTED: compiles.
!
!     g % branch(i)                      public component
!     g % branch(i) % status()           branch query
!     g % branch(i) % known()            branch reference
!     p => g % branch(i) % known()       pointer binding, for depth
!     associate (x => ... % known())     ASSOCIATE name, for depth
program encapsulated_navigation
  use fractal_graph
  implicit none
  type(graph), target  :: a, b
  type(graph), pointer :: p
  call a % declare(); call b % declare()
  a % branch(1) = known_branch(b)
  b % branch(1) = known_branch(a)
  if (a % branch(1) % status() /= GRAPH_KNOWN) error stop 1
  if (.not. associated(a % branch(1) % known(), b)) error stop 2
  p => a % branch(1) % known()
  if (.not. associated(p % branch(1) % known(), a)) error stop 3
  associate (x => a % branch(1) % known())
    if (x % branch(1) % status() /= GRAPH_KNOWN) error stop 4
  end associate
end program encapsulated_navigation
