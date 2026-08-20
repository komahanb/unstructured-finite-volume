! Q8.E, against the shipped kernel. A graph reached through known() is
! mutated by the caller. EXPECTED: compiles and runs; the mutation is
! visible on the referenced graph, whose identity does not change.
!
! Immutability is therefore not a property one graph can hold: it would
! have to be refused independently by every graph in the reachable set.
program transitive_mutation
  use graph_fractal
  implicit none
  type(graph), target  :: a, b, c
  type(graph), pointer :: p
  integer              :: before, after

  call a % declare(); call b % declare(); call c % declare()
  a % branch(1) = known_branch(b)
  b % branch(1) = unknown_branch()

  p => a % branch(1) % known()
  before = p % branch(1) % status()
  p % branch(1) = known_branch(c)
  after  = b % branch(1) % status()

  if (before /= GRAPH_UNKNOWN) error stop 'b did not start UNKNOWN'
  if (after  /= GRAPH_KNOWN)   error stop 'the transitive mutation did not take'
  if (.not. b % same_as(b))    error stop 'b lost its identity'
  print *, ' E: b mutated through a pointer obtained from a; identity intact'
end program transitive_mutation
