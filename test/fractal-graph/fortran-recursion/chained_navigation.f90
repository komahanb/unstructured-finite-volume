! NAVIGATION PROBE. A function reference in the middle of a traversal.
! EXPECTED: rejected - 'Unclassifiable statement'; the same rule as
! function_result_base, reported by the parser.
!
! This is a Fortran restriction on function results, not a restriction
! on the ontology: branch(2) remains a public component, and
! g % branch(i) % status() and g % branch(i) % known() remain direct.
! Depth is reached by a pointer or an ASSOCIATE name; see
! encapsulated_navigation.f90.
program chained_navigation
  use graph_fractal
  implicit none
  type(graph), target :: a, b
  integer :: s
  call a % declare(); call b % declare()
  a % branch(1) = known_branch(b)
  s = a % branch(1) % known() % branch(1) % status()
  print *, s
end program chained_navigation
