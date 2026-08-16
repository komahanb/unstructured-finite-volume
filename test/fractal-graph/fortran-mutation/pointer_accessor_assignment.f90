! CANDIDATE D, mutation. Assignment through a pointer-returning
! accessor. EXPECTED: compiles - a pointer function result is a
! definable variable. D therefore enforces nothing, while costing the
! same navigation as C.
program pointer_accessor_assignment
  use accessor_candidates
  implicit none
  type(graph_t), target   :: g
  type(branch_t), pointer :: p
  g % by_pointer(1) = branch_t(9)
  p => g % by_pointer(1)
  if (p % status() /= 9) error stop 'assignment through the accessor had no effect'
  print *, ' D: private branch, mutated through the accessor: s =', p % status()
end program pointer_accessor_assignment
