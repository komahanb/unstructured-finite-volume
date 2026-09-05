! CANDIDATE C, mutation. Assignment through a value-returning accessor.
! EXPECTED: rejected - 'The function result on the lhs of the
! assignment at (1) must have the pointer attribute'. C therefore does
! enforce immutability of the branch slot.
program value_accessor_assignment
  use accessor_candidates
  implicit none
  type(graph_t) :: g
  g % by_value(1) = branch_t(9)
end program value_accessor_assignment
