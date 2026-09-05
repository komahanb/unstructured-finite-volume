! CANDIDATE B, mechanism 1. A defined assignment on the parent type,
! intended to intercept assignment to one of its components and refuse
! it when the parent is sealed. EXPECTED: compiles and runs, and the
! interception DOES NOT HAPPEN - component assignment does not invoke
! the parent's defined assignment. Printed evidence, not a refusal.
module seal_parent_assignment_m
  implicit none
  type :: branch_t
     integer :: s = 0
  end type branch_t
  type :: graph_t
     type(branch_t) :: branch(2)
     logical        :: sealed = .false.
  end type graph_t
  interface assignment(=)
     module procedure assign_graph
  end interface assignment(=)
contains
  subroutine assign_graph(lhs, rhs)
    type(graph_t), intent(out) :: lhs
    type(graph_t), intent(in)  :: rhs
    lhs % branch = rhs % branch
    lhs % sealed = rhs % sealed
    print *, '   graph-level defined assignment ran (whole-object only)'
  end subroutine assign_graph
end module seal_parent_assignment_m

program seal_parent_assignment
  use seal_parent_assignment_m
  implicit none
  type(graph_t) :: g, h
  g % sealed = .true.
  g % branch(1) = branch_t(9)
  if (g % branch(1) % s /= 9) error stop 'component assignment was intercepted'
  print *, ' A: sealed parent, component assigned anyway: s =', g % branch(1) % s
  h = g
end program seal_parent_assignment
