! CANDIDATE B, mechanism 2. A defined assignment on the BRANCH type.
! EXPECTED: compiles and runs, and the assignment IS intercepted - but
! the procedure's only arguments are two branches. No argument, host
! association or inquiry yields the owning graph, so the seal state it
! would have to consult is not reachable from inside.
module seal_branch_assignment_m
  implicit none
  type :: branch_t
     integer :: s = 0
  end type branch_t
  type :: graph_t
     type(branch_t) :: branch(2)
     logical        :: sealed = .false.
  end type graph_t
  interface assignment(=)
     module procedure assign_branch
  end interface assignment(=)
contains
  subroutine assign_branch(lhs, rhs)
    type(branch_t), intent(out) :: lhs
    type(branch_t), intent(in)  :: rhs
    lhs % s = rhs % s
    print *, '   branch-level defined assignment ran; owner not reachable'
  end subroutine assign_branch
end module seal_branch_assignment_m

program seal_branch_assignment
  use seal_branch_assignment_m
  implicit none
  type(graph_t) :: g
  g % sealed = .true.
  g % branch(1) = branch_t(9)
  print *, ' B: intercepted, but the seal could not be consulted: s =', &
       & g % branch(1) % s
end program seal_branch_assignment
