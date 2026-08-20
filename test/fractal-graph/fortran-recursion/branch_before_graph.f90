! DECLARATION ORDER PROBE. branch declared before graph, with a
! forward reference by pointer. EXPECTED: compiles. This is the
! representation the kernel uses.
module branch_before_graph
  implicit none
  type :: branch
     integer              :: status = 0
     type(graph), pointer :: known  => null()
  end type branch
  type :: graph
     type(branch) :: branch(2)
  end type graph
end module branch_before_graph
