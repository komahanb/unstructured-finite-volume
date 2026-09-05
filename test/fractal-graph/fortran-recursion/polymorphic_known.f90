! DECLARATION ORDER PROBE. As branch_before_graph, with a polymorphic
! forward reference. EXPECTED: compiles. Not used: no type extends
! graph, so class(graph) adds no admissible value.
module polymorphic_known
  implicit none
  type :: branch
     integer               :: status = 0
     class(graph), pointer :: known  => null()
  end type branch
  type :: graph
     type(branch) :: branch(2)
  end type graph
end module polymorphic_known
