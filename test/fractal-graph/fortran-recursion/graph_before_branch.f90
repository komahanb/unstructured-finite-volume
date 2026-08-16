! DECLARATION ORDER PROBE. graph declared before graph_branch.
! EXPECTED: rejected - 'Derived type at (1) has not been previously
! defined and so cannot appear in a derived type definition'. The
! constraint is on declaration order only; the two-type ontology is
! unaffected.
module graph_before_branch
  implicit none
  type :: graph
     type(graph_branch) :: branch(2)
  end type graph
  type :: graph_branch
     integer              :: status = 0
     type(graph), pointer :: known  => null()
  end type graph_branch
end module graph_before_branch
