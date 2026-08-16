! CANDIDATE 1. The sequence is a graph_branch.
!
!     NULL          the empty sequence
!     KNOWN -> cell a nonempty sequence
!     UNKNOWN       the sequence is not known
!
! The empty sequence needs no graph, because it has no cell.
module branch_form
  use fractal_graph, only : graph, graph_branch, GRAPH_NULL, GRAPH_KNOWN
  implicit none
  private
  public :: size_of
contains
  integer function size_of(b) result(n)
    type(graph_branch), intent(in) :: b
    type(graph), pointer :: cell
    n = 0
    if (b % status() .eq. GRAPH_NULL) return
    cell => b % known()
    do
       n = n + 1
       if (cell % branch(2) % status() .eq. GRAPH_NULL) exit
       cell => cell % branch(2) % known()
    end do
  end function size_of
end module branch_form
