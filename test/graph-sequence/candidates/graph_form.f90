! CANDIDATE 2. The sequence is a graph: the first cell.
!
! There is no cell for the empty sequence, so the caller must decide
! emptiness before calling, and every consumer repeats that decision.
module graph_form
  use graph_fractal, only : graph, branch, BRANCH_NULL, BRANCH_KNOWN
  implicit none
  private
  public :: num_members, is_empty
contains
  logical function is_empty(b) result(e)
    type(branch), intent(in) :: b       ! still needs the branch
    e = b % status() .eq. BRANCH_NULL
  end function is_empty
  integer function num_members(first) result(n)
    type(graph), target, intent(in) :: first
    type(graph), pointer :: cell
    n = 1
    cell => first
    do
       if (cell % branch(2) % status() .eq. BRANCH_NULL) exit
       cell => cell % branch(2) % known()
       n = n + 1
    end do
  end function num_members
end module graph_form
