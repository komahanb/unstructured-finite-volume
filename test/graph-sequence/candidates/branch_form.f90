! CANDIDATE 1. The sequence is a branch.
!
!     NULL          the empty sequence
!     KNOWN -> cell a nonempty sequence
!     UNKNOWN       the sequence is not known
!
! The empty sequence needs no graph, because it has no cell.
module branch_form
  use graph_fractal, only : graph, branch, BRANCH_NULL, BRANCH_KNOWN
  implicit none
  private
  public :: num_members
contains
  integer function num_members(b) result(n)
    type(branch), intent(in) :: b
    type(graph), pointer :: cell
    n = 0
    if (b % status() .eq. BRANCH_NULL) return
    cell => b % known()
    do
       n = n + 1
       if (cell % branch(2) % status() .eq. BRANCH_NULL) exit
       cell => cell % branch(2) % known()
    end do
  end function num_members
end module branch_form
