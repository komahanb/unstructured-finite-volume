! CANDIDATES C and D. A private branch component reached through an
! accessor. C returns the branch by value, D by pointer.
!
! This module compiles; what the two accessors admit is established by
! value_accessor_assignment.f90, pointer_accessor_assignment.f90 and
! accessor_navigation.f90.
module accessor_candidates
  implicit none
  public

  type :: branch_t
     integer :: s = 0
   contains
     procedure :: status
  end type branch_t

  type :: graph_t
     private
     type(branch_t) :: branch_(2)
   contains
     procedure :: by_value          ! candidate C
     procedure :: by_pointer        ! candidate D
  end type graph_t

contains

  pure integer function status(this)
    class(branch_t), intent(in) :: this
    status = this % s
  end function status

  type(branch_t) function by_value(this, i) result(b)
    class(graph_t), intent(in) :: this
    integer       , intent(in) :: i
    b = this % branch_(i)
  end function by_value

  function by_pointer(this, i) result(b)
    class(graph_t), target, intent(in) :: this
    integer               , intent(in) :: i
    type(branch_t), pointer            :: b
    b => this % branch_(i)
  end function by_pointer

end module accessor_candidates
