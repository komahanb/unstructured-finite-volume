!=====================================================================!
! A private generic assignment binding does not prohibit assignment.
!
! Outside the module the binding is inaccessible, so intrinsic
! assignment substitutes for it - silently, with no diagnostic. For a
! type whose rows hold pointers that is a shallow copy and a double
! free: worse than no mechanism at all.
!=====================================================================!

module private_generic_m

  implicit none

  private
  public :: bind_t, defined_assignment_ran

  logical :: defined_assignment_ran = .false.

  type :: bind_t
     integer, pointer, private :: object => null()
   contains
     procedure, private :: copy
     generic, private :: assignment(=) => copy
  end type bind_t

contains

  subroutine copy(lhs, rhs)

    class(bind_t), intent(out) :: lhs
    type(bind_t) , intent(in)  :: rhs

    defined_assignment_ran = .true.
    if (associated(rhs % object)) allocate(lhs % object, source=rhs % object)

  end subroutine copy

end module private_generic_m

program private_generic

  use private_generic_m

  implicit none

  type(bind_t) :: a, b

  a = b

  if (defined_assignment_ran) then
     print *, '   the private binding was used'
  else
     print *, '   intrinsic assignment substituted, with no diagnostic'
  end if

end program private_generic
