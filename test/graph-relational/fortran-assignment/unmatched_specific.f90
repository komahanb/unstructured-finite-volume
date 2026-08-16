!=====================================================================!
! A public generic assignment with no matching specific does not
! prohibit assignment either.
!
! There is no way to declare an operation deleted: where no specific
! matches, intrinsic assignment substitutes, again with no diagnostic.
!=====================================================================!

module unmatched_specific_m

  implicit none

  private
  public :: bind_t, defined_assignment_ran

  logical :: defined_assignment_ran = .false.

  type :: other_t
     integer :: k = 0
  end type other_t

  type :: bind_t
     integer, pointer, private :: object => null()
   contains
     procedure, private :: copy_other
     generic :: assignment(=) => copy_other
  end type bind_t

contains

  subroutine copy_other(lhs, rhs)

    class(bind_t), intent(out) :: lhs
    type(other_t), intent(in)  :: rhs

    defined_assignment_ran = .true.
    lhs % object => null()

  end subroutine copy_other

end module unmatched_specific_m

program unmatched_specific

  use unmatched_specific_m

  implicit none

  type(bind_t) :: a, b

  a = b

  if (defined_assignment_ran) then
     print *, '   the declared specific was used'
  else
     print *, '   intrinsic assignment substituted, with no diagnostic'
  end if

end program unmatched_specific
