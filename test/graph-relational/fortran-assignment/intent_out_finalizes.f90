!=====================================================================!
! Why the refusing assignment takes INTENT(INOUT).
!
! An INTENT(OUT) dummy of a finalizable type is finalized on entry, so
! a refusal declared that way destroys the lender before it refuses to
! destroy it. INTENT(INOUT) is not finalized, and the left-hand side is
! still intact when the body runs.
!=====================================================================!

module intent_out_finalizes_m

  implicit none

  private
  public :: out_t, inout_t

  type :: out_t
     integer :: k = 0
   contains
     procedure, private :: assign_out
     generic :: assignment(=) => assign_out
     final :: release_out
  end type out_t

  type :: inout_t
     integer :: k = 0
   contains
     procedure, private :: assign_inout
     generic :: assignment(=) => assign_inout
     final :: release_inout
  end type inout_t

  logical, public :: finalized = .false.

contains

  subroutine assign_out(lhs, rhs)
    class(out_t), intent(out) :: lhs
    type(out_t) , intent(in)  :: rhs
    print *, '   INTENT(OUT)   body: the left-hand side was finalized first =', finalized
  end subroutine assign_out

  subroutine release_out(this)
    type(out_t), intent(inout) :: this
    finalized = .true.
  end subroutine release_out

  subroutine assign_inout(lhs, rhs)
    class(inout_t), intent(inout) :: lhs
    type(inout_t) , intent(in)    :: rhs
    print *, '   INTENT(INOUT) body: the left-hand side survives, k =', lhs % k
  end subroutine assign_inout

  subroutine release_inout(this)
    type(inout_t), intent(inout) :: this
  end subroutine release_inout

end module intent_out_finalizes_m

program intent_out_finalizes

  use intent_out_finalizes_m

  implicit none

  block
    type(out_t) :: a, b
    a % k = 7
    a = b
  end block

  block
    type(inout_t) :: a, b
    a % k = 7
    a = b
  end block

end program intent_out_finalizes
