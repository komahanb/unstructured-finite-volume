!=====================================================================!
! The refusals that must die at the change-protocol seat, one per
! invocation:
!
!      noapply      an apply that neither marks applied nor fails
!      nocheck      a check that neither marks checked nor fails
!      nokeep       an accepted keep that does not mark kept
!      norevert     a rejected revert that does not mark reverted
!      badterminal  a terminal record kept AND reverted at once
!
! A change may lawfully fail; it may not lie about its own
! lifecycle. Every case must error stop; a case that returns is a
! failure of the suite.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_liar_changes

  use gti_change_protocols, only : gti_reversible_change, gti_change_result

  implicit none

  private
  public :: liar_change

  !===================================================================!
  ! One configurable liar: each verb does its marking only if its
  ! knob says so - and never marks failure, which is the lie.
  !===================================================================!

  type, extends(gti_reversible_change) :: liar_change
     logical :: says_applied  = .true.
     logical :: says_checked  = .true.
     logical :: says_kept     = .true.
     logical :: says_reverted = .true.
   contains
     procedure :: apply  => liar_apply
     procedure :: check  => liar_check
     procedure :: keep   => liar_keep
     procedure :: revert => liar_revert
  end type liar_change

contains

  subroutine liar_apply(this, result)
    class(liar_change)     , intent(inout) :: this
    type(gti_change_result), intent(inout) :: result
    if (this % says_applied) call result % mark_applied()
  end subroutine liar_apply

  subroutine liar_check(this, result)
    class(liar_change)     , intent(inout) :: this
    type(gti_change_result), intent(inout) :: result
    if (this % says_checked) call result % mark_checked(.true.)
  end subroutine liar_check

  subroutine liar_keep(this, result)
    class(liar_change)     , intent(inout) :: this
    type(gti_change_result), intent(inout) :: result
    if (this % says_kept) call result % mark_kept()
  end subroutine liar_keep

  subroutine liar_revert(this, result)
    class(liar_change)     , intent(inout) :: this
    type(gti_change_result), intent(inout) :: result
    if (this % says_reverted) call result % mark_reverted()
  end subroutine liar_revert

end module gti_liar_changes

program refusal

  use gti_change_protocols, only : gti_change_controller, gti_change_result
  use gti_liar_changes    , only : liar_change

  implicit none

  type(gti_change_controller) :: controller
  type(gti_change_result)     :: result
  type(liar_change)           :: liar

  character(len=32) :: which

  call get_command_argument(1, which)

  select case (trim(which))

  case ('noapply')

     liar % says_applied = .false.
     call controller % run(liar, .true., result)

  case ('nocheck')

     liar % says_checked = .false.
     call controller % run(liar, .true., result)

  case ('nokeep')

     liar % says_kept = .false.
     call controller % run(liar, .true., result)

  case ('norevert')

     liar % says_reverted = .false.
     call controller % run(liar, .false., result)

  case ('badterminal')

     result % kept     = .true.
     result % reverted = .true.
     call result % validate_terminal()

  case default

     error stop 'refusal: unknown case'

  end select

  write(*,*) 'refusal case survived: ', trim(which)

end program refusal
