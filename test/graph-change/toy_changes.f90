!=====================================================================!
! THE TOY CHANGES: concrete members of the reversible change
! family, each owning exactly what the protocol leaves to the
! change - the mutation, its memory, and its undoing.
!
!      counting_change   touches structure alone: one room grown
!                        on apply, ungrown on revert; a veto flag
!                        and a failing-apply flag make one type
!                        serve the accept, reject, veto, and
!                        failure walks
!      mixed_change      touches structure AND value in one apply,
!                        and revert restores both - the change the
!                        controller must not care about
!      liar_apply        does its work and never marks it
!      liar_revert       fails, then reverts without saying so
!
! The liars exist to die: a change that lies about its own
! lifecycle is refused by the controller, and the refusal suite
! proves it.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module toy_changes

  use iso_fortran_env      , only : dp => REAL64
  use graph_change_protocol, only : reversible_change, change_result

  implicit none

  private
  public :: counting_change, mixed_change, liar_apply, liar_revert

  !===================================================================!
  ! The structural specimen: the state is one room count.
  !===================================================================!

  type, extends(reversible_change) :: counting_change

     integer :: rooms = 0

     logical :: check_passes = .true.
     logical :: fail_apply   = .false.

     logical, private :: grown = .false.

   contains

     procedure :: apply  => counting_apply
     procedure :: check  => counting_check
     procedure :: keep   => counting_keep
     procedure :: revert => counting_revert

  end type counting_change

  !===================================================================!
  ! The mixed specimen: one apply moves a count and a number, one
  ! revert restores both.
  !===================================================================!

  type, extends(reversible_change) :: mixed_change

     integer  :: rooms = 0
     real(dp) :: value = 1.0_dp

     logical :: check_passes = .true.

     logical, private :: grown = .false.

   contains

     procedure :: apply  => mixed_apply
     procedure :: check  => mixed_check
     procedure :: keep   => mixed_keep
     procedure :: revert => mixed_revert

  end type mixed_change

  !===================================================================!
  ! The liars.
  !===================================================================!

  type, extends(reversible_change) :: liar_apply
   contains
     procedure :: apply  => liar_apply_apply
     procedure :: check  => liar_apply_check
     procedure :: keep   => liar_apply_keep
     procedure :: revert => liar_apply_revert
  end type liar_apply

  type, extends(reversible_change) :: liar_revert
   contains
     procedure :: apply  => liar_revert_apply
     procedure :: check  => liar_revert_check
     procedure :: keep   => liar_revert_keep
     procedure :: revert => liar_revert_revert
  end type liar_revert

contains

  !===================================================================!
  ! The counting change: grow one room, or fail before touching
  ! anything - so revert has exactly as much to undo as apply did.
  !===================================================================!

  subroutine counting_apply(this, result)

    class(counting_change), intent(inout) :: this
    type(change_result)   , intent(inout) :: result

    result % touches_structure = .true.

    if (this % fail_apply) then
       call result % mark_failed()
       return
    end if

    this % rooms = this % rooms + 1
    this % grown = .true.
    call result % mark_applied()

  end subroutine counting_apply

  subroutine counting_check(this, result)

    class(counting_change), intent(inout) :: this
    type(change_result)   , intent(inout) :: result

    call result % mark_checked(this % check_passes)

  end subroutine counting_check

  subroutine counting_keep(this, result)

    class(counting_change), intent(inout) :: this
    type(change_result)   , intent(inout) :: result

    this % grown = .false.
    call result % mark_kept()

  end subroutine counting_keep

  subroutine counting_revert(this, result)

    class(counting_change), intent(inout) :: this
    type(change_result)   , intent(inout) :: result

    if (this % grown) then
       this % rooms = this % rooms - 1
       this % grown = .false.
    end if
    call result % mark_reverted()

  end subroutine counting_revert

  !===================================================================!
  ! The mixed change: structure and value move together, and are
  ! restored together.
  !===================================================================!

  subroutine mixed_apply(this, result)

    class(mixed_change), intent(inout) :: this
    type(change_result), intent(inout) :: result

    result % touches_structure = .true.
    result % touches_value     = .true.

    this % rooms = this % rooms + 1
    this % value = 2.0_dp * this % value
    this % grown = .true.
    call result % mark_applied()

  end subroutine mixed_apply

  subroutine mixed_check(this, result)

    class(mixed_change), intent(inout) :: this
    type(change_result), intent(inout) :: result

    call result % mark_checked(this % check_passes)

  end subroutine mixed_check

  subroutine mixed_keep(this, result)

    class(mixed_change), intent(inout) :: this
    type(change_result), intent(inout) :: result

    this % grown = .false.
    call result % mark_kept()

  end subroutine mixed_keep

  subroutine mixed_revert(this, result)

    class(mixed_change), intent(inout) :: this
    type(change_result), intent(inout) :: result

    if (this % grown) then
       this % rooms = this % rooms - 1
       this % value = 0.5_dp * this % value
       this % grown = .false.
    end if
    call result % mark_reverted()

  end subroutine mixed_revert

  !===================================================================!
  ! The apply liar: returns from apply having marked neither
  ! applied nor failed.
  !===================================================================!

  subroutine liar_apply_apply(this, result)
    class(liar_apply)  , intent(inout) :: this
    type(change_result), intent(inout) :: result
    associate(u1 => this, u2 => result); end associate
  end subroutine liar_apply_apply

  subroutine liar_apply_check(this, result)
    class(liar_apply)  , intent(inout) :: this
    type(change_result), intent(inout) :: result
    associate(u1 => this); end associate
    call result % mark_checked(.true.)
  end subroutine liar_apply_check

  subroutine liar_apply_keep(this, result)
    class(liar_apply)  , intent(inout) :: this
    type(change_result), intent(inout) :: result
    associate(u1 => this); end associate
    call result % mark_kept()
  end subroutine liar_apply_keep

  subroutine liar_apply_revert(this, result)
    class(liar_apply)  , intent(inout) :: this
    type(change_result), intent(inout) :: result
    associate(u1 => this); end associate
    call result % mark_reverted()
  end subroutine liar_apply_revert

  !===================================================================!
  ! The revert liar: fails at apply, then returns from revert
  ! without marking it.
  !===================================================================!

  subroutine liar_revert_apply(this, result)
    class(liar_revert) , intent(inout) :: this
    type(change_result), intent(inout) :: result
    associate(u1 => this); end associate
    call result % mark_failed()
  end subroutine liar_revert_apply

  subroutine liar_revert_check(this, result)
    class(liar_revert) , intent(inout) :: this
    type(change_result), intent(inout) :: result
    associate(u1 => this); end associate
    call result % mark_checked(.true.)
  end subroutine liar_revert_check

  subroutine liar_revert_keep(this, result)
    class(liar_revert) , intent(inout) :: this
    type(change_result), intent(inout) :: result
    associate(u1 => this); end associate
    call result % mark_kept()
  end subroutine liar_revert_keep

  subroutine liar_revert_revert(this, result)
    class(liar_revert) , intent(inout) :: this
    type(change_result), intent(inout) :: result
    associate(u1 => this, u2 => result); end associate
  end subroutine liar_revert_revert

end module toy_changes
