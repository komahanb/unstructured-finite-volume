!=====================================================================!
! Test fixtures shared by test.f90 and refusal.f90: concrete
! extensions of reversible_change.
!
!      counting_change   increments an integer counter on apply
!                        and decrements it on revert. check_passes
!                        sets the check verdict; fail_apply makes
!                        apply mark failure and skip the mutation.
!                        One type covers the accept, reject, veto,
!                        and failed-apply paths.
!      mixed_change      mutates a counter and a real value in one
!                        apply and restores both on revert;
!                        reports touches_structure and
!                        touches_value.
!      silent_apply_change        returns from apply without marking
!                        applied or failed; the controller must
!                        error stop on it.
!      silent_revert_change       marks failure in apply, then returns from
!                        revert without marking reverted; the
!                        controller must error stop on it.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module toy_changes

  use iso_fortran_env      , only : dp => REAL64
  use graph_change_protocol, only : reversible_change, change_result

  implicit none

  private
  public :: counting_change, mixed_change, silent_apply_change, silent_revert_change

  !===================================================================!
  ! The structure-only change: its state is one integer counter.
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
  ! The mixed change: one apply mutates a counter and a real
  ! value; one revert restores both.
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
  ! Changes that omit a required lifecycle mark; used by the
  ! refusal cases.
  !===================================================================!

  type, extends(reversible_change) :: silent_apply_change
   contains
     procedure :: apply  => silent_apply_apply
     procedure :: check  => silent_apply_check
     procedure :: keep   => silent_apply_keep
     procedure :: revert => silent_apply_revert
  end type silent_apply_change

  type, extends(reversible_change) :: silent_revert_change
   contains
     procedure :: apply  => silent_revert_apply
     procedure :: check  => silent_revert_check
     procedure :: keep   => silent_revert_keep
     procedure :: revert => silent_revert_revert
  end type silent_revert_change

contains

  !===================================================================!
  ! counting_change. When fail_apply is set, apply marks failure
  ! and returns without mutating; otherwise it increments rooms
  ! and records the mutation in grown, so that revert decrements
  ! only when a mutation actually happened.
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
  ! mixed_change. apply increments rooms and doubles value; revert
  ! undoes both, guarded by grown as above.
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
  ! silent_apply_change: apply returns having marked neither applied nor
  ! failed; the controller must refuse this.
  !===================================================================!

  subroutine silent_apply_apply(this, result)
    class(silent_apply_change)  , intent(inout) :: this
    type(change_result), intent(inout) :: result
    associate(u1 => this, u2 => result); end associate
  end subroutine silent_apply_apply

  subroutine silent_apply_check(this, result)
    class(silent_apply_change)  , intent(inout) :: this
    type(change_result), intent(inout) :: result
    associate(u1 => this); end associate
    call result % mark_checked(.true.)
  end subroutine silent_apply_check

  subroutine silent_apply_keep(this, result)
    class(silent_apply_change)  , intent(inout) :: this
    type(change_result), intent(inout) :: result
    associate(u1 => this); end associate
    call result % mark_kept()
  end subroutine silent_apply_keep

  subroutine silent_apply_revert(this, result)
    class(silent_apply_change)  , intent(inout) :: this
    type(change_result), intent(inout) :: result
    associate(u1 => this); end associate
    call result % mark_reverted()
  end subroutine silent_apply_revert

  !===================================================================!
  ! silent_revert_change: apply marks failure; revert returns without
  ! marking reverted; the controller must refuse this.
  !===================================================================!

  subroutine silent_revert_apply(this, result)
    class(silent_revert_change) , intent(inout) :: this
    type(change_result), intent(inout) :: result
    associate(u1 => this); end associate
    call result % mark_failed()
  end subroutine silent_revert_apply

  subroutine silent_revert_check(this, result)
    class(silent_revert_change) , intent(inout) :: this
    type(change_result), intent(inout) :: result
    associate(u1 => this); end associate
    call result % mark_checked(.true.)
  end subroutine silent_revert_check

  subroutine silent_revert_keep(this, result)
    class(silent_revert_change) , intent(inout) :: this
    type(change_result), intent(inout) :: result
    associate(u1 => this); end associate
    call result % mark_kept()
  end subroutine silent_revert_keep

  subroutine silent_revert_revert(this, result)
    class(silent_revert_change) , intent(inout) :: this
    type(change_result), intent(inout) :: result
    associate(u1 => this, u2 => result); end associate
  end subroutine silent_revert_revert

end module toy_changes
