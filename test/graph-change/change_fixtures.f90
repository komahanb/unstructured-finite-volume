!=====================================================================!
! Test fixtures shared by test.f90 and refusal.f90: concrete
! extensions of reversible_change.
!
!      counting_change   increments an integer counter on apply
!                        and decrements it on revert. check_passes
!                        sets the check result; fail_apply makes
!                        apply mark failure and skip the mutation.
!                        One type covers the accept, reject, veto,
!                        and failed-apply paths.
!      mixed_change      mutates a counter and a real value in one
!                        apply and restores both on revert;
!                        reports touches_structure and
!                        touches_value.
!      unreported_apply_change        returns from apply without marking
!                        applied or failed; run_change must
!                        error stop on it.
!      unreported_revert_change       marks failure in apply, then returns from
!                        revert without marking reverted; the
!                        run_change must error stop on it.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module change_fixtures

  use iso_fortran_env      , only : dp => REAL64
  use map_change_protocol, only : reversible_change, change_record

  implicit none

  private
  public :: counting_change, mixed_change, unreported_apply_change, unreported_revert_change

  !===================================================================!
  ! The structure-only change: its state is one integer counter.
  !===================================================================!

  type, extends(reversible_change) :: counting_change

     integer :: counter = 0

     logical :: check_passes = .true.
     logical :: fail_apply   = .false.

     logical, private :: increment_applied = .false.

   contains

     procedure :: apply  => counting_apply
     procedure :: check  => counting_check
     procedure :: commit => counting_commit
     procedure :: revert => counting_revert

  end type counting_change

  !===================================================================!
  ! The mixed change: one apply mutates a counter and a real
  ! value; one revert restores both.
  !===================================================================!

  type, extends(reversible_change) :: mixed_change

     integer  :: counter = 0
     real(dp) :: value = 1.0_dp

     logical :: check_passes = .true.

     logical, private :: increment_applied = .false.

   contains

     procedure :: apply  => mixed_apply
     procedure :: check  => mixed_check
     procedure :: commit => mixed_commit
     procedure :: revert => mixed_revert

  end type mixed_change

  !===================================================================!
  ! Changes that omit a required lifecycle mark; used by the
  ! refusal cases.
  !===================================================================!

  type, extends(reversible_change) :: unreported_apply_change
   contains
     procedure :: apply  => unreported_apply_apply
     procedure :: check  => unreported_apply_check
     procedure :: commit => unreported_apply_commit
     procedure :: revert => unreported_apply_revert
  end type unreported_apply_change

  type, extends(reversible_change) :: unreported_revert_change
   contains
     procedure :: apply  => unreported_revert_apply
     procedure :: check  => unreported_revert_check
     procedure :: commit => unreported_revert_commit
     procedure :: revert => unreported_revert_revert
  end type unreported_revert_change

contains

  !===================================================================!
  ! counting_change. When fail_apply is set, apply marks failure
  ! and returns without mutating; otherwise it increments counter
  ! and records the mutation in increment_applied, so that revert decrements
  ! only when a mutation occurred.
  !===================================================================!

  subroutine counting_apply(this, result)

    class(counting_change), intent(inout) :: this
    type(change_record)   , intent(inout) :: result

    result % touches_structure = .true.

    if (this % fail_apply) then
       call result % mark_failed()
       return
    end if

    this % counter = this % counter + 1
    this % increment_applied = .true.
    call result % mark_applied()

  end subroutine counting_apply

  subroutine counting_check(this, result)

    class(counting_change), intent(inout) :: this
    type(change_record)   , intent(inout) :: result

    call result % mark_checked(this % check_passes)

  end subroutine counting_check

  subroutine counting_commit(this, result)

    class(counting_change), intent(inout) :: this
    type(change_record)   , intent(inout) :: result

    this % increment_applied = .false.
    call result % mark_committed()

  end subroutine counting_commit

  subroutine counting_revert(this, result)

    class(counting_change), intent(inout) :: this
    type(change_record)   , intent(inout) :: result

    if (this % increment_applied) then
       this % counter = this % counter - 1
       this % increment_applied = .false.
    end if
    call result % mark_reverted()

  end subroutine counting_revert

  !===================================================================!
  ! mixed_change. apply increments counter and doubles value; revert
  ! undoes both, conditioned on increment_applied as above.
  !===================================================================!

  subroutine mixed_apply(this, result)

    class(mixed_change), intent(inout) :: this
    type(change_record), intent(inout) :: result

    result % touches_structure = .true.
    result % touches_value     = .true.

    this % counter = this % counter + 1
    this % value = 2.0_dp * this % value
    this % increment_applied = .true.
    call result % mark_applied()

  end subroutine mixed_apply

  subroutine mixed_check(this, result)

    class(mixed_change), intent(inout) :: this
    type(change_record), intent(inout) :: result

    call result % mark_checked(this % check_passes)

  end subroutine mixed_check

  subroutine mixed_commit(this, result)

    class(mixed_change), intent(inout) :: this
    type(change_record), intent(inout) :: result

    this % increment_applied = .false.
    call result % mark_committed()

  end subroutine mixed_commit

  subroutine mixed_revert(this, result)

    class(mixed_change), intent(inout) :: this
    type(change_record), intent(inout) :: result

    if (this % increment_applied) then
       this % counter = this % counter - 1
       this % value = 0.5_dp * this % value
       this % increment_applied = .false.
    end if
    call result % mark_reverted()

  end subroutine mixed_revert

  !===================================================================!
  ! unreported_apply_change: apply returns having marked neither applied nor
  ! failed; run_change must refuse this.
  !===================================================================!

  subroutine unreported_apply_apply(this, result)
    class(unreported_apply_change)  , intent(inout) :: this
    type(change_record), intent(inout) :: result
    associate(u1 => this, u2 => result); end associate
  end subroutine unreported_apply_apply

  subroutine unreported_apply_check(this, result)
    class(unreported_apply_change)  , intent(inout) :: this
    type(change_record), intent(inout) :: result
    associate(u1 => this); end associate
    call result % mark_checked(.true.)
  end subroutine unreported_apply_check

  subroutine unreported_apply_commit(this, result)
    class(unreported_apply_change)  , intent(inout) :: this
    type(change_record), intent(inout) :: result
    associate(u1 => this); end associate
    call result % mark_committed()
  end subroutine unreported_apply_commit

  subroutine unreported_apply_revert(this, result)
    class(unreported_apply_change)  , intent(inout) :: this
    type(change_record), intent(inout) :: result
    associate(u1 => this); end associate
    call result % mark_reverted()
  end subroutine unreported_apply_revert

  !===================================================================!
  ! unreported_revert_change: apply marks failure; revert returns without
  ! marking reverted; run_change must refuse this.
  !===================================================================!

  subroutine unreported_revert_apply(this, result)
    class(unreported_revert_change) , intent(inout) :: this
    type(change_record), intent(inout) :: result
    associate(u1 => this); end associate
    call result % mark_failed()
  end subroutine unreported_revert_apply

  subroutine unreported_revert_check(this, result)
    class(unreported_revert_change) , intent(inout) :: this
    type(change_record), intent(inout) :: result
    associate(u1 => this); end associate
    call result % mark_checked(.true.)
  end subroutine unreported_revert_check

  subroutine unreported_revert_commit(this, result)
    class(unreported_revert_change) , intent(inout) :: this
    type(change_record), intent(inout) :: result
    associate(u1 => this); end associate
    call result % mark_committed()
  end subroutine unreported_revert_commit

  subroutine unreported_revert_revert(this, result)
    class(unreported_revert_change) , intent(inout) :: this
    type(change_record), intent(inout) :: result
    associate(u1 => this, u2 => result); end associate
  end subroutine unreported_revert_revert

end module change_fixtures
