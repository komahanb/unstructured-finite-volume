!=====================================================================!
! The reversible change protocol suite: one controller, one
! lifecycle -
!
!      apply -> check -> keep | revert
!
! - proven over five toy changes: a pure attached-value change, a
! pure structure change, a mixed change, a change whose check
! judges against, and a change whose apply fails. Accepted changes
! keep; rejected changes revert; failed steps revert lawfully; and
! the controller never learns which kind it ran. Identity lives in
! the graph. Numbers live in attached values. Updates are
! reversible changes.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_toy_changes

  use gti_change_protocols, only : gti_reversible_change, gti_change_result

  implicit none

  private
  public :: toy_value_change
  public :: toy_structure_change
  public :: toy_mixed_change
  public :: toy_vetoed_change
  public :: toy_broken_apply_change

  !===================================================================!
  ! A: the pure attached-value change - one number, no structure.
  !===================================================================!

  type, extends(gti_reversible_change) :: toy_value_change
     integer :: value = 0
   contains
     procedure :: apply  => value_apply
     procedure :: check  => value_check
     procedure :: keep   => value_keep
     procedure :: revert => value_revert
  end type toy_value_change

  !===================================================================!
  ! B: the pure structure change - one seat count, no numbers.
  !===================================================================!

  type, extends(gti_reversible_change) :: toy_structure_change
     integer :: structure_count = 0
   contains
     procedure :: apply  => structure_apply
     procedure :: check  => structure_check
     procedure :: keep   => structure_keep
     procedure :: revert => structure_revert
  end type toy_structure_change

  !===================================================================!
  ! C: the mixed change - structure and attached values at once.
  !===================================================================!

  type, extends(gti_reversible_change) :: toy_mixed_change
     integer :: structure_count = 0
     integer :: value_count = 0
   contains
     procedure :: apply  => mixed_apply
     procedure :: check  => mixed_check
     procedure :: keep   => mixed_keep
     procedure :: revert => mixed_revert
  end type toy_mixed_change

  !===================================================================!
  ! D: the vetoed change - its apply works, its check judges
  ! against it, and the controller must revert.
  !===================================================================!

  type, extends(gti_reversible_change) :: toy_vetoed_change
     integer :: value = 0
   contains
     procedure :: apply  => vetoed_apply
     procedure :: check  => vetoed_check
     procedure :: keep   => vetoed_keep
     procedure :: revert => vetoed_revert
  end type toy_vetoed_change

  !===================================================================!
  ! E: the broken apply - it marks failed and builds nothing, and
  ! the controller must revert lawfully, never refuse.
  !===================================================================!

  type, extends(gti_reversible_change) :: toy_broken_apply_change
     integer :: value = 0
   contains
     procedure :: apply  => broken_apply
     procedure :: check  => broken_check
     procedure :: keep   => broken_keep
     procedure :: revert => broken_revert
  end type toy_broken_apply_change

contains

  subroutine value_apply(this, result)
    class(toy_value_change), intent(inout) :: this
    type(gti_change_result), intent(inout) :: result
    result % touches_value = .true.
    this % value = 1
    call result % mark_applied()
  end subroutine value_apply

  subroutine value_check(this, result)
    class(toy_value_change), intent(inout) :: this
    type(gti_change_result), intent(inout) :: result
    associate(unread => this)
    end associate
    call result % mark_checked(.true.)
  end subroutine value_check

  subroutine value_keep(this, result)
    class(toy_value_change), intent(inout) :: this
    type(gti_change_result), intent(inout) :: result
    associate(unread => this)
    end associate
    call result % mark_kept()
  end subroutine value_keep

  subroutine value_revert(this, result)
    class(toy_value_change), intent(inout) :: this
    type(gti_change_result), intent(inout) :: result
    this % value = 0
    call result % mark_reverted()
  end subroutine value_revert

  subroutine structure_apply(this, result)
    class(toy_structure_change), intent(inout) :: this
    type(gti_change_result)    , intent(inout) :: result
    result % touches_structure = .true.
    this % structure_count = this % structure_count + 1
    call result % mark_applied()
  end subroutine structure_apply

  subroutine structure_check(this, result)
    class(toy_structure_change), intent(inout) :: this
    type(gti_change_result)    , intent(inout) :: result
    associate(unread => this)
    end associate
    call result % mark_checked(.true.)
  end subroutine structure_check

  subroutine structure_keep(this, result)
    class(toy_structure_change), intent(inout) :: this
    type(gti_change_result)    , intent(inout) :: result
    associate(unread => this)
    end associate
    call result % mark_kept()
  end subroutine structure_keep

  subroutine structure_revert(this, result)
    class(toy_structure_change), intent(inout) :: this
    type(gti_change_result)    , intent(inout) :: result
    this % structure_count = this % structure_count - 1
    call result % mark_reverted()
  end subroutine structure_revert

  subroutine mixed_apply(this, result)
    class(toy_mixed_change), intent(inout) :: this
    type(gti_change_result), intent(inout) :: result
    result % touches_structure = .true.
    result % touches_value     = .true.
    this % structure_count = this % structure_count + 1
    this % value_count     = this % value_count + 1
    call result % mark_applied()
  end subroutine mixed_apply

  subroutine mixed_check(this, result)
    class(toy_mixed_change), intent(inout) :: this
    type(gti_change_result), intent(inout) :: result
    associate(unread => this)
    end associate
    call result % mark_checked(.true.)
  end subroutine mixed_check

  subroutine mixed_keep(this, result)
    class(toy_mixed_change), intent(inout) :: this
    type(gti_change_result), intent(inout) :: result
    associate(unread => this)
    end associate
    call result % mark_kept()
  end subroutine mixed_keep

  subroutine mixed_revert(this, result)
    class(toy_mixed_change), intent(inout) :: this
    type(gti_change_result), intent(inout) :: result
    this % structure_count = this % structure_count - 1
    this % value_count     = this % value_count - 1
    call result % mark_reverted()
  end subroutine mixed_revert

  subroutine vetoed_apply(this, result)
    class(toy_vetoed_change), intent(inout) :: this
    type(gti_change_result) , intent(inout) :: result
    result % touches_value = .true.
    this % value = 1
    call result % mark_applied()
  end subroutine vetoed_apply

  subroutine vetoed_check(this, result)
    class(toy_vetoed_change), intent(inout) :: this
    type(gti_change_result) , intent(inout) :: result
    associate(unread => this)
    end associate
    call result % mark_checked(.false.)
  end subroutine vetoed_check

  subroutine vetoed_keep(this, result)
    class(toy_vetoed_change), intent(inout) :: this
    type(gti_change_result) , intent(inout) :: result
    associate(unread => this)
    end associate
    call result % mark_kept()
  end subroutine vetoed_keep

  subroutine vetoed_revert(this, result)
    class(toy_vetoed_change), intent(inout) :: this
    type(gti_change_result) , intent(inout) :: result
    this % value = 0
    call result % mark_reverted()
  end subroutine vetoed_revert

  subroutine broken_apply(this, result)
    class(toy_broken_apply_change), intent(inout) :: this
    type(gti_change_result)       , intent(inout) :: result
    associate(unread => this)
    end associate
    call result % mark_failed()
  end subroutine broken_apply

  subroutine broken_check(this, result)
    class(toy_broken_apply_change), intent(inout) :: this
    type(gti_change_result)       , intent(inout) :: result
    associate(unread => this)
    end associate
    call result % mark_checked(.true.)
  end subroutine broken_check

  subroutine broken_keep(this, result)
    class(toy_broken_apply_change), intent(inout) :: this
    type(gti_change_result)       , intent(inout) :: result
    associate(unread => this)
    end associate
    call result % mark_kept()
  end subroutine broken_keep

  subroutine broken_revert(this, result)
    class(toy_broken_apply_change), intent(inout) :: this
    type(gti_change_result)       , intent(inout) :: result
    associate(unread => this)
    end associate
    call result % mark_reverted()
  end subroutine broken_revert

end module gti_toy_changes

program test_gti_change_protocol

  use gti_change_protocols, only : gti_change_controller, gti_change_result
  use gti_toy_changes     , only : toy_value_change, toy_structure_change, &
       & toy_mixed_change, toy_vetoed_change, toy_broken_apply_change

  implicit none

  type(gti_change_controller) :: controller
  type(gti_change_result)     :: result

  type(toy_value_change)        :: value_change
  type(toy_structure_change)    :: structure_change
  type(toy_mixed_change)        :: mixed_change
  type(toy_vetoed_change)       :: vetoed_change
  type(toy_broken_apply_change) :: broken_change

  integer :: nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "gti reversible change protocol suite"
  write(*,'(1x,a)') "============================================="

  !-------------------------------------------------------------------!
  ! A: the pure attached-value change, accepted then rejected.
  !-------------------------------------------------------------------!

  call controller % run(value_change, .true., result)
  call report(value_change % value == 1, &
       & "an accepted pure value change keeps the value", nfail)
  call report(result % kept .and. .not. result % reverted, &
       & "the accepted result reads kept and not reverted", nfail)
  call report(result % accepted .and. result % attempted .and. &
       & result % applied .and. result % checked .and. result % check_passed, &
       & "the accepted lifecycle is fully recorded", nfail)
  call report(result % touches_value .and. .not. result % touches_structure, &
       & "the value change declares itself: values, not structure", nfail)

  value_change % value = 0
  call controller % run(value_change, .false., result)
  call report(value_change % value == 0, &
       & "a rejected pure value change reverts the value", nfail)
  call report(.not. result % kept .and. result % reverted .and. &
       & .not. result % accepted, &
       & "the rejected result reads reverted and not kept", nfail)

  !-------------------------------------------------------------------!
  ! B: the pure structure change, accepted then rejected.
  !-------------------------------------------------------------------!

  call controller % run(structure_change, .true., result)
  call report(structure_change % structure_count == 1, &
       & "an accepted pure structure change keeps the structure", nfail)
  call report(result % touches_structure .and. .not. result % touches_value, &
       & "the structure change declares itself: structure, not values", nfail)

  call controller % run(structure_change, .false., result)
  call report(structure_change % structure_count == 1, &
       & "a rejected pure structure change reverts the structure", nfail)

  !-------------------------------------------------------------------!
  ! C: the mixed change - the adaptive-growth shape - accepted
  ! then rejected.
  !-------------------------------------------------------------------!

  call controller % run(mixed_change, .true., result)
  call report(mixed_change % structure_count == 1 .and. &
       & mixed_change % value_count == 1, &
       & "an accepted mixed change keeps structure and values", nfail)
  call report(result % touches_structure .and. result % touches_value, &
       & "the mixed change declares itself: both at once", nfail)

  call controller % run(mixed_change, .false., result)
  call report(mixed_change % structure_count == 1 .and. &
       & mixed_change % value_count == 1, &
       & "a rejected mixed change reverts structure and values", nfail)

  !-------------------------------------------------------------------!
  ! D: the vetoed change - check judges against it, and even an
  ! accepting caller cannot keep it.
  !-------------------------------------------------------------------!

  call controller % run(vetoed_change, .true., result)
  call report(vetoed_change % value == 0, &
       & "a change whose check judges against it is reverted", nfail)
  call report(result % attempted .and. result % applied .and. &
       & result % checked .and. .not. result % check_passed .and. &
       & .not. result % accepted .and. .not. result % kept .and. &
       & result % reverted .and. .not. result % failed, &
       & "the vetoed lifecycle reads exactly: checked, not passed, reverted", nfail)

  !-------------------------------------------------------------------!
  ! E: the broken apply - failure is a lawful answer, reverted and
  ! reported, never refused.
  !-------------------------------------------------------------------!

  call controller % run(broken_change, .true., result)
  call report(result % attempted .and. .not. result % applied .and. &
       & result % failed .and. result % reverted .and. .not. result % kept, &
       & "a failed apply reverts and reports: attempted, failed, reverted", nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all reversible change checks passed"
  else
     error stop
  end if

contains

  subroutine report(ok, label, nfail)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: label
    integer, intent(inout) :: nfail
    if (ok) then
       write(*,'(1x,a,a)') "PASS : ", label
    else
       write(*,'(1x,a,a)') "FAIL : ", label
       nfail = nfail + 1
    end if
  end subroutine report

end program test_gti_change_protocol
