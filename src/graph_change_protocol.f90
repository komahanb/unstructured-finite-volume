!=====================================================================!
! THE REVERSIBLE CHANGE PROTOCOL
!
! One lifecycle for every reversible mutation of a graph and its
! attached values:
!
!      apply -> check -> keep | revert
!
! The controller owns only the lifecycle; the change object owns
! the mutation, its memory, and its undoing. A change may touch
! structure, attached values, or both; the controller never reads
! which.
!
! A failure reported by apply or check is reverted and returned,
! not refused. What stops the program: a change that returns from
! a step without marking either the step or failure, because its
! record would then misreport what happened; and a terminal record
! with contradictory flags. This module imports nothing; concrete
! changes live next to the data they mutate.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module graph_change_protocol

  implicit none

  private
  public :: change_result
  public :: reversible_change
  public :: change_controller

  !===================================================================!
  ! The record of one change's lifecycle: what was attempted, what
  ! each step reported, and what the change declares it touches.
  ! The touch flags are the change's self-description; the
  ! controller never reads them.
  !===================================================================!

  type :: change_result

     logical :: attempted    = .false.
     logical :: applied      = .false.
     logical :: checked      = .false.
     logical :: check_passed = .false.
     logical :: accepted     = .false.
     logical :: kept         = .false.
     logical :: reverted     = .false.
     logical :: failed       = .false.

     logical :: touches_structure = .false.
     logical :: touches_value     = .false.

   contains

     procedure :: reset
     procedure :: mark_attempted
     procedure :: mark_applied
     procedure :: mark_checked
     procedure :: mark_kept
     procedure :: mark_reverted
     procedure :: mark_failed
     procedure :: validate_terminal

  end type change_result

  !===================================================================!
  ! The abstract reversible change: four deferred steps. What
  ! apply builds, revert must undo; keep makes the applied state
  ! permanent.
  !===================================================================!

  type, abstract :: reversible_change

   contains

     procedure(change_apply) , deferred :: apply
     procedure(change_check) , deferred :: check
     procedure(change_keep)  , deferred :: keep
     procedure(change_revert), deferred :: revert

  end type reversible_change

  abstract interface

     subroutine change_apply(this, result)
       import :: reversible_change, change_result
       class(reversible_change), intent(inout) :: this
       type(change_result)     , intent(inout) :: result
     end subroutine change_apply

     subroutine change_check(this, result)
       import :: reversible_change, change_result
       class(reversible_change), intent(inout) :: this
       type(change_result)     , intent(inout) :: result
     end subroutine change_check

     subroutine change_keep(this, result)
       import :: reversible_change, change_result
       class(reversible_change), intent(inout) :: this
       type(change_result)     , intent(inout) :: result
     end subroutine change_keep

     subroutine change_revert(this, result)
       import :: reversible_change, change_result
       class(reversible_change), intent(inout) :: this
       type(change_result)     , intent(inout) :: result
     end subroutine change_revert

  end interface

  !===================================================================!
  ! The stateless lifecycle runner.
  !===================================================================!

  type :: change_controller

   contains

     procedure :: run

  end type change_controller

contains

  !===================================================================!
  ! The result verbs: each marks exactly what its name says, and
  ! reset returns the whole record to false.
  !===================================================================!

  pure subroutine reset(this)

    class(change_result), intent(inout) :: this

    this % attempted    = .false.
    this % applied      = .false.
    this % checked      = .false.
    this % check_passed = .false.
    this % accepted     = .false.
    this % kept         = .false.
    this % reverted     = .false.
    this % failed       = .false.

    this % touches_structure = .false.
    this % touches_value     = .false.

  end subroutine reset

  pure subroutine mark_attempted(this)
    class(change_result), intent(inout) :: this
    this % attempted = .true.
  end subroutine mark_attempted

  pure subroutine mark_applied(this)
    class(change_result), intent(inout) :: this
    this % applied = .true.
  end subroutine mark_applied

  pure subroutine mark_checked(this, pass)
    class(change_result), intent(inout) :: this
    logical                 , intent(in)    :: pass
    this % checked      = .true.
    this % check_passed = pass
  end subroutine mark_checked

  pure subroutine mark_kept(this)
    class(change_result), intent(inout) :: this
    this % kept = .true.
  end subroutine mark_kept

  pure subroutine mark_reverted(this)
    class(change_result), intent(inout) :: this
    this % reverted = .true.
  end subroutine mark_reverted

  pure subroutine mark_failed(this)
    class(change_result), intent(inout) :: this
    this % failed = .true.
  end subroutine mark_failed

  !===================================================================!
  ! Check a terminal record for contradictions, each stopping the
  ! program: kept and reverted cannot both be set; accepted
  ! requires kept; failed excludes kept.
  !===================================================================!

  pure subroutine validate_terminal(this)

    class(change_result), intent(in) :: this

    if (this % kept .and. this % reverted) then
       error stop 'change_result: terminal state is consistent'
    end if

    if (this % accepted .and. .not. this % kept) then
       error stop 'change_result: terminal state is consistent'
    end if

    if (this % failed .and. this % kept) then
       error stop 'change_result: terminal state is consistent'
    end if

  end subroutine validate_terminal

  !===================================================================!
  ! Run the lifecycle: apply, then check, then keep on acceptance
  ! or revert otherwise. A failure reported by apply or check is
  ! reverted and returned. A step that returns without marking its
  ! work stops the program - including every controller-called
  ! revert, on failure and reject paths alike - so no terminal
  ! record can be silently incomplete.
  !===================================================================!

  subroutine run(this, change, accept, result)

    class(change_controller), intent(in)    :: this
    class(reversible_change), intent(inout) :: change
    logical                     , intent(in)    :: accept
    type(change_result)     , intent(inout) :: result

    call result % reset()
    call result % mark_attempted()

    call change % apply(result)

    if (result % failed) then
       call change % revert(result)
       if (.not. result % reverted) then
          error stop 'change_controller: reverted change reports reverted'
       end if
       call result % validate_terminal()
       return
    end if

    if (.not. result % applied) then
       error stop 'change_controller: applied change reports applied'
    end if

    call change % check(result)

    if (.not. result % failed) then
       if (.not. result % checked) then
          error stop 'change_controller: checked change reports checked'
       end if
    end if

    if (result % failed .or. .not. result % check_passed) then
       call change % revert(result)
       if (.not. result % reverted) then
          error stop 'change_controller: reverted change reports reverted'
       end if
       call result % validate_terminal()
       return
    end if

    if (accept) then
       result % accepted = .true.
       call change % keep(result)
       if (.not. result % kept) then
          error stop 'change_controller: kept change reports kept'
       end if
    else
       call change % revert(result)
       if (.not. result % reverted) then
          error stop 'change_controller: reverted change reports reverted'
       end if
    end if

    call result % validate_terminal()

  end subroutine run

end module graph_change_protocol
