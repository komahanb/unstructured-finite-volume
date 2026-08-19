!=====================================================================!
! THE REVERSIBLE CHANGE PROTOCOL
!
! One lifecycle for every reversible mutation of a graph and its
! attached values:
!
!      (G, V) --apply/check--> (G*, V*)
!             --keep---------> (G_new, V_new)
!             --revert-------> (G, V)
!
! where G is graph structure and V is attached values. A change may
! touch structure only, attached values only, or both, and the
! controller does not care which: it owns only the lifecycle
!
!      apply -> check -> keep | revert
!
! while the change object owns every detail. Identity lives in the
! graph. Numbers live in attached values. Updates are reversible
! changes.
!
! Failure is a lawful answer at every step: a change whose apply or
! check marks failed is reverted and reported, never refused. What
! dies loudly is a change that lies about its own lifecycle - one
! that returns from apply, check, keep, or revert without either
! marking the step or marking failure - and a terminal record that
! claims impossible things at once.
!
! The protocol is lifecycle only. It names no graph, no values, no
! forms, no time; it imports nothing at all. Concrete changes live
! with the seats they mutate - the first was transactional
! graph growth - and every future structure amendment or
! attached-value update is one more extension of the same abstract
! change, run by the same controller.
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
  ! The abstract reversible change: four deferred verbs, no state
  ! demanded. What apply builds, revert unbuilds; what check
  ! judges, keep makes permanent. The change owns the details.
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
  ! The stateless lifecycle verb. The types keep their public
  ! singular names; Fortran denies a type its host module's name,
  ! so the module speaks in the plural.
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
  ! A terminal record must not claim impossible things at once: a
  ! change is kept or reverted, never both; an accepted change was
  ! kept; a failed change was not kept.
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
  ! The lifecycle, owned whole: apply, then check, then keep on
  ! acceptance or revert otherwise. A failure reported by apply or
  ! check reverts and returns lawfully; a step that did its work
  ! without marking it is refused - and EVERY revert the controller
  ! calls, on failure paths and reject paths alike, must report
  ! reverted, so no terminal record can be silently incomplete.
  ! The controller never asks what the change touched - structural,
  ! numerical, or mixed is the change's own affair.
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
