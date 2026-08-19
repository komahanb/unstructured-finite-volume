!=====================================================================!
! The refusals that must die on the change seam, one per
! invocation:
!
!      attachtwice   a second seat for one identity
!      updatefree    vouching a value nobody attached
!      readunknown   reading a number nobody vouched for
!      detachfree    detaching a seat that does not exist
!      emptyknown    a KNOWN claim with no numbers behind it
!      undeclared    keying a writer on unassigned identity
!      liarapply     a change that works without marking applied
!      liarrevert    a reverted change that never says reverted
!      impossible    a terminal record kept AND reverted at once
!      unbound       a value change applied with no map bound
!
! Every case must error stop; a case that returns is a failure of
! the suite.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program refusal

  use iso_fortran_env      , only : dp => REAL64
  use fractal_graph        , only : graph
  use graph_change_protocol, only : change_controller, change_result
  use graph_value_map      , only : value_map, VALUE_UNKNOWN
  use graph_value_change   , only : value_change
  use toy_changes          , only : liar_apply, liar_revert

  implicit none

  type(change_controller) :: controller
  type(change_result)     :: result
  type(value_map)         :: map
  type(value_change)      :: loose
  type(liar_apply)        :: worker
  type(liar_revert)       :: mute

  type(graph) :: a, ghost
  real(dp), allocatable :: rv(:)

  character(len=32) :: which

  call get_command_argument(1, which)

  call a % declare()

  select case (trim(which))

  case ('attachtwice')

     call map % attach_unknown(a)
     call map % attach_unknown(a)

  case ('updatefree')

     call map % mark_known(a, [1.0_dp])

  case ('readunknown')

     call map % attach_unknown(a)
     call map % value_of(a, rv)

  case ('detachfree')

     call map % detach(a)

  case ('emptyknown')

     call map % attach_unknown(a)
     call map % mark_known(a, [real(dp) ::])

  case ('undeclared')

     call map % attach_unknown(ghost)

  case ('liarapply')

     call controller % run(worker, .true., result)

  case ('liarrevert')

     call controller % run(mute, .true., result)

  case ('impossible')

     call result % mark_kept()
     call result % mark_reverted()
     call result % validate_terminal()

  case ('unbound')

     call controller % run(loose, .true., result)

  case default

     error stop 'refusal: unknown case'

  end select

  write(*,*) 'refusal case survived: ', trim(which)

end program refusal
