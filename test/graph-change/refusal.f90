!=====================================================================!
! Invalid-input cases for the change modules, one per invocation.
! The case is selected by the first command-line argument; each
! must terminate in error stop with the message run.sh expects,
! and a case that returns normally is reported as a failure by
! run.sh.
!
!      attachtwice   attach_unknown called twice for one graph
!      updatefree    mark_known on a graph with no map entry
!      readunknown   value_of on an entry whose status is UNKNOWN
!      detachfree    detach on a graph with no map entry
!      emptyknown    mark_known with a zero-length value array
!      undeclared    a map writer keyed on a graph with no
!                    assigned identity
!      silentapply   a change whose apply marks neither applied
!                    nor failed
!      silentrevert  a change whose revert does not mark reverted
!      impossible    validate_terminal on a record marked both
!                    committed and reverted
!      unbound       a value_change applied with no map bound
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program refusal

  use iso_fortran_env      , only : dp => REAL64
  use graph_fractal        , only : graph
  use map_change_protocol, only : run_change, change_record
  use map_value      , only : value_map, VALUE_UNKNOWN
  use map_value_change   , only : value_change
  use change_fixtures          , only : unreported_apply_change, unreported_revert_change

  implicit none
  type(change_record)         :: result
  type(value_map)             :: map
  type(value_change)          :: unbound_change
  type(unreported_apply_change)   :: unreported_apply
  type(unreported_revert_change)  :: unreported_revert

  type(graph) :: a, undeclared_graph
  real(dp), allocatable :: rv(:)

  character(len=32) :: case_name

  call get_command_argument(1, case_name)

  call a % declare()

  select case (trim(case_name))

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

     call map % attach_unknown(undeclared_graph)

  case ('silentapply')

     call run_change(unreported_apply, .true., result)

  case ('silentrevert')

     call run_change(unreported_revert, .true., result)

  case ('impossible')

     call result % mark_committed()
     call result % mark_reverted()
     call result % validate_terminal()

  case ('unbound')

     call run_change(unbound_change, .true., result)

  case default

     error stop 'refusal: unknown case'

  end select

  write(*,*) 'refusal case survived: ', trim(case_name)

end program refusal
