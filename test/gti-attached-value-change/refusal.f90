!=====================================================================!
! The refusals that must die at the attached-value change seat,
! one per invocation:
!
!      unbound          apply on a change nobody bound
!      unboundrun       the controller running an unbound change
!      emptyvalue       a bound update carrying no numbers - dies
!                       at the map's own known-value law
!      undeclaredgraph  a bound update on an unassigned identity -
!                       binding does not die; the map's identity
!                       gate does, when apply reaches it
!
! The change adds exactly one refusal of its own; every other law
! stays where it lives, in the map. Every case must error stop; a
! case that returns is a failure of the suite.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program refusal

  use iso_fortran_env      , only : dp => REAL64
  use fractal_graph        , only : graph
  use gti_value_buffers    , only : gti_value_buffer
  use gti_attached_value_maps, only : gti_attached_value_map
  use gti_attached_value_changes, only : gti_attached_value_change
  use gti_change_protocols , only : gti_change_controller, gti_change_result

  implicit none

  type(gti_attached_value_map)    :: map
  type(gti_attached_value_change) :: change
  type(gti_change_controller)     :: controller
  type(gti_change_result)         :: result

  type(graph) :: declared_graph
  type(graph) :: bare

  type(gti_value_buffer) :: values, empty
  character(len=32) :: which

  call get_command_argument(1, which)

  call declared_graph % declare()
  call values % set_real([1.0_dp])

  select case (trim(which))

  case ('unbound')

     call change % apply(result)

  case ('unboundrun')

     call controller % run(change, .true., result)

  case ('emptyvalue')

     call change % bind(map, declared_graph, empty)
     call controller % run(change, .true., result)

  case ('undeclaredgraph')

     call change % bind(map, bare, values)
     call controller % run(change, .true., result)

  case default

     error stop 'refusal: unknown case'

  end select

  write(*,*) 'refusal case survived: ', trim(which)

end program refusal
