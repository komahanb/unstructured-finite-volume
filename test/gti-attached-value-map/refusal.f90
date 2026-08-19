!=====================================================================!
! The refusals that must die at the attached-value seat, one per
! invocation:
!
!      attachundeclared      attaching on an unassigned identity
!      attachtwice           attaching one seat twice
!      markknownundeclared   vouching on an unassigned identity
!      markknownunattached   vouching where no seat exists
!      markunknownunattached withdrawing trust where no seat exists
!      detachundeclared      detaching on an unassigned identity
!      detachunattached      detaching where no seat exists
!      readunattached        reading where no seat exists
!      readunknown           reading a number nobody vouched for
!      knownempty            vouching KNOWN with no numbers behind it
!
! Every case must error stop; a case that returns is a failure of
! the suite.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program refusal

  use iso_fortran_env      , only : dp => REAL64
  use fractal_graph        , only : graph
  use gti_value_buffers    , only : gti_value_buffer
  use gti_attached_value_maps, only : gti_attached_value_map

  implicit none

  type(gti_attached_value_map) :: map
  type(graph) :: declared_graph
  type(graph) :: bare

  type(gti_value_buffer) :: values, empty, stored
  character(len=32) :: which

  call get_command_argument(1, which)

  call declared_graph % declare()
  call values % set_real([1.0_dp])

  select case (trim(which))

  case ('attachundeclared')

     call map % attach_unknown(bare)

  case ('attachtwice')

     call map % attach_unknown(declared_graph)
     call map % attach_unknown(declared_graph)

  case ('markknownundeclared')

     call map % mark_known(bare, values)

  case ('markknownunattached')

     call map % mark_known(declared_graph, values)

  case ('markunknownunattached')

     call map % mark_unknown(declared_graph)

  case ('detachundeclared')

     call map % detach(bare)

  case ('detachunattached')

     call map % detach(declared_graph)

  case ('readunattached')

     call map % value_of(declared_graph, stored)

  case ('readunknown')

     call map % attach_unknown(declared_graph)
     call map % value_of(declared_graph, stored)

  case ('knownempty')

     call map % attach_unknown(declared_graph)
     call map % mark_known(declared_graph, empty)

  case default

     error stop 'refusal: unknown case'

  end select

  write(*,*) 'refusal case survived: ', trim(which)

end program refusal
