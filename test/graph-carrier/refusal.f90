!=====================================================================!
! The carrier refusals, each EXPECTED TO DIE for its stated reason:
!
!      twice        a second signing of a declared domain
!      outsider     a subset member from beyond the ambient
!      unhosted     a subset carved from an unsigned ambient
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program carrier_refusal

  use graph_carrier, only : counted_set, subset_set

  implicit none

  type(counted_set) :: cells, raw
  type(subset_set)  :: bad
  character(len=32) :: which

  which = ''
  call get_command_argument(1, which)

  cells = counted_set('cells', 3)

  select case (trim(which))

  case ('twice')
     call cells % declare('cells-again')

  case ('outsider')
     bad = subset_set('bad', cells, [2, 7])

  case ('unhosted')
     bad = subset_set('bad', raw, [1])

  case default
     write(*,'(1x,a)') "usage: refusal twice|outsider|unhosted"
     error stop 'no case chosen'

  end select

  ! Reaching this line is the failure.
  write(*,'(1x,a,a)') "REACHED PAST THE REFUSAL: ", trim(which)

end program carrier_refusal
