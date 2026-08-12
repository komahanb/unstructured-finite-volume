!=====================================================================!
! The carrier refusal: a domain never signs twice.
!
! This program is EXPECTED TO DIE. It declares a domain and then
! asks it to sign again; the second signing must stop the program
! loudly, because a silent second stamp would let one question
! answer two ways in one lifetime. run.sh asserts the death.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program carrier_refusal

  use graph_carrier, only : counted_set

  implicit none

  type(counted_set) :: cells

  cells = counted_set('cells', 3)

  call cells % declare('cells-again')

  ! Reaching this line is the failure.
  write(*,'(1x,a)') "REACHED PAST A SECOND DECLARE"

end program carrier_refusal
