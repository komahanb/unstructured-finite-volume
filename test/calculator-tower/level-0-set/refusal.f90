!=====================================================================!
! CALCULATOR TOWER . LEVEL 0 . THE REFUSAL
!
! The one invalid construction the set contract promises to
! refuse (CALCULATOR.md 23): a graph never signs twice. This
! program is EXPECTED TO DIE.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program calculator_level_0_refusal

  use graph_set, only : index_set

  implicit none

  type(index_set) :: x

  x = index_set('value-slots', 5)

  call x % declare('value-slots-again')

  ! Reaching this line is the failure.
  write(*,'(1x,a)') "REACHED PAST A SECOND DECLARE"

end program calculator_level_0_refusal
