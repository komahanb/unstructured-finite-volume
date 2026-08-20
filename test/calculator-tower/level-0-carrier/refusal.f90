!=====================================================================!
! CALCULATOR TOWER . LEVEL 0 . THE REFUSAL
!
! The one invalid construction the carrier contract promises to
! refuse (CALCULATOR.md 23): a domain never signs twice. This
! program is EXPECTED TO DIE.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program calculator_level_0_refusal

  use graph_fractal        , only : graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set        , only : set_map

  implicit none

  type(graph) :: x
  type(set_map)     :: sets

  call x % declare()
  call sets % bind(x, counted_set_representation(5))

  call x % declare()

  ! Reaching this line is the failure.
  write(*,'(1x,a)') "REACHED PAST A SECOND DECLARE"

end program calculator_level_0_refusal
