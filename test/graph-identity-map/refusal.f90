!=====================================================================!
! REFUSALS OF THE LABEL MAP
!
! Two storage laws. Neither states what a set may be called; both state
! what may be BOUND.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program refusal

  use graph_fractal  , only : graph
  use map_label, only : label_map

  implicit none

  character(len=32) :: which

  call get_command_argument(1, which)

  select case (trim(which))

     !================================================================!
     ! An undeclared token names nothing, so it cannot be the key of
     ! anything either.
     !================================================================!

  case ('unsigned')
     block
       type(graph)     :: g
       type(label_map) :: m
       call m % bind(g, 'cells')
     end block

     !================================================================!
     ! A set is named once. A second binding would leave two labels
     ! for one key, and the later would replace the first silently.
     !================================================================!

  case ('twice')
     block
       type(graph)     :: g
       type(label_map) :: m
       call g % declare()
       call m % bind(g, 'cells')
       call m % bind(g, 'faces')
     end block

  case default
     error stop 'refusal: unknown case'

  end select

end program refusal
