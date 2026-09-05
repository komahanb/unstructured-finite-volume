!=====================================================================!
! REFUSALS OF THE SET FOUNDATION
!
! Three, each of which must terminate the program with its own message.
! All three are STORAGE laws: they say what the map may key on and what
! it may evaluate, never what a set may be.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program refusal

  use graph_fractal          , only : graph
  use map_set_representation, only : counted_set_representation
  use map_set          , only : set_map

  implicit none

  character(len=32) :: which

  call get_command_argument(1, which)

  select case (trim(which))

     !================================================================!
     ! The map is keyed on identity, so an undeclared set cannot be a
     ! key: it does not match itself, and nothing could ever find it.
     !================================================================!

  case ('unsigned')
     block
       type(graph), target :: s
       type(set_map)       :: m
       call m % bind(s, counted_set_representation(4))
     end block

     !================================================================!
     ! A set is described once. Two representations of one set would
     ! be two extents for one key.
     !================================================================!

  case ('twice')
     block
       type(graph), target :: s
       type(set_map)       :: m
       call s % declare()
       call m % bind(s, counted_set_representation(4))
       call m % bind(s, counted_set_representation(9))
     end block

     !================================================================!
     ! A set with no representation has no extension in this map, and
     ! the map refuses rather than inventing one.
     !================================================================!

  case ('undescribed')
     block
       type(graph), target :: s, other
       type(set_map)       :: m
       call s % declare(); call other % declare()
       call m % bind(other, counted_set_representation(4))
       print *, m % num_members_of(s)
     end block

  case default
     error stop 'refusal: unknown case'

  end select

end program refusal
