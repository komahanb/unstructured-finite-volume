!=====================================================================!
! REFUSALS OF THE INCLUSION MAP
!
! Three storage laws. None of them says what a set may be; all three
! say what may be DECLARED.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program refusal

  use fractal_graph      , only : graph
  use graph_inclusion_map, only : inclusion_map, declared_subobject

  implicit none

  character(len=32) :: which

  call get_command_argument(1, which)

  select case (trim(which))

  case ('unsigned')
     block
       type(graph), target :: s, a
       type(inclusion_map) :: m
       call a % declare()
       call m % include_in(s, a)
     end block

  case ('selfsame')
     block
       type(graph), target :: s
       type(inclusion_map) :: m
       call s % declare()
       call m % include_in(s, s)
     end block

  case ('twoambients')
     block
       type(graph), target :: s, a, b
       type(inclusion_map) :: m
       call s % declare(); call a % declare(); call b % declare()
       call m % include_in(s, a)
       call m % include_in(s, b)
     end block

     !================================================================!
     ! A cycle is not a chain. Declaring A into B and B into A leaves
     ! a walk with no end, and the closure refuses rather than looping.
     !================================================================!

  case ('cycle')
     block
       type(graph), target :: a, b, outside
       type(inclusion_map) :: m
       call a % declare(); call b % declare(); call outside % declare()
       call m % include_in(a, b)
       call m % include_in(b, a)
       print *, declared_subobject(a, outside, m)
     end block

  case default
     error stop 'refusal: unknown case'

  end select

end program refusal
