!=====================================================================!
! THE EPISTEMIC VIEW REFUSALS
!
! Three laws, each of which must terminate the program with its own
! message. The retired suite also refused a second identity assignment
! and a fifth state; those belong to the kernel and to its own suite
! now, and are not restated here.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program refusal

  use graph_fractal       , only : graph, null_branch, unknown_branch, known_branch
  use graph_epistemic_view, only : epistemic_name, data_of, residual_of

  implicit none

  character(len=32)    :: which
  type(graph), pointer :: p

  call get_command_argument(1, which)

  select case (trim(which))

     !================================================================!
     ! Neither NULL nor UNKNOWN is a value; no accessor manufactures
     ! one.
     !================================================================!

  case ('dataunknown')
     block
       type(graph), target :: g
       call g % declare()
       g % branch = [unknown_branch(), unknown_branch()]
       p => data_of(g)
     end block

  case ('residualunknown')
     block
       type(graph), target :: g, q
       call g % declare(); call q % declare()
       g % branch = [known_branch(q), unknown_branch()]
       p => residual_of(g)
     end block

     !================================================================!
     ! NULL is outside the reading's domain. The 3x3 state space is
     ! not forced back into the 2x2 classification.
     !================================================================!

  case ('nullname')
     block
       type(graph), target :: g, q
       character(len=:), allocatable :: said
       call g % declare(); call q % declare()
       g % branch = [known_branch(q), null_branch()]
       said = epistemic_name(g)
       print *, said
     end block

  case default
     error stop 'refusal: no such case'

  end select

  print *, 'refusal: unreachable'
  error stop 1

end program refusal
