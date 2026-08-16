!=====================================================================!
! The fractal refusals, each EXPECTED TO DIE for its stated reason:
!
!      unknownanswer   known() asked of an UNKNOWN branch
!      nullanswer      known() asked of a NULL branch
!      realizenull     realizing an absence
!      realizetwice    realizing knowledge again
!      foreign         realizing across universes
!      unhatched       questioning a graph that was never minted
!      novalue         asking a bare atom for a number
!      fourth          naming a status that is not one of three
!      cyclelength     asking a cyclic spine for its extent
!      unknowntail     asking a half-learned sequence for its extent
!      boundaryhead    asking a boundary occurrence for its head
!      unbound         asking an environment for a stranger's value
!      lawless         evaluating an operation no reading binds
!      cyclicexpr      evaluating an expression that never bottoms out
!      unlearnedfar    reading a relation through an unlearned far side
!      numeralmember   binding a compressed numeral as a member
!      blockzero       asking a block for its zeroth value
!      staleuniverse   questioning a universe that was signed over
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program fractal_refusal

  use iso_fortran_env, only : dp => REAL64
  use fractal_graph  , only : graph_arena, graph, null_branch, &
       &                      unknown_branch, known_branch
  use fractal_views  , only : seq_len, head_of, env_of, env_value, &
       &                      evaluate, pair_of, seq_of, branch_name, &
       &                      rel_image_count, block_value

  implicit none

  type(graph_arena), target     :: a
  type(graph)                   :: g, h, r, w, items(1)
  type(graph_arena), target     :: b
  character(len=:), allocatable :: said
  character(len=32)             :: which
  real(dp)                      :: v

  which = ''
  call get_command_argument(1, which)

  a = graph_arena()

  select case (trim(which))

  case ('unknownanswer')
     g = a % node(unknown_branch(), null_branch())
     r = g % known(1)

  case ('nullanswer')
     g = a % node(null_branch(), null_branch())
     r = g % known(1)

  case ('realizenull')
     g = a % node(null_branch(), null_branch())
     h = a % node(null_branch(), null_branch())
     call g % realize(1, h)

  case ('realizetwice')
     h = a % node(null_branch(), null_branch())
     g = a % node(known_branch(h), null_branch())
     call g % realize(1, h)

  case ('foreign')
     b = graph_arena()
     g = a % node(unknown_branch(), null_branch())
     h = b % node(null_branch(), null_branch())
     call g % realize(1, h)

  case ('unhatched')
     write(*,'(1x,i0)') w % status(1)

  case ('novalue')
     g = a % node(null_branch(), null_branch())
     v = g % value()

  case ('fourth')
     said = branch_name(7)

  case ('cyclelength')
     g = a % node(unknown_branch(), unknown_branch())
     h = a % node(known_branch(g), known_branch(g))
     call g % realize(1, h)
     call g % realize(2, h)
     write(*,'(1x,i0)') seq_len(g)

  case ('unknowntail')
     h = a % node(null_branch(), null_branch())
     g = a % node(known_branch(h), unknown_branch())
     write(*,'(1x,i0)') seq_len(g)

  case ('boundaryhead')
     h = a % node(null_branch(), null_branch())
     g = a % node(known_branch(h), null_branch())
     items(1) = g
     g = seq_of(a, items)
     r = head_of(g, 1)

  case ('unbound')
     g = a % node(null_branch(), null_branch())
     h = a % node(null_branch(), null_branch())
     r = env_of(a, g, 1.0_dp)
     v = env_value(r, h)

  case ('lawless')
     g = a % node(null_branch(), null_branch())   ! an op nobody binds
     h = a % node(null_branch(), null_branch())   ! plus
     r = a % node(null_branch(), null_branch())   ! times
     items(1) = a % value_atom(1.0_dp)
     w = seq_of(a, items)
     w = pair_of(a, g, w)
     v = evaluate(w, h, r)

  case ('cyclicexpr')
     h = a % node(null_branch(), null_branch())          ! plus
     r = a % node(null_branch(), null_branch())          ! times
     g = a % node(known_branch(r), unknown_branch())     ! e = times[bottom]
     w = a % node(known_branch(g), null_branch())        ! args = [e]
     call g % realize(2, w)                              ! e = times[e]
     v = evaluate(g, h, r)

  case ('unlearnedfar')
     g = a % node(null_branch(), null_branch())
     h = a % node(known_branch(g), unknown_branch())
     items(1) = h
     w = seq_of(a, items)
     write(*,'(1x,i0)') rel_image_count(w, g)

  case ('numeralmember')
     g = a % value_atom(1.0_dp)
     r = env_of(a, g, 2.0_dp)

  case ('blockzero')
     g = a % value_atom(1.0_dp)
     v = block_value(a, g, 0)

  case ('staleuniverse')
     g = a % node(null_branch(), null_branch())
     a = graph_arena()
     write(*,'(1x,i0)') g % status(1)

  case default
     write(*,'(1x,a)') "usage: refusal <case>"
     error stop 'no case chosen'

  end select

  ! Reaching this line is the failure.
  write(*,'(1x,a,a)') "REACHED PAST THE REFUSAL: ", trim(which)

end program fractal_refusal
