!=====================================================================!
! The set refusals, each EXPECTED TO DIE for its stated reason:
!
!      twice        a second signing of a declared domain
!      outsider     a subset member from beyond the ambient
!      unhosted     a subset carved from an unsigned ambient
!      undeclared   a graph seat holding an unsigned domain
!      dupset       one domain seated twice in one graph
!      slot         set_at asked outside {1 .. num_sets()}
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program set_refusal

  use graph_set, only : set, index_set, subset
  use graph_set, only : unrelated_graph, declared_set

  implicit none

  type(index_set)       :: cells, raw
  type(subset)          :: bad
  type(unrelated_graph), target :: g
  class(set), pointer   :: p
  character(len=32)     :: which

  which = ''
  call get_command_argument(1, which)

  cells = index_set('cells', 3)

  select case (trim(which))

  case ('twice')
     call cells % declare('cells-again')

  case ('outsider')
     bad = subset('bad', cells, [2, 7])

  case ('unhosted')
     bad = subset('bad', raw, [1])

  case ('undeclared')
     ! raw never signed, so no graph may seat it.
     g = unrelated_graph('bad', [declared_set(cells), declared_set(raw)])

  case ('dupset')
     g = unrelated_graph('bad', [declared_set(cells), declared_set(cells)])

  case ('slot')
     g = unrelated_graph('good', [declared_set(cells)])
     p => g % set_at(2)

  case default
     write(*,'(1x,a)') "usage: refusal twice|outsider|unhosted|undeclared|dupset|slot"
     error stop 'no case chosen'

  end select

  ! Reaching this line is the failure.
  write(*,'(1x,a,a)') "REACHED PAST THE REFUSAL: ", trim(which)

end program set_refusal
