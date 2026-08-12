!=====================================================================!
! The structure refusals, each EXPECTED TO DIE for its stated
! reason:
!
!      foreign      a relation over a domain the graph does not hold
!      twice        a second signing of a declared graph
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program structure_refusal

  use graph_carrier        , only : counted_set
  use graph_binary_relation, only : csr_relation
  use graph_structure      , only : relational_graph, held_set, &
       &                            held_relation

  implicit none

  type(counted_set)      :: ops, foreign
  type(csr_relation)     :: dep
  type(relational_graph) :: g
  type(held_relation)    :: none(0)
  character(len=32)      :: which

  which = ''
  call get_command_argument(1, which)

  ops     = counted_set('operations', 3)
  foreign = counted_set('foreign'   , 3)

  select case (trim(which))

  case ('foreign')
     dep = csr_relation('feeds', ops, foreign, reshape([1, 1], [2, 1]))
     g = relational_graph('bad', [held_set(ops)], [held_relation(dep)])

  case ('twice')
     g = relational_graph('good', [held_set(ops)], none)
     call g % declare('good-again')

  case default
     write(*,'(1x,a)') "usage: refusal foreign|twice"
     error stop 'no case chosen'

  end select

  ! Reaching this line is the failure.
  write(*,'(1x,a,a)') "REACHED PAST THE REFUSAL: ", trim(which)

end program structure_refusal
