!=====================================================================!
! The structure refusals, each EXPECTED TO DIE for its stated
! reason:
!
!      foreign      a relation over a domain the graph does not hold
!      twice        a second signing of a declared graph
!      undeclared   a seat holding a domain that never signed
!      dupset       one domain seated twice
!      duprel       one relation seated twice
!      view         a borrowing view offered for ownership
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program structure_refusal

  use graph_carrier        , only : counted_set
  use graph_binary_relation, only : csr_relation, transposed_view, &
       &                            transpose_of
  use graph_structure      , only : relational_graph, held_set, &
       &                            held_relation

  implicit none

  type(counted_set)          :: ops, foreign, raw
  type(csr_relation), target :: dep
  type(transposed_view)      :: flipped
  type(relational_graph)     :: g
  type(held_relation)        :: none(0)
  character(len=32)          :: which

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

  case ('undeclared')
     g = relational_graph('bad', [held_set(ops), held_set(raw)], none)

  case ('dupset')
     g = relational_graph('bad', [held_set(ops), held_set(ops)], none)

  case ('duprel')
     dep = csr_relation('feeds', ops, ops, reshape([1, 2], [2, 1]))
     g = relational_graph('bad', [held_set(ops)], &
          & [held_relation(dep), held_relation(dep)])

  case ('view')
     dep = csr_relation('feeds', ops, ops, reshape([1, 2], [2, 1]))
     flipped = transpose_of(dep)
     g = relational_graph('bad', [held_set(ops)], &
          & [held_relation(flipped)])

  case default
     write(*,'(1x,a)') &
          & "usage: refusal foreign|twice|undeclared|dupset|duprel|view"
     error stop 'no case chosen'

  end select

  ! Reaching this line is the failure.
  write(*,'(1x,a,a)') "REACHED PAST THE REFUSAL: ", trim(which)

end program structure_refusal
