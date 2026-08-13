!=====================================================================!
! LEARNING TOWER . LEVEL 3 . THE REFUSALS
!
! The container refusals that define what owned graph structure
! means, each EXPECTED TO DIE for its stated reason:
!
!      foreign   a relation over a carrier the graph does not hold
!      dupset    one domain seated twice
!      duprel    one relation seated twice
!      view      a borrowing view offered for ownership - refused
!                for OWNERSHIP/LIFETIME, not for being binary or
!                transposed; this case alone touches the binary
!                module, and only to build the view
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program learning_level_3_refusal

  use learning_assert      , only : OP_PREDICT, OP_ERROR, SLOT_W
  use graph_carrier        , only : counted_set
  use graph_binary_relation, only : csr_relation, transposed_view, &
       &                            transpose_of
  use graph_structure      , only : relational_graph, held_set, &
       &                            held_relation

  implicit none

  type(counted_set)          :: v, o
  type(csr_relation), target :: dep
  type(transposed_view)      :: flipped
  type(relational_graph)     :: g
  character(len=32)          :: which

  which = ''
  call get_command_argument(1, which)

  v = counted_set('value-slots', 5)
  o = counted_set('operations' , 2)

  select case (trim(which))

  case ('foreign')
     dep = csr_relation('uses', o, v, reshape([OP_PREDICT, SLOT_W], [2, 1]))
     g = relational_graph('bad', [held_set(o)], [held_relation(dep)])

  case ('dupset')
     g = relational_graph('bad', [held_set(o), held_set(o)], &
          & [held_relation ::])

  case ('duprel')
     dep = csr_relation('precedes', o, o, &
          & reshape([OP_PREDICT, OP_ERROR], [2, 1]))
     g = relational_graph('bad', [held_set(o)], &
          & [held_relation(dep), held_relation(dep)])

  case ('view')
     dep = csr_relation('precedes', o, o, &
          & reshape([OP_PREDICT, OP_ERROR], [2, 1]))
     flipped = transpose_of(dep)
     g = relational_graph('bad', [held_set(o)], [held_relation(flipped)])

  case default
     write(*,'(1x,a)') "usage: refusal foreign|dupset|duprel|view"
     error stop 'no case chosen'

  end select

  ! Reaching this line is the failure.
  write(*,'(1x,a,a)') "REACHED PAST THE REFUSAL: ", trim(which)

end program learning_level_3_refusal
