!=====================================================================!
! LEARNING TOWER . LEVEL 3 . THE REFUSALS
!
! The container refusals that define what owned graph structure
! means, each EXPECTED TO DIE for its stated reason:
!
!      unheld   a relation over a set the graph does not hold
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
  use graph_set        , only : index_set
  use graph_binary_relation, only : csr_relation, transposed_view, &
       &                            transpose_of
  use graph_structure      , only : related_graph, declared_set, &
       &                            declared_relation

  implicit none

  type(index_set)          :: v, o
  type(csr_relation), target :: dep
  type(transposed_view)      :: flipped
  type(related_graph)     :: g
  character(len=32)          :: which

  which = ''
  call get_command_argument(1, which)

  v = index_set('value-slots', 5)
  o = index_set('operations' , 2)

  select case (trim(which))

  case ('unheld-domain')
     dep = csr_relation('uses', o, v, reshape([OP_PREDICT, SLOT_W], [2, 1]))
     g = related_graph('bad', [declared_set(o)], [declared_relation(dep)])

  case ('dupset')
     g = related_graph('bad', [declared_set(o), declared_set(o)], &
          & [declared_relation ::])

  case ('duprel')
     dep = csr_relation('precedes', o, o, &
          & reshape([OP_PREDICT, OP_ERROR], [2, 1]))
     g = related_graph('bad', [declared_set(o)], &
          & [declared_relation(dep), declared_relation(dep)])

  case ('view')
     dep = csr_relation('precedes', o, o, &
          & reshape([OP_PREDICT, OP_ERROR], [2, 1]))
     flipped = transpose_of(dep)
     g = related_graph('bad', [declared_set(o)], [declared_relation(flipped)])

  case default
     write(*,'(1x,a)') "usage: refusal unheld-domain|dupset|duprel|view"
     error stop 'no case chosen'

  end select

  ! Reaching this line is the failure.
  write(*,'(1x,a,a)') "REACHED PAST THE REFUSAL: ", trim(which)

end program learning_level_3_refusal
