!=====================================================================!
! THE SEQUENCE VIEW REFUSALS
!
! Two distinct conditions:
!
!     malformed   a cell whose branch(1) is not KNOWN. The
!                 representation is wrong, and every traversal refuses.
!     unknown     the representation is right and the result is not
!                 determined. Refused only where the result depends on
!                 the unknown part.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program refusal

  use graph_fractal      , only : graph, &
       & null_branch, unknown_branch, known_branch
  use view_sequence, only : sequence_num_elements, sequence_element, &
       & sequence_has

  implicit none

  character(len=32)    :: case_name
  type(graph), pointer :: e
  integer              :: n
  logical              :: contains_member

  call get_command_argument(1, case_name)

  select case (trim(case_name))

     !================================================================!
     ! Malformed representation: a cell must hold a KNOWN element.
     !================================================================!

  case ('cellnull')
     block
       type(graph), target :: container, cell
       call container % declare(); call cell % declare()
       cell % branch(1) = null_branch()          ! no element
       cell % branch(2) = null_branch()
       container % branch(1) = known_branch(cell)
       n = sequence_num_elements(container % branch(1))
       print *, n
     end block

  case ('cellunknown')
     block
       type(graph), target :: container, cell
       call container % declare(); call cell % declare()
       cell % branch(1) = unknown_branch()       ! element not known
       cell % branch(2) = null_branch()
       container % branch(1) = known_branch(cell)
       n = sequence_num_elements(container % branch(1))
       print *, n
     end block

     !================================================================!
     ! The extent depends on the unknown part.
     !================================================================!

  case ('sizeunknownholder')
     block
       type(graph), target :: container
       call container % declare()
       container % branch(1) = unknown_branch()
       n = sequence_num_elements(container % branch(1))
       print *, n
     end block

  case ('sizeunknowntail')
     block
       type(graph), target :: container, cell, elem
       call container % declare(); call cell % declare(); call elem % declare()
       cell % branch(1) = known_branch(elem)
       cell % branch(2) = unknown_branch()
       container % branch(1) = known_branch(cell)
       n = sequence_num_elements(container % branch(1))
       print *, n
     end block

  case ('containsunknowntail')
     block
       type(graph), target :: container, cell, elem, nonmember
       call container % declare(); call cell % declare()
       call elem % declare(); call nonmember % declare()
       cell % branch(1) = known_branch(elem)
       cell % branch(2) = unknown_branch()
       container % branch(1) = known_branch(cell)
       contains_member = sequence_has(container % branch(1), nonmember)
       print *, contains_member
     end block

     !================================================================!
     ! Indexing.
     !================================================================!

  case ('indexzero')
     block
       type(graph), target :: container, cell, elem
       call container % declare(); call cell % declare(); call elem % declare()
       cell % branch(1) = known_branch(elem)
       cell % branch(2) = null_branch()
       container % branch(1) = known_branch(cell)
       e => sequence_element(container % branch(1), 0)
     end block

  case ('pastend')
     block
       type(graph), target :: container, cell, elem
       call container % declare(); call cell % declare(); call elem % declare()
       cell % branch(1) = known_branch(elem)
       cell % branch(2) = null_branch()
       container % branch(1) = known_branch(cell)
       e => sequence_element(container % branch(1), 2)
     end block

  case ('pastunknown')
     block
       type(graph), target :: container, cell, elem
       call container % declare(); call cell % declare(); call elem % declare()
       cell % branch(1) = known_branch(elem)
       cell % branch(2) = unknown_branch()
       container % branch(1) = known_branch(cell)
       e => sequence_element(container % branch(1), 2)
     end block

  case ('emptyindexed')
     block
       type(graph), target :: container
       call container % declare()
       container % branch(1) = null_branch()
       e => sequence_element(container % branch(1), 1)
     end block

  case default
     error stop 'refusal: no such case'

  end select

  print *, 'refusal: unreachable'
  error stop 1

end program refusal
