!=====================================================================!
! THE SEQUENCE VIEW REFUSALS
!
! Two kinds, kept apart:
!
!     malformed   a cell whose branch(1) is not KNOWN. The
!                 representation is wrong, and every traversal refuses.
!     unknown     the representation is right and the answer is not
!                 determined. Refused only where the answer depends on
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

  character(len=32)    :: which
  type(graph), pointer :: e
  integer              :: n
  logical              :: yes

  call get_command_argument(1, which)

  select case (trim(which))

     !================================================================!
     ! Malformed representation: a cell must hold a KNOWN element.
     !================================================================!

  case ('cellnull')
     block
       type(graph), target :: holder, cell
       call holder % declare(); call cell % declare()
       cell % branch(1) = null_branch()          ! no element
       cell % branch(2) = null_branch()
       holder % branch(1) = known_branch(cell)
       n = sequence_num_elements(holder % branch(1))
       print *, n
     end block

  case ('cellunknown')
     block
       type(graph), target :: holder, cell
       call holder % declare(); call cell % declare()
       cell % branch(1) = unknown_branch()       ! element not known
       cell % branch(2) = null_branch()
       holder % branch(1) = known_branch(cell)
       n = sequence_num_elements(holder % branch(1))
       print *, n
     end block

     !================================================================!
     ! The extent depends on the unknown part.
     !================================================================!

  case ('sizeunknownholder')
     block
       type(graph), target :: holder
       call holder % declare()
       holder % branch(1) = unknown_branch()
       n = sequence_num_elements(holder % branch(1))
       print *, n
     end block

  case ('sizeunknowntail')
     block
       type(graph), target :: holder, cell, elem
       call holder % declare(); call cell % declare(); call elem % declare()
       cell % branch(1) = known_branch(elem)
       cell % branch(2) = unknown_branch()
       holder % branch(1) = known_branch(cell)
       n = sequence_num_elements(holder % branch(1))
       print *, n
     end block

  case ('containsunknowntail')
     block
       type(graph), target :: holder, cell, elem, stranger
       call holder % declare(); call cell % declare()
       call elem % declare(); call stranger % declare()
       cell % branch(1) = known_branch(elem)
       cell % branch(2) = unknown_branch()
       holder % branch(1) = known_branch(cell)
       yes = sequence_has(holder % branch(1), stranger)
       print *, yes
     end block

     !================================================================!
     ! Indexing.
     !================================================================!

  case ('indexzero')
     block
       type(graph), target :: holder, cell, elem
       call holder % declare(); call cell % declare(); call elem % declare()
       cell % branch(1) = known_branch(elem)
       cell % branch(2) = null_branch()
       holder % branch(1) = known_branch(cell)
       e => sequence_element(holder % branch(1), 0)
     end block

  case ('pastend')
     block
       type(graph), target :: holder, cell, elem
       call holder % declare(); call cell % declare(); call elem % declare()
       cell % branch(1) = known_branch(elem)
       cell % branch(2) = null_branch()
       holder % branch(1) = known_branch(cell)
       e => sequence_element(holder % branch(1), 2)
     end block

  case ('pastunknown')
     block
       type(graph), target :: holder, cell, elem
       call holder % declare(); call cell % declare(); call elem % declare()
       cell % branch(1) = known_branch(elem)
       cell % branch(2) = unknown_branch()
       holder % branch(1) = known_branch(cell)
       e => sequence_element(holder % branch(1), 2)
     end block

  case ('emptyindexed')
     block
       type(graph), target :: holder
       call holder % declare()
       holder % branch(1) = null_branch()
       e => sequence_element(holder % branch(1), 1)
     end block

  case default
     error stop 'refusal: no such case'

  end select

  print *, 'refusal: unreachable'
  error stop 1

end program refusal
