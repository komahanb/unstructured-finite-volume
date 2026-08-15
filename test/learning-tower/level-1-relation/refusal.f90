!=====================================================================!
! LEARNING TOWER . LEVEL 1 . THE REFUSALS
!
! The invalid constructions the relation contract promises to
! refuse, each EXPECTED TO DIE for its stated reason:
!
!      arity        a two-row table for a three-slot signature
!      member       the port value 3 offered in the operation slot,
!                   and 3 is not in O = {1, 2}: the ordered
!                   signature supplies the typing
!      undeclared   a set that never signed
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program learning_level_1_refusal

  use learning_assert, only : SLOT_W, OP_PREDICT, PORT_IN1, PORT_OUT
  use graph_set  , only : index_set
  use graph_relation , only : stored_relation

  implicit none

  type(index_set)     :: v, o, p, raw
  type(stored_relation) :: flow
  character(len=32)     :: which

  which = ''
  call get_command_argument(1, which)

  v = index_set('value-slots', 5)
  o = index_set('operations' , 2)
  p = index_set('ports'      , 3)

  select case (trim(which))

  case ('arity')
     flow = stored_relation('bad', [o, v, p], &
          & reshape([OP_PREDICT, SLOT_W], [2, 1]))

  case ('member')
     ! PORT_OUT = 3 in the operation slot: 3 is not in O = {1, 2}.
     flow = stored_relation('bad', [o, v, p], &
          & reshape([PORT_OUT, SLOT_W, PORT_IN1], [3, 1]))

  case ('undeclared')
     flow = stored_relation('bad', [o, v, raw], &
          & reshape([OP_PREDICT, SLOT_W, PORT_IN1], [3, 1]))

  case default
     write(*,'(1x,a)') "usage: refusal arity|member|undeclared"
     error stop 'no case chosen'

  end select

  ! Reaching this line is the failure.
  write(*,'(1x,a,a)') "REACHED PAST THE REFUSAL: ", trim(which)

end program learning_level_1_refusal
