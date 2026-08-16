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
!      undeclared   a carrier that never signed
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program learning_level_1_refusal

  use learning_assert, only : SLOT_W, OP_PREDICT, PORT_IN1, PORT_OUT
  use fractal_graph        , only : set_graph => graph
  use graph_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use graph_set_map        , only : set_map
  use graph_relation , only : stored_relation

  implicit none

  type(set_graph)     :: v, o, p, raw
  type(stored_relation) :: flow
  character(len=32)     :: which
  type(set_map)     :: sets

  which = ''
  call get_command_argument(1, which)

  call v % declare()
  call sets % bind(v, counted_set_representation(5))
  call o % declare()
  call sets % bind(o, counted_set_representation(2))
  call p % declare()
  call sets % bind(p, counted_set_representation(3))

  select case (trim(which))

  case ('arity')
     flow = stored_relation('bad', [o, v, p], &
          & reshape([OP_PREDICT, SLOT_W], [2, 1]), sets)

  case ('member')
     ! PORT_OUT = 3 in the operation slot: 3 is not in O = {1, 2}.
     flow = stored_relation('bad', [o, v, p], &
          & reshape([PORT_OUT, SLOT_W, PORT_IN1], [3, 1]), sets)

  case ('undeclared')
     flow = stored_relation('bad', [o, v, raw], &
          & reshape([OP_PREDICT, SLOT_W, PORT_IN1], [3, 1]), sets)

  case default
     write(*,'(1x,a)') "usage: refusal arity|member|undeclared"
     error stop 'no case chosen'

  end select

  ! Reaching this line is the failure.
  write(*,'(1x,a,a)') "REACHED PAST THE REFUSAL: ", trim(which)

end program learning_level_1_refusal
