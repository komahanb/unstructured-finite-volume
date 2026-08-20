!=====================================================================!
! CALCULATOR TOWER . LEVEL 1 . THE REFUSALS
!
! The invalid constructions the relation contract promises to
! refuse, each EXPECTED TO DIE for its stated reason:
!
!      arity        a tuple table whose rows are not the signature
!      member       a member its slot's domain does not hold - the
!                   port value 3 offered in the operation slot,
!                   and 3 is not in O = {1, 2}: the ordered
!                   signature supplies the typing
!      undeclared   a carrier that never signed
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program calculator_level_1_refusal

  use calculator_assert, only : SLOT_A, OP_PLUS, PORT_IN1, PORT_OUT
  use graph_fractal        , only : set_graph => graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set        , only : set_map
  use relation_finitary   , only : stored_relation

  implicit none

  type(set_graph)     :: x, o, p, raw
  type(stored_relation) :: flow
  character(len=32)     :: which
  type(set_map)     :: sets

  which = ''
  call get_command_argument(1, which)

  call x % declare()
  call sets % bind(x, counted_set_representation(5))
  call o % declare()
  call sets % bind(o, counted_set_representation(2))
  call p % declare()
  call sets % bind(p, counted_set_representation(3))

  select case (trim(which))

  case ('arity')
     ! Two rows for a three-slot signature.
     flow = stored_relation('bad', [o, x, p], &
          & reshape([OP_PLUS, SLOT_A], [2, 1]), sets)

  case ('member')
     ! PORT_OUT = 3 in the operation slot: 3 is not in O = {1, 2}.
     flow = stored_relation('bad', [o, x, p], &
          & reshape([PORT_OUT, SLOT_A, PORT_IN1], [3, 1]), sets)

  case ('undeclared')
     ! raw never signed; a signature refers to declared domains only.
     flow = stored_relation('bad', [o, x, raw], &
          & reshape([OP_PLUS, SLOT_A, PORT_IN1], [3, 1]), sets)

  case default
     write(*,'(1x,a)') "usage: refusal arity|member|undeclared"
     error stop 'no case chosen'

  end select

  ! Reaching this line is the failure.
  write(*,'(1x,a,a)') "REACHED PAST THE REFUSAL: ", trim(which)

end program calculator_level_1_refusal
