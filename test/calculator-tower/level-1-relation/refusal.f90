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
!      undeclared   a set that never signed
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program calculator_level_1_refusal

  use calculator_assert, only : SLOT_A, OP_PLUS, PORT_IN1, PORT_OUT
  use graph_set    , only : index_set
  use graph_relation   , only : stored_relation

  implicit none

  type(index_set)     :: x, o, p, raw
  type(stored_relation) :: flow
  character(len=32)     :: which

  which = ''
  call get_command_argument(1, which)

  x = index_set('value-slots', 5)
  o = index_set('operations' , 2)
  p = index_set('ports'      , 3)

  select case (trim(which))

  case ('arity')
     ! Two rows for a three-slot signature.
     flow = stored_relation('bad', [o, x, p], &
          & reshape([OP_PLUS, SLOT_A], [2, 1]))

  case ('member')
     ! PORT_OUT = 3 in the operation slot: 3 is not in O = {1, 2}.
     flow = stored_relation('bad', [o, x, p], &
          & reshape([PORT_OUT, SLOT_A, PORT_IN1], [3, 1]))

  case ('undeclared')
     ! raw never signed; a signature refers to declared domains only.
     flow = stored_relation('bad', [o, x, raw], &
          & reshape([OP_PLUS, SLOT_A, PORT_IN1], [3, 1]))

  case default
     write(*,'(1x,a)') "usage: refusal arity|member|undeclared"
     error stop 'no case chosen'

  end select

  ! Reaching this line is the failure.
  write(*,'(1x,a,a)') "REACHED PAST THE REFUSAL: ", trim(which)

end program calculator_level_1_refusal
