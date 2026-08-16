!=====================================================================!
! REFUSALS OF THE RELATIONAL VIEW
!
! Three laws, each of which must terminate the program with its own
! message. Driven by run.sh.
!
! The binding is storage, and its laws are storage laws: an element it
! does not hold cannot be resolved, and a binding that has lent
! pointers cannot be replaced.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program refusal

  use fractal_graph        , only : graph
  use graph_carrier        , only : member_set, counted_set
  use graph_relation       , only : relation, stored_relation, slot
  use graph_relational_view, only : relational_binding

  implicit none

  character(len=32) :: which

  call get_command_argument(1, which)

  select case (trim(which))

     !================================================================!
     ! The binding is partial: an unbound element has no object.
     !================================================================!

  case ('nosetbound')
     block
       type(graph), target        :: e1, e2
       type(relational_binding)   :: b
       type(counted_set)          :: s
       class(member_set), pointer :: p
       call e1 % declare(); call e2 % declare()
       s = counted_set('E', 4)
       call b % bind_set(e1, s)
       p => b % set_for(e2)
       print *, p % name()
     end block

  case ('norelationbound')
     block
       type(graph), target      :: e1, e2
       type(relational_binding) :: b
       type(counted_set)        :: s
       type(stored_relation)    :: r1
       class(relation), pointer :: p
       call e1 % declare(); call e2 % declare()
       s  = counted_set('E', 4)
       r1 = stored_relation('R1', [slot(s)], reshape([1, 2], [1, 2]))
       call b % bind_relation(e1, r1)
       p => b % relation_for(e2)
       print *, p % name()
     end block

     !================================================================!
     ! Extension preserves what a binding lent; replacement cannot, so
     ! replacement is not an operation. Measured before it was refused:
     ! a deep copy finalizes its left-hand side first, and the lender's
     ! borrowers then read storage it had already freed.
     !================================================================!

  case ('assign')
     block
       type(graph), target      :: e1, e2
       type(relational_binding) :: b, d
       type(counted_set)        :: s
       type(stored_relation)    :: r1, r2
       class(relation), pointer :: p
       call e1 % declare(); call e2 % declare()
       s  = counted_set('E', 4)
       r1 = stored_relation('R1', [slot(s)], reshape([1, 2], [1, 2]))
       r2 = stored_relation('R2', [slot(s)], reshape([3, 4], [1, 2]))
       call b % bind_relation(e1, r1)
       call d % bind_relation(e2, r2)
       p => b % relation_for(e1)                      ! b has lent
       b = d
       print *, p % name()
     end block

  case default
     error stop 'refusal: unknown case'

  end select

end program refusal
