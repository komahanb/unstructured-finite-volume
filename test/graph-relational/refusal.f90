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

  use graph_fractal        , only : graph
  use graph_fractal        , only : set_graph => graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set        , only : set_map
  use map_label      , only : label_map
  use relation_finitary       , only : relation, stored_relation
  use relation_binary, only : csr_relation, transposed_relation, &
       &                            transpose_of
  use view_relational, only : relational_binding

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
       type(set_graph)          :: s
       type(set_graph), pointer :: p
       type(set_map)     :: sets
       type(label_map)     :: labels
       call e1 % declare(); call e2 % declare()
       call s % declare()
       call sets % bind(s, counted_set_representation(4))
       call b % bind_set(e1, s)
       p => b % set_for(e2)
       print *, labels % label_of(p)
     end block

  case ('norelationbound')
     block
       type(graph), target      :: e1, e2
       type(relational_binding) :: b
       type(set_graph)        :: s
       type(stored_relation)    :: r1
       class(relation), pointer :: p
       type(set_map)     :: sets
       type(label_map)     :: labels
       call e1 % declare(); call e2 % declare()
       call s % declare()
       call sets % bind(s, counted_set_representation(4))
       r1 = stored_relation('R1', [s], reshape([1, 2], [1, 2]), sets)
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
       type(set_graph)        :: s
       type(stored_relation)    :: r1, r2
       class(relation), pointer :: p
       type(set_map)     :: sets
       type(label_map)     :: labels
       call e1 % declare(); call e2 % declare()
       call s % declare()
       call sets % bind(s, counted_set_representation(4))
       r1 = stored_relation('R1', [s], reshape([1, 2], [1, 2]), sets)
       r2 = stored_relation('R2', [s], reshape([3, 4], [1, 2]), sets)
       call b % bind_relation(e1, r1)
       call d % bind_relation(e2, r2)
       p => b % relation_for(e1)                      ! b has lent
       b = d
       print *, p % name()
     end block

     !================================================================!
     ! The binding stores identified objects. An undeclared object
     ! does not match itself, so nothing the view asks about it can be
     ! answered - not membership, not distinctness.
     !================================================================!

  case ('unsignedset')
     block
       type(graph), target      :: e1
       type(relational_binding) :: b
       type(set_graph)        :: raw
       call e1 % declare()
       call b % bind_set(e1, raw)
     end block

  case ('unsignedrelation')
     block
       type(graph), target      :: e1
       type(relational_binding) :: b
       type(stored_relation)    :: raw
       call e1 % declare()
       call b % bind_relation(e1, raw)
     end block

     !================================================================!
     ! The binding owns whole relations. Copying a borrowing view into
     ! owned storage copies a reference to a base the binding does not
     ! keep alive: it would own the view and not what makes it true.
     ! A view rides above a bound relation, never inside one.
     !================================================================!

  case ('boundview')
     block
       type(graph), target        :: e1
       type(relational_binding)   :: b
       type(set_graph)          :: ops
       type(csr_relation), target :: dep
       type(transposed_relation)      :: flipped
       type(set_map)     :: sets
       type(label_map)     :: labels
       call e1 % declare()
       call ops % declare()
       call sets % bind(ops, counted_set_representation(3))
       dep     = csr_relation('feeds', ops, ops, reshape([1, 2], [2, 1]), sets)
       flipped = transpose_of(dep)
       call b % bind_relation(e1, flipped)
     end block

  case default
     error stop 'refusal: unknown case'

  end select

end program refusal
