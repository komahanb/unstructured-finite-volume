!=====================================================================!
! REFUSALS
!
! Two kernel laws and six view laws, each of which must terminate the
! program with its own message. Driven by run.sh.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program refusal

  use fractal_graph, only : graph, null_branch, unknown_branch, known_branch
  use graph_views  , only : dp, attribute_map, evaluate, relation_view

  implicit none

  character(len=32) :: which

  call get_command_argument(1, which)

  select case (trim(which))

     !================================================================!
     ! Kernel: identity is assigned once.
     !================================================================!

  case ('twice')
     block
       type(graph), target :: a
       call a % declare()
       call a % declare()
     end block

     !================================================================!
     ! Kernel: KNOWN requires a graph with assigned identity.
     !================================================================!

  case ('undeclared')
     block
       type(graph), target :: a, b
       call a % declare()
       a % branch(1) = known_branch(b)
     end block

     !================================================================!
     ! Views: the attribute map is partial and refuses outside its
     ! domain.
     !================================================================!

  case ('nonumber')
     block
       type(graph), target :: a
       type(attribute_map) :: m
       real(dp)            :: v
       call a % declare()
       v = evaluate(a, m)
       print *, v
     end block

  case ('nosymbol')
     block
       type(graph), target :: a, x, y
       type(attribute_map) :: m
       real(dp)            :: v
       call a % declare(); call x % declare(); call y % declare()
       a % branch(1) = known_branch(x)
       a % branch(2) = known_branch(y)
       call m % bind(x, number = 1.0_dp)
       call m % bind(y, number = 2.0_dp)
       v = evaluate(a, m)
       print *, v
     end block

  case ('noindex')
     block
       type(graph), target  :: cell, face, node
       type(attribute_map)  :: m
       integer, allocatable :: left(:), right(:)
       call cell % declare(); call face % declare(); call node % declare()
       face % branch(1) = null_branch()
       face % branch(2) = known_branch(cell)      ! cell not in the domain
       node % branch(1) = known_branch(face)
       node % branch(2) = null_branch()
       call relation_view(node, m, left, right)
       print *, left, right
     end block

     !================================================================!
     ! Views: the expression view has no operator for this symbol, and
     ! no value on an UNKNOWN branch.
     !================================================================!

  case ('nooperator')
     block
       type(graph), target :: a, x, y
       type(attribute_map) :: m
       real(dp)            :: v
       call a % declare(); call x % declare(); call y % declare()
       a % branch(1) = known_branch(x)
       a % branch(2) = known_branch(y)
       call m % bind(x, number = 1.0_dp)
       call m % bind(y, number = 2.0_dp)
       call m % bind(a, symbol = '-')
       v = evaluate(a, m)
       print *, v
     end block

  case ('unknownvalue')
     block
       type(graph), target :: a, x
       type(attribute_map) :: m
       real(dp)            :: v
       call a % declare(); call x % declare()
       a % branch(1) = known_branch(x)
       a % branch(2) = unknown_branch()
       call m % bind(x, number = 1.0_dp)
       call m % bind(a, symbol = '+')
       v = evaluate(a, m)
       print *, v
     end block

     !================================================================!
     ! Views: only NULL denotes absence. An UNKNOWN member is not a
     ! boundary, and the relation view must not merge the two.
     !================================================================!

  case ('unknownmember')
     block
       type(graph), target  :: cell, face, node
       type(attribute_map)  :: m
       integer, allocatable :: left(:), right(:)
       call cell % declare(); call face % declare(); call node % declare()
       call m % bind(cell, index = 1)
       face % branch(1) = unknown_branch()
       face % branch(2) = known_branch(cell)
       node % branch(1) = known_branch(face)
       node % branch(2) = null_branch()
       call relation_view(node, m, left, right)
       print *, left, right
     end block

  case default
     error stop 'refusal: no such case'

  end select

  print *, 'refusal: unreachable'
  error stop 1

end program refusal
