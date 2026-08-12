!=====================================================================!
! The interpretation and algorithm refusals, each EXPECTED TO DIE
! for its stated reason:
!
!      notowned     the selector names a relation the graph lacks
!      notbinary    a ternary relation offered as adjacency
!      notsquare    a binary relation over two different domains
!      cycle        a lawful cyclic view whose topological order
!                   is rightly refused - the VIEW is valid; only
!                   the ordering is undefined
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program algorithms_refusal

  use graph_carrier        , only : counted_set
  use graph_relation       , only : stored_relation, slot
  use graph_binary_relation, only : csr_relation
  use graph_structure      , only : relational_graph, held_set, &
       &                            held_relation
  use graph_profile        , only : directed_adjacency_view
  use graph_algorithms     , only : topological_order

  implicit none

  type(counted_set)              :: a, b, c
  type(csr_relation)             :: adj, stray, lopsided, ring
  type(stored_relation)          :: fat
  type(relational_graph), target :: g
  type(directed_adjacency_view)  :: view
  integer, allocatable           :: order(:)
  character(len=32)              :: which

  which = ''
  call get_command_argument(1, which)

  a = counted_set('a-things', 3)
  b = counted_set('b-things', 2)

  select case (trim(which))

  case ('notowned')
     adj   = csr_relation('inside' , a, a, reshape([1, 2], [2, 1]))
     stray = csr_relation('outside', a, a, reshape([2, 3], [2, 1]))
     g = relational_graph('bad', [held_set(a)], [held_relation(adj)])
     view = directed_adjacency_view(g, stray)

  case ('notbinary')
     c   = counted_set('c-things', 2)
     fat = stored_relation('fat', [slot(a), slot(a), slot(c)], &
          & reshape([1, 2, 1], [3, 1]))
     g = relational_graph('bad', [held_set(a), held_set(c)], &
          & [held_relation(fat)])
     view = directed_adjacency_view(g, fat)

  case ('notsquare')
     lopsided = csr_relation('lopsided', a, b, reshape([1, 1], [2, 1]))
     g = relational_graph('bad', [held_set(a), held_set(b)], &
          & [held_relation(lopsided)])
     view = directed_adjacency_view(g, lopsided)

  case ('cycle')
     ring = csr_relation('ring', a, a, &
          & reshape([1,2,  2,3,  3,1], [2, 3]))
     g = relational_graph('cyclic', [held_set(a)], [held_relation(ring)])
     ! The view is lawful: a cycle is a directed graph.
     view = directed_adjacency_view(g, ring)
     ! Only the ordering is undefined.
     call topological_order(view, order)

  case default
     write(*,'(1x,a)') "usage: refusal notowned|notbinary|notsquare|cycle"
     error stop 'no case chosen'

  end select

  ! Reaching this line is the failure.
  write(*,'(1x,a,a)') "REACHED PAST THE REFUSAL: ", trim(which)

end program algorithms_refusal
