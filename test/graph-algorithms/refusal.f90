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

  use graph_set        , only : index_set
  use graph_relation       , only : stored_relation, declared_domain
  use graph_binary_relation, only : csr_relation
  use graph_structure      , only : related_graph, declared_set, &
       &                            declared_relation
  use graph_interpretation        , only : directed_adjacency_view
  use graph_algorithms     , only : topological_order

  implicit none

  type(index_set)              :: a, b, c
  type(csr_relation)             :: adj, stray, lopsided, ring
  type(stored_relation)          :: fat
  type(related_graph), target :: g
  type(directed_adjacency_view)  :: view
  integer, allocatable           :: order(:)
  character(len=32)              :: which

  which = ''
  call get_command_argument(1, which)

  a = index_set('a-things', 3)
  b = index_set('b-things', 2)

  select case (trim(which))

  case ('notowned')
     adj   = csr_relation('inside' , a, a, reshape([1, 2], [2, 1]))
     stray = csr_relation('outside', a, a, reshape([2, 3], [2, 1]))
     g = related_graph('bad', [declared_set(a)], [declared_relation(adj)])
     view = directed_adjacency_view(g, stray)

  case ('notbinary')
     c   = index_set('c-things', 2)
     fat = stored_relation('fat', [declared_domain(a), declared_domain(a), declared_domain(c)], &
          & reshape([1, 2, 1], [3, 1]))
     g = related_graph('bad', [declared_set(a), declared_set(c)], &
          & [declared_relation(fat)])
     view = directed_adjacency_view(g, fat)

  case ('notsquare')
     lopsided = csr_relation('lopsided', a, b, reshape([1, 1], [2, 1]))
     g = related_graph('bad', [declared_set(a), declared_set(b)], &
          & [declared_relation(lopsided)])
     view = directed_adjacency_view(g, lopsided)

  case ('cycle')
     ring = csr_relation('ring', a, a, &
          & reshape([1,2,  2,3,  3,1], [2, 3]))
     g = related_graph('cyclic', [declared_set(a)], [declared_relation(ring)])
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
