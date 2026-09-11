!=====================================================================!
! Invalid-input cases for the graph algorithms, one per invocation.
! Each must terminate in error stop with the message run.sh expects;
! a case that returns normally is reported as a failure by run.sh.
!
!      notbinary    a ternary relation offered as the adjacency
!      notsquare    a binary relation over two different domains
!      cycle        a lawful cyclic relation whose topological
!                   order is refused - the relation is valid; only
!                   the ordering is undefined on it
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program algorithms_refusal

  use graph_fractal           , only : graph
  use map_set_representation, only : counted_set_representation
  use map_set_store     , only : set_store
  use relation_finitary          , only : stored_relation
  use relation_binary   , only : csr_relation
  use relation_algorithms        , only : topological_order

  implicit none

  type(graph)           :: a, b, c
  type(set_store)       :: sets
  type(csr_relation)    :: unequal_domains, ring
  type(stored_relation) :: ternary_relation
  integer, allocatable  :: order(:)
  character(len=32)     :: case_name

  case_name = ''
  call get_command_argument(1, case_name)

  call a % declare()
  call sets % bind(a, counted_set_representation(3))
  call sets % name(a, 'a-domain')
  call b % declare()
  call sets % bind(b, counted_set_representation(2))
  call sets % name(b, 'b-domain')

  select case (trim(case_name))

  case ('notbinary')

     call c % declare()
     call sets % bind(c, counted_set_representation(2))
     call sets % name(c, 'c-domain')
     ternary_relation = stored_relation('fat', [a, a, c], &
          & reshape([1, 2, 1], [3, 1]), sets % set_map)
     call topological_order(ternary_relation, sets % set_map, order)

  case ('notsquare')

     unequal_domains = csr_relation('lopsided', a, b, reshape([1, 1], [2, 1]), sets % set_map)
     call topological_order(unequal_domains, sets % set_map, order)

  case ('cycle')

     ring = csr_relation('ring', a, a, &
          & reshape([1,2,  2,3,  3,1], [2, 3]), sets % set_map)
     call topological_order(ring, sets % set_map, order)

  case default
     write(*,'(1x,a)') "usage: refusal notbinary|notsquare|cycle"
     error stop 'no case chosen'

  end select

  ! Reaching this line is the failure.
  write(*,'(1x,a,a)') "REACHED PAST THE REFUSAL: ", trim(case_name)

end program algorithms_refusal
