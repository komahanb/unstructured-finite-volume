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
  use map_set           , only : set_map
  use map_label         , only : label_map
  use relation_finitary          , only : stored_relation
  use relation_binary   , only : csr_relation
  use graph_algorithms        , only : topological_order

  implicit none

  type(graph)           :: a, b, c
  type(set_map)         :: sets
  type(label_map)       :: labels
  type(csr_relation)    :: lopsided, ring
  type(stored_relation) :: fat
  integer, allocatable  :: order(:)
  character(len=32)     :: which

  which = ''
  call get_command_argument(1, which)

  call a % declare()
  call sets % bind(a, counted_set_representation(3))
  call labels % bind(a, 'a-things')
  call b % declare()
  call sets % bind(b, counted_set_representation(2))
  call labels % bind(b, 'b-things')

  select case (trim(which))

  case ('notbinary')

     call c % declare()
     call sets % bind(c, counted_set_representation(2))
     call labels % bind(c, 'c-things')
     fat = stored_relation('fat', [a, a, c], &
          & reshape([1, 2, 1], [3, 1]), sets)
     call topological_order(fat, sets, order)

  case ('notsquare')

     lopsided = csr_relation('lopsided', a, b, reshape([1, 1], [2, 1]), sets)
     call topological_order(lopsided, sets, order)

  case ('cycle')

     ring = csr_relation('ring', a, a, &
          & reshape([1,2,  2,3,  3,1], [2, 3]), sets)
     call topological_order(ring, sets, order)

  case default
     write(*,'(1x,a)') "usage: refusal notbinary|notsquare|cycle"
     error stop 'no case chosen'

  end select

  ! Reaching this line is the failure.
  write(*,'(1x,a,a)') "REACHED PAST THE REFUSAL: ", trim(which)

end program algorithms_refusal
