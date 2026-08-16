!=====================================================================!
! The binary refusals: three ways to ask for a CSR relation that
! cannot exist, each EXPECTED TO DIE with its own stated reason.
!
!      member       a tuple names a member outside its slot's domain
!      arity        the tuple table's rows are not two
!      undeclared   a carrier that never signed
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program binary_refusal

  use fractal_graph        , only : set_graph => graph
  use graph_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use graph_set_map        , only : set_map
  use graph_binary_relation, only : csr_relation

  implicit none

  type(set_graph)  :: cells, faces, raw
  type(csr_relation) :: r
  character(len=32)  :: which
  type(set_map)     :: sets

  which = ''
  call get_command_argument(1, which)

  call cells % declare()
  call sets % bind(cells, counted_set_representation(4))
  call faces % declare()
  call sets % bind(faces, counted_set_representation(5))

  select case (trim(which))

  case ('member')
     r = csr_relation('bad', cells, faces, reshape([5, 1], [2, 1]), sets)

  case ('arity')
     r = csr_relation('bad', cells, faces, reshape([1, 1, 1], [3, 1]), sets)

  case ('undeclared')
     r = csr_relation('bad', cells, raw, reshape([1, 1], [2, 1]), sets)

  case default
     write(*,'(1x,a)') "usage: refusal member|arity|undeclared"
     error stop 'no case chosen'

  end select

  ! Reaching this line is the failure.
  write(*,'(1x,a,a)') "REACHED PAST THE REFUSAL: ", trim(which)

end program binary_refusal
