!=====================================================================!
! The binary refusals: three ways to ask for a CSR relation that
! cannot exist, each EXPECTED TO DIE with its own stated reason.
!
!      member       a tuple names a member outside its slot's domain
!      arity        the tuple table's rows are not two
!      undeclared   a set that never signed
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program binary_refusal

  use graph_set        , only : index_set
  use graph_binary_relation, only : csr_relation

  implicit none

  type(index_set)  :: cells, faces, raw
  type(csr_relation) :: r
  character(len=32)  :: which

  which = ''
  call get_command_argument(1, which)

  cells = index_set('cells', 4)
  faces = index_set('faces', 5)

  select case (trim(which))

  case ('member')
     r = csr_relation('bad', cells, faces, reshape([5, 1], [2, 1]))

  case ('arity')
     r = csr_relation('bad', cells, faces, reshape([1, 1, 1], [3, 1]))

  case ('undeclared')
     r = csr_relation('bad', cells, raw, reshape([1, 1], [2, 1]))

  case default
     write(*,'(1x,a)') "usage: refusal member|arity|undeclared"
     error stop 'no case chosen'

  end select

  ! Reaching this line is the failure.
  write(*,'(1x,a,a)') "REACHED PAST THE REFUSAL: ", trim(which)

end program binary_refusal
