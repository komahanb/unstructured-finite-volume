!=====================================================================!
! The relation refusals: five ways to ask for a relation that cannot
! exist, each EXPECTED TO DIE with its own stated reason. run.sh
! picks the case by argument and asserts both the death and the
! message:
!
!      member       a tuple names a member outside its slot's domain
!      arity        the tuple table's rows do not match the slots
!      undeclared   a signature slot that never signed
!      empty        a signature with no slots at all
!      twice        a second signing of a declared relation
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program relation_refusal

  use fractal_graph        , only : set_graph => graph
  use graph_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use graph_set_map        , only : set_map
  use graph_relation, only : stored_relation

  implicit none

  type(set_graph)              :: cells, faces, raw
  type(set_graph), allocatable :: nothing(:)
  type(stored_relation)          :: r
  character(len=32)              :: which
  type(set_map)     :: sets

  which = ''
  call get_command_argument(1, which)

  call cells % declare()
  call sets % bind(cells, counted_set_representation(4))
  call faces % declare()
  call sets % bind(faces, counted_set_representation(5))

  select case (trim(which))

  case ('member')
     ! Cell 5 does not exist; the domain must refuse it.
     r = stored_relation('bad', [cells, faces], reshape([5, 1], [2, 1]), sets)

  case ('arity')
     ! Three rows for two slots.
     r = stored_relation('bad', [cells, faces], reshape([1, 1, 1], [3, 1]), sets)

  case ('undeclared')
     ! raw never signed; a signature refers to declared domains only.
     r = stored_relation('bad', [cells, raw], reshape([1, 1], [2, 1]), sets)

  case ('empty')
     ! No slots at all.
     allocate(nothing(0))
     r = stored_relation('bad', nothing, reshape([integer ::], [0, 0]), sets)

  case ('twice')
     r = stored_relation('good', [cells, faces], reshape([1, 1], [2, 1]), sets)
     call r % declare('good-again')

  case default
     write(*,'(1x,a)') "usage: refusal member|arity|undeclared|empty|twice"
     error stop 'no case chosen'

  end select

  ! Reaching this line is the failure.
  write(*,'(1x,a,a)') "REACHED PAST THE REFUSAL: ", trim(which)

end program relation_refusal
