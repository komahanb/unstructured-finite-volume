!=====================================================================!
! LEARNING TOWER . LEVEL 3 . THE REFUSALS
!
! One law is a refusal, and it is a STORAGE law:
!
!      view      a borrowing view offered for binding - refused for
!                OWNERSHIP/LIFETIME, not for being binary or
!                transposed. The binding owns what it stores, and a
!                copied view carries a reference to a base the
!                binding does not keep alive.
!
! A foreign carrier, a domain seated twice and a relation seated
! twice are properties of a representation already built, so the
! view ANSWERS them rather than refusing them. They are checked in
! test.f90.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program learning_level_3_refusal

  use learning_assert      , only : OP_PREDICT, OP_ERROR
  use fractal_graph        , only : set_graph => graph
  use graph_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use graph_set_map        , only : set_map
  use graph_binary_relation, only : csr_relation, transposed_view, &
       &                            transpose_of
  use fractal_graph        , only : graph
  use graph_relational_view, only : relational_binding

  implicit none

  type(set_graph)          :: o
  type(csr_relation), target :: dep
  type(transposed_view)      :: flipped
  type(graph)       , target :: elem
  type(relational_binding)   :: bnd
  character(len=32)          :: which
  type(set_map)     :: sets

  which = ''
  call get_command_argument(1, which)

  call o % declare()
  call sets % bind(o, counted_set_representation(2))
  call elem % declare()

  select case (trim(which))

  case ('view')
     dep = csr_relation('precedes', o, o, &
          & reshape([OP_PREDICT, OP_ERROR], [2, 1]), sets)
     flipped = transpose_of(dep)
     call bnd % bind_relation(elem, flipped)

  case default
     write(*,'(1x,a)') "usage: refusal view"
     error stop 'no case chosen'

  end select

  ! Reaching this line is the failure.
  write(*,'(1x,a,a)') "REACHED PAST THE REFUSAL: ", trim(which)

end program learning_level_3_refusal
