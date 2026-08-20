!=====================================================================!
! CALCULATOR TOWER . LEVEL 3 . THE REFUSALS
!
! One law is a refusal, and it is a STORAGE law:
!
!      view      a borrowing view offered for binding - built with
!                the lower-level binary infrastructure, never with an
!                ordinary-graph profile
!
! The binding owns what it stores, and copying a borrowing view into
! owned storage copies a reference to a base the binding does not keep
! alive. It would own the view and not what makes it true.
!
! The three laws that used to stand here - a foreign carrier, a domain
! seated twice, a relation seated twice - are now properties of a
! representation already built, so the view ANSWERS them rather than
! refusing them. They are checked in test.f90.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program calculator_level_3_refusal

  use calculator_assert     , only : OP_PLUS, OP_TIMES
  use graph_fractal        , only : set_graph => graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set        , only : set_map
  use relation_binary , only : csr_relation, transposed_view, &
       &                             transpose_of
  use graph_fractal         , only : graph
  use view_relational , only : relational_binding

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
          & reshape([OP_PLUS, OP_TIMES], [2, 1]), sets)
     flipped = transpose_of(dep)
     call bnd % bind_relation(elem, flipped)

  case default
     write(*,'(1x,a)') "usage: refusal view"
     error stop 'no case chosen'

  end select

  ! Reaching this line is the failure.
  write(*,'(1x,a,a)') "REACHED PAST THE REFUSAL: ", trim(which)

end program calculator_level_3_refusal
