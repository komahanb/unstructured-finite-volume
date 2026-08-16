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
  use graph_carrier         , only : counted_set
  use graph_binary_relation , only : csr_relation, transposed_view, &
       &                             transpose_of
  use fractal_graph         , only : graph
  use graph_relational_view , only : relational_binding

  implicit none

  type(counted_set)          :: o
  type(csr_relation), target :: dep
  type(transposed_view)      :: flipped
  type(graph)       , target :: elem
  type(relational_binding)   :: bnd
  character(len=32)          :: which

  which = ''
  call get_command_argument(1, which)

  o = counted_set('operations', 2)
  call elem % declare()

  select case (trim(which))

  case ('view')
     dep = csr_relation('precedes', o, o, &
          & reshape([OP_PLUS, OP_TIMES], [2, 1]))
     flipped = transpose_of(dep)
     call bnd % bind_relation(elem, flipped)

  case default
     write(*,'(1x,a)') "usage: refusal view"
     error stop 'no case chosen'

  end select

  ! Reaching this line is the failure.
  write(*,'(1x,a,a)') "REACHED PAST THE REFUSAL: ", trim(which)

end program calculator_level_3_refusal
