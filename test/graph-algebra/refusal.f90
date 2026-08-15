!=====================================================================!
! The algebra refusals, each EXPECTED TO DIE for its stated reason:
!
!      slot        a slot index outside the signature
!      embed       a restriction domain from another world
!      none        a projection choosing no slots
!      range       a projection choosing a slot that is not there
!      repeat      a projection choosing one slot twice
!      binary      composition offered a ternary relation
!      middle      composition across two different middle domains
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program algebra_refusal

  use graph_set         , only : index_set, subset
  use graph_relation        , only : stored_relation
  use graph_relation_algebra, only : restrict_slot, project_slots, &
       &                             compose_binary
  use graph_binary_relation , only : csr_relation

  implicit none

  type(index_set)     :: a, b, c, other
  type(subset)      :: outside_slot_domain
  type(stored_relation) :: r, r_ab, r_bc, fat
  type(stored_relation) :: out1
  type(csr_relation)    :: out2
  character(len=32)     :: which

  which = ''
  call get_command_argument(1, which)

  a = index_set('a-things', 3)
  b = index_set('b-things', 4)
  c = index_set('c-things', 2)

  r = stored_relation('r', [a, b, c], &
       & reshape([1,1,1,  2,2,2], [3, 2]))

  select case (trim(which))

  case ('slot')
     out1 = restrict_slot(r, 4, b)

  case ('embed')
     other   = index_set('elsewhere', 4)
     outside_slot_domain = subset('not a subobject of slot 2', other, [2])
     out1 = restrict_slot(r, 2, outside_slot_domain)

  case ('none')
     out1 = project_slots(r, [integer ::])

  case ('range')
     out1 = project_slots(r, [1, 5])

  case ('repeat')
     out1 = project_slots(r, [1, 1])

  case ('binary')
     r_bc = stored_relation('bc', [b, c], reshape([1, 1], [2, 1]))
     out2 = compose_binary(r, r_bc)

  case ('middle')
     r_ab = stored_relation('ab', [a, b], reshape([1, 1], [2, 1]))
     fat  = stored_relation('cb', [c, b], reshape([1, 1], [2, 1]))
     out2 = compose_binary(r_ab, fat)

  case default
     write(*,'(1x,a)') &
          & "usage: refusal slot|embed|none|range|repeat|binary|middle"
     error stop 'no case chosen'

  end select

  ! Reaching this line is the failure.
  write(*,'(1x,a,a)') "REACHED PAST THE REFUSAL: ", trim(which)

end program algebra_refusal
