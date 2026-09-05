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

  use graph_fractal        , only : graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set        , only : set_map
  use map_inclusion  , only : inclusion_map, declared_subobject
  use relation_finitary        , only : stored_relation
  use relation_algebra, only : restrict_slot, project_slots, &
       &                             compose_binary
  use relation_binary , only : csr_relation

  implicit none

  type(graph)     :: a, b, c, other
  type(graph)      :: foreign
  type(stored_relation) :: r, r_ab, r_bc, fat
  type(stored_relation) :: out1
  type(csr_relation)    :: out2
  character(len=32)     :: which
  type(set_map)     :: sets
  type(inclusion_map)     :: inclusions

  which = ''
  call get_command_argument(1, which)

  call a % declare()
  call sets % bind(a, counted_set_representation(3))
  call b % declare()
  call sets % bind(b, counted_set_representation(4))
  call c % declare()
  call sets % bind(c, counted_set_representation(2))

  r = stored_relation('r', [a, b, c], &
       & reshape([1,1,1,  2,2,2], [3, 2]), sets)

  select case (trim(which))

  case ('slot')
     out1 = restrict_slot(r, 4, b, sets, inclusions)

  case ('embed')
     call other % declare()
     call sets % bind(other, counted_set_representation(4))
     call foreign % declare()
     call sets       % bind(foreign, listed_set_representation([2]))
     call inclusions % include_in(foreign, other)
     out1 = restrict_slot(r, 2, foreign, sets, inclusions)

  case ('none')
     out1 = project_slots(r, [integer ::], sets)

  case ('range')
     out1 = project_slots(r, [1, 5], sets)

  case ('repeat')
     out1 = project_slots(r, [1, 1], sets)

  case ('binary')
     r_bc = stored_relation('bc', [b, c], reshape([1, 1], [2, 1]), sets)
     out2 = compose_binary(r, r_bc, sets)

  case ('middle')
     r_ab = stored_relation('ab', [a, b], reshape([1, 1], [2, 1]), sets)
     fat  = stored_relation('cb', [c, b], reshape([1, 1], [2, 1]), sets)
     out2 = compose_binary(r_ab, fat, sets)

  case default
     write(*,'(1x,a)') &
          & "usage: refusal slot|embed|none|range|repeat|binary|middle"
     error stop 'no case chosen'

  end select

  ! Reaching this line is the failure.
  write(*,'(1x,a,a)') "REACHED PAST THE REFUSAL: ", trim(which)

end program algebra_refusal
