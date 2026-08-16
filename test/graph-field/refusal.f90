!=====================================================================!
! The field refusals, each EXPECTED TO DIE for its stated reason:
!      shape        a value vector that does not fill its domain
!      unsigned     a field on a domain that never signed
!=====================================================================!
program field_refusal
  use iso_fortran_env  , only : dp => REAL64
  use fractal_graph           , only : graph
  use graph_set_representation, only : counted_set_representation
  use graph_set_map           , only : set_map
  use class_graph_field, only : field
  implicit none
  type(graph)       :: cells, raw
  type(set_map)     :: sets
  type(field)       :: q
  character(len=32) :: which
  which=''; call get_command_argument(1, which)
  call cells % declare()
  call sets % bind(cells, counted_set_representation(4))
  select case (trim(which))
  case ('ishape')
     q = field('q', cells, 4)
     call q % set_integer_vector([1, 2])
  case ('rshape')
     q = field('q', cells, 4)
     call q % set_real_vector([1.0_dp, 2.0_dp])
  case ('cshape')
     q = field('q', cells, 4)
     call q % set_complex_vector([(1.0_dp, 0.0_dp)])
  case ('lshape')
     q = field('q', cells, 4)
     call q % set_logical_vector([.true., .false., .true.])
  case ('sshape')
     q = field('q', cells, 4)
     call q % set_character_vector(['a', 'b'])
  case ('unsigned')
     q = field('q', raw, 4)
  case default
     error stop 'no case chosen'
  end select
  write(*,'(1x,a,a)') "REACHED PAST THE REFUSAL: ", trim(which)
end program field_refusal
