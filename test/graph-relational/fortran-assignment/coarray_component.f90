!=====================================================================!
! A coarray ultimate component does not prohibit assignment.
!
! The assignment is not diagnosed at all, at -Wall -pedantic. What the
! component does prohibit is the pointer form the callers need: see
! coarray_pointer.
!=====================================================================!

module coarray_component_m

  implicit none

  type :: bind_t
     integer, allocatable :: guard[:]
     integer              :: n = 0
  end type bind_t

end module coarray_component_m

program coarray_component

  use coarray_component_m

  implicit none

  type(bind_t) :: a, b

  a = b

  print *, '   the assignment is not diagnosed, at -Wall -pedantic'

end program coarray_component
