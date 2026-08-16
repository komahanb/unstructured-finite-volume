!=====================================================================!
! What a coarray ultimate component does prohibit is the API.
!
! An object of such a type must be a nonpointer nonallocatable scalar,
! so it can be neither the pointer a fixture hands back nor an element
! of an array of bindings. The mechanism removes the storage shape and
! leaves the assignment.
!=====================================================================!

module coarray_pointer_m

  implicit none

  type :: bind_t
     integer, allocatable :: guard[:]
     integer              :: n = 0
  end type bind_t

end module coarray_pointer_m

program coarray_pointer

  use coarray_pointer_m

  implicit none

  type(bind_t), pointer :: p => null()

  print *, associated(p)

end program coarray_pointer
