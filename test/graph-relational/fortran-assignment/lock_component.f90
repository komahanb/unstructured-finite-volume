!=====================================================================!
! A LOCK_TYPE component rejects the declaration, not the assignment.
!
! It forces every variable of the type to be a coarray, which forbids
! the pointer component that a fixture and a caller need. And once the
! codimension is paid for, assignment is still accepted - so the
! mechanism removes the API and leaves the operation.
!=====================================================================!

module lock_component_m

  use iso_fortran_env, only : lock_type

  implicit none

  type :: bind_t
     type(lock_type) :: guard
     integer         :: n = 0
  end type bind_t

  type :: holder
     type(bind_t), pointer :: object => null()
  end type holder

end module lock_component_m
