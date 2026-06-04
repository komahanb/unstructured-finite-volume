!=====================================================================!
! Abstract preconditioner: applies an approximate inverse z = M^-1 r.
! A linear solver (e.g. CG) calls apply once per iteration; concrete
! preconditioners (algebraic multigrid, jacobi, ...) extend this.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module interface_preconditioner

  use iso_fortran_env, only : dp => REAL64

  implicit none

  private
  public :: preconditioner

  type, abstract :: preconditioner
   contains
     procedure(apply_interface), deferred :: apply
  end type preconditioner

  abstract interface

     !================================================================!
     ! z = M^-1 r  (the approximate-inverse action)
     !================================================================!

     subroutine apply_interface(this, r, z)
       import :: preconditioner, dp
       class(preconditioner), intent(in)  :: this
       real(dp)             , intent(in)  :: r(:)
       real(dp)             , intent(out) :: z(:)
     end subroutine apply_interface

  end interface

end module interface_preconditioner
