!=====================================================================!
! Interface for linear solvers to implement.
!=====================================================================!

module interface_linear_solver

  use iso_fortran_env , only : dp => REAL64

  implicit none
  
  ! Expose only the linear solver interface
  private
  public :: linear_solver

  !===================================================================!
  ! Linear solver datatype
  !===================================================================!
  
  type, abstract :: linear_solver

     real(dp) :: max_tol
     integer  :: max_it

   contains

     ! type bound procedures
     procedure(solve_interface), deferred :: solve

  end type linear_solver
  
  !===================================================================!
  ! Interface for multiple constructors
  !===================================================================!
  
  interface
     subroutine solve_interface(this, x)
       import linear_solver
       import dp
       class(linear_solver)  , intent(in)  :: this
       real(dp), allocatable , intent(out) :: x(:)
     end subroutine solve_interface
  end interface
  
contains

end module interface_linear_solver
