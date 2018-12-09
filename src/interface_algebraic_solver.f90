!=====================================================================!
! Interface for linear and nonlinear solvers to implement.
! 
! Author : Komahan Boopathy
!=====================================================================!

module interface_algebraic_solver

  use iso_fortran_env , only : dp => REAL64
  use class_assembler , only : assembler

  implicit none
  
  private
  public :: algebraic_solver

  !===================================================================!
  ! Any solver for algebraic set of equations resulting after
  ! discretization
  !===================================================================!
  
  type, abstract :: algebraic_solver

     ! Governing System to be solved until R = 0
     class(assembler), allocatable :: system

   contains

     ! Deferred behavior for system that is solved until get_residual
     ! = 0
     procedure(solve_interface), deferred :: solve

  end type algebraic_solver
  
  !===================================================================!
  ! Interface for multiple constructors
  !===================================================================!
  
  interface
     subroutine solve_interface(this, x)
       import algebraic_solver
       import dp
       class(algebraic_solver) , intent(in)  :: this
       real(dp), allocatable   , intent(out) :: x(:)
     end subroutine solve_interface
  end interface
  
contains

end module interface_algebraic_solver
