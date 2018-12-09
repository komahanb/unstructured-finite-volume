!=====================================================================!
! Interface for linear and nonlinear solvers to implement.
! 
! Author : Komahan Boopathy
!=====================================================================!

module interface_algebraic_solver

  use iso_fortran_env     , only : dp => REAL64
  use interface_assembler , only : assembler

  implicit none
  
  private
  public :: algebraic_solver

  !===================================================================!
  ! Any solver for algebraic set of equations resulting after
  ! discretization
  !===================================================================!
  
  type, abstract :: algebraic_solver

   contains

     ! Deferred behavior for system that is solved until R = 0
     procedure(solve_interface), deferred :: solve

  end type algebraic_solver
  
  !===================================================================!
  ! Solve function to be implemented by inherting types
  !===================================================================!
  
  interface
     subroutine solve_interface(this, system)
       import algebraic_solver
       import assembler
       import dp
       class(algebraic_solver), intent(in)  :: this
       class(assembler)       , intent(in)  :: system
       !real(dp), allocatable  , intent(out) :: x(:) ! maybe 2D for qdot
     end subroutine solve_interface
  end interface
  
contains

end module interface_algebraic_solver
