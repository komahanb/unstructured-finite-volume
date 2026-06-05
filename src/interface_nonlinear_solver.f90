#include "scalar.fpp"

!=====================================================================!
! Abstract interface for nonlinear solvers. A nonlinear solver drives a
! system's residual R(U) = 0 by repeatedly linearising and solving with
! an (inner) linear solver. It extends the common algebraic_solver base
! and carries the linear solver used for those linearised steps.
!
! The deferred solve takes the assembler (the system), the integrator's
! linearisation coefficients coeff, and the multi-order state U it
! updates in place - matching the condensed Newton step the BDF marcher
! and the adjoint forward solve drive.
!
! Author : Komahan Boopathy
!=====================================================================!

module interface_nonlinear_solver

  use iso_fortran_env           , only : dp => REAL64
  use interface_algebraic_solver, only : algebraic_solver
  use interface_linear_solver   , only : linear_solver
  use interface_assembler       , only : assembler

  implicit none

  private
  public :: nonlinear_solver

  !===================================================================!
  ! Abstract nonlinear solver datatype
  !===================================================================!

  type, abstract, extends(algebraic_solver) :: nonlinear_solver

     ! the inner linear solver for each linearised step (optional; a
     ! matrix-free inner iteration may be used instead)
     class(linear_solver), allocatable :: linear_solver

     integer  :: print_level = 0
     real(dp) :: max_tol     = 1.0d-12
     integer  :: max_it      = 25

   contains

     procedure(nonlinear_solve_interface), deferred :: solve

  end type nonlinear_solver

  !===================================================================!
  ! Deferred solve: drive R(U) = 0 for the system at the given
  ! linearisation coefficients, updating the multi-order state U.
  !===================================================================!

  abstract interface

     subroutine nonlinear_solve_interface(this, system, coeff, U)
       import nonlinear_solver
       import assembler
       class(nonlinear_solver), intent(inout) :: this   ! configures its held linear solver
       class(assembler)       , intent(inout) :: system
       type(scalar)           , intent(in)    :: coeff(:)
       type(scalar)           , intent(inout) :: U(:,:)
     end subroutine nonlinear_solve_interface

  end interface

contains

end module interface_nonlinear_solver
