#include "scalar.fpp"

!=====================================================================!
! Abstract interface for nonlinear solvers. A nonlinear solver drives a
! system's residual R(U) = 0 by repeatedly linearising and solving with
! an (inner) linear solver. It extends the common marcher base and
! carries the linear solver used for those linearised steps; its march
! is bound to the context name solve through a generic.
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
  use interface_marcher         , only : marcher
  use interface_state           , only : state
  use class_differential_state  , only : differential_state
  use module_solve_mode         , only : REVERSE
  use interface_linear_solver   , only : linear_solver
  use interface_assembler       , only : assembler

  implicit none

  private
  public :: nonlinear_solver

  !===================================================================!
  ! Abstract nonlinear solver datatype
  !===================================================================!

  type, abstract, extends(marcher) :: nonlinear_solver

     ! the inner linear solver for each linearised step (optional; a
     ! matrix-free inner iteration may be used instead)
     class(linear_solver), allocatable :: linear_solver

     ! max_tol (1.0d-12), max_it (25) and print_level are inherited
     ! from marcher with the same defaults

   contains

     ! the family's per-step deferred entry (bdf drives this directly)
     procedure(nonlinear_solve_interface), deferred :: solve

     ! the marcher contract, provided once for the family: drive the
     ! deferred solve at the state's linearization
     procedure :: march

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

  !===================================================================!
  ! The marcher contract for the nonlinear family: drive the deferred
  ! solve at the state's linearization. Requires a differential_state
  ! (the linearization coefficients live on the state) - a stated
  ! precondition with a clear error. REVERSE is refused until the
  ! functional supplies the adjoint right-hand side.
  !===================================================================!

  impure subroutine march(this, system, s, mode)

    class(nonlinear_solver), intent(inout) :: this
    class(assembler)       , intent(inout) :: system
    class(state)           , intent(inout) :: s
    integer                , intent(in), optional :: mode

    if (present(mode)) then
       if (mode .eq. REVERSE) then
          error stop "nonlinear_solver % march: REVERSE requires the functional's " // &
               & "adjoint right-hand side - use the adjoint driver"
       end if
    end if

    select type (s)
    type is (differential_state)
       call this % solve(system, s % coeff, s % U)
    class default
       error stop "nonlinear_solver % march: requires a differential_state"
    end select

  end subroutine march

end module interface_nonlinear_solver
