#include "scalar.fpp"

!=====================================================================!
! A scalar function of interest  J = f(u, x)  evaluated on the state u
! held by an assembler (and, through it, the design variables x). The
! adjoint needs its partials: df/du (the adjoint right-hand side) and
! df/dx (an explicit design dependence, usually zero).
!
! Concrete functionals (state energy, mean field, ...) extend this. The
! state is read vector-opaquely from system % S(:,1); a functional that
! genuinely owns its layout may index, like the physics does.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module interface_function

  use iso_fortran_env    , only : dp => REAL64
  use interface_assembler, only : assembler

  implicit none

  private
  public :: functional

  !-------------------------------------------------------------------!
  ! Abstract function of interest
  !-------------------------------------------------------------------!

  type, abstract :: functional

     character(len=:), allocatable :: description

   contains

     ! Deferred: the value and its state derivative
     procedure(eval_interface)   , deferred :: eval
     procedure(add_dfdu_interface), deferred :: add_dfdu

     ! Provided: explicit design derivative (default: none)
     procedure :: add_dfdx

  end type functional

  !-------------------------------------------------------------------!
  ! Interfaces for the deferred procedures
  !-------------------------------------------------------------------!

  abstract interface

     !================================================================!
     ! Scalar value  J = f(u, x)  at the assembler's current state
     !================================================================!

     impure subroutine eval_interface(this, system, fval)

       import :: functional, assembler

       class(functional), intent(in)  :: this
       class(assembler) , intent(in)  :: system
       type(scalar)     , intent(out) :: fval

     end subroutine eval_interface

     !================================================================!
     ! Accumulate the state derivative  dfdu += df/du  (adjoint rhs)
     !================================================================!

     impure subroutine add_dfdu_interface(this, system, dfdu)

       import :: functional, assembler

       class(functional), intent(in)    :: this
       class(assembler) , intent(in)    :: system
       type(scalar)     , intent(inout) :: dfdu(:)

     end subroutine add_dfdu_interface

  end interface

contains

  !===================================================================!
  ! Explicit design derivative  dfdx += df/dx. Default: the functional
  ! has no direct design dependence (it sees x only through u), so this
  ! adds nothing. Functionals with an explicit x term override it.
  !===================================================================!

  pure subroutine add_dfdx(this, system, dfdx)

    class(functional), intent(in)    :: this
    class(assembler) , intent(in)    :: system
    real(dp)         , intent(inout) :: dfdx(:)

  end subroutine add_dfdx

end module interface_function
