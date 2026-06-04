#include "scalar.fpp"

!=====================================================================!
! State-energy functional  J = 1/2 sum_i u_i^2 = 1/2 ||u||^2. A simple
! differentiable objective for verifying the adjoint: its state
! derivative df/du = u is a non-trivial vector, and it has no explicit
! design dependence (df/dx = 0, inherited default). Everything is
! vector-level, so it stays valid as flat arrays become distributed.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_state_energy

  use iso_fortran_env    , only : dp => REAL64
  use interface_assembler, only : assembler
  use interface_function , only : functional

  implicit none

  private
  public :: state_energy

  type, extends(functional) :: state_energy

   contains

     procedure :: eval
     procedure :: add_dfdu

  end type state_energy

  interface state_energy
     module procedure create
  end interface state_energy

contains

  !===================================================================!
  ! Construct the state-energy functional
  !===================================================================!

  type(state_energy) function create() result(this)

    this % description = "state energy J = 1/2 ||u||^2"

  end function create

  !===================================================================!
  ! Value  J = 1/2 ||u||^2,  u = S(:,1)
  !===================================================================!

  subroutine eval(this, system, fval)

    class(state_energy), intent(in)  :: this
    class(assembler)   , intent(in)  :: system
    type(scalar)       , intent(out) :: fval

    fval = 0.5_dp*sum(system % S(:,1)**2)

  end subroutine eval

  !===================================================================!
  ! State derivative  df/du = u
  !===================================================================!

  subroutine add_dfdu(this, system, dfdu)

    class(state_energy), intent(in)    :: this
    class(assembler)   , intent(in)    :: system
    type(scalar)       , intent(inout) :: dfdu(:)

    dfdu = dfdu + system % S(:,1)

  end subroutine add_dfdu

end module class_state_energy
