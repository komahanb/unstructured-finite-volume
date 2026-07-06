!=====================================================================!
! The common ancestor of every iterative process in the framework: the
! linear solver, the nonlinear solver and the time integrator are one
! concept - advance a state until a termination criterion is met. The
! deferred contract is march:
!
!     call marcher % march(system, s, mode)
!
! where the system supplies residuals and products, s is the state the
! caller owns (class(state): the concrete knows how a correction moves
! it), and mode selects the forward or the adjoint direction - one
! march plus a mode, never a second procedure.
!
! Each family implements march and keeps its context entry point as a
! thin provided wrapper that delegates to it:
!
!   family            | march implementation      | wrapper
!   ------------------+---------------------------+--------------------
!   linear solver     | converge (residual        | solve(system, x,
!                     | minimization around the   | mode)
!                     | deferred iterate)         |
!   nonlinear solver  | drive the deferred        | (bdf calls the
!                     | solve(system, coeff, U)   | per-step solve
!                     | at the state's            | directly)
!                     | linearization             |
!   integrator        | advance the state window  | integrate(mode)
!                     | along the step chain,     |
!                     | recording the trajectory  |
!
! What also lives here: the shared attributes (tolerance, iteration
! budget, verbosity) previously duplicated across the families.
!
! In graph terms: a marcher traverses a chain - the iterate sequence in
! solver iterations, the step sequence in physical time. The chain is a
! graph (class_chain), the states are its vertex payloads, and the
! discrete adjoint traverses the chain in reverse
! (interface_graph % accumulate_adjoint).
!
! Author : Komahan Boopathy
!=====================================================================!

module interface_marcher

  use iso_fortran_env    , only : dp => REAL64
  use interface_assembler, only : assembler
  use interface_state    , only : state

  implicit none

  private
  public :: marcher

  !===================================================================!
  ! The abstract marcher
  !===================================================================!

  type, abstract :: marcher

     ! termination tolerance, iteration budget, verbosity. defaults
     ! reproduce the nonlinear family's previous values; the linear
     ! solvers and the integrator always set their own.
     real(dp) :: max_tol     = 1.0d-12
     integer  :: max_it      = 25
     integer  :: print_level = 0

   contains

     ! the contract: advance the state until termination
     procedure(march_interface), deferred :: march

  end type marcher

  !===================================================================!
  ! Deferred interface
  !===================================================================!

  abstract interface

     impure subroutine march_interface(this, system, s, mode)
       import :: marcher, assembler, state
       class(marcher)  , intent(inout) :: this
       class(assembler), intent(inout) :: system
       class(state)    , intent(inout) :: s
       integer         , intent(in), optional :: mode  ! FORWARD (default) / REVERSE
     end subroutine march_interface

  end interface

contains

end module interface_marcher
