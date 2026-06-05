!=====================================================================!
! Common abstract base for every algebraic-equation solver. Both a linear
! solver (interface_linear_solver) and a nonlinear solver
! (interface_nonlinear_solver) extend this.
!
! It carries no deferred procedure of its own on purpose: a linear and a
! nonlinear solve have different contracts (solve(this, x) versus
! solve(this, system, coeff, U)), so the specific solve interface is
! declared by each child. This type is the shared ancestor that lets the
! rest of the code refer to "some algebraic solver".
!
! Author : Komahan Boopathy
!=====================================================================!

module interface_algebraic_solver

  implicit none

  private
  public :: algebraic_solver

  !===================================================================!
  ! Any solver for the algebraic set of equations resulting after
  ! discretization (solved until R = 0).
  !===================================================================!

  type, abstract :: algebraic_solver
  end type algebraic_solver

contains

end module interface_algebraic_solver
