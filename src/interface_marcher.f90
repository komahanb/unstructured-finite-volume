!=====================================================================!
! The common ancestor of every iterative process in the framework: the
! linear solver, the nonlinear solver and the time integrator are one
! concept - advance a deferred one-step map along a chain until a
! termination criterion is met. Each family binds its own context name
! to the common operation:
!
!   common (marcher)  | linear solver | newton            | integrator
!   ------------------+---------------+-------------------+-------------
!   march             | solve         | solve             | integrate
!   step (deferred    | iterate       | linearize + inner | bdf step
!    per family)      |               | solve             |
!   termination       | residual tol  | residual tol      | t >= t_final
!   budget (max_it)   | max_it        | max_it            | num_steps
!   reverse march     | adjoint       | steady adjoint    | backward
!                     | accumulation  |                   | adjoint
!
! What lives here: the shared attributes (tolerance, iteration budget,
! verbosity) that were previously duplicated across the solver families.
! What stays per family, deliberately: the step signature - linear steps
! move plain vectors, newton steps move multi-order states with
! coefficients, time steps move trajectory windows. Fortran has no
! generics, and a common state type would be heavy machinery for a
! ten-line loop; each family declares its own step contract with the
! shared name and shape, and its own short march bound to its context
! name through a generic (solve => march, integrate => march).
!
! In graph terms: a marcher traverses a chain - the iterate sequence in
! solver iterations, the step sequence in physical time. The chain is a
! graph (class_chain); the marcher advances along it, and the discrete
! adjoint traverses it in reverse (interface_graph % accumulate_adjoint).
!
! Author : Komahan Boopathy
!=====================================================================!

module interface_marcher

  use iso_fortran_env, only : dp => REAL64

  implicit none

  private
  public :: marcher

  !===================================================================!
  ! The abstract marcher. It carries no deferred procedure of its own
  ! on purpose: the march signatures differ per family (solve(this,
  ! system, x, mode) versus solve(this, system, coeff, U) versus
  ! integrate(this, mode)), so each family declares its own contract.
  ! This type is the shared ancestor that lets the rest of the code
  ! refer to "some marcher", and the one home of the shared attributes.
  !===================================================================!

  type, abstract :: marcher

     ! termination tolerance, iteration budget, verbosity. defaults
     ! reproduce the nonlinear family's previous values; the linear
     ! solvers and the integrator always set their own.
     real(dp) :: max_tol     = 1.0d-12
     integer  :: max_it      = 25
     integer  :: print_level = 0

  end type marcher

contains

end module interface_marcher
