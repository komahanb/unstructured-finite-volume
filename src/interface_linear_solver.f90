!=====================================================================!
! Interface for linear solvers, plus the preconditioner contract they
! consume. A linear solver solves A x = b for x; a preconditioner applies
! an approximate inverse z = M^-1 r once per iteration. Both abstract
! types live here so a solver and its preconditioner share one home.
!
! linear_solver extends the common algebraic_solver base.
!
! Author : Komahan Boopathy
!=====================================================================!

module interface_linear_solver

  use iso_fortran_env           , only : dp => REAL64
  use interface_algebraic_solver, only : algebraic_solver

  implicit none

  private
  public :: linear_solver
  public :: preconditioner

  !===================================================================!
  ! Linear solver datatype
  !===================================================================!

  type, abstract, extends(algebraic_solver) :: linear_solver

     real(dp) :: max_tol
     integer  :: max_it

   contains

     ! type bound procedures
     procedure(solve_interface), deferred :: solve

  end type linear_solver

  !===================================================================!
  ! Abstract preconditioner: applies an approximate inverse z = M^-1 r.
  ! A linear solver (e.g. CG) calls apply once per iteration; concrete
  ! preconditioners (algebraic multigrid, jacobi, ...) extend this.
  !===================================================================!

  type, abstract :: preconditioner
   contains
     procedure(apply_interface), deferred :: apply
  end type preconditioner

  !===================================================================!
  ! Deferred interfaces
  !===================================================================!

  interface

     subroutine solve_interface(this, x)
       import linear_solver
       import dp
       class(linear_solver)  , intent(in)  :: this
       real(dp), allocatable , intent(out) :: x(:)
     end subroutine solve_interface

     ! z = M^-1 r  (the approximate-inverse action)
     subroutine apply_interface(this, r, z)
       import preconditioner
       import dp
       class(preconditioner), intent(in)  :: this
       real(dp)             , intent(in)  :: r(:)
       real(dp)             , intent(out) :: z(:)
     end subroutine apply_interface

  end interface

contains

end module interface_linear_solver
