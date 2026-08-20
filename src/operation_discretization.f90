!=====================================================================!
! LEVEL 6 OF THE TOWER: DISCRETIZATION CALCULUS
!
! The two abstract contracts by which a continuous statement
! becomes discrete algebra:
!
!      discretization   an operation built from an operation by
!                       binding it to a graph's arithmetic; owes
!                       its dependency pattern as a graph
!      linearization    an operation built from an operation by
!                       freezing a state: the tangent J v at that
!                       state, behind the ordinary operation
!                       interface; freeze moves the state
!
! Concrete members: stencil and scheme discretize on the dependent
! and independent axes; difference_linearization and
! exact_linearization are the two tangents, selected by tangent_of
! from the operation's max_degree; chain_rule composes partial
! actions to any degree.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_discretization

  use iso_fortran_env     , only : dp => REAL64
  use view_directed , only : directed_graph
  use operation_action, only : operation

  implicit none

  private
  public :: discretization
  public :: linearization

  !===================================================================!
  ! DISCRETIZATION. Owes by contract its dependency
  ! pattern, as a graph, on the axis the concrete type represents:
  ! a stencil's pattern is on the dependent variable
  ! (which unknown feeds which), a scheme's on the
  ! independent variable (which instants the residual reads).
  !===================================================================!

  type, abstract, extends(operation) :: discretization

   contains

     procedure(discretization_pattern_interface), deferred :: dependencies

  end type discretization

  !===================================================================!
  ! LINEARIZATION. The tangent of S at a frozen state,
  ! behind the operation interface, so a minimizer sees an
  ! ordinary linear operation. One deferred routine beyond the
  ! interface: freeze, which moves the state (and may cache the
  ! base residual) between a governor's steps.
  !===================================================================!

  type, abstract, extends(operation) :: linearization

   contains

     procedure(linearization_freeze_interface), deferred :: freeze

  end type linearization

  abstract interface

     !===============================================================!
     ! The dependency pattern a discretization owes: a graph, on
     ! the concrete type's own axis.
     !===============================================================!

     subroutine discretization_pattern_interface(this, pattern)
       import :: discretization, directed_graph
       class(discretization), intent(in) :: this
       class(directed_graph), allocatable, intent(out)     :: pattern
     end subroutine discretization_pattern_interface

     !===============================================================!
     ! Move the frozen state; base optionally carries the residual
     ! at that state for members that cache it.
     !===============================================================!

     subroutine linearization_freeze_interface(this, at, base)
       import :: linearization, dp
       class(linearization), intent(inout) :: this
       real(dp), intent(in)           :: at(:)
       real(dp), intent(in), optional :: base(:)
     end subroutine linearization_freeze_interface

  end interface

end module operation_discretization
