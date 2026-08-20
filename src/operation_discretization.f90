!=====================================================================!
! LEVEL 6 OF THE TOWER: DISCRETIZATION CALCULUS
!
! The abstract contract by which a continuous statement becomes
! discrete algebra: a discretization is an operation built from an
! operation by binding it to a graph's arithmetic, and it owes its
! dependency pattern as a graph. stencil and scheme discretize on the
! dependent and independent axes; the tangent of any statement is the
! derived operation in operation_linearization; chain_rule composes
! partial actions to any degree.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_discretization

  use view_directed , only : directed_graph
  use operation_action, only : operation

  implicit none

  private
  public :: discretization

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

  end interface

end module operation_discretization
