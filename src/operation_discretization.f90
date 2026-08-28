!=====================================================================!
! LEVEL 6 OF THE TOWER: DISCRETIZATION CALCULUS
!
! The abstract contract by which a continuous statement becomes
! discrete algebra: a discretization is an operation built from an
! operation by binding it to a graph's arithmetic. stencil and scheme
! discretize on the dependent and independent axes; the tangent of
! any statement is the derived operation in operation_linearization;
! chain_rule composes partial actions to any degree.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_discretization

  use operation_action, only : operation

  implicit none

  private
  public :: discretization

  !===================================================================!
  ! DISCRETIZATION. The common marker for operators created by
  ! binding a continuous statement to discrete arithmetic.
  !===================================================================!

  type, abstract, extends(operation) :: discretization

  end type discretization

end module operation_discretization
