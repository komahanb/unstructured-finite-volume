!=====================================================================!
! LEVEL 6 OF THE TOWER: DISCRETIZATION CALCULUS
!
! The three abstract contracts by which a continuous statement
! becomes discrete algebra and exposes its derivatives:
!
!      discretization_operator   an operation built from an
!                                operation by binding it to a
!                                graph's arithmetic; owes its
!                                dependency pattern as a graph
!      linearization_operator    an operation built from an
!                                operation by freezing a state:
!                                the tangent J v at that state,
!                                behind the ordinary operation
!                                interface; freeze moves the state
!      differentiable_operation  an operation that computes exact
!                                partial actions - multilinear
!                                contractions of its derivatives
!                                against direction fields - up to
!                                a declared maximum order
!
! Concrete members: stencil_operator and step_operator discretize
! on the dependent and independent axes; difference_linearization
! and exact_linearization are the two tangents, selected by
! tangent_of; chain_rule composes partial actions to any degree.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_discretization

  use iso_fortran_env     , only : dp => REAL64
  use view_directed , only : directed_graph
  use field_calculus, only : graph_field
  use operation_action, only : graph_operation

  implicit none

  private
  public :: discretization_operator
  public :: linearization_operator
  public :: differentiable_operation

  !===================================================================!
  ! DISCRETIZATION_OPERATOR. Owes by contract its dependency
  ! pattern, as a graph, on the axis the concrete type represents:
  ! a stencil_operator's pattern is on the dependent variable
  ! (which unknown feeds which), a step_operator's on the
  ! independent variable (which instants the residual reads).
  !===================================================================!

  type, abstract, extends(graph_operation) :: discretization_operator

   contains

     procedure(discretization_pattern_interface), deferred :: dependencies

  end type discretization_operator

  !===================================================================!
  ! LINEARIZATION_OPERATOR. The tangent of S at a frozen state,
  ! behind the operation interface, so a minimizer sees an
  ! ordinary linear operation. One deferred routine beyond the
  ! interface: freeze, which moves the state (and may cache the
  ! base residual) between a governor's steps.
  !===================================================================!

  type, abstract, extends(graph_operation) :: linearization_operator

   contains

     procedure(linearization_freeze_interface), deferred :: freeze

  end type linearization_operator

  !===================================================================!
  ! DIFFERENTIABLE_OPERATION. An operation that, beyond apply,
  ! computes exact partial actions: multilinear contractions of
  ! its derivatives against direction fields, one input slot named
  ! per derivative, up to a declared maximum order. No derivative
  ! tensor is stored; an action is computed, contracted, and
  ! discarded. The exact linearization and the chain rule are
  ! built on this contract.
  !===================================================================!

  type, abstract, extends(graph_operation) :: differentiable_operation

   contains

     procedure(differentiable_max_degree_interface)    , deferred :: max_degree
     procedure(differentiable_partial_action_interface), deferred :: partial_action

  end type differentiable_operation

  abstract interface

     !===============================================================!
     ! The dependency pattern a discretization owes: a graph, on
     ! the concrete type's own axis.
     !===============================================================!

     subroutine discretization_pattern_interface(this, pattern)
       import :: discretization_operator, directed_graph
       class(discretization_operator), intent(in) :: this
       class(directed_graph), allocatable, intent(out)     :: pattern
     end subroutine discretization_pattern_interface

     !===============================================================!
     ! Move the frozen state; base optionally carries the residual
     ! at that state for members that cache it.
     !===============================================================!

     subroutine linearization_freeze_interface(this, at, base)
       import :: linearization_operator, dp
       class(linearization_operator), intent(inout) :: this
       real(dp), intent(in)           :: at(:)
       real(dp), intent(in), optional :: base(:)
     end subroutine linearization_freeze_interface

     !===============================================================!
     ! The two deferred routines: max_degree reports how deep the
     ! exact calculus goes; partial_action computes one mixed
     ! partial - the statement differentiated once per entry of
     ! slots(:), contracted against the matching direction field,
     ! returned on the statement's own domain.
     !===============================================================!

     pure function differentiable_max_degree_interface(this) result(degree)
       import :: differentiable_operation
       class(differentiable_operation), intent(in) :: this
       integer :: degree
     end function differentiable_max_degree_interface

     subroutine differentiable_partial_action_interface(this, input_graph, &
          & input_data, slots, directions, output)
       import :: differentiable_operation, directed_graph, graph_field
       class(differentiable_operation), intent(in)    :: this
       class(directed_graph)          , intent(in)    :: input_graph
       class(graph_field)             , intent(in)    :: input_data(:)
       integer                        , intent(in)    :: slots(:)
       class(graph_field)             , intent(in)    :: directions(:)
       class(graph_field), allocatable, intent(inout) :: output
     end subroutine differentiable_partial_action_interface

  end interface

end module operation_discretization
