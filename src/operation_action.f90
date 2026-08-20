!=====================================================================!
! The operation prime: the verb within a graph, (graph, fields) ->
! field. Three symbols. name says what it is; domain says which
! member set the answer lives on and how many entries it has; apply
! does the work.
!
! apply writes its result into the output argument and never adds to
! what was there. The argument is intent(inout) only so a caller
! already holding a buffer of the right shape can lend it and save
! an allocation; lending changes the cost of the call, not its
! meaning.
!
! A concrete operation is handed the fields it reads when it is
! constructed - a coefficient, a measure, a geometry field arrives as
! an argument the compiler checks - so apply fetches nothing by name.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_action

  use view_directed , only : directed_graph
  use graph_fractal       , only : set_graph => graph
  use field_calculus, only : graph_field

  implicit none

  private

  public :: graph_operation

  type, abstract :: graph_operation

   contains

     procedure(operation_name_interface)  , deferred :: name
     procedure(operation_domain_interface), deferred :: domain
     procedure(operation_apply_interface) , deferred :: apply

  end type graph_operation

  abstract interface

     pure function operation_name_interface(this) result(name)
       import :: graph_operation
       class(graph_operation), intent(in) :: this
       character(len=:), allocatable :: name
     end function operation_name_interface

     !---------------------------------------------------------------!
     ! Where the answer lives: WHICH set, and HOW MANY entries it
     ! has. The count travels beside the identity because every
     ! caller wants exactly those two things - to check the domain
     ! matches, and to size a field.
     !---------------------------------------------------------------!

     subroutine operation_domain_interface(this, input_graph, domain, &
          & nentries)
       import :: graph_operation, directed_graph, set_graph
       class(graph_operation), intent(in)  :: this
       class(directed_graph)          , intent(in)  :: input_graph
       type(set_graph)       , intent(out) :: domain
       integer               , intent(out) :: nentries
     end subroutine operation_domain_interface

     subroutine operation_apply_interface(this, input_graph, input_data, output)
       import :: graph_operation, directed_graph, graph_field
       class(graph_operation), intent(in) :: this
       class(directed_graph), intent(in) :: input_graph
       class(graph_field), intent(in), optional :: input_data(:)
       class(graph_field), allocatable, intent(inout) :: output
     end subroutine operation_apply_interface

  end interface

end module operation_action
