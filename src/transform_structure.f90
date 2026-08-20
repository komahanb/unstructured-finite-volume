!=====================================================================!
! The transform prime: the verb between graphs, graph -> graph with
! data carried along. Two symbols at this level, both admissibility
! questions: may this transform act on that graph, and on that data
! riding on it. The concretions - partition, assemble, coarsen,
! refine - are each judged by a round-trip law:
!
!      exact        assemble(partition(G)) = G      both ways
!      one-sided    coarsen(refine(G)) = G          one way only
!
! The interfaces are impure for a historical language constraint that
! no longer binds (F2018 C1594 on copying a pointer-carrying set
! graph inside a pure subprogram); an impure interface permits a pure
! override, and every implementation is pure.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module transform_structure

  use view_directed , only : directed_graph
  use field_calculus, only : field

  implicit none

  private

  public :: transform

  type, abstract :: transform

   contains

     procedure(transform_on_graph_interface), deferred :: defined_on_graph
     procedure(transform_on_data_interface) , deferred :: defined_on_data

  end type transform

  abstract interface

     logical function transform_on_graph_interface(this, input_graph)
       import :: transform, directed_graph
       class(transform), intent(in) :: this
       class(directed_graph), intent(in) :: input_graph
     end function transform_on_graph_interface

     logical function transform_on_data_interface(this, input_graph, input_data)
       import :: transform, directed_graph, field
       class(transform), intent(in) :: this
       class(directed_graph), intent(in) :: input_graph
       class(field), intent(in) :: input_data
     end function transform_on_data_interface

  end interface

end module transform_structure
