!=====================================================================!
! Things that remake a graph.
!
! A transform maps one graph to another and carries the values across
! by the same map, so structure and data never drift apart. Four verbs
! are declared here: cut the whole into parts, put the parts back,
! glue cells into blocks, and split cells into children.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module transform_graph_transform

  use structure_graph, only : graph
  use data_graph_field, only : graph_field

  implicit none

  private

  public :: graph_transform
  public :: graph_partitioner
  public :: graph_assembler
  public :: graph_coarsener
  public :: graph_refiner

  !===================================================================!
  ! GRAPH_TRANSFORM. The verb between graphs.
  !
  ! Two symbols at this level, both admissibility questions: may this
  ! transform act on that graph, and on that data riding on it. The
  ! maps themselves are verbs, and verbs are level-1 concretions -
  ! partition, assemble, coarsen, refine - each pair judged by its
  ! round-trip law:
  !
  !      exact        assemble(partition(G)) = G      both ways
  !      one-sided    coarsen(refine(G)) = G          one way only
  !
  ! The central law of the whole tower is a sentence about this role:
  ! split a graph into parts, work on the parts, put the answer back
  ! together.
  !
  !        G'  =  P^-1 ( A ( P ( G ) ) )
  !
  !             G            P(G)           A(P(G))          G'
  !             |             |                |              |
  !             +-----P------>+-------A------->+----P^-1----->+
  !             |             |                |              |
  !        whole graph    the parts      worked-on parts  whole again
  !
  ! P is a transform and nothing else. P^-1 is a transform and
  ! nothing else. A is an operation and does neither.
  !===================================================================!

  type, abstract :: graph_transform

   contains

     procedure(transform_on_graph_interface), deferred :: defined_on_graph

  end type graph_transform

  !===================================================================!
  ! GRAPH_PARTITIONER. P: the whole becomes parts.
  !
  !         o---o---o---o---o---o
  !                     :                 cut where few edges cross,
  !         o---o---o   :   o---o---o     so the parts talk little
  !                part 1     part 2
  !
  ! What comes out is still a graph, and it holds the record of how
  ! it relates back to the whole - the frame the grammar already
  ! knows how to read.
  !===================================================================!

  type, abstract, extends(graph_transform) :: graph_partitioner

   contains

     procedure(partition_graph_interface), deferred :: partition_graph
     procedure(partition_data_interface) , deferred :: partition_data

  end type graph_partitioner

  !===================================================================!
  ! GRAPH_ASSEMBLER. P^-1, and only that: the parts become the whole
  ! again, in whole-graph order. The law:
  !
  !      assemble( partition( G ) )  =  G
  !
  ! Only owned values are collected. A borrowed value is a copy of a
  ! value another part owns; counting both copies violates
  ! conservation. No physics, no matrices, no solver behaviour
  ! belongs in here.
  !===================================================================!

  type, abstract, extends(graph_transform) :: graph_assembler

   contains

     procedure(assemble_graph_interface), deferred :: assemble_graph
     procedure(assemble_data_interface) , deferred :: assemble_data

  end type graph_assembler

  !===================================================================!
  ! GRAPH_COARSENER. C: fewer, larger cells.
  !
  !      o o o o                O   O
  !      o o o o    ------>              a multigrid level, where the
  !      o o o o                O   O    slow, smooth part of the
  !                                      error is cheap to kill
  !===================================================================!

  type, abstract, extends(graph_transform) :: graph_coarsener

   contains

     procedure(coarsen_graph_interface), deferred :: coarsen_graph
     procedure(coarsen_data_interface) , deferred :: coarsen_data

  end type graph_coarsener

  !===================================================================!
  ! GRAPH_REFINER. R: the other way, one cell becomes several. The
  ! pair's law is one-sided:
  !
  !      coarsen( refine( G ) )  =  G
  !
  ! and the other direction loses, because refinement invents detail
  ! that coarsening cannot recover.
  !===================================================================!

  type, abstract, extends(graph_transform) :: graph_refiner

   contains

     procedure(refine_graph_interface), deferred :: refine_graph
     procedure(refine_data_interface) , deferred :: refine_data

  end type graph_refiner

  abstract interface

     pure logical function transform_on_graph_interface(this, input_graph)
       import :: graph_transform, graph
       class(graph_transform), intent(in) :: this
       class(graph), intent(in) :: input_graph
     end function transform_on_graph_interface

     subroutine partition_graph_interface(this, global_graph, part_graph)
       import :: graph_partitioner, graph
       class(graph_partitioner), intent(in) :: this
       class(graph), intent(in) :: global_graph
       class(graph), allocatable, intent(out) :: part_graph
     end subroutine partition_graph_interface

     subroutine partition_data_interface(this, global_graph, global_data, &
          & part_graph, part_data)
       import :: graph_partitioner, graph, graph_field
       class(graph_partitioner), intent(in) :: this
       class(graph), intent(in) :: global_graph
       class(graph_field), intent(in) :: global_data
       class(graph), intent(in) :: part_graph
       class(graph_field), allocatable, intent(out) :: part_data
     end subroutine partition_data_interface

     subroutine assemble_graph_interface(this, part_graph, global_graph)
       import :: graph_assembler, graph
       class(graph_assembler), intent(in) :: this
       class(graph), intent(in) :: part_graph
       class(graph), allocatable, intent(out) :: global_graph
     end subroutine assemble_graph_interface

     subroutine assemble_data_interface(this, part_graph, part_data, &
          & global_graph, global_data)
       import :: graph_assembler, graph, graph_field
       class(graph_assembler), intent(in) :: this
       class(graph), intent(in) :: part_graph
       class(graph_field), intent(in) :: part_data
       class(graph), intent(in) :: global_graph
       class(graph_field), allocatable, intent(out) :: global_data
     end subroutine assemble_data_interface

     subroutine coarsen_graph_interface(this, fine_graph, coarse_graph)
       import :: graph_coarsener, graph
       class(graph_coarsener), intent(in) :: this
       class(graph), intent(in) :: fine_graph
       class(graph), allocatable, intent(out) :: coarse_graph
     end subroutine coarsen_graph_interface

     subroutine coarsen_data_interface(this, fine_graph, fine_data, &
          & coarse_graph, coarse_data)
       import :: graph_coarsener, graph, graph_field
       class(graph_coarsener), intent(in) :: this
       class(graph), intent(in) :: fine_graph
       class(graph_field), intent(in) :: fine_data
       class(graph), intent(in) :: coarse_graph
       class(graph_field), allocatable, intent(out) :: coarse_data
     end subroutine coarsen_data_interface

     subroutine refine_graph_interface(this, coarse_graph, fine_graph)
       import :: graph_refiner, graph
       class(graph_refiner), intent(in) :: this
       class(graph), intent(in) :: coarse_graph
       class(graph), allocatable, intent(out) :: fine_graph
     end subroutine refine_graph_interface

     subroutine refine_data_interface(this, coarse_graph, coarse_data, &
          & fine_graph, fine_data)
       import :: graph_refiner, graph, graph_field
       class(graph_refiner), intent(in) :: this
       class(graph), intent(in) :: coarse_graph
       class(graph_field), intent(in) :: coarse_data
       class(graph), intent(in) :: fine_graph
       class(graph_field), allocatable, intent(out) :: fine_data
     end subroutine refine_data_interface

  end interface

end module transform_graph_transform
