!=====================================================================!
! LEVEL 1 OF THE STRATIFICATION . THE CALCULUS
!
! The first level above the ground. This module answers one question
! and no other: HOW QUANTITIES RELATE on a graph. It holds the eight
! named citizens of graph mathematics, each one a partial concretion
! of a grammar role - one bound parameter per name, and nothing else
! new.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!
!
!                      THE EIGHT CONCRETIONS
!
!    graph -----[edge count = 0]----------> graph_support
!    field -----[domain size = 1]---------> graph_functional
!    transform -[verb = partition]--------> graph_partitioner
!    transform -[verb = assemble]---------> graph_assembler
!    transform -[verb = coarsen]----------> graph_coarsener
!    transform -[verb = refine]-----------> graph_refiner
!
! and, standing BESIDE the operation role rather than under it:
!
!    graph_reduction                        field      -> functional
!    graph_broadcast                        functional -> field
!
! ONE DELIBERATE DEVIATION. The reduction and the broadcast are the
! maps that touch the one-entry domain, so by the letter of the
! tower they would be concretions of graph_operation. They stand
! apart for a stated reason: the reduction's measure argument and
! its staged parallel combine are settled designs that the operation
! signature cannot carry without reopening them. A deviation chosen
! and recorded beats a symmetry restored at the price of a won
! argument.
!
! The census of this level: eight types, fifteen operation symbols.
! With the grammar beneath, the tower stands at twelve types and
! seventy symbols.
!
!=====================================================================!
!
!                       THE ROUND-TRIP LAWS
!
! Every pair on this level is judged by what a round trip preserves:
!
!    exact        assemble( partition( G ) ) = G       both ways
!    one-sided    coarsen( refine( G ) ) = G           one way only
!    transpose    average( copy( c ) ) = c             the pair of
!                 sum( share( c ) ) = c                one-to-many
!
! The laws live in the test suite, because the types cannot state
! them; a law the sorts stopped enforcing is an equation the suite
! must carry.
!
!=====================================================================!

module graph_calculus

  use graph_grammar, only : graph, graph_field, graph_transform

  implicit none

  private

  public :: graph_support
  public :: graph_functional
  public :: graph_reduction
  public :: graph_broadcast
  public :: graph_partitioner
  public :: graph_assembler
  public :: graph_coarsener
  public :: graph_refiner

  public :: GRAPH_SIDE_VERTEX
  public :: GRAPH_SIDE_EDGE

  ! Which side of its host a support references. A finite axis,
  ! absorbed: it rides as an answer, never as a pair of types.
  integer, parameter :: GRAPH_SIDE_VERTEX = 1
  integer, parameter :: GRAPH_SIDE_EDGE   = 2

  !===================================================================!
  ! GRAPH_SUPPORT. The grammar's graph at edge count zero: pure
  ! membership, no incidence.
  !
  !      all_vertices           tagged_edges('wall')
  !      +-------------+        +-------------+
  !      | 1 2 3 4 5 6 |        |  11  14  19 |
  !      +-------------+        +-------------+
  !
  ! Almost everything a membership list is asked, the inherited
  ! grammar already answers: its size is num_vertices, its members
  ! are its global indices, its identity is its id. One question has
  ! no inherited answer, and it is the one symbol added here: which
  ! SIDE of the host graph the members reference - its vertices or
  ! its edges. An index list cannot say this about itself, and the
  ! graph no longer carries data that could.
  !===================================================================!

  type, abstract, extends(graph) :: graph_support

   contains

     procedure(support_side_interface), deferred :: side

  end type graph_support

  !===================================================================!
  ! GRAPH_FUNCTIONAL. The grammar's field at domain size one: a
  ! single value, with the whole inherited battery intact.
  !
  !      a field    [ v1 v2 v3 v4 v5 v6 ]     values over members
  !      a functional      [ J ]              one value: the answer
  !
  ! Nothing is declared here because nothing new can be asked: name,
  ! units, kind, and the get/set roads all mean exactly what they
  ! meant one level down, read at one entry. The type exists so that
  ! an argument may DEMAND the one-entry case at compile time - a
  ! reduction returns a functional, not a field that happens to be
  ! small. The suite holds the law the compiler cannot:
  !
  !      num_entries() == 1
  !===================================================================!

  type, abstract, extends(graph_field) :: graph_functional
  end type graph_functional

  !===================================================================!
  ! GRAPH_REDUCTION. Many values become one: field -> functional.
  !
  ! The whole four-step dance, and why it is four steps:
  !
  !   part 1  [q q q]  --accumulate-->  (sum 6, count 3) --+
  !                                                        +--> combine
  !   part 2  [q q]    --accumulate-->  (sum 14, count 2) -+       |
  !                                                                v
  !                                                   (sum 20, count 5)
  !                                                                |
  !                                                            finalize
  !                                                                |
  !                                                                v
  !                                                           J = 4.0
  !
  ! Means of 2 and 7 do not average to 4.5. They average to 4,
  ! because the sum and the count travel together and the division
  ! happens once, at the very end. A reduction that finished early
  ! on each part would get this wrong.
  !
  ! The measure argument is what turns a bare sum into an integral -
  ! weight each cell by its volume, or each face by its area, and
  ! the answer stops depending on how finely the mesh was cut. The
  ! measure seat is also the inner product's second field:
  !
  !      sum       J = sum q_i
  !      integral  J = sum q_i V_i          <- measure is the volume
  !      average   J = sum q_i V_i / sum V_i
  !      norm      J = sqrt( sum q_i^2 V_i )
  !      product   J = sum u_i v_i          <- measure is the field v
  !===================================================================!

  type, abstract :: graph_reduction

   contains

     ! The four steps, for a caller that owns the parallel dance.
     procedure(reduction_initialize_interface), deferred :: initialize
     procedure(reduction_accumulate_interface), deferred :: accumulate
     procedure(reduction_combine_interface)   , deferred :: combine
     procedure(reduction_finalize_interface)  , deferred :: finalize

     ! All four in one call, for a caller holding the whole thing.
     procedure(reduction_reduce_interface), deferred :: reduce

  end type graph_reduction

  !===================================================================!
  ! GRAPH_BROADCAST. One value becomes many: functional -> field.
  ! The transpose of the reduction, and the round trip is the law:
  !
  !      average( copy( c ) )  = c
  !      sum(     share( c ) ) = c
  !
  ! One step, not four. A reduction walks parts and must combine
  ! their partial answers; a broadcast writes the same fill on every
  ! part, so there is nothing to communicate and nothing to stage.
  !===================================================================!

  type, abstract :: graph_broadcast

   contains

     procedure(broadcast_interface), deferred :: broadcast

  end type graph_broadcast

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

     !===============================================================!
     ! The support's one question.
     !===============================================================!

     pure integer function support_side_interface(this)
       import :: graph_support
       class(graph_support), intent(in) :: this
     end function support_side_interface

     !===============================================================!
     ! The reduction's four steps and its one-call form. The state is
     ! an unfinished functional, which is what keeps every numeric
     ! kind out of this contract.
     !===============================================================!

     ! Start empty. For a sum that is zero, for a product one, for a
     ! minimum the largest number there is.
     pure subroutine reduction_initialize_interface(this, state)
       import :: graph_reduction, graph_functional
       class(graph_reduction), intent(in) :: this
       class(graph_functional), allocatable, intent(inout) :: state
     end subroutine reduction_initialize_interface

     ! Add one part's values into the running state.
     pure subroutine reduction_accumulate_interface(this, field, state, measure)
       import :: graph_reduction, graph_field, graph_functional
       class(graph_reduction), intent(in) :: this
       class(graph_field), intent(in) :: field
       class(graph_functional), intent(inout) :: state
       class(graph_field), intent(in), optional :: measure
     end subroutine reduction_accumulate_interface

     ! Join two part answers into one. This is what lets the same
     ! reduction serve one image or a thousand, and it must not care
     ! which order the parts arrive in.
     pure subroutine reduction_combine_interface(this, left, right, combined)
       import :: graph_reduction, graph_functional
       class(graph_reduction), intent(in) :: this
       class(graph_functional), intent(in) :: left
       class(graph_functional), intent(in) :: right
       class(graph_functional), allocatable, intent(inout) :: combined
     end subroutine reduction_combine_interface

     ! Finish, once, after every part has been joined in. An average
     ! divides here; a norm takes its root here. Doing either any
     ! earlier is how a parallel run quietly gets a different answer
     ! from a serial one.
     pure subroutine reduction_finalize_interface(this, state, functional)
       import :: graph_reduction, graph_functional
       class(graph_reduction), intent(in) :: this
       class(graph_functional), intent(in) :: state
       class(graph_functional), allocatable, intent(inout) :: functional
     end subroutine reduction_finalize_interface

     ! Not pure: a reduction whose field lies across images must
     ! communicate somewhere, and this is the one place in the tower
     ! where that is allowed.
     subroutine reduction_reduce_interface(this, field, functional, measure)
       import :: graph_reduction, graph_field, graph_functional
       class(graph_reduction), intent(in) :: this
       class(graph_field), intent(in) :: field
       class(graph_functional), allocatable, intent(inout) :: functional
       class(graph_field), intent(in), optional :: measure
     end subroutine reduction_reduce_interface

     !===============================================================!
     ! The broadcast's one step. The field arrives constructed on its
     ! domain, so it knows its own entry count and its own side; only
     ! its values change.
     !===============================================================!

     pure subroutine broadcast_interface(this, functional, field)
       import :: graph_broadcast, graph_functional, graph_field
       class(graph_broadcast) , intent(in)    :: this
       class(graph_functional), intent(in)    :: functional
       class(graph_field)     , intent(inout) :: field
     end subroutine broadcast_interface

     !===============================================================!
     ! The four verbs. Each carries its graph map and its data map,
     ! and the data travels by the very same cut, merge, or split as
     ! the structure, so the values cannot drift out of step with it.
     !===============================================================!

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

     ! Several fine values land on one coarse cell, so the concrete
     ! says how they merge - added for a residual, averaged for a
     ! state, volume-weighted when the cells are unequal.
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

     ! One coarse value feeds several fine ones, so the concrete says
     ! how it lands - copied straight down, or interpolated so the
     ! result stays smooth across the new cells.
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

end module graph_calculus
