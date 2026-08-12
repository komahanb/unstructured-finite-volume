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
!                      THE TEN CONCRETIONS
!
!    field -----[domain size = 1]---------> graph_functional
!    operation -[answer domain = 1]-------> graph_reduction
!    operation -[source domain = 1]-------> graph_broadcast
!    operation -[built from an operation
!                by binding a graph]------> discretization_operator
!    operation -[built from an operation
!                by freezing a state]-----> linearization_operator
!    transform -[verb = partition]--------> graph_partitioner
!    transform -[verb = assemble]---------> graph_assembler
!    transform -[verb = coarsen]----------> graph_coarsener
!    transform -[verb = refine]-----------> graph_refiner
!
! THE FOLD OF THE DUAL PAIR. The reduction and the broadcast once
! stood beside the operation role, for fear that the measure
! argument and the staged parallel combine would not survive the
! apply signature. They survive untouched: a reduction is an
! operation whose ANSWER lives on the one-entry domain - and a
! functional IS a field, so apply carries it; the measure rides as
! the second input field; the staged combine remains the abstract's
! own discipline. The broadcast is the mirror, its SOURCE on the
! one-entry domain. Every extends-chain in the tower now ends at
! one of the grammar's four roots.
!
! THE DERIVED OPERATORS. Two families of operations are built FROM
! operations. The discretization binds a continuous statement to a
! graph's arithmetic - the act that turns pde and ode into algebra -
! and its contract exposes the algebraic structure: the dependency
! pattern is a graph, visible by law, because that is what the
! minimizers one level up interrogate. The linearization freezes a
! statement at a state, so a nonlinear question becomes the linear
! one newton governs. Both bind one thing and defer the rest; a
! parent over the pair would own no symbol, and generation refuses
! it a name.
!
! The census of this level: ten types, seventeen operation symbols.
! With the grammar beneath, the tower stands at fourteen types and
! seventy-two symbols.
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

  use iso_fortran_env, only : dp => REAL64
  use graph_grammar  , only : graph, graph_field, graph_operation, &
       &                      graph_transform

  implicit none

  private

  public :: graph_functional
  public :: graph_reduction
  public :: graph_broadcast
  public :: discretization_operator
  public :: linearization_operator
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

  type, abstract, extends(graph_operation) :: graph_reduction

   contains

     ! The four steps, for a caller that owns the parallel dance.
     procedure(reduction_initialize_interface), deferred :: initialize
     procedure(reduction_accumulate_interface), deferred :: accumulate
     procedure(reduction_combine_interface)   , deferred :: combine
     procedure(reduction_finalize_interface)  , deferred :: finalize

     ! All four in one call, for a caller holding the whole thing.
     ! The inherited apply is this verb's operation face: the field
     ! reduced over its own domain, the measure riding as the second
     ! input field, the functional leaving as the output - for a
     ! functional IS a field, one entry wide.
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

  type, abstract, extends(graph_operation) :: graph_broadcast

   contains

     procedure(broadcast_interface), deferred :: broadcast

  end type graph_broadcast

  !===================================================================!
  ! DISCRETIZATION_OPERATOR. An operation built from an operation by
  ! binding it to a graph's arithmetic - the act that turns pde and
  ! ode into algebra. A scheme is a MOTIF stamped along the domain
  ! graph:
  !
  !      backward euler    o<--o          reach 1, weights [1, -1]
  !      bdf-k             o<--o<-..<--o  reach k, the k-step table
  !      two-point flux    one mesh edge, weights -+ 1/delta
  !      fitted, degree p  the p-ring, weights solved by the fit
  !
  ! What every member owes by contract: its dependency PATTERN, as a
  ! graph. The support is to a field what this pattern is to a
  ! derived operator - values sit on members, arithmetic flows on
  ! pairs. The minimizers one level up interrogate the pattern - the
  ! diagonal, the colouring, the triangularity, the Galerkin road -
  ! so it is exposed by law, never by inspection.
  !===================================================================!

  type, abstract, extends(graph_operation) :: discretization_operator

   contains

     procedure(discretization_pattern_interface), deferred :: dependencies

  end type discretization_operator

  !===================================================================!
  ! LINEARIZATION_OPERATOR. An operation built from an operation by
  ! freezing a state: the tangent of S at q, wearing the operation
  ! face, so a governed minimizer sees an ordinary linear question,
  !
  !      J v  at  q        the derivative of S along v, at the
  !                        standing state
  !
  ! One deferred verb beyond the face: freeze, which moves the
  ! standing state (and may cache the base residual) between
  ! newton's steps. The difference road divides two residuals; the
  ! exact road, when a statement speaks its own tangent, joins as a
  ! second concretion with no change here.
  !===================================================================!

  type, abstract, extends(graph_operation) :: linearization_operator

   contains

     procedure(linearization_freeze_interface), deferred :: freeze

  end type linearization_operator

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
          !===============================================================!


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
     ! The derived operators' two verbs: the pattern a discretization
     ! owes, and the freeze a linearization answers.
     !===============================================================!

     subroutine discretization_pattern_interface(this, pattern)
       import :: discretization_operator, graph
       class(discretization_operator), intent(in) :: this
       class(graph), allocatable, intent(out)     :: pattern
     end subroutine discretization_pattern_interface

     subroutine linearization_freeze_interface(this, at, base)
       import :: linearization_operator, dp
       class(linearization_operator), intent(inout) :: this
       real(dp), intent(in)           :: at(:)
       real(dp), intent(in), optional :: base(:)
     end subroutine linearization_freeze_interface

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

! NOTE (phase 5B): GRAPH_SIDE_VERTEX / GRAPH_SIDE_EDGE survive only
! as the legacy differential operator's finite LANDING choice - which
! side its answer lands on. They no longer describe any field's
! domain: domains are member_set identities, and the old
! support-as-edgeless-graph ontology is retired.
