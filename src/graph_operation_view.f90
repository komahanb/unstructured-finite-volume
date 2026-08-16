!=====================================================================!
! THE VERBS: OPERATION AND TRANSFORM
!
! Two abstract contracts, and the reason they share one module is the
! import graph rather than tidiness:
!
!     graph_operation    the verb WITHIN a graph: data in, data out
!     graph_transform    the verb BETWEEN graphs: one becomes another
!
! A transform IS a graph-to-graph operation contract. They are named
! in the same files, written against the same three imported names,
! and after the ordinary view moved out there was nothing left to
! disentangle them from. Two modules here would be ceremony, not
! structure - so the measurement was accepted and one module written
! (doc/final-codebase-cutover-plan.md, PR2).
!
! This is the last of graph_grammar, which is deleted with this file's
! arrival. What that module held is now in three places:
!
!     graph                 -> graph_ordinary_view
!     graph_field, kinds    -> graph_field_calculus
!     set_graph             -> fractal_graph
!     `-- the two verbs     -> here
!
! WHAT IT LENDS. Two names, both defined here. graph, set_graph and
! graph_field are IMPORTED to spell the signatures below and are not
! re-exported: a verb within a graph must name the structure it reads,
! the domain its answer lands on, and the data it moves, and it
! defines none of the three.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!
!
!                           THE FOUR ROLES
!
! Every citizen of every level is one of these, and each answers a
! different question. The list outlived the module that first stated
! it; what changed is that each role now lives where it belongs:
!
!    graph ............ structure     graph_ordinary_view
!    graph_field ...... value         graph_field_calculus
!    graph_operation .. verb within   here
!    graph_transform .. verb between  here
!
! A graph here is the ordinary reading: two finite domains joined by
! tail and head. A domain is a set graph, never an edgeless graph; a
! field's domain is an identity and a count, full stop.
!
!=====================================================================!
!
!                         THE ADMISSION LAW
!
! Nothing enters this contract, at any level, except through four
! rules:
!
!    ABSORPTION     a choice from a finite list (vertex or edge, sum
!                   or max, real or complex) rides as an argument or
!                   a constant, never as a type. A type is earned
!                   only when the ROLE of an argument changes.
!
!    GENERATION     a question whose answer composes from other
!                   answers earns no procedure. Only generators
!                   enter; compositions are written at call sites.
!
!    COMMUTATION    two choices are independent axes only if applying
!                   them in either order gives one result. If the
!                   order matters they are one coupled concern and
!                   are handled in one place.
!
!    INHABITATION   no abstract type without a concrete citizen that
!                   answers every symbol meaningfully. Empty types
!                   and dead procedures are structural defects.
!
! Under these rules the minimum is exactly the four roles below. A
! fifth role always either absorbs into one of the four or forces
! dead procedures somewhere. The census of this level: four roles,
! fifty-five operation symbols.
!
!
!=====================================================================!
!
!                  WHAT apply DOES TO A LENT BUFFER
!
! WHAT apply DOES TO A LENT BUFFER. It writes the result into it and
! never adds to what was there. The output argument is intent(inout)
! for one reason only: a caller that already holds a buffer of the
! right shape can lend it and save an allocation. Lending changes the
! cost of the call, not its meaning.
!
!
!=====================================================================!

module graph_operation_view

  !===================================================================!
  ! THREE NAMES ARRIVE, AND ALL THREE ARE SOMEONE ELSE'S.
  !
  !     graph        graph_ordinary_view   the structure it reads
  !     set_graph    fractal_graph         the domain it answers on,
  !     |                                  renamed because the kernel
  !     |                                  calls it `graph` too
  !     graph_field  graph_field_calculus  the values it moves
  !
  ! The rename survives here as a convenience, not a disambiguation:
  ! graph_ordinary_view's abstract type is `ordinary_graph` now, so
  ! nothing this module imports competes for the word `graph`. The
  ! kernel's type is still renamed at the door because a signature
  ! that says set_graph says WHICH question it is asking.
  !===================================================================!

  use graph_ordinary_view , only : ordinary_graph
  use fractal_graph       , only : set_graph => graph
  use graph_field_calculus, only : graph_field

  implicit none

  private

  public :: graph_operation
  public :: graph_transform

  !===================================================================!
  ! GRAPH_OPERATION. The verb within a graph: data in, data out.
  !
  !      input_graph, input_data(:)  --- apply --->  output
  !
  ! Three symbols. name says what it is. domain says where the answer
  ! lives: a member set of the input graph. apply does the work,
  ! under the lent-buffer rule of the header.
  !
  ! A concrete operation is handed the fields it reads when it is
  ! constructed - a coefficient, a measure, a geometry field arrives
  ! as an argument the compiler checks, and apply fetches nothing by
  ! name. This is what keeps the call structure visible from
  ! construction to result.
  !
  ! The concretions of this role, arriving on level 1, split by
  ! order: first-order kernels act on fields (a differential
  ! operator, a balance), and the one higher-order citizen acts on
  ! operators (the walk, a traversal whose kernel is not yet bound).
  ! The maps that touch the one-entry domain - the reduction and the
  ! broadcast - stand beside this role rather than under it, one
  ! deliberate deviation, recorded where they are declared.
  !===================================================================!

  type, abstract :: graph_operation

   contains

     procedure(operation_name_interface)  , deferred :: name
     procedure(operation_domain_interface), deferred :: domain
     procedure(operation_apply_interface) , deferred :: apply

  end type graph_operation

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
     procedure(transform_on_data_interface) , deferred :: defined_on_data

  end type graph_transform

  abstract interface

     !===============================================================!
     ! Verb within: name, domain, apply.
     !===============================================================!

     pure function operation_name_interface(this) result(name)
       import :: graph_operation
       class(graph_operation), intent(in) :: this
       character(len=:), allocatable :: name
     end function operation_name_interface

     !---------------------------------------------------------------!
     ! Where the answer lives: WHICH set, and HOW MANY entries it has.
     !
     ! The count travels beside the identity because every caller of
     ! this symbol wants exactly those two things - to check the
     ! domain matches, and to size a field. Neither is a question
     ! about membership, so neither needs a map, and an operation that
     ! carries a domain carries an integer rather than an extension.
     !---------------------------------------------------------------!

     subroutine operation_domain_interface(this, input_graph, domain, &
          & nentries)
       import :: graph_operation, ordinary_graph, set_graph
       class(graph_operation), intent(in)  :: this
       class(ordinary_graph)          , intent(in)  :: input_graph
       type(set_graph)       , intent(out) :: domain
       integer               , intent(out) :: nentries
     end subroutine operation_domain_interface

     subroutine operation_apply_interface(this, input_graph, input_data, output)
       import :: graph_operation, ordinary_graph, graph_field
       class(graph_operation), intent(in) :: this
       class(ordinary_graph), intent(in) :: input_graph
       class(graph_field), intent(in), optional :: input_data(:)
       class(graph_field), allocatable, intent(inout) :: output
     end subroutine operation_apply_interface

     !===============================================================!
     ! Verb between: the two admissibility questions.
     !===============================================================!

     pure logical function transform_on_graph_interface(this, input_graph)
       import :: graph_transform, ordinary_graph
       class(graph_transform), intent(in) :: this
       class(ordinary_graph), intent(in) :: input_graph
     end function transform_on_graph_interface

     logical function transform_on_data_interface(this, input_graph, input_data)
       import :: graph_transform, ordinary_graph, graph_field
       class(graph_transform), intent(in) :: this
       class(ordinary_graph), intent(in) :: input_graph
       class(graph_field), intent(in) :: input_data
     end function transform_on_data_interface

  end interface

end module graph_operation_view
