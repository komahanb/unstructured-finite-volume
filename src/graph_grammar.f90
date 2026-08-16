!=====================================================================!
! THE LEGACY VERB CONTRACTS: OPERATION AND TRANSFORM
!
! ONCE the ground level of the old stratification, then the one place
! everything legacy met. What is left is two abstract verbs:
!
!     graph_operation    the verb WITHIN a graph: data in, data out
!     graph_transform    the verb BETWEEN graphs
!
! Everything else has gone to the module that owns it. This one is
! being drained, and two draughts are taken
! (doc/final-codebase-cutover-plan.md, PR2):
!
!     graph_field, GRAPH_FIELD_*  ->  graph_field_calculus, which
!     |                               defines them
!     set_graph                   ->  fractal_graph, which mints it
!     `-- graph                   ->  graph_ordinary_view, which is
!                                     what that contract always was
!
! It lends no name it did not define. graph, set_graph and graph_field
! are still IMPORTED, because the two interfaces below are written in
! them - a verb within a graph must name the graph, the domain it
! answers on, and the data it moves. Importing a name to spell a
! signature is not lending it onward, and none of the three is public
! here.
!
! The module is deleted when `use graph_grammar` is empty - not
! before, and not by renaming it to something that sounds newer. Two
! contracts remain, in 25 files, and both have a home named for them.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
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
!=====================================================================!
!
!                           THE FOUR ROLES
!
! Every citizen of every level is one of these, and each answers a
! different question:
!
!    graph ............ structure     what is joined to what
!    graph_field ...... value         what the members carry
!    graph_operation .. verb within   how data becomes other data
!    graph_transform .. verb between  how one graph becomes another
!
! A graph here is the ordinary reading: two finite domains joined by
! tail and head. A domain is a set graph now, never an edgeless
! graph; a field's domain is an identity and a count, full stop.
!
!=====================================================================!
!
!                      WHAT A GRAPH IS MADE OF
!
! A vertex is a thing. An edge joins two of them, tail to head.
!
!                            e
!                     i ----------> j       edge_tail(e) = i
!                                           edge_head(e) = j
!                                           edge_has_head(e) = .true.
!
!                            b
!                     i ----------o         edge_tail(b) = i
!                                           edge_has_head(b) = .false.
!
! The second edge is attached to vertex i alone. That is a boundary
! face, and this is how it is written without inventing an imaginary
! cell on the far side of the wall.
!
!=====================================================================!
!
!                           THE THREE RULES
!
! Three questions come up over and over in the procedures below:
! where does a value sit in a field, can a graph change after it is
! built, and what does apply do to a buffer it is handed? Each has
! one answer, fixed here, and the whole tower is written assuming it.
!
! WHERE A VALUE SITS. Suppose a domain lists three cells in the
! order 7, 3, 11, and a field on that domain holds two components
! per entry. The field keeps its values in that same order, with the
! components of each entry side by side:
!
!        cell            7        7        3        3      ...
!        component       1        2        1        2
!                     +--------+--------+--------+--------+--
!        values       |  v(1)  |  v(2)  |  v(3)  |  v(4)  |
!                     +--------+--------+--------+--------+--
!
! So if a cell is the p-th entry of the domain, its component c is
! the number at position
!
!        (p - 1) * num_components + c
!
! and anyone holding the flat vector finds any value by this formula
! alone. Once the layout is fixed, degree counting and vector
! numbering are one-line formulas, and a formula needs no procedure.
!
! CAN A GRAPH CHANGE? No. Everything a graph holds - structure,
! tags, its relation to the whole it came from - goes in at
! construction, and no procedure below accepts data afterwards.
! When an operation computes something new, the result leaves through
! that operation's output argument. The reason is repeatability: ask
! a graph the same question twice and it gives the same answer twice,
! no matter what ran in between.
!
! WHAT apply DOES TO A LENT BUFFER. It writes the result into it and
! never adds to what was there. The output argument is intent(inout)
! for one reason only: a caller that already holds a buffer of the
! right shape can lend it and save an allocation. Lending changes the
! cost of the call, not its meaning.
!
!=====================================================================!

module graph_grammar


  !===================================================================!
  ! THREE NAMES ARRIVE, AND ALL THREE ARE SOMEONE ELSE'S.
  !
  ! A verb within a graph names three things in its signature and
  ! defines none of them: the graph it acts on, the domain its answer
  ! lands on, and the data it moves. Each comes from its owner.
  !
  !     graph        graph_ordinary_view   the structure it reads
  !     set_graph    fractal_graph         the domain it answers on,
  !     |                                  renamed because the kernel
  !     |                                  calls it `graph` too
  !     graph_field  graph_field_calculus  the values it moves
  !
  ! The collision that forced the set_graph rename is now in
  ! graph_ordinary_view, where the type named `graph` actually lives;
  ! here it survives only because this module still spells signatures
  ! in both.
  !===================================================================!

  use graph_ordinary_view , only : graph
  use fractal_graph       , only : set_graph => graph
  use graph_field_calculus, only : graph_field

  implicit none

  private

  !===================================================================!
  ! WHAT THIS MODULE OWNS, AND IT IS ONLY THIS.
  !
  ! Two abstract types. Nothing imported above is public here, and
  ! that is the rule the drain is following: a module lends only what
  ! it defines. Both survivors have a home named for them -
  ! graph_operation_view and graph_transform_view - and the split
  ! between them is now a two-line decision rather than two phases,
  ! because after the ordinary view left there is nothing else here to
  ! disentangle them from.
  !===================================================================!

  public :: graph_operation
  public :: graph_transform


  !===================================================================!
  ! GRAPH_FIELD. The carrier of values.
  !
  ! A field is a function: values over a domain, and the domain is a
  ! graph - a member set of some host. One value kind per field,
  ! num_components values per entry, laid out by the formula in the
  ! header.
  !
  !      field  ---get--->  [ v1 v2 v3 v4 ]  ---> a solver, a file
  !                                               writer, an outside
  !      field  <--set----  [ v1 v2 v3 v4 ]  <--- library; the answer
  !                                               coming back
  !
  ! Fetch once, work in plain arrays where the arithmetic is free,
  ! write back once. Scaling and adding are not graph theory and have
  ! no procedures here on purpose. One get/set pair per value kind:
  ! the five roads of one absorbed axis.
  !
  ! A field whose domain has ONE entry is a single number. Level 1
  ! names that case the functional; nothing in this role depends on
  ! the size of the domain, which is why the name can wait.
  !
  ! num_entries repeats the size of the domain. That repetition is a
  ! priced convenience for call sites, recorded here so nobody
  ! mistakes it for a generator.
  !===================================================================!

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
       import :: graph_operation, graph, set_graph
       class(graph_operation), intent(in)  :: this
       class(graph)          , intent(in)  :: input_graph
       type(set_graph)       , intent(out) :: domain
       integer               , intent(out) :: nentries
     end subroutine operation_domain_interface

     subroutine operation_apply_interface(this, input_graph, input_data, output)
       import :: graph_operation, graph, graph_field
       class(graph_operation), intent(in) :: this
       class(graph), intent(in) :: input_graph
       class(graph_field), intent(in), optional :: input_data(:)
       class(graph_field), allocatable, intent(inout) :: output
     end subroutine operation_apply_interface

     !===============================================================!
     ! Verb between: the two admissibility questions.
     !===============================================================!

     pure logical function transform_on_graph_interface(this, input_graph)
       import :: graph_transform, graph
       class(graph_transform), intent(in) :: this
       class(graph), intent(in) :: input_graph
     end function transform_on_graph_interface

     logical function transform_on_data_interface(this, input_graph, input_data)
       import :: graph_transform, graph, graph_field
       class(graph_transform), intent(in) :: this
       class(graph), intent(in) :: input_graph
       class(graph_field), intent(in) :: input_data
     end function transform_on_data_interface

  end interface

end module graph_grammar
