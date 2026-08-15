!=====================================================================!
! LEVEL 3 OF THE NEW TOWER . THE RELATED GRAPH
!
! The container the whole redesign aims at (AGENTS.md 14):
!
!      G = ( S, R )
!
! a structured collection of member sets and typed relations over
! them - and NOTHING ELSE. No vertex, no edge, no tail, no head:
! those words belong to interpretations, one level up, that interpret a
! particular schema of relations. This container holds a mesh's
! cells and faces exactly as gladly as a calculator's operations,
! values and ports.
!
!      Graph CONTAINS Relations - composition, never inheritance.
!      A graph is not a relation; a relation needs no graph.
!
!                     THE SIGNATURE VALIDITY LAW
!
! Every relation admitted must relate the graph's own member sets:
! each slot of each relation's signature answers equals against
! one of the owned sets, or the graph refuses to exist
! (AGENTS.md 62). Two relations of one signature coexist freely -
! identity is the address, never the signature.
!
!                  OWNERSHIP, DECLARED AND KEPT
!
! The graph OWNS stable relations; views and fibre borrows BORROW
! them. So the accessors answer POINTERS into owned storage, never
! copies: an interpretation that borrows a relation through relation_at may
! hold fibre views into it for as long as the graph lives, and the
! graph, immutable after construction, never pulls that storage out
! from under anyone.
!
! And ownership means WHOLENESS: only a MATERIALIZED relation - one
! whole unto itself - may be owned. A borrowing view is refused at
! the door, because copying a view into owned storage would copy a
! pointer to a base the graph does not keep alive: the graph would
! own the view and not what makes it true. Views ride ABOVE
! graph-owned relations, never inside them.
!
! The graph is the third citizen to sign the identity roll, after
! the sets and the relations - one law, one roll, three hands.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module graph_structure

  use graph_set , only : graph, set, declared_set, unrelated_graph
  use graph_relation, only : relation

  implicit none

  private
  public :: related_graph, declared_set, declared_relation
  public :: declare_relations, forget_relations

  !===================================================================!
  ! The declared elements. A Fortran array carries one dynamic type, so the
  ! graph's heterogeneous collections sit in declared-element wrappers - any
  ! set concretion, any relation concretion, side by side.
  !===================================================================!

  type :: declared_relation

     class(relation), allocatable :: contents

  end type declared_relation

  interface declared_relation
     module procedure create_declared_relation
  end interface declared_relation

  !===================================================================!
  ! The related graph: sets, relations, an identity, and the
  ! questions a container answers. Immutable after construction.
  !===================================================================!

  type, extends(graph) :: related_graph

     type(declared_set)     , allocatable, private :: sets(:)
     type(declared_relation), allocatable, private :: relations(:)

   contains

     procedure :: num_sets
     procedure :: set_at
     procedure :: num_relations
     procedure :: relation_at
     procedure :: holds_set

  end type related_graph

  interface related_graph
     module procedure create_graph
  end interface related_graph

contains

  !===================================================================!
  ! Declare one element of S, or one element of R.
  !===================================================================!

  type(declared_relation) function create_declared_relation(contents) result(this)

    class(relation), intent(in) :: contents

    allocate(this % contents, source=contents)

  end function create_declared_relation

  !===================================================================!
  ! Declare a graph: a name, the member sets, the relations. The
  ! refusals that guard the level, in the order they are checked:
  !
  !     an empty or unsigned element  a graph holds declared domains
  !                                   only
  !     the same domain twice         S_i /= S_j: a collection of
  !                                   identified objects holds each
  !                                   once
  !     an empty relation family      |R| > 0 IS what related_graph
  !                                   means; a graph over which no
  !                                   relation is declared is an
  !                                   unrelated_graph, and this is
  !                                   not that type
  !     a hollow or unsigned relation a graph holds declared
  !                                   relations only
  !     a borrowing view              a graph owns whole relations; a
  !                                   view rides above, never inside
  !     the same relation twice       R_i /= R_j, by identity
  !     an unheld domain              every slot of every relation
  !                                   names one of the graph's own
  !                                   member sets
  !===================================================================!

  type(related_graph) function create_graph(name, sets, relations) &
       & result(this)

    character(len=*)   , intent(in) :: name
    type(declared_set)     , intent(in) :: sets(:)
    type(declared_relation), intent(in) :: relations(:)

    class(set), allocatable :: d
    integer                        :: r, k, s
    logical                        :: found

    do s = 1, size(sets)
       if (.not. allocated(sets(s) % set)) then
          error stop 'graph_structure: a graph holds declared domains only'
       end if
       if (.not. sets(s) % set % equals(sets(s) % set)) then
          error stop 'graph_structure: a graph holds declared domains only'
       end if
       do k = 1, s - 1
          if (sets(s) % set % equals(sets(k) % set)) then
             error stop 'graph_structure: a graph holds each domain once'
          end if
       end do
    end do

    ! |R| > 0 is the defining law of this type, checked after the
    ! domains because S is prior: a graph's sets are what its
    ! relations are declared OVER. The count is of DECLARED RELATION
    ! OBJECTS, never of tuples - a relation holding no tuples is a
    ! relation, and one of them is enough.
    if (size(relations) == 0) then
       error stop 'graph_structure: a related graph declares at least one relation'
    end if

    do r = 1, size(relations)
       if (.not. allocated(relations(r) % contents)) then
          error stop 'graph_structure: a graph holds declared relations only'
       end if
       if (.not. relations(r) % contents % equals(relations(r) % contents)) then
          error stop 'graph_structure: a graph holds declared relations only'
       end if
       if (.not. relations(r) % contents % materialized()) then
          error stop 'graph_structure: a graph owns whole relations; a view cannot be owned'
       end if
       do k = 1, r - 1
          if (relations(r) % contents % equals(relations(k) % contents)) then
             error stop 'graph_structure: a graph holds each relation once'
          end if
       end do
    end do

    do r = 1, size(relations)
       do k = 1, relations(r) % contents % arity()
          d = relations(r) % contents % domain(k)
          found = .false.
          do s = 1, size(sets)
             found = found .or. d % equals(sets(s) % set)
          end do
          if (.not. found) then
             error stop 'graph_structure: a relation must relate the graph''s own member sets'
          end if
       end do
    end do

    allocate(this % sets(size(sets)))
    do s = 1, size(sets)
       allocate(this % sets(s) % set, source=sets(s) % set)
    end do

    allocate(this % relations(size(relations)))
    do r = 1, size(relations)
       allocate(this % relations(r) % contents, &
            &   source=relations(r) % contents)
    end do

    call this % declare(name)

  end function create_graph

  pure integer function num_sets(this)

    class(related_graph), intent(in) :: this

    num_sets = size(this % sets)

  end function num_sets

  pure integer function num_relations(this)

    class(related_graph), intent(in) :: this

    num_relations = size(this % relations)

  end function num_relations

  !===================================================================!
  ! The accessors: references into owned, stable storage - the
  ! ownership policy made flesh. A borrower's fibres stay good for
  ! the graph's whole life.
  !===================================================================!

  function set_at(this, slot) result(domain)

    class(related_graph), target, intent(in) :: this
    integer             ,         intent(in) :: slot
    class(set)          , pointer            :: domain

    ! set_at : {1 .. num_sets()} -> S, and nowhere else. The range is
    ! the root's contract, so a slot outside it is refused here as it
    ! is on a set, rather than read off the end of the collection.
    if (slot < 1 .or. slot > size(this % sets)) then
       error stop 'graph_structure: set_at is asked outside {1 .. num_sets()}'
    end if

    select type (element => this % sets(slot) % set)
    class is (set)
       domain => element
    class default
       error stop 'graph_structure: a declared element of S is a set, and this one is not'
    end select

  end function set_at

  function relation_at(this, k) result(contents)

    class(related_graph), target, intent(in) :: this
    integer                        , intent(in) :: k
    class(relation), pointer                    :: contents

    contents => this % relations(k) % contents

  end function relation_at

  !===================================================================!
  ! Does this declared domain sit in the graph's collection.
  !===================================================================!

  pure logical function holds_set(this, domain)

    class(related_graph), intent(in) :: this
    class(set)             , intent(in) :: domain

    integer :: s

    holds_set = .false.
    do s = 1, size(this % sets)
       holds_set = holds_set .or. domain % equals(this % sets(s) % set)
    end do

  end function holds_set


  !===================================================================!
  ! THE TRANSITION MAPS between the two branches of the tree.
  !
  !      declare_relations : (S, empty), R  ->  (S, R),  |R| > 0
  !      forget_relations  : (S, R)         ->  (S, empty)
  !
  ! MODULE PROCEDURES, AND NOT BY PREFERENCE. declare_relations
  ! cannot be bound to unrelated_graph: its result is a level-3 type
  ! and unrelated_graph is level 0, so the binding would make the
  ! ground level name the container that is built over it. Its
  ! partner is written the same way rather than bound to
  ! related_graph, because a complementary pair that reads two
  ! different ways at the call site has stopped being a pair.
  !
  ! NEITHER MUTATES. Both build a NEW declared graph and leave the
  ! source untouched, so a borrow already taken from the source stays
  ! good; structural objects are immutable and a transformation is a
  ! new truth, not a rewritten one (AGENTS.md 23, 54).
  !
  ! WHAT SURVIVES, AND WHAT DOES NOT.
  !
  !      S    survives BY IDENTITY - each set of the result equals
  !           the set it came from, because a whole-object copy IS
  !           the declared domain
  !      R    is exactly what the map says: handed in, or emptied
  !      G    does NOT survive: the result signs its own name, and
  !           equals(source, result) is false
  !
  ! NO VALIDATION LIVES HERE. Both maps rebuild the collection from
  ! the source's own root contract and hand it to the constructor
  ! that already owns the laws - so the empty-family refusal, the
  ! view-ownership refusal and the signature-closure refusal are the
  ! SAME refusals, raised in the same place, and cannot drift from
  ! the direct construction.
  !===================================================================!

  type(related_graph) function declare_relations(over, relations, name) &
       & result(this)

    class(unrelated_graph), target, intent(in) :: over
    type(declared_relation)      , intent(in) :: relations(:)
    character(len=*)             , intent(in) :: name

    type(declared_set), allocatable :: domains(:)
    integer                         :: s

    allocate(domains(over % num_sets()))
    do s = 1, over % num_sets()
       domains(s) = declared_set(over % set_at(s))
    end do

    this = related_graph(name, domains, relations)

  end function declare_relations

  type(unrelated_graph) function forget_relations(over, name) result(this)

    class(related_graph), target, intent(in) :: over
    character(len=*)           , intent(in) :: name

    type(declared_set), allocatable :: domains(:)
    integer                         :: s

    allocate(domains(over % num_sets()))
    do s = 1, over % num_sets()
       domains(s) = declared_set(over % set_at(s))
    end do

    this = unrelated_graph(name, domains)

  end function forget_relations

end module graph_structure
