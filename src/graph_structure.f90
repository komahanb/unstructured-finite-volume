!=====================================================================!
! LEVEL 3 OF THE NEW TOWER . THE RELATIONAL GRAPH
!
! The container the whole redesign aims at (AGENTS.md 14):
!
!      G = ( S, R )
!
! a structured collection of member sets and typed relations over
! them - and NOTHING ELSE. No vertex, no edge, no tail, no head:
! those words belong to profiles, one level up, that interpret a
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
! each slot of each relation's signature answers same_as against
! one of the owned carriers, or the graph refuses to exist
! (AGENTS.md 62). Two relations of one signature coexist freely -
! identity is the address, never the signature.
!
!                  OWNERSHIP, DECLARED AND KEPT
!
! The graph OWNS stable relations; views and fibre borrows BORROW
! them. So the accessors answer POINTERS into owned storage, never
! copies: a profile that borrows a relation through relation_at may
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
! the carriers and the relations - one law, one roll, three hands.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module graph_structure

  use graph_identity, only : token, mint_token
  use graph_carrier , only : member_set
  use graph_relation, only : relation

  implicit none

  private
  public :: relational_graph, held_set, held_relation

  !===================================================================!
  ! The seats. A Fortran array carries one dynamic type, so the
  ! graph's heterogeneous collections sit in wrapper seats - any
  ! carrier concretion, any relation concretion, side by side.
  !===================================================================!

  type :: held_set

     class(member_set), allocatable :: carrier

  end type held_set

  interface held_set
     module procedure hold_set
  end interface held_set

  type :: held_relation

     class(relation), allocatable :: contents

  end type held_relation

  interface held_relation
     module procedure hold_relation
  end interface held_relation

  !===================================================================!
  ! The relational graph: sets, relations, an identity, and the
  ! questions a container answers. Immutable after construction.
  !===================================================================!

  type :: relational_graph

     type(token)                  , private :: identity
     character(len=:), allocatable, private :: label

     type(held_set)     , allocatable, private :: sets(:)
     type(held_relation), allocatable, private :: relations(:)

   contains

     procedure :: num_member_sets
     procedure :: member_set_at
     procedure :: num_relations
     procedure :: relation_at
     procedure :: holds_set

     procedure :: declare
     procedure :: id
     procedure :: same_as
     procedure :: name

  end type relational_graph

  interface relational_graph
     module procedure create_graph
  end interface relational_graph

contains

  !===================================================================!
  ! Wrap one carrier, or one relation, for a seat.
  !===================================================================!

  type(held_set) function hold_set(carrier) result(this)

    class(member_set), intent(in) :: carrier

    allocate(this % carrier, source=carrier)

  end function hold_set

  type(held_relation) function hold_relation(contents) result(this)

    class(relation), intent(in) :: contents

    allocate(this % contents, source=contents)

  end function hold_relation

  !===================================================================!
  ! Declare a graph: a name, the member sets, the relations. The
  ! refusals that guard the level, in the order they are checked:
  !
  !     an empty or unsigned seat     a graph holds declared domains
  !                                   only
  !     the same domain twice         S_i /= S_j: a collection of
  !                                   identified objects holds each
  !                                   once
  !     a hollow or unsigned relation a graph holds declared
  !                                   relations only
  !     a borrowing view              a graph owns whole relations; a
  !                                   view rides above, never inside
  !     the same relation twice       R_i /= R_j, by identity
  !     a foreign domain              every slot of every relation
  !                                   names one of the graph's own
  !                                   member sets
  !===================================================================!

  type(relational_graph) function create_graph(name, sets, relations) &
       & result(this)

    character(len=*)   , intent(in) :: name
    type(held_set)     , intent(in) :: sets(:)
    type(held_relation), intent(in) :: relations(:)

    class(member_set), allocatable :: d
    integer                        :: r, k, s
    logical                        :: found

    do s = 1, size(sets)
       if (.not. allocated(sets(s) % carrier)) then
          error stop 'graph_structure: a graph holds declared domains only'
       end if
       if (.not. sets(s) % carrier % same_as(sets(s) % carrier)) then
          error stop 'graph_structure: a graph holds declared domains only'
       end if
       do k = 1, s - 1
          if (sets(s) % carrier % same_as(sets(k) % carrier)) then
             error stop 'graph_structure: a graph holds each domain once'
          end if
       end do
    end do

    do r = 1, size(relations)
       if (.not. allocated(relations(r) % contents)) then
          error stop 'graph_structure: a graph holds declared relations only'
       end if
       if (.not. relations(r) % contents % same_as(relations(r) % contents)) then
          error stop 'graph_structure: a graph holds declared relations only'
       end if
       if (.not. relations(r) % contents % materialized()) then
          error stop 'graph_structure: a graph owns whole relations; a view cannot be owned'
       end if
       do k = 1, r - 1
          if (relations(r) % contents % same_as(relations(k) % contents)) then
             error stop 'graph_structure: a graph holds each relation once'
          end if
       end do
    end do

    do r = 1, size(relations)
       do k = 1, relations(r) % contents % arity()
          d = relations(r) % contents % domain(k)
          found = .false.
          do s = 1, size(sets)
             found = found .or. d % same_as(sets(s) % carrier)
          end do
          if (.not. found) then
             error stop 'graph_structure: a relation must relate the graph''s own member sets'
          end if
       end do
    end do

    allocate(this % sets(size(sets)))
    do s = 1, size(sets)
       allocate(this % sets(s) % carrier, source=sets(s) % carrier)
    end do

    allocate(this % relations(size(relations)))
    do r = 1, size(relations)
       allocate(this % relations(r) % contents, &
            &   source=relations(r) % contents)
    end do

    call this % declare(name)

  end function create_graph

  pure integer function num_member_sets(this)

    class(relational_graph), intent(in) :: this

    num_member_sets = size(this % sets)

  end function num_member_sets

  pure integer function num_relations(this)

    class(relational_graph), intent(in) :: this

    num_relations = size(this % relations)

  end function num_relations

  !===================================================================!
  ! The accessors: references into owned, stable storage - the
  ! ownership policy made flesh. A borrower's fibres stay good for
  ! the graph's whole life.
  !===================================================================!

  function member_set_at(this, k) result(carrier)

    class(relational_graph), target, intent(in) :: this
    integer                        , intent(in) :: k
    class(member_set), pointer                  :: carrier

    carrier => this % sets(k) % carrier

  end function member_set_at

  function relation_at(this, k) result(contents)

    class(relational_graph), target, intent(in) :: this
    integer                        , intent(in) :: k
    class(relation), pointer                    :: contents

    contents => this % relations(k) % contents

  end function relation_at

  !===================================================================!
  ! Does this declared domain sit in the graph's collection.
  !===================================================================!

  pure logical function holds_set(this, carrier)

    class(relational_graph), intent(in) :: this
    class(member_set)      , intent(in) :: carrier

    integer :: s

    holds_set = .false.
    do s = 1, size(this % sets)
       holds_set = holds_set .or. carrier % same_as(this % sets(s) % carrier)
    end do

  end function holds_set

  !===================================================================!
  ! The identity block, one law with the carriers' and the
  ! relations'.
  !===================================================================!

  subroutine declare(this, name)

    class(relational_graph), intent(inout)        :: this
    character(len=*)       , intent(in), optional :: name

    if (this % identity % declared()) then
       error stop 'graph_structure: a graph never signs twice'
    end if

    this % identity = mint_token()
    if (present(name)) this % label = name

  end subroutine declare

  pure type(token) function id(this)

    class(relational_graph), intent(in) :: this

    id = this % identity

  end function id

  pure logical function same_as(this, other)

    class(relational_graph), intent(in) :: this
    class(relational_graph), intent(in) :: other

    same_as = this % identity % matches(other % identity)

  end function same_as

  function name(this)

    class(relational_graph), intent(in) :: this
    character(len=:), allocatable       :: name

    if (allocated(this % label)) then
       name = this % label
    else
       name = ''
    end if

  end function name

end module graph_structure
