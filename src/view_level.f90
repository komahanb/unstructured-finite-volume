!=====================================================================!
! LEVEL VIEW
!
! One view of a graph, in which the two branches are interpreted as
! one level of a nested hierarchy:
!
!     branch(1) = the members of this level, as a sequence
!     branch(2) = the coupling among them, as a relational pair
!
! A member is itself a graph read under this same view, so the
! hierarchy is the branch recursion and no level introduces a new
! kind of object. The three views divide as follows: this one
! names the two branches, view_sequence traverses branch(1), and
! view_relational reads branch(2) as its carriers and its relations.
!
!                        WHERE A LEVEL ENDS
!
! A leaf spends no branch(1): its members are values rather than
! graphs, and their extent is held as a counted set representation.
! That boundary keeps a domain of N freedoms at O(1) semantic
! objects instead of N.
!
! A leaf may still carry branch(2). A component of an ordinary
! differential equation has no coupling and leaves it NULL; the same
! component of a field problem carries the spatial coupling there.
!
!                       THE ONE CONSISTENCY CHECK
!
! A coupling's carriers are graphs: this level's members first, in
! the order the spine lists them, and after them the constraint
! instances the relations run into. So the member spine must be a
! prefix of the carrier spine, compared by identity. Equal counts
! are not that claim, and two levels whose members were built
! separately are indistinguishable by count. level_consistent is
! what refuses a coupling that belongs to another level.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module view_level

  use graph_fractal , only : graph, branch, BRANCH_NULL, BRANCH_KNOWN, &
       & known_branch, null_branch
  use view_sequence , only : sequence_num_elements, sequence_element, &
       & sequence_empty, sequence_first, sequence_rest

  implicit none

  private
  public :: level_is_leaf, level_num_members, level_member, level_members
  public :: level_couples, level_coupling, level_consistent
  public :: level_storage

  !===================================================================!
  ! THE OWNER OF A HIERARCHY
  !
  ! Branch references do not own their targets, so every graph a
  ! hierarchy is built from - the level nodes and the spine cells
  ! alike - has to outlive the branches pointing at it. One storage
  ! owns them all, and a level is addressed by its index in it.
  !
  ! Each node is allocated on its own and held by pointer, never as
  ! an element of the array: the array grows by move_alloc and its
  ! elements move, while separately allocated targets do not. That is
  ! the same arrangement relational_binding uses, and for the same
  ! measured reason.
  !
  ! Assignment is refused at run time. A copy would carry the
  ! pointers of the original, and freeing either would leave the
  ! other referencing released storage. No mechanism in the language
  ! forbids the copy at compile time, so it is stopped when attempted.
  !===================================================================!

  type :: node_holder
     type(graph), pointer :: node => null()
  end type node_holder

  type :: level_storage

     type(node_holder), allocatable, private :: nodes(:)
     integer          , private              :: filled = 0

   contains

     procedure :: fresh
     procedure :: node
     procedure :: num_nodes
     procedure :: spine
     procedure :: assemble
     procedure :: couple

     procedure, private :: refuse_assignment
     generic :: assignment(=) => refuse_assignment

     final :: release

  end type level_storage

contains

  !===================================================================!
  ! A level is a leaf when it spends no branch(1). Its members are
  ! then values held in a map rather than graphs on a spine.
  !===================================================================!

  logical function level_is_leaf(g) result(leaf)

    type(graph), intent(in) :: g

    leaf = g % branch(1) % status() == BRANCH_NULL

  end function level_is_leaf

  !===================================================================!
  ! How many members this level holds. A leaf holds none. A branch(1)
  ! that reaches UNKNOWN stops the program inside view_sequence,
  ! because the count is not determined and reporting zero would be
  ! indistinguishable from a leaf.
  !===================================================================!

  integer function level_num_members(g) result(n)

    type(graph), intent(in) :: g

    if (level_is_leaf(g)) then
       n = 0
       return
    end if

    n = sequence_num_elements(g % branch(1))

  end function level_num_members

  !===================================================================!
  ! Member k of this level, in the order the spine lists them. A leaf
  ! stops the program: it has no members to index. An index outside
  ! the spine is refused by view_sequence.
  !===================================================================!

  function level_member(g, k) result(member)

    type(graph), intent(in) :: g
    integer    , intent(in) :: k
    type(graph), pointer    :: member

    if (level_is_leaf(g)) then
       error stop 'view_level: a leaf has no members to index'
    end if

    member => sequence_element(g % branch(1), k)

  end function level_member

  !===================================================================!
  ! The members as a sequence, for a traversal that reaches them once
  ! through sequence_first and sequence_rest. A leaf's sequence is
  ! empty.
  !===================================================================!

  function level_members(g) result(members)

    type(graph), intent(in) :: g
    type(branch)            :: members

    members = g % branch(1)

  end function level_members

  !===================================================================!
  ! Whether a coupling is present. A level whose members do not read
  ! one another leaves branch(2) NULL, which is a different statement
  ! from a coupling that has not yet been built.
  !===================================================================!

  logical function level_couples(g) result(couples)

    type(graph), intent(in) :: g

    couples = g % branch(2) % status() == BRANCH_KNOWN

  end function level_couples

  !===================================================================!
  ! The coupling of this level. Asking for one that is absent stops
  ! the program, because a caller that gathers along edges cannot
  ! proceed on a disassociated reference.
  !===================================================================!

  function level_coupling(g) result(coupling)

    type(graph), intent(in) :: g
    type(graph), pointer    :: coupling

    if (.not. level_couples(g)) then
       error stop 'view_level: this level carries no coupling'
    end if

    coupling => g % branch(2) % known()

  end function level_coupling

  !===================================================================!
  ! Whether a sequence begins with another, by identity: an empty
  ! prefix is a prefix of anything; a nonempty prefix needs a first
  ! element that is the same graph and a rest that is again a prefix.
  ! A carrier spine that runs out first is a disagreement, not an
  ! error.
  !===================================================================!

  recursive logical function begins_with(carriers, members) result(agrees)

    type(branch), intent(in) :: carriers, members

    type(graph), pointer :: x, y

    if (sequence_empty(members)) then
       agrees = .true.
       return
    end if

    if (sequence_empty(carriers)) then
       agrees = .false.
       return
    end if

    x => sequence_first(carriers)
    y => sequence_first(members)

    agrees = x % same_as(y)
    if (agrees) agrees = begins_with(sequence_rest(carriers), sequence_rest(members))

  end function begins_with

  !===================================================================!
  ! Whether the coupling's carriers begin with this level's members,
  ! compared by identity. A level with no coupling is consistent,
  ! there being nothing to disagree with.
  !===================================================================!

  logical function level_consistent(g) result(agrees)

    type(graph), intent(in) :: g

    type(graph), pointer :: coupling

    agrees = .true.
    if (.not. level_couples(g)) return

    coupling => level_coupling(g)
    agrees   =  begins_with(coupling % branch(1), g % branch(1))

  end function level_consistent

  !===================================================================!
  ! A fresh graph, allocated on its own, identity assigned, owned by
  ! this storage. The index it is named by is how every other
  ! procedure here names it.
  !===================================================================!

  integer function fresh(this) result(at)

    class(level_storage), intent(inout) :: this

    type(node_holder), allocatable :: grown(:)

    if (.not. allocated(this % nodes)) allocate(this % nodes(8))

    if (this % filled == size(this % nodes)) then
       allocate(grown(2 * this % filled))
       grown(1:this % filled) = this % nodes
       call move_alloc(grown, this % nodes)
    end if

    this % filled = this % filled + 1
    at = this % filled

    allocate(this % nodes(at) % node)
    call this % nodes(at) % node % declare()

  end function fresh

  !===================================================================!
  ! The graph an index names. An index outside 1 .. num_nodes stops
  ! the program: a branch built on a disassociated reference would
  ! break the kernel's iff between status and association.
  !===================================================================!

  function node(this, at) result(g)

    class(level_storage), intent(in) :: this
    integer             , intent(in) :: at
    type(graph), pointer :: g

    if (at < 1 .or. at > this % filled) then
       error stop 'view_level: the index names a node this storage owns'
    end if

    g => this % nodes(at) % node

  end function node

  pure integer function num_nodes(this)

    class(level_storage), intent(in) :: this

    num_nodes = this % filled

  end function num_nodes

  !===================================================================!
  ! A spine over the given members: no members is the empty spine,
  ! otherwise a cell holding the first member followed by the spine
  ! over the rest. The rest is built before the cell that points at
  ! it, so every target of a KNOWN branch exists when it is named.
  ! The result is the index of the head cell, or zero for the empty spine,
  ! which is what a leaf's branch(1) is built from.
  !===================================================================!

  recursive integer function spine(this, members) result(head)

    class(level_storage), intent(inout) :: this
    integer             , intent(in)    :: members(:)

    type(graph), pointer :: element, rest
    integer :: tail

    if (size(members) == 0) then
       head = 0
       return
    end if

    tail    =  this % spine(members(2:))
    head    =  this % fresh()
    element => this % nodes(members(1)) % node

    this % nodes(head) % node % branch(1) = known_branch(element)

    if (tail == 0) then
       this % nodes(head) % node % branch(2) = null_branch()
    else
       rest => this % nodes(tail) % node
       this % nodes(head) % node % branch(2) = known_branch(rest)
    end if

  end function spine

  !===================================================================!
  ! One level: its members on a spine in branch(1), its coupling in
  ! branch(2). A coupling index of zero leaves branch(2) NULL, which
  ! is the level whose members do not read one another. A coupling
  ! whose carriers do not begin with this level's members, by
  ! identity, stops the program - the two are indistinguishable by
  ! count.
  !===================================================================!

  integer function assemble(this, members, coupling) result(at)

    class(level_storage), intent(inout) :: this
    integer             , intent(in)    :: members(:)
    integer             , intent(in)    :: coupling

    type(graph), pointer :: first, pairing
    integer :: head

    head = this % spine(members)
    at   = this % fresh()

    if (head == 0) then
       this % nodes(at) % node % branch(1) = null_branch()
    else
       first => this % nodes(head) % node
       this % nodes(at) % node % branch(1) = known_branch(first)
    end if

    if (coupling == 0) then
       this % nodes(at) % node % branch(2) = null_branch()
    else
       pairing => this % nodes(coupling) % node
       this % nodes(at) % node % branch(2) = known_branch(pairing)
    end if

    if (.not. level_consistent(this % nodes(at) % node)) then
       error stop 'view_level: the coupling carries this level''s own members'
    end if

  end function assemble

  !===================================================================!
  ! A relational node: its carriers on a spine in branch(1), its
  ! relations on a spine in branch(2). Both branches are spines,
  ! which is what separates this node from a level, whose branch(2)
  ! holds one coupling graph and whose carriers are checked against
  ! its members. Nothing is checked here: the agreement a coupling
  ! owes its level is checked by that level, and the agreement it
  ! owes its relations is checked by relational_valid.
  !===================================================================!

  integer function couple(this, carriers, relations) result(at)

    class(level_storage), intent(inout) :: this
    integer             , intent(in)    :: carriers(:), relations(:)

    type(graph), pointer :: first_carrier, first_relation
    integer :: carrier_head, relation_head

    carrier_head  = this % spine(carriers)
    relation_head = this % spine(relations)
    at            = this % fresh()

    if (carrier_head == 0) then
       this % nodes(at) % node % branch(1) = null_branch()
    else
       first_carrier => this % nodes(carrier_head) % node
       this % nodes(at) % node % branch(1) = known_branch(first_carrier)
    end if

    if (relation_head == 0) then
       this % nodes(at) % node % branch(2) = null_branch()
    else
       first_relation => this % nodes(relation_head) % node
       this % nodes(at) % node % branch(2) = known_branch(first_relation)
    end if

  end function couple

  !===================================================================!
  ! A storage lends pointers into its own nodes, so a copy would
  ! share them and either release would strand the other.
  !===================================================================!

  subroutine refuse_assignment(lhs, rhs)

    class(level_storage), intent(out) :: lhs
    class(level_storage), intent(in)  :: rhs

    associate (u1 => lhs, u2 => rhs); end associate

    error stop 'view_level: a level storage is not assignable'

  end subroutine refuse_assignment

  !===================================================================!
  ! Release every node this storage allocated.
  !===================================================================!

  subroutine release(this)

    type(level_storage), intent(inout) :: this

    integer :: k

    if (.not. allocated(this % nodes)) return

    do k = 1, this % filled
       if (associated(this % nodes(k) % node)) deallocate(this % nodes(k) % node)
    end do

    deallocate(this % nodes)
    this % filled = 0

  end subroutine release

end module view_level
