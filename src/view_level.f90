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
! A leaf uses no branch(1): its members are values rather than
! graphs, and their extent is stored as a counted set representation.
! That boundary keeps a domain of N freedoms at O(1) semantic
! objects instead of N.
!
! A leaf may still have branch(2). A component of an ordinary
! differential equation has no coupling and leaves it NULL; the same
! component of a field problem stores the spatial coupling there.
!
!                       THE ONE CONSISTENCY CHECK
!
! A coupling's carriers are graphs: this level's members first, in
! the order of the member list, and after them the constraint
! instances the relations map into. So the member list must be a
! prefix of the carrier list, compared by identity. Equal counts
! are not that claim, and two levels whose members were built
! separately are indistinguishable by count. level_consistent is
! what rejects a coupling that belongs to another level.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module view_level

  use graph_fractal , only : graph, branch, BRANCH_NULL, BRANCH_KNOWN, &
       & known_branch, null_branch
  use view_sequence , only : sequence_num_elements, sequence_element, &
       & sequence_empty, sequence_first, sequence_rest
  use util_counted_storage, only : counted_storage, counted_reference

  implicit none

  private
  public :: level_is_leaf, level_num_members, level_member, level_members
  public :: level_couples, level_coupling, level_consistent
  public :: level_storage

  !===================================================================!
  ! THE OWNERS OF A HIERARCHY
  !
  ! Branch references do not own their targets, so every graph a
  ! hierarchy is built from - the level nodes and the list cells
  ! alike - has to outlive the branches pointing at it. One cell
  ! owns them all, and a level is addressed by its index in it.
  !
  ! Each node is allocated separately and referenced by pointer, never as
  ! an element of the array: the array grows by move_alloc and its
  ! elements move, while separately allocated targets do not. That is
  ! the same arrangement relational_binding uses, and for the same
  ! measured reason.
  !
  ! A storage value is a counted reference to the cell
  ! (util_counted_storage). Assignment binds one more owner of the
  ! same nodes; the last owner's finalization deallocates them. A
  ! cell with more than one owner is immutable: extension is refused,
  ! so every owner reads the hierarchy it was bound to. A copy made
  ! without defined assignment (allocate with source=, a structure
  ! constructor, polymorphic assignment) is not an owner: the first
  ! of the two to be finalized releases the binding and the other is
  ! refused at its next access, never read from released storage.
  !===================================================================!

  type :: node_pointer
     type(graph), pointer :: node => null()
  end type node_pointer

  type, extends(counted_storage) :: level_nodes

     type(node_pointer), allocatable :: nodes(:)
     integer                         :: filled = 0

   contains

     procedure :: clear => release_nodes

  end type level_nodes

  type :: level_storage

     type(counted_reference), private :: reference

   contains

     procedure :: allocate_node
     procedure :: node
     procedure :: num_nodes
     procedure :: num_owners
     procedure :: member_list
     procedure :: assemble
     procedure :: couple

     procedure, private :: branch_to
     procedure, private :: cell

  end type level_storage

contains

  !===================================================================!
  ! A level is a leaf when it uses no branch(1). Its members are
  ! then values stored in a map rather than graphs in a list.
  !===================================================================!

  logical function level_is_leaf(g) result(leaf)

    type(graph), intent(in) :: g

    leaf = g % branch(1) % status() == BRANCH_NULL

  end function level_is_leaf

  !===================================================================!
  ! How many members this level contains. A leaf contains none. A branch(1)
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
  ! Member k of this level, in the order of the member list. A leaf
  ! stops the program: it has no members to index. An index outside
  ! the list is rejected by view_sequence.
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
  ! The coupling of this level. Reading one that is absent stops
  ! the program, because a caller that gathers along edges cannot
  ! proceed on a disassociated reference.
  !===================================================================!

  function level_coupling(g) result(coupling)

    type(graph), intent(in) :: g
    type(graph), pointer    :: coupling

    if (.not. level_couples(g)) then
       error stop 'view_level: this level has no coupling'
    end if

    coupling => g % branch(2) % known()

  end function level_coupling

  !===================================================================!
  ! Whether a sequence begins with another, by identity: an empty
  ! prefix is a prefix of anything; a nonempty prefix needs a first
  ! element that is the same graph and a rest that is again a prefix.
  ! A carrier list that ends first is a disagreement, not an
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
  ! The cell this storage is bound to. A storage that was bound and
  ! whose binding is no longer live stops the program: its nodes were
  ! deallocated by the last owner, or this value is a bitwise copy
  ! whose twin released the binding. A storage never bound has no
  ! cell.
  !===================================================================!

  function cell(this) result(nodes)

    class(level_storage), intent(in) :: this
    type(level_nodes), pointer :: nodes

    class(counted_storage), pointer :: storage

    nodes => null()
    if (this % reference % released()) then
       error stop 'view_level: this storage''s hierarchy has been released'
    end if
    storage => this % reference % storage()
    if (.not. associated(storage)) return
    select type (storage)
    type is (level_nodes)
       nodes => storage
    end select

  end function cell

  !===================================================================!
  ! A new graph, allocated separately, identity assigned, owned by
  ! this storage. Every other procedure here refers to the graph by
  ! the index returned. The first node binds this storage to a cell
  ! of its own. A hierarchy with more than one owner is immutable:
  ! extending it stops the program.
  !===================================================================!

  integer function allocate_node(this) result(at)

    class(level_storage), intent(inout) :: this

    type(level_nodes) :: template
    type(level_nodes), pointer :: nodes
    type(node_pointer), allocatable :: expanded_nodes(:)

    if (this % reference % num_owners() > 1) then
       error stop 'view_level: a hierarchy is extended by its sole owner'
    end if
    if (.not. this % reference % live()) call this % reference % acquire(template)
    nodes => this % cell()

    if (.not. allocated(nodes % nodes)) allocate(nodes % nodes(8))

    if (nodes % filled == size(nodes % nodes)) then
       allocate(expanded_nodes(2 * nodes % filled))
       expanded_nodes(1:nodes % filled) = nodes % nodes
       call move_alloc(expanded_nodes, nodes % nodes)
    end if

    nodes % filled = nodes % filled + 1
    at = nodes % filled

    allocate(nodes % nodes(at) % node)
    call nodes % nodes(at) % node % declare()

  end function allocate_node

  !===================================================================!
  ! The graph an index names. An index outside 1 .. num_nodes stops
  ! the program: a branch built on a disassociated reference would
  ! break the kernel's iff between status and association.
  !===================================================================!

  function node(this, at) result(g)

    class(level_storage), intent(in) :: this
    integer             , intent(in) :: at
    type(graph), pointer :: g

    type(level_nodes), pointer :: nodes

    nodes => this % cell()
    if (at < 1 .or. .not. associated(nodes)) then
       error stop 'view_level: the index names a node this storage owns'
    end if
    if (at > nodes % filled) then
       error stop 'view_level: the index names a node this storage owns'
    end if

    g => nodes % nodes(at) % node

  end function node

  ! Zero for a storage without a hierarchy.
  integer function num_nodes(this)

    class(level_storage), intent(in) :: this

    type(level_nodes), pointer :: nodes

    num_nodes = 0
    nodes => this % cell()
    if (associated(nodes)) num_nodes = nodes % filled

  end function num_nodes

  ! The number of storage values bound to this hierarchy; zero for a
  ! storage without one.
  pure integer function num_owners(this)

    class(level_storage), intent(in) :: this

    num_owners = this % reference % num_owners()

  end function num_owners

  !===================================================================!
  ! A list over the given members: no members is the empty list,
  ! otherwise a cell storing the first member followed by the list
  ! over the rest. The rest is built before the cell that points at
  ! it, so every target of a KNOWN branch exists when it is named.
  ! The result is the index of the head cell, or zero for the empty list,
  ! which is what a leaf's branch(1) is built from.
  !===================================================================!

  recursive integer function member_list(this, members) result(head)

    class(level_storage), intent(inout) :: this
    integer             , intent(in)    :: members(:)

    type(graph), pointer :: g
    integer :: tail

    if (size(members) == 0) then
       head = 0
       return
    end if

    tail = this % member_list(members(2:))
    head = this % allocate_node()

    g => this % node(head)
    g % branch(1) = this % branch_to(members(1))
    g % branch(2) = this % branch_to(tail)

  end function member_list

  !===================================================================!
  ! The branch that references the node at an index: NULL for index
  ! zero, KNOWN -> node otherwise.
  !===================================================================!

  function branch_to(this, at) result(b)

    class(level_storage), intent(in) :: this
    integer             , intent(in) :: at
    type(branch)                     :: b

    type(graph), pointer :: g

    if (at == 0) then
       b = null_branch()
    else
       g => this % node(at)
       b = known_branch(g)
    end if

  end function branch_to

  !===================================================================!
  ! One level: its members as a list in branch(1), its coupling in
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

    type(graph), pointer :: g
    integer :: head

    head = this % member_list(members)
    at   = this % allocate_node()

    g => this % node(at)
    g % branch(1) = this % branch_to(head)
    g % branch(2) = this % branch_to(coupling)

    if (.not. level_consistent(g)) then
       error stop 'view_level: the coupling''s carriers begin with this level''s own members'
    end if

  end function assemble

  !===================================================================!
  ! A relational node: its carriers as a list in branch(1), its
  ! relations as a list in branch(2). Both branches are lists,
  ! which is what separates this node from a level, whose branch(2)
  ! stores one coupling graph and whose carriers are checked against
  ! its members. Nothing is checked here: the agreement between a
  ! coupling and its level is checked by that level, and the agreement
  ! between the coupling and its relations is checked by relational_valid.
  !===================================================================!

  integer function couple(this, carriers, relations) result(at)

    class(level_storage), intent(inout) :: this
    integer             , intent(in)    :: carriers(:), relations(:)

    type(graph), pointer :: g
    integer :: carrier_head, relation_head

    carrier_head  = this % member_list(carriers)
    relation_head = this % member_list(relations)
    at            = this % allocate_node()

    g => this % node(at)
    g % branch(1) = this % branch_to(carrier_head)
    g % branch(2) = this % branch_to(relation_head)

  end function couple

  !===================================================================!
  ! Release every node of a cell, when its last owner is finalized.
  !===================================================================!

  subroutine release_nodes(this)

    class(level_nodes), intent(inout) :: this

    integer :: k

    if (.not. allocated(this % nodes)) return

    do k = 1, this % filled
       if (associated(this % nodes(k) % node)) deallocate(this % nodes(k) % node)
    end do

    deallocate(this % nodes)
    this % filled = 0

  end subroutine release_nodes

end module view_level
