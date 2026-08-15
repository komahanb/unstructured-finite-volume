!=====================================================================!
! LEVEL 0 OF THE NEW TOWER . THE SETS
!
! The ground object of the relation-centered tower (AGENTS.md): a
! MEMBER SET,
!
!      A = { a_1, ..., a_n }
!
! It answers what members may exist, and nothing else. No structure,
! no relations, no graph words - vertices, edges, cells and faces
! are what higher levels call the sets they declare here.
!
! The structural contract is five questions: how many (size), which
! one (member), all of them (members), IS THIS ONE OF YOURS (has),
! and WHERE DOES IT STAND (local_index - the inverse of member,
! zero for an outsider). Membership and position are primitives,
! not searches: a relation signature validates tuples against has,
! and an indexed relation finds its row through local_index, so
! neither ever has to enumerate a domain. The two enumeration laws
! bind them:
!
!      member(local_index(a)) = a        for every a in A
!      local_index(member(i)) = i        for i = 1 .. size
!
! which forces enumeration to be injective: a member set holds each
! member once - set semantics on the sets themselves. The name
! a set was declared with is metadata for the reader - carried,
! printable, and no part of the mathematics.
!
!                       STRUCTURAL IDENTITY
!
! Two member sets are the same only when they are the SAME DECLARED
! DOMAIN - never because they hold equal integers, or count the same,
! or happen to enumerate identical values,
!
!      cells  = { 1 2 3 4 }        four cells
!      faces  = { 1 2 3 4 }        four faces - a different world
!
! The identity machinery - the opaque token, its minting, its one
! law - lives beneath the tower in graph_identity, because relations
! sign the same way and the related graph will be the third to
! do so. Here the law is only applied: a set signs once, at
! declaration; a second signing is refused loudly; id() answers the
! whole token, so identity is never mistaken for an image-local
! integer.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module graph_set

  use graph_identity, only : token, mint_token

  implicit none

  private
  public :: graph, unrelated_graph, declared_set
  public :: set, index_set, subset

  !===================================================================!
  ! THE CANONICAL ROOT. A graph is a collection of sets and a
  ! collection of declared relation objects,
  !
  !      G = (S, R)
  !
  ! and that pair is the whole of what a graph IS. The root asks
  ! three questions about it - how many sets, which one, how many
  ! declared relation objects - and answers three more about itself,
  ! because a graph is a declared object and signs like every other.
  !
  ! WHAT IS NOT HERE, AND WHY. relation_at is absent. Naming a
  ! relation needs the relation type, which is level 1 and is built
  ! over this module; a root that could name one would close the
  ! cycle graph -> set -> relation -> graph and no compiler would
  ! accept it. Counting needs no such type, so num_relations stands
  ! and relation_at waits for related_graph, which is exactly the
  ! structure that level adds. The vertex/edge vocabulary - tail,
  ! head, incident, adjacent, outgoing, incoming - is not here
  ! either, and never will be: those are readings of declared
  ! relations, and they live where readings live (AGENTS.md 7, 16).
  !
  ! THE INDEX RANGE. set_at is declared over {1 .. num_sets()}. The
  ! empty graph G = (empty, empty) is admissible, and for it that
  ! range is empty and set_at is the empty map - so admitting it
  ! costs no axiom about S being occupied. Distinct from it, and not
  ! to be confused: an EMPTY SET is one domain with zero members,
  ! which is a set like any other and answers num_sets() = 1.
  !===================================================================!

  type, abstract :: graph

     type(token)                  , private :: identity
     character(len=:), allocatable, private :: label

   contains

     !----------------------------------------------------------------!
     ! The pair, deferred: how many sets, which one, how many
     ! declared relation objects.
     !----------------------------------------------------------------!

     procedure(graph_num_sets_interface)     , deferred :: num_sets
     procedure(graph_set_at_interface)       , deferred :: set_at
     procedure(graph_num_relations_interface), deferred :: num_relations

     !----------------------------------------------------------------!
     ! Identity, answered once for every graph beneath this root.
     !----------------------------------------------------------------!

     procedure :: declare
     procedure :: id
     procedure :: equals

     !----------------------------------------------------------------!
     ! Metadata, not mathematics: the declared name, or ''.
     !----------------------------------------------------------------!

     procedure :: name

  end type graph

  !===================================================================!
  ! THE DECLARED SET WRAPPER. A Fortran array carries one dynamic
  ! type, so a graph's
  ! heterogeneous collection of sets sits in declared-element wrappers - any set
  ! realization beside any other.
  !
  ! WHY THE COMPONENT IS class(graph) AND NOT class(set). The wrapper
  ! must be declared before unrelated_graph, which must be declared
  ! before set, which extends it: a component typed class(set) would
  ! name a type not yet defined. The declared type is therefore the
  ! root, and the INVARIANT IS KEPT AT THE DOOR - create_declared_set admits
  ! only a set, so every wrapper carries one, and set_at re-establishes
  ! that in the only place it is read.
  !
  ! The recursion - a graph holding graphs - is legal because it
  ! runs through allocatables, and it is honest: a set IS a graph, so
  ! a collection of sets IS a collection of graphs.
  !===================================================================!

  type :: declared_set

     class(graph), allocatable :: set

  end type declared_set

  interface declared_set
     module procedure create_declared_set
  end interface declared_set

  !===================================================================!
  ! THE RELATION-FREE GRAPH. G = (S, empty): a graph over which no
  ! relation object has been declared.
  !
  ! Its law is one fixed answer: num_relations() is 0 and cannot
  ! vary, so a consumer may DEMAND the relation-free state in a
  ! signature instead of testing for it at run time. That is what
  ! earns the type its existence (AGENTS.md 50, clauses 3 and 5).
  !
  ! THE DISTINCTION IS THE OBJECT, NOT THE TUPLE. A declared
  ! relation holding zero tuples is still a declared relation, and a
  ! graph holding one is a related_graph. What is empty here is the
  ! COLLECTION R, never a count of tuples inside it.
  !
  ! CONSTRUCTIBLE FOR ANY S. This type is concrete, so
  !
  !      (empty, empty)        the empty graph
  !      ({A}, empty)          one domain
  !      ({A,B,C}, empty)      three domains, no structure between them
  !
  ! are all declarable objects. `set` is the |S| = 1 specialization
  ! below, and NOTHING SITS BETWEEN THEM.
  !
  ! THE STORAGE COMPROMISE, STATED. Because set extends this
  ! concrete type, every set - and so every index_set and every
  ! subset - carries the collection component. It is never allocated
  ! there, since set overrides num_sets and set_at to answer from
  ! itself, so the cost is one unallocated array descriptor per set
  ! and no heap traffic. That is the price of keeping `set` a
  ! subtype of the object it specializes rather than a sibling of
  ! it, and it is paid deliberately.
  !===================================================================!

  type, extends(graph) :: unrelated_graph

     type(declared_set), allocatable, private :: sets(:)

   contains

     procedure                  :: num_sets      => unrelated_num_sets
     procedure                  :: set_at        => unrelated_set_at
     procedure, non_overridable :: num_relations => unrelated_num_relations
     procedure                  :: holds_set     => unrelated_holds_set

  end type unrelated_graph

  interface unrelated_graph
     module procedure create_unrelated
  end interface unrelated_graph

  !===================================================================!
  ! THE SET. An identity, a count, its members, and membership.
  ! The contract is deliberately this small (AGENTS.md
  ! 4.2); subsets, tags, parts and every other refinement arrive on
  ! higher levels as predicates and relations, never here.
  !
  ! IT IS A GRAPH, AND THIS IS THE CONSTRUCTION. A set is the
  ! one-domain relation-free graph,
  !
  !      A  =  ({A}, empty)
  !
  ! so num_sets() is 1 and set_at(1) is the set itself. Every answer
  ! is total and every answer is the mathematics; nothing degenerate
  ! is fabricated to satisfy the root, which is the whole difference
  ! between this and the edgeless-graph support that AGENTS.md 6 and
  ! 37 refused. That support had to invent edge_tail. This invents
  ! nothing.
  !===================================================================!

  type, abstract, extends(unrelated_graph) :: set

   contains

     !----------------------------------------------------------------!
     ! The graph face: one domain, and it is this one.
     !----------------------------------------------------------------!

     procedure :: num_sets => set_num_sets
     procedure :: set_at   => set_set_at

     !----------------------------------------------------------------!
     ! The structural questions, deferred to each concretion.
     !----------------------------------------------------------------!

     procedure(set_size_interface)       , deferred :: size
     procedure(set_member_interface)     , deferred :: member
     procedure(set_members_interface)    , deferred :: members
     procedure(set_has_interface)        , deferred :: has
     procedure(set_local_index_interface), deferred :: local_index

     !----------------------------------------------------------------!
     ! Embedding, transitive: A is a subobject of itself, and a
     ! subset is a subobject of everything its ambient chain
     ! reaches. This query - never a side flag, never a select type
     ! - is how a consumer asks where a domain ultimately lives.
     !
     ! Identity, naming and the graph pair are inherited from the
     ! root and are not restated here.
     !----------------------------------------------------------------!

     procedure :: is_subobject_of

  end type set

  abstract interface

     !----------------------------------------------------------------!
     ! The root's three questions about the pair G = (S, R).
     !----------------------------------------------------------------!

     pure integer function graph_num_sets_interface(this)
       import graph
       class(graph), intent(in) :: this
     end function graph_num_sets_interface

     function graph_set_at_interface(this, slot) result(domain)
       import graph, set
       class(graph), target, intent(in) :: this
       integer             , intent(in) :: slot
       class(set) , pointer             :: domain
     end function graph_set_at_interface

     pure integer function graph_num_relations_interface(this)
       import graph
       class(graph), intent(in) :: this
     end function graph_num_relations_interface

     pure integer function set_size_interface(this)
       import set
       class(set), intent(in) :: this
     end function set_size_interface

     pure integer function set_member_interface(this, local_index)
       import set
       class(set), intent(in) :: this
       integer          , intent(in) :: local_index
     end function set_member_interface

     pure subroutine set_members_interface(this, indices)
       import set
       class(set)   , intent(in)  :: this
       integer, allocatable, intent(out) :: indices(:)
     end subroutine set_members_interface

     pure logical function set_has_interface(this, member)
       import set
       class(set), intent(in) :: this
       integer          , intent(in) :: member
     end function set_has_interface

     pure integer function set_local_index_interface(this, member)
       import set
       class(set), intent(in) :: this
       integer          , intent(in) :: member
     end function set_local_index_interface

  end interface

  !===================================================================!
  ! THE INDEX SET: the first concrete realization of a set,
  !
  !      A = { 1, ..., n }
  !
  ! A REALIZATION, NOT A RUNG. The taxonomy is graph ->
  ! unrelated_graph -> set -> subset, and index_set is no part of
  ! it: it is one WAY OF BEING a set, chosen because every domain
  ! the repository names today - cells, faces, points, parts,
  ! instants - enumerates exactly so, which makes one integer its
  ! whole storage and membership one comparison. A domain that must
  ! list its members arrives as a second realization the day
  ! something needs it, BESIDE this one and neither above nor below.
  !===================================================================!

  type, extends(set) :: index_set

     integer, private :: n = 0

   contains

     procedure :: size        => index_size
     procedure :: member      => index_member
     procedure :: members     => index_members
     procedure :: has         => index_has
     procedure :: local_index => index_local_index

  end type index_set

  interface index_set
     module procedure create_index
  end interface index_set

  !===================================================================!
  ! The subobject: a subset that IS a member set,
  !
  !      S c--> A          s in S  =>  s in A
  !
  ! This is what the old edgeless-graph "support" was reaching for
  ! (AGENTS.md 6, 37): a chosen family of an ambient domain's
  ! members, itself a declared domain - so a field lives on a
  ! set, ambient or subset, and never learns two domain
  ! kinds. The inclusion law is sealed at construction - every
  ! member must already belong to the ambient - and immutability
  ! keeps it sealed for life. The subset signs its own identity: a
  ! subobject is its own declared domain, never a disguise of its
  ! host. The unary-predicate face is has(); the relational face,
  ! I_S <= S x A, is built by inclusion_of in the binary level.
  !
  ! Membership and standing are honest scans here, as in any listed
  ! roster; a subset that matters at million-member scale owes
  ! itself an indexed concretion, and says so.
  !===================================================================!

  type, extends(set) :: subset

     class(set), allocatable, private :: host
     integer          , allocatable, private :: roll(:)

   contains

     procedure :: size        => subset_size
     procedure :: member      => subset_member
     procedure :: members     => subset_members
     procedure :: has         => subset_has
     procedure :: local_index => subset_local_index

     procedure :: ambient

     procedure :: is_subobject_of => subset_is_subobject_of

  end type subset

  interface subset
     module procedure create_subset
  end interface subset

contains

  !===================================================================!
  ! declare stamps a set with a fresh identity and, if given, its
  ! name. A set signs once; a second signing stops the program,
  ! because a silent second stamp would leave the caller believing a
  ! domain it never made.
  !===================================================================!

  subroutine declare(this, name)

    class(graph)    , intent(inout)         :: this
    character(len=*), intent(in) , optional :: name

    if (this % identity % declared()) then
       error stop 'graph_set: a graph never signs twice'
    end if

    this % identity = mint_token()
    if (present(name)) this % label = name

  end subroutine declare

  !===================================================================!
  ! id answers the whole opaque token - the identity itself, honest
  ! across images, never a bare local integer.
  !===================================================================!

  pure type(token) function id(this)

    class(graph), intent(in) :: this

    id = this % identity

  end function id

  !===================================================================!
  ! The relation-free answer, fixed and non-overridable: no relation
  ! object has been declared over this graph. Not a tuple count - a
  ! declared relation holding no tuples is still declared, and puts
  ! its graph on the other branch of the tree.
  !===================================================================!

  pure integer function unrelated_num_relations(this)

    class(unrelated_graph), intent(in) :: this

    associate (u1 => this); end associate

    unrelated_num_relations = 0

  end function unrelated_num_relations

  !===================================================================!
  ! Declare one element of S. The declared component type is the
  ! root; this door is where it is narrowed, and every element that
  ! exists came through here.
  !===================================================================!

  type(declared_set) function create_declared_set(domain) result(this)

    class(set), intent(in) :: domain

    allocate(this % set, source=domain)

  end function create_declared_set

  !===================================================================!
  ! Declare a relation-free graph: a name and its domains. The
  ! refusals, in the order they are checked:
  !
  !     an empty or unsigned element  a graph holds declared domains
  !                                   only
  !     the same domain twice         S_i /= S_j: a collection of
  !                                   identified objects holds each
  !                                   once
  !
  ! |S| = 0 is admitted. The empty graph is a declared object like
  ! any other, and set_at over an empty range is the empty map.
  !===================================================================!

  type(unrelated_graph) function create_unrelated(name, sets) result(this)

    character(len=*), intent(in) :: name
    type(declared_set)  , intent(in) :: sets(:)

    integer :: s, k

    do s = 1, size(sets)
       if (.not. allocated(sets(s) % set)) then
          error stop 'graph_set: a graph holds declared domains only'
       end if
       if (.not. sets(s) % set % equals(sets(s) % set)) then
          error stop 'graph_set: a graph holds declared domains only'
       end if
       do k = 1, s - 1
          if (sets(s) % set % equals(sets(k) % set)) then
             error stop 'graph_set: a graph holds each domain once'
          end if
       end do
    end do

    allocate(this % sets(size(sets)))
    do s = 1, size(sets)
       allocate(this % sets(s) % set, source=sets(s) % set)
    end do

    call this % declare(name)

  end function create_unrelated

  pure integer function unrelated_num_sets(this)

    class(unrelated_graph), intent(in) :: this

    if (allocated(this % sets)) then
       unrelated_num_sets = size(this % sets)
    else
       unrelated_num_sets = 0
    end if

  end function unrelated_num_sets

  !===================================================================!
  ! set_at : {1 .. num_sets()} -> S, borrowed, refused outside the
  ! range. The select type is where the element's declared root type is
  ! narrowed back to what create_declared_set admitted; the default arm cannot
  ! be reached through the constructor and says so if it ever is.
  !===================================================================!

  function unrelated_set_at(this, slot) result(domain)

    class(unrelated_graph), target, intent(in) :: this
    integer               ,         intent(in) :: slot
    class(set)            , pointer            :: domain

    if (slot < 1 .or. slot > this % num_sets()) then
       error stop 'graph_set: set_at is asked outside {1 .. num_sets()}'
    end if

    select type (element => this % sets(slot) % set)
    class is (set)
       domain => element
    class default
       error stop 'graph_set: a declared element of S is a set, and this one is not'
    end select

  end function unrelated_set_at

  pure logical function unrelated_holds_set(this, domain)

    class(unrelated_graph), intent(in) :: this
    class(set)            , intent(in) :: domain

    integer :: s

    unrelated_holds_set = .false.
    if (.not. allocated(this % sets)) return
    do s = 1, size(this % sets)
       unrelated_holds_set = unrelated_holds_set .or. &
            & domain % equals(this % sets(s) % set)
    end do

  end function unrelated_holds_set

  !===================================================================!
  ! The set's graph face: A = ({A}, empty). One domain, and it is
  ! this one - so set_at hands back THE SET ITSELF, by pointer.
  !
  ! IDENTITY-PRESERVING, WHICH IS THE WHOLE POINT. A copy would
  ! carry the same token and so would answer equals, but it would be
  ! a second object with the graph's storage duplicated behind it. A
  ! borrow is the declared domain, not a likeness of it, and one
  ! rule now governs every graph: set_at borrows, never copies, so a
  ! set and a related_graph answer the same contract the same way.
  !
  ! THE BORROWER'S OBLIGATION, AS EVERYWHERE. The pointer is good
  ! for as long as the graph it came from lives - the ownership
  ! policy graph_structure and graph_interpretation already state.
  !
  ! The slot is checked because {1 .. 1} is the whole of set_at's
  ! index range here.
  !===================================================================!

  pure integer function set_num_sets(this)

    class(set), intent(in) :: this

    associate (u1 => this); end associate

    set_num_sets = 1

  end function set_num_sets

  function set_set_at(this, slot) result(domain)

    class(set), target, intent(in) :: this
    integer           , intent(in) :: slot
    class(set), pointer            :: domain

    if (slot /= 1) then
       error stop 'graph_set: a set holds one domain, and it stands at slot 1'
    end if

    domain => this

  end function set_set_at

  !===================================================================!
  ! equals: one declared domain, or not. Equal contents prove
  ! nothing; an undeclared set equals nothing, itself included.
  !===================================================================!

  pure logical function equals(this, other)

    class(graph), intent(in) :: this
    class(graph), intent(in) :: other

    equals = this % identity % matches(other % identity)

  end function equals

  !===================================================================!
  ! The embedding order: A precedes A, and nothing else at ground
  ! level - the subset overrides this to walk its ambient chain.
  !===================================================================!

  pure logical function is_subobject_of(this, ancestor)

    class(set), intent(in) :: this
    class(set), intent(in) :: ancestor

    is_subobject_of = this % equals(ancestor)

  end function is_subobject_of

  !===================================================================!
  ! name answers the declared label, or '' for a set declared
  ! nameless. Metadata only: no law reads it.
  !===================================================================!

  function name(this)

    class(graph)    , intent(in)  :: this
    character(len=:), allocatable :: name

    if (allocated(this % label)) then
       name = this % label
    else
       name = ''
    end if

  end function name

  !===================================================================!
  ! Declare a counted domain: a name for the reader, a count for the
  ! mathematics, a fresh stamp for the identity.
  !===================================================================!

  type(index_set) function create_index(name, n) result(this)

    character(len=*), intent(in) :: name
    integer         , intent(in) :: n

    this % n = max(0, n)
    call this % declare(name)

  end function create_index

  pure integer function index_size(this)

    class(index_set), intent(in) :: this

    index_size = this % n

  end function index_size

  pure integer function index_member(this, local_index)

    class(index_set), intent(in) :: this
    integer           , intent(in) :: local_index

    index_member = local_index

  end function index_member

  pure subroutine index_members(this, indices)

    class(index_set)  , intent(in)  :: this
    integer, allocatable, intent(out) :: indices(:)

    integer :: k

    allocate(indices(this % n))
    do k = 1, this % n
       indices(k) = k
    end do

  end subroutine index_members

  !===================================================================!
  ! Membership in one comparison - the primitive a relation
  ! signature leans on, never an enumeration and a search.
  !===================================================================!

  pure logical function index_has(this, member)

    class(index_set), intent(in) :: this
    integer           , intent(in) :: member

    index_has = (member >= 1) .and. (member <= this % n)

  end function index_has

  !===================================================================!
  ! Where a member stands - for the counted set, where it says. Zero
  ! for an outsider, so an indexed relation reads absence without a
  ! second question.
  !===================================================================!

  pure integer function index_local_index(this, member)

    class(index_set), intent(in) :: this
    integer           , intent(in) :: member

    if (member >= 1 .and. member <= this % n) then
       index_local_index = member
    else
       index_local_index = 0
    end if

  end function index_local_index

  !===================================================================!
  ! Declare a subobject: a name, the ambient domain, the chosen
  ! members. Refusals guard the inclusion:
  !
  !      an unsigned ambient       a subset needs a declared ambient
  !                                domain
  !      a member from elsewhere   a subset holds members of its
  !                                ambient domain only
  !
  ! and a member handed in twice is in the subset once, first
  ! appearance keeping its place - sets are sets everywhere.
  !===================================================================!

  type(subset) function create_subset(name, ambient, members) &
       & result(this)

    character(len=*) , intent(in) :: name
    class(set), intent(in) :: ambient
    integer          , intent(in) :: members(:)

    integer :: keep(size(members))
    integer :: j, nkept

    if (.not. ambient % equals(ambient)) then
       error stop 'graph_set: a subset needs a declared ambient domain'
    end if

    nkept = 0
    do j = 1, size(members)
       if (.not. ambient % has(members(j))) then
          error stop 'graph_set: a subset holds members of its ambient domain only'
       end if
       if (.not. any(keep(1:nkept) == members(j))) then
          nkept       = nkept + 1
          keep(nkept) = members(j)
       end if
    end do

    allocate(this % host, source=ambient)
    this % roll = keep(1:nkept)
    call this % declare(name)

  end function create_subset

  !===================================================================!
  ! The ambient domain, as a copy - which is to say, as the same
  ! declared domain the subset was carved from.
  !===================================================================!

  function ambient(this) result(host)

    class(subset), intent(in)  :: this
    class(set), allocatable :: host

    allocate(host, source=this % host)

  end function ambient

  pure integer function subset_size(this)

    class(subset), intent(in) :: this

    subset_size = size(this % roll)

  end function subset_size

  pure integer function subset_member(this, local_index)

    class(subset), intent(in) :: this
    integer          , intent(in) :: local_index

    subset_member = this % roll(local_index)

  end function subset_member

  pure subroutine subset_members(this, indices)

    class(subset)   , intent(in)  :: this
    integer, allocatable, intent(out) :: indices(:)

    indices = this % roll

  end subroutine subset_members

  pure logical function subset_has(this, member)

    class(subset), intent(in) :: this
    integer          , intent(in) :: member

    subset_has = any(this % roll == member)

  end function subset_has

  pure integer function subset_local_index(this, member)

    class(subset), intent(in) :: this
    integer          , intent(in) :: member

    integer :: k

    subset_local_index = 0
    do k = 1, size(this % roll)
       if (this % roll(k) == member) then
          subset_local_index = k
          return
       end if
    end do

  end function subset_local_index

  !===================================================================!
  ! A subset precedes itself, its ambient, and everything its
  ! ambient precedes - the chain walked to any depth.
  !===================================================================!

  pure logical function subset_is_subobject_of(this, ancestor)

    class(subset) , intent(in) :: this
    class(set) , intent(in) :: ancestor

    subset_is_subobject_of = this % equals(ancestor)
    if (.not. subset_is_subobject_of) then
       subset_is_subobject_of = this % host % is_subobject_of(ancestor)
    end if

  end function subset_is_subobject_of

end module graph_set
