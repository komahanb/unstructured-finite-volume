!=====================================================================!
! LEVEL 1 OF THE NEW TOWER . THE RELATIONS
!
! The second object of the relation-centered tower (AGENTS.md): a
! named finite-arity subset of a cartesian product,
!
!      R  <=  A_1 x A_2 x ... x A_k ,        k >= 1
!
! A relation has identity, arity, an ordered signature, and a
! membership law - and it is FIRST-CLASS: constructible, queryable
! and testable with no graph anywhere in sight. A graph, when it
! arrives on a higher level, will CONTAIN relations; a relation
! never needs a graph in order to exist.
!
!                          THE SIGNATURE
!
! Each slot of the signature holds one set - ANY concretion of
! set, and different concretions may stand in different domains
! of one relation. The signature holds copies, and a copy IS the
! declared domain - so two relations built over one set answer
! equals across each other's domains, and no relation ever assumes
! its domains are owned by one graph. A slot may repeat a set:
!
!      R_CC  <=  cells x cells         adjacency
!      R_CF  <=  cells x faces         incidence
!
! adjacency and incidence are interpretations of the signature, not
! separate primitives. Higher arity is a slot, not a special case:
!
!      R_end <=  edges x vertices x roles
!
! holds (e, v, tail) and (e, v, head) for an interior edge and one
! lone (e, v, tail) for a boundary face - no imaginary far-side
! member, the same wall the old grammar drew with a headless edge.
!
!                         A SET, NOT A BAG
!
! A relation is a set of tuples: no tuple is in it twice. The
! constructor collapses duplicate columns to the first occurrence,
! order kept, so num_tuples, tuples and has all answer set
! semantics and nothing else. Multiplicity that MEANS something -
! two parallel edges between one pair of cells - is already carried
! by distinct members of an edge domain; the day counted repetition
! itself is the mathematics, that is a distinct multirelation
! abstraction, not a quiet flag here.
!
!                       WHAT IS VALIDATED
!
! Construction refuses, loudly: a tuple table whose row count is not
! the arity; a signature slot that was never declared; a tuple
! component its slot's set does not hold (checked through the
! set's own has - membership is a primitive, never an
! enumeration and a search). After construction the relation is
! immutable, so the laws hold for life.
!
!                     CAPABILITY, NOT FICTION
!
! This stored relation carries NO per-slot index yet: has() is a
! linear scan over the tuple table, and the constructor's duplicate
! collapse is quadratic in the tuple count - both honest, both
! stated (AGENTS.md 53), both construction-or-query costs no hot
! loop should sit on. The indexed, CSR-backed binary specialization
! is the next phase's business; nothing here pretends to be
! O(degree).
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module graph_relation

  use graph_identity, only : token, mint_token
  use graph_set , only : set, index_set

  implicit none

  private
  public :: relation, stored_relation, declared_domain

  !===================================================================!
  ! The abstract relation: identity, arity, ordered signature,
  ! membership, and tuples as non-hot generic access (AGENTS.md
  ! 5.1). The tuple convention everywhere: a table t(arity, count),
  ! column j the j-th tuple.
  !===================================================================!

  type, abstract :: relation

     type(token)                  , private :: identity
     character(len=:), allocatable, private :: label

   contains

     !----------------------------------------------------------------!
     ! The structural questions, deferred to each concretion.
     !----------------------------------------------------------------!

     procedure(relation_arity_interface)  , deferred :: arity
     procedure(relation_domain_interface) , deferred :: domain
     procedure(relation_has_interface)    , deferred :: has
     procedure(relation_count_interface)  , deferred :: num_tuples
     procedure(relation_tuples_interface) , deferred :: tuples

     !----------------------------------------------------------------!
     ! Identity, answered once for every concretion - the one token
     ! law of graph_identity: sign once, refuse twice, copies carry
     ! the stamp, the undeclared equal nothing.
     !----------------------------------------------------------------!

     procedure :: declare
     procedure :: id
     procedure :: equals

     !----------------------------------------------------------------!
     ! Self-containment, FAILING CLOSED. A relation is assumed a
     ! borrower - unsafe to own - until a concretion CLAIMS to be
     ! materialized: whole unto itself, safe to copy and to own.
     ! Stored citizens claim it; views never do, and a view an
     ! author forgets to mark stays unownable by default, which is
     ! the only safe direction for the mistake.
     !----------------------------------------------------------------!

     procedure :: materialized

     !----------------------------------------------------------------!
     ! Metadata, not mathematics.
     !----------------------------------------------------------------!

     procedure :: name

  end type relation

  abstract interface

     pure integer function relation_arity_interface(this)
       import relation
       class(relation), intent(in) :: this
     end function relation_arity_interface

     function relation_domain_interface(this, slot_index) result(domain)
       import relation, set
       class(relation), intent(in)    :: this
       integer        , intent(in)    :: slot_index
       class(set), allocatable :: domain
     end function relation_domain_interface

     pure logical function relation_has_interface(this, tuple)
       import relation
       class(relation), intent(in) :: this
       integer        , intent(in) :: tuple(:)
     end function relation_has_interface

     pure integer function relation_count_interface(this)
       import relation
       class(relation), intent(in) :: this
     end function relation_count_interface

     pure subroutine relation_tuples_interface(this, table)
       import relation
       class(relation)     , intent(in)  :: this
       integer, allocatable, intent(out) :: table(:,:)
     end subroutine relation_tuples_interface

  end interface

  !===================================================================!
  ! ONE DECLARED DOMAIN, standing at one position of the ordered
  ! signature. The wrapper exists because a Fortran array carries
  ! one dynamic type: an array of domains is how different set
  ! realizations stand side by side in one signature.
  !
  ! The POSITION is the slot, and slot_index names it. THIS is what
  ! occupies the position, and it is a domain - so a signature is an
  ! ordered tuple of declared domains, which is what the mathematics
  ! says it is. Build one with declared_domain(A).
  !
  ! A domain MAY REPEAT across positions: R <= A x A is ordinary,
  ! and that is the law separating a signature from the collection S
  ! of a graph, which holds each domain once.
  !===================================================================!

  type :: declared_domain

     class(set), allocatable :: set

  end type declared_domain

  interface declared_domain
     module procedure create_declared_domain
  end interface declared_domain

  !===================================================================!
  ! The stored relation: the deduplicated tuple table, the signature
  ! as set copies. The first inhabitant of the contract, and the
  ! validation gate of the level.
  !===================================================================!

  type, extends(relation) :: stored_relation

     type(declared_domain), allocatable, private :: signature(:)
     integer              , allocatable, private :: tuple_table(:,:)

   contains

     procedure :: arity        => stored_arity
     procedure :: domain       => stored_domain
     procedure :: has          => stored_has
     procedure :: num_tuples   => stored_num_tuples
     procedure :: tuples       => stored_tuples
     procedure :: materialized => stored_materialized

  end type stored_relation

  interface stored_relation
     module procedure create_stored
     module procedure create_stored_index
  end interface stored_relation

contains

  !===================================================================!
  ! The identity block, one law with the sets'.
  !===================================================================!

  subroutine declare(this, name)

    class(relation) , intent(inout)        :: this
    character(len=*), intent(in), optional :: name

    if (this % identity % declared()) then
       error stop 'graph_relation: a relation never signs twice'
    end if

    this % identity = mint_token()
    if (present(name)) this % label = name

  end subroutine declare

  !===================================================================!
  ! id answers the whole opaque token - the identity itself, honest
  ! across images, never a bare local integer.
  !===================================================================!

  pure type(token) function id(this)

    class(relation), intent(in) :: this

    id = this % identity

  end function id

  pure logical function equals(this, other)

    class(relation), intent(in) :: this
    class(relation), intent(in) :: other

    equals = this % identity % matches(other % identity)

  end function equals

  pure logical function materialized(this)

    class(relation), intent(in) :: this

    materialized = .false.

  end function materialized

  pure logical function stored_materialized(this)

    class(stored_relation), intent(in) :: this

    stored_materialized = .true.

  end function stored_materialized

  function name(this)

    class(relation), intent(in)   :: this
    character(len=:), allocatable :: name

    if (allocated(this % label)) then
       name = this % label
    else
       name = ''
    end if

  end function name

  !===================================================================!
  ! Wrap one set - any concretion - as a signature slot.
  !===================================================================!

  type(declared_domain) function create_declared_domain(domain) result(this)

    class(set), intent(in) :: domain

    allocate(this % set, source=domain)

  end function create_declared_domain

  !===================================================================!
  ! Declare a stored relation: a name, the ordered domains, and the
  ! tuple table, one column per tuple. Refusals, in the order they
  ! are checked:
  !
  !     an empty signature            k >= 1 is the definition
  !     an undeclared set         a signature refers to declared
  !                                   domains only
  !     a row count off the arity     each tuple has exactly k parts
  !     a member no slot holds        domain validity, through has
  !
  ! Duplicate tuple columns then collapse to the first occurrence,
  ! order kept: a relation is a set.
  !===================================================================!

  type(stored_relation) function create_stored(name, domains, table) &
       & result(this)

    character(len=*), intent(in) :: name
    type(declared_domain)      , intent(in) :: domains(:)
    integer         , intent(in) :: table(:,:)

    integer, allocatable :: kept(:)
    integer              :: k, j, i, nkept
    logical              :: fresh

    if (size(domains) < 1) then
       error stop 'graph_relation: a relation relates at least one domain'
    end if

    do k = 1, size(domains)
       if (.not. allocated(domains(k) % set)) then
          error stop 'graph_relation: a signature refers to declared domains only'
       end if
       if (.not. domains(k) % set % equals(domains(k) % set)) then
          error stop 'graph_relation: a signature refers to declared domains only'
       end if
    end do

    if (size(table, 1) /= size(domains)) then
       error stop 'graph_relation: each tuple has exactly one part per slot'
    end if

    do j = 1, size(table, 2)
       do k = 1, size(domains)
          if (.not. domains(k) % set % has(table(k, j))) then
             error stop 'graph_relation: a tuple names a member its domain does not hold'
          end if
       end do
    end do

    ! A set, not a bag: keep each tuple's first appearance, in order.
    allocate(kept(size(table, 2)))
    nkept = 0
    do j = 1, size(table, 2)
       fresh = .true.
       do i = 1, nkept
          if (all(table(:, kept(i)) == table(:, j))) then
             fresh = .false.
             exit
          end if
       end do
       if (fresh) then
          nkept       = nkept + 1
          kept(nkept) = j
       end if
    end do

    allocate(this % signature(size(domains)))
    do k = 1, size(domains)
       allocate(this % signature(k) % set, source=domains(k) % set)
    end do

    allocate(this % tuple_table(size(domains), nkept))
    do i = 1, nkept
       this % tuple_table(:, i) = table(:, kept(i))
    end do

    call this % declare(name)

  end function create_stored

  !===================================================================!
  ! The counted convenience: the common case, a signature of counted
  ! domains, wrapped into domains and handed to the one gate above.
  !===================================================================!

  type(stored_relation) function create_stored_index(name, domains, &
       & table) result(this)

    character(len=*)  , intent(in) :: name
    type(index_set)   , intent(in) :: domains(:)
    integer           , intent(in) :: table(:,:)

    type(declared_domain), allocatable :: signature(:)
    integer                            :: k

    allocate(signature(size(domains)))
    do k = 1, size(domains)
       allocate(signature(k) % set, source=domains(k))
    end do

    this = create_stored(name, signature, table)

  end function create_stored_index

  pure integer function stored_arity(this)

    class(stored_relation), intent(in) :: this

    stored_arity = size(this % signature)

  end function stored_arity

  !===================================================================!
  ! The slot's set, as a copy - which is to say, as the same
  ! declared domain.
  !===================================================================!

  function stored_domain(this, slot_index) result(domain)

    class(stored_relation), intent(in) :: this
    integer               , intent(in) :: slot_index
    class(set), allocatable     :: domain

    allocate(domain, source=this % signature(slot_index) % set)

  end function stored_domain

  !===================================================================!
  ! Membership by linear scan - stated capability, no fiction. The
  ! indexed lookup belongs to the binary specialization.
  !===================================================================!

  pure logical function stored_has(this, tuple)

    class(stored_relation), intent(in) :: this
    integer               , intent(in) :: tuple(:)

    integer :: j

    stored_has = .false.

    if (size(tuple) /= size(this % signature)) return

    do j = 1, size(this % tuple_table, 2)
       if (all(this % tuple_table(:, j) == tuple)) then
          stored_has = .true.
          return
       end if
    end do

  end function stored_has

  pure integer function stored_num_tuples(this)

    class(stored_relation), intent(in) :: this

    stored_num_tuples = size(this % tuple_table, 2)

  end function stored_num_tuples

  pure subroutine stored_tuples(this, table)

    class(stored_relation), intent(in)  :: this
    integer, allocatable  , intent(out) :: table(:,:)

    allocate(table, source=this % tuple_table)

  end subroutine stored_tuples

end module graph_relation
