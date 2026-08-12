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
! Each slot of the signature names one carrier, by the carrier's own
! structural identity. The signature holds copies, and a copy IS the
! declared domain - so two relations built over one carrier answer
! same_as across each other's slots, and no relation ever assumes
! its domains are owned by one graph. A slot may repeat a carrier:
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
!                       WHAT IS VALIDATED
!
! Construction refuses, loudly: a tuple table whose row count is not
! the arity; a signature slot that was never declared; a tuple
! component its slot's carrier does not hold (checked through the
! carrier's own has - membership is a primitive, never an
! enumeration and a search). After construction the relation is
! immutable, so the laws hold for life.
!
!                     CAPABILITY, NOT FICTION
!
! This stored relation keeps its tuples verbatim: order kept,
! duplicates kept - a deliberate departure from pure set semantics,
! stated here (AGENTS.md 11.5); has() answers membership regardless
! of multiplicity. It carries NO per-slot index yet: has() is a
! linear scan over the tuple table, honest and documented
! (AGENTS.md 53). The indexed, CSR-backed binary specialization is
! the next phase's business; nothing here pretends to be O(degree).
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module graph_relation

  use graph_carrier, only : member_set, counted_set, token, mint_token

  implicit none

  private
  public :: relation, stored_relation

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
     ! Identity, answered once for every concretion - the same
     ! opaque-token law the carriers keep: sign once, refuse twice,
     ! copies carry the stamp, the undeclared equal nothing.
     !----------------------------------------------------------------!

     procedure :: declare
     procedure :: id
     procedure :: same_as

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

     function relation_domain_interface(this, slot) result(domain)
       import relation, member_set
       class(relation), intent(in)   :: this
       integer        , intent(in)   :: slot
       class(member_set), allocatable :: domain
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
  ! One signature slot. The wrapper exists because Fortran arrays
  ! hold one dynamic type: each slot carries its own copy of a
  ! declared carrier, whatever concretion that carrier is.
  !===================================================================!

  type :: slot_domain

     class(member_set), allocatable :: carrier

  end type slot_domain

  !===================================================================!
  ! The stored relation: the tuple table verbatim, the signature as
  ! carrier copies. The first inhabitant of the contract, and the
  ! validation gate of the level.
  !===================================================================!

  type, extends(relation) :: stored_relation

     type(slot_domain), allocatable, private :: signature(:)
     integer          , allocatable, private :: entry(:,:)

   contains

     procedure :: arity      => stored_arity
     procedure :: domain     => stored_domain
     procedure :: has        => stored_has
     procedure :: num_tuples => stored_num_tuples
     procedure :: tuples     => stored_tuples

  end type stored_relation

  interface stored_relation
     module procedure create_stored
  end interface stored_relation

contains

  !===================================================================!
  ! The identity block, one law with the carriers'.
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

  pure integer function id(this)

    class(relation), intent(in) :: this

    id = this % identity % serial_number()

  end function id

  pure logical function same_as(this, other)

    class(relation), intent(in) :: this
    class(relation), intent(in) :: other

    same_as = this % identity % matches(other % identity)

  end function same_as

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
  ! Declare a stored relation: a name, the ordered carriers, and the
  ! tuple table, one column per tuple. Refusals, in the order they
  ! are checked:
  !
  !     an empty signature            k >= 1 is the definition
  !     an undeclared carrier         a signature refers to declared
  !                                   domains only
  !     a row count off the arity     each tuple has exactly k parts
  !     a member no slot holds        domain validity, through has
  !===================================================================!

  type(stored_relation) function create_stored(name, domains, table) &
       & result(this)

    character(len=*)  , intent(in) :: name
    type(counted_set) , intent(in) :: domains(:)
    integer           , intent(in) :: table(:,:)

    integer :: k, j

    if (size(domains) < 1) then
       error stop 'graph_relation: a relation relates at least one domain'
    end if

    do k = 1, size(domains)
       if (.not. domains(k) % same_as(domains(k))) then
          error stop 'graph_relation: a signature refers to declared domains only'
       end if
    end do

    if (size(table, 1) /= size(domains)) then
       error stop 'graph_relation: each tuple has exactly one part per slot'
    end if

    do j = 1, size(table, 2)
       do k = 1, size(domains)
          if (.not. domains(k) % has(table(k, j))) then
             error stop 'graph_relation: a tuple names a member its domain does not hold'
          end if
       end do
    end do

    allocate(this % signature(size(domains)))
    do k = 1, size(domains)
       allocate(this % signature(k) % carrier, source=domains(k))
    end do

    allocate(this % entry, source=table)

    call this % declare(name)

  end function create_stored

  pure integer function stored_arity(this)

    class(stored_relation), intent(in) :: this

    stored_arity = size(this % signature)

  end function stored_arity

  !===================================================================!
  ! The slot's carrier, as a copy - which is to say, as the same
  ! declared domain.
  !===================================================================!

  function stored_domain(this, slot) result(domain)

    class(stored_relation), intent(in) :: this
    integer               , intent(in) :: slot
    class(member_set), allocatable     :: domain

    allocate(domain, source=this % signature(slot) % carrier)

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

    do j = 1, size(this % entry, 2)
       if (all(this % entry(:, j) == tuple)) then
          stored_has = .true.
          return
       end if
    end do

  end function stored_has

  pure integer function stored_num_tuples(this)

    class(stored_relation), intent(in) :: this

    stored_num_tuples = size(this % entry, 2)

  end function stored_num_tuples

  pure subroutine stored_tuples(this, table)

    class(stored_relation), intent(in)  :: this
    integer, allocatable  , intent(out) :: table(:,:)

    allocate(table, source=this % entry)

  end subroutine stored_tuples

end module graph_relation
