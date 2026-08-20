!=====================================================================!
! LEVEL 1 OF THE NEW TOWER . THE RELATIONS
!
! The second object of the relation-centered tower (AGENTS.md): a
! named finite-arity subset of a cartesian product,
!
!      P  <=  A_1 x A_2 x ... x A_k ,        k >= 1
!
! A relation has identity, arity, an ordered signature, and a
! membership law - and it is FIRST-CLASS: constructible, queryable
! and testable with no graph anywhere in sight. A graph, when it
! arrives on a higher level, will CONTAIN relations; a relation
! never needs a graph in order to exist.
!
!                          THE SIGNATURE
!
! An ORDERED SEQUENCE OF SET GRAPHS, and nothing else:
!
!      sig(P) = (A_1, ..., A_k)
!
! It is small control-plane data - k identities - so it is stored
! contiguously, as an array of graph values. That IS the ordered graph
! sequence; view_sequence represents the same mathematics as a
! branch spine and says so itself: indexed access there is O(k) and
! hands back a POINTER into cells the holder must keep alive. A
! signature is indexed constantly and copied freely, so it takes the
! view's own advice and compiles to the contiguous form. No borrow, no
! spine to own, and the map law of the identity maps is kept: a
! signature owns its identities by value.
!
! THE SLOT IS GONE. It existed because an array carries one dynamic
! type and the old carriers were a class hierarchy; every domain is a
! type(graph) now, so the wrapper wrapped nothing. It carried no
! mathematical information and left with the type that needed it.
!
! The signature holds copies, and a copy IS the declared domain - so
! two relations built over one domain answer same_as across each
! other's positions, and no relation ever assumes its domains are
! owned by one graph. A position may repeat a domain:
!
!      P_CC  <=  cells x cells         adjacency
!      P_CF  <=  cells x faces         incidence
!
! adjacency and incidence are interpretations of the signature, not
! separate primitives. Higher arity is one more position, not a
! special case:
!
!      T_end <=  edges x vertices x roles
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
! the arity; a signature position that was never declared; a tuple
! component its domain does not hold.
!
! That last check is a MEMBERSHIP question, and membership belongs to
! a representation - so construction takes the caller's set map and
! asks it. The map is used HERE and nowhere else: it is a compilation
! input, not a stored dependency, and the relation keeps no reference
! to it. After construction the relation is immutable, so the laws
! hold for life without it.
!
!                     CAPABILITY, NOT FICTION
!
! This stored relation carries NO per-domain index yet: has() is a
! linear scan over the tuple table, and the constructor's duplicate
! collapse is quadratic in the tuple count - both honest, both
! stated (AGENTS.md 53), both construction-or-query costs no hot
! loop should sit on. The indexed, CSR-backed binary specialization
! is the next phase's business; nothing here pretends to be
! O(degree).
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module relation_finitary

  use token_identity, only : token, mint_token
  use graph_fractal , only : set_graph => graph
  use map_set , only : set_map

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
     ! Identity, answered once for every concretion - the one token
     ! law of token_identity: sign once, refuse twice, copies carry
     ! the stamp, the undeclared equal nothing.
     !----------------------------------------------------------------!

     procedure :: declare
     procedure :: id
     procedure :: same_as

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

     !--------------------------------------------------------------!
     ! WHICH domain stands at this position, by value. Not pure: a set
     ! graph carries a pointer component, so copying one out of an
     ! INTENT(IN) dummy is barred from a pure subprogram (F2018
     ! C1594). This is a control-plane question and no hot path asks
     ! it - the numbering the hot path needs lives in the
     ! representation, never here.
     !--------------------------------------------------------------!

     type(set_graph) function relation_domain_interface(this, position)
       import relation, set_graph
       class(relation), intent(in) :: this
       integer        , intent(in) :: position
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
  ! The stored relation: the deduplicated tuple table, and the
  ! signature as domain identities. The first inhabitant of the
  ! contract, and the validation gate of the level.
  !
  ! It keeps NO representation. A generic table relation answers has()
  ! by scanning its own tuples, and answers domain(k) by identity, so
  ! nothing it does after construction is a membership question. The
  ! set map validates it into existence and is then not its business.
  !===================================================================!

  type, extends(relation) :: stored_relation

     type(set_graph), allocatable, private :: signature(:)
     integer        , allocatable, private :: entry(:,:)

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
  end interface stored_relation

contains

  !===================================================================!
  ! The identity block, one law with the carriers'.
  !===================================================================!

  subroutine declare(this, name)

    class(relation) , intent(inout)        :: this
    character(len=*), intent(in), optional :: name

    if (this % identity % declared()) then
       error stop 'relation_finitary: a relation never signs twice'
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

  pure logical function same_as(this, other)

    class(relation), intent(in) :: this
    class(relation), intent(in) :: other

    same_as = this % identity % matches(other % identity)

  end function same_as

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
  ! Declare a stored relation: a name, the ordered domains, the tuple
  ! table one column per tuple, and the set map that says what those
  ! domains contain. Refusals, in the order they are checked:
  !
  !     an empty signature            k >= 1 is the definition
  !     an undeclared domain          a signature names declared sets
  !     a row count off the arity     each tuple has exactly k parts
  !     a member no domain holds      domain validity, through the map
  !
  ! The map is a COMPILATION INPUT. It is read here, to decide what may
  ! exist, and never stored - the finished relation holds identities
  ! and integers, so it copies freely and outlives the map that judged
  ! it.
  !
  ! Duplicate tuple columns then collapse to the first occurrence,
  ! order kept: a relation is a set.
  !===================================================================!

  type(stored_relation) function create_stored(name, domains, table, sets) &
       & result(this)

    character(len=*), intent(in) :: name
    type(set_graph) , intent(in) :: domains(:)
    integer         , intent(in) :: table(:,:)
    type(set_map)   , intent(in) :: sets

    integer, allocatable :: kept(:)
    integer              :: k, j, i, nkept
    logical              :: fresh

    if (size(domains) < 1) then
       error stop 'relation_finitary: a relation relates at least one domain'
    end if

    do k = 1, size(domains)
       if (.not. domains(k) % same_as(domains(k))) then
          error stop 'relation_finitary: a signature refers to declared domains only'
       end if
       if (.not. sets % describes(domains(k))) then
          error stop 'relation_finitary: a signature refers to described domains only'
       end if
    end do

    if (size(table, 1) /= size(domains)) then
       error stop 'relation_finitary: each tuple has exactly one part per domain'
    end if

    do j = 1, size(table, 2)
       do k = 1, size(domains)
          if (.not. sets % has_in(domains(k), table(k, j))) then
             error stop 'relation_finitary: a tuple names a member its domain does not hold'
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
       this % signature(k) = domains(k)
    end do

    allocate(this % entry(size(domains), nkept))
    do i = 1, nkept
       this % entry(:, i) = table(:, kept(i))
    end do

    call this % declare(name)

  end function create_stored

  pure integer function stored_arity(this)

    class(stored_relation), intent(in) :: this

    stored_arity = size(this % signature)

  end function stored_arity

  !===================================================================!
  ! The domain standing at this position, as a copy - which is to say,
  ! as the same declared domain.
  !===================================================================!

  type(set_graph) function stored_domain(this, position) result(domain)

    class(stored_relation), intent(in) :: this
    integer               , intent(in) :: position

    domain = this % signature(position)

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

end module relation_finitary
