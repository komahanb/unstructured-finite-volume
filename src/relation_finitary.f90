!=====================================================================!
! LEVEL 1 OF THE NEW TOWER . THE RELATIONS
!
! The second object of the relation-centred tower (AGENTS.md): a
! named finite-arity subset of a cartesian product,
!
!      P  <=  A_1 x A_2 x ... x A_k ,        k >= 1
!
! A relation has identity, arity, an ordered signature, and a
! membership law - and it is FIRST-CLASS: constructible, queryable
! and testable with no graph required. A graph, defined on a higher
! level, CONTAINS relations; a relation does not require a graph in
! order to exist.
!
!                          THE SIGNATURE
!
! An ORDERED SEQUENCE OF SET GRAPHS, and nothing else:
!
!      sig(P) = (A_1, ..., A_k)
!
! It is small control data - k identities - so it is stored
! contiguously, as an array of graph values. That IS the ordered graph
! sequence; view_sequence represents the same mathematics as a linked
! chain of cells and states so: indexed access there is O(k) and
! returns a POINTER into cells the owner must retain. A signature is
! indexed frequently and copied freely, so it is compiled to the
! contiguous form. No non-owning reference, no chain of cells to own,
! and the map law of the identity maps is preserved: a signature owns
! its identities by value.
!
! THE SLOT WRAPPER IS REMOVED. It existed because an array stores one
! dynamic type and the previous domain types were a class hierarchy;
! every domain is a type(graph) now, so the wrapper contained nothing.
! It stored no mathematical information and was removed with the type
! that required it.
!
! The signature stores copies, and a copy IS the declared domain - so
! two relations built over one domain return true from same_as across
! each other's positions, and no relation ever assumes its domains are
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
! contains (e, v, tail) and (e, v, head) for an interior edge and one
! single (e, v, tail) for a boundary face - no fictitious far-side
! member, the same boundary the previous grammar represented with a
! headless edge.
!
!                       A SET, NOT A MULTISET
!
! A relation is a set of tuples: no tuple is in it twice. The
! constructor collapses duplicate columns to the first occurrence,
! order preserved, so num_tuples, tuples and has all follow set
! semantics and nothing else. Multiplicity that denotes something -
! two parallel edges between one pair of cells - is already
! represented by distinct members of an edge domain; when counted
! repetition itself is the mathematics required, that is a distinct
! multirelation abstraction, not an undocumented flag here.
!
!                       WHAT IS VALIDATED
!
! Construction stops the program on: a tuple table whose row count is
! not the arity; a signature position that was never declared; a
! tuple component its domain does not contain.
!
! That last check is a MEMBERSHIP query, and membership belongs to
! a representation - so construction takes the caller's set map and
! queries it. The map is used HERE and nowhere else: it is a
! compilation input, not a stored dependency, and the relation
! retains no reference to it. After construction the relation is
! immutable, so the laws hold for its whole lifetime without it.
!
!                       STATED CAPABILITY
!
! This stored relation stores NO per-domain index yet: has() is a
! linear scan over the tuple table, and the constructor's duplicate
! collapse is quadratic in the tuple count - both exact, both
! stated (AGENTS.md 53), both construction-or-query costs no hot
! loop should depend on. The indexed, CSR-backed binary
! specialisation is a later phase; nothing here is O(degree).
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module relation_finitary

  use token_identity, only : token, next_token
  use graph_fractal , only : graph
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
     type(graph)     , allocatable, private :: signature(:)

   contains

     !----------------------------------------------------------------!
     ! The membership queries, deferred to each concrete type.
     !----------------------------------------------------------------!

     procedure(relation_has_interface)    , deferred :: has
     procedure(relation_count_interface)  , deferred :: num_tuples
     procedure(relation_tuples_interface) , deferred :: tuples

     !----------------------------------------------------------------!
     ! Identity and signature, implemented once for every concrete
     ! type - the one token law of token_identity: declare once,
     ! reject a second declaration, copies retain the token, the
     ! undeclared equal nothing - and the ordered domains, checked
     ! declared at the one declaration.
     !----------------------------------------------------------------!

     procedure :: declare
     procedure :: id
     procedure :: same_as
     procedure :: arity
     procedure :: domain

     !----------------------------------------------------------------!
     ! Self-containment, DEFAULTING TO FALSE. A relation is assumed
     ! to store non-owning references - unsafe to own - until a concrete
     ! type DECLARES itself materialized: self-contained, safe to copy
     ! and to own. Stored types declare it; views never do, and a view
     ! an author omits to mark stays unownable by default, which is
     ! the only safe direction for the omission.
     !----------------------------------------------------------------!

     procedure :: materialized

     !----------------------------------------------------------------!
     ! Metadata, not mathematics.
     !----------------------------------------------------------------!

     procedure :: name

  end type relation

  abstract interface

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
  ! The stored relation: the deduplicated tuple table over the
  ! signature the root stores. The first implementation of the
  ! contract, and the validation check of the level.
  !
  ! It stores NO representation. A generic table relation evaluates
  ! has() by scanning its own tuples, so nothing it does after
  ! construction is a membership query. The set map validates it at
  ! construction and is not referenced afterwards.
  !===================================================================!

  type, extends(relation) :: stored_relation

     integer, allocatable, private :: entry(:,:)

   contains

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
  ! The identity block, the same law as the graph types', and the
  ! signature. Invalid inputs: a second declaration; an empty
  ! signature, since k >= 1 is the definition; an undeclared domain,
  ! since a signature names declared sets (an undeclared token does
  ! not match itself).
  !===================================================================!

  subroutine declare(this, name, domains)

    class(relation) , intent(inout) :: this
    character(len=*), intent(in)    :: name
    type(graph)     , intent(in)    :: domains(:)

    integer :: k

    if (this % identity % declared()) then
       error stop 'relation_finitary: a relation is declared at most once'
    end if

    if (size(domains) < 1) then
       error stop 'relation_finitary: a relation relates at least one domain'
    end if

    do k = 1, size(domains)
       if (.not. domains(k) % same_as(domains(k))) then
          error stop 'relation_finitary: a signature refers to declared domains only'
       end if
    end do

    this % identity  = next_token()
    this % label     = name
    this % signature = domains

  end subroutine declare

  pure integer function arity(this)

    class(relation), intent(in) :: this

    arity = size(this % signature)

  end function arity

  !===================================================================!
  ! WHICH domain is at this position, by value - that is, the same
  ! declared domain. Not pure: a set graph contains a pointer
  ! component, so copying one out of an INTENT(IN) dummy is barred
  ! from a pure subprogram (F2018 C1594). This is a control query and
  ! no hot path calls it - the numbering the hot path requires is
  ! stored in the representation, never here.
  !===================================================================!

  type(graph) function domain(this, position)

    class(relation), intent(in) :: this
    integer        , intent(in) :: position

    domain = this % signature(position)

  end function domain

  !===================================================================!
  ! id returns the whole opaque token - the identity itself,
  ! consistent across images, never a bare local integer.
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
  ! table one column per tuple, and the set map that states what those
  ! domains contain. Invalid inputs, in the order they are checked:
  !
  !     an empty signature            k >= 1, checked by declare
  !     an undeclared domain          a signature names declared sets
  !     an undescribed domain         the map describes every domain
  !     a row count off the arity     each tuple has exactly k parts
  !     a member no domain contains   domain validity, through the map
  !
  ! The map is a COMPILATION INPUT. It is read here, to decide what may
  ! exist, and never stored - the finished relation stores identities
  ! and integers, so it copies freely and outlives the map that
  ! validated it.
  !
  ! Duplicate tuple columns then collapse to the first occurrence,
  ! order preserved: a relation is a set.
  !===================================================================!

  type(stored_relation) function create_stored(name, domains, table, sets) &
       & result(this)

    character(len=*), intent(in) :: name
    type(graph) , intent(in) :: domains(:)
    integer         , intent(in) :: table(:,:)
    type(set_map)   , intent(in) :: sets

    integer, allocatable :: kept(:)
    integer              :: k, j, i, nkept
    logical              :: first_occurrence

    call this % declare(name, domains)

    do k = 1, size(domains)
       if (.not. sets % describes(domains(k))) then
          error stop 'relation_finitary: a signature refers to described domains only'
       end if
    end do

    if (size(table, 1) /= size(domains)) then
       error stop 'relation_finitary: each tuple has exactly one part per domain'
    end if

    do j = 1, size(table, 2)
       do k = 1, size(domains)
          if (.not. sets % has(domains(k), table(k, j))) then
             error stop 'relation_finitary: a tuple names a member its domain does not hold'
          end if
       end do
    end do

    ! A set, not a multiset: retain each tuple's first occurrence, in order.
    allocate(kept(size(table, 2)))
    nkept = 0
    do j = 1, size(table, 2)
       first_occurrence = .true.
       do i = 1, nkept
          if (all(table(:, kept(i)) == table(:, j))) then
             first_occurrence = .false.
             exit
          end if
       end do
       if (first_occurrence) then
          nkept       = nkept + 1
          kept(nkept) = j
       end if
    end do

    allocate(this % entry(size(domains), nkept))
    do i = 1, nkept
       this % entry(:, i) = table(:, kept(i))
    end do

  end function create_stored

  !===================================================================!
  ! Membership by linear scan, as stated. The indexed lookup belongs
  ! to the binary specialisation.
  !===================================================================!

  pure logical function stored_has(this, tuple)

    class(stored_relation), intent(in) :: this
    integer               , intent(in) :: tuple(:)

    integer :: j

    stored_has = .false.

    if (size(tuple) /= this % arity()) return

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
