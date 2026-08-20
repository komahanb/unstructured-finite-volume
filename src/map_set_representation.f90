!=====================================================================!
! SET REPRESENTATION
!
! How a finite set's members are stored and enumerated HERE. Five
! questions, and nothing else:
!
!     size          how many
!     member(k)     which one stands at position k
!     members       all of them
!     has(v)        is v one of them
!     local_index(v) where does v stand, 0 for an outsider
!
! bound by the enumeration laws
!
!     member(local_index(v)) = v      for every member v
!     local_index(member(k)) = k      for k = 1 .. size
!
! which force enumeration to be injective: a representation lists each
! member once.
!
!                       NO SEMANTIC IDENTITY
!
! A representation carries NO token, NO declare, NO id, NO same_as and
! NO name. WHICH set this describes is not its question - that is
! answered by a type(graph), and by nothing else. Two representations
! with equal extensions are interchangeable descriptions; two graphs
! with equal extensions are two sets.
!
! This is the whole of the split: identity above, storage here.
!
!                        SCALE, AND WHY
!
!     |G|                    = O(N_semantic)
!     |representation|       = O(N_extent)
!
! A bulk extension must not require O(N_extent) semantic graph objects
! merely in order to exist. A counted representation of 10^9 members is
! one integer; a genuinely million-node semantic task graph may
! legitimately hold 10^6 graphs. Semantic complexity and extensional
! cardinality are different quantities.
!
!                     THE NAMESPACE CONVENTION
!
! Fortran forbids a module and its type sharing one name. The
! convention, recorded once and applied throughout: the MODULE carries
! the systematic name graph_<noun>_<role>; the TYPES inside drop the
! graph_ prefix. map_set/set_map is the precedent.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module map_set_representation

  implicit none

  private
  public :: set_representation
  public :: counted_set_representation, listed_set_representation

  !===================================================================!
  ! The contract. Deferred everywhere, identity nowhere.
  !===================================================================!

  type, abstract :: set_representation

   contains

     procedure(representation_size_interface)       , deferred :: size
     procedure(representation_member_interface)     , deferred :: member
     procedure(representation_members_interface)    , deferred :: members
     procedure(representation_has_interface)        , deferred :: has
     procedure(representation_local_index_interface), deferred :: local_index

  end type set_representation

  abstract interface

     pure integer function representation_size_interface(this)
       import set_representation
       class(set_representation), intent(in) :: this
     end function representation_size_interface

     pure integer function representation_member_interface(this, position)
       import set_representation
       class(set_representation), intent(in) :: this
       integer                  , intent(in) :: position
     end function representation_member_interface

     pure subroutine representation_members_interface(this, values)
       import set_representation
       class(set_representation), intent(in)  :: this
       integer, allocatable     , intent(out) :: values(:)
     end subroutine representation_members_interface

     pure logical function representation_has_interface(this, value)
       import set_representation
       class(set_representation), intent(in) :: this
       integer                  , intent(in) :: value
     end function representation_has_interface

     pure integer function representation_local_index_interface(this, value)
       import set_representation
       class(set_representation), intent(in) :: this
       integer                  , intent(in) :: value
     end function representation_local_index_interface

  end interface

  !===================================================================!
  ! The counted representation: members 1..n, stored as n. O(1),
  ! whatever n is, and every query is one comparison.
  !===================================================================!

  type, extends(set_representation) :: counted_set_representation

     integer, private :: n = 0

   contains

     procedure :: size        => counted_size
     procedure :: member      => counted_member
     procedure :: members     => counted_members
     procedure :: has         => counted_has
     procedure :: local_index => counted_local_index

  end type counted_set_representation

  interface counted_set_representation
     module procedure create_counted
  end interface counted_set_representation

  !===================================================================!
  ! The listed representation: an explicit roll of member values, in
  ! declaration order. It describes a set with or without any declared
  ! ambient - the embedding is not its business.
  !===================================================================!

  type, extends(set_representation) :: listed_set_representation

     integer, allocatable, private :: roll(:)

   contains

     procedure :: size        => listed_size
     procedure :: member      => listed_member
     procedure :: members     => listed_members
     procedure :: has         => listed_has
     procedure :: local_index => listed_local_index

  end type listed_set_representation

  interface listed_set_representation
     module procedure create_listed
  end interface listed_set_representation

contains

  !===================================================================!
  ! Counted: n, clamped at zero. The empty set is a set.
  !===================================================================!

  pure type(counted_set_representation) function create_counted(n) &
       & result(this)

    integer, intent(in) :: n

    this % n = max(0, n)

  end function create_counted

  pure integer function counted_size(this)

    class(counted_set_representation), intent(in) :: this

    counted_size = this % n

  end function counted_size

  pure integer function counted_member(this, position)

    class(counted_set_representation), intent(in) :: this
    integer                          , intent(in) :: position

    counted_member = position

  end function counted_member

  pure subroutine counted_members(this, values)

    class(counted_set_representation), intent(in)  :: this
    integer, allocatable             , intent(out) :: values(:)

    integer :: k

    allocate(values(this % n))
    do k = 1, this % n
       values(k) = k
    end do

  end subroutine counted_members

  pure logical function counted_has(this, value)

    class(counted_set_representation), intent(in) :: this
    integer                          , intent(in) :: value

    counted_has = (value >= 1) .and. (value <= this % n)

  end function counted_has

  pure integer function counted_local_index(this, value)

    class(counted_set_representation), intent(in) :: this
    integer                          , intent(in) :: value

    if (this % has(value)) then
       counted_local_index = value
    else
       counted_local_index = 0
    end if

  end function counted_local_index

  !===================================================================!
  ! Listed: the values handed in, each kept once, first appearance
  ! keeping its place - the enumeration law needs injectivity, and
  ! declaration order is the representation's own.
  !===================================================================!

  pure type(listed_set_representation) function create_listed(values) &
       & result(this)

    integer, intent(in) :: values(:)

    integer :: keep(size(values))
    integer :: j, nkept

    nkept = 0
    do j = 1, size(values)
       if (.not. any(keep(1:nkept) == values(j))) then
          nkept       = nkept + 1
          keep(nkept) = values(j)
       end if
    end do

    this % roll = keep(1:nkept)

  end function create_listed

  pure integer function listed_size(this)

    class(listed_set_representation), intent(in) :: this

    if (allocated(this % roll)) then
       listed_size = size(this % roll)
    else
       listed_size = 0
    end if

  end function listed_size

  pure integer function listed_member(this, position)

    class(listed_set_representation), intent(in) :: this
    integer                         , intent(in) :: position

    listed_member = this % roll(position)

  end function listed_member

  pure subroutine listed_members(this, values)

    class(listed_set_representation), intent(in)  :: this
    integer, allocatable            , intent(out) :: values(:)

    if (allocated(this % roll)) then
       values = this % roll
    else
       allocate(values(0))
    end if

  end subroutine listed_members

  pure logical function listed_has(this, value)

    class(listed_set_representation), intent(in) :: this
    integer                         , intent(in) :: value

    listed_has = .false.
    if (allocated(this % roll)) listed_has = any(this % roll == value)

  end function listed_has

  pure integer function listed_local_index(this, value)

    class(listed_set_representation), intent(in) :: this
    integer                         , intent(in) :: value

    integer :: k

    listed_local_index = 0
    if (.not. allocated(this % roll)) return

    do k = 1, size(this % roll)
       if (this % roll(k) == value) then
          listed_local_index = k
          return
       end if
    end do

  end function listed_local_index

end module map_set_representation
