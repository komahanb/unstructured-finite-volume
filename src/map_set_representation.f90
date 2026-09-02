!=====================================================================!
! SET REPRESENTATION
!
! How a finite set's members are stored and enumerated HERE. Three
! primitive operations, and nothing else:
!
!     num_members    how many
!     member(k)      the member at position k
!     local_index(v) the position of v, 0 for a non-member
!
! bound by the enumeration laws
!
!     member(local_index(v)) = v      for every member v
!     local_index(member(k)) = k      for k = 1 .. num_members
!
! which force enumeration to be injective: a representation lists each
! member once. Two further operations are consequences of those laws
! and are implemented once, on the abstract type, for every concretion:
!
!     has(v)    = (local_index(v) /= 0)
!     members   = [ member(k), k = 1 .. num_members ]
!
! A concretion that overrides either must return the same; in
! particular a non-member's local_index is 0 and no other sentinel.
!
!                       NO SEMANTIC IDENTITY
!
! A representation stores NO token, NO declare, NO id, NO same_as and
! NO name. WHICH set this describes is not the representation's
! concern - that is determined by a type(graph), and by nothing else.
! Two representations
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
! one integer; a million-node semantic task graph may validly
! contain 10^6 graphs. Semantic complexity and extensional
! cardinality are different quantities.
!
!                     THE NAMESPACE CONVENTION
!
! Fortran forbids a module and its type sharing one name. The
! convention, recorded once and applied throughout: the MODULE has
! the systematic name graph_<noun>_<role>; the TYPES inside omit the
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
  ! The contract: the three primitives deferred, the two theorems
  ! concrete, identity nowhere.
  !===================================================================!

  type, abstract :: set_representation

   contains

     procedure(representation_num_members_interface), deferred :: num_members
     procedure(representation_member_interface)     , deferred :: member
     procedure(representation_local_index_interface), deferred :: local_index

     procedure :: members => representation_members
     procedure :: has     => representation_has

  end type set_representation

  abstract interface

     pure integer function representation_num_members_interface(this)
       import set_representation
       class(set_representation), intent(in) :: this
     end function representation_num_members_interface

     pure integer function representation_member_interface(this, position)
       import set_representation
       class(set_representation), intent(in) :: this
       integer                  , intent(in) :: position
     end function representation_member_interface

     pure integer function representation_local_index_interface(this, value)
       import set_representation
       class(set_representation), intent(in) :: this
       integer                  , intent(in) :: value
     end function representation_local_index_interface

  end interface

  !===================================================================!
  ! The counted representation: members 1..n, stored as n. O(1),
  ! whatever n is: every primitive is a comparison or an identity.
  !===================================================================!

  type, extends(set_representation) :: counted_set_representation

     integer, private :: n = 0

   contains

     procedure :: num_members => counted_num_members
     procedure :: member      => counted_member
     procedure :: local_index => counted_local_index

  end type counted_set_representation

  interface counted_set_representation
     module procedure create_counted
  end interface counted_set_representation

  !===================================================================!
  ! The listed representation: an explicit array of member values, in
  ! declaration order. The representation describes a set with or
  ! without any declared ambient - the embedding is not stored here.
  !===================================================================!

  type, extends(set_representation) :: listed_set_representation

     integer, allocatable, private :: listed(:)

   contains

     procedure :: num_members => listed_num_members
     procedure :: member      => listed_member
     procedure :: local_index => listed_local_index

  end type listed_set_representation

  interface listed_set_representation
     module procedure create_listed
  end interface listed_set_representation

contains

  !===================================================================!
  ! The two theorems, for every concretion: membership is a nonzero
  ! position, and the members are the enumeration.
  !===================================================================!

  pure logical function representation_has(this, value)

    class(set_representation), intent(in) :: this
    integer                  , intent(in) :: value

    representation_has = this % local_index(value) /= 0

  end function representation_has

  pure subroutine representation_members(this, values)

    class(set_representation), intent(in)  :: this
    integer, allocatable     , intent(out) :: values(:)

    integer :: k

    values = [(this % member(k), k = 1, this % num_members())]

  end subroutine representation_members

  !===================================================================!
  ! Counted: n, clamped at zero. The empty set is a set.
  !===================================================================!

  pure type(counted_set_representation) function create_counted(n) &
       & result(this)

    integer, intent(in) :: n

    this % n = max(0, n)

  end function create_counted

  pure integer function counted_num_members(this)

    class(counted_set_representation), intent(in) :: this

    counted_num_members = this % n

  end function counted_num_members

  pure integer function counted_member(this, position)

    class(counted_set_representation), intent(in) :: this
    integer                          , intent(in) :: position

    counted_member = position

  end function counted_member

  pure integer function counted_local_index(this, value)

    class(counted_set_representation), intent(in) :: this
    integer                          , intent(in) :: value

    counted_local_index = 0
    if (value >= 1 .and. value <= this % n) counted_local_index = value

  end function counted_local_index

  !===================================================================!
  ! Listed: the values passed in, each retained once, the first
  ! occurrence retaining its position - the enumeration law needs
  ! injectivity, and declaration order is the representation's own.
  !===================================================================!

  pure type(listed_set_representation) function create_listed(values) &
       & result(this)

    integer, intent(in) :: values(:)

    integer :: retained(size(values))
    integer :: j, num_retained

    num_retained = 0
    do j = 1, size(values)
       if (.not. any(retained(1:num_retained) == values(j))) then
          num_retained           = num_retained + 1
          retained(num_retained) = values(j)
       end if
    end do

    this % listed = retained(1:num_retained)

  end function create_listed

  pure integer function listed_num_members(this)

    class(listed_set_representation), intent(in) :: this

    if (allocated(this % listed)) then
       listed_num_members = size(this % listed)
    else
       listed_num_members = 0
    end if

  end function listed_num_members

  pure integer function listed_member(this, position)

    class(listed_set_representation), intent(in) :: this
    integer                         , intent(in) :: position

    listed_member = this % listed(position)

  end function listed_member

  pure integer function listed_local_index(this, value)

    class(listed_set_representation), intent(in) :: this
    integer                         , intent(in) :: value

    integer :: k

    listed_local_index = 0
    if (.not. allocated(this % listed)) return

    do k = 1, size(this % listed)
       if (this % listed(k) == value) then
          listed_local_index = k
          return
       end if
    end do

  end function listed_local_index

end module map_set_representation
