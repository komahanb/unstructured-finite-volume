!=====================================================================!
! Module that states the set representation: how a finite set's
! members are stored and enumerated here. Three primitive operations,
! and nothing else:
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
!         position k                      member v
!      1  2  3  ...  n    -- member -->   values of the set
!                         <- local_index --
!
! which force the enumeration to be injective: a representation lists
! each member once. Two further operations are consequences of those
! laws and are implemented once, on the abstract type, for every
! concretion:
!
!     has(v)    = (local_index(v) /= 0)
!     members   = [ member(k), k = 1 .. num_members ]
!
! A concretion that overrides either must return the same; in
! particular a non-member's local_index is 0 and no other sentinel.
!
! No semantic identity. A representation stores no token, no declare,
! no id, no same_as and no name. Which set it describes is not the
! representation's concern: that is determined by a type(graph), and
! by nothing else. Two representations with equal extensions are
! interchangeable descriptions; two graphs with equal extensions are
! two sets. This is the whole of the split: identity above, storage
! here.
!
! Scale, and why:
!
!     |G|                    = O(N_semantic)
!     |representation|       = O(N_extent)
!
! A bulk extension must not require O(N_extent) semantic graph objects
! merely in order to exist. A counted representation of 10^9 members
! is one integer; a million-node semantic task graph may validly
! contain 10^6 graphs. Semantic complexity and extensional cardinality
! are different quantities.
!
! The namespace convention. Fortran forbids a module and its type
! sharing one name. The convention, recorded once and applied
! throughout: the module has the systematic name graph_<noun>_<role>;
! the types inside omit the graph_ prefix. map_set/set_map is the
! precedent.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module map_set_representation

  implicit none

  private

  public :: set_representation
  public :: counted_set_representation, listed_set_representation

  !-------------------------------------------------------------------!
  ! Type that models the contract: the three primitives deferred, the
  ! two theorems concrete, identity nowhere
  !-------------------------------------------------------------------!

  type, abstract :: set_representation

   contains

     ! Deferred procedures
     procedure(representation_num_members_interface), deferred :: num_members
     procedure(representation_member_interface)     , deferred :: member
     procedure(representation_local_index_interface), deferred :: local_index

     ! Provided procedures
     procedure :: members => representation_members
     procedure :: has     => representation_has

  end type set_representation

  abstract interface

     !----------------------------------------------------------------!
     ! Interface for the number of members
     !----------------------------------------------------------------!

     pure integer function representation_num_members_interface(this)
       import set_representation
       class(set_representation), intent(in) :: this
     end function representation_num_members_interface

     !----------------------------------------------------------------!
     ! Interface for the member at a position
     !----------------------------------------------------------------!

     pure integer function representation_member_interface(this, position)
       import set_representation
       class(set_representation), intent(in) :: this
       integer                  , intent(in) :: position
     end function representation_member_interface

     !----------------------------------------------------------------!
     ! Interface for the position of a value, 0 for a non-member
     !----------------------------------------------------------------!

     pure integer function representation_local_index_interface(this, value)
       import set_representation
       class(set_representation), intent(in) :: this
       integer                  , intent(in) :: value
     end function representation_local_index_interface

  end interface

  !-------------------------------------------------------------------!
  ! Type that models the counted representation: members 1..n, stored
  ! as n. O(1), whatever n is: every primitive is a comparison or an
  ! identity.
  !-------------------------------------------------------------------!

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

  !-------------------------------------------------------------------!
  ! Type that models the listed representation: an explicit array of
  ! member values, in declaration order. The representation describes
  ! a set with or without any declared ambient; the embedding is not
  ! stored here. A sorted permutation indexes that same array:
  ! construction takes O(n log n), membership O(log n), and the
  ! enumeration retains the first occurrence of every value. Storage is
  ! two integers per member:
  !
  !   values given   [ 7, 3, 7, 5 ]
  !   listed         [ 7, 3, 5 ]      first occurrences, declaration order
  !   inverse        [ 2, 3, 1 ]      positions in listed, by increasing value:
  !                                   listed(2) = 3 < listed(3) = 5 < listed(1) = 7
  !   local_index(5)                  a binary search over inverse: 3
  !-------------------------------------------------------------------!

  type, extends(set_representation) :: listed_set_representation

     integer, allocatable, private :: listed(:)
     integer, allocatable, private :: inverse(:)

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
  ! Returns whether a value is a member: a nonzero position, for every
  ! concretion
  !===================================================================!

  pure logical function representation_has(this, value)

    class(set_representation), intent(in) :: this
    integer                  , intent(in) :: value

    representation_has = this % local_index(value) /= 0

  end function representation_has

  !===================================================================!
  ! Returns the members: the enumeration, for every concretion
  !===================================================================!

  pure subroutine representation_members(this, values)

    class(set_representation), intent(in)  :: this
    integer, allocatable     , intent(out) :: values(:)

    integer :: k

    values = [(this % member(k), k = 1, this % num_members())]

  end subroutine representation_members

  !===================================================================!
  ! Creates the counted representation of n members, n clamped at
  ! zero: the empty set is a set
  !===================================================================!

  pure type(counted_set_representation) function create_counted(n) &
       & result(this)

    integer, intent(in) :: n

    this % n = max(0, n)

  end function create_counted

  !===================================================================!
  ! Returns n, the number of counted members
  !===================================================================!

  pure integer function counted_num_members(this)

    class(counted_set_representation), intent(in) :: this

    counted_num_members = this % n

  end function counted_num_members

  !===================================================================!
  ! Returns the counted member at a position, the position itself
  !===================================================================!

  pure integer function counted_member(this, position)

    class(counted_set_representation), intent(in) :: this
    integer                          , intent(in) :: position

    counted_member = position

  end function counted_member

  !===================================================================!
  ! Returns the position of a value among 1..n, or 0
  !===================================================================!

  pure integer function counted_local_index(this, value)

    class(counted_set_representation), intent(in) :: this
    integer                          , intent(in) :: value

    counted_local_index = 0
    if (value >= 1 .and. value <= this % n) counted_local_index = value

  end function counted_local_index

  !===================================================================!
  ! Creates the listed representation of the values passed in, each
  ! retained once, the first occurrence keeping its position: the
  ! enumeration law needs injectivity, and the declaration order is the
  ! representation's own
  !===================================================================!

  pure type(listed_set_representation) function create_listed(values) &
       & result(this)

    integer, intent(in) :: values(:)

    integer, allocatable :: order(:), work(:)
    integer              :: j, num_retained, last_value, at

    ! The positions sorted by value, stably
    allocate(order(size(values)), work(size(values)))
    order = [(j, j = 1, size(values))]
    call order_by_member(values, order, work, 1, size(values))

    ! The first position of each distinct value
    num_retained = 0
    do j = 1, size(values)
       if (num_retained > 0) then
          if (values(order(j)) == last_value) cycle
       end if
       last_value = values(order(j))
       num_retained = num_retained + 1
       order(num_retained) = order(j)
    end do

    ! The sort is stable, so each retained index is the member's first
    ! occurrence. Number those indices in declaration order, then read
    ! the resulting positions in sorted member order for inverse lookup.
    work = 0
    work(order(1:num_retained)) = 1
    at = 0
    allocate(this % listed(num_retained))
    do j = 1, size(values)
       if (work(j) == 0) cycle
       at = at + 1
       work(j) = at
       this % listed(at) = values(j)
    end do
    this % inverse = work(order(1:num_retained))

  end function create_listed

  !===================================================================!
  ! Sorts the positions order(first:last) by their values, a stable
  ! merge sort: equal values keep their order
  !===================================================================!

  recursive pure subroutine order_by_member(values, order, work, first, last)

    integer, intent(in)    :: values(:), first, last
    integer, intent(inout) :: order(:), work(:)

    integer :: middle, left, right, k

    if (last <= first) return
    middle = first + (last - first) / 2
    call order_by_member(values, order, work, first, middle)
    call order_by_member(values, order, work, middle + 1, last)
    left = first
    right = middle + 1
    do k = first, last
       if (left > middle) then
          work(k) = order(right)
          right = right + 1
       else if (right > last) then
          work(k) = order(left)
          left = left + 1
       else if (values(order(left)) <= values(order(right))) then
          work(k) = order(left)
          left = left + 1
       else
          work(k) = order(right)
          right = right + 1
       end if
    end do
    order(first:last) = work(first:last)

  end subroutine order_by_member

  !===================================================================!
  ! Returns the number of listed members
  !===================================================================!

  pure integer function listed_num_members(this)

    class(listed_set_representation), intent(in) :: this

    if (allocated(this % listed)) then
       listed_num_members = size(this % listed)
    else
       listed_num_members = 0
    end if

  end function listed_num_members

  !===================================================================!
  ! Returns the listed member at a position
  !===================================================================!

  pure integer function listed_member(this, position)

    class(listed_set_representation), intent(in) :: this
    integer                         , intent(in) :: position

    listed_member = this % listed(position)

  end function listed_member

  !===================================================================!
  ! Returns the position of a value, by a binary search over the
  ! listed members in increasing order, or 0 for a non-member
  !===================================================================!

  pure integer function listed_local_index(this, value)

    class(listed_set_representation), intent(in) :: this
    integer                         , intent(in) :: value

    integer :: first, last, middle, at

    listed_local_index = 0
    if (.not. allocated(this % listed)) return

    first = 1
    last = size(this % listed)
    do while (first <= last)
       middle = first + (last - first) / 2
       at = this % inverse(middle)
       if (this % listed(at) < value) then
          first = middle + 1
       else if (this % listed(at) > value) then
          last = middle - 1
       else
          listed_local_index = at
          return
       end if
    end do

  end function listed_local_index

end module map_set_representation
