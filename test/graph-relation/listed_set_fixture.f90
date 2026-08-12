!=====================================================================!
! A second carrier concretion, living in the tests on purpose: the
! library keeps one citizen until real code needs another, but the
! relation's claim to be generic over member_set must be PROVED, not
! asserted. This listed set stores its members outright - sparse,
! unordered-world indices like sensor numbers - and stands in a
! signature slot beside a counted set.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module listed_set_fixture

  use graph_carrier, only : member_set

  implicit none

  private
  public :: listed_set

  type, extends(member_set) :: listed_set

     integer, allocatable, private :: roll(:)

   contains

     procedure :: size    => listed_size
     procedure :: member  => listed_member
     procedure :: members => listed_members
     procedure :: has     => listed_has

  end type listed_set

  interface listed_set
     module procedure create_listed
  end interface listed_set

contains

  type(listed_set) function create_listed(name, members) result(this)

    character(len=*), intent(in) :: name
    integer         , intent(in) :: members(:)

    this % roll = members
    call this % declare(name)

  end function create_listed

  pure integer function listed_size(this)

    class(listed_set), intent(in) :: this

    listed_size = size(this % roll)

  end function listed_size

  pure integer function listed_member(this, local_index)

    class(listed_set), intent(in) :: this
    integer          , intent(in) :: local_index

    listed_member = this % roll(local_index)

  end function listed_member

  pure subroutine listed_members(this, indices)

    class(listed_set)   , intent(in)  :: this
    integer, allocatable, intent(out) :: indices(:)

    indices = this % roll

  end subroutine listed_members

  pure logical function listed_has(this, member)

    class(listed_set), intent(in) :: this
    integer          , intent(in) :: member

    listed_has = any(this % roll == member)

  end function listed_has

end module listed_set_fixture
