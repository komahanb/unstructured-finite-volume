!=====================================================================!
! THE IDENTITY-ROW TABLE beneath every token-keyed map.
!
! A map keyed on identity is a set of rows, each a copied identity
! token beside a payload the map owns. The token column and its
! mechanics are the same in every such map, and they are defined here
! once:
!
!      position      the row of a key, or zero when no row stores it;
!                    unallocated storage returns zero
!      row           the row of a key; no row stops the program with
!                    the map's own message
!      append        one more row, its key copied by value, the new
!                    row's position returned; an undeclared key or a
!                    key already present stops the program with the
!                    map's own message
!      key           the token stored at a row
!      removal       one row removed, the rest retaining their order
!
! A key is a declared token: a row keyed on an undeclared token could
! never be found again, so append rejects one. This module owns the
! token column ALONE. It has no dependency on set representations,
! labels, inclusions, fields, or graph mathematics: a map indexes its
! own payload by the row positions this table returns, and reads that
! payload itself. So the different laws the maps enforce - one ambient
! per part, a status that must be known before it is read, a
! representation copied whole - stay in the maps, and only the
! copied-token bookkeeping is shared.
!
!             THE LIFETIME LAW
!
! Keys are type(token), copied by value at append. The table references
! no object it was built from, so a map may outlive every variable
! that populated it, and append requires no TARGET. Nothing here stores
! a pointer or a graph.
!
!             WHAT IS REJECTED
!
! An append of an undeclared key or of a key already present, a row
! query for a key with no row, and a removal of a position outside the
! rows, each stop the program. The map passes the message that states
! its own law (named once, attached once, one ambient), so the
! rejection is reported in the map's terms and checked here once.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module map_token_rows

  use token_identity, only : token, index_of

  implicit none

  private
  public :: identity_rows

  type :: identity_rows

     type(token), allocatable, private :: keys(:)

   contains

     procedure :: position
     procedure :: row
     procedure :: key
     procedure :: num_rows
     procedure :: append
     procedure :: remove

  end type identity_rows

contains

  !===================================================================!
  ! The row of a key, or zero. Unallocated storage returns zero, so a
  ! map need not check for an empty table itself.
  !===================================================================!

  pure integer function position(this, key) result(at)

    class(identity_rows), intent(in) :: this
    type(token)      , intent(in) :: key

    at = 0
    if (.not. allocated(this % keys)) return
    at = index_of(this % keys, key)

  end function position

  !===================================================================!
  ! The row of a key. No row stops the program with the caller's
  ! message, so a map's read is one call.
  !===================================================================!

  integer function row(this, key, message) result(at)

    class(identity_rows), intent(in) :: this
    type(token)         , intent(in) :: key
    character(len=*)    , intent(in) :: message

    at = this % position(key)
    if (at == 0) error stop message

  end function row

  pure type(token) function key(this, at)

    class(identity_rows), intent(in) :: this
    integer             , intent(in) :: at

    key = this % keys(at)

  end function key

  pure integer function num_rows(this)

    class(identity_rows), intent(in) :: this

    num_rows = 0
    if (allocated(this % keys)) num_rows = size(this % keys)

  end function num_rows

  !===================================================================!
  ! Append one row and return its position. The key is copied by
  ! value. An undeclared key stops the program with the first
  ! message, a key already present with the second.
  !===================================================================!

  function append(this, key, undeclared, duplicate) result(at)

    class(identity_rows), intent(inout) :: this
    type(token)         , intent(in)    :: key
    character(len=*)    , intent(in)    :: undeclared, duplicate
    integer :: at

    if (.not. key % declared())     error stop undeclared
    if (this % position(key) /= 0) error stop duplicate

    if (.not. allocated(this % keys)) allocate(this % keys(0))
    this % keys = [this % keys, key]
    at = size(this % keys)

  end function append

  !===================================================================!
  ! Remove the row at a position; the rest retain their order, because a
  ! map's parallel payload is compacted the same way. A position
  ! outside the rows stops the program.
  !===================================================================!

  subroutine remove(this, at)

    class(identity_rows), intent(inout) :: this
    integer          , intent(in)    :: at

    if (at < 1 .or. at > this % num_rows()) then
       error stop 'map_token_rows: a removal requires an existing row'
    end if

    this % keys = [this % keys(1:at - 1), this % keys(at + 1:)]

  end subroutine remove

end module map_token_rows
