!=====================================================================!
! Module that states the identity-row table beneath every token-keyed
! map.
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
!      removal       one row removed, the rest keeping their order
!
!   this table               a map's payload, parallel to it
!   row 1   token  ------>   payload(1)
!   row 2   token  ------>   payload(2)
!   ...
!
! A key is a declared token: a row keyed on an undeclared token could
! never be found again, so append rejects one. This module owns the
! token column alone. It does not depend on set representations,
! labels, inclusions, fields or graph mathematics: a map indexes its
! own payload by the row positions this table returns, and reads that
! payload itself. So the different laws the maps enforce, one ambient
! per part, a status that must be known before it is read, a
! representation copied whole, stay in the maps, and only the
! copied-token bookkeeping is shared.
!
! The lifetime law. Keys are type(token), copied by value at append.
! The table references no object it was built from, so a map may
! outlive every variable that populated it, and append requires no
! TARGET. Nothing here stores a pointer or a graph.
!
! What is rejected. An append of an undeclared key or of a key already
! present, a row query for a key with no row, and a removal of a
! position outside the rows each stop the program. The map passes the
! message that states its own law (named once, attached once, one
! ambient), so the rejection is reported in the map's terms and
! checked here once.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module map_token_rows

  use iso_fortran_env, only : int64
  use token_identity , only : token

  implicit none

  private

  public :: identity_rows

  !-------------------------------------------------------------------!
  ! Type that models the table: a column of keys and an open-addressing
  ! table of slots on their serial numbers. For example,
  !
  !   keys                  slot      1   2   3   4   5   6   7   8
  !   row 1   (1, 7)        slot_row  0   2   0   1   0   3   0   0
  !   row 2   (1, 3)
  !   row 3   (1, 12)       a slot stores the row of a key, or 0
  !-------------------------------------------------------------------!

  type :: identity_rows

     ! The rows 1..num_stored of a column with spare capacity, so that
     ! an append costs amortised constant time rather than a whole copy
     type(token), allocatable, private :: keys(:)
     integer                 , private :: num_stored = 0

     ! The row at each slot of an open-addressing table on the serial
     ! number, zero for an empty slot. The number of slots is a power of
     ! two of at least twice the rows, so a position costs constant
     ! time.
     integer, allocatable, private :: slot_row(:)

   contains

     ! Queries
     procedure :: position
     procedure :: row
     procedure :: key
     procedure :: num_rows

     ! Writers
     procedure :: append
     procedure :: remove

  end type identity_rows

contains

  !===================================================================!
  ! Returns the row of a key, or zero. Unallocated storage returns
  ! zero, so a map need not check for an empty table itself. The slots
  ! are read in turn from the key's own slot, wrapping at the end,
  ! until its row or an empty slot:
  !
  !   slot_of(key)
  !      |
  !      v
  !   [ r ][ r ][ 0 ]        a row whose key matches: found
  !                          an empty slot: not present, zero
  !===================================================================!

  pure integer function position(this, key) result(at)

    class(identity_rows), intent(in) :: this
    type(token)         , intent(in) :: key

    integer :: slot

    at = 0
    if (this % num_stored == 0 .or. .not. key % declared()) return
    slot = slot_of(key, size(this % slot_row))
    do
       at = this % slot_row(slot)
       if (at == 0) return
       if (this % keys(at) % matches(key)) return
       slot = mod(slot, size(this % slot_row)) + 1
    end do

  end function position

  !===================================================================!
  ! Returns the slot a key starts from: the high bits of the product of
  ! its serial with an odd constant, masked to the power-of-two number
  ! of slots,
  !
  !   serial * 2654435761  -->  shifted right by 16  -->  its low bits,
  !                                                       num_slots - 1
  !                                                       as the mask, + 1
  !===================================================================!

  pure integer function slot_of(key, num_slots)

    type(token), intent(in) :: key
    integer    , intent(in) :: num_slots

    slot_of = int(iand(ishft(int(key % serial_number(), int64) * 2654435761_int64, -16), &
         & int(num_slots - 1, int64))) + 1

  end function slot_of

  !===================================================================!
  ! Places one row into the first empty slot from its key's slot
  !===================================================================!

  subroutine place(this, at)

    class(identity_rows), intent(inout) :: this
    integer             , intent(in)    :: at

    integer :: slot

    slot = slot_of(this % keys(at), size(this % slot_row))
    do while (this % slot_row(slot) /= 0)
       slot = mod(slot, size(this % slot_row)) + 1
    end do
    this % slot_row(slot) = at

  end subroutine place

  !===================================================================!
  ! Reserves capacity for the required rows: the column doubles when
  ! full, and the slots are rebuilt at twice the rows, a power of two
  ! from 16, when they fall below that
  !===================================================================!

  subroutine reserve(this, required)

    class(identity_rows), intent(inout) :: this
    integer             , intent(in)    :: required

    type(token), allocatable :: larger(:)
    integer                  :: num_slots, at

    ! The column
    if (.not. allocated(this % keys)) allocate(this % keys(max(required, 8)))
    if (required > size(this % keys)) then
       allocate(larger(max(required, 2 * size(this % keys))))
       larger(1:this % num_stored) = this % keys(1:this % num_stored)
       call move_alloc(larger, this % keys)
    end if

    ! The slots, rebuilt only when too few
    num_slots = 16
    if (allocated(this % slot_row)) num_slots = size(this % slot_row)
    if (allocated(this % slot_row) .and. 2 * required <= num_slots) return
    do while (num_slots < 2 * required)
       num_slots = 2 * num_slots
    end do
    if (allocated(this % slot_row)) deallocate(this % slot_row)
    allocate(this % slot_row(num_slots), source=0)
    do at = 1, this % num_stored
       call place(this, at)
    end do

  end subroutine reserve

  !===================================================================!
  ! Returns the row of a key. No row stops the program with the
  ! caller's message, so that a map's read is one call.
  !===================================================================!

  integer function row(this, key, message) result(at)

    class(identity_rows), intent(in) :: this
    type(token)         , intent(in) :: key
    character(len=*)    , intent(in) :: message

    at = this % position(key)
    if (at == 0) error stop message

  end function row

  !===================================================================!
  ! Returns the token stored at a row
  !===================================================================!

  pure type(token) function key(this, at)

    class(identity_rows), intent(in) :: this
    integer             , intent(in) :: at

    key = this % keys(at)

  end function key

  !===================================================================!
  ! Returns the number of rows
  !===================================================================!

  pure integer function num_rows(this)

    class(identity_rows), intent(in) :: this

    num_rows = this % num_stored

  end function num_rows

  !===================================================================!
  ! Appends one row and returns its position. The key is copied by
  ! value. An undeclared key stops the program with the first message,
  ! a key already present with the second.
  !===================================================================!

  function append(this, key, undeclared, duplicate) result(at)

    class(identity_rows), intent(inout) :: this
    type(token)         , intent(in)    :: key
    character(len=*)    , intent(in)    :: undeclared, duplicate
    integer                             :: at

    if (.not. key % declared())     error stop undeclared
    if (this % position(key) /= 0) error stop duplicate

    call reserve(this, this % num_stored + 1)
    this % num_stored = this % num_stored + 1
    at = this % num_stored
    this % keys(at) = key
    call place(this, at)

  end function append

  !===================================================================!
  ! Removes the row at a position. The rest keep their order, because
  ! a map's parallel payload is compacted the same way, and the slots
  ! are placed again:
  !
  !   before    [ k1 ][ k2 ][ k3 ][ k4 ]        remove row 2
  !   after     [ k1 ][ k3 ][ k4 ]              k3 is row 2, k4 row 3
  !
  ! A position outside the rows stops the program.
  !===================================================================!

  subroutine remove(this, at)

    class(identity_rows), intent(inout) :: this
    integer             , intent(in)    :: at

    integer            :: slot
    character(len=100) :: message

    if (at < 1 .or. at > this % num_rows()) then
       write(message,'(a,i0,a,i0)') 'map_token_rows: remove requires an existing row; at = ', &
            & at, ', num_rows = ', this % num_rows()
       error stop trim(message)
    end if

    this % keys(at:this % num_stored - 1) = this % keys(at + 1:this % num_stored)
    this % num_stored = this % num_stored - 1
    this % slot_row = 0
    do slot = 1, this % num_stored
       call place(this, slot)
    end do

  end subroutine remove

end module map_token_rows
