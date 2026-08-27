!=====================================================================!
! THE IDENTITY-ROW TABLE beneath every token-keyed map.
!
! A map keyed on identity is a set of rows, each a copied identity
! token beside a payload the map owns. The token column and its
! mechanics are the same in every such map, and they live here once:
!
!      position      the row of a key, or zero when no row holds it;
!                    unallocated storage answers zero
!      append        one more row, its key copied by value, the new
!                    row's position returned
!      removal       one row dropped, the rest keeping their order
!
! A key is a declared token: a map checks token % declared() before a
! write, since a row keyed on an undeclared token could never be found
! again, and states its own refusal. This module owns the token column ALONE. It knows nothing of set
! representations, labels, inclusions, fields, or graph mathematics: a
! map hangs its own payload off the row positions this table hands
! back, and reads that payload itself. So the different laws the maps
! enforce - one ambient per part, a status that must be known before
! it is read, a representation copied whole - stay in the maps, and
! only the copied-token bookkeeping is shared.
!
!             THE LIFETIME LAW
!
! Keys are type(token), copied by value at append. The table borrows
! no object it was built from, so a map may outlive every variable
! that populated it, and append demands no TARGET. Nothing here stores
! a pointer or a graph.
!
!             WHAT IS REFUSED
!
! An append of a key already present, and a removal of a position
! outside the rows, each stop the program: a caller that has not first
! asked position would otherwise duplicate or misplace a row. The
! message is generic here; a map states its own law before it calls,
! so the map-specific refusal (named once, attached once, one ambient)
! is reported there and this is only the backstop.
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
     procedure :: num_rows
     procedure :: append
     procedure :: drop

  end type identity_rows

contains

  !===================================================================!
  ! The row of a key, or zero. Unallocated storage answers zero, so a
  ! map need not guard an empty table itself.
  !===================================================================!

  pure integer function position(this, key) result(at)

    class(identity_rows), intent(in) :: this
    type(token)      , intent(in) :: key

    at = 0
    if (.not. allocated(this % keys)) return
    at = index_of(this % keys, key)

  end function position

  pure integer function num_rows(this)

    class(identity_rows), intent(in) :: this

    num_rows = 0
    if (allocated(this % keys)) num_rows = size(this % keys)

  end function num_rows

  !===================================================================!
  ! Append one row and return its position. The key is copied by
  ! value. A key already present stops the program; the caller refuses
  ! it first with its own message, so this is the backstop.
  !===================================================================!

  function append(this, key) result(at)

    class(identity_rows), intent(inout) :: this
    type(token)      , intent(in)    :: key
    integer :: at

    type(token), allocatable :: grown(:)
    integer :: n

    if (this % position(key) /= 0) then
       error stop 'map_token_rows: a row is appended for a key not already present'
    end if

    if (.not. allocated(this % keys)) allocate(this % keys(0))

    n = size(this % keys)
    allocate(grown(n + 1))
    grown(1:n)   = this % keys
    grown(n + 1) = key
    call move_alloc(grown, this % keys)

    at = n + 1

  end function append

  !===================================================================!
  ! Drop the row at a position; the rest keep their order, because a
  ! map's parallel payload is compacted the same way. A position
  ! outside the rows stops the program.
  !===================================================================!

  subroutine drop(this, at)

    class(identity_rows), intent(inout) :: this
    integer          , intent(in)    :: at

    type(token), allocatable :: kept(:)
    integer :: n, k, m

    n = this % num_rows()
    if (at < 1 .or. at > n) then
       error stop 'map_token_rows: a removal touches an existing row'
    end if

    allocate(kept(n - 1))
    m = 0
    do k = 1, n
       if (k == at) cycle
       m = m + 1
       kept(m) = this % keys(k)
    end do
    call move_alloc(kept, this % keys)

  end subroutine drop

end module map_token_rows
