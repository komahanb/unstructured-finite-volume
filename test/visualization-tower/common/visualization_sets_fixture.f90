!=====================================================================!
! THE SPECIMEN'S SETS - earned at LEVEL 0, and the sole home of
! the seven declared domains.
!
!      X0 = { a  b  c  d }
!      X1 = { p  q  r }
!      X2 = { u  v  w }
!      X3 = { m  n }
!
!      E1 = { e11 e12 e13 e14 e15 }
!      E2 = { e21 e22 e23 e24 }
!      E3 = { e31 e32 e33 }
!
! Its framework dependency is graph_set and nothing else. Level 0
! knows about members; it does not know that a relation exists, or an
! operator, or a graph, or a picture.
!
! WHY THE SIZES ARE WRITTEN OUT HERE. They are literals in this file
! and named constants in visualization_assert, so a level that checks
! |X0| = NX0 is comparing two independent statements rather than
! reading one number back to itself.
!
!                    LABELS ARE METADATA, NOT MATHEMATICS
!
! label_for answers the reader's name for a member - 'a' for X0's
! first, 'e23' for E2's third - and every level above uses it for
! display only. Three properties keep it honest:
!
!   1. It is keyed by the set's DECLARED NAME, which graph_set
!      already states is metadata carried for the reader and no part
!      of the mathematics. No law reads it.
!
!   2. It is keyed by MEMBER VALUE, never by position. A caller must
!      already hold the member - which it can only have got by asking
!      the set - before it can ask what to call it.
!
!   3. A set it does not recognise gets its member printed as the
!      integer it is. Nothing is invented, and the renderer above can
!      be pointed at a set this file has never heard of.
!
! So a picture's ORDER can never come from here: order is whatever
! member(1), member(2), ... say it is, and this file is asked only
! afterwards, one member at a time.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module visualization_sets_fixture

  use graph_set, only : index_set, set

  implicit none

  private
  public :: structural_sets, label_for

contains

  !===================================================================!
  ! The seven declared domains, each signing its own identity once.
  ! Four state sets and three occurrence sets: |X1| = |X2|,
  ! and X1 is still not X2.
  !===================================================================!

  subroutine structural_sets(x0, x1, x2, x3, e1, e2, e3)

    type(index_set), intent(out) :: x0, x1, x2, x3
    type(index_set), intent(out) :: e1, e2, e3

    x0 = index_set('X0', 4)
    x1 = index_set('X1', 3)
    x2 = index_set('X2', 3)
    x3 = index_set('X3', 2)

    e1 = index_set('E1', 5)
    e2 = index_set('E2', 4)
    e3 = index_set('E3', 3)

  end subroutine structural_sets

  !===================================================================!
  ! What the reader calls one member of one set. Display only.
  !===================================================================!

  function label_for(set, member) result(text)

    class(set), intent(in) :: set
    integer          , intent(in) :: member

    character(len=:), allocatable :: text
    character(len=12)             :: buf

    select case (set % name())
    case ('X0'); text = pick(['a', 'b', 'c', 'd'], member)
    case ('X1'); text = pick(['p', 'q', 'r'], member)
    case ('X2'); text = pick(['u', 'v', 'w'], member)
    case ('X3'); text = pick(['m', 'n'], member)
    case ('E1'); text = pick(['e11', 'e12', 'e13', 'e14', 'e15'], member)
    case ('E2'); text = pick(['e21', 'e22', 'e23', 'e24'], member)
    case ('E3'); text = pick(['e31', 'e32', 'e33'], member)
    case default; text = ''
    end select

    ! An unrecognised set - or a member outside the alphabet -
    ! gets the integer it actually is. Nothing is invented.
    if (len(text) .eq. 0) then
       write(buf, '(i0)') member
       text = trim(buf)
    end if

  end function label_for

  !===================================================================!
  ! Read the k-th entry of an alphabet, or answer the empty name.
  !===================================================================!

  function pick(alphabet, k) result(text)

    character(len=*), intent(in) :: alphabet(:)
    integer         , intent(in) :: k

    character(len=:), allocatable :: text

    if (k .ge. 1 .and. k .le. size(alphabet)) then
       text = trim(alphabet(k))
    else
       text = ''
    end if

  end function pick

end module visualization_sets_fixture
