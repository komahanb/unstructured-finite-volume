!=====================================================================!
! THE SPECIMEN'S CARRIERS - earned at LEVEL 0, and the sole home of
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
! Its framework dependency is graph_carrier and nothing else. Level 0
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
!   1. It is keyed by the carrier's DECLARED NAME, which graph_carrier
!      already states is metadata carried for the reader and no part
!      of the mathematics. No law reads it.
!
!   2. It is keyed by MEMBER VALUE, never by position. A caller must
!      already hold the member - which it can only have got by asking
!      the carrier - before it can ask what to call it.
!
!   3. A carrier it does not recognise gets its member printed as the
!      integer it is. Nothing is invented, and the renderer above can
!      be pointed at a carrier this file has never heard of.
!
! So a picture's ORDER can never come from here: order is whatever
! member(1), member(2), ... say it is, and this file is asked only
! afterwards, one member at a time.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module visualization_carriers_fixture

  use fractal_graph        , only : set_graph => graph
  use graph_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use graph_set_map        , only : set_map
  use graph_label_map      , only : label_map

  implicit none

  private
  public :: structural_carriers, label_for

contains

  !===================================================================!
  ! The seven declared domains, each signing its own identity once.
  ! Four state carriers and three occurrence carriers: |X1| = |X2|,
  ! and X1 is still not X2.
  !===================================================================!

  subroutine structural_carriers(x0, x1, x2, x3, e1, e2, e3, sets, labels)

    type(set_graph), intent(out)   :: x0, x1, x2, x3
    type(set_graph), intent(out)   :: e1, e2, e3
    type(set_map)  , intent(inout) :: sets
    type(label_map), intent(inout) :: labels

    call x0 % declare()
    call sets   % bind(x0, counted_set_representation(4))
    call labels % bind(x0, 'X0')
    call x1 % declare()
    call sets   % bind(x1, counted_set_representation(3))
    call labels % bind(x1, 'X1')
    call x2 % declare()
    call sets   % bind(x2, counted_set_representation(3))
    call labels % bind(x2, 'X2')
    call x3 % declare()
    call sets   % bind(x3, counted_set_representation(2))
    call labels % bind(x3, 'X3')

    call e1 % declare()
    call sets   % bind(e1, counted_set_representation(5))
    call labels % bind(e1, 'E1')
    call e2 % declare()
    call sets   % bind(e2, counted_set_representation(4))
    call labels % bind(e2, 'E2')
    call e3 % declare()
    call sets   % bind(e3, counted_set_representation(3))
    call labels % bind(e3, 'E3')

  end subroutine structural_carriers

  !===================================================================!
  ! What the reader calls one member of one carrier. Display only.
  !===================================================================!

  function label_for(carrier, member, labels) result(text)

    type(set_graph), intent(in) :: carrier
    integer        , intent(in) :: member
    type(label_map), intent(in) :: labels

    character(len=:), allocatable :: text
    character(len=12)             :: buf

    select case (labels % label_of(carrier))
    case ('X0'); text = pick(['a', 'b', 'c', 'd'], member)
    case ('X1'); text = pick(['p', 'q', 'r'], member)
    case ('X2'); text = pick(['u', 'v', 'w'], member)
    case ('X3'); text = pick(['m', 'n'], member)
    case ('E1'); text = pick(['e11', 'e12', 'e13', 'e14', 'e15'], member)
    case ('E2'); text = pick(['e21', 'e22', 'e23', 'e24'], member)
    case ('E3'); text = pick(['e31', 'e32', 'e33'], member)
    case default; text = ''
    end select

    ! An unrecognised carrier - or a member outside the alphabet -
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

end module visualization_carriers_fixture
