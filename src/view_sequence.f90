!=====================================================================!
! SEQUENCE VIEW
!
! One view of a graph. The kernel supplies possibility - graph, branch,
! graph - and never the word sequence. A finite sequence is
! REPRESENTABLE over that possibility; it is not encoded by it.
!
! THE REPRESENTATION LAW. A sequence is represented by a BRANCH, not by
! a graph:
!
!     NULL            the empty sequence
!     KNOWN -> cell   a nonempty sequence, beginning at that cell
!     UNKNOWN         the sequence is not known
!
! and a cell is a graph:
!
!     branch(1) = KNOWN -> element      the element, always KNOWN
!     branch(2)                         the rest, again a sequence branch
!
! So [a,b,c] is (a, (b, (c, NULL))), and the empty sequence has no cell
! and needs no graph. A holder writes
!
!     holder % branch(i) = null_branch()
!
! for the empty sequence, which is why these procedures take a branch.
!
! TWO KINDS OF FAILURE, KEPT APART.
!
!     malformed     a cell whose branch(1) is NULL or UNKNOWN. The
!                   representation is wrong. Refused, always.
!     unknown       a spine that reaches UNKNOWN. The representation is
!                   right and the answer is not yet determined. Refused
!                   only when the answer depends on the unknown part -
!                   sequence_element(b, 1) succeeds on a sequence whose
!                   tail is UNKNOWN.
!
! NULL is not UNKNOWN here either: NULL ends a sequence, UNKNOWN
! withholds it.
!
! COMPLEXITY. size is O(n) and element(k) is O(k), and this is the
! semantic view, not the hot representation. Where repeated indexed
! access matters, compile to a contiguous representation.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module view_sequence

  use graph_fractal, only : graph, graph_branch, &
       & GRAPH_NULL, GRAPH_UNKNOWN, GRAPH_KNOWN

  implicit none

  private
  public :: sequence_defined, sequence_size
  public :: sequence_element, sequence_contains

contains

  !===================================================================!
  ! Is the whole extent known: does the spine reach NULL. False when
  ! it reaches UNKNOWN. A malformed cell is refused rather than
  ! reported, because it is not a state of the sequence.
  !===================================================================!

  logical function sequence_defined(b) result(known_extent)

    type(graph_branch), intent(in) :: b

    type(graph), pointer :: cell

    known_extent = .false.
    if (b % status() .eq. GRAPH_UNKNOWN) return
    if (b % status() .eq. GRAPH_NULL) then
       known_extent = .true.
       return
    end if

    cell => b % known()
    do
       call require_cell(cell)
       if (cell % branch(2) % status() .eq. GRAPH_NULL) exit
       if (cell % branch(2) % status() .eq. GRAPH_UNKNOWN) return
       cell => cell % branch(2) % known()
    end do
    known_extent = .true.

  end function sequence_defined

  !===================================================================!
  ! The number of elements. The answer depends on the whole spine, so
  ! an unknown extent is refused.
  !===================================================================!

  integer function sequence_size(b) result(n)

    type(graph_branch), intent(in) :: b

    type(graph), pointer :: cell

    if (.not. sequence_defined(b)) then
       error stop 'view_sequence: the extent depends on an unknown tail'
    end if

    n = 0
    if (b % status() .eq. GRAPH_NULL) return

    cell => b % known()
    do
       n = n + 1
       if (cell % branch(2) % status() .eq. GRAPH_NULL) exit
       cell => cell % branch(2) % known()
    end do

  end function sequence_size

  !===================================================================!
  ! The k-th element, counting from one. Only the first k cells are
  ! traversed, so an unknown tail beyond k is no obstacle.
  !===================================================================!

  function sequence_element(b, k) result(element)

    type(graph_branch), intent(in) :: b
    integer           , intent(in) :: k
    type(graph), pointer           :: element

    type(graph), pointer :: cell
    integer              :: i

    if (k .lt. 1) then
       error stop 'view_sequence: a sequence is indexed from one'
    end if
    call require_reachable(b)

    cell => b % known()
    do i = 1, k - 1
       call require_cell(cell)
       call require_reachable(cell % branch(2))
       cell => cell % branch(2) % known()
    end do

    call require_cell(cell)
    element => cell % branch(1) % known()

  end function sequence_element

  !===================================================================!
  ! Does the sequence hold this graph, by identity. Answered as soon
  ! as it is found; refused if the spine runs into UNKNOWN first,
  ! because then the answer depends on the unknown part.
  !
  ! This is one traversal. A caller looping sequence_element instead
  ! would be O(n^2), which is why membership is here and not derived
  ! outside.
  !===================================================================!

  logical function sequence_contains(b, g) result(found)

    type(graph_branch), intent(in) :: b
    type(graph)       , intent(in) :: g

    type(graph), pointer :: cell, element

    found = .false.
    if (b % status() .eq. GRAPH_NULL) return
    if (b % status() .eq. GRAPH_UNKNOWN) then
       error stop 'view_sequence: membership depends on an unknown sequence'
    end if

    cell => b % known()
    do
       call require_cell(cell)
       element => cell % branch(1) % known()
       if (element % same_as(g)) then
          found = .true.
          return
       end if
       if (cell % branch(2) % status() .eq. GRAPH_NULL) return
       if (cell % branch(2) % status() .eq. GRAPH_UNKNOWN) then
          error stop 'view_sequence: membership depends on an unknown tail'
       end if
       cell => cell % branch(2) % known()
    end do

  end function sequence_contains

  !===================================================================!
  ! The two guards. require_cell refuses a malformed representation;
  ! require_reachable refuses a step that cannot be taken.
  !===================================================================!

  subroutine require_cell(cell)

    type(graph), intent(in) :: cell

    if (cell % branch(1) % status() .ne. GRAPH_KNOWN) then
       error stop 'view_sequence: a sequence cell holds a KNOWN element'
    end if

  end subroutine require_cell

  subroutine require_reachable(b)

    type(graph_branch), intent(in) :: b

    if (b % status() .eq. GRAPH_NULL) then
       error stop 'view_sequence: the sequence has no such element'
    end if
    if (b % status() .eq. GRAPH_UNKNOWN) then
       error stop 'view_sequence: that element lies beyond an unknown tail'
    end if

  end subroutine require_reachable

end module view_sequence
