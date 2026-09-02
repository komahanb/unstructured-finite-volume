!=====================================================================!
! SEQUENCE VIEW
!
! One view of a graph. The kernel supplies the structure - graph, branch,
! graph - and never the word sequence. A finite sequence is
! REPRESENTABLE over that structure; it is not encoded by it.
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
! and needs no graph. A containing graph writes
!
!     holder % branch(i) = null_branch()
!
! for the empty sequence, which is why these procedures take a branch.
!
! TWO KINDS OF FAILURE, KEPT APART.
!
!     malformed     a cell whose branch(1) is NULL or UNKNOWN. The
!                   representation is wrong. Refused, always.
!     unknown       a chain of cells that reaches UNKNOWN. The
!                   representation is correct and the result is not yet
!                   determined. Refused only when the result depends on
!                   the unknown part -
!                   sequence_element(b, 1) succeeds on a sequence whose
!                   tail is UNKNOWN.
!
! NULL is not UNKNOWN here either: NULL ends a sequence, UNKNOWN
! withholds it.
!
! COMPLEXITY. size is O(n) and element(k) is O(k), and this is the
! semantic view, not the low-cost representation. Where repeated indexed
! access matters, compile to a contiguous representation.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module view_sequence

  use graph_fractal, only : graph, branch, &
       & BRANCH_NULL, BRANCH_UNKNOWN, BRANCH_KNOWN

  implicit none

  private
  public :: sequence_defined, sequence_num_elements
  public :: sequence_element, sequence_has
  public :: sequence_empty, sequence_first, sequence_rest

contains

  !===================================================================!
  ! Is the whole extent known: does the chain of cells reach NULL. False when
  ! it reaches UNKNOWN. A malformed cell is refused rather than
  ! reported, because it is not a state of the sequence.
  !===================================================================!

  logical function sequence_defined(b) result(known_extent)

    type(branch), intent(in) :: b

    type(graph), pointer :: cell

    known_extent = .false.
    if (b % status() .eq. BRANCH_UNKNOWN) return
    if (b % status() .eq. BRANCH_NULL) then
       known_extent = .true.
       return
    end if

    cell => b % known()
    do
       call require_cell(cell)
       if (cell % branch(2) % status() .eq. BRANCH_NULL) exit
       if (cell % branch(2) % status() .eq. BRANCH_UNKNOWN) return
       cell => cell % branch(2) % known()
    end do
    known_extent = .true.

  end function sequence_defined

  !===================================================================!
  ! The number of elements. The result depends on the whole chain of
  ! cells, so an unknown extent is refused.
  !===================================================================!

  integer function sequence_num_elements(b) result(n)

    type(branch), intent(in) :: b

    type(graph), pointer :: cell

    if (.not. sequence_defined(b)) then
       error stop 'view_sequence: the extent depends on an unknown tail'
    end if

    n = 0
    if (b % status() .eq. BRANCH_NULL) return

    cell => b % known()
    do
       n = n + 1
       if (cell % branch(2) % status() .eq. BRANCH_NULL) exit
       cell => cell % branch(2) % known()
    end do

  end function sequence_num_elements

  !===================================================================!
  ! The k-th element, counting from one. Only the first k cells are
  ! traversed, so an unknown tail beyond k is not an error.
  !===================================================================!

  function sequence_element(b, k) result(element)

    type(branch), intent(in) :: b
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
  ! Does the sequence contain this graph, by identity. Returned as
  ! soon as it is found; refused if the chain of cells reaches UNKNOWN
  ! first, because then the result depends on the unknown part.
  !
  ! This is one traversal. A caller looping sequence_element instead
  ! would be O(n^2), which is why membership is here and not derived
  ! outside.
  !===================================================================!

  logical function sequence_has(b, g) result(found)

    type(branch), intent(in) :: b
    type(graph)       , intent(in) :: g

    type(graph), pointer :: cell, element

    found = .false.
    if (b % status() .eq. BRANCH_NULL) return
    if (b % status() .eq. BRANCH_UNKNOWN) then
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
       if (cell % branch(2) % status() .eq. BRANCH_NULL) return
       if (cell % branch(2) % status() .eq. BRANCH_UNKNOWN) then
          error stop 'view_sequence: membership depends on an unknown tail'
       end if
       cell => cell % branch(2) % known()
    end do

  end function sequence_has

  !===================================================================!
  ! The recursive form: a sequence is empty, or it is a first
  ! element followed by the rest. A traversal written on these two
  ! traverses the chain of cells once, where indexing by position restarts from
  ! the head at every step.
  !===================================================================!

  logical function sequence_empty(b) result(empty)

    type(branch), intent(in) :: b

    empty = b % status() .eq. BRANCH_NULL

  end function sequence_empty

  !===================================================================!
  ! The first element. An empty or unknown sequence has none and
  ! stops the program, as does a malformed cell.
  !===================================================================!

  function sequence_first(b) result(element)

    type(branch), intent(in) :: b
    type(graph), pointer     :: element

    type(graph), pointer :: cell

    call require_reachable(b)
    cell => b % known()
    call require_cell(cell)
    element => cell % branch(1) % known()

  end function sequence_first

  !===================================================================!
  ! The sequence after the first element, refused on the same
  ! grounds. A copy of the branch stores the reference and owns
  ! nothing, which is what lets the result be returned by value.
  !===================================================================!

  function sequence_rest(b) result(rest)

    type(branch), intent(in) :: b
    type(branch)             :: rest

    type(graph), pointer :: cell

    call require_reachable(b)
    cell => b % known()
    call require_cell(cell)
    rest = cell % branch(2)

  end function sequence_rest

  !===================================================================!
  ! The two guards. require_cell refuses a malformed representation;
  ! require_reachable refuses a step that cannot be taken.
  !===================================================================!

  subroutine require_cell(cell)

    type(graph), intent(in) :: cell

    if (cell % branch(1) % status() .ne. BRANCH_KNOWN) then
       error stop 'view_sequence: a sequence cell contains a KNOWN element'
    end if

  end subroutine require_cell

  subroutine require_reachable(b)

    type(branch), intent(in) :: b

    if (b % status() .eq. BRANCH_NULL) then
       error stop 'view_sequence: the sequence has no such element'
    end if
    if (b % status() .eq. BRANCH_UNKNOWN) then
       error stop 'view_sequence: that element lies beyond an unknown tail'
    end if

  end subroutine require_reachable

end module view_sequence
