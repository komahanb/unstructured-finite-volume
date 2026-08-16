!=====================================================================!
! THE SEQUENCE VIEW SUITE
!
! The representation under test:
!
!     branch    NULL          the empty sequence
!               KNOWN -> cell a nonempty sequence
!               UNKNOWN       the sequence is not known
!
!     cell      branch(1) = KNOWN -> element
!               branch(2) = the rest, again a sequence branch
!
! The kernel is unchanged and knows none of these words.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test

  use fractal_graph      , only : graph, GRAPH_NULL, GRAPH_UNKNOWN, &
       & null_branch, unknown_branch, known_branch
  use graph_sequence_view, only : sequence_defined, sequence_size, &
       & sequence_element, sequence_contains

  implicit none

  integer :: failures = 0

  write(*,'(1x,a)') "sequence view suite"

  !===================================================================!
  ! The empty sequence. It has no cell, and needs no graph.
  !===================================================================!

  empty_block: block

    type(graph), target :: holder, stranger

    call holder % declare(); call stranger % declare()
    holder % branch(1) = null_branch()

    call check('empty: NULL is the empty sequence, not a missing one', &
         & sequence_defined(holder % branch(1)))
    call check('empty: size = 0', sequence_size(holder % branch(1)) .eq. 0)
    call check('empty: contains nothing', &
         & .not. sequence_contains(holder % branch(1), stranger))

  end block empty_block

  !===================================================================!
  ! Singleton, length two, length five.
  !===================================================================!

  lengths_block: block

    type(graph), target  :: holder
    type(graph), target  :: cell(5), elem(5)
    type(graph), pointer :: e
    integer              :: i
    logical              :: ok

    call holder % declare()
    do i = 1, 5
       call cell(i) % declare(); call elem(i) % declare()
    end do

    call wire(cell, elem, 1, .true.)
    holder % branch(1) = known_branch(cell(1))
    e => sequence_element(holder % branch(1), 1)
    call check('singleton: size = 1, element(1) is the element', &
         & sequence_size(holder % branch(1)) .eq. 1 .and. e % same_as(elem(1)))

    call wire(cell, elem, 2, .true.)
    holder % branch(1) = known_branch(cell(1))
    e => sequence_element(holder % branch(1), 2)
    call check('length 2: size = 2, element(2) is the second element', &
         & sequence_size(holder % branch(1)) .eq. 2 .and. e % same_as(elem(2)))

    call wire(cell, elem, 5, .true.)
    holder % branch(1) = known_branch(cell(1))
    call check('length n: size = 5', sequence_size(holder % branch(1)) .eq. 5)

    ok = .true.
    do i = 1, 5
       e => sequence_element(holder % branch(1), i)
       ok = ok .and. e % same_as(elem(i))
    end do
    call check('length n: element(k) is the k-th element, k = 1..5', ok)

    call check('length n: membership holds for a member', &
         & sequence_contains(holder % branch(1), elem(3)))
    call check('length n: and fails for a non-member', &
         & .not. sequence_contains(holder % branch(1), holder))

  end block lengths_block

  !===================================================================!
  ! An UNKNOWN holder: the sequence itself is not known. That is not
  ! the empty sequence, and it is not a malformed one.
  !===================================================================!

  unknown_holder_block: block

    type(graph), target :: holder

    call holder % declare()
    holder % branch(1) = unknown_branch()

    call check('unknown holder: the extent is not defined', &
         & .not. sequence_defined(holder % branch(1)))
    call check('unknown holder: UNKNOWN is not NULL', &
         & holder % branch(1) % status() .ne. GRAPH_NULL .and. &
         & holder % branch(1) % status() .eq. GRAPH_UNKNOWN)

  end block unknown_holder_block

  !===================================================================!
  ! An UNKNOWN tail. The known prefix is still readable: only an
  ! answer that depends on the unknown part is refused.
  !===================================================================!

  unknown_tail_block: block

    type(graph), target  :: holder
    type(graph), target  :: cell(3), elem(3)
    type(graph), pointer :: e
    integer              :: i

    call holder % declare()
    do i = 1, 3
       call cell(i) % declare(); call elem(i) % declare()
    end do

    call wire(cell, elem, 3, .false.)          ! last tail is UNKNOWN
    holder % branch(1) = known_branch(cell(1))

    call check('unknown tail: the extent is not defined', &
         & .not. sequence_defined(holder % branch(1)))

    e => sequence_element(holder % branch(1), 1)
    call check('unknown tail: element(1) is still answered', e % same_as(elem(1)))
    e => sequence_element(holder % branch(1), 3)
    call check('unknown tail: so is element(3), the last known one', &
         & e % same_as(elem(3)))

    call check('unknown tail: membership holds for a member of the prefix', &
         & sequence_contains(holder % branch(1), elem(2)))

  end block unknown_tail_block

  !===================================================================!

  if (failures .eq. 0) then
     print *, ''
     print *, ' ALL PROPOSITIONS HOLD'
  else
     print *, ''
     print *, ' FAILURES :', failures
     error stop 'test: a proposition failed'
  end if

contains

  !===================================================================!
  ! Wire cell(1..n) into a sequence of elem(1..n). Closed sequences
  ! end in NULL; open ones end in UNKNOWN.
  !===================================================================!

  subroutine wire(cell, elem, n, closed)

    type(graph), target, intent(inout) :: cell(:)
    type(graph), target, intent(inout) :: elem(:)
    integer            , intent(in)    :: n
    logical            , intent(in)    :: closed

    integer :: i

    do i = 1, n
       cell(i) % branch(1) = known_branch(elem(i))
       if (i .lt. n) then
          cell(i) % branch(2) = known_branch(cell(i + 1))
       else if (closed) then
          cell(i) % branch(2) = null_branch()
       else
          cell(i) % branch(2) = unknown_branch()
       end if
    end do

  end subroutine wire

  subroutine check(label, ok)

    character(len=*), intent(in) :: label
    logical         , intent(in) :: ok

    if (ok) then
       print *, ' PASS : ', label
    else
       print *, ' FAIL : ', label
       failures = failures + 1
    end if

  end subroutine check

end program test
