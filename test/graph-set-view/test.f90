!=====================================================================!
! THE SET FOUNDATION SUITE
!
! A finite set is a GRAPH and a REPRESENTATION, associated by a MAP and
! read through a VIEW. Four roles, and the laws that keep them apart:
!
!     identity is the graph's, and only the graph's
!     the extension is the representation's, and costs O(N_extent)
!     the branches stay NULL, at any cardinality
!     nothing is lent, so nothing can dangle
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program set_foundation

  use graph_fractal          , only : graph, null_branch, GRAPH_NULL
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set          , only : set_map
  use view_set         , only : set_defined, set_size, set_member, &
       & set_members, set_has, set_local_index, set_equivalent

  implicit none

  integer :: failures = 0

  write(*,'(1x,a)') "set foundation suite"

  !===================================================================!
  ! 1 . IDENTITY IS NOT EXTENSIONAL EQUALITY.
  !
  ! Two sets over 1..8, and two empty sets. Equal members, different
  ! sets. The representations take no part: they carry no identity to
  ! take part with.
  !===================================================================!

  identity_block: block

    type(graph), target :: a, b, empty_a, empty_b
    type(set_map)       :: m

    call a % declare(); call b % declare()
    call empty_a % declare(); call empty_b % declare()

    call m % bind(a, counted_set_representation(8))
    call m % bind(b, counted_set_representation(8))
    call m % bind(empty_a, counted_set_representation(0))
    call m % bind(empty_b, listed_set_representation([integer ::]))

    call check('1  two sets over 1..8 are extensionally equal', &
         & set_equivalent(a, b, m))
    call check('1  and are NOT the same set', .not. a % same_as(b))

    call check('1  two empty sets are extensionally equal', &
         & set_equivalent(empty_a, empty_b, m))
    call check('1  and are NOT the same set', &
         & .not. empty_a % same_as(empty_b))
    call check('1  even when described by different representations', &
         & set_size(empty_a, m) .eq. 0 .and. set_size(empty_b, m) .eq. 0)

    call check('1  identity is the graph''s alone', &
         & a % same_as(a) .and. b % same_as(b))

  end block identity_block

  !===================================================================!
  ! 2 . COUNTED, AT 10^9.
  !
  ! One graph, one integer. No member object, and both branches NULL -
  ! asserted, because the whole architecture rests on it.
  !===================================================================!

  counted_block: block

    type(graph), target  :: cells
    type(set_map)        :: m
    integer, allocatable :: v(:)

    call cells % declare()
    cells % branch(1) = null_branch()
    cells % branch(2) = null_branch()

    call m % bind(cells, counted_set_representation(1000000000))

    call check('2  |cells| = 10^9, and no member object exists', &
         & set_size(cells, m) .eq. 1000000000)
    call check('2  membership is one comparison', &
         & set_has(cells, m, 999999999) .and. &
         &  .not. set_has(cells, m, 1000000001))
    call check('2  and position is the representation''s numbering', &
         & set_local_index(cells, m, 7) .eq. 7 .and. &
         &  set_member(cells, m, 7) .eq. 7)
    call check('2  the graph carries no extension: both branches NULL', &
         & cells % branch(1) % status() .eq. GRAPH_NULL .and. &
         &  cells % branch(2) % status() .eq. GRAPH_NULL)

    ! And at a size a test may enumerate, the roster is the identity map.
    block
      type(graph), target :: small
      call small % declare()
      call m % bind(small, counted_set_representation(4))
      call set_members(small, m, v)
      call check('2  members enumerate 1..n', &
           & size(v) .eq. 4 .and. all(v .eq. [1, 2, 3, 4]))
    end block

  end block counted_block

  !===================================================================!
  ! 3 . LISTED, AND THE TWO COORDINATE SYSTEMS.
  !
  ! S = {2,5,6} and A = 1..8. The member VALUE is shared; the POSITION
  ! is each representation's own. Both enumeration laws hold inside
  ! each, independently.
  !===================================================================!

  listed_block: block

    type(graph), target  :: s, a
    type(set_map)        :: m
    integer, allocatable :: v(:)
    integer              :: k
    logical              :: both_ways

    call s % declare(); call a % declare()
    call m % bind(a, counted_set_representation(8))
    call m % bind(s, listed_set_representation([2, 5, 6]))

    call check('3  |S| = 3, and its members are its own values', &
         & set_size(s, m) .eq. 3)
    call set_members(s, m, v)
    call check('3  enumerated in declaration order', &
         & all(v .eq. [2, 5, 6]))

    call check('3  the value 5 stands 2nd in S', &
         & set_member(s, m, 2) .eq. 5 .and. set_local_index(s, m, 5) .eq. 2)
    call check('3  and 5th in A - one value, two positions', &
         & set_member(a, m, 5) .eq. 5 .and. set_local_index(a, m, 5) .eq. 5)

    both_ways = .true.
    do k = 1, set_size(s, m)
       both_ways = both_ways .and. &
            & set_local_index(s, m, set_member(s, m, k)) .eq. k
       both_ways = both_ways .and. &
            & set_member(a, m, set_local_index(a, m, set_member(s, m, k))) &
            &   .eq. set_member(s, m, k)
    end do
    call check('3  member(local_index(v)) = v holds inside EACH', both_ways)

    call check('3  an outsider stands nowhere', &
         & set_local_index(s, m, 3) .eq. 0 .and. .not. set_has(s, m, 3))

    ! A listed representation describes a set with no ambient at all.
    block
      type(graph), target :: loose
      call loose % declare()
      call m % bind(loose, listed_set_representation([30, 10, 20]))
      call check('3  and a listed set needs no ambient to exist', &
           & set_size(loose, m) .eq. 3 .and. &
           &  set_member(loose, m, 1) .eq. 30 .and. &
           &  set_local_index(loose, m, 20) .eq. 3)
    end block

    ! Repetition collapses: a representation lists each member once.
    block
      type(graph), target :: dup
      call dup % declare()
      call m % bind(dup, listed_set_representation([4, 9, 4, 9, 1]))
      call set_members(dup, m, v)
      call check('3  a repeated value is listed once, first place kept', &
           & size(v) .eq. 3 .and. all(v .eq. [4, 9, 1]))
    end block

  end block listed_block

  !===================================================================!
  ! 4 . THE MAP LENDS NOTHING.
  !
  ! Growth relocates its rows freely, because no caller holds anything
  ! that could point into them; and a copy is a deep copy, because the
  ! rows own their representations rather than referencing them.
  !===================================================================!

  storage_block: block

    type(graph), target :: a, b, c
    type(set_map)       :: m, copy

    call a % declare(); call b % declare(); call c % declare()

    call m % bind(a, counted_set_representation(4))
    call check('4  a is described', set_size(a, m) .eq. 4)

    call m % bind(b, counted_set_representation(9))          ! growth
    call check('4  and still is after the row array grows', &
         & set_size(a, m) .eq. 4 .and. set_size(b, m) .eq. 9)

    copy = m
    call copy % bind(c, listed_set_representation([7]))
    call check('4  a copy is independent: c is in the copy', &
         & set_defined(c, copy) .and. .not. set_defined(c, m))
    call check('4  and the original''s answers are unchanged', &
         & set_size(a, m) .eq. 4 .and. set_size(a, copy) .eq. 4)

  end block storage_block

  !===================================================================!

  if (failures .eq. 0) then
     print *, ''
     print *, ' ALL PROPOSITIONS HOLD'
  else
     print *, ''
     print *, ' FAILURES :', failures
     error stop 'set_foundation: a proposition failed'
  end if

contains

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

end program set_foundation
