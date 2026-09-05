!=====================================================================!
! THE SET FOUNDATION SUITE
!
! A finite set is a GRAPH and a REPRESENTATION, associated by a MAP.
! Every extent query is the map's own procedure. Three roles, and the
! laws that keep them apart:
!
!     identity is the graph's, and only the graph's
!     the extension is the representation's, and costs O(N_extent)
!     the branches stay NULL, at any cardinality
!     no reference is returned, so no reference can dangle
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program set_foundation

  use graph_fractal          , only : graph, null_branch, BRANCH_NULL
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set          , only : set_map

  implicit none

  integer :: failures = 0

  write(*,'(1x,a)') "set foundation suite"

  !===================================================================!
  ! 1 . IDENTITY IS NOT EXTENSIONAL EQUALITY.
  !
  ! Two sets over 1..8, and two empty sets. Equal members, different
  ! sets. The representations take no part: they store no identity.
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
         & equivalent(a, b, m))
    call check('1  and are NOT the same set', .not. a % same_as(b))

    call check('1  two empty sets are extensionally equal', &
         & equivalent(empty_a, empty_b, m))
    call check('1  and are NOT the same set', &
         & .not. empty_a % same_as(empty_b))
    call check('1  even when described by different representations', &
         & m % num_members_of(empty_a) .eq. 0 .and. m % num_members_of(empty_b) .eq. 0)

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
         & m % num_members_of(cells) .eq. 1000000000)
    call check('2  membership is one comparison', &
         & m % has(cells, 999999999) .and. &
         &  .not. m % has(cells, 1000000001))
    call check('2  and position is the representation''s numbering', &
         & m % index_in(cells, 7) .eq. 7 .and. &
         &  m % member_of(cells, 7) .eq. 7)
    call check('2  the graph stores no extension: both branches NULL', &
         & cells % branch(1) % status() .eq. BRANCH_NULL .and. &
         &  cells % branch(2) % status() .eq. BRANCH_NULL)

    ! At a size a test may enumerate, the member list is 1..n.
    block
      type(graph), target :: small
      call small % declare()
      call m % bind(small, counted_set_representation(4))
      call m % members_of(small, v)
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
         & m % num_members_of(s) .eq. 3)
    call m % members_of(s, v)
    call check('3  enumerated in declaration order', &
         & all(v .eq. [2, 5, 6]))

    call check('3  the value 5 is 2nd in S', &
         & m % member_of(s, 2) .eq. 5 .and. m % index_in(s, 5) .eq. 2)
    call check('3  and 5th in A - one value, two positions', &
         & m % member_of(a, 5) .eq. 5 .and. m % index_in(a, 5) .eq. 5)

    both_ways = .true.
    do k = 1, m % num_members_of(s)
       both_ways = both_ways .and. &
            & m % index_in(s, m % member_of(s, k)) .eq. k
       both_ways = both_ways .and. &
            & m % member_of(a, m % index_in(a, m % member_of(s, k))) &
            &   .eq. m % member_of(s, k)
    end do
    call check('3  member(local_index(v)) = v is true inside EACH', both_ways)

    call check('3  an outsider is at no position', &
         & m % index_in(s, 3) .eq. 0 .and. .not. m % has(s, 3))

    ! A listed representation describes a set with no ambient at all.
    block
      type(graph), target :: loose
      call loose % declare()
      call m % bind(loose, listed_set_representation([30, 10, 20]))
      call check('3  and a listed set needs no ambient to exist', &
           & m % num_members_of(loose) .eq. 3 .and. &
           &  m % member_of(loose, 1) .eq. 30 .and. &
           &  m % index_in(loose, 20) .eq. 3)
    end block

    ! Repetition collapses: a representation lists each member once.
    block
      type(graph), target :: dup
      call dup % declare()
      call m % bind(dup, listed_set_representation([4, 9, 4, 9, 1]))
      call m % members_of(dup, v)
      call check('3  a repeated value is listed once, at its first position', &
           & size(v) .eq. 3 .and. all(v .eq. [4, 9, 1]))
    end block

  end block listed_block

  !===================================================================!
  ! 4 . THE MAP RETURNS NO REFERENCE.
  !
  ! Growth relocates its rows freely, because no caller has a pointer
  ! into them; and a copy is a deep copy, because the rows own their
  ! representations rather than referencing them.
  !===================================================================!

  storage_block: block

    type(graph), target :: a, b, c
    type(set_map)       :: m, copy

    call a % declare(); call b % declare(); call c % declare()

    call m % bind(a, counted_set_representation(4))
    call check('4  a is described', m % num_members_of(a) .eq. 4)

    call m % bind(b, counted_set_representation(9))          ! growth
    call check('4  and still is after the row array grows', &
         & m % num_members_of(a) .eq. 4 .and. m % num_members_of(b) .eq. 9)

    copy = m
    call copy % bind(c, listed_set_representation([7]))
    call check('4  a copy is independent: c is in the copy', &
         & copy % describes(c) .and. .not. m % describes(c))
    call check('4  and the original''s extents are unchanged', &
         & m % num_members_of(a) .eq. 4 .and. copy % num_members_of(a) .eq. 4)

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

  !===================================================================!
  ! Extensional equality: the same members, in the same enumeration
  ! order, under one map.
  !===================================================================!

  logical function equivalent(a, b, m) result(equal)

    type(graph)  , intent(in) :: a, b
    type(set_map), intent(in) :: m

    integer, allocatable :: va(:), vb(:)

    call m % members_of(a, va)
    call m % members_of(b, vb)

    equal = size(va) .eq. size(vb)
    if (equal) equal = all(va .eq. vb)

  end function equivalent

  subroutine check(label, passes)

    character(len=*), intent(in) :: label
    logical         , intent(in) :: passes

    if (passes) then
       print *, ' PASS : ', label
    else
       print *, ' FAIL : ', label
       failures = failures + 1
    end if

  end subroutine check

end program set_foundation
