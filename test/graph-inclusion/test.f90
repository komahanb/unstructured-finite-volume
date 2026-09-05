!=====================================================================!
! THE INCLUSION SUITE
!
! The declared embedding, keyed on graph identity, over the graph +
! map foundation.
!
!     A and B have equal extensions and are two sets
!     S is declared into A, and into A alone
!     S <= A, and S is NOT <= B, though every member of S fits in B
!     T c--> S c--> A, so T <= A, transitively, and T is NOT <= B
!
! and the distinction the whole map exists to preserve:
!
!     extensional subset   a question about two extents
!     declared subobject   a question about what was declared
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program inclusion_suite

  use graph_fractal          , only : graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set          , only : set_map
  use map_inclusion    , only : inclusion_map, declared_subobject

  implicit none

  integer :: failures = 0

  write(*,'(1x,a)') "inclusion suite"

  !===================================================================!
  ! 1 . DECLARED PROVENANCE. The law, on the new foundation.
  !===================================================================!

  provenance_block: block

    type(graph), target :: a, b, s, t
    type(set_map)       :: sets
    type(inclusion_map) :: inc
    integer             :: k
    logical             :: fits

    call a % declare(); call b % declare()
    call s % declare(); call t % declare()

    call sets % bind(a, counted_set_representation(8))
    call sets % bind(b, counted_set_representation(8))
    call sets % bind(s, listed_set_representation([2, 5, 6]))
    call sets % bind(t, listed_set_representation([5]))

    call inc % include_in(s, a)
    call inc % include_in(t, s)

    call check('1  A and B are two sets, equal extensions notwithstanding', &
         & equivalent(a, b, sets) .and. .not. a % same_as(b))

    fits = .true.
    do k = 1, sets % num_members_of(s)
       fits = fits .and. sets % has(b, sets % member_of(s, k))
    end do
    call check('1  every member of S fits inside B', fits)

    call check('1  S <= S', declared_subobject(s, s, inc))
    call check('1  S <= A, where it was declared', &
         & declared_subobject(s, a, inc))
    call check('1  S is NOT <= B, though every member fits', &
         & .not. declared_subobject(s, b, inc))

    call check('1  T <= T', declared_subobject(t, t, inc))
    call check('1  T <= S', declared_subobject(t, s, inc))
    call check('1  T <= A, transitively', declared_subobject(t, a, inc))
    call check('1  T is NOT <= B: no chain reaches it', &
         & .not. declared_subobject(t, b, inc))

    call check('1  and the order points one way: A is not <= S', &
         & .not. declared_subobject(a, s, inc))

  end block provenance_block

  !===================================================================!
  ! 2 . THE TWO NOTIONS, SIDE BY SIDE.
  !
  ! S = {2,5,6} is extensionally inside every set that contains 2, 5 and
  ! 6, and is a declared subobject of exactly one of them.
  !===================================================================!

  distinction_block: block

    type(graph), target :: a, b, s
    type(set_map)       :: sets
    type(inclusion_map) :: inc
    integer             :: k
    logical             :: inside_a, inside_b

    call a % declare(); call b % declare(); call s % declare()

    call sets % bind(a, counted_set_representation(8))
    call sets % bind(b, listed_set_representation([2, 5, 6, 7]))
    call sets % bind(s, listed_set_representation([2, 5, 6]))

    call inc % include_in(s, a)

    inside_a = .true.; inside_b = .true.
    do k = 1, sets % num_members_of(s)
       inside_a = inside_a .and. sets % has(a, sets % member_of(s, k))
       inside_b = inside_b .and. sets % has(b, sets % member_of(s, k))
    end do

    call check('2  S is extensionally inside A and inside B', &
         & inside_a .and. inside_b)
    call check('2  but declared into A only', &
         & declared_subobject(s, a, inc) .and. &
         &  .not. declared_subobject(s, b, inc))
    call check('2  so containment does not imply an inclusion edge', &
         & inside_b .and. .not. declared_subobject(b, a, inc) &
         &          .and. .not. declared_subobject(b, s, inc))

  end block distinction_block

  !===================================================================!
  ! 3 . THE VALUE MAP IS THE IDENTITY, AND IS NOT STORED.
  !
  ! A member of S is the same value in A. What changes is where it
  ! is placed. So the inclusion stores an association and no table.
  !===================================================================!

  coordinates_block: block

    type(graph), target :: a, s
    type(set_map)       :: sets
    type(inclusion_map) :: inc
    integer             :: k
    logical             :: same_value, position_differs

    call a % declare(); call s % declare()
    call sets % bind(a, counted_set_representation(8))
    call sets % bind(s, listed_set_representation([2, 5, 6]))
    call inc % include_in(s, a)

    same_value       = .true.
    position_differs = .false.
    do k = 1, sets % num_members_of(s)
       associate (v => sets % member_of(s, k))
         same_value = same_value .and. sets % has(a, v)
         if (sets % index_in(a, v) /= sets % index_in(s, v)) &
              & position_differs = .true.
       end associate
    end do

    call check('3  every member of S is the same VALUE in A', same_value)
    call check('3  and is at a different POSITION', position_differs)
    call check('3  5 is 2nd in S and 5th in A', &
         & sets % index_in(s, 5) .eq. 2 .and. &
         &  sets % index_in(a, 5) .eq. 5)

    !----------------------------------------------------------------!
    ! The map stores the directed pair id(S) -> id(A) and evaluates
    ! the order predicate; it never rebuilds a graph object from a
    ! stored identity. The order is directed: S <= A and not A <= S.
    !----------------------------------------------------------------!

    call check('3  A is the declared ambient of S, and S is not of A', &
         & declared_subobject(s, a, inc) .and. .not. declared_subobject(a, s, inc))

  end block coordinates_block

  !===================================================================!

  if (failures .eq. 0) then
     print *, ''
     print *, ' ALL PROPOSITIONS HOLD'
  else
     print *, ''
     print *, ' FAILURES :', failures
     error stop 'inclusion: a proposition failed'
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

end program inclusion_suite
