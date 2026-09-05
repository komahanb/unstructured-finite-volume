!=====================================================================!
! THE IDENTITY-MAP LIFETIME SUITE
!
! Two maps are keyed on graph identity - set_map and inclusion_map -
! and this suite checks the one storage law they share:
!
!     an identity map OWNS ITS KEYS BY VALUE.
!     It references no graph object in order to recognise it.
!
! A map may return only values and still store a pointer INWARD, to
! the caller's graph variable; then the map outlives its own key and
! every lookup reads storage the binder already deallocated.
!
!                   WHAT WAS MEASURED, AND WHEN
!
! Before the storage law, both maps stored
!
!     type(graph), pointer :: element        set_map
!     type(graph), pointer :: part, ambient  inclusion_map
!
! and bind/include_in took their graphs with the TARGET attribute. The
! check below - bind inside a scope, deallocate the binder's graph,
! then query the map with a COPY storing the same token - measured:
!
!     native    describes = T, num_members_of = 10   read off freed storage
!     valgrind  Invalid read of size 4                in the row scan
!
! The native run returning CORRECT values is the failure: a freed page
! that has not yet been reused still reads as it did. The pointer
! itself established the dependency; the measurement states the cost.
!
! Under the storage law a row stores type(token), copied at bind. The
! graph dummies lost TARGET, which is why the binder below declares
! none: this file would not compile against the pointer-keyed maps,
! and that is the compile-time half of the proof.
!
!                    WHY A COPY IS THE RIGHT KEY
!
! Looking a map up with the very variable that bound it cannot
! distinguish the two designs - that variable is alive by
! construction. A copy stores the token and nothing else, so the map
! must recognise an identity rather than an address. That is the
! question the law is about.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program lifetime

  use graph_fractal           , only : graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set           , only : set_map
  use map_inclusion     , only : inclusion_map, declared_subobject

  implicit none

  integer :: failures = 0

  write(*,'(1x,a)') "identity map lifetime suite"

  !===================================================================!
  ! A . THE SET MAP OUTLIVES THE GRAPH THAT BOUND IT.
  !
  ! bind_and_die allocates its graphs and deallocates them on return,
  ! so the storage is released: a stack frame out of scope may remain
  ! intact, a deallocated block is released.
  !===================================================================!

  set_map_block: block

    type(set_map) :: sets
    type(graph)   :: a, s
    integer, allocatable :: v(:)

    call bind_and_die(sets, a, s)

    call check('A  the map still describes a set whose binder is gone', &
         & sets % describes(a) .and. sets % describes(s))

    call check('A  and returns its size', sets % num_members_of(a) .eq. 10)

    call check('A  and its members, positions and membership', &
         & sets % member_of(s, 2) .eq. 5 .and. &
         &  sets % index_in(s, 6) .eq. 3   .and. &
         &  sets % has(s, 5)            .and. &
         &  .not. sets % has(s, 4))

    call sets % members_of(s, v)
    call check('A  and enumerates them', &
         & size(v) .eq. 3 .and. all(v .eq. [2, 5, 6]))

    !----------------------------------------------------------------!
    ! A set the map never bound is still undescribed - the copy
    ! stores a token, not a binding.
    !----------------------------------------------------------------!

    outsider: block
      type(graph) :: outsider
      call outsider % declare()
      call check('A  an unbound set is still undescribed', &
           & .not. sets % describes(outsider))
    end block outsider

  end block set_map_block

  !===================================================================!
  ! B . THE INCLUSION MAP OUTLIVES THE GRAPHS THAT DECLARED IT.
  !
  ! The chain T c--> S c--> A is declared inside a scope that then
  ! deallocates every graph. The subobject order, traversed entirely
  ! on identities, must still evaluate - including the negative,
  ! which is the value a dangling scan is most likely to get wrong.
  ! The order predicate is the map's one query: a part is included
  ! iff it is a declared subobject of some other identity.
  !===================================================================!

  inclusion_block: block

    type(inclusion_map) :: inc
    type(graph)         :: a, b, s, t

    call declare_and_die(inc, a, b, s, t)

    call check('B  the declared ambients are still recognised', &
         & declared_subobject(s, a, inc) .and. declared_subobject(t, s, inc))

    call check('B  S <= A, after every declaring variable is gone', &
         & declared_subobject(s, a, inc))

    call check('B  T <= A, transitively', declared_subobject(t, a, inc))

    call check('B  T is NOT <= B: the negative survives too', &
         & .not. declared_subobject(t, b, inc))

    call check('B  and the direct edge is evaluated by identity', &
         & declared_subobject(s, a, inc) .and. .not. declared_subobject(s, b, inc))

    call check('B  while an undeclared part has no ambient', &
         & .not. declared_subobject(a, b, inc) .and. &
         &  .not. declared_subobject(a, s, inc) .and. &
         &  .not. declared_subobject(a, t, inc))

  end block inclusion_block

  !===================================================================!
  ! C . A MAP IS A VALUE. Intrinsic assignment deep-copies the rows,
  ! keys included, so a copy evaluates the same queries and the
  ! original is untouched by the copy's growth.
  !===================================================================!

  value_block: block

    type(set_map) :: sets, twin
    type(graph)   :: a, s, extra

    call bind_and_die(sets, a, s)

    twin = sets
    call extra % declare()
    call twin % bind(extra, counted_set_representation(4))

    call check('C  a copied map describes the keys it copied', &
         & twin % describes(a) .and. twin % num_members_of(a) .eq. 10)

    call check('C  and growth of the copy does not reach the original', &
         & twin % describes(extra) .and. .not. sets % describes(extra))

  end block value_block

  !===================================================================!

  if (failures .eq. 0) then
     print *, ''
     print *, ' ALL PROPOSITIONS HOLD'
  else
     print *, ''
     print *, ' FAILURES :', failures
     error stop 'lifetime: a proposition failed'
  end if

contains

  !===================================================================!
  ! Bind two sets, return copies of their identities, and deallocate
  ! every graph object the map was built from. No TARGET anywhere:
  ! under the storage law, binding requires none.
  !===================================================================!

  subroutine bind_and_die(m, a_key, s_key)

    type(set_map), intent(out) :: m
    type(graph)  , intent(out) :: a_key, s_key

    type(graph), allocatable :: a, s

    allocate(a, s)
    call a % declare(); call s % declare()

    call m % bind(a, counted_set_representation(10))
    call m % bind(s, listed_set_representation([2, 5, 6]))

    a_key = a
    s_key = s

    deallocate(a, s)

  end subroutine bind_and_die

  !===================================================================!
  ! Declare T c--> S c--> A, return the four identities, and
  ! deallocate every graph the declarations named.
  !===================================================================!

  subroutine declare_and_die(m, a_key, b_key, s_key, t_key)

    type(inclusion_map), intent(out) :: m
    type(graph)        , intent(out) :: a_key, b_key, s_key, t_key

    type(graph), allocatable :: a, b, s, t

    allocate(a, b, s, t)
    call a % declare(); call b % declare()
    call s % declare(); call t % declare()

    call m % include_in(s, a)
    call m % include_in(t, s)

    a_key = a; b_key = b; s_key = s; t_key = t

    deallocate(a, b, s, t)

  end subroutine declare_and_die

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

end program lifetime
