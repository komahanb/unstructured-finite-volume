!=====================================================================!
! THE IDENTITY-MAP LIFETIME SUITE
!
! Two maps are keyed on graph identity - set_map and inclusion_map -
! and this suite holds the one storage law they share:
!
!     an identity map OWNS ITS KEYS BY VALUE.
!     It borrows no graph object merely to recognize it.
!
! That is stronger than "the API lends nothing". A map may hand back
! only values and still hold a pointer INWARD, to the caller's graph
! variable; then the map outlives its own key and every lookup reads
! storage the binder already gave back.
!
!                   WHAT WAS MEASURED, AND WHEN
!
! Before the gate, both maps stored
!
!     type(graph), pointer :: element        set_map
!     type(graph), pointer :: part, ambient  inclusion_map
!
! and bind/include_in took their graphs with the TARGET attribute. The
! probe below - bind inside a scope, let the binder's graph die, then
! ask the map about a COPY carrying the same token - measured:
!
!     native    describes = T, size_of = 10   answers, off freed store
!     valgrind  Invalid read of size 4        in row_of, the scan
!
! The native run answering CORRECTLY is the danger, not the comfort:
! a freed page that has not yet been reused still reads as it did. No
! crash was needed to establish the dependency; the pointer itself
! established it. The measurement only says what the dependency costs
! when it is paid.
!
! After the gate a row stores type(token), copied at bind. The graph
! dummies lost TARGET, which is why the binder below declares none:
! this file would not compile against the pointer-keyed maps, and that
! is the compile-time half of the proof.
!
!                    WHY A COPY IS THE RIGHT PROBE
!
! Looking a map up with the very variable that bound it cannot tell
! the two designs apart - that variable is alive by construction. A
! copy carries the token and nothing else, so it asks the map to
! recognize an identity rather than to recognize an address. That is
! precisely the question the law is about.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program lifetime

  use fractal_graph           , only : graph
  use graph_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use graph_set_map           , only : set_map
  use graph_inclusion_map     , only : inclusion_map, declared_subobject

  implicit none

  integer :: failures = 0

  write(*,'(1x,a)') "identity map lifetime suite"

  !===================================================================!
  ! A . THE SET MAP OUTLIVES THE GRAPH THAT BOUND IT.
  !
  ! bind_and_die allocates its graph on the heap and lets it be
  ! deallocated on return, so the storage is genuinely given back -
  ! a dead stack frame merely looks dead, a freed block is dead.
  !===================================================================!

  set_map_block: block

    type(set_map) :: sets
    type(graph)   :: a, s
    integer, allocatable :: v(:)

    call bind_and_die(sets, a, s)

    call check('A  the map still describes a set whose binder is gone', &
         & sets % describes(a) .and. sets % describes(s))

    call check('A  and answers its size', sets % size_of(a) .eq. 10)

    call check('A  and its members, positions and membership', &
         & sets % member_of(s, 2) .eq. 5 .and. &
         &  sets % index_in(s, 6) .eq. 3   .and. &
         &  sets % has_in(s, 5)            .and. &
         &  .not. sets % has_in(s, 4))

    call sets % members_of(s, v)
    call check('A  and enumerates them', &
         & size(v) .eq. 3 .and. all(v .eq. [2, 5, 6]))

    !----------------------------------------------------------------!
    ! And a set the map never heard of is still refused - the copy
    ! carries a token, not a licence.
    !----------------------------------------------------------------!

    stranger: block
      type(graph) :: outsider
      call outsider % declare()
      call check('A  an unbound set is still undescribed', &
           & .not. sets % describes(outsider))
    end block stranger

  end block set_map_block

  !===================================================================!
  ! B . THE INCLUSION MAP OUTLIVES THE GRAPHS THAT DECLARED IT.
  !
  ! The chain T c--> S c--> A is declared inside a scope that then
  ! dies. The subobject order, walked entirely on identities, must
  ! still answer - including the negative, which is the answer a
  ! dangling scan is most likely to get wrong.
  !===================================================================!

  inclusion_block: block

    type(inclusion_map) :: inc
    type(graph)         :: a, b, s, t

    call declare_and_die(inc, a, b, s, t)

    call check('B  the declared ambient is still recognized', &
         & inc % included(s) .and. inc % included(t))

    call check('B  S <= A, after every declaring variable is gone', &
         & declared_subobject(s, a, inc))

    call check('B  T <= A, transitively', declared_subobject(t, a, inc))

    call check('B  T is NOT <= B: the negative survives too', &
         & .not. declared_subobject(t, b, inc))

    call check('B  and the direct edge answers by identity', &
         & inc % declared_into(s, a) .and. .not. inc % declared_into(s, b))

    call check('B  while an undeclared part has no ambient', &
         & .not. inc % included(a))

  end block inclusion_block

  !===================================================================!
  ! C . A MAP IS A VALUE. Intrinsic assignment deep-copies the rows,
  ! keys included, so a copy answers the same questions and the
  ! original is untouched by the copy's growth.
  !===================================================================!

  value_block: block

    type(set_map) :: sets, twin
    type(graph)   :: a, s, extra

    call bind_and_die(sets, a, s)

    twin = sets
    call extra % declare()
    call twin % bind(extra, counted_set_representation(4))

    call check('C  a copied map answers for the keys it copied', &
         & twin % describes(a) .and. twin % size_of(a) .eq. 10)

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
  ! Bind two sets, hand back copies of their identities, and destroy
  ! every graph object the map was built from. No TARGET anywhere:
  ! after the gate, binding does not ask for one.
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
  ! Declare T c--> S c--> A, hand back the four identities, and
  ! destroy every graph the declarations named.
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

end program lifetime
