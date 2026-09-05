!=====================================================================!
! PROTOTYPE . WHAT IS A RELATION IN THE FRACTAL SYSTEM
!
! Three things are separated and tested apart:
!
!     R                  the semantic relation, an identity
!     (A_1,...,A_k)      the ordered signature, small
!     extension(R)       the tuples, large
!
! The signature is built with view_sequence over set graphs -
! repeated domains, mixed domains, arity three. The extension is left
! in the production representations, which already store integers.
!
! The load-bearing experiment is TUPLE EQUALITY: if a tuple were ever
! given a graph of its own, tuple-holder identity would NOT be tuple
! equality, and a relation would stop being a set. That trap is built
! here so it is visible rather than argued.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module relation_prototype

  use graph_fractal, only : graph, null_branch, known_branch
  use view_sequence, only : sequence_num_elements, sequence_element

  implicit none

  private
  public :: same_tuple, transpose_role, apply_role

  !===================================================================!
  ! A transpose expressed as a ROLE PERMUTATION over one base, rather
  ! than as a new relation. Composing two of them returns the identity
  ! permutation, so the pair (base, role) is again the base itself -
  ! involution at the level of the SEMANTIC OBJECT, not merely of the
  ! extension.
  !===================================================================!

  type, public :: role_view
     type(graph), pointer :: base => null()
     integer              :: role(2) = [1, 2]
  end type role_view

contains

  !===================================================================!
  ! Tuple equality is COMPONENT IDENTITY, slot by slot. Never the
  ! identity of whatever object happens to hold the tuple.
  !===================================================================!

  logical function same_tuple(x, y) result(equal)

    type(graph), intent(in) :: x, y        ! two tuple-holder graphs

    type(graph), pointer :: cx, cy
    integer              :: k, n

    equal = .false.
    n = sequence_num_elements(x % branch(1))
    if (n .ne. sequence_num_elements(y % branch(1))) return

    equal = .true.
    do k = 1, n
       cx => sequence_element(x % branch(1), k)
       cy => sequence_element(y % branch(1), k)
       equal = equal .and. cx % same_as(cy)
    end do

  end function same_tuple

  type(role_view) function transpose_role(v) result(t)

    type(role_view), intent(in) :: v

    t % base => v % base
    t % role = [v % role(2), v % role(1)]

  end function transpose_role

  pure integer function apply_role(v, slot_index) result(actual)

    type(role_view), intent(in) :: v
    integer        , intent(in) :: slot_index

    actual = v % role(slot_index)

  end function apply_role

end module relation_prototype

program relation_map

  use graph_fractal        , only : graph, null_branch, known_branch
  use view_sequence  , only : sequence_num_elements, sequence_element
  use graph_fractal           , only : graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set           , only : set_map
  use map_label         , only : label_map
  use map_inclusion     , only : inclusion_map, declared_subobject
  use relation_finitary       , only : relation, stored_relation
  use relation_binary, only : binary_relation, csr_relation, &
       & transposed_relation, transpose_of
  use view_relational, only : relational_binding
  use relation_prototype   , only : same_tuple, role_view, &
       & transpose_role, apply_role

  implicit none

  integer :: failures = 0

  write(*,'(1x,a)') "relation prototype"

  !===================================================================!
  ! 1 . THE SIGNATURE IS A SEQUENCE, and arity is small.
  !
  !     R <= A x A     repeated domain
  !     R <= A x B     mixed domains
  !     R <= A x B x C arity three
  !
  ! One mechanism serves all three. Arity two earns no ontology.
  !===================================================================!

  signature_block: block

    type(graph), target  :: a, b, c
    type(graph), target  :: raa, rab, rabc            ! three relations
    type(graph), target  :: saa(2), sab(2), sabc(3)   ! their signatures
    type(graph), pointer :: d
    integer              :: k
    logical              :: ok

    call a % declare(); call b % declare(); call c % declare()
    call raa % declare(); call rab % declare(); call rabc % declare()
    do k = 1, 2
       call saa(k) % declare(); call sab(k) % declare()
    end do
    do k = 1, 3
       call sabc(k) % declare()
    end do

    ! A x A: the SAME set graph stands in both slots.
    saa(1) % branch(1) = known_branch(a)
    saa(1) % branch(2) = known_branch(saa(2))
    saa(2) % branch(1) = known_branch(a)
    raa % branch(1) = known_branch(saa(1))

    call check('1  A x A: arity 2, one repeated domain', &
         & sequence_num_elements(raa % branch(1)) .eq. 2)
    d => sequence_element(raa % branch(1), 1)
    ok = d % same_as(a)
    d => sequence_element(raa % branch(1), 2)
    call check('1  and both slots name the same set, by identity', &
         & ok .and. d % same_as(a))

    ! A x B.
    sab(1) % branch(1) = known_branch(a)
    sab(1) % branch(2) = known_branch(sab(2))
    sab(2) % branch(1) = known_branch(b)
    rab % branch(1) = known_branch(sab(1))

    d => sequence_element(rab % branch(1), 2)
    call check('1  A x B: mixed domains, same mechanism', &
         & sequence_num_elements(rab % branch(1)) .eq. 2 .and. &
         &  d % same_as(b) .and. .not. d % same_as(a))

    ! A x B x C.
    sabc(1) % branch(1) = known_branch(a)
    sabc(1) % branch(2) = known_branch(sabc(2))
    sabc(2) % branch(1) = known_branch(b)
    sabc(2) % branch(2) = known_branch(sabc(3))
    sabc(3) % branch(1) = known_branch(c)
    rabc % branch(1) = known_branch(sabc(1))

    ok = sequence_num_elements(rabc % branch(1)) .eq. 3
    d => sequence_element(rabc % branch(1), 2); ok = ok .and. d % same_as(b)
    d => sequence_element(rabc % branch(1), 3); ok = ok .and. d % same_as(c)
    call check('1  A x B x C: arity 3 is a longer sequence, not a case', ok)
    call check('1  the signature costs O(k) cells, and k is small', &
         & size(saa) + size(sab) + size(sabc) .eq. 7)

  end block signature_block

  !===================================================================!
  ! 2 . TUPLE EQUALITY IS COMPONENT IDENTITY.
  !
  ! Two tuple-holder graphs are built independently, both spelling
  ! (a, b). They are DIFFERENT graphs. If a relation deduped on holder
  ! identity it would hold (a,b) twice and stop being a set.
  !===================================================================!

  tuple_block: block

    type(graph), target :: a, b
    type(graph), target :: h1, c1(2), h2, c2(2), h3, c3(2)
    integer             :: k

    call a % declare(); call b % declare()
    call h1 % declare(); call h2 % declare(); call h3 % declare()
    do k = 1, 2
       call c1(k) % declare(); call c2(k) % declare(); call c3(k) % declare()
    end do

    ! h1 = (a, b)
    c1(1) % branch(1) = known_branch(a); c1(1) % branch(2) = known_branch(c1(2))
    c1(2) % branch(1) = known_branch(b)
    h1 % branch(1) = known_branch(c1(1))

    ! h2 = (a, b), built independently
    c2(1) % branch(1) = known_branch(a); c2(1) % branch(2) = known_branch(c2(2))
    c2(2) % branch(1) = known_branch(b)
    h2 % branch(1) = known_branch(c2(1))

    ! h3 = (b, a), the other order
    c3(1) % branch(1) = known_branch(b); c3(1) % branch(2) = known_branch(c3(2))
    c3(2) % branch(1) = known_branch(a)
    h3 % branch(1) = known_branch(c3(1))

    call check('2  two holders of (a,b) are NOT the same graph', &
         & .not. h1 % same_as(h2))
    call check('2  yet they are the same TUPLE, by component identity', &
         & same_tuple(h1, h2))
    call check('2  and (b,a) is a different tuple: order is part of it', &
         & .not. same_tuple(h1, h3))
    call check('2  so holder identity must never be tuple equality', &
         & .not. h1 % same_as(h2) .and. same_tuple(h1, h2))

  end block tuple_block

  !===================================================================!
  ! 3 . ONE SEMANTIC RELATION, TWO REPRESENTATIONS.
  !
  ! R is a graph. A table and a CSR both describe R's extension. They
  ! answer the same membership questions and they are NOT the same
  ! object - representation identity is storage identity, and it is
  ! not the relation's.
  !===================================================================!

  representation_block: block

    type(graph), target      :: r
    type(relational_binding) :: table_of, csr_of
    type(graph)        :: e, v
    type(stored_relation)    :: t
    type(csr_relation)       :: c
    class(relation), pointer :: pt, pc
    integer                  :: pairs(2, 3)
    logical                  :: agree
    integer                  :: k
    type(set_map)       :: sets
    type(label_map)     :: labels

    call r % declare()

    call e % declare()
    call sets % bind(e, counted_set_representation(3))
    call labels % bind(e, 'E')
    call v % declare()
    call sets % bind(v, counted_set_representation(4))
    call labels % bind(v, 'V')
    pairs = reshape([1, 2,  2, 3,  3, 1], [2, 3])

    t = stored_relation('R', [e, v], pairs, sets)
    c = csr_relation('R', e, v, pairs, sets)

    call table_of % bind_relation(r, t)      ! one binding per storage
    call csr_of   % bind_relation(r, c)

    pt => table_of % relation_for(r)
    pc => csr_of   % relation_for(r)

    agree = pt % num_tuples() .eq. pc % num_tuples()
    do k = 1, 3
       agree = agree .and. (pt % has(pairs(:, k)) .eqv. pc % has(pairs(:, k)))
    end do
    agree = agree .and. (pt % has([1, 4]) .eqv. pc % has([1, 4]))

    call check('3  table and CSR answer the same extension', agree)
    call check('3  and are NOT the same object', .not. pt % same_as(pc))
    call check('3  while R is one relation, whichever storage answers', &
         & r % same_as(r))
    call check('3  the extension is integers, never one graph per tuple', &
         & pt % num_tuples() .eq. 3)

  end block representation_block

  !===================================================================!
  ! 4 . TRANSPOSE, AND WHAT ITS INVOLUTION PRESERVES.
  !
  ! Measured on the production view first: transpose_of mints a fresh
  ! identity, so T(T(R)) is a THIRD object and the law holds only
  ! extensionally. Then the role-permutation prototype, where T(T(R))
  ! is the base itself and the law holds at graph identity.
  !===================================================================!

  transpose_block: block

    type(graph), target        :: r
    type(graph)          :: e, v
    type(csr_relation), target :: c
    type(transposed_relation), target :: ct
    type(transposed_relation)      :: ctt
    type(role_view)            :: rv, rvt, rvtt
    integer                    :: pairs(2, 3)
    integer, allocatable       :: back(:,:)
    logical                    :: same_extension
    integer                    :: k
    type(set_map)       :: sets
    type(label_map)     :: labels

    call r % declare()
    call e % declare()
    call sets % bind(e, counted_set_representation(3))
    call labels % bind(e, 'E')
    call v % declare()
    call sets % bind(v, counted_set_representation(4))
    call labels % bind(v, 'V')
    pairs = reshape([1, 2,  2, 3,  3, 1], [2, 3])
    c = csr_relation('R', e, v, pairs, sets)

    ct  = transpose_of(c)
    ctt = transpose_of(ct)

    call ctt % tuples(back)
    same_extension = size(back, 2) .eq. 3
    do k = 1, 3
       same_extension = same_extension .and. c % has(back(:, k))
    end do

    call check('4  production: T(T(R)) has R''s extension', same_extension)
    call check('4  but is a THIRD identity, not R', &
         & .not. ctt % same_as(c) .and. .not. ct % same_as(c))

    ! The role permutation: nothing is constructed, the roles swap.
    rv % base => r
    rvt  = transpose_role(rv)
    rvtt = transpose_role(rvt)

    call check('4  role view: T swaps the slots', &
         & apply_role(rvt, 1) .eq. 2 .and. apply_role(rvt, 2) .eq. 1)
    call check('4  T(T(.)) is the identity permutation', &
         & apply_role(rvtt, 1) .eq. 1 .and. apply_role(rvtt, 2) .eq. 2)
    call check('4  and denotes the SAME semantic relation R', &
         & rvtt % base % same_as(r) .and. rvt % base % same_as(r))

  end block transpose_block

  !===================================================================!

  if (failures .eq. 0) then
     print *, ''
     print *, ' ALL PROPOSITIONS HOLD'
  else
     print *, ''
     print *, ' FAILURES :', failures
     error stop 'relation: a proposition failed'
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

end program relation_map
