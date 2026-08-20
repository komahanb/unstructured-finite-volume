!=====================================================================!
! VISUALIZATION TOWER . LEVEL 3 . RELATIONAL GRAPH
!
! The level answers one question: CAN THE COMPLETE OPERATOR CHAIN
! EXIST AS ONE TYPED RELATIONAL STRUCTURE.
!
!      GAMMA = ( { X0 X1 X2 X3 E1 E2 E3 },
!                { T1 H1 T2 H2 T3 H3 } )
!
! Seven carriers, six primitive relations, and NOTHING ELSE. What is
! absent from that list is as much a claim as what is present.
!
!                       DERIVED STRUCTURE STAYS DERIVED
!
! D1, D2, D3, D2:1, D3:1 and every transpose are absent from the
! graph. They are true of it, and Level 2 showed how to obtain them,
! but a graph that stored them would hold one mathematical truth in
! two places and would have to keep them agreeing forever. The graph
! owns the twelve occurrences' incidence; everything else is a
! question the algebra answers on demand, and this level proves that
! the answer can still be reached with nothing but graph-owned
! relations in hand.
!
!                     NO UNION CARRIER IS MANUFACTURED
!
! There is no V = X0 u X1 u X2 u X3 here. The temptation to make one
! is real - it is what an ordinary graph would want, and Level 4 will
! examine the request properly - but a union would erase exactly the
! typing that makes D1 a rectangular X0 -> X1 rather than an
! adjacency on twelve anonymous vertices. The graph holds seven
! domains because the mathematics has seven domains.
!
!                        OWNERSHIP AND WHOLENESS
!
! The graph owns MATERIALIZED relations only. A transposed view
! borrows its base by pointer, so owning one would mean owning a
! pointer to something the graph does not keep alive. Views ride
! above graph-owned relations; that is why Level 2 materialized its
! transposes before passing them anywhere.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program visualization_level_3

  use visualization_assert , only : report, verdict
  use visualization_assert , only : ND1, ND2, ND3, ND21, ND31
  use visualization_assert , only : NX0, NX1, NX2, NX3, NE1, NE2, NE3
  use visualization_assert , only : X0_A, X0_B, X0_C, X0_D
  use visualization_assert , only : X1_P, X1_Q, X1_R
  use visualization_assert , only : X2_U, X2_V, X2_W
  use visualization_assert , only : X3_M, X3_N
  use graph_fractal        , only : graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set        , only : set_map
  use map_label      , only : label_map
  use map_inclusion  , only : inclusion_map, declared_subobject
  use relation_finitary       , only : relation
  use relation_binary, only : csr_relation, binary_relation
  use relation_binary, only : transposed_relation, transpose_of
  use graph_fractal        , only : graph, known_branch, null_branch
  use view_relational, only : relational_binding, &
       & num_member_sets, member_set_at, num_relations, relation_at, &
       & has_set
  use visualization_carriers_fixture , only : structural_carriers
  use visualization_relations_fixture, only : occurrences_of_a1
  use visualization_relations_fixture, only : occurrences_of_a2
  use visualization_relations_fixture, only : occurrences_of_a3
  use visualization_algebra_fixture  , only : derive_dependency
  use visualization_algebra_fixture  , only : derive_composition
  use visualization_algebra_fixture  , only : materialized_transpose
  use visualization_algebra_fixture  , only : same_extension

  implicit none

  type(graph)              :: x0, x1, x2, x3, e1, e2, e3
  type(set_map)     :: sets
  type(label_map)     :: labels
  type(csr_relation)     , target :: t1, h1, t2, h2, t3, h3
  type(graph)             , target :: g
  type(graph)             , target :: scell(7), selem(7)
  type(graph)             , target :: rcell(6), relem(6)
  type(relational_binding)         :: bnd
  integer                          :: k
  integer                         :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "visualization tower . level 3 . relational graph"
  write(*,'(1x,a)') "============================================="

  call structural_carriers(x0, x1, x2, x3, e1, e2, e3, sets, labels)
  call occurrences_of_a1(e1, x0, x1, t1, h1, sets)
  call occurrences_of_a2(e2, x1, x2, t2, h2, sets)
  call occurrences_of_a3(e3, x2, x3, t3, h3, sets)

  ! 'the operator chain A3 o A2 o A1': (S, P) as one sequence on each branch.
  call g % declare()
  do k = 1, 7
     call scell(k) % declare()
     call selem(k) % declare()
  end do
  do k = 1, 6
     call rcell(k) % declare()
     call relem(k) % declare()
  end do

  call bnd % bind_set(selem(1), x0)
  call bnd % bind_set(selem(2), x1)
  call bnd % bind_set(selem(3), x2)
  call bnd % bind_set(selem(4), x3)
  call bnd % bind_set(selem(5), e1)
  call bnd % bind_set(selem(6), e2)
  call bnd % bind_set(selem(7), e3)
  call bnd % bind_relation(relem(1), t1)
  call bnd % bind_relation(relem(2), h1)
  call bnd % bind_relation(relem(3), t2)
  call bnd % bind_relation(relem(4), h2)
  call bnd % bind_relation(relem(5), t3)
  call bnd % bind_relation(relem(6), h3)

  do k = 1, 7
     scell(k) % branch(1) = known_branch(selem(k))
     if (k .lt. 7) scell(k) % branch(2) = &
          & known_branch(scell(k + 1))
  end do
  do k = 1, 6
     rcell(k) % branch(1) = known_branch(relem(k))
     if (k .lt. 6) rcell(k) % branch(2) = &
          & known_branch(rcell(k + 1))
  end do

  g % branch(1) = known_branch(scell(1))
  g % branch(2) = known_branch(rcell(1))

  call check_the_seven_carriers_are_owned(nfail)
  call check_the_six_primitives_are_owned(nfail)
  call check_signature_closure(nfail)
  call check_the_carriers_stay_apart(nfail)
  call check_no_union_carrier_exists(nfail)
  call check_nothing_derived_is_stored(nfail)
  call check_the_chain_is_recoverable(nfail)

  call verdict(nfail, "level 3")

contains

  !===================================================================!
  ! Seven, held by structural identity - not by name, not by size,
  ! and not by having been handed to the constructor.
  !===================================================================!

  subroutine check_the_seven_carriers_are_owned(nfail)

    integer, intent(inout) :: nfail

    call report(num_member_sets(g) .eq. 7, &
         & "the graph holds SEVEN member sets - four state carriers " // &
         & "and three occurrence carriers", nfail)

    call report(has_set(g, bnd, x0) .and. has_set(g, bnd, x1) .and. &
         &      has_set(g, bnd, x2) .and. has_set(g, bnd, x3) .and. &
         &      has_set(g, bnd, e1) .and. has_set(g, bnd, e2) .and. &
         &      has_set(g, bnd, e3), &
         & "and each of X0 X1 X2 X3 E1 E2 E3 answers same_as against " // &
         & "one of them: OWNERSHIP IS IDENTITY", nfail)

    call report(sizes_kept(), &
         & "each owned carrier still holds what it declared: " // &
         & "4 3 3 2 5 4 3", nfail)

  end subroutine check_the_seven_carriers_are_owned

  !===================================================================!
  ! Six primitive relations, each found by identity rather than by
  ! the seat it happens to occupy.
  !===================================================================!

  subroutine check_the_six_primitives_are_owned(nfail)

    integer, intent(inout) :: nfail

    call report(num_relations(g) .eq. 6, &
         & "the graph holds SIX relations - the two ends of each of " // &
         & "the three operators' occurrences", nfail)

    call report(seat_of(t1) .gt. 0 .and. seat_of(h1) .gt. 0 .and. &
         &      seat_of(t2) .gt. 0 .and. seat_of(h2) .gt. 0 .and. &
         &      seat_of(t3) .gt. 0 .and. seat_of(h3) .gt. 0, &
         & "T1 H1 T2 H2 T3 H3 are all found by identity - a selector " // &
         & "is an address, never a shape", nfail)

    call report(counted_occurrences() .eq. NE1 + NE2 + NE3, &
         & "and the twelve occurrences survive the admission: 5 + 4 " // &
         & "+ 3 tails owned", nfail)

  end subroutine check_the_six_primitives_are_owned

  !===================================================================!
  ! The signature validity law, checked from outside rather than
  ! trusted: every slot of every owned relation names an owned
  ! carrier.
  !===================================================================!

  subroutine check_signature_closure(nfail)

    integer, intent(inout) :: nfail

    class(relation)  , pointer     :: r
    type(graph) :: d
    integer                        :: k, s
    logical                        :: closed

    closed = .true.
    do k = 1, num_relations(g)
       r => relation_at(g, bnd, k)
       do s = 1, r % arity()
          d = r % domain(s)
          closed = closed .and. has_set(g, bnd, d)
       end do
    end do

    call report(closed, &
         & "EVERY SLOT OF EVERY OWNED RELATION RESOLVES TO AN OWNED " // &
         & "CARRIER - the graph is closed under its own signatures", nfail)

    call report(num_relations(g) .eq. 6 .and. all_binary(), &
         & "all twelve slots come from six binary relations, so the " // &
         & "closure covers twelve slot resolutions", nfail)

  end subroutine check_signature_closure

  !===================================================================!
  ! Admission changed nothing about identity: the four state carriers
  ! are still four, and the three occurrence carriers still three.
  !===================================================================!

  subroutine check_the_carriers_stay_apart(nfail)

    integer, intent(inout) :: nfail

    type(graph), pointer :: a, b
    integer                    :: i, j
    logical                    :: apart

    apart = .true.
    do i = 1, num_member_sets(g)
       a => member_set_at(g, bnd, i)
       do j = 1, num_member_sets(g)
          if (i .eq. j) cycle
          b => member_set_at(g, bnd, j)
          apart = apart .and. (.not. a % same_as(b))
       end do
    end do

    call report(apart, &
         & "no two OWNED carriers are the same domain - X0 X1 X2 X3 " // &
         & "stay distinct, and so do E1 E2 E3", nfail)

    call report(.not. x1 % same_as(x2) .and. has_set(g, bnd, x1) .and. &
         &      has_set(g, bnd, x2), &
         & "the graph holds both three-member state carriers, and " // &
         & "still tells X1 from X2", nfail)

  end subroutine check_the_carriers_stay_apart

  !===================================================================!
  ! NO FAKE VERTEX CARRIER. A twelve-member union of the four state
  ! carriers would be exactly what an ordinary graph asks for, and
  ! the graph does not hold one.
  !===================================================================!

  subroutine check_no_union_carrier_exists(nfail)

    integer, intent(inout) :: nfail

    type(graph) :: union_like

    ! A manufactured carrier: declared and described here so the
    ! test can ask its size, and NOT held by the graph.
    call union_like % declare()
    call sets % bind(union_like, &
         & counted_set_representation(NX0 + NX1 + NX2 + NX3))

    call report(.not. has_set(g, bnd, union_like), &
         & "a manufactured twelve-member vertex carrier IS NOT HELD - " // &
         & "the chain was not collapsed to make it renderable", nfail)

    call report(num_member_sets(g) .eq. 7 .and. sets % num_members_of(union_like) .eq. 12, &
         & "the graph keeps seven typed domains where an ordinary " // &
         & "graph would want one untyped set of twelve", nfail)

  end subroutine check_no_union_carrier_exists

  !===================================================================!
  ! Derived structure stays derived. No owned relation runs between
  ! two state carriers, and a view - the other thing that must never
  ! be owned - fails the wholeness test the constructor applies.
  !===================================================================!

  subroutine check_nothing_derived_is_stored(nfail)

    integer, intent(inout) :: nfail

    type(transposed_relation) :: view
    class(relation), pointer :: r
    integer :: k
    logical :: any_dependency

    any_dependency = .false.
    do k = 1, num_relations(g)
       r => relation_at(g, bnd, k)
       any_dependency = any_dependency .or. runs_between_states(r)
    end do

    call report(.not. any_dependency, &
         & "NO OWNED RELATION RUNS X_(k-1) -> X_k: D1 D2 D3 D2:1 D3:1 " // &
         & "and every transpose are derived, and stay derived", nfail)

    view = transpose_of(t1)
    call report(.not. view % materialized() .and. t1 % materialized(), &
         & "and a transposed view is not whole unto itself, so the " // &
         & "graph would refuse to own one - views ride above", nfail)

  end subroutine check_nothing_derived_is_stored

  !===================================================================!
  ! THE LEVEL'S POINT. Hand back nothing but the graph, and the whole
  ! forward chain still follows - derived from graph-owned primitive
  ! relations reached through relation_at, never from the local
  ! copies this program happens to be holding.
  !===================================================================!

  subroutine check_the_chain_is_recoverable(nfail)

    integer, intent(inout) :: nfail

    class(binary_relation), pointer :: gt1, gh1, gt2, gh2, gt3, gh3
    type(csr_relation), target      :: gd1, gd2, gd3, gd21, gd31
    type(csr_relation)              :: gd31t, grev
    type(csr_relation), target      :: gd1t, gd2t, gd3t
    type(csr_relation)              :: gmid

    gt1 => owned_binary(t1); gh1 => owned_binary(h1)
    gt2 => owned_binary(t2); gh2 => owned_binary(h2)
    gt3 => owned_binary(t3); gh3 => owned_binary(h3)

    gd1 = derive_dependency('D1', gt1, gh1, sets)
    gd2 = derive_dependency('D2', gt2, gh2, sets)
    gd3 = derive_dependency('D3', gt3, gh3, sets)

    gd21 = derive_composition('D2:1', gd1,  gd2, sets)
    gd31 = derive_composition('D3:1', gd21, gd3, sets)

    call report(gd1 % num_tuples() .eq. ND1 .and. &
         &      gd2 % num_tuples() .eq. ND2 .and. &
         &      gd3 % num_tuples() .eq. ND3, &
         & "D1, D2, D3 re-derived from GRAPH-OWNED incidence alone", nfail)

    call report(gd21 % num_tuples() .eq. ND21 .and. &
         &      gd21 % has([X0_B, X2_U]) .and. gd21 % has([X0_C, X2_V]) .and. &
         &      .not. gd21 % has([X0_A, X2_W]), &
         & "D2:1 follows from them, six tuples, the witness collapse " // &
         & "intact", nfail)

    call report(gd31 % num_tuples() .eq. ND31 .and. &
         &      gd31 % has([X0_A, X3_M]) .and. gd31 % has([X0_D, X3_N]) .and. &
         &      .not. gd31 % has([X0_A, X3_N]) .and. &
         &      .not. gd31 % has([X0_D, X3_M]), &
         & "AND SO DOES THE WHOLE SKELETON D3:1 - the graph carries " // &
         & "the chain without storing it", nfail)

    ! The reverse chain too, from the same owned primitives.
    gd1t = materialized_transpose('D1^T', gd1, sets)
    gd2t = materialized_transpose('D2^T', gd2, sets)
    gd3t = materialized_transpose('D3^T', gd3, sets)
    gmid = derive_composition('D2^T o D3^T', gd3t, gd2t, sets)
    grev = derive_composition('Drev', gmid, gd1t, sets)
    gd31t = materialized_transpose('D3:1^T', gd31, sets)

    call report(same_extension(grev, gd31t), &
         & "and the reverse chain closes on the graph's own relations " // &
         & "as well: (D3 o D2 o D1)^T = D1^T o D2^T o D3^T", nfail)

    call report(gd1 % has([X0_A, X1_P]) .and. gd3 % has([X2_W, X3_N]) .and. &
         &      .not. same_extension(gd31, gd31t), &
         & "the recovered chain is the same mathematics Level 2 " // &
         & "derived, orientation and all", nfail)

  end subroutine check_the_chain_is_recoverable

  !===================================================================!
  ! Helpers - all of them asking the graph, never the local copies.
  !===================================================================!

  integer function seat_of(selector)

    class(relation), intent(in) :: selector

    class(relation), pointer :: r
    integer                  :: k

    seat_of = 0
    do k = 1, num_relations(g)
       r => relation_at(g, bnd, k)
       if (r % same_as(selector)) then
          seat_of = k
          return
       end if
    end do

  end function seat_of

  function owned_binary(selector) result(rp)

    class(relation), intent(in)     :: selector
    class(binary_relation), pointer :: rp

    class(relation), pointer :: r

    rp => null()
    r  => relation_at(g, bnd, seat_of(selector))
    select type (r)
    class is (binary_relation)
       rp => r
    class default
       error stop 'level 3: the graph owns a primitive that is not binary'
    end select

  end function owned_binary

  logical function sizes_kept()

    type(graph), pointer :: c
    integer                    :: want(7), k

    want = [NX0, NX1, NX2, NX3, NE1, NE2, NE3]
    sizes_kept = .true.
    do k = 1, 7
       c => member_set_at(g, bnd, k)
       sizes_kept = sizes_kept .and. (sets % num_members_of(c) .eq. want(k))
    end do

  end function sizes_kept

  logical function all_binary()

    class(relation), pointer :: r
    integer                  :: k

    all_binary = .true.
    do k = 1, num_relations(g)
       r => relation_at(g, bnd, k)
       all_binary = all_binary .and. (r % arity() .eq. 2)
    end do

  end function all_binary

  integer function counted_occurrences()

    class(relation), pointer :: r
    integer                  :: k

    counted_occurrences = 0
    do k = 1, num_relations(g)
       r => relation_at(g, bnd, k)
       if (r % name() .eq. 'T1' .or. r % name() .eq. 'T2' .or. &
            & r % name() .eq. 'T3') then
          counted_occurrences = counted_occurrences + r % num_tuples()
       end if
    end do

  end function counted_occurrences

  !-------------------------------------------------------------------!
  ! Does this relation run from one state carrier to the next - which
  ! is to say, is it a stored dependency wearing a primitive's seat.
  !-------------------------------------------------------------------!

  logical function runs_between_states(r)

    class(relation), intent(in) :: r

    type(graph) :: first, second

    runs_between_states = .false.
    if (r % arity() .ne. 2) return

    first  = r % domain(1)
    second = r % domain(2)

    runs_between_states = &
         & (first % same_as(x0) .and. second % same_as(x1)) .or. &
         & (first % same_as(x1) .and. second % same_as(x2)) .or. &
         & (first % same_as(x2) .and. second % same_as(x3)) .or. &
         & (first % same_as(x0) .and. second % same_as(x2)) .or. &
         & (first % same_as(x0) .and. second % same_as(x3))

  end function runs_between_states

end program visualization_level_3
