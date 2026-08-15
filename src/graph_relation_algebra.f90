!=====================================================================!
! LEVEL 2 OF THE NEW TOWER . THE RELATION ALGEBRA
!
! The level answers one question: HOW RELATIONS GENERATE RELATIONS.
! It holds exactly the three primitives its first real caller - the
! calculator tower's dependency derivation - has earned, and not
! one more (AGENTS.md 9, CALCULATOR.md 9):
!
!      restrict_slot     P|_S       keep the tuples whose i-th part
!                                   the subobject S admits; the
!                                   signature stands unchanged
!
!      project_slots     pi(P)      keep the chosen parts, in the
!                                   chosen order; the signature is
!                                   exactly the selection
!
!      compose_binary    T o P      (a, c) wherever some b carries
!                                   (a, b) in P and (b, c) in T
!
! No natural join, no union, no intersection, no identity relation,
! no general permutation: each waits for the caller that earns it.
! The calculator's derivation reads
!
!      T_flow restricted to the output port, projected to O x X,
!      composed with
!      T_flow restricted to the input ports, projected to X x O,
!
! and answers the one dependency (+, x) - the join-then-project of
! the textbook, factored through a smaller algebra.
!
!                      SEMANTICS, THEN SPEED
!
! Every result here is MATERIALIZED, and says so: restriction and
! projection answer stored_relations, composition answers the
! established csr_relation - no second binary storage, no lazy view
! hierarchy built on speculation.
!
! The costs, parametrically and honestly - the semantic pass is only
! the first half of each bill, because materialization inherits the
! constructors' own validation and collapse:
!
!      restrict    tuple filtering O(|R| * T_has(allowed)), then
!                  stored_relation materialization: carrier
!                  membership validation per slot, and its current
!                  QUADRATIC worst-case duplicate collapse
!
!      project     slot extraction O(m |R|), then the same generic
!                  stored_relation materialization costs - and here
!                  the collapse genuinely works, since projection
!                  makes tuples indistinct
!
!      compose     witness search O(|R| |S|), then csr_relation
!                  materialization, whose cost rides the carriers'
!                  local_index complexity and the number of witness
!                  pairs produced before collapse
!
! None of this stands on a hot path, and none of it is optimized
! ahead of the caller that would earn it; the day a large caller
! composes large relations, an indexed composition can be earned
! beside this one.
!
! Set semantics ride on the constructors: a projection that
! collapses many tuples to one, or a composition reached along two
! witnesses, holds the tuple once - relations are sets everywhere.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module graph_relation_algebra

  use graph_carrier        , only : member_set
  use graph_relation       , only : relation, stored_relation, slot
  use graph_binary_relation, only : csr_relation

  implicit none

  private
  public :: restrict_slot, project_slots, compose_binary

contains

  !===================================================================!
  ! Restriction: R|_S at one slot. The allowed domain must EMBED in
  ! the slot's own domain - is_subobject_of, never a cardinality
  ! coincidence - so restricting by the full domain itself is the
  ! lawful identity restriction. Refusals:
  !
  !      a slot index outside the signature
  !      an allowed domain that does not embed in the slot
  !===================================================================!

  type(stored_relation) function restrict_slot(r, slot_index, allowed) &
       & result(narrowed)

    class(relation)  , intent(in) :: r
    integer          , intent(in) :: slot_index
    class(member_set), intent(in) :: allowed

    type(slot), allocatable        :: seats(:)
    class(member_set), allocatable :: d
    integer, allocatable           :: table(:,:), kept(:,:)
    integer                        :: k, j, n

    if (slot_index < 1 .or. slot_index > r % arity()) then
       error stop 'graph_relation_algebra: a slot index must name a slot of the relation'
    end if

    d = r % domain(slot_index)
    if (.not. allowed % is_subobject_of(d)) then
       error stop 'graph_relation_algebra: a restriction domain must embed in the slot it restricts'
    end if

    allocate(seats(r % arity()))
    do k = 1, r % arity()
       d = r % domain(k)
       seats(k) = slot(d)
    end do

    call r % tuples(table)
    allocate(kept(size(table, 1), size(table, 2)))
    n = 0
    do j = 1, size(table, 2)
       if (allowed % has(table(slot_index, j))) then
          n = n + 1
          kept(:, n) = table(:, j)
       end if
    end do

    narrowed = stored_relation(r % name() // ' restricted', &
         &                     seats, kept(:, 1:n))

  end function restrict_slot

  !===================================================================!
  ! Projection: pi(R) onto the chosen slots, in the chosen order -
  ! the order is structural, so [2,1] answers the reversed
  ! signature. Duplicates born of forgetting slots collapse in the
  ! constructor: the image is a set. Refusals:
  !
  !      no slots chosen
  !      a slot index outside the signature
  !      one slot chosen twice - repeated-slot projection would be
  !      a different operation, and is not quietly interpreted
  !===================================================================!

  type(stored_relation) function project_slots(r, slot_indices) &
       & result(image)

    class(relation), intent(in) :: r
    integer        , intent(in) :: slot_indices(:)

    type(slot), allocatable        :: seats(:)
    class(member_set), allocatable :: d
    integer, allocatable           :: table(:,:), proj(:,:)
    integer                        :: k, l, j, m

    m = size(slot_indices)

    if (m < 1) then
       error stop 'graph_relation_algebra: a projection selects at least one slot'
    end if
    do k = 1, m
       if (slot_indices(k) < 1 .or. slot_indices(k) > r % arity()) then
          error stop 'graph_relation_algebra: a slot index must name a slot of the relation'
       end if
       do l = 1, k - 1
          if (slot_indices(l) == slot_indices(k)) then
             error stop 'graph_relation_algebra: a projection selects each slot at most once'
          end if
       end do
    end do

    allocate(seats(m))
    do k = 1, m
       d = r % domain(slot_indices(k))
       seats(k) = slot(d)
    end do

    call r % tuples(table)
    allocate(proj(m, size(table, 2)))
    do j = 1, size(table, 2)
       do k = 1, m
          proj(k, j) = table(slot_indices(k), j)
       end do
    end do

    image = stored_relation(r % name() // ' projected', &
         &                  seats, proj)

  end function project_slots

  !===================================================================!
  ! Binary composition, argument order and formula locked together:
  !
  !      compose_binary(P_AB, P_BC)  =  P_BC o P_AB
  !          =  { (a, c) : exists b, (a,b) in P_AB and (b,c) in P_BC }
  !
  ! The middle domains must be the SAME declared domain -
  ! structural identity, never a size coincidence. The result is
  ! binary, and lands in the established binary citizen. Refusals:
  !
  !      an argument that is not binary
  !      middle domains that are not one domain
  !===================================================================!

  type(csr_relation) function compose_binary(r_ab, r_bc) result(chained)

    class(relation), intent(in) :: r_ab
    class(relation), intent(in) :: r_bc

    class(member_set), allocatable :: da, db, db2, dc
    integer, allocatable           :: tab(:,:), tbc(:,:), pairs(:,:)
    integer                        :: i, j, n

    if (r_ab % arity() /= 2 .or. r_bc % arity() /= 2) then
       error stop 'graph_relation_algebra: composition takes two binary relations'
    end if

    db  = r_ab % domain(2)
    db2 = r_bc % domain(1)
    if (.not. db % same_as(db2)) then
       error stop 'graph_relation_algebra: composition requires one shared middle domain'
    end if

    call r_ab % tuples(tab)
    call r_bc % tuples(tbc)

    allocate(pairs(2, size(tab, 2) * size(tbc, 2)))
    n = 0
    do i = 1, size(tab, 2)
       do j = 1, size(tbc, 2)
          if (tab(2, i) == tbc(1, j)) then
             n = n + 1
             pairs(:, n) = [tab(1, i), tbc(2, j)]
          end if
       end do
    end do

    da = r_ab % domain(1)
    dc = r_bc % domain(2)

    chained = csr_relation(r_ab % name() // ' then ' // r_bc % name(), &
         &                 da, dc, pairs(:, 1:n))

  end function compose_binary

end module graph_relation_algebra
