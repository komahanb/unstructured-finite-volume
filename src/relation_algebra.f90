!=====================================================================!
! LEVEL 2 OF THE TOWER. THE RELATION ALGEBRA
!
! The level defines one operation family: HOW RELATIONS GENERATE
! RELATIONS. It contains exactly the three primitives its first
! caller - the calculator tower's dependency derivation - requires,
! and no more (AGENTS.md 9, CALCULATOR.md 9):
!
!      restrict_slot     P|_S       retain the tuples whose i-th part
!                                   is a member of the subobject S;
!                                   the signature is unchanged
!
!      project_slots     pi(P)      retain the selected parts, in the
!                                   selected order; the signature is
!                                   exactly the selection
!
!      compose_binary    T o P      (a, c) wherever some b has
!                                   (a, b) in P and (b, c) in T
!
! No natural join, no union, no intersection, no identity relation,
! no general permutation: each is added when a caller requires it.
! The calculator's derivation reads
!
!      T_flow restricted to the output port, projected to O x X,
!      composed with
!      T_flow restricted to the input ports, projected to X x O,
!
! and returns the one dependency (+, x) - the standard
! join-then-project, factored through a smaller algebra.
!
!                      SEMANTICS, THEN COST
!
! Every result here is MATERIALIZED: restriction and projection
! return stored_relations, composition returns the existing
! csr_relation - no second binary storage, no deferred view hierarchy
! added before a caller requires it.
!
! The costs, parametrically - the semantic pass is only the first
! part of each cost, because materialization inherits the
! constructors' own validation and duplicate collapse:
!
!      restrict    tuple filtering O(|R| * T_has(allowed)), then
!                  stored_relation materialization: carrier
!                  membership validation per slot, and its current
!                  QUADRATIC upper bound, duplicate collapse
!
!      project     slot extraction O(m |R|), then the same generic
!                  stored_relation materialization costs - and here
!                  the collapse has effect, since projection can
!                  make tuples equal
!
!      compose     witness search O(|R| |S|), then csr_relation
!                  materialization, whose cost depends on the
!                  carriers' local_index complexity and the number
!                  of witness pairs produced before collapse
!
! None of this is on a hot path, and none of it is optimized before
! a caller requires it; when a caller composes large relations, an
! indexed composition can be added beside this one.
!
! Set semantics are enforced by the constructors: a projection that
! collapses many tuples to one, or a composition reached through two
! witnesses, stores the tuple once - relations are sets everywhere.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module relation_algebra

  use graph_fractal      , only : graph
  use relation_finitary       , only : relation, stored_relation
  use map_set        , only : set_map
  use map_inclusion  , only : inclusion_map, declared_subobject
  use relation_binary, only : group_by_key, csr_relation

  implicit none

  private
  public :: restrict_slot, project_slots, compose_binary

contains

  !===================================================================!
  ! Restriction: R|_S at one slot. The allowed domain must EMBED in
  ! the slot's own domain - is_subobject_of, never a cardinality
  ! coincidence - so restricting by the full domain itself is the
  ! valid identity restriction. Rejected inputs:
  !
  !      a slot index outside the signature
  !      an allowed domain that does not embed in the slot
  !===================================================================!

  type(stored_relation) function restrict_slot(r, slot_index, allowed, &
       & sets, inclusions) result(restricted)

    class(relation)    , intent(in) :: r
    integer            , intent(in) :: slot_index
    type(graph)    , intent(in) :: allowed
    type(set_map)      , intent(in) :: sets
    type(inclusion_map), intent(in) :: inclusions

    type(graph), allocatable :: domains(:)
    type(graph)              :: d
    integer, allocatable         :: table(:,:), restricted_tuples(:,:)
    integer                      :: k, j, n

    if (slot_index < 1 .or. slot_index > r % arity()) then
       error stop 'relation_algebra: a slot index must name a slot of the relation'
    end if

    !----------------------------------------------------------------!
    ! The embedding is a DECLARED predicate, so the inclusion map
    ! evaluates it; the membership below is an extensional predicate,
    ! so the set map evaluates it. Two predicates, two maps, and
    ! neither is stored in the graph.
    !----------------------------------------------------------------!

    d = r % domain(slot_index)
    if (.not. declared_subobject(allowed, d, inclusions)) then
       error stop 'relation_algebra: a restriction domain must embed in the slot it restricts'
    end if

    allocate(domains(r % arity()))
    do k = 1, r % arity()
       domains(k) = r % domain(k)
    end do

    call r % tuples(table)
    allocate(restricted_tuples(size(table, 1), size(table, 2)))
    n = 0
    do j = 1, size(table, 2)
       if (sets % has(allowed, table(slot_index, j))) then
          n = n + 1
          restricted_tuples(:, n) = table(:, j)
       end if
    end do

    restricted = stored_relation(r % name() // ' restricted', &
         &                     domains, restricted_tuples(:, 1:n), sets)

  end function restrict_slot

  !===================================================================!
  ! Projection: pi(R) onto the chosen slots, in the chosen order -
  ! the order is structural, so [2,1] returns the reversed
  ! signature. Duplicates created by discarding slots collapse in the
  ! constructor: the image is a set. Rejected inputs:
  !
  !      no slots selected
  !      a slot index outside the signature
  !      one slot selected twice - repeated-slot projection would be
  !      a different operation, and is not interpreted implicitly
  !===================================================================!

  type(stored_relation) function project_slots(r, slot_indices, sets) &
       & result(image)

    class(relation), intent(in) :: r
    integer        , intent(in) :: slot_indices(:)
    type(set_map)  , intent(in) :: sets

    type(graph), allocatable :: domains(:)
    integer, allocatable         :: table(:,:), proj(:,:)
    integer                      :: k, l, j, m

    m = size(slot_indices)

    if (m < 1) then
       error stop 'relation_algebra: a projection selects at least one slot'
    end if
    do k = 1, m
       if (slot_indices(k) < 1 .or. slot_indices(k) > r % arity()) then
          error stop 'relation_algebra: a slot index must name a slot of the relation'
       end if
       do l = 1, k - 1
          if (slot_indices(l) == slot_indices(k)) then
             error stop 'relation_algebra: a projection selects each slot at most once'
          end if
       end do
    end do

    allocate(domains(m))
    do k = 1, m
       domains(k) = r % domain(slot_indices(k))
    end do

    call r % tuples(table)
    allocate(proj(m, size(table, 2)))
    do j = 1, size(table, 2)
       do k = 1, m
          proj(k, j) = table(slot_indices(k), j)
       end do
    end do

    image = stored_relation(r % name() // ' projected', &
         &                  domains, proj, sets)

  end function project_slots

  !===================================================================!
  ! Binary composition, argument order and formula fixed together:
  !
  !      compose_binary(P_AB, P_BC)  =  P_BC o P_AB
  !          =  { (a, c) : exists b, (a,b) in P_AB and (b,c) in P_BC }
  !
  ! The middle domains must be the SAME declared domain -
  ! structural identity, never a size coincidence. The result is
  ! binary, and is stored in the existing binary representation,
  ! csr_relation. Rejected inputs:
  !
  !      an argument that is not binary
  !      middle domains that are not one domain
  !===================================================================!

  type(csr_relation) function compose_binary(r_ab, r_bc, sets) result(chained)

    class(relation), intent(in) :: r_ab
    class(relation), intent(in) :: r_bc
    type(set_map)  , intent(in) :: sets

    type(graph)      :: da, db, db2, dc
    integer, allocatable :: tab(:,:), tbc(:,:), pairs(:,:)
    integer              :: i, j, n

    if (r_ab % arity() /= 2 .or. r_bc % arity() /= 2) then
       error stop 'relation_algebra: composition takes two binary relations'
    end if

    db  = r_ab % domain(2)
    db2 = r_bc % domain(1)
    if (.not. db % same_as(db2)) then
       error stop 'relation_algebra: composition requires one shared middle domain'
    end if

    call r_ab % tuples(tab)
    call r_bc % tuples(tbc)

    ! The sparse product: group the right factor's tuples by the
    ! local index of their first slot once, then traverse each left
    ! tuple's matching fibre - linear in the tuples plus the
    ! output, not the all-pairs scan. Output order is unchanged:
    ! left tuples in order, matches in the right factor's tuple
    ! order.
    sparse_product : block

      integer, allocatable :: keys(:), identity(:), ptr(:), grouped(:)
      integer :: nmid, k, m

      nmid = sets % num_members_of(db)

      allocate(keys(size(tbc, 2)), identity(size(tbc, 2)))
      do j = 1, size(tbc, 2)
         keys(j)     = sets % index_in(db, tbc(1, j))
         identity(j) = j
      end do
      call group_by_key(nmid, keys, identity, ptr, grouped)

      n = 0
      do i = 1, size(tab, 2)
         m = sets % index_in(db, tab(2, i))
         n = n + ptr(m + 1) - ptr(m)
      end do
      allocate(pairs(2, n))

      n = 0
      do i = 1, size(tab, 2)
         m = sets % index_in(db, tab(2, i))
         do k = ptr(m), ptr(m + 1) - 1
            n = n + 1
            pairs(:, n) = [tab(1, i), tbc(2, grouped(k))]
         end do
      end do

    end block sparse_product

    da = r_ab % domain(1)
    dc = r_bc % domain(2)

    chained = csr_relation(r_ab % name() // ' then ' // r_bc % name(), &
         &                 da, dc, pairs(:, 1:n), sets)

  end function compose_binary

end module relation_algebra
