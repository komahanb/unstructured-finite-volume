!=====================================================================!
! PARTITIONED IMPLICIT PDE TOWER . LEVEL 4 . GRAPH CALCULUS
!
! The level answers one question: WHAT HAPPENS WHEN THE WHOLE GRAPH
! IS INTERPRETED THROUGH THE PARTITION TRANSFORM. This is where the
! production partitioner first appears, and where the words OWNED,
! BORROWED and OVERLAP first mean something:
!
!      GLOBAL   1 -- 2 -- 3 -- 4 -- 5 -- 6
!                          |
!                         cut
!
!      G1  1 -- 2 -- 3 -- (4)        globals [1,2,3,4], borrows 4
!      G2 (3) -- 4 -- 5 -- 6         globals [4,5,6,3], borrows 3
!
! The parentheses are VISIBILITY, not ownership: a borrowed vertex
! is present because a stencil will need to see it, and it belongs
! to the other part.
!
! And the rung's sharpest truth: production's edge ownership is
! read as a CHOICE BETWEEN TWO NAMED POLICIES, not as agreement with
! a single forced law. Level 2 derived both, from the primitive
! facts alone, before any partitioner existed:
!
!      TailOwner = Own^T o Tail : E -> K
!      HeadOwner = Own^T o Head : E -> K
!
! Both are total functions, so both satisfy the reconstruction law
! and neither is singled out by it. They differ on exactly one edge
! here - the crossing edge e3 - and that one edge is what lets this
! level say which policy production implements:
!
!      production(e3) = part1 = TailOwner(e3) /= HeadOwner(e3)
!
! so the conclusion is empirical and specific: PRODUCTION SELECTS
! THE TAIL-BASED POLICY. Not "production satisfies the uniquely
! forced ownership theorem" - there is no such theorem.
!
! This also names an old prose defect exactly. The partitioner's
! comment once described a hybrid - tail owner normally, head owner
! when the tail is borrowed - while the executable behaviour was
! uniformly TailOwner. With both policies in hand that sentence is
! no longer merely wrong; it is wrong in a way this level can spell.
!
! This level is graph-to-graph interpretation ONLY. No field is
! transported here; that is Level 5's mathematics.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program partitioned_pde_level_4

  use partitioned_pde_assert , only : report, verdict
  use graph_fractal        , only : set_graph => graph
  use map_set_representation, only : counted_set_representation
  use map_set        , only : set_map
  use map_inclusion  , only : inclusion_map, declared_subobject
  use map_label      , only : label_map
  use view_directed    , only : directed_graph
  use relation_binary  , only : csr_relation
  use view_directed_stored            , only : stored_directed_graph
  use transform_partitioner, only : partitioner, PARTITION_LINEAR
  use chain_carriers_fixture , only : chain_carriers
  use chain_relations_fixture, only : tail_relation, head_relation, &
       &                              own_relation
  use chain_algebra_fixture  , only : derive_tail_owner, derive_head_owner

  use relation_partition, only : partition_relation
  implicit none
  type(partition_relation) :: rel

  type(set_graph)          :: v, e, k
  type(set_map)          :: sets
  type(csr_relation), target :: tail, head, own
  type(csr_relation)         :: tail_owner, head_owner
  type(stored_directed_graph)         :: g
  class(directed_graph), allocatable  :: g1, g2
  integer                    :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "partitioned pde tower . level 4 . partition"
  write(*,'(1x,a)') "============================================="

  ! BOTH Level-2 candidate policies, re-derived here so this level
  ! selects between mathematics rather than recalling an observation.
  call chain_carriers(sets, v, e, k)
  tail = tail_relation(e, v, sets)
  head = head_relation(e, v, sets)
  own  = own_relation(k, v, sets)
  tail_owner = derive_tail_owner(tail, own, sets)
  head_owner = derive_head_owner(head, own, sets)

  g = stored_directed_graph(6, tails=[1,2,3,4,5], heads=[2,3,4,5,6])
  call sets % bind(g % vertex_set(), &
       & counted_set_representation(g % num_vertices()))
  call sets % bind(g % edge_set(), &
       & counted_set_representation(g % num_edges()))
  call cut(g1, 1)
  call cut(g2, 2)

  call check_part_identity(g1, 1, nfail)
  call check_part_identity(g2, 2, nfail)
  call check_maps_and_ownership(g1, 1, [1,2,3,4], 4, nfail)
  call check_maps_and_ownership(g2, 2, [4,5,6,3], 3, nfail)
  call check_crossing_edge_presence(nfail)
  call check_production_selects_tail_ownership(nfail)

  call verdict(nfail, "level 4")

contains

  subroutine cut(part, kpart)

    class(directed_graph), allocatable, intent(out) :: part
    integer                  , intent(in)  :: kpart

    type(partitioner) :: p

    p = partitioner(PARTITION_LINEAR, nparts=2, part=kpart)
    call p % partition_graph(g, part, rel)

  end subroutine cut

  subroutine check_part_identity(part, kpart, nfail)

    class(directed_graph), intent(in)    :: part
    integer     , intent(in)    :: kpart
    integer     , intent(inout) :: nfail

    character(len=1) :: tag
    type(partition_relation) :: rel

    write(tag,'(i1)') kpart

    select type (part)
    type is (stored_directed_graph)
       rel = part % whole_relation()
       call report(rel % has_part_relation() .and. &
            &      rel % num_parts() .eq. 2 .and. &
            &      part % id() .eq. kpart, &
            & "G" // tag // " knows it is part " // tag // " of two", &
            & nfail)
    end select

  end subroutine check_part_identity

  !===================================================================!
  ! The maps, and the distinction the whole tower rests on: a part
  ! HOLDS members it does not OWN.
  !===================================================================!

  subroutine check_maps_and_ownership(part, kpart, globals, &
       & borrowed_global, nfail)

    class(directed_graph), intent(in)    :: part
    integer     , intent(in)    :: kpart, globals(:), borrowed_global
    integer     , intent(inout) :: nfail

    type(set_graph) :: owned, borrowed, overlap
    type(label_map)     :: labels
    type(inclusion_map)     :: inclusions
    character(len=1) :: tag
    integer          :: i
    logical          :: ok
    type(partition_relation) :: rel

    write(tag,'(i1)') kpart

    select type (part)
    type is (stored_directed_graph)
       rel = part % whole_relation()

       ok = part % num_vertices() .eq. size(globals)
       do i = 1, min(part % num_vertices(), size(globals))
          ok = ok .and. (rel % global_vertex_index(i) .eq. globals(i))
       end do
       call report(ok, &
            & "G" // tag // "'s local-to-global map is exactly as " // &
            & "declared - local order is not global order", nfail)

       ok = .true.
       do i = 1, part % num_vertices()
          if (rel % global_vertex_index(i) .eq. borrowed_global) then
             ok = ok .and. (rel % vertex_owner_part(i) .ne. kpart)
          else
             ok = ok .and. (rel % vertex_owner_part(i) .eq. kpart)
          end if
       end do
       call report(ok, &
            & "G" // tag // " owns every local vertex but the one it " // &
            & "borrows: PRESENCE IS NOT OWNERSHIP", nfail)

       call part % owned_vertices(kpart, sets, labels, inclusions, owned)
       call part % borrowed_vertices(kpart, sets, labels, inclusions, borrowed)
       call part % overlap_vertices(kpart, sets, labels, inclusions, overlap)
       call report(sets % size_of(owned) .eq. 3 .and. sets % size_of(borrowed) .eq. 1 &
            & .and. sets % size_of(overlap) .eq. part % num_vertices(), &
            & "G" // tag // ": three owned, one borrowed, and the " // &
            & "overlap is the whole local carrier", nfail)

       ok = .false.
       do i = 1, part % num_vertices()
          if (sets % has_in(borrowed, i)) then
             ok = rel % global_vertex_index(i) .eq. borrowed_global
          end if
       end do
       call report(ok, &
            & "G" // tag // " borrows global vertex " // &
            & achar(48 + borrowed_global) // ", and only that one", &
            & nfail)

    end select

  end subroutine check_maps_and_ownership

  subroutine check_crossing_edge_presence(nfail)

    integer, intent(inout) :: nfail

    call report(holds_global_edge(g1, 3) .and. holds_global_edge(g2, 3), &
         & "the crossing edge e3 = 3->4 is PRESENT in both parts", &
         & nfail)
    call report(.not. holds_global_edge(g2, 1) .and. &
         &      .not. holds_global_edge(g1, 5), &
         & "while e1 and e5 each live in one part only", nfail)

  end subroutine check_crossing_edge_presence

  !===================================================================!
  ! THE rung's selection: production's ownership read against BOTH
  ! Level-2 candidate policies, edge for edge. Agreement with
  ! TailOwner alone is what identifies the policy; agreement with
  ! "some total function" would identify nothing.
  !===================================================================!

  subroutine check_production_selects_tail_ownership(nfail)

    integer, intent(inout) :: nfail

    integer :: ge, produced
    logical :: ok

    ok = .true.
    do ge = 1, 5
       produced = owner_of_global_edge(g1, ge)
       if (produced .eq. 0) produced = owner_of_global_edge(g2, ge)
       ok = ok .and. (produced .ne. 0)
       ok = ok .and. tail_owner % has([ge, produced])
    end do
    call report(ok, &
         & "every edge_owner_part production reports agrees with " // &
         & "TailOwner = Own^T o Tail, edge for edge", nfail)

    ! The crossing edge is the whole discriminator: it is the only
    ! place the two policies differ, so it is the only place
    ! production can be caught choosing.
    call report(owner_of_global_edge(g1, 3) .eq. 1 .and. &
         &      owner_of_global_edge(g2, 3) .eq. 1 .and. &
         &      tail_owner % has([3, 1]), &
         & "both parts report the crossing edge as part1's, and so " // &
         & "does TailOwner(e3)", nfail)

    call report(head_owner % has([3, 2]) .and. &
         &      .not. head_owner % has([3, 1]) .and. &
         &      owner_of_global_edge(g1, 3) .ne. 2 .and. &
         &      owner_of_global_edge(g2, 3) .ne. 2, &
         & "while HeadOwner(e3) = part2, which production does NOT " // &
         & "report: PRODUCTION SELECTS THE TAIL-BASED POLICY", nfail)

  end subroutine check_production_selects_tail_ownership

  logical function holds_global_edge(part, ge)

    class(directed_graph), intent(in) :: part
    integer     , intent(in) :: ge

    integer :: i
    type(partition_relation) :: rel

    holds_global_edge = .false.
    select type (part)
    type is (stored_directed_graph)
       rel = part % whole_relation()
       do i = 1, part % num_edges()
          if (rel % global_edge_index(i) .eq. ge) holds_global_edge = .true.
       end do
    end select

  end function holds_global_edge

  integer function owner_of_global_edge(part, ge)

    class(directed_graph), intent(in) :: part
    integer     , intent(in) :: ge

    integer :: i
    type(partition_relation) :: rel

    owner_of_global_edge = 0
    select type (part)
    type is (stored_directed_graph)
       rel = part % whole_relation()
       do i = 1, part % num_edges()
          if (rel % global_edge_index(i) .eq. ge) then
             owner_of_global_edge = rel % edge_owner_part(i)
          end if
       end do
    end select

  end function owner_of_global_edge

end program partitioned_pde_level_4
