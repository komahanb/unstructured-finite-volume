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
! checked against the LEVEL-2 RELATIONAL LAW, not against a
! previous observation of production. Level 2 derived
!
!      EdgeOwner = Own^T o Tail : E -> K
!
! from the primitive facts alone, before any partitioner existed.
! Here every edge_owner_part the machinery reports must agree with
! it - which turns tail-ownership from a convention someone noticed
! into a theorem production happens to satisfy.
!
! This level is graph-to-graph interpretation ONLY. No field is
! transported here; that is Level 5's mathematics.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program partitioned_pde_level_4

  use partitioned_pde_assert , only : report, verdict
  use graph_carrier          , only : counted_set, member_set
  use graph_grammar          , only : graph
  use graph_binary_relation  , only : csr_relation
  use class_graph            , only : stored_graph
  use class_graph_partitioner, only : partitioner, PARTITION_LINEAR
  use chain_relations_fixture, only : chain_carriers, tail_relation, &
       &                              own_relation
  use chain_algebra_fixture  , only : derive_edge_owner

  implicit none

  type(counted_set)          :: v, e, k
  type(csr_relation), target :: tail, own
  type(csr_relation)         :: edge_owner
  type(stored_graph)         :: g
  class(graph), allocatable  :: g1, g2
  integer                    :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "partitioned pde tower . level 4 . partition"
  write(*,'(1x,a)') "============================================="

  ! The Level-2 law, re-derived here so this level checks production
  ! against mathematics rather than against memory.
  call chain_carriers(v, e, k)
  tail = tail_relation(e, v)
  own  = own_relation(k, v)
  edge_owner = derive_edge_owner(tail, own)

  g = stored_graph(6, tails=[1,2,3,4,5], heads=[2,3,4,5,6])
  call cut(g1, 1)
  call cut(g2, 2)

  call check_part_identity(g1, 1, nfail)
  call check_part_identity(g2, 2, nfail)
  call check_maps_and_ownership(g1, 1, [1,2,3,4], 4, nfail)
  call check_maps_and_ownership(g2, 2, [4,5,6,3], 3, nfail)
  call check_crossing_edge_presence(nfail)
  call check_production_matches_the_law(nfail)

  call verdict(nfail, "level 4")

contains

  subroutine cut(part, kpart)

    class(graph), allocatable, intent(out) :: part
    integer                  , intent(in)  :: kpart

    type(partitioner) :: p

    p = partitioner(PARTITION_LINEAR, nparts=2, part=kpart)
    call p % partition_graph(g, part)

  end subroutine cut

  subroutine check_part_identity(part, kpart, nfail)

    class(graph), intent(in)    :: part
    integer     , intent(in)    :: kpart
    integer     , intent(inout) :: nfail

    character(len=1) :: tag

    write(tag,'(i1)') kpart

    select type (part)
    type is (stored_graph)
       call report(part % has_part_relation() .and. &
            &      part % num_parts() .eq. 2 .and. &
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

    class(graph), intent(in)    :: part
    integer     , intent(in)    :: kpart, globals(:), borrowed_global
    integer     , intent(inout) :: nfail

    class(member_set), allocatable :: owned, borrowed, overlap
    character(len=1) :: tag
    integer          :: i
    logical          :: ok

    write(tag,'(i1)') kpart

    select type (part)
    type is (stored_graph)

       ok = part % num_vertices() .eq. size(globals)
       do i = 1, min(part % num_vertices(), size(globals))
          ok = ok .and. (part % global_vertex_index(i) .eq. globals(i))
       end do
       call report(ok, &
            & "G" // tag // "'s local-to-global map is exactly as " // &
            & "declared - local order is not global order", nfail)

       ok = .true.
       do i = 1, part % num_vertices()
          if (part % global_vertex_index(i) .eq. borrowed_global) then
             ok = ok .and. (part % vertex_owner_part(i) .ne. kpart)
          else
             ok = ok .and. (part % vertex_owner_part(i) .eq. kpart)
          end if
       end do
       call report(ok, &
            & "G" // tag // " owns every local vertex but the one it " // &
            & "borrows: PRESENCE IS NOT OWNERSHIP", nfail)

       call part % owned_vertices(kpart, owned)
       call part % borrowed_vertices(kpart, borrowed)
       call part % overlap_vertices(kpart, overlap)
       call report(owned % size() .eq. 3 .and. borrowed % size() .eq. 1 &
            & .and. overlap % size() .eq. part % num_vertices(), &
            & "G" // tag // ": three owned, one borrowed, and the " // &
            & "overlap is the whole local carrier", nfail)

       ok = .false.
       do i = 1, part % num_vertices()
          if (borrowed % has(i)) then
             ok = part % global_vertex_index(i) .eq. borrowed_global
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
  ! THE rung's theorem-check: production's ownership against the
  ! Level-2 derivation, edge for edge, in both directions.
  !===================================================================!

  subroutine check_production_matches_the_law(nfail)

    integer, intent(inout) :: nfail

    integer :: ge, produced
    logical :: ok

    ok = .true.
    do ge = 1, 5
       produced = owner_of_global_edge(g1, ge)
       if (produced .eq. 0) produced = owner_of_global_edge(g2, ge)
       ok = ok .and. (produced .ne. 0)
       ok = ok .and. edge_owner % has([ge, produced])
    end do
    call report(ok, &
         & "every edge_owner_part production reports satisfies the " // &
         & "Level-2 law EdgeOwner = Own^T o Tail", nfail)

    call report(owner_of_global_edge(g1, 3) .eq. 1 .and. &
         &      owner_of_global_edge(g2, 3) .eq. 1 .and. &
         &      edge_owner % has([3, 1]), &
         & "and both parts agree the crossing edge is part1's - as " // &
         & "the law derived before any partitioner existed", nfail)

  end subroutine check_production_matches_the_law

  logical function holds_global_edge(part, ge)

    class(graph), intent(in) :: part
    integer     , intent(in) :: ge

    integer :: i

    holds_global_edge = .false.
    select type (part)
    type is (stored_graph)
       do i = 1, part % num_edges()
          if (part % global_edge_index(i) .eq. ge) holds_global_edge = .true.
       end do
    end select

  end function holds_global_edge

  integer function owner_of_global_edge(part, ge)

    class(graph), intent(in) :: part
    integer     , intent(in) :: ge

    integer :: i

    owner_of_global_edge = 0
    select type (part)
    type is (stored_graph)
       do i = 1, part % num_edges()
          if (part % global_edge_index(i) .eq. ge) then
             owner_of_global_edge = part % edge_owner_part(i)
          end if
       end do
    end select

  end function owner_of_global_edge

end program partitioned_pde_level_4
