!=====================================================================!
! PARTITIONED IMPLICIT PDE TOWER . GATE A . PARTITION
!
! The gate answers one question: DOES THE PARTITION MACHINERY
! PRESERVE ENOUGH TRUTH FOR A TOPOLOGY-CONSUMING OPERATION TO RUN
! LOCALLY. No operator and no solver appear here - only structure,
! ownership, visibility, transport and reconstruction.
!
!      GLOBAL   1 -- 2 -- 3 -- 4 -- 5 -- 6
!                          |
!                         cut
!
!      PART 1   1 -- 2 -- 3 -- (4)        (4) borrowed
!      PART 2  (3) -- 4 -- 5 -- 6         (3) borrowed
!
! The parentheses are VISIBILITY, not ownership, and the whole gate
! turns on keeping those apart:
!
!      borrowed INPUTS  are necessary - a stencil at owned vertex 3
!                       cannot be evaluated without seeing q(4)
!      borrowed OUTPUTS are disposable - only owned members may
!                       contribute back to the whole
!
! The crossing edge e3 = 3->4 lives in BOTH parts, and this gate
! does NOT guess which part owns it. It imposes the assembly law
! first -
!
!      one global entity  ->  one assembled contribution
!
! - with an edge probe whose reconstruction must be exact, and only
! then reads edge_owner_part to DOCUMENT the canonical owner. A law
! decides; an inspection records.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program partitioned_pde_gate_a

  use iso_fortran_env  , only : dp => REAL64
  use partitioned_pde_assert, only : report, verdict
  use partitioned_pde_assert, only : NV, NE, Q_EXACT
  use graph_carrier    , only : counted_set, subset_set, member_set
  use graph_grammar    , only : graph, graph_field
  use class_graph      , only : stored_graph
  use class_graph_field, only : field
  use class_graph_partitioner, only : partitioner, PARTITION_LINEAR
  use class_graph_assembler  , only : assembler

  implicit none

  type(stored_graph)        :: g
  type(assembler)           :: a
  class(graph), allocatable :: g1, g2
  integer                   :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "partitioned pde tower . gate A . partition"
  write(*,'(1x,a)') "============================================="

  g = stored_graph(NV, tails=[1,2,3,4,5], heads=[2,3,4,5,6])
  a = assembler()

  call cut_the_graph(g1, 1)
  call cut_the_graph(g2, 2)

  call check_global_topology(nfail)
  call check_part_structure(g1, 1, [1,2,3,4], 4, nfail)
  call check_part_structure(g2, 2, [4,5,6,3], 3, nfail)
  call check_crossing_edge_is_present_twice(nfail)
  call check_edge_assembly_law(nfail)
  call document_canonical_edge_owner(nfail)
  call check_vertex_transport(nfail)
  call check_vertex_round_trip(nfail)
  call check_proper_subset_transport(nfail)

  call verdict(nfail, "gate A")

contains

  !===================================================================!
  ! One part, cut from the global graph by the production rule.
  !===================================================================!

  subroutine cut_the_graph(part, k)

    class(graph), allocatable, intent(out) :: part
    integer                  , intent(in)  :: k

    type(partitioner) :: p

    p = partitioner(PARTITION_LINEAR, nparts=2, part=k)
    call p % partition_graph(g, part)

  end subroutine cut_the_graph

  !===================================================================!
  ! The global chain, exactly - six vertices, five edges, and each
  ! edge running from its stated tail to its stated head.
  !===================================================================!

  subroutine check_global_topology(nfail)

    integer, intent(inout) :: nfail

    type(counted_set) :: vs, es
    integer           :: e
    logical           :: ok

    call report(g % num_vertices() .eq. NV .and. &
         &      g % num_edges() .eq. NE, &
         & "G is the six-vertex chain with five edges", nfail)

    ok = .true.
    do e = 1, NE
       ok = ok .and. (g % edge_tail(e) .eq. e)
       ok = ok .and. g % edge_has_head(e)
       ok = ok .and. (g % edge_head(e) .eq. e + 1)
    end do
    call report(ok, &
         & "every edge e runs from vertex e to vertex e+1", nfail)

    vs = g % vertex_set()
    es = g % edge_set()
    call report(vs % size() .eq. NV .and. es % size() .eq. NE .and. &
         &      .not. vs % same_as(es), &
         & "its vertex and edge carriers are distinct domains", nfail)

    call report(.not. g % has_part_relation(), &
         & "and the whole graph is not itself a part", nfail)

  end subroutine check_global_topology

  !===================================================================!
  ! One part's structure: its identity, its global map, and the
  ! agreement between owner-per-local-member and the owned/borrowed/
  ! overlap subsets. Local member numbers are NOT global numbers.
  !===================================================================!

  subroutine check_part_structure(part, k, globals, borrowed_global, &
       & nfail)

    class(graph)      , intent(in)    :: part
    integer           , intent(in)    :: k
    integer           , intent(in)    :: globals(:)
    integer           , intent(in)    :: borrowed_global
    integer           , intent(inout) :: nfail

    class(member_set), allocatable :: owned, borrowed, overlap
    character(len=1)  :: tag
    integer           :: i
    logical           :: ok

    write(tag,'(i1)') k

    select type (part)
    type is (stored_graph)

       call report(part % has_part_relation() .and. &
            &      part % num_parts() .eq. 2 .and. &
            &      part % id() .eq. k, &
            & "G" // tag // " knows it is part " // tag // " of two", &
            & nfail)

       call report(part % num_vertices() .eq. size(globals), &
            & "G" // tag // " holds four vertices: three owned and " // &
            & "one borrowed", nfail)

       ok = .true.
       do i = 1, min(part % num_vertices(), size(globals))
          ok = ok .and. (part % global_vertex_index(i) .eq. globals(i))
       end do
       call report(ok, &
            & "G" // tag // "'s global map is exactly as declared - " // &
            & "local order is not global order", nfail)

       ! Owner per local member, against the owned/borrowed subsets.
       ok = .true.
       do i = 1, part % num_vertices()
          if (part % global_vertex_index(i) .eq. borrowed_global) then
             ok = ok .and. (part % vertex_owner_part(i) .ne. k)
          else
             ok = ok .and. (part % vertex_owner_part(i) .eq. k)
          end if
       end do
       call report(ok, &
            & "G" // tag // " owns every local vertex but the one it " // &
            & "borrows", nfail)

       call part % owned_vertices(k, owned)
       call part % borrowed_vertices(k, borrowed)
       call part % overlap_vertices(k, overlap)

       call report(owned % size() .eq. 3 .and. &
            &      borrowed % size() .eq. 1 .and. &
            &      overlap % size() .eq. 4, &
            & "G" // tag // ": |owned| = 3, |borrowed| = 1, " // &
            & "|overlap| = 4", nfail)

       ! The borrowed member, named by its GLOBAL index.
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

       ! Overlap is the carrier a stencil may read.
       call report(overlap % size() .eq. part % num_vertices(), &
            & "G" // tag // "'s overlap is its whole local carrier: " // &
            & "owned union borrowed", nfail)

    class default
       call report(.false., "G" // tag // " is a stored graph", nfail)
    end select

  end subroutine check_part_structure

  !===================================================================!
  ! The crossing edge is PRESENT in both parts - presence is not
  ! ownership, and the next two checks are about telling them apart.
  !===================================================================!

  subroutine check_crossing_edge_is_present_twice(nfail)

    integer, intent(inout) :: nfail

    call report(holds_global_edge(g1, 3) .and. holds_global_edge(g2, 3), &
         & "global e3 = 3->4 is present in BOTH parts", nfail)
    call report(.not. holds_global_edge(g2, 1) .and. &
         &      .not. holds_global_edge(g1, 5), &
         & "while e1 and e5 each live in one part only", nfail)

  end subroutine check_crossing_edge_is_present_twice

  !===================================================================!
  ! THE LAW, imposed before any inspection: partition an edge probe
  ! to both parts, assemble each home, and the sum must be the probe
  ! EXACTLY. The crossing edge must contribute 30 once - not twice,
  ! not never. This is what makes an edge owner canonical; the
  ! ownership arrays are consulted afterwards, to document.
  !===================================================================!

  subroutine check_edge_assembly_law(nfail)

    integer, intent(inout) :: nfail

    real(dp), parameter :: PROBE(NE) = &
         & [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 50.0_dp]
    type(field)           :: z
    real(dp), allocatable :: total(:), got(:)

    z = field('edge probe', g % edge_set())
    call z % set_real_vector(PROBE)

    allocate(total(NE))
    total = 0.0_dp
    call add_round_trip(z, 1, total)
    call add_round_trip(z, 2, total)

    call report(maxval(abs(total - PROBE)) < 1.0d-13, &
         & "every global edge is assembled exactly once: the probe " // &
         & "returns unchanged", nfail)
    call report(abs(total(3) - 30.0_dp) < 1.0d-13, &
         & "and the CROSSING edge contributes 30 exactly once - " // &
         & "not doubled, not lost", nfail)

    ! The same law for a full vertex probe.
    deallocate(total)
    allocate(total(NV))
    total = 0.0_dp
    z = field('vertex probe', g % vertex_set())
    call z % set_real_vector(Q_EXACT)
    call add_round_trip(z, 1, total)
    call add_round_trip(z, 2, total)
    call report(maxval(abs(total - Q_EXACT)) < 1.0d-13, &
         & "and every global vertex likewise, borrowed copies " // &
         & "notwithstanding", nfail)

    if (allocated(got)) deallocate(got)

  end subroutine check_edge_assembly_law

  !===================================================================!
  ! Only NOW: which part does the machinery consider canonical for
  ! the crossing edge? This is documentation of a law already
  ! satisfied, never a decision.
  !
  ! The production comment beside the assignment says an edge is
  ! owned by its tail's part "unless the tail is borrowed, in which
  ! case the head's owner answers for it" - but both branches of
  ! that if/else assign the tail's owner, so the code is
  ! unconditionally tail-owned. The law above holds regardless,
  ! because tail-ownership is still a single well-defined owner per
  ! global edge. What is wrong is the PROSE, not the behaviour.
  !===================================================================!

  subroutine document_canonical_edge_owner(nfail)

    integer, intent(inout) :: nfail

    integer :: o1, o2

    o1 = owner_of_global_edge(g1, 3)
    o2 = owner_of_global_edge(g2, 3)

    call report(o1 .eq. o2 .and. o1 .ne. 0, &
         & "both parts agree on who owns the crossing edge", nfail)
    call report(o1 .eq. 1, &
         & "and it is part 1 - the owner of its TAIL, vertex 3", &
         & nfail)

    call report(owner_of_global_edge(g1, 1) .eq. 1 .and. &
         &      owner_of_global_edge(g2, 4) .eq. 2 .and. &
         &      owner_of_global_edge(g2, 5) .eq. 2, &
         & "every other edge is owned by its tail's part too: the " // &
         & "rule is uniform", nfail)

  end subroutine document_canonical_edge_owner

  !===================================================================!
  ! A FULL global vertex field becomes a FULL overlap field on each
  ! part - values read by GLOBAL member, never by position.
  !===================================================================!

  subroutine check_vertex_transport(nfail)

    integer, intent(inout) :: nfail

    type(field) :: q

    q = field('q star', g % vertex_set())
    call q % set_real_vector(Q_EXACT)

    call check_one_transport(q, g1, 1, [1,2,3,4], &
         & [1.0_dp, 2.0_dp, 4.0_dp, 7.0_dp], nfail)
    call check_one_transport(q, g2, 2, [4,5,6,3], &
         & [7.0_dp, 11.0_dp, 16.0_dp, 4.0_dp], nfail)

  end subroutine check_vertex_transport

  subroutine check_one_transport(q, part, k, globals, expect, nfail)

    type(field)      , intent(in)    :: q
    class(graph)     , intent(in)    :: part
    integer          , intent(in)    :: k, globals(:)
    real(dp)         , intent(in)    :: expect(:)
    integer          , intent(inout) :: nfail

    type(partitioner)               :: p
    class(graph_field), allocatable :: pd
    class(member_set), allocatable  :: dom
    type(counted_set)               :: pvs
    real(dp), allocatable           :: v(:)
    character(len=1)                :: tag
    integer                         :: i
    logical                         :: ok

    write(tag,'(i1)') k

    ! The data is transported onto THIS part - not onto a second cut
    ! of the same graph, which would be a different carrier with a
    ! different identity.
    p = partitioner(PARTITION_LINEAR, nparts=2, part=k)
    call p % partition_data(g, q, part, pd)

    call pd % domain(dom)
    select type (part)
    type is (stored_graph)
       pvs = part % vertex_set()
       call report(dom % same_as(pvs) .and. dom % size() .eq. 4, &
            & "q" // tag // " lives on G" // tag // "'s whole vertex " // &
            & "carrier: a full field becomes a full OVERLAP field", &
            & nfail)
    end select

    call pd % get_real_vector(v)
    ok = .true.
    select type (part)
    type is (stored_graph)
       do i = 1, size(globals)
          ok = ok .and. (part % global_vertex_index(i) .eq. globals(i))
          ok = ok .and. (abs(v(i) - expect(i)) < 1.0d-13)
       end do
    end select
    call report(ok, &
         & "and carries the right value at every GLOBAL member it " // &
         & "holds", nfail)

  end subroutine check_one_transport

  !===================================================================!
  ! Assemble each part's field home: each contributes only what it
  ! OWNS, and the two together are exactly q* - no borrowed value
  ! counted twice, none dropped.
  !===================================================================!

  subroutine check_vertex_round_trip(nfail)

    integer, intent(inout) :: nfail

    type(field)           :: q
    real(dp)              :: from1(NV), from2(NV), total(NV)

    q = field('q star', g % vertex_set())
    call q % set_real_vector(Q_EXACT)

    from1 = 0.0_dp
    from2 = 0.0_dp
    call add_round_trip(q, 1, from1)
    call add_round_trip(q, 2, from2)

    call report(abs(from1(1) - 1.0_dp) < 1.0d-13 .and. &
         &      abs(from1(3) - 4.0_dp) < 1.0d-13 .and. &
         &      abs(from1(4)) < 1.0d-13, &
         & "part 1 contributes globals 1,2,3 and NOTHING at 4 - its " // &
         & "borrowed copy stays home", nfail)
    call report(abs(from2(4) - 7.0_dp) < 1.0d-13 .and. &
         &      abs(from2(6) - 16.0_dp) < 1.0d-13 .and. &
         &      abs(from2(3)) < 1.0d-13, &
         & "part 2 contributes globals 4,5,6 and nothing at 3", nfail)

    total = from1 + from2
    call report(maxval(abs(total - Q_EXACT)) < 1.0d-13, &
         & "and their sum is q* exactly: owned contributions tile " // &
         & "the whole", nfail)

  end subroutine check_vertex_round_trip

  !===================================================================!
  ! A PROPER global subset, declared out of global order, survives
  ! the same road: each part receives only the members it can see,
  ! read through local_index, and the assembled pieces recover S
  ! extensionally. Transformed subsets are not asked to keep S's
  ! identity token - only its members and values.
  !===================================================================!

  subroutine check_proper_subset_transport(nfail)

    integer, intent(inout) :: nfail

    type(counted_set)               :: vs
    type(subset_set)                :: s
    type(field)                     :: d
    type(partitioner)               :: p
    class(graph), allocatable       :: part
    class(graph_field), allocatable :: pd, fd
    class(member_set), allocatable  :: dom
    real(dp), allocatable           :: v(:)
    integer , allocatable           :: mem(:)
    real(dp)                        :: total(NV)
    integer                         :: k, i, gm
    logical                         :: ok, seen(NV)

    vs = g % vertex_set()
    s = subset_set('probe', vs, [6, 3, 4])      ! non-global order
    d = field('subset probe', s)
    call d % set_real_vector([600.0_dp, 300.0_dp, 400.0_dp])

    call report(abs(600.0_dp - value_at(d, s, 6)) < 1.0d-13 .and. &
         &      abs(300.0_dp - value_at(d, s, 3)) < 1.0d-13, &
         & "S = {6,3,4} carries its values in DECLARATION order", &
         & nfail)

    total = 0.0_dp
    seen  = .false.
    do k = 1, 2
       p = partitioner(PARTITION_LINEAR, nparts=2, part=k)
       call p % partition_graph(g, part)
       call p % partition_data(g, d, part, pd)

       ! Each part receives only the members it can see.
       call pd % domain(dom)
       select type (pp => part)
       type is (stored_graph)
          ok = .true.
          do i = 1, pp % num_vertices()
             if (dom % has(i)) then
                gm = pp % global_vertex_index(i)
                ok = ok .and. (gm .eq. 6 .or. gm .eq. 3 .or. gm .eq. 4)
                seen(gm) = .true.
             end if
          end do
       end select
       call report(ok, &
            & "the transported piece holds only members of S that " // &
            & "this part can see", nfail)

       ! The assembled piece is a proper SUBOBJECT of the global
       ! carrier, so its storage is its own - read every value
       ! through the assembled domain's own local_index.
       call a % assemble_data(part, pd, g, fd)
       call fd % domain(dom)
       call fd % get_real_vector(v)
       select type (dom)
       type is (subset_set)
          call report(dom % is_subobject_of(vs), &
               & "the assembled piece embeds in the global vertex " // &
               & "carrier", nfail)
          call dom % members(mem)
          do i = 1, size(mem)
             total(mem(i)) = total(mem(i)) + v(dom % local_index(mem(i)))
          end do
       class default
          call report(.false., &
               & "the assembled piece embeds in the global vertex " // &
               & "carrier", nfail)
       end select
    end do

    call report(seen(3) .and. seen(4) .and. seen(6), &
         & "and between them the parts saw all three members", nfail)

    call report(abs(total(6) - 600.0_dp) < 1.0d-13 .and. &
         &      abs(total(3) - 300.0_dp) < 1.0d-13 .and. &
         &      abs(total(4) - 400.0_dp) < 1.0d-13 .and. &
         &      abs(total(1)) < 1.0d-13 .and. &
         &      abs(total(2)) < 1.0d-13 .and. &
         &      abs(total(5)) < 1.0d-13, &
         & "the assembled pieces recover S extensionally, values " // &
         & "and absences alike", nfail)

  end subroutine check_proper_subset_transport

  !===================================================================!
  ! Helpers - each composed locally from the production API.
  !===================================================================!

  ! Partition d to part k, assemble it home, add it into total.
  subroutine add_round_trip(d, k, total)

    type(field), intent(in)    :: d
    integer    , intent(in)    :: k
    real(dp)   , intent(inout) :: total(:)

    type(partitioner)               :: p
    class(graph), allocatable       :: part
    class(graph_field), allocatable :: pd, fd
    real(dp), allocatable           :: v(:)

    p = partitioner(PARTITION_LINEAR, nparts=2, part=k)
    call p % partition_graph(g, part)
    call p % partition_data(g, d, part, pd)
    call a % assemble_data(part, pd, g, fd)
    call fd % get_real_vector(v)
    total = total + v(1:size(total))

  end subroutine add_round_trip

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

  real(dp) function value_at(d, dom, member)

    type(field)      , intent(in) :: d
    class(member_set), intent(in) :: dom
    integer          , intent(in) :: member

    real(dp), allocatable :: v(:)

    call d % get_real_vector(v)
    value_at = v(dom % local_index(member))

  end function value_at

end program partitioned_pde_gate_a
