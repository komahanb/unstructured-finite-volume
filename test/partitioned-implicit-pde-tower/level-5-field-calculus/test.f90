!=====================================================================!
! PARTITIONED IMPLICIT PDE TOWER . LEVEL 5 . FIELD CALCULUS
!
! The level answers one question: HOW DO VALUES MOVE BETWEEN THE
! WHOLE AND ITS PART CARRIERS. Level 4 established what the parts
! ARE; this level establishes what may be read from them and what
! may be written back:
!
!      READ DOMAIN           = overlap  (owned union borrowed)
!      WRITE-BACK AUTHORITY  = owned
!
! A full global field becomes a FULL field on each part's whole
! vertex carrier - borrowed member included - because that is the
! shape a stencil will need one level above. But assembling a part
! home contributes only what that part OWNS, so the two
! contributions tile the whole exactly and no borrowed copy is
! counted twice.
!
! The same law is checked on edges, where it bites hardest: the
! crossing edge e3 lives in both parts, and an edge probe must
! come home unchanged - contributing 30 once, not twice, not never.
!
! Level 2 derived WHY that works, and the reason is narrower than
! it looks: the selected ownership policy is a TOTAL, SINGLE-VALUED
! map E -> K, so each edge answers exactly once. That property, not
! the tail anchor specifically, is what reconstruction needs -
! HeadOwner has it too. Level 4 then found which policy production
! actually selected (the tail-based one), and this level confirms
! the field machinery honours production's selection operationally.
!
! No operator appears here. Values move; nothing is computed.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program partitioned_pde_level_5

  use iso_fortran_env  , only : dp => REAL64
  use partitioned_pde_assert, only : report, verdict
  use partitioned_pde_assert, only : NV, NE, Q_EXACT
  use fractal_graph        , only : set_graph => graph
  use graph_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use graph_set_map        , only : set_map
  use graph_inclusion_map  , only : inclusion_map, declared_subobject
  use graph_label_map      , only : label_map
  use graph_grammar    , only : graph
  use graph_field_calculus, only : graph_field
  use class_graph      , only : stored_graph
  use class_graph_field, only : field
  use class_graph_partitioner, only : partitioner, PARTITION_LINEAR
  use class_graph_assembler  , only : assembler

  implicit none

  type(stored_graph)        :: g
  type(assembler)           :: a
  class(graph), allocatable :: g1, g2
  integer                   :: nfail
  type(set_map)     :: sets

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "partitioned pde tower . level 5 . fields"
  write(*,'(1x,a)') "============================================="

  g = stored_graph(NV, tails=[1,2,3,4,5], heads=[2,3,4,5,6])
  call sets % bind(g % vertex_set(), &
       & counted_set_representation(g % num_vertices()))
  call sets % bind(g % edge_set(), &
       & counted_set_representation(g % num_edges()))
  a = assembler()
  call cut(g1, 1)
  call cut(g2, 2)

  call check_edge_assembly_law(nfail)
  call check_vertex_transport(nfail)
  call check_vertex_round_trip(nfail)
  call check_proper_subset_transport(nfail)

  call verdict(nfail, "level 5")

contains
  !===================================================================!
  ! One part, cut from the global graph by the production rule.
  !===================================================================!

  subroutine cut(part, kpart)

    class(graph), allocatable, intent(out) :: part
    integer                  , intent(in)  :: kpart

    type(partitioner) :: p

    p = partitioner(PARTITION_LINEAR, nparts=2, part=kpart)
    call p % partition_graph(g, part)

    ! A part is a NEW graph, so its carriers are new declared domains
    ! and must be described before anything is seated on them.
    call sets % bind(part % vertex_set(), &
         & counted_set_representation(part % num_vertices()))
    call sets % bind(part % edge_set(), &
         & counted_set_representation(part % num_edges()))

  end subroutine cut
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
    type(label_map)     :: labels
    type(inclusion_map)     :: inclusions

    z = field('edge probe', g % edge_set(), g % num_edges())
    call z % set_real_vector(PROBE)

    allocate(total(NE))
    total = 0.0_dp
    call add_round_trip(z, 1, sets, labels, inclusions, total)
    call add_round_trip(z, 2, sets, labels, inclusions, total)

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
    z = field('vertex probe', g % vertex_set(), g % num_vertices())
    call z % set_real_vector(Q_EXACT)
    call add_round_trip(z, 1, sets, labels, inclusions, total)
    call add_round_trip(z, 2, sets, labels, inclusions, total)
    call report(maxval(abs(total - Q_EXACT)) < 1.0d-13, &
         & "and every global vertex likewise, borrowed copies " // &
         & "notwithstanding", nfail)

    if (allocated(got)) deallocate(got)

  end subroutine check_edge_assembly_law
  !===================================================================!
  ! A FULL global vertex field becomes a FULL overlap field on each
  ! part - values read by GLOBAL member, never by position.
  !===================================================================!

  subroutine check_vertex_transport(nfail)

    integer, intent(inout) :: nfail

    type(field) :: q

    q = field('q star', g % vertex_set(), g % num_vertices())
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
    type(set_graph)  :: dom
    type(inclusion_map)  :: inclusions
    type(label_map)  :: labels
    type(set_graph)               :: pvs
    real(dp), allocatable           :: v(:)
    character(len=1)                :: tag
    integer                         :: i
    logical                         :: ok

    write(tag,'(i1)') k

    ! The data is transported onto THIS part - not onto a second cut
    ! of the same graph, which would be a different carrier with a
    ! different identity.
    p = partitioner(PARTITION_LINEAR, nparts=2, part=k)
    call p % partition_data(g, q, part, sets, labels, inclusions, pd)

    dom = pd % domain()
    select type (part)
    type is (stored_graph)
       pvs = part % vertex_set()
       call report(dom % same_as(pvs) .and. sets % size_of(dom) .eq. 4, &
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
    type(label_map)     :: labels
    type(inclusion_map)     :: inclusions

    q = field('q star', g % vertex_set(), g % num_vertices())
    call q % set_real_vector(Q_EXACT)

    from1 = 0.0_dp
    from2 = 0.0_dp
    call add_round_trip(q, 1, sets, labels, inclusions, from1)
    call add_round_trip(q, 2, sets, labels, inclusions, from2)

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

    type(set_graph)               :: vs
    type(set_graph)                :: s
    type(field)                     :: d
    type(partitioner)               :: p
    class(graph), allocatable       :: part
    class(graph_field), allocatable :: pd, fd
    type(set_graph)  :: dom
    type(label_map)     :: labels
    type(inclusion_map)     :: inclusions
    real(dp), allocatable           :: v(:)
    integer , allocatable           :: mem(:)
    real(dp)                        :: total(NV)
    integer                         :: k, i, gm
    logical                         :: ok, seen(NV)

    vs = g % vertex_set()
    call s % declare()      ! non-global order
    call sets       % bind(s, listed_set_representation([6, 3, 4]))
    call inclusions % include_in(s, vs)
    d = field('subset probe', s, 3)
    call d % set_real_vector([600.0_dp, 300.0_dp, 400.0_dp])

    call report(abs(600.0_dp - value_at(sets, d, s, 6)) < 1.0d-13 .and. &
         &      abs(300.0_dp - value_at(sets, d, s, 3)) < 1.0d-13, &
         & "S = {6,3,4} carries its values in DECLARATION order", &
         & nfail)

    total = 0.0_dp
    seen  = .false.
    do k = 1, 2
       p = partitioner(PARTITION_LINEAR, nparts=2, part=k)
       call p % partition_graph(g, part)
       call p % partition_data(g, d, part, sets, labels, inclusions, pd)

       ! Each part receives only the members it can see.
       dom = pd % domain()
       select type (pp => part)
       type is (stored_graph)
          ok = .true.
          do i = 1, pp % num_vertices()
             if (sets % has_in(dom, i)) then
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
       call a % assemble_data(part, pd, g, sets, labels, inclusions, fd)
       dom = fd % domain()
       call fd % get_real_vector(v)
          call report(declared_subobject(dom, vs, inclusions), &
               & "the assembled piece embeds in the global vertex " // &
               & "carrier", nfail)
          call sets % members_of(dom, mem)
          do i = 1, size(mem)
             total(mem(i)) = total(mem(i)) + v(sets % index_in(dom, mem(i)))
          end do
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
  subroutine add_round_trip(d, k, sets, labels, inclusions, total)

    type(field)        , intent(in)    :: d
    integer            , intent(in)    :: k
    type(set_map)      , intent(inout) :: sets
    type(label_map)    , intent(inout) :: labels
    type(inclusion_map), intent(inout) :: inclusions
    real(dp)           , intent(inout) :: total(:)

    type(partitioner)               :: p
    class(graph), allocatable       :: part
    class(graph_field), allocatable :: pd, fd
    real(dp), allocatable           :: v(:)

    p = partitioner(PARTITION_LINEAR, nparts=2, part=k)
    call p % partition_graph(g, part)
    call sets % bind(part % vertex_set(), &
         & counted_set_representation(part % num_vertices()))
    call sets % bind(part % edge_set(), &
         & counted_set_representation(part % num_edges()))
    call p % partition_data(g, d, part, sets, labels, inclusions, pd)
    call a % assemble_data(part, pd, g, sets, labels, inclusions, fd)
    call fd % get_real_vector(v)
    total = total + v(1:size(total))

  end subroutine add_round_trip
  real(dp) function value_at(sets, d, dom, member)

    type(set_map)  , intent(in) :: sets
    type(field)    , intent(in) :: d
    type(set_graph), intent(in) :: dom
    integer        , intent(in) :: member

    real(dp), allocatable :: v(:)

    call d % get_real_vector(v)
    value_at = v(sets % index_in(dom, member))

  end function value_at
end program partitioned_pde_level_5
