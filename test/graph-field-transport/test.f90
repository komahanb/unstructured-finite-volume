!=====================================================================!
! The domain-transport suite (phase 5B acceptance): partition and
! assembly must carry FIELD DOMAINS, not merely values. A six-cell
! chain is cut linearly in two, so carrier identities differ, part
! numbering differs from global numbering, and borrowed members
! stand across the cut. Five laws, both families:
!
!      full vertex . full edge . proper vertex subset {6,3,1} .
!      proper edge subset {5,3,2} . empty subsets
!
! Proper subsets are declared in NONNUMERIC order with num_components=2, so
! domain indexing and component indexing are tested at once, and
! every value is fetched through local_index - never by assuming
! member equals position. Transported subobjects are NEW
! declarations: the round trip is extensional - same global
! ambient, same members, same values by member - never same token.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_graph_field_transport

  use iso_fortran_env        , only : dp => REAL64
  use graph_fractal           , only : graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set           , only : set_map
  use map_label         , only : label_map
  use map_inclusion     , only : inclusion_map, declared_subobject
  use view_directed    , only : directed_graph
  use field_calculus   , only : field
  use view_directed_stored            , only : stored_directed_graph
  use field_stored      , only : stored_field
  use transform_partitioner, only : partitioner, PARTITION_LINEAR
  use transform_assembler  , only : assembler

  use relation_partition, only : partition_relation
  implicit none
  type(partition_relation) :: rel

  type(stored_directed_graph) :: g
  type(assembler)    :: a
  integer            :: nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "graph field transport suite (phase 5B)"
  write(*,'(1x,a)') "============================================="

  g = stored_directed_graph(6, tails=[1,2,3,4,5], heads=[2,3,4,5,6])
  a = assembler()

  call check_full(.true. , 6, nfail)
  call check_full(.false., 5, nfail)
  call check_proper(.true. , [6, 3, 1], nfail)
  call check_proper(.false., [5, 3, 2], nfail)
  call check_empty(.true. , nfail)
  call check_empty(.false., nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all transport checks passed"
  else
     error stop
  end if

contains

  subroutine report(ok, label, nfail)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: label
    integer, intent(inout) :: nfail
    if (ok) then
       write(*,'(1x,a,a)') "PASS : ", label
    else
       write(*,'(1x,a,a)') "FAIL : ", label
       nfail = nfail + 1
    end if
  end subroutine report

  function family(v) result(s)
    logical, intent(in) :: v
    character(len=6) :: s
    if (v) then; s = 'vertex'; else; s = 'edge  '; end if
  end function family

  !===================================================================!
  ! Full fields: the established law - the two assembled owned
  ! contributions sum back to the whole.
  !===================================================================!

  subroutine check_full(verts, n, nfail)

    logical, intent(in)    :: verts
    integer, intent(in)    :: n
    integer, intent(inout) :: nfail

    type(stored_field)                     :: d
    type(partitioner)               :: p
    class(directed_graph), allocatable       :: part
    class(field), allocatable :: pd, fd
    type(set_map)                   :: sets
    type(label_map)                 :: labels
    type(inclusion_map)             :: inclusions
    real(dp), allocatable           :: v(:)
    real(dp)                        :: total(n)
    integer                         :: k, i

    if (verts) then
       call sets % bind(g % vertex_set(), &
            & counted_set_representation(g % num_vertices()))
       d = stored_field('q', g % vertex_set(), g % num_vertices())
    else
       call sets % bind(g % edge_set(), &
            & counted_set_representation(g % num_edges()))
       d = stored_field('q', g % edge_set(), g % num_edges())
    end if
    call d % set_real_vector([(10.0_dp * i, i = 1, n)])

    total = 0.0_dp
    do k = 1, 2
       p = partitioner(PARTITION_LINEAR, num_parts=2, part=k)
       call p % partition_graph(g, part, rel)
       call sets % bind(part % vertex_set(), &
            & counted_set_representation(part % num_vertices()))
       call sets % bind(part % edge_set(), &
            & counted_set_representation(part % num_edges()))
       call p % partition_data(rel, g, d, part, sets, labels, inclusions, pd)
       call a % assemble_data(rel, part, pd, g, sets, labels, inclusions, fd)
       call fd % real_vector(v)
       total = total + v(1:n)
    end do
    call report(all(abs(total - [(10.0_dp * i, i = 1, n)]) < 1.0d-13), &
         & "a full " // family(verts) // " field survives the round trip", &
         & nfail)

  end subroutine check_full

  !===================================================================!
  ! Proper subsets, nonnumeric order, two components: partitioned
  ! pieces are subobjects of the part carriers; assembled pieces
  ! are subobjects of the global carrier; their disjoint union by
  ! member reconstructs the source, values through local_index.
  !===================================================================!

  subroutine check_proper(verts, chosen, nfail)

    logical, intent(in)    :: verts
    integer, intent(in)    :: chosen(:)
    integer, intent(inout) :: nfail

    type(graph)                 :: carrier, s
    type(set_map)                   :: sets
    type(label_map)                 :: labels
    type(inclusion_map)             :: inclusions
    type(stored_field)                     :: d
    type(partitioner)               :: p
    class(directed_graph), allocatable       :: part
    class(field), allocatable :: pd, fd
    type(graph)                 :: dp_, dg
    real(dp), allocatable           :: sv(:), v(:)
    integer, allocatable            :: mem(:)
    integer                         :: k, i, c, m, seen(size(chosen))
    logical                         :: ok, okp

    if (verts) then
       carrier = g % vertex_set()
       call sets % bind(carrier, counted_set_representation(g % num_vertices()))
    else
       carrier = g % edge_set()
       call sets % bind(carrier, counted_set_representation(g % num_edges()))
    end if

    ! The chosen members are a new declared set, carved from the
    ! carrier: identity, extension and embedding together.
    call s % declare()
    call sets       % bind(s, listed_set_representation(chosen))
    call labels     % bind(s, 'chosen')
    call inclusions % include_in(s, carrier)

    d = stored_field('q', s, size(chosen), num_components=2)

    allocate(sv(2 * size(chosen)))
    do i = 1, size(chosen)
       sv(2 * i - 1) = 10.0_dp * chosen(i) + 1.0_dp
       sv(2 * i)     = 10.0_dp * chosen(i) + 2.0_dp
    end do
    call d % set_real_vector(sv)

    seen = 0
    ok   = .true.
    okp  = .true.
    do k = 1, 2
       p = partitioner(PARTITION_LINEAR, num_parts=2, part=k)
       call p % partition_graph(g, part, rel)
       call sets % bind(part % vertex_set(), &
            & counted_set_representation(part % num_vertices()))
       call sets % bind(part % edge_set(), &
            & counted_set_representation(part % num_edges()))

       call p % partition_data(rel, g, d, part, sets, labels, inclusions, pd)

       !-------------------------------------------------------------!
       ! The select type asked whether the domain was a subset TYPE.
       ! There is one domain type now, so the question that survives
       ! is the one that always mattered: was it CARVED from the part
       ! carrier. Provenance, asked of the map that records it.
       !-------------------------------------------------------------!

       dp_ = pd % domain()
       if (verts) then
          okp = okp .and. declared_subobject(dp_, part % vertex_set(), inclusions)
       else
          okp = okp .and. declared_subobject(dp_, part % edge_set(), inclusions)
       end if

       call a % assemble_data(rel, part, pd, g, sets, labels, inclusions, fd)

       dg = fd % domain()
       ok = ok .and. declared_subobject(dg, carrier, inclusions)
       call sets % members_of(dg, mem)
       call fd % real_vector(v)
       do i = 1, size(mem)
          m = mem(i)
          ! member identity: find m in the source declaration
          do c = 1, size(chosen)
             if (chosen(c) == m) seen(c) = seen(c) + 1
          end do
          ok = ok .and. &
               & abs(v((sets % index_in(dg, m) - 1) * 2 + 1) &
               &     - (10.0_dp * m + 1.0_dp)) < 1.0d-13 .and. &
               & abs(v((sets % index_in(dg, m) - 1) * 2 + 2) &
               &     - (10.0_dp * m + 2.0_dp)) < 1.0d-13
       end do
    end do

    call report(okp, &
         & "partitioned proper " // family(verts) // &
         & " subsets embed in the part carriers", nfail)
    call report(ok, &
         & "assembled contributions embed globally, values by member", &
         & nfail)
    call report(all(seen .eq. 1), &
         & "and their disjoint union rebuilds the source exactly once", &
         & nfail)

  end subroutine check_proper

  !===================================================================!
  ! The empty subset travels as itself: zero members, zero values,
  ! nothing manufactured, in either family.
  !===================================================================!

  subroutine check_empty(verts, nfail)

    logical, intent(in)    :: verts
    integer, intent(inout) :: nfail

    type(graph)                 :: carrier, s
    type(set_map)                   :: sets
    type(label_map)                 :: labels
    type(inclusion_map)             :: inclusions
    type(stored_field)                     :: d
    type(partitioner)               :: p
    class(directed_graph), allocatable       :: part
    class(field), allocatable :: pd, fd
    real(dp), allocatable           :: v(:)
    logical                         :: ok

    if (verts) then
       carrier = g % vertex_set()
       call sets % bind(carrier, counted_set_representation(g % num_vertices()))
    else
       carrier = g % edge_set()
       call sets % bind(carrier, counted_set_representation(g % num_edges()))
    end if

    call s % declare()
    call sets       % bind(s, listed_set_representation([integer ::]))
    call labels     % bind(s, 'none')
    call inclusions % include_in(s, carrier)

    d = stored_field('q', s, 0)
    call d % set_real_vector([real(dp) ::])

    p = partitioner(PARTITION_LINEAR, num_parts=2, part=1)
    call p % partition_graph(g, part, rel)
    call sets % bind(part % vertex_set(), &
         & counted_set_representation(part % num_vertices()))
    call sets % bind(part % edge_set(), &
         & counted_set_representation(part % num_edges()))
    call p % partition_data(rel, g, d, part, sets, labels, inclusions, pd)
    call a % assemble_data(rel, part, pd, g, sets, labels, inclusions, fd)
    call fd % real_vector(v)

    ok = pd % num_entries() .eq. 0 .and. fd % num_entries() .eq. 0 &
         & .and. size(v) .eq. 0
    call report(ok, &
         & "the empty " // family(verts) // &
         & " subset travels as itself: no fake member anywhere", nfail)

  end subroutine check_empty

end program test_graph_field_transport
