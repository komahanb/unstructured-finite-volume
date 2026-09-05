!=====================================================================!
! The domain-transport suite (phase 5B acceptance): partition and
! assembly must transport FIELD DOMAINS, not merely values. A six-cell
! chain is cut linearly in two, so carrier identities differ, part
! numbering differs from global numbering, and shared members
! exist across the cut. Five laws, both families:
!
!      full vertex . full edge . proper vertex subset {6,3,1} .
!      proper edge subset {5,3,2} . empty subsets
!
! Proper subsets are declared in NONNUMERIC order with num_components=2, so
! domain indexing and component indexing are tested at once, and
! every value is read through index_in - never by assuming
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
  use map_set_store     , only : set_store
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

  subroutine report(passed, label, nfail)
    logical, intent(in) :: passed
    character(len=*), intent(in) :: label
    integer, intent(inout) :: nfail
    if (passed) then
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
    type(set_store)                 :: sets
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
       call p % partition_data(rel, g, d, part, sets, pd)
       call a % assemble_data(rel, part, pd, g, sets, fd)
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
  ! member reconstructs the source, values through index_in.
  !===================================================================!

  subroutine check_proper(verts, chosen, nfail)

    logical, intent(in)    :: verts
    integer, intent(in)    :: chosen(:)
    integer, intent(inout) :: nfail

    type(graph)                 :: carrier, s
    type(set_store)                 :: sets
    type(stored_field)                     :: d
    type(partitioner)               :: p
    class(directed_graph), allocatable       :: part
    class(field), allocatable :: pd, fd
    type(graph)                 :: dp_, dg
    real(dp), allocatable           :: sv(:), v(:)
    integer, allocatable            :: mem(:)
    integer                         :: k, i, c, m, counted(size(chosen))
    logical                         :: passed, passed_part

    if (verts) then
       carrier = g % vertex_set()
       call sets % bind(carrier, counted_set_representation(g % num_vertices()))
    else
       carrier = g % edge_set()
       call sets % bind(carrier, counted_set_representation(g % num_edges()))
    end if

    ! The chosen members are a new declared set, a subobject of the
    ! carrier: identity, extension and embedding together.
    call sets % declare_subobject(s, chosen, 'chosen', carrier)

    d = stored_field('q', s, size(chosen), num_components=2)

    allocate(sv(2 * size(chosen)))
    do i = 1, size(chosen)
       sv(2 * i - 1) = 10.0_dp * chosen(i) + 1.0_dp
       sv(2 * i)     = 10.0_dp * chosen(i) + 2.0_dp
    end do
    call d % set_real_vector(sv)

    counted = 0
    passed   = .true.
    passed_part  = .true.
    do k = 1, 2
       p = partitioner(PARTITION_LINEAR, num_parts=2, part=k)
       call p % partition_graph(g, part, rel)
       call sets % bind(part % vertex_set(), &
            & counted_set_representation(part % num_vertices()))
       call sets % bind(part % edge_set(), &
            & counted_set_representation(part % num_edges()))

       call p % partition_data(rel, g, d, part, sets, pd)

       !-------------------------------------------------------------!
       ! There is one domain type, so the question is whether the
       ! domain was declared a subobject of the part carrier: the
       ! provenance the set store records.
       !-------------------------------------------------------------!

       dp_ = pd % domain()
       if (verts) then
          passed_part = passed_part .and. sets % subobject_of(dp_, part % vertex_set())
       else
          passed_part = passed_part .and. sets % subobject_of(dp_, part % edge_set())
       end if

       call a % assemble_data(rel, part, pd, g, sets, fd)

       dg = fd % domain()
       passed = passed .and. sets % subobject_of(dg, carrier)
       call sets % members_of(dg, mem)
       call fd % real_vector(v)
       do i = 1, size(mem)
          m = mem(i)
          ! member identity: find m in the source declaration
          do c = 1, size(chosen)
             if (chosen(c) == m) counted(c) = counted(c) + 1
          end do
          passed = passed .and. &
               & abs(v((sets % index_in(dg, m) - 1) * 2 + 1) &
               &     - (10.0_dp * m + 1.0_dp)) < 1.0d-13 .and. &
               & abs(v((sets % index_in(dg, m) - 1) * 2 + 2) &
               &     - (10.0_dp * m + 2.0_dp)) < 1.0d-13
       end do
    end do

    call report(passed_part, &
         & "partitioned proper " // family(verts) // &
         & " subsets embed in the part carriers", nfail)
    call report(passed, &
         & "assembled contributions embed globally, values by member", &
         & nfail)
    call report(all(counted .eq. 1), &
         & "and their disjoint union rebuilds the source exactly once", &
         & nfail)

  end subroutine check_proper

  !===================================================================!
  ! The empty subset is transported as itself: zero members, zero
  ! values, no member manufactured, in either family.
  !===================================================================!

  subroutine check_empty(verts, nfail)

    logical, intent(in)    :: verts
    integer, intent(inout) :: nfail

    type(graph)                 :: carrier, s
    type(set_store)                 :: sets
    type(stored_field)                     :: d
    type(partitioner)               :: p
    class(directed_graph), allocatable       :: part
    class(field), allocatable :: pd, fd
    real(dp), allocatable           :: v(:)
    logical                         :: passed

    if (verts) then
       carrier = g % vertex_set()
       call sets % bind(carrier, counted_set_representation(g % num_vertices()))
    else
       carrier = g % edge_set()
       call sets % bind(carrier, counted_set_representation(g % num_edges()))
    end if

    call sets % declare_subobject(s, [integer ::], 'none', carrier)

    d = stored_field('q', s, 0)
    call d % set_real_vector([real(dp) ::])

    p = partitioner(PARTITION_LINEAR, num_parts=2, part=1)
    call p % partition_graph(g, part, rel)
    call sets % bind(part % vertex_set(), &
         & counted_set_representation(part % num_vertices()))
    call sets % bind(part % edge_set(), &
         & counted_set_representation(part % num_edges()))
    call p % partition_data(rel, g, d, part, sets, pd)
    call a % assemble_data(rel, part, pd, g, sets, fd)
    call fd % real_vector(v)

    passed = pd % num_entries() .eq. 0 .and. fd % num_entries() .eq. 0 &
         & .and. size(v) .eq. 0
    call report(passed, &
         & "the empty " // family(verts) // &
         & " subset is transported as itself: no member manufactured", nfail)

  end subroutine check_empty

end program test_graph_field_transport
