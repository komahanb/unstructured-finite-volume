!=====================================================================!
! The domain-transport suite (phase 5B acceptance): partition and
! assembly must carry FIELD DOMAINS, not merely values. A six-cell
! chain is cut linearly in two, so set identities differ, part
! numbering differs from global numbering, and borrowed members
! stand across the cut. Five laws, both families:
!
!      full vertex . full edge . proper vertex subset {6,3,1} .
!      proper edge subset {5,3,2} . empty subsets
!
! Proper subsets are declared in NONNUMERIC order with ncomp=2, so
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
  use graph_set          , only : set, index_set, subset
  use graph_grammar          , only : ordinary_graph, graph_field
  use class_graph            , only : stored_graph
  use class_graph_field      , only : field
  use class_graph_partitioner, only : partitioner, PARTITION_LINEAR
  use class_graph_assembler  , only : assembler

  implicit none

  type(stored_graph) :: g
  type(assembler)    :: a
  integer            :: nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "graph field transport suite (phase 5B)"
  write(*,'(1x,a)') "============================================="

  g = stored_graph(6, tails=[1,2,3,4,5], heads=[2,3,4,5,6])
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

    type(field)                     :: d
    type(partitioner)               :: p
    class(ordinary_graph), allocatable       :: part
    class(graph_field), allocatable :: pd, fd
    real(dp), allocatable           :: v(:)
    real(dp)                        :: total(n)
    integer                         :: k, i

    if (verts) then
       d = field('q', g % vertex_set())
    else
       d = field('q', g % edge_set())
    end if
    call d % set_real_vector([(10.0_dp * i, i = 1, n)])

    total = 0.0_dp
    do k = 1, 2
       p = partitioner(PARTITION_LINEAR, nparts=2, part=k)
       call p % partition_graph(g, part)
       call p % partition_data(g, d, part, pd)
       call a % assemble_data(part, pd, g, fd)
       call fd % get_real_vector(v)
       total = total + v(1:n)
    end do
    call report(all(abs(total - [(10.0_dp * i, i = 1, n)]) < 1.0d-13), &
         & "a full " // family(verts) // " field survives the round trip", &
         & nfail)

  end subroutine check_full

  !===================================================================!
  ! Proper subsets, nonnumeric order, two components: partitioned
  ! pieces are subobjects of the part sets; assembled pieces
  ! are subobjects of the global set; their disjoint union by
  ! member reconstructs the source, values through local_index.
  !===================================================================!

  subroutine check_proper(verts, chosen, nfail)

    logical, intent(in)    :: verts
    integer, intent(in)    :: chosen(:)
    integer, intent(inout) :: nfail

    type(index_set)               :: set
    type(subset)                :: s
    type(field)                     :: d
    type(partitioner)               :: p
    class(ordinary_graph), allocatable       :: part
    class(graph_field), allocatable :: pd, fd
    class(set), allocatable  :: dp_, dg
    real(dp), allocatable           :: sv(:), v(:)
    integer, allocatable            :: mem(:)
    integer                         :: k, i, c, m, seen(size(chosen))
    logical                         :: ok, okp

    if (verts) then
       set = g % vertex_set()
    else
       set = g % edge_set()
    end if
    s = subset('chosen', set, chosen)
    d = field('q', s, ncomp=2)

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
       p = partitioner(PARTITION_LINEAR, nparts=2, part=k)
       call p % partition_graph(g, part)
       call p % partition_data(g, d, part, pd)

       call pd % domain(dp_)
       select type (dp_)
       type is (subset)
          if (verts) then
             okp = okp .and. dp_ % is_subobject_of(part % vertex_set())
          else
             okp = okp .and. dp_ % is_subobject_of(part % edge_set())
          end if
       class default
          okp = .false.
       end select

       call a % assemble_data(part, pd, g, fd)
       call fd % domain(dg)
       select type (dg)
       type is (subset)
          ok = ok .and. dg % is_subobject_of(set)
          call dg % members(mem)
          call fd % get_real_vector(v)
          do i = 1, size(mem)
             m = mem(i)
             ! member identity: find m in the source declaration
             do c = 1, size(chosen)
                if (chosen(c) == m) seen(c) = seen(c) + 1
             end do
             ok = ok .and. &
                  & abs(v((dg % local_index(m) - 1) * 2 + 1) &
                  &     - (10.0_dp * m + 1.0_dp)) < 1.0d-13 .and. &
                  & abs(v((dg % local_index(m) - 1) * 2 + 2) &
                  &     - (10.0_dp * m + 2.0_dp)) < 1.0d-13
          end do
       class default
          ok = .false.
       end select
    end do

    call report(okp, &
         & "partitioned proper " // family(verts) // &
         & " subsets embed in the part sets", nfail)
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

    type(index_set)               :: set
    type(subset)                :: s
    type(field)                     :: d
    type(partitioner)               :: p
    class(ordinary_graph), allocatable       :: part
    class(graph_field), allocatable :: pd, fd
    real(dp), allocatable           :: v(:)
    logical                         :: ok

    if (verts) then
       set = g % vertex_set()
    else
       set = g % edge_set()
    end if
    s = subset('none', set, [integer ::])
    d = field('q', s)
    call d % set_real_vector([real(dp) ::])

    p = partitioner(PARTITION_LINEAR, nparts=2, part=1)
    call p % partition_graph(g, part)
    call p % partition_data(g, d, part, pd)
    call a % assemble_data(part, pd, g, fd)
    call fd % get_real_vector(v)

    ok = pd % num_entries() .eq. 0 .and. fd % num_entries() .eq. 0 &
         & .and. size(v) .eq. 0
    call report(ok, &
         & "the empty " // family(verts) // &
         & " subset travels as itself: no fake member anywhere", nfail)

  end subroutine check_empty

end program test_graph_field_transport
