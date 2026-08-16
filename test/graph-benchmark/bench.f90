!=====================================================================!
! The hot-traversal benchmark: the current stack, timed before the
! relational refactor moves in (AGENTS.md, Phase 0 / section 66).
!
! One structured grid graph, big enough that per-query overheads
! show: nx x ny cells as vertices, right and down neighbours as
! edges. The measured acts are the ones section 66 names:
!
!      incident query          incident_edges over every vertex
!      adjacent query          adjacent_vertices over every vertex
!      incoming/outgoing       both edge queries over every vertex
!      differential apply      divergence and laplacian, whole graph
!      field assembly          partition_data + assemble_data, 4 ways
!      partition construction  partition_graph, 4 ways
!
! The output is prose-free on purpose: one act, one time, one rate,
! so the same program rerun after any migration step reads as a
! before/after table. Numbers land in doc/phase0-benchmark.md.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program bench_graph_traversal

  use iso_fortran_env        , only : dp => REAL64, int64
  use graph_calculus         , only : GRAPH_SIDE_VERTEX, GRAPH_SIDE_EDGE
  use graph_grammar          , only : graph, graph_field
  use fractal_graph          , only : set_graph => graph
  use graph_set_representation, only : counted_set_representation
  use graph_set_map          , only : set_map
  use graph_binary_relation  , only : csr_relation
  use class_graph            , only : stored_graph
  use class_graph_support    , only : support
  use class_graph_field      , only : field
  use class_graph_partitioner, only : partitioner, PARTITION_LINEAR
  use class_graph_assembler  , only : assembler
  use class_graph_differential_operator, only : divergence, laplacian, &
       &                                        differential_operator

  implicit none

  integer, parameter :: nx = 700, ny = 700
  integer, parameter :: nv = nx * ny
  integer, parameter :: ne = (nx - 1) * ny + nx * (ny - 1)

  type(stored_graph)              :: g
  type(set_graph)                 :: vcarrier, ecarrier
  type(set_map)                   :: sets
  type(csr_relation), target      :: rel
  integer, pointer                :: fp(:)
  integer, allocatable            :: tab(:,:)
  type(differential_operator)     :: op
  type(partitioner)               :: p
  type(assembler)                 :: a
  class(graph), allocatable       :: part
  class(graph_field), allocatable :: pd, fd, yf
  type(support)                   :: von, eon
  type(field)                     :: q, z
  integer, allocatable            :: tails(:), heads(:), idx(:)
  real(dp), allocatable           :: qv(:), zv(:)
  integer(int64)                  :: t0, t1, rate
  integer(int64)                  :: touched
  integer                         :: i, j, e, v, k, rep
  real(dp)                        :: seconds

  write(*,'(1x,a)')       "============================================="
  write(*,'(1x,a)')       "hot traversal benchmark (AGENTS phase 0)"
  write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a)') "grid ", nx, " x ", ny, &
       & "  ->  ", nv, " vertices, ", ne, " edges"
  write(*,'(1x,a)')       "============================================="

  ! -- build the edge list: right neighbour, then down neighbour.
  allocate(tails(ne), heads(ne))
  e = 0
  do j = 1, ny
     do i = 1, nx
        v = (j - 1) * nx + i
        if (i < nx) then
           e = e + 1
           tails(e) = v
           heads(e) = v + 1
        end if
        if (j < ny) then
           e = e + 1
           tails(e) = v
           heads(e) = v + nx
        end if
     end do
  end do

  call system_clock(t0, rate)
  g = stored_graph(nv, tails=tails, heads=heads)
  call system_clock(t1)
  call line("stored_graph construction", t0, t1, rate, int(nv, int64))

  ! -- traversal sweeps: every vertex once per repetition.
  touched = 0
  call system_clock(t0)
  do rep = 1, 3
     do v = 1, nv
        call g % incident_edges(v, idx)
        touched = touched + size(idx)
     end do
  end do
  call system_clock(t1)
  call line("incident_edges sweep (x3)", t0, t1, rate, int(3, int64) * nv)

  call system_clock(t0)
  do rep = 1, 3
     do v = 1, nv
        call g % adjacent_vertices(v, idx)
        touched = touched + size(idx)
     end do
  end do
  call system_clock(t1)
  call line("adjacent_vertices sweep (x3)", t0, t1, rate, int(3, int64) * nv)

  call system_clock(t0)
  do rep = 1, 3
     do v = 1, nv
        call g % outgoing_edges(v, idx)
        touched = touched + size(idx)
     end do
  end do
  call system_clock(t1)
  call line("outgoing_edges sweep (x3)", t0, t1, rate, int(3, int64) * nv)

  call system_clock(t0)
  do rep = 1, 3
     do v = 1, nv
        call g % incoming_edges(v, idx)
        touched = touched + size(idx)
     end do
  end do
  call system_clock(t1)
  call line("incoming_edges sweep (x3)", t0, t1, rate, int(3, int64) * nv)

  ! -- the phase-3 citizen on the same topology: the tail relation
  !    V x E, image(v) = outgoing edges, preimage(e) = tail vertex.
  call vcarrier % declare()
  call ecarrier % declare()
  call sets % bind(vcarrier, counted_set_representation(nv))
  call sets % bind(ecarrier, counted_set_representation(ne))
  allocate(tab(2, ne))
  do e = 1, ne
     tab(1, e) = tails(e)
     tab(2, e) = e
  end do

  call system_clock(t0)
  rel = csr_relation('tail', vcarrier, ecarrier, tab, sets)
  call system_clock(t1)
  call line("csr_relation construction", t0, t1, rate, int(nv, int64))

  call system_clock(t0)
  do rep = 1, 3
     do v = 1, nv
        call rel % image(v, idx)
        touched = touched + size(idx)
     end do
  end do
  call system_clock(t1)
  call line("csr image sweep (x3)", t0, t1, rate, int(3, int64) * nv)

  call system_clock(t0)
  do rep = 1, 3
     do e = 1, ne
        call rel % preimage(e, idx)
        touched = touched + size(idx)
     end do
  end do
  call system_clock(t1)
  call line("csr preimage sweep (x3)", t0, t1, rate, int(3, int64) * ne)

  ! -- the borrows: the same fibres with no allocation anywhere.
  call system_clock(t0)
  do rep = 1, 3
     do v = 1, nv
        fp => rel % image_view(v)
        touched = touched + size(fp)
     end do
  end do
  call system_clock(t1)
  call line("csr image_view sweep (x3)", t0, t1, rate, int(3, int64) * nv)

  call system_clock(t0)
  do rep = 1, 3
     do e = 1, ne
        fp => rel % preimage_view(e)
        touched = touched + size(fp)
     end do
  end do
  call system_clock(t1)
  call line("csr preimage_view sweep (x3)", t0, t1, rate, int(3, int64) * ne)

  ! -- differential operators over the whole graph.
  eon = support(GRAPH_SIDE_EDGE, [(e, e = 1, ne)])
  z   = field('z', eon)
  allocate(zv(ne))
  do e = 1, ne
     zv(e) = real(mod(e, 97), dp) - 48.0_dp
  end do
  call z % set_real_vector(zv)

  op = divergence()
  call system_clock(t0)
  do rep = 1, 3
     call op % apply(g, [z], yf)
  end do
  call system_clock(t1)
  call line("divergence apply (x3)", t0, t1, rate, int(3, int64) * ne)

  von = support(GRAPH_SIDE_VERTEX, [(v, v = 1, nv)])
  q   = field('q', von)
  allocate(qv(nv))
  do v = 1, nv
     qv(v) = real(mod(v, 89), dp) - 44.0_dp
  end do
  call q % set_real_vector(qv)

  op = laplacian()
  call system_clock(t0)
  do rep = 1, 3
     call op % apply(g, [q], yf)
  end do
  call system_clock(t1)
  call line("laplacian apply (x3)", t0, t1, rate, int(3, int64) * ne)

  ! -- partition construction and the field round trip, four ways.
  a = assembler()
  call system_clock(t0)
  do k = 1, 4
     p = partitioner(PARTITION_LINEAR, nparts=4, part=k)
     call p % partition_graph(g, part)
  end do
  call system_clock(t1)
  call line("partition_graph, 4 parts", t0, t1, rate, int(nv, int64))

  call system_clock(t0)
  do k = 1, 4
     p = partitioner(PARTITION_LINEAR, nparts=4, part=k)
     call p % partition_graph(g, part)
     call p % partition_data(g, q, part, pd)
     call a % assemble_data(part, pd, g, fd)
  end do
  call system_clock(t1)
  call line("carry + assemble field, 4 parts", t0, t1, rate, int(nv, int64))

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a,i0)') "touched (keeps sweeps honest): ", touched

contains

  !===================================================================!
  ! line prints one benchmark row: what, how long, and the per-item
  ! rate in nanoseconds.
  !===================================================================!

  subroutine line(label, t0, t1, rate, items)

    character(len=*), intent(in) :: label
    integer(int64)  , intent(in) :: t0, t1, rate, items

    real(dp) :: s

    s = real(t1 - t0, dp) / real(rate, dp)
    write(*,'(1x,a32,f10.3,a,f10.1,a)') label, s, " s  ", &
         & 1.0d9 * s / real(items, dp), " ns/item"

  end subroutine line

end program bench_graph_traversal
