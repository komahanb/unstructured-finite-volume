!=====================================================================!
! The characterization suite: the present, photographed before the
! relational refactor moves in (AGENTS.md, Phase 0).
!
! Nothing here states how the code OUGHT to behave. Every claim below
! was measured against the current stack, and the suite exists so
! that the migration to member sets and relations can prove, at every
! step, that these corners did not shift:
!
!      parallel edges         two edges, same endpoints, two identities
!      boundary half-edges    a wall needs no imaginary neighbour
!      directed traversal     out and in, edges and vertices
!      supports               membership, side, emptiness  (TRANSITIONAL)
!      partition round trip   vertex AND edge fields rebuild exactly
!      differential adjoints  the pairing holds on awkward topology
!
! One block is marked TRANSITIONAL: it pins vocabulary the tower has
! already sentenced (the support as an edgeless graph). Such checks
! guard the migration while it runs and are REWRITTEN, not obeyed,
! when their phase arrives - they must never veto the redesign they
! chaperone.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_graph_characterization

  use iso_fortran_env        , only : dp => REAL64
  use graph_calculus         , only : GRAPH_SIDE_VERTEX, GRAPH_SIDE_EDGE
  use graph_grammar          , only : graph, graph_field
  use class_graph            , only : stored_graph
  use class_graph_support    , only : support
  use class_graph_field      , only : field
  use class_graph_partitioner, only : partitioner, PARTITION_LINEAR
  use class_graph_assembler  , only : assembler
  use class_graph_differential_operator, only : divergence, &
       &                                        differential_operator, &
       &                                        vertex_differential_operator

  implicit none

  integer :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "graph characterization suite (AGENTS phase 0)"
  write(*,'(1x,a)') "============================================="

  call check_parallel_edges(nfail)
  call check_boundary_half_edges(nfail)
  call check_directed_traversal(nfail)
  call check_supports(nfail)
  call check_partition_round_trip(nfail)
  call check_adjoints_on_awkward_topology(nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all characterization checks passed"
  else
     write(*,'(1x,a,i0,a)') "FAILED: ", nfail, " characterization check(s)"
     error stop
  end if

contains

  subroutine report(ok, label, nfail)

    logical         , intent(in)    :: ok
    character(len=*), intent(in)    :: label
    integer         , intent(inout) :: nfail

    if (ok) then
       write(*,'(1x,a,a)') "PASS : ", label
    else
       write(*,'(1x,a,a)') "FAIL : ", label
       nfail = nfail + 1
    end if

  end subroutine report

  !===================================================================!
  ! Two edges joining the same pair of vertices are two edges. The
  ! multigraph identity: traversal reports both, the divergence
  ! counts both, and only the derived adjacency collapses them.
  !===================================================================!

  subroutine check_parallel_edges(nfail)

    integer, intent(inout) :: nfail

    type(stored_graph)              :: g
    type(differential_operator)     :: div
    class(graph_field), allocatable :: yf
    type(support)                   :: eon
    type(field)                     :: zf
    integer, allocatable            :: idx(:)
    real(dp), allocatable           :: y(:)

    g = stored_graph(2, tails=[1, 1], heads=[2, 2])

    call report(g % num_edges() .eq. 2, &
         & "parallel edges keep two identities", nfail)

    call g % incident_edges(1, idx)
    call report(size(idx) .eq. 2 .and. all(idx .eq. [1, 2]), &
         & "the shared tail is incident to both", nfail)

    call g % outgoing_edges(1, idx)
    call report(size(idx) .eq. 2 .and. all(idx .eq. [1, 2]), &
         & "both leave the tail", nfail)

    call g % incoming_edges(2, idx)
    call report(size(idx) .eq. 2 .and. all(idx .eq. [1, 2]), &
         & "both arrive at the head", nfail)

    ! The derived adjacency answers the neighbour once, however many
    ! edges carry the connection.
    call g % adjacent_vertices(1, idx)
    call report(size(idx) .eq. 1 .and. idx(1) .eq. 2, &
         & "derived adjacency names the neighbour once", nfail)

    ! The divergence spends each edge separately: out minus in.
    eon = support(GRAPH_SIDE_EDGE, [1, 2])
    zf  = field('z', eon)
    call zf % set_real_vector([3.0_dp, 5.0_dp])
    div = divergence()
    call div % apply(g, [zf], yf)
    call yf % get_real_vector(y)
    call report(abs(y(1) - 8.0_dp) < 1.0d-12 .and. &
         &      abs(y(2) + 8.0_dp) < 1.0d-12, &
         & "the divergence counts each parallel edge once", nfail)

  end subroutine check_parallel_edges

  !===================================================================!
  ! A boundary face is an edge with a tail and no head. It needs no
  ! imaginary member on the far side, it traverses like any edge from
  ! its one vertex, and its value leaves the graph exactly once.
  !===================================================================!

  subroutine check_boundary_half_edges(nfail)

    integer, intent(inout) :: nfail

    type(stored_graph)              :: g
    type(differential_operator)     :: div
    class(graph_field), allocatable :: yf
    type(support)                   :: eon
    type(field)                     :: zf
    class(graph), allocatable       :: sset
    integer, allocatable            :: idx(:)
    real(dp), allocatable           :: y(:)

    ! 1 --> 2 --> 3 --> wall
    g = stored_graph(3, tails=[1, 2, 3], heads=[2, 3, 0])

    call report(g % edge_has_head(1) .and. .not. g % edge_has_head(3), &
         & "the wall edge has a tail and no head", nfail)

    call g % interior_edges(sset)
    call members(sset, idx)
    call report(size(idx) .eq. 2 .and. all(idx .eq. [1, 2]), &
         & "interior edges are the two with heads", nfail)

    call g % boundary_edges(sset)
    call members(sset, idx)
    call report(size(idx) .eq. 1 .and. idx(1) .eq. 3, &
         & "the boundary set holds exactly the headless edge", nfail)

    call g % adjacent_vertices(3, idx)
    call report(size(idx) .eq. 1 .and. idx(1) .eq. 2, &
         & "no imaginary neighbour stands across the wall", nfail)

    call g % incident_edges(3, idx)
    call report(size(idx) .eq. 2 .and. all(idx .eq. [2, 3]), &
         & "the wall edge is incident to its one vertex", nfail)

    call g % outgoing_edges(3, idx)
    call report(size(idx) .eq. 1 .and. idx(1) .eq. 3, &
         & "and leaves it, as every edge leaves its tail", nfail)

    ! Divergence with samples 2, 4, 7: out minus in at each vertex
    ! gives 2, 2, 3 - and the total, 7, is the wall edge's sample:
    ! the half-edge contributes exactly once, to its tail alone.
    eon = support(GRAPH_SIDE_EDGE, [1, 2, 3])
    zf  = field('z', eon)
    call zf % set_real_vector([2.0_dp, 4.0_dp, 7.0_dp])
    div = divergence()
    call div % apply(g, [zf], yf)
    call yf % get_real_vector(y)
    call report(all(abs(y - [2.0_dp, 2.0_dp, 3.0_dp]) < 1.0d-12), &
         & "the half-edge contributes once, to its tail alone", nfail)
    call report(abs(sum(y) - 7.0_dp) < 1.0d-12, &
         & "what does not cancel is the flux through the wall", nfail)

  end subroutine check_boundary_half_edges

  !===================================================================!
  ! Directed traversal on a diamond: 1 --> 2 --> 4 and 1 --> 3 --> 4.
  ! Out and in, edges and vertices, sources and sinks - the queries
  ! the compatibility profile must go on answering.
  !===================================================================!

  subroutine check_directed_traversal(nfail)

    integer, intent(inout) :: nfail

    type(stored_graph)   :: g
    integer, allocatable :: idx(:)

    g = stored_graph(4, tails=[1, 1, 2, 3], heads=[2, 3, 4, 4])

    call g % outgoing_edges(1, idx)
    call report(size(idx) .eq. 2 .and. all(idx .eq. [1, 2]), &
         & "the source's outgoing edges", nfail)

    call g % incoming_edges(1, idx)
    call report(size(idx) .eq. 0, &
         & "a source has no incoming edge", nfail)

    call g % outgoing_vertices(1, idx)
    call report(size(idx) .eq. 2 .and. all(idx .eq. [2, 3]), &
         & "the source's downstream vertices", nfail)

    call g % incoming_edges(4, idx)
    call report(size(idx) .eq. 2 .and. all(idx .eq. [3, 4]), &
         & "the sink's incoming edges", nfail)

    call g % incoming_vertices(4, idx)
    call report(size(idx) .eq. 2 .and. all(idx .eq. [2, 3]), &
         & "the sink's upstream vertices", nfail)

    call g % outgoing_edges(4, idx)
    call report(size(idx) .eq. 0, &
         & "a sink has no outgoing edge", nfail)

  end subroutine check_directed_traversal

  !===================================================================!
  ! TRANSITIONAL - LEGACY COMPATIBILITY, NOT DESTINATION.
  !
  ! A support is a chosen set of members on one side of its host.
  ! Today it answers through the graph vocabulary - side constants,
  ! num_vertices, even num_edges == 0 - and these checks pin that
  ! TODAY, so nothing shifts unseen while the ground moves. The
  ! destination is different on purpose: support becomes a SUBOBJECT
  ! S c--> A (AGENTS.md sections 6 and 37, refined by review) - and
  ! THE DESTINATION NOW STANDS: subset_set in graph_carrier, with
  ! its own laws in test/graph-carrier. The pins below keep guarding
  ! the OLD support until the old fields retire onto the new
  ! carriers; then these checks are rewritten, not obeyed - what
  ! survives is membership, the host domain, order, and emptiness;
  ! the graph-flavoured spelling does not.
  !===================================================================!

  subroutine check_supports(nfail)

    integer, intent(inout) :: nfail

    type(support)        :: es, vs, none
    integer, allocatable :: idx(:)
    integer              :: k

    es = support(GRAPH_SIDE_EDGE, [11, 14, 19])
    call report(es % side() .eq. GRAPH_SIDE_EDGE, &
         & "a support knows which side its members reference", nfail)
    call report(es % num_vertices() .eq. 3, &
         & "and how many it holds", nfail)

    allocate(idx(es % num_vertices()))
    do k = 1, size(idx)
       idx(k) = es % global_vertex_index(k)
    end do
    call report(all(idx .eq. [11, 14, 19]), &
         & "members return exactly as given, in order", nfail)
    deallocate(idx)

    vs = support(GRAPH_SIDE_VERTEX, [2, 4, 6])
    call report(vs % side() .eq. GRAPH_SIDE_VERTEX, &
         & "the vertex side is the other side", nfail)

    call report(es % num_edges() .eq. 0 .and. vs % num_edges() .eq. 0, &
         & "a support carries members and no incidence", nfail)

    none = support(GRAPH_SIDE_VERTEX, [integer ::])
    call report(none % num_vertices() .eq. 0, &
         & "the empty support is a support", nfail)

  end subroutine check_supports

  !===================================================================!
  ! Partition then assemble rebuilds the whole - for a field on
  ! either side, through the same two calls. The transport symmetry
  ! the relational migration must keep while it collapses the
  ! vertex/edge branches into one domain-parametric act.
  !===================================================================!

  subroutine check_partition_round_trip(nfail)

    integer, intent(inout) :: nfail

    type(stored_graph)              :: g
    type(partitioner)               :: p
    type(assembler)                 :: a
    class(graph), allocatable       :: part
    class(graph_field), allocatable :: pd, fd
    type(support)                   :: von, eon
    type(field)                     :: q, w
    real(dp), allocatable           :: v(:)
    real(dp)                        :: vtotal(6), etotal(5)
    integer                         :: k

    g = stored_graph(6, tails=[1, 2, 3, 4, 5], heads=[2, 3, 4, 5, 6])
    a = assembler()

    von = support(GRAPH_SIDE_VERTEX, [1, 2, 3, 4, 5, 6])
    q   = field('q', von)
    call q % set_real_vector([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp])

    eon = support(GRAPH_SIDE_EDGE, [1, 2, 3, 4, 5])
    w   = field('w', eon)
    call w % set_real_vector([10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 50.0_dp])

    vtotal = 0.0_dp
    etotal = 0.0_dp

    do k = 1, 2
       p = partitioner(PARTITION_LINEAR, nparts=2, part=k)
       call p % partition_graph(g, part)

       call p % partition_data(g, q, part, pd)
       call a % assemble_data(part, pd, g, fd)
       select type (fd)
       class is (field)
          call fd % get_real_vector(v)
          vtotal = vtotal + v(1:6)
       end select

       call p % partition_data(g, w, part, pd)
       call a % assemble_data(part, pd, g, fd)
       select type (fd)
       class is (field)
          call fd % get_real_vector(v)
          etotal = etotal + v(1:5)
       end select
    end do

    call report(all(abs(vtotal - [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
         &                        5.0_dp, 6.0_dp]) < 1.0d-13), &
         & "a vertex field survives the round trip exactly", nfail)
    call report(all(abs(etotal - [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, &
         &                        50.0_dp]) < 1.0d-13), &
         & "and so does an edge field, through the same two calls", nfail)

  end subroutine check_partition_round_trip

  !===================================================================!
  ! The adjoint pairing on the topology this suite exists for: a
  ! graph holding a parallel pair AND a boundary half-edge at once.
  ! For every order, <A q, p> = <q, A* p> to machine precision.
  !===================================================================!

  subroutine check_adjoints_on_awkward_topology(nfail)

    integer, intent(inout) :: nfail

    type(stored_graph)              :: g
    type(differential_operator)     :: fwd, rev
    class(graph_field), allocatable :: yf
    type(support)                   :: on
    type(field)                     :: qf, pf
    real(dp), allocatable           :: aq(:), ap(:)
    real(dp)                        :: q(4), p(4), cs(5)
    integer                         :: v, order
    logical                         :: ok

    ! 1 ==> 2 (twice, in parallel), 2 --> 3 --> 4 --> wall.
    g  = stored_graph(4, tails=[1, 1, 2, 3, 4], heads=[2, 2, 3, 4, 0])
    on = support(GRAPH_SIDE_VERTEX, [1, 2, 3, 4])
    qf = field('q', on)
    pf = field('p', on)

    do v = 1, 4
       q(v) = real(v, dp)**2 - 3.0_dp * v
       p(v) = 2.0_dp * v + real(4 - v, dp)**2
    end do
    call qf % set_real_vector(q)
    call pf % set_real_vector(p)

    cs = [2.0_dp, 0.5_dp, 3.0_dp, 1.5_dp, 4.0_dp]

    ok = .true.
    do order = 1, 2
       fwd = vertex_differential_operator(order=order, coefficients=cs)
       rev = vertex_differential_operator(order=order, coefficients=cs, &
            &                             adjoint=.true.)

       call fwd % apply(g, [qf], yf)
       call yf % get_real_vector(aq)
       call rev % apply(g, [pf], yf)
       call yf % get_real_vector(ap)

       ok = ok .and. abs(sum(aq * p) - sum(q * ap)) < 1.0d-11
    end do

    call report(ok, &
         & "the pairing holds beside parallel edges and a wall", nfail)

  end subroutine check_adjoints_on_awkward_topology

  !===================================================================!
  ! members reads a named set back as plain indices.
  !===================================================================!

  subroutine members(g, indices)

    class(graph), intent(in)          :: g
    integer, allocatable, intent(out) :: indices(:)
    integer :: k

    allocate(indices(g % num_vertices()))
    do k = 1, size(indices)
       indices(k) = g % global_vertex_index(k)
    end do

  end subroutine members

end program test_graph_characterization
