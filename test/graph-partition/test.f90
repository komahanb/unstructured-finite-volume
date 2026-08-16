!=====================================================================!
! THE PARTITION / ASSEMBLY LAW, FOR MORE THAN ONE PART
!
! What was proven before this suite: assemble(partition(G)) = G at
! nparts = 1, and the test that proves it calls that the weakest case.
! Nothing was proven for nparts > 1, and the maps therefore refused to
! call the pair a section.
!
! This suite asks what the CURRENT PUBLIC INTERFACES can actually
! express at nparts = 2 and 3. Nothing is migrated and no merge
! operation is invented so that an equation can be written; where the
! API cannot state a law, the suite records that it cannot.
!
! Four things are kept apart, because they have different laws:
!
!     A  graph structure   what assemble gives back from ONE part
!     B  vertex data       sum_i A_i(P_i(D)) = D , over owned entries
!     C  edge data         the same question one dimension over
!     D  ownership         every member owned exactly once
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program partition_law

  use iso_fortran_env      , only : dp => REAL64
  use class_graph          , only : directed_stored_graph
  use graph_directed_view  , only : directed_graph
  use graph_field_calculus , only : graph_field
  use class_graph_field    , only : field
  use fractal_graph      , only : set_graph => graph
  use graph_set_representation, only : counted_set_representation
  use graph_set_map      , only : set_map
  use graph_label_map    , only : label_map
  use graph_inclusion_map, only : inclusion_map
  use class_graph_partitioner, only : partitioner, PARTITION_LINEAR, &
       &                              PARTITION_BREADTH_FIRST
  use class_graph_assembler, only : assembler

  implicit none

  integer :: failures = 0

  write(*,'(1x,a)') "partition law suite"

  call structure_law(1)
  call structure_law(2)
  call structure_law(3)

  call vertex_data_law(1)
  call vertex_data_law(2)
  call vertex_data_law(3)

  call edge_data_law(2)
  call edge_data_law(3)

  call ownership_law(2)
  call ownership_law(3)

  call one_part_is_restricted_identity()

  if (failures .eq. 0) then
     print *, ''
     print *, ' ALL PROPOSITIONS HOLD'
  else
     print *, ''
     print *, ' FAILURES :', failures
     error stop 'partition_law: a proposition failed'
  end if

contains

  type(directed_stored_graph) function chain_of_six() result(g)

    g = directed_stored_graph(6, tails=[1, 2, 3, 4, 5], heads=[2, 3, 4, 5, 6])

  end function chain_of_six

  !===================================================================!
  ! A . STRUCTURE. assemble(partition(G, k)) renames one part's edges
  ! to whole-graph names. For k parts it is therefore a SUBGRAPH of G,
  ! not G - and the union over parts is not a public operation, so the
  ! statable law is COVERAGE, both ways:
  !
  !     every assembled edge is an edge of G          soundness
  !     every edge of G is assembled by some part     completeness
  !
  ! At nparts = 1 coverage collapses to equality, which is the law that
  ! was already proven.
  !===================================================================!

  subroutine structure_law(nparts)

    integer, intent(in) :: nparts

    type(directed_stored_graph)        :: g
    type(partitioner)         :: p
    type(assembler)           :: a
    class(directed_graph), allocatable :: part, back
    integer                   :: k, e, f
    integer                   :: seen(5)
    logical                   :: sound
    character(len=96)         :: what

    g = chain_of_six()
    a = assembler()
    seen  = 0
    sound = .true.

    do k = 1, nparts
       p = partitioner(PARTITION_LINEAR, nparts=nparts, part=k)
       call p % partition_graph(g, part)
       call a % assemble_graph(part, back)

       do e = 1, back % num_edges()
          f = whole_edge(g, back % edge_tail(e), back % edge_head(e))
          if (f .eq. 0) then
             sound = .false.
          else
             seen(f) = seen(f) + 1
          end if
       end do
    end do

    write(what, '(a,i0,a)') "A  nparts = ", nparts, &
         & ": every assembled edge is an edge of G"
    call check(trim(what), sound)

    write(what, '(a,i0,a)') "A  nparts = ", nparts, &
         & ": every edge of G is assembled by some part"
    call check(trim(what), all(seen .ge. 1))

  end subroutine structure_law

  !===================================================================!
  ! Which edge of G runs between these two whole-graph cells, or 0.
  !===================================================================!

  integer function whole_edge(g, t, h) result(which)

    type(directed_stored_graph), intent(in) :: g
    integer           , intent(in) :: t, h

    integer :: e

    which = 0
    do e = 1, g % num_edges()
       if (g % edge_tail(e) .eq. t .and. g % edge_head(e) .eq. h) then
          which = e
          return
       end if
    end do

  end function whole_edge

  !===================================================================!
  ! B . VERTEX DATA. The assembler writes only what a part OWNS and
  ! leaves the rest at zero, so the parts' answers add:
  !
  !     sum_k  A_k( P_k(D) )  =  D
  !
  ! This is a genuine law over more than one part, and it is the one
  ! the current API does express.
  !===================================================================!

  subroutine vertex_data_law(nparts)

    integer, intent(in) :: nparts

    type(directed_stored_graph)              :: g
    type(partitioner)               :: p
    type(assembler)                 :: a
    class(directed_graph), allocatable       :: part
    class(graph_field), allocatable :: pd, fd
    type(set_graph)                 :: on
    type(field)                     :: d

    !----------------------------------------------------------------!
    ! The transport CARVES: a part field lands on a domain that did
    ! not exist before the call, and the round trip carves again on
    ! the way home. These maps outlive both calls because the answers
    ! do.
    !----------------------------------------------------------------!

    type(set_map)                   :: sets
    type(label_map)                 :: labels
    type(inclusion_map)             :: inclusions
    real(dp)                        :: want(6), total(6)
    real(dp), allocatable           :: v(:)
    integer                         :: k
    character(len=96)               :: what

    g    = chain_of_six()
    a    = assembler()
    want = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 50.0_dp, 60.0_dp]

    !----------------------------------------------------------------!
    ! The whole graph's domains, described before anything transports
    ! data across them: the transform asks WHERE a global member
    ! stands, and only a representation can say.
    !----------------------------------------------------------------!

    on = g % vertex_set()
    call sets % bind(on, counted_set_representation(g % num_vertices()))

    d  = field('q', on, g % num_vertices())
    call d % set_real_vector(want)

    total = 0.0_dp
    do k = 1, nparts
       p = partitioner(PARTITION_LINEAR, nparts=nparts, part=k)
       call p % partition_graph(g, part)

       ! Each part is a new graph, so its own carriers are new domains
       ! and must be described before a field can be seated on them.
       call sets % bind(part % vertex_set(), &
            & counted_set_representation(part % num_vertices()))
       call sets % bind(part % edge_set(), &
            & counted_set_representation(part % num_edges()))

       call p % partition_data(g, d, part, sets, labels, inclusions, pd)
       call a % assemble_data(part, pd, g, sets, labels, inclusions, fd)

       select type (fd)
       class is (field)
          call fd % get_real_vector(v)
          if (size(v) .eq. 6) total = total + v
       end select
    end do

    write(what, '(a,i0,a)') "B  nparts = ", nparts, &
         & ": sum_k A_k(P_k(D)) = D on vertex data"
    call check(trim(what), all(abs(total - want) .lt. 1.0e-13_dp))

  end subroutine vertex_data_law

  !===================================================================!
  ! C . EDGE DATA. The same question one dimension over. Edges are not
  ! owned the way cells are - a cut edge stands on both sides - so this
  ! reports what the API answers rather than assuming the sum law
  ! carries over.
  !===================================================================!

  subroutine edge_data_law(nparts)

    integer, intent(in) :: nparts

    type(directed_stored_graph)              :: g
    type(partitioner)               :: p
    type(assembler)                 :: a
    class(directed_graph), allocatable       :: part
    class(graph_field), allocatable :: pd, fd
    type(set_graph)                 :: on
    type(field)                     :: d

    !----------------------------------------------------------------!
    ! The transport CARVES: a part field lands on a domain that did
    ! not exist before the call, and the round trip carves again on
    ! the way home. These maps outlive both calls because the answers
    ! do.
    !----------------------------------------------------------------!

    type(set_map)                   :: sets
    type(label_map)                 :: labels
    type(inclusion_map)             :: inclusions
    real(dp)                        :: want(5), total(5)
    real(dp), allocatable           :: v(:)
    integer                         :: k
    logical                         :: expressible
    character(len=96)               :: what

    g    = chain_of_six()
    a    = assembler()
    want = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]

    !----------------------------------------------------------------!
    ! The whole graph's domains, described before anything transports
    ! data across them: the transform asks WHERE a global member
    ! stands, and only a representation can say.
    !----------------------------------------------------------------!

    on = g % edge_set()
    call sets % bind(on, counted_set_representation(g % num_edges()))

    d  = field('w', on, g % num_edges())
    call d % set_real_vector(want)

    total       = 0.0_dp
    expressible = .true.
    do k = 1, nparts
       p = partitioner(PARTITION_LINEAR, nparts=nparts, part=k)
       call p % partition_graph(g, part)

       ! Each part is a new graph, so its own carriers are new domains
       ! and must be described before a field can be seated on them.
       call sets % bind(part % vertex_set(), &
            & counted_set_representation(part % num_vertices()))
       call sets % bind(part % edge_set(), &
            & counted_set_representation(part % num_edges()))

       call p % partition_data(g, d, part, sets, labels, inclusions, pd)
       call a % assemble_data(part, pd, g, sets, labels, inclusions, fd)

       select type (fd)
       class is (field)
          call fd % get_real_vector(v)
          if (size(v) .eq. 5) then
             total = total + v
          else
             expressible = .false.
          end if
       end select
    end do

    write(what, '(a,i0,a)') "C  nparts = ", nparts, &
         & ": edge data round trip is expressible"
    call check(trim(what), expressible)

    if (expressible) then
       write(what, '(a,i0,a)') "C  nparts = ", nparts, &
            & ": sum_k A_k(P_k(D)) = D on edge data"
       call check(trim(what), all(abs(total - want) .lt. 1.0e-13_dp))
    end if

  end subroutine edge_data_law

  !===================================================================!
  ! D . OWNERSHIP. Every whole-graph cell is owned by exactly one part.
  ! This is what actually holds at nparts > 1, and it is weaker than
  ! any round trip: it is a statement about the partition alone.
  !===================================================================!

  subroutine ownership_law(nparts)

    integer, intent(in) :: nparts

    type(directed_stored_graph)             :: g
    type(partitioner)              :: p
    class(directed_graph), allocatable      :: part
    type(set_graph)                :: owned, borrowed
    type(set_map)                  :: sets
    type(label_map)                :: labels
    type(inclusion_map)            :: inclusions
    integer                        :: times(6), k, l
    integer, allocatable           :: idx(:)
    logical                        :: borrows
    character(len=96)              :: what

    g     = chain_of_six()
    times = 0
    borrows = .false.

    do k = 1, nparts
       p = partitioner(PARTITION_LINEAR, nparts=nparts, part=k)
       call p % partition_graph(g, part)

       call part % owned_vertices(k, sets, labels, inclusions, owned)
       call roster(sets, owned, idx)
       do l = 1, size(idx)
          times(part % global_vertex_index(idx(l))) = &
               & times(part % global_vertex_index(idx(l))) + 1
       end do

       call part % borrowed_vertices(k, sets, labels, inclusions, borrowed)
       if (sets % size_of(borrowed) .gt. 0) borrows = .true.
    end do

    write(what, '(a,i0,a)') "D  nparts = ", nparts, &
         & ": every cell is owned exactly once"
    call check(trim(what), all(times .eq. 1))

    write(what, '(a,i0,a)') "D  nparts = ", nparts, &
         & ": and cells across a cut are borrowed, not owned twice"
    call check(trim(what), borrows)

  end subroutine ownership_law

  subroutine roster(sets, s, idx)

    type(set_map)       , intent(in)  :: sets
    type(set_graph)     , intent(in)  :: s
    integer, allocatable, intent(out) :: idx(:)

    integer :: k

    allocate(idx(sets % size_of(s)))
    do k = 1, sets % size_of(s)
       idx(k) = sets % member_of(s, k)
    end do

  end subroutine roster

  !===================================================================!
  ! The precise name for the one-part case. It is not a global section:
  ! it is the identity RESTRICTED to the one-part subdomain, and the
  ! restriction is where its whole force lies.
  !===================================================================!

  subroutine one_part_is_restricted_identity()

    type(directed_stored_graph)        :: g
    type(partitioner)         :: p
    type(assembler)           :: a
    class(directed_graph), allocatable :: part, back
    integer                   :: e
    logical                   :: same, some_part_is_short

    g = chain_of_six()
    a = assembler()
    p = partitioner(PARTITION_LINEAR, nparts=1, part=1)

    call p % partition_graph(g, part)
    call a % assemble_graph(part, back)

    same = back % num_vertices() .eq. g % num_vertices() .and. &
         & back % num_edges() .eq. g % num_edges()
    do e = 1, g % num_edges()
       same = same .and. back % edge_tail(e) .eq. g % edge_tail(e)
       same = same .and. back % edge_head(e) .eq. g % edge_head(e)
    end do

    call check('E  at nparts = 1 the round trip is the identity', same)

    ! And at nparts = 2 it is not: at least one part assembles back to
    ! strictly fewer edges than G has. That is what RESTRICTED means,
    ! and it is why the one-part law is not a global section.
    some_part_is_short = .false.
    do e = 1, 2
       p = partitioner(PARTITION_LINEAR, nparts=2, part=e)
       call p % partition_graph(g, part)
       call a % assemble_graph(part, back)
       if (back % num_edges() .lt. g % num_edges()) some_part_is_short = .true.
    end do

    call check('E  at nparts = 2 one part alone does NOT give G back', &
         & some_part_is_short)

  end subroutine one_part_is_restricted_identity

  subroutine check(label, ok)

    character(len=*), intent(in) :: label
    logical         , intent(in) :: ok

    if (ok) then
       print *, ' PASS : ', label
    else
       print *, ' FAIL : ', label
       failures = failures + 1
    end if

  end subroutine check

end program partition_law
