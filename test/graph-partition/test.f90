!=====================================================================!
! THE PARTITION / ASSEMBLY LAW, FOR MORE THAN ONE PART
!
! What was proven before this suite: assemble(partition(G)) = G at
! num_parts = 1, and the test that proves it calls that the weakest case.
! Nothing was proven for num_parts > 1, and the maps therefore refused to
! call the pair a section.
!
! This suite tests the mathematical laws that the CURRENT PUBLIC INTERFACES can
! express at num_parts = 2 and 3. Nothing is migrated and no merge
! operation is invented so that an equation can be written; where the
! API cannot state a law, the suite records that it cannot.
!
! Four objects are distinguished, because they have different laws:
!
!     A  graph structure   what assemble gives back from ONE part
!     B  vertex data       sum_i A_i(P_i(D)) = D , over owned entries
!     C  edge data         the same question one dimension over
!     D  ownership         every member owned exactly once
!     F  one relation      all four verbs driven by the SAME r_p
!     G  the wrong one     a relation from another part is refused
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program partition_law

  use iso_fortran_env      , only : dp => REAL64
  use view_directed_stored          , only : stored_directed_graph
  use view_directed  , only : directed_graph
  use field_calculus , only : field
  use field_stored    , only : stored_field
  use graph_fractal      , only : graph
  use map_set_representation, only : counted_set_representation
  use map_set_store, only : set_store
  use transform_partitioner, only : partitioner, PARTITION_LINEAR, &
       &                              PARTITION_BREADTH_FIRST
  use transform_assembler, only : assembler

  use relation_partition, only : partition_relation
  implicit none
  type(partition_relation) :: rel

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

  call check_partition_relation_consistency()
  call check_relation_identity()

  if (failures .eq. 0) then
     print *, ''
     print *, ' ALL PROPOSITIONS ARE TRUE'
  else
     print *, ''
     print *, ' FAILURES :', failures
     error stop 'partition_law: a proposition failed'
  end if

contains

  type(stored_directed_graph) function chain_of_six() result(g)

    g = stored_directed_graph(6, tails=[1, 2, 3, 4, 5], heads=[2, 3, 4, 5, 6])

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
  ! At num_parts = 1 coverage collapses to equality, which is the law that
  ! was already proven.
  !===================================================================!

  subroutine structure_law(num_parts)

    integer, intent(in) :: num_parts

    type(stored_directed_graph)        :: g
    type(partitioner)         :: p
    type(assembler)           :: a
    class(directed_graph), allocatable :: part, back
    integer                   :: k, e, f
    integer                   :: edge_counts(5)
    logical                   :: edges_valid
    character(len=96)         :: assertion_label

    g = chain_of_six()
    a = assembler()
    edge_counts  = 0
    edges_valid = .true.

    do k = 1, num_parts
       p = partitioner(PARTITION_LINEAR, num_parts=num_parts, part=k)
       call p % partition_graph(g, part, rel)
       call a % assemble_graph(rel, part, back)

       do e = 1, back % num_edges()
          f = whole_edge(g, back % edge_tail(e), back % edge_head(e))
          if (f .eq. 0) then
             edges_valid = .false.
          else
             edge_counts(f) = edge_counts(f) + 1
          end if
       end do
    end do

    write(assertion_label, '(a,i0,a)') "A  num_parts = ", num_parts, &
         & ": every assembled edge is an edge of G"
    call check(trim(assertion_label), edges_valid)

    write(assertion_label, '(a,i0,a)') "A  num_parts = ", num_parts, &
         & ": every edge of G is assembled by some part"
    call check(trim(assertion_label), all(edge_counts .ge. 1))

  end subroutine structure_law

  !===================================================================!
  ! Which edge of G runs between these two whole-graph cells, or 0.
  !===================================================================!

  integer function whole_edge(g, t, h) result(edge_index)

    type(stored_directed_graph), intent(in) :: g
    integer           , intent(in) :: t, h

    integer :: e

    edge_index = 0
    do e = 1, g % num_edges()
       if (g % edge_tail(e) .eq. t .and. g % edge_head(e) .eq. h) then
          edge_index = e
          return
       end if
    end do

  end function whole_edge

  !===================================================================!
  ! B . VERTEX DATA. The assembler writes only what a part OWNS and
  ! leaves the rest at zero, so the parts' contributions add:
  !
  !     sum_k  A_k( P_k(D) )  =  D
  !
  ! This is a genuine law over more than one part, and it is the one
  ! the current API does express.
  !===================================================================!

  subroutine vertex_data_law(num_parts)

    integer, intent(in) :: num_parts

    type(stored_directed_graph)              :: g
    type(partitioner)               :: p
    type(assembler)                 :: a
    class(directed_graph), allocatable       :: part
    class(field), allocatable :: pd, fd
    type(graph)                 :: on
    type(stored_field)                     :: d

    !----------------------------------------------------------------!
    ! The transport declares: a part field is defined on a domain
    ! that did not exist before the call, and the return declares
    ! again. The store outlives both calls because the results do.
    !----------------------------------------------------------------!

    type(set_store)                 :: sets
    real(dp)                        :: reference_values(6), total(6)
    real(dp), allocatable           :: v(:)
    integer                         :: k
    character(len=96)               :: assertion_label

    g    = chain_of_six()
    a    = assembler()
    reference_values = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 50.0_dp, 60.0_dp]

    !----------------------------------------------------------------!
    ! The whole graph's domains, described before anything transports
    ! data across them: the transform queries a global member
    ! position through the representation.
    !----------------------------------------------------------------!

    on = g % vertex_set()
    call sets % bind(on, counted_set_representation(g % num_vertices()))

    d  = stored_field('q', on, g % num_vertices())
    call d % set_real_vector(reference_values)

    total = 0.0_dp
    do k = 1, num_parts
       p = partitioner(PARTITION_LINEAR, num_parts=num_parts, part=k)
       call p % partition_graph(g, part, rel)

       ! Each part is a new graph, so its own carriers are new domains
       ! and must be described before a field can be defined on them.
       call sets % bind(part % vertex_set(), &
            & counted_set_representation(part % num_vertices()))
       call sets % bind(part % edge_set(), &
            & counted_set_representation(part % num_edges()))

       call p % partition_data(rel, g, d, part, sets, pd)
       call a % assemble_data(rel, part, pd, g, sets, fd)

       select type (fd)
       class is (stored_field)
          call fd % real_vector(v)
          if (size(v) .eq. 6) total = total + v
       end select
    end do

    write(assertion_label, '(a,i0,a)') "B  num_parts = ", num_parts, &
         & ": sum_k A_k(P_k(D)) = D on vertex data"
    call check(trim(assertion_label), all(abs(total - reference_values) .lt. 1.0e-13_dp))

  end subroutine vertex_data_law

  !===================================================================!
  ! C . EDGE DATA. The same question one dimension over. Edges are not
  ! owned the way cells are - a cut edge belongs to both parts - so this
  ! reports what the API returns rather than assuming the sum law
  ! transfers.
  !===================================================================!

  subroutine edge_data_law(num_parts)

    integer, intent(in) :: num_parts

    type(stored_directed_graph)              :: g
    type(partitioner)               :: p
    type(assembler)                 :: a
    class(directed_graph), allocatable       :: part
    class(field), allocatable :: pd, fd
    type(graph)                 :: on
    type(stored_field)                     :: d

    !----------------------------------------------------------------!
    ! The transport declares: a part field is defined on a domain
    ! that did not exist before the call, and the return declares
    ! again. The store outlives both calls because the results do.
    !----------------------------------------------------------------!

    type(set_store)                 :: sets
    real(dp)                        :: reference_values(5), total(5)
    real(dp), allocatable           :: v(:)
    integer                         :: k
    logical                         :: expressible
    character(len=96)               :: assertion_label

    g    = chain_of_six()
    a    = assembler()
    reference_values = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]

    !----------------------------------------------------------------!
    ! The whole graph's domains, described before anything transports
    ! data across them: the transform queries a global member
    ! position through the representation.
    !----------------------------------------------------------------!

    on = g % edge_set()
    call sets % bind(on, counted_set_representation(g % num_edges()))

    d  = stored_field('w', on, g % num_edges())
    call d % set_real_vector(reference_values)

    total       = 0.0_dp
    expressible = .true.
    do k = 1, num_parts
       p = partitioner(PARTITION_LINEAR, num_parts=num_parts, part=k)
       call p % partition_graph(g, part, rel)

       ! Each part is a new graph, so its own carriers are new domains
       ! and must be described before a field can be defined on them.
       call sets % bind(part % vertex_set(), &
            & counted_set_representation(part % num_vertices()))
       call sets % bind(part % edge_set(), &
            & counted_set_representation(part % num_edges()))

       call p % partition_data(rel, g, d, part, sets, pd)
       call a % assemble_data(rel, part, pd, g, sets, fd)

       select type (fd)
       class is (stored_field)
          call fd % real_vector(v)
          if (size(v) .eq. 5) then
             total = total + v
          else
             expressible = .false.
          end if
       end select
    end do

    write(assertion_label, '(a,i0,a)') "C  num_parts = ", num_parts, &
         & ": edge data round trip is expressible"
    call check(trim(assertion_label), expressible)

    if (expressible) then
       write(assertion_label, '(a,i0,a)') "C  num_parts = ", num_parts, &
            & ": sum_k A_k(P_k(D)) = D on edge data"
       call check(trim(assertion_label), all(abs(total - reference_values) .lt. 1.0e-13_dp))
    end if

  end subroutine edge_data_law

  !===================================================================!
  ! D . OWNERSHIP. Every whole-graph cell is owned by exactly one part.
  ! This is the ownership law at num_parts > 1, and it is weaker than any
  ! round trip: it is a statement about the partition alone. The
  ! owner of each part member is read from the partition relation;
  ! a member whose owner is another part is a halo member across a cut.
  !===================================================================!

  subroutine ownership_law(num_parts)

    integer, intent(in) :: num_parts

    type(stored_directed_graph)             :: g
    type(partitioner)              :: p
    class(directed_graph), allocatable      :: part
    integer                        :: times(6), k, l, gm
    logical                        :: has_halo_members
    character(len=96)              :: assertion_label

    g     = chain_of_six()
    times = 0
    has_halo_members = .false.

    do k = 1, num_parts
       p = partitioner(PARTITION_LINEAR, num_parts=num_parts, part=k)
       call p % partition_graph(g, part, rel)

       do l = 1, part % num_vertices()
          gm = rel % global_vertex_index(l)
          if (rel % vertex_owner_part(l) .eq. rel % part_id()) then
             times(gm) = times(gm) + 1
          else
             has_halo_members = .true.
          end if
       end do
    end do

    write(assertion_label, '(a,i0,a)') "D  num_parts = ", num_parts, &
         & ": every cell is owned exactly once"
    call check(trim(assertion_label), all(times .eq. 1))

    write(assertion_label, '(a,i0,a)') "D  num_parts = ", num_parts, &
         & ": and cells across a cut are halo members, each with one owner"
    call check(trim(assertion_label), has_halo_members)

  end subroutine ownership_law

  !===================================================================!
  ! The precise name for the one-part case. It is not a global section:
  ! it is the identity RESTRICTED to the one-part subdomain, and the
  ! restriction is where its whole force lies.
  !===================================================================!

  subroutine one_part_is_restricted_identity()

    type(stored_directed_graph)        :: g
    type(partitioner)         :: p
    type(assembler)           :: a
    class(directed_graph), allocatable :: part, back
    integer                   :: e
    logical                   :: same, has_proper_edge_subset

    g = chain_of_six()
    a = assembler()
    p = partitioner(PARTITION_LINEAR, num_parts=1, part=1)

    call p % partition_graph(g, part, rel)
    call a % assemble_graph(rel, part, back)

    same = back % num_vertices() .eq. g % num_vertices() .and. &
         & back % num_edges() .eq. g % num_edges()
    do e = 1, g % num_edges()
       same = same .and. back % edge_tail(e) .eq. g % edge_tail(e)
       same = same .and. back % edge_head(e) .eq. g % edge_head(e)
    end do

    call check('E  at num_parts = 1 the round trip is the identity', same)

    ! And at num_parts = 2 it is not: at least one part assembles back to
    ! strictly fewer edges than G has. That is what RESTRICTED means,
    ! and it is why the one-part law is not a global section.
    has_proper_edge_subset = .false.
    do e = 1, 2
       p = partitioner(PARTITION_LINEAR, num_parts=2, part=e)
       call p % partition_graph(g, part, rel)
       call a % assemble_graph(rel, part, back)
       if (back % num_edges() .lt. g % num_edges()) has_proper_edge_subset = .true.
    end do

    call check('E  at num_parts = 2 one part alone does NOT give G back', &
         & has_proper_edge_subset)

  end subroutine one_part_is_restricted_identity

  !===================================================================!
  ! F . ONE RELATION, FOUR VERBS.
  !
  ! r_p <= S_part x S_whole is written ONCE, by the cut, and the three
  ! verbs that follow receive it:
  !
  !     partition_graph(G)      -> G_p , r        r written here
  !     partition_data (r, D)   -> D_p            read forward
  !     assemble_graph (r, G_p) -> G contribution read backward
  !     assemble_data  (r, D_p) -> D contribution read backward
  !
  ! The claim this makes is not that the four agree - it is that there
  ! is nothing for them to disagree ABOUT, because exactly one value
  ! reaches all four. So the test is written to make a second relation
  ! impossible: `rel_p` is set by partition_graph and is never
  ! assigned again, and every later verb names it.
  !
  ! What it checks afterwards is the consequence: the values that come
  ! map to the global members specified by r_p.
  !===================================================================!

  subroutine check_partition_relation_consistency()

    type(stored_directed_graph)        :: g
    type(partitioner)         :: p
    type(assembler)           :: a
    class(directed_graph), allocatable :: part, back
    class(field), allocatable    :: pd, fd
    type(partition_relation)  :: rel_p
    type(set_store)           :: sets
    type(stored_field)               :: d
    type(graph)           :: on
    real(dp), allocatable     :: assembled_values(:)
    integer                   :: l, gm
    logical                   :: satisfied

    g = chain_of_six()
    call sets % bind(g % vertex_set(), &
         & counted_set_representation(g % num_vertices()))

    ! 1. THE CUT WRITES r. Nothing else in this routine writes rel_p.
    p = partitioner(PARTITION_LINEAR, num_parts=2, part=2)
    call p % partition_graph(g, part, rel_p)
    call sets % bind(part % vertex_set(), &
         & counted_set_representation(part % num_vertices()))

    on = g % vertex_set()
    d  = stored_field('q', on, sets % num_members_of(on))
    call d % set_real_vector([10.0_dp, 20.0_dp, 30.0_dp, &
         &                    40.0_dp, 50.0_dp, 60.0_dp])

    ! 2, 3, 4. THE OTHER THREE RECEIVE THE SAME r.
    a = assembler()
    call check('F  the partition relation describes the part', &
         & a % defined_on_relation(rel_p, part))

    call p % partition_data(rel_p, g, d, part, sets, pd)
    call a % assemble_graph(rel_p, part, back)
    call a % assemble_data(rel_p, part, pd, g, sets, fd)

    ! The inverse relation maps the structure to the whole.
    call check('F  assemble_graph put the part back in whole names', &
         & back % num_vertices() .eq. g % num_vertices())

    ! And so did the values: every entry this part OWNS is at
    ! the member r_p names, and every entry it does not own is zero -
    ! which is what stops the other part''s copy being counted twice.
    call fd % real_vector(assembled_values)
    satisfied = .true.
    do l = 1, part % num_vertices()
       gm = rel_p % global_vertex_index(l)
       if (rel_p % vertex_owner_part(l) .eq. rel_p % part_id()) then
          satisfied = satisfied .and. abs(assembled_values(gm) - 10.0_dp * gm) .lt. 1.0e-12_dp
       end if
    end do
    call check('F  every owned value maps to its specified global member', satisfied)

  end subroutine check_partition_relation_consistency

  !===================================================================!
  ! G . THE WRONG RELATION IS REFUSED.
  !
  ! ONE RELATION PER PART. NEVER ONE ACROSS TWO. A caller cutting a
  ! graph in two has two relations, and passing part 2 the relation
  ! written for part 1 is the exact mistake this design exists to make
  ! impossible - it was a real defect, found when a single relation was
  ! reused across two parts and the second part assembled to the whole through
  ! the first part''s numbering.
  !
  ! So it must be REPORTABLE, not merely wrong. defined_on_relation
  ! returns .false., and it decides on identity: r1 and r2 here have
  ! equal shape in every respect a count could see, and are still not
  ! interchangeable, because they relate different sets.
  !===================================================================!

  subroutine check_relation_identity()

    type(stored_directed_graph)        :: g
    type(partitioner)         :: p1, p2
    type(assembler)           :: a
    class(directed_graph), allocatable :: part1, part2
    type(partition_relation)  :: r1, r2

    g = chain_of_six()
    a = assembler()

    p1 = partitioner(PARTITION_LINEAR, num_parts=2, part=1)
    p2 = partitioner(PARTITION_LINEAR, num_parts=2, part=2)
    call p1 % partition_graph(g, part1, r1)
    call p2 % partition_graph(g, part2, r2)

    call check('G  each relation describes its part', &
         & a % defined_on_relation(r1, part1) .and. &
         & a % defined_on_relation(r2, part2))

    call check('G  and in neither of the other''s: the wrong relation ' // &
         & 'is REFUSED, not obeyed', &
         & .not. a % defined_on_relation(r1, part2) .and. &
         & .not. a % defined_on_relation(r2, part1))

    ! The refusal is not a size argument. Both parts have four cells
    ! and three edges here, so a check on counts alone would pass the
    ! wrong relation straight through.
    call check('G  and the two are indistinguishable by size, so ' // &
         & 'IDENTITY is what refused them', &
         & part1 % num_vertices() .eq. part2 % num_vertices() .and. &
         & part1 % num_edges()    .eq. part2 % num_edges())

  end subroutine check_relation_identity

  subroutine check(label, satisfied)

    character(len=*), intent(in) :: label
    logical         , intent(in) :: satisfied

    if (satisfied) then
       print *, ' PASS : ', label
    else
       print *, ' FAIL : ', label
       failures = failures + 1
    end if

  end subroutine check

end program partition_law
