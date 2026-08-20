!=====================================================================!
! PARTITIONED IMPLICIT PDE TOWER . LEVEL 6 . DISCRETIZATION
!
! The level answers one question: DOES THE SAME DISCRETE SPATIAL LAW
! APPLIED ON THE PARTS RECONSTRUCT THE LAW APPLIED ON THE WHOLE.
!
!      L   the production vertex Laplacian
!      A   the test-local shift,  A(q) = 2q - L(q)
!
! The Laplacian genuinely traverses whatever graph it is handed, so
! each part answers over its OWN incidence - and its answers at
! BORROWED members are wrong on purpose, because a borrowed vertex
! sits at the part's edge with half its stencil missing:
!
!      G1's copy of global 4 says 17;  the whole says 13
!      G2's copy of global 3 says  5;  the whole says  7
!
! Keep only what each part owns and the pieces reconstruct the
! global action exactly:
!
!      A_G q  =  sum_p  P_p^T O_p A_{G_p} P_p q
!
! Level 5 stated the read/write law structurally. This level proves
! it NUMERICALLY: perturbing only a borrowed copy moves an owned
! result across the cut, so visibility is not a courtesy - the
! stencil genuinely depends on it.
!
!      VISIBILITY governs what a local calculation may READ.
!      OWNERSHIP governs what it may authoritatively WRITE back.
!
! No solver appears here.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program partitioned_pde_level_6

  use iso_fortran_env  , only : dp => REAL64
  use partitioned_pde_assert, only : report, verdict
  use partitioned_pde_assert, only : NV, Q_EXACT, B_EXACT, L_EXACT
  use graph_fractal        , only : graph
  use map_set_representation, only : counted_set_representation
  use map_set        , only : set_map
  use map_inclusion  , only : inclusion_map, declared_subobject
  use map_label      , only : label_map
  use view_directed, only : directed_graph
  use field_calculus, only : field
  use view_directed_stored      , only : stored_directed_graph
  use field_stored, only : stored_field
  use transform_partitioner, only : partitioner, PARTITION_LINEAR
  use transform_assembler  , only : assembler
  use operation_differential, only : differential_operator, &
       &                                        laplacian
  use shifted_laplacian_fixture, only : shifted_laplacian

  use relation_partition, only : partition_relation
  implicit none
  type(partition_relation) :: rel

  type(stored_directed_graph)        :: g
  type(assembler)           :: a
  type(shifted_laplacian)   :: shifted
  class(directed_graph), allocatable :: g1, g2
  integer                   :: nfail
  type(set_map)     :: sets

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "partitioned pde tower . level 6 . operator"
  write(*,'(1x,a)') "============================================="

  g = stored_directed_graph(NV, tails=[1,2,3,4,5], heads=[2,3,4,5,6])
  call sets % bind(g % vertex_set(), &
       & counted_set_representation(g % num_vertices()))
  call sets % bind(g % edge_set(), &
       & counted_set_representation(g % num_edges()))
  a = assembler()
  call cut(g1, 1)
  call cut(g2, 2)

  call check_global_action(nfail)
  call check_transported_state(nfail)
  call check_local_actions(nfail)
  call check_assembled_equivalence(nfail)
  call check_borrowed_input_matters(nfail)

  call verdict(nfail, "level 6")

contains
  subroutine cut(part, k)

    class(directed_graph), allocatable, intent(out) :: part
    integer                  , intent(in)  :: k

    type(partitioner) :: p

    p = partitioner(PARTITION_LINEAR, num_parts=2, part=k)
    call p % partition_graph(g, part, rel)

  end subroutine cut
  !===================================================================!
  ! The global oracles, from the production Laplacian traversing G -
  ! never from a matrix.
  !===================================================================!

  subroutine check_global_action(nfail)

    integer, intent(inout) :: nfail

    type(differential_operator)     :: lap
    type(stored_field)                     :: q
    class(field), allocatable :: lq, aq
    type(graph)  :: dom
    real(dp), allocatable           :: v(:)
    type(graph)               :: vs
    integer         :: n_vs

    vs = g % vertex_set()
    n_vs = g % num_vertices()
    q = stored_field('q star', vs, n_vs)
    call q % set_real_vector(Q_EXACT)

    lap = laplacian(coefficient=1.0_dp, spacing=1.0_dp, measure=1.0_dp)
    call lap % apply(g, [q], lq)
    call lq % real_vector(v)
    call report(by_member(sets, v, vs, L_EXACT), &
         & "the production Laplacian on G gives " // &
         & "[1,1,1,1,1,-5], by member", nfail)

    call shifted % apply(g, [q], aq)
    dom = aq % domain()
    call aq % real_vector(v)
    call report(dom % same_as(vs) .and. by_member(sets, v, vs, B_EXACT), &
         & "and A q* = 2q* - Lq* = b = [1,3,7,13,21,37] on V(G)", &
         & nfail)

  end subroutine check_global_action
  !===================================================================!
  ! The minimal re-check Gate B needs before it may trust its own
  ! inputs: the transported state lives on the PART's whole vertex
  ! carrier, by identity, and carries the right values in the
  ! part's own enumeration. Gate A owns the full ownership story;
  ! this is only the handshake.
  !===================================================================!

  subroutine check_transported_state(nfail)

    integer, intent(inout) :: nfail

    call one_transported_state(g1, 1, &
         & [1.0_dp, 2.0_dp, 4.0_dp, 7.0_dp], nfail)
    call one_transported_state(g2, 2, &
         & [7.0_dp, 11.0_dp, 16.0_dp, 4.0_dp], nfail)

  end subroutine check_transported_state
  subroutine one_transported_state(part, k, expect, nfail)

    class(directed_graph), intent(in)    :: part
    integer     , intent(in)    :: k
    real(dp)    , intent(in)    :: expect(:)
    integer     , intent(inout) :: nfail

    type(partitioner)               :: p
    type(stored_field)                     :: q
    class(field), allocatable :: pd
    type(graph)  :: dom
    real(dp), allocatable           :: v(:)
    type(graph)               :: pvs
    type(label_map)     :: labels
    type(inclusion_map)     :: inclusions
    character(len=1)                :: tag
    type(partition_relation) :: rel

    write(tag,'(i1)') k

    q = stored_field('q star', g % vertex_set(), g % num_vertices())
    call q % set_real_vector(Q_EXACT)

    p = partitioner(PARTITION_LINEAR, num_parts=2, part=k)

    ! The relation of THIS part, not whichever cut ran last.
    select type (part)
    type is (stored_directed_graph)
       rel = part % whole_relation()
    end select

    call p % partition_data(rel, g, q, part, sets, labels, inclusions, pd)

    dom = pd % domain()
    call pd % real_vector(v)

    select type (part)
    type is (stored_directed_graph)
       pvs = part % vertex_set()
       call report(dom % same_as(pvs), &
            & "q" // tag // " lives on G" // tag // "'s whole vertex " // &
            & "carrier - the overlap a stencil will read", nfail)
    end select

    call report(size(v) .eq. size(expect) .and. &
         &      maxval(abs(v - expect)) < 1.0d-13, &
         & "and carries q* in G" // tag // "'s own enumeration", nfail)

  end subroutine one_transported_state
  !===================================================================!
  ! The same operation on each part, over the part's own topology.
  ! The complete local answers include BORROWED entries, and those
  ! are deliberately not the global ones.
  !===================================================================!

  subroutine check_local_actions(nfail)

    integer, intent(inout) :: nfail

    call one_local_action(g1, 1, [1,2,3,4], &
         & [1.0_dp, 1.0_dp, 1.0_dp, -3.0_dp], &
         & [1.0_dp, 3.0_dp, 7.0_dp, 17.0_dp], 4, 17.0_dp, 13.0_dp, nfail)
    call one_local_action(g2, 2, [4,5,6,3], &
         & [1.0_dp, 1.0_dp, -5.0_dp, 3.0_dp], &
         & [13.0_dp, 21.0_dp, 37.0_dp, 5.0_dp], 3, 5.0_dp, 7.0_dp, nfail)

  end subroutine check_local_actions
  subroutine one_local_action(part, k, globals, expect_l, expect_a, &
       & borrowed_global, borrowed_says, global_says, nfail)

    class(directed_graph), intent(in)    :: part
    integer     , intent(in)    :: k, globals(:), borrowed_global
    real(dp)    , intent(in)    :: expect_l(:), expect_a(:)
    real(dp)    , intent(in)    :: borrowed_says, global_says
    integer     , intent(inout) :: nfail

    type(differential_operator)     :: lap
    type(stored_field)                     :: qp
    class(field), allocatable :: lq, aq
    real(dp), allocatable           :: v(:)
    character(len=1)                :: tag
    integer                         :: i, seat
    logical                         :: ok
    type(label_map)     :: labels
    type(inclusion_map)     :: inclusions

    write(tag,'(i1)') k

    qp = local_state(part, k, sets, labels, inclusions)

    ! Every value is located through the part's OWN global map -
    ! never by array position - so a map inconsistency would show
    ! up here rather than hide behind a coincidental ordering.
    lap = laplacian(coefficient=1.0_dp, spacing=1.0_dp, measure=1.0_dp)
    call lap % apply(part, [qp], lq)
    call lq % real_vector(v)
    ok = .true.
    do i = 1, size(globals)
       seat = seat_of_global(part, globals(i))
       ok = ok .and. (seat .gt. 0)
       if (seat .gt. 0) then
          ok = ok .and. (abs(v(seat) - expect_l(i)) < 1.0d-12)
       end if
    end do
    call report(ok, &
         & "L on G" // tag // " traverses the PART's topology, read " // &
         & "by global member", nfail)

    call shifted % apply(part, [qp], aq)
    call aq % real_vector(v)
    ok = .true.
    do i = 1, size(globals)
       seat = seat_of_global(part, globals(i))
       ok = ok .and. (seat .gt. 0)
       if (seat .gt. 0) then
          ok = ok .and. (abs(v(seat) - expect_a(i)) < 1.0d-12)
       end if
    end do
    call report(ok, &
         & "and A on G" // tag // " answers on all four local " // &
         & "members, borrowed one included", nfail)

    ! The borrowed seat, located through the global map.
    seat = seat_of_global(part, borrowed_global)
    call report(seat .gt. 0 .and. &
         &      abs(v(seat) - borrowed_says) < 1.0d-12 .and. &
         &      abs(borrowed_says - global_says) > 1.0_dp, &
         & "G" // tag // "'s BORROWED output disagrees with the " // &
         & "global action - it is a copy, not an answer", nfail)

  end subroutine one_local_action
  !===================================================================!
  ! THE GATE-B THEOREM: assemble the two local actions, keeping only
  ! owned members, and the sum is the global action exactly.
  !===================================================================!

  subroutine check_assembled_equivalence(nfail)

    integer, intent(inout) :: nfail

    real(dp) :: total(NV)

    total = 0.0_dp
    call add_local_action(g1, 1, total)
    call add_local_action(g2, 2, total)

    call report(maxval(abs(total - B_EXACT)) < 1.0d-12, &
         & "assembled owned contributions reproduce A_G q* exactly: " // &
         & "the partitioned action IS the global action", nfail)
    call report(abs(total(3) - 7.0_dp) < 1.0d-12 .and. &
         &      abs(total(4) - 13.0_dp) < 1.0d-12, &
         & "and at the two interface vertices - where the borrowed " // &
         & "copies said 5 and 17 - the owned answers stand", nfail)

  end subroutine check_assembled_equivalence
  !===================================================================!
  ! The negative control: borrowed INPUT is numerically load-bearing.
  ! Perturb only the borrowed copy and an OWNED result must move.
  !===================================================================!

  subroutine check_borrowed_input_matters(nfail)

    integer, intent(inout) :: nfail

    call perturb_and_watch(g1, 1, 4, 3, 7.0_dp, -3.0_dp, nfail)
    call perturb_and_watch(g2, 2, 3, 4, 13.0_dp, 3.0_dp, nfail)

  end subroutine check_borrowed_input_matters
  subroutine perturb_and_watch(part, k, borrowed_global, watched_global, &
       & before, after, nfail)

    class(directed_graph), intent(in)    :: part
    integer     , intent(in)    :: k, borrowed_global, watched_global
    real(dp)    , intent(in)    :: before, after
    integer     , intent(inout) :: nfail

    type(stored_field)                     :: qp
    class(field), allocatable :: aq
    real(dp), allocatable           :: v(:)
    real(dp), allocatable           :: state(:)
    character(len=1)                :: tag
    integer                         :: bseat, wseat
    type(label_map)     :: labels
    type(inclusion_map)     :: inclusions

    write(tag,'(i1)') k

    bseat = seat_of_global(part, borrowed_global)
    wseat = seat_of_global(part, watched_global)

    ! Unperturbed: the owned answer stands where it should.
    qp = local_state(part, k, sets, labels, inclusions)
    call shifted % apply(part, [qp], aq)
    call aq % real_vector(v)
    call report(abs(v(wseat) - before) < 1.0d-12, &
         & "G" // tag // ": the owned result at global " // &
         & achar(48 + watched_global) // " starts correct", nfail)

    ! Perturb ONLY the borrowed copy, by its global identity.
    call qp % real_vector(state)
    state(bseat) = state(bseat) + 10.0_dp
    call qp % set_real_vector(state)

    call shifted % apply(part, [qp], aq)
    call aq % real_vector(v)
    call report(abs(v(wseat) - after) < 1.0d-12, &
         & "G" // tag // ": +10 on the BORROWED copy of global " // &
         & achar(48 + borrowed_global) // " moves the OWNED result " // &
         & "at global " // achar(48 + watched_global) // " by -10", &
         & nfail)

    ! Restore, and the correct answer returns.
    qp = local_state(part, k, sets, labels, inclusions)
    call shifted % apply(part, [qp], aq)
    call aq % real_vector(v)
    call report(abs(v(wseat) - before) < 1.0d-12, &
         & "G" // tag // ": restore the halo and equivalence returns", &
         & nfail)

  end subroutine perturb_and_watch
  !===================================================================!
  ! Helpers.
  !===================================================================!

  ! q* transported onto this part - the overlap-complete local state.
  type(stored_field) function local_state(part, k, sets, labels, inclusions) &
       & result(qp)

    class(directed_graph)       , intent(in)    :: part
    integer            , intent(in)    :: k
    type(set_map)      , intent(inout) :: sets
    type(label_map)    , intent(inout) :: labels
    type(inclusion_map), intent(inout) :: inclusions

    type(partitioner)               :: p
    type(stored_field)                     :: q
    class(field), allocatable :: pd
    real(dp), allocatable           :: v(:)
    type(partition_relation) :: rel

    q = stored_field('q star', g % vertex_set(), g % num_vertices())
    call q % set_real_vector(Q_EXACT)

    p = partitioner(PARTITION_LINEAR, num_parts=2, part=k)

    ! The relation of THIS part, not whichever cut ran last.
    select type (part)
    type is (stored_directed_graph)
       rel = part % whole_relation()
    end select

    call p % partition_data(rel, g, q, part, sets, labels, inclusions, pd)

    call pd % real_vector(v)
    qp = stored_field('local q', part % vertex_set(), part % num_vertices())
    call qp % set_real_vector(v)

  end function local_state
  ! Apply A on a part and add its OWNED contribution into total.
  subroutine add_local_action(part, k, total)

    class(directed_graph), intent(in)    :: part
    integer     , intent(in)    :: k
    real(dp)    , intent(inout) :: total(:)

    type(stored_field)                     :: qp, aq_local
    class(field), allocatable :: aq, fd
    type(graph)  :: dom
    type(label_map)     :: labels
    type(inclusion_map)     :: inclusions
    real(dp), allocatable           :: v(:)
    integer , allocatable           :: mem(:)
    integer                         :: i
    type(partition_relation) :: rel

    select type (part)
    type is (stored_directed_graph)
       rel = part % whole_relation()
    end select

    qp = local_state(part, k, sets, labels, inclusions)
    call shifted % apply(part, [qp], aq)
    call aq % real_vector(v)

    aq_local = stored_field('A local', part % vertex_set(), part % num_vertices())
    call aq_local % set_real_vector(v)

    ! The assembler is handed the relation of the part it gathers
    ! from. It binds nothing.
    call a % assemble_data(rel, part, aq_local, g, sets, labels, inclusions, fd)
    dom = fd % domain()
    call fd % real_vector(v)

       call sets % members_of(dom, mem)
       do i = 1, size(mem)
          total(mem(i)) = total(mem(i)) + v(sets % index_in(dom, mem(i)))
       end do

  end subroutine add_local_action
  ! The local seat holding this global member, or 0.
  integer function seat_of_global(part, gm)

    class(directed_graph), intent(in) :: part
    integer     , intent(in) :: gm

    integer :: i
    type(partition_relation) :: rel

    seat_of_global = 0
    select type (part)
    type is (stored_directed_graph)
       rel = part % whole_relation()
       do i = 1, part % num_vertices()
          if (rel % global_vertex_index(i) .eq. gm) seat_of_global = i
       end do
    end select

  end function seat_of_global
  ! Values compared through the domain's own map, never by position.
  logical function by_member(sets, v, dom, expect)

    type(set_map)  , intent(in) :: sets
    real(dp)       , intent(in) :: v(:)
    type(graph), intent(in) :: dom
    real(dp)       , intent(in) :: expect(:)

    integer :: i, m

    by_member = .true.
    do i = 1, sets % num_members_of(dom)
       m = sets % member_of(dom, i)
       by_member = by_member .and. &
            & (abs(v(sets % index_in(dom, m)) - expect(m)) < 1.0d-12)
    end do

  end function by_member
end program partitioned_pde_level_6
