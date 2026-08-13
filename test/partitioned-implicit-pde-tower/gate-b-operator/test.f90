!=====================================================================!
! PARTITIONED IMPLICIT PDE TOWER . GATE B . TOPOLOGY-CONSUMING ACTION
!
! The gate answers TWO different questions, and proves them apart.
!
! FIRST: does a real Class-1 operation, evaluated on
! overlap-complete part fields and assembled through OWNED outputs,
! reproduce the global operator?
!
!      A_G q  =  sum_p  P_p^T O_p A_{G_p} P_p q
!
! It does, exactly - and the borrowed outputs it discards are
! visibly wrong: G1's copy of global vertex 4 answers 17 where the
! global action answers 13, and G2's copy of global vertex 3
! answers 5 where the global action answers 7. A part holds enough
! topology to answer for what it OWNS and no more. Borrowed inputs
! are necessary; borrowed outputs are not authoritative.
!
! SECOND: is the graph a minimizer carries actually load-bearing?
! Four earlier towers reported it as scenery - correctly, about
! their own actions, which had no topology to traverse. Here the
! attached action does. The same shifted_laplacian type, which
! stores no graph at all, is attached to two solvers over two
! DIFFERENT six-vertex five-edge topologies, and handed the same
! probe:
!
!      solver on the chain -> b
!      solver on the star  -> something else entirely
!
! Nothing changed but the host. The minimizer never reads topology;
! it CARRIES the graph to an operation that does. Those are two
! roles, and this gate is where the distinction stops being an
! inference from production code and becomes a measurement.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program partitioned_pde_gate_b

  use iso_fortran_env  , only : dp => REAL64
  use partitioned_pde_assert, only : report, verdict
  use partitioned_pde_assert, only : NV, NE, Q_EXACT, B_EXACT, L_EXACT
  use graph_carrier    , only : counted_set, subset_set, member_set
  use graph_grammar    , only : graph, graph_field
  use class_graph      , only : stored_graph
  use class_graph_field, only : field
  use class_graph_partitioner, only : partitioner, PARTITION_LINEAR
  use class_graph_assembler  , only : assembler
  use class_graph_differential_operator, only : differential_operator, &
       &                                        laplacian
  use class_graph_gmres, only : gmres
  use shifted_laplacian_fixture, only : shifted_laplacian

  implicit none

  type(stored_graph)        :: g, g_alt
  type(assembler)           :: a
  type(shifted_laplacian)   :: shifted
  class(graph), allocatable :: g1, g2
  integer                   :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "partitioned pde tower . gate B . the action"
  write(*,'(1x,a)') "============================================="

  g = stored_graph(NV, tails=[1,2,3,4,5], heads=[2,3,4,5,6])
  a = assembler()

  call cut(g1, 1)
  call cut(g2, 2)

  call check_global_action(nfail)
  call check_local_actions(nfail)
  call check_assembled_equivalence(nfail)
  call check_borrowed_input_matters(nfail)
  call check_global_solve(nfail)
  call check_host_is_load_bearing(nfail)

  call verdict(nfail, "gate B")

contains

  subroutine cut(part, k)

    class(graph), allocatable, intent(out) :: part
    integer                  , intent(in)  :: k

    type(partitioner) :: p

    p = partitioner(PARTITION_LINEAR, nparts=2, part=k)
    call p % partition_graph(g, part)

  end subroutine cut

  !===================================================================!
  ! The global oracles, from the production Laplacian traversing G -
  ! never from a matrix.
  !===================================================================!

  subroutine check_global_action(nfail)

    integer, intent(inout) :: nfail

    type(differential_operator)     :: lap
    type(field)                     :: q
    class(graph_field), allocatable :: lq, aq
    class(member_set), allocatable  :: dom
    real(dp), allocatable           :: v(:)
    type(counted_set)               :: vs

    vs = g % vertex_set()
    q = field('q star', vs)
    call q % set_real_vector(Q_EXACT)

    lap = laplacian(coefficient=1.0_dp, spacing=1.0_dp, measure=1.0_dp)
    call lap % apply(g, [q], lq)
    call lq % get_real_vector(v)
    call report(by_member(v, vs, L_EXACT), &
         & "the production Laplacian on G gives " // &
         & "[1,1,1,1,1,-5], by member", nfail)

    call shifted % apply(g, [q], aq)
    call aq % domain(dom)
    call aq % get_real_vector(v)
    call report(dom % same_as(vs) .and. by_member(v, vs, B_EXACT), &
         & "and A q* = 2q* - Lq* = b = [1,3,7,13,21,37] on V(G)", &
         & nfail)

  end subroutine check_global_action

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

    class(graph), intent(in)    :: part
    integer     , intent(in)    :: k, globals(:), borrowed_global
    real(dp)    , intent(in)    :: expect_l(:), expect_a(:)
    real(dp)    , intent(in)    :: borrowed_says, global_says
    integer     , intent(inout) :: nfail

    type(differential_operator)     :: lap
    type(field)                     :: qp
    class(graph_field), allocatable :: lq, aq
    real(dp), allocatable           :: v(:)
    character(len=1)                :: tag
    integer                         :: i, seat
    logical                         :: ok

    write(tag,'(i1)') k

    qp = local_state(part, k)

    lap = laplacian(coefficient=1.0_dp, spacing=1.0_dp, measure=1.0_dp)
    call lap % apply(part, [qp], lq)
    call lq % get_real_vector(v)
    ok = .true.
    do i = 1, size(expect_l)
       ok = ok .and. (abs(v(i) - expect_l(i)) < 1.0d-12)
    end do
    call report(ok, &
         & "L on G" // tag // " traverses the PART's topology", nfail)

    call shifted % apply(part, [qp], aq)
    call aq % get_real_vector(v)
    ok = .true.
    do i = 1, size(expect_a)
       ok = ok .and. (abs(v(i) - expect_a(i)) < 1.0d-12)
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

    class(graph), intent(in)    :: part
    integer     , intent(in)    :: k, borrowed_global, watched_global
    real(dp)    , intent(in)    :: before, after
    integer     , intent(inout) :: nfail

    type(field)                     :: qp
    class(graph_field), allocatable :: aq
    real(dp), allocatable           :: v(:)
    real(dp), allocatable           :: state(:)
    character(len=1)                :: tag
    integer                         :: bseat, wseat

    write(tag,'(i1)') k

    bseat = seat_of_global(part, borrowed_global)
    wseat = seat_of_global(part, watched_global)

    ! Unperturbed: the owned answer stands where it should.
    qp = local_state(part, k)
    call shifted % apply(part, [qp], aq)
    call aq % get_real_vector(v)
    call report(abs(v(wseat) - before) < 1.0d-12, &
         & "G" // tag // ": the owned result at global " // &
         & achar(48 + watched_global) // " starts correct", nfail)

    ! Perturb ONLY the borrowed copy, by its global identity.
    call qp % get_real_vector(state)
    state(bseat) = state(bseat) + 10.0_dp
    call qp % set_real_vector(state)

    call shifted % apply(part, [qp], aq)
    call aq % get_real_vector(v)
    call report(abs(v(wseat) - after) < 1.0d-12, &
         & "G" // tag // ": +10 on the BORROWED copy of global " // &
         & achar(48 + borrowed_global) // " moves the OWNED result " // &
         & "at global " // achar(48 + watched_global) // " by -10", &
         & nfail)

    ! Restore, and the correct answer returns.
    qp = local_state(part, k)
    call shifted % apply(part, [qp], aq)
    call aq % get_real_vector(v)
    call report(abs(v(wseat) - before) < 1.0d-12, &
         & "G" // tag // ": restore the halo and equivalence returns", &
         & nfail)

  end subroutine perturb_and_watch

  !===================================================================!
  ! The Class-1 operation behind production minimization, on G.
  !===================================================================!

  subroutine check_global_solve(nfail)

    integer, intent(inout) :: nfail

    type(gmres)                     :: solver
    type(field)                     :: rhs
    class(member_set), allocatable  :: dom
    class(graph_field), allocatable :: sol
    real(dp), allocatable           :: gv(:), v(:)
    type(counted_set)               :: vs

    vs = g % vertex_set()

    call solver % attach(shifted, g, vs)
    solver % tolerance      = 1.0d-12
    solver % max_iterations = 200

    call solver % domain(g, dom)
    call report(dom % same_as(vs), &
         & "the solver answers on V(G)", nfail)
    call shifted % domain(g, dom)
    call report(dom % same_as(vs), &
         & "and so does its action: unknown and residual domains " // &
         & "coincide here, both being V(G)", nfail)

    call solver % constant(gv)
    call report(maxval(abs(gv)) < 1.0d-12, &
         & "the affine constant is zero: A is linear", nfail)

    rhs = field('b', vs)
    call rhs % set_real_vector(B_EXACT)
    call solver % apply(g, [rhs], sol)

    call sol % domain(dom)
    call sol % get_real_vector(v)
    call report(dom % same_as(vs), &
         & "the solution is a field on V(G)", nfail)
    call report(by_member(v, vs, Q_EXACT), &
         & "and A q = b solves to q* = [1,2,4,7,11,16], by member", &
         & nfail)

  end subroutine check_global_solve

  !===================================================================!
  ! THE CONDUIT TEST, behavioural. One operation type storing no
  ! graph; two solvers over two different six-vertex five-edge
  ! topologies; one probe. If the host were scenery the two matvecs
  ! would agree. They do not.
  !===================================================================!

  subroutine check_host_is_load_bearing(nfail)

    integer, intent(inout) :: nfail

    type(gmres)           :: solver_g, solver_alt
    type(counted_set)     :: vs, vs_alt
    real(dp), allocatable :: y(:), y_alt(:)

    ! Same counts, different shape: a star, not a chain.
    g_alt = stored_graph(NV, tails=[1,1,1,1,1], heads=[2,3,4,5,6])
    vs     = g % vertex_set()
    vs_alt = g_alt % vertex_set()

    call report(g_alt % num_vertices() .eq. g % num_vertices() .and. &
         &      g_alt % num_edges() .eq. g % num_edges() .and. &
         &      .not. vs_alt % same_as(vs), &
         & "G_alt is a star: same six vertices, same five edges, " // &
         & "different topology", nfail)

    call solver_g % attach(shifted, g, vs)
    call solver_alt % attach(shifted, g_alt, vs_alt)

    ! The SAME operation object, the SAME probe, two hosts.
    call solver_g % matvec(Q_EXACT, y)
    call solver_alt % matvec(Q_EXACT, y_alt)

    call report(maxval(abs(y - B_EXACT)) < 1.0d-12, &
         & "on the chain host the solver's own matvec gives b", nfail)
    call report(maxval(abs(y_alt - B_EXACT)) > 1.0_dp, &
         & "on the star host it does NOT - nothing changed but the " // &
         & "graph the solver carries", nfail)
    call report(maxval(abs(y_alt - y)) > 1.0_dp, &
         & "so the host reaches the mathematics: a minimizer that " // &
         & "reads no topology still CARRIES it to one that does", &
         & nfail)

  end subroutine check_host_is_load_bearing

  !===================================================================!
  ! Helpers.
  !===================================================================!

  ! q* transported onto this part - the overlap-complete local state.
  type(field) function local_state(part, k) result(qp)

    class(graph), intent(in) :: part
    integer     , intent(in) :: k

    type(partitioner)               :: p
    type(field)                     :: q
    class(graph_field), allocatable :: pd
    real(dp), allocatable           :: v(:)

    q = field('q star', g % vertex_set())
    call q % set_real_vector(Q_EXACT)

    p = partitioner(PARTITION_LINEAR, nparts=2, part=k)
    call p % partition_data(g, q, part, pd)

    call pd % get_real_vector(v)
    qp = field('local q', part % vertex_set())
    call qp % set_real_vector(v)

  end function local_state

  ! Apply A on a part and add its OWNED contribution into total.
  subroutine add_local_action(part, k, total)

    class(graph), intent(in)    :: part
    integer     , intent(in)    :: k
    real(dp)    , intent(inout) :: total(:)

    type(field)                     :: qp, aq_local
    class(graph_field), allocatable :: aq, fd
    class(member_set), allocatable  :: dom
    real(dp), allocatable           :: v(:)
    integer , allocatable           :: mem(:)
    integer                         :: i

    qp = local_state(part, k)
    call shifted % apply(part, [qp], aq)
    call aq % get_real_vector(v)

    aq_local = field('A local', part % vertex_set())
    call aq_local % set_real_vector(v)

    call a % assemble_data(part, aq_local, g, fd)
    call fd % domain(dom)
    call fd % get_real_vector(v)

    select type (dom)
    type is (subset_set)
       call dom % members(mem)
       do i = 1, size(mem)
          total(mem(i)) = total(mem(i)) + v(dom % local_index(mem(i)))
       end do
    class default
       do i = 1, size(total)
          total(i) = total(i) + v(i)
       end do
    end select

  end subroutine add_local_action

  ! The local seat holding this global member, or 0.
  integer function seat_of_global(part, gm)

    class(graph), intent(in) :: part
    integer     , intent(in) :: gm

    integer :: i

    seat_of_global = 0
    select type (part)
    type is (stored_graph)
       do i = 1, part % num_vertices()
          if (part % global_vertex_index(i) .eq. gm) seat_of_global = i
       end do
    end select

  end function seat_of_global

  ! Values compared through the domain's own map, never by position.
  logical function by_member(v, dom, expect)

    real(dp)         , intent(in) :: v(:)
    class(member_set), intent(in) :: dom
    real(dp)         , intent(in) :: expect(:)

    integer :: i, m

    by_member = .true.
    do i = 1, dom % size()
       m = dom % member(i)
       by_member = by_member .and. &
            & (abs(v(dom % local_index(m)) - expect(m)) < 1.0d-12)
    end do

  end function by_member

end program partitioned_pde_gate_b
