!=====================================================================!
! PARTITIONED IMPLICIT PDE TOWER . LEVEL 7 . MINIMIZATION
!
! The level answers two questions. Can the discrete global problem
! be solved by ordinary minimization - and what does the graph the
! minimizer carries actually DO?
!
! The first is quickly settled: production GMRES, attached to the
! Level-6 shifted Laplacian on G, solves A q = b to q*.
!
! The second is the one four earlier towers could not answer. They
! all reported the minimizer's graph host as scenery, correctly,
! about their own actions - which had no topology to traverse.
! Here the attached action does. So: the SAME shifted_laplacian
! type, which stores no graph at all, is attached to two solvers
! over two topologies with the same six vertices and the same five
! edges - the chain, and a star - and handed the same probe:
!
!      solver on the chain -> b
!      solver on the star  -> something else entirely
!
! Nothing changed but the host. The conclusion is precise, and the
! precision matters:
!
!      GMRES itself does NOT traverse topology.
!      It CARRIES the graph to an operation that does.
!
!      minimizer                graph as conduit / context
!      differential operation   graph as numerical topology operand
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program partitioned_pde_level_7

  use iso_fortran_env  , only : dp => REAL64
  use partitioned_pde_assert, only : report, verdict
  use partitioned_pde_assert, only : NV, Q_EXACT, B_EXACT
  use graph_fractal        , only : set_graph => graph
  use map_set_representation, only : counted_set_representation
  use map_set        , only : set_map
  use graph_directed_view, only : directed_graph
  use graph_field_calculus, only : graph_field
  use class_graph      , only : directed_stored_graph
  use class_graph_field, only : field
  use class_graph_gmres, only : gmres
  use shifted_laplacian_fixture, only : shifted_laplacian

  implicit none

  type(directed_stored_graph)      :: g, g_alt
  type(shifted_laplacian) :: shifted
  integer                 :: nfail
  type(set_map)     :: sets

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "partitioned pde tower . level 7 . minimization"
  write(*,'(1x,a)') "============================================="

  g = directed_stored_graph(NV, tails=[1,2,3,4,5], heads=[2,3,4,5,6])
  call sets % bind(g % vertex_set(), &
       & counted_set_representation(g % num_vertices()))
  call sets % bind(g % edge_set(), &
       & counted_set_representation(g % num_edges()))

  call check_global_solve(nfail)
  call check_host_is_load_bearing(nfail)

  call verdict(nfail, "level 7")

contains
  !===================================================================!
  ! The Class-1 operation behind production minimization, on G.
  !===================================================================!

  subroutine check_global_solve(nfail)

    integer, intent(inout) :: nfail

    type(gmres)                     :: solver
    type(field)                     :: rhs
    type(set_graph)  :: dom
    class(graph_field), allocatable :: sol
    real(dp), allocatable           :: gv(:), v(:)
    type(set_graph)               :: vs
    integer         :: n_dom
    integer         :: n_vs

    vs = g % vertex_set()
    n_vs = g % num_vertices()

    call solver % attach(shifted, g, vs, n_vs)
    solver % tolerance      = 1.0d-12
    solver % max_iterations = 200

    call solver % domain(g, dom, n_dom)
    call report(dom % same_as(vs), &
         & "the solver answers on V(G)", nfail)
    call shifted % domain(g, dom, n_dom)
    call report(dom % same_as(vs), &
         & "and so does its action: unknown and residual domains " // &
         & "coincide here, both being V(G)", nfail)

    call solver % constant(gv)
    call report(maxval(abs(gv)) < 1.0d-12, &
         & "the affine constant is zero: A is linear", nfail)

    rhs = field('b', vs, n_vs)
    call rhs % set_real_vector(B_EXACT)
    call solver % apply(g, [rhs], sol)

    dom = sol % domain()
    call sol % get_real_vector(v)
    call report(dom % same_as(vs), &
         & "the solution is a field on V(G)", nfail)
    call report(by_member(sets, v, vs, Q_EXACT), &
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
    type(set_graph)     :: vs, vs_alt
    real(dp), allocatable :: y(:), y_alt(:)
    integer         :: n_vs
    integer         :: n_vs_alt

    ! Same counts, different shape: a star, not a chain.
    g_alt = directed_stored_graph(NV, tails=[1,1,1,1,1], heads=[2,3,4,5,6])
    call sets % bind(g_alt % vertex_set(), &
         & counted_set_representation(g_alt % num_vertices()))
    call sets % bind(g_alt % edge_set(), &
         & counted_set_representation(g_alt % num_edges()))
    vs     = g % vertex_set()
    n_vs = g % num_vertices()
    vs_alt = g_alt % vertex_set()
    n_vs_alt = g_alt % num_vertices()

    call report(g_alt % num_vertices() .eq. g % num_vertices() .and. &
         &      g_alt % num_edges() .eq. g % num_edges() .and. &
         &      .not. vs_alt % same_as(vs), &
         & "G_alt is a star: same six vertices, same five edges, " // &
         & "different topology", nfail)

    call solver_g % attach(shifted, g, vs, n_vs)
    call solver_alt % attach(shifted, g_alt, vs_alt, n_vs_alt)

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
  ! Values compared through the domain's own map, never by position.
  logical function by_member(sets, v, dom, expect)

    type(set_map)  , intent(in) :: sets
    real(dp)       , intent(in) :: v(:)
    type(set_graph), intent(in) :: dom
    real(dp)       , intent(in) :: expect(:)

    integer :: i, m

    by_member = .true.
    do i = 1, sets % size_of(dom)
       m = sets % member_of(dom, i)
       by_member = by_member .and. &
            & (abs(v(sets % index_in(dom, m)) - expect(m)) < 1.0d-12)
    end do

  end function by_member
end program partitioned_pde_level_7
