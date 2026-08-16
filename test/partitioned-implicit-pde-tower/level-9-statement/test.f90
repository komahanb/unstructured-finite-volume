!=====================================================================!
! PARTITIONED IMPLICIT PDE TOWER . LEVEL 9 . THE STATEMENT
!
! The sealing level. It answers one question: CAN A USER STATE THE
! COMPLETE PROBLEM -
!
!      solve  A q = b
!
! - while the implementation chooses the partitioned matrix-free
! realization? The user's sentence says nothing about parts. The
! machinery decomposes, refreshes overlaps, acts locally, assembles
! owned outputs and hands GMRES an operator that behaves exactly
! like the global one.
!
! Equivalence is required first at solver % matvec - the exact seat
! GMRES consumes, not merely between fixtures - on all four Level-8
! probes. Then both roads are solved INDEPENDENTLY from the same b,
! and must meet:
!
!      q_partitioned = q_global = q* = [1,2,4,7,11,16].
!
! The decomposition changed the road, not the answer.
!
! And the honest boundary, stated plainly because it would be easy
! to overclaim: this is a PARTITIONED MATRIX-FREE SOLVE, not a
! distributed solver. One process; one global trial vector;
! partition_data reads global state; the parts run sequentially in
! one address space; assembly happens in-process; inner products
! and norms reduce over global arrays; the Krylov state is global.
! What a genuinely distributed road would need is derived - not
! implemented - in NUCLEUS-OBSERVATIONS.md under PIP-8.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program partitioned_pde_level_9

  use iso_fortran_env  , only : dp => REAL64
  use partitioned_pde_assert, only : report, verdict
  use partitioned_pde_assert, only : NV, Q_EXACT, B_EXACT
  use fractal_graph        , only : set_graph => graph
  use graph_set_representation, only : counted_set_representation
  use graph_set_map        , only : set_map
  use graph_directed_view, only : directed_graph
  use graph_field_calculus, only : graph_field
  use class_graph      , only : directed_stored_graph
  use class_graph_field, only : field
  use class_graph_gmres, only : gmres
  use shifted_laplacian_fixture, only : shifted_laplacian
  use partitioned_shifted_laplacian_fixture, only : &
       & partitioned_shifted_laplacian

  implicit none

  real(dp), parameter :: MIXED(NV) = &
       & [3.0_dp, -1.0_dp, 4.0_dp, 1.0_dp, 5.0_dp, -9.0_dp]
  real(dp), parameter :: E3(NV) = &
       & [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
  real(dp), parameter :: E4(NV) = &
       & [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]

  type(directed_stored_graph)                  :: g
  type(shifted_laplacian)             :: direct
  type(partitioned_shifted_laplacian) :: composite
  type(gmres)                         :: solver_global, solver_part
  real(dp)                            :: q_part(NV), q_global(NV)
  integer                             :: nfail
  type(set_map)     :: sets

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "partitioned pde tower . level 9 . statement"
  write(*,'(1x,a)') "============================================="

  g = directed_stored_graph(NV, tails=[1,2,3,4,5], heads=[2,3,4,5,6])
  call sets % bind(g % vertex_set(), &
       & counted_set_representation(g % num_vertices()))
  call sets % bind(g % edge_set(), &
       & counted_set_representation(g % num_edges()))
  composite = partitioned_shifted_laplacian(g)

  call check_solver_matvec_equivalence(nfail)
  call check_the_two_solves(nfail)

  ! The tower's answer: the partitioned solution, as it stands.
  write(*,'(1x,a)', advance='no') "PARTITIONED_PDE_RESULT ="
  call write_field(q_part)

  call verdict(nfail, "level 9")

contains

  !===================================================================!
  ! Equivalence where GMRES actually consumes it.
  !===================================================================!

  subroutine check_solver_matvec_equivalence(nfail)

    integer, intent(inout) :: nfail

    type(set_graph)     :: vs
    real(dp), allocatable :: gv(:)
    integer         :: n_vs

    vs = g % vertex_set()
    n_vs = g % num_vertices()

    call solver_global % attach(direct, g, vs, n_vs)
    call solver_part % attach(composite, g, vs, n_vs)
    solver_global % tolerance      = 1.0d-12
    solver_global % max_iterations = 200
    solver_part % tolerance        = 1.0d-12
    solver_part % max_iterations   = 200

    call solver_global % constant(gv)
    call report(maxval(abs(gv)) < 1.0d-12, &
         & "the global action's affine constant is zero", nfail)
    call solver_part % constant(gv)
    call report(maxval(abs(gv)) < 1.0d-12, &
         & "and so is the partitioned action's: both are linear", nfail)

    call one_matvec(Q_EXACT, "q*", nfail)
    call one_matvec(MIXED  , "the mixed-sign probe", nfail)
    call one_matvec(E3     , "the interface basis e3", nfail)
    call one_matvec(E4     , "the interface basis e4", nfail)

  end subroutine check_solver_matvec_equivalence

  subroutine one_matvec(v, tag, nfail)

    real(dp)        , intent(in)    :: v(:)
    character(len=*), intent(in)    :: tag
    integer         , intent(inout) :: nfail

    real(dp), allocatable :: yg(:), yp(:)

    call solver_global % matvec(v, yg)
    call solver_part % matvec(v, yp)

    call report(maxval(abs(yg - yp)) < 1.0d-12, &
         & "solver matvec agrees on " // tag // ": the composite is " // &
         & "the operator GMRES sees", nfail)

  end subroutine one_matvec

  !===================================================================!
  ! The statement itself: one sentence, two roads, one answer.
  !===================================================================!

  subroutine check_the_two_solves(nfail)

    integer, intent(inout) :: nfail

    type(field)                     :: rhs
    type(set_graph)               :: vs
    class(graph_field), allocatable :: sol
    type(set_graph)  :: dom
    real(dp), allocatable           :: v(:)
    integer         :: n_vs

    vs = g % vertex_set()
    n_vs = g % num_vertices()
    rhs = field('b', vs, n_vs)
    call rhs % set_real_vector(B_EXACT)

    ! The partitioned road - the tower's own.
    call solver_part % apply(g, [rhs], sol)
    dom = sol % domain()
    call sol % get_real_vector(v)
    q_part = v(1:NV)

    call report(dom % same_as(vs), &
         & "the partitioned solution is a field on V(G)", nfail)
    call report(by_member(sets, q_part, vs, Q_EXACT), &
         & "and it is q* = [1,2,4,7,11,16], by global member", nfail)

    ! The global baseline, run independently.
    call solver_global % apply(g, [rhs], sol)
    call sol % get_real_vector(v)
    q_global = v(1:NV)

    call report(by_member(sets, q_global, vs, Q_EXACT), &
         & "the global solve independently reaches q* as well", nfail)
    call report(maxval(abs(q_part - q_global)) < 1.0d-9, &
         & "q_partitioned = q_global = q*: the decomposition changed " // &
         & "the road, not the answer", nfail)

  end subroutine check_the_two_solves

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
            & (abs(v(sets % index_in(dom, m)) - expect(m)) < 1.0d-9)
    end do

  end function by_member

  ! One real per global vertex, in global enumeration order,
  ! unrounded - the honest image of the arithmetic.
  subroutine write_field(v)

    real(dp), intent(in) :: v(:)

    integer :: i

    do i = 1, size(v)
       write(*,'(es24.16)', advance='no') v(i)
    end do
    write(*,'(a)') ""

  end subroutine write_field

end program partitioned_pde_level_9
