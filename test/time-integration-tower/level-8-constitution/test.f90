!=====================================================================!
! TIME INTEGRATION TOWER . LEVEL 8 . CONSTITUTION
!
! The level answers one question:
!
!      CAN EVERYTHING THE TOWER HAS EARNED BE CONSTITUTED INTO AN
!      ACTUAL MULTI-STEP TIME MARCH?
!
! This is where the production driver (operation_driver) is earned:
! the march is a bipartite digraph of steps and states the driver
! evaluates in topological order. Structural reach from Level 2,
! scheme coefficients from Level 6, implicit governance from
! Level 7, repeated along a causal chain.
!
!                    THREE FIVE-ELEMENT OBJECTS
!
! The specimen puts three different five-element things in one
! program, and the level's first duty is to keep them apart:
!
!      T                       the tower's instant CARRIER (L0)
!      V(H_context)            the operation host's vertices
!      the states' projection  the march digraph's CONTROL CHAIN
!
! plus the two-member state carrier Q, which is none of them. That
! all three have five elements is a coincidence of THIS specimen,
! and the assertions below refuse to rely on it: no two are the
! same carrier, by identity.
!
!      H_context is the compatibility conduit the operation
!      contract requires. It is NOT the time graph.
!
!      The march digraph projected onto its states is the driver's
!      own control chain. It is NOT T either - it is a second
!      realization of the same one-step structure.
!
! The bridge is therefore EXTENSIONAL, never identity: the chain's
! incidence says what A1 says, step for step. Two parties who agree
! need not be the same party.
!
!                    TWO SPECIALIZATIONS, NOT DEFECTS
!
! The march digraph is generated from nsteps rather than read from
! G_time, and the step stores ONE scalar h rather than the field
! h : E -> reals. For this specimen both are EXACT
! specializations - the time graph is a simple chain and h is
! uniform - so no defect is established here. They are recorded as
! frontier, for clients that would supply a nonuniform or
! nonlinear time structure. One tower cannot decide that.
!
!                    THE IMPLICIT GOVERNOR IS NEWTON
!
! Not bare GMRES. Each implicit step drives the FULL stated residual
! to zero via newton % solve(zeros, q, ...), and a bare GMRES matvec
! subtracts the affine constant - which would solve a different
! question. Newton reaches the difference linearization on its own.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program time_level_8

  use iso_fortran_env       , only : dp => REAL64
  use time_assert           , only : report, verdict
  use time_assert           , only : NQ, NT, NE, NSTEPS, TOL, TOL_MARCH
  use time_assert           , only : T0, T4, H_STEP, Q0
  use time_assert           , only : FE_TRAJECTORY, BE_TRAJECTORY, &
       &                             BDF2_TRAJECTORY
  use graph_fractal        , only : graph
  use map_set        , only : set_map
  use relation_binary , only : csr_relation
  use field_calculus  , only : field
  use view_directed_stored           , only : stored_directed_graph
  use view_read_write , only : bipartite_digraph, SECOND_PART
  use field_stored     , only : stored_field
  use time_carriers_fixture , only : time_carriers
  use time_relations_fixture, only : tail_relation, head_relation
  use time_algebra_fixture  , only : derive_one_step_reach
  use time_fields_fixture   , only : step_sizes, state_field
  use triangular_decay_fixture, only : triangular_decay
  use temporal_step_fixture , only : temporal_step, backward_euler
  use temporal_march_fixture, only : marched, march_incidence, &
       &                             MARCH_FORWARD, MARCH_BACKWARD, MARCH_BDF2

  implicit none

  type(graph)          :: q, t, e
  type(set_map)          :: sets
  type(csr_relation), target :: tail, head
  type(csr_relation)         :: a1
  type(stored_directed_graph)         :: hcontext
  type(triangular_decay)     :: decay
  type(stored_field)         :: q_initial
  integer                    :: nfail

  ! The implicit solve's tolerances: newton's and its inner gmres's.
  real(dp), parameter :: NEWTON_TOLERANCE = 1.0e-13_dp
  real(dp), parameter :: LINEAR_TOLERANCE = 1.0e-14_dp

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "time integration tower . level 8 . march"
  write(*,'(1x,a)') "============================================="

  call time_carriers(sets, q, t, e)
  tail = tail_relation(e, t, sets)
  head = head_relation(e, t, sets)
  a1   = derive_one_step_reach(tail, head, sets)

  ! The OPERATION HOST - the conduit, not the clock.
  hcontext = stored_directed_graph(NT, tails=[1,2,3,4], heads=[2,3,4,5])

  decay     = triangular_decay(q, NQ)
  q_initial = state_field(q)

  call check_control_chain_realizes_a1(nfail)
  call check_the_three_carriers_stay_apart(nfail)
  call check_scalar_step_specializes_the_field(nfail)
  call check_forward_march(nfail)
  call check_backward_march(nfail)
  call check_bdf2_march(nfail)

  call verdict(nfail, "level 8")

contains

  !===================================================================!
  ! THE first Rosetta bridge: the driver's control chain says what
  ! A1 says. Extensionally, step for step - and NOT by identity,
  ! which is not required and would be wrong to demand.
  !===================================================================!

  subroutine check_control_chain_realizes_a1(nfail)

    integer, intent(inout) :: nfail

    type(bipartite_digraph)     :: incidence
    type(stored_directed_graph) :: chain
    integer            :: i
    logical            :: ok

    incidence = march_incidence(NSTEPS, 1)
    chain     = incidence % projection(SECOND_PART)

    call report(chain % num_vertices() .eq. NT .and. &
         &      chain % num_edges() .eq. NE, &
         & "the march digraph's projection onto the states has five " // &
         & "instants and four steps, as T and E do", nfail)

    ! Each control step joins the instants A1 joins, read through the
    ! TIME carrier's own members rather than through the integers the
    ! chain happens to use.
    ok = .true.
    do i = 1, chain % num_edges()
       ok = ok .and. (chain % edge_tail(i) .eq. i)
       ok = ok .and. (chain % edge_head(i) .eq. i + 1)
       ok = ok .and. a1 % has([sets % member_of(t, i), sets % member_of(t, i + 1)])
    end do
    call report(ok, &
         & "and step i joins instant i to instant i+1 - exactly the " // &
         & "pairs A1 holds: THE DRIVER'S CONTROL CHAIN REALIZES THE " // &
         & "RELATIONAL TIME STRUCTURE, extensionally", nfail)

    call report(a1 % num_tuples() .eq. chain % num_edges(), &
         & "one control step per one-step reach, and no more", nfail)

  end subroutine check_control_chain_realizes_a1

  !===================================================================!
  ! THREE five-element objects and one two-element one, and no two
  ! of them are the same carrier. The coincidence of sizes is the
  ! specimen's, not the mathematics'.
  !===================================================================!

  subroutine check_the_three_carriers_stay_apart(nfail)

    integer, intent(inout) :: nfail

    type(bipartite_digraph)        :: incidence
    type(stored_directed_graph)             :: chain
    type(graph) :: cv, hv

    incidence = march_incidence(NSTEPS, 1)
    chain     = incidence % projection(SECOND_PART)
    cv = chain % vertex_set()
    hv = hcontext % vertex_set()

    call report(.not. cv % same_as(q), &
         & "the control chain's vertices are NOT Q: the clock is " // &
         & "not the state", nfail)

    call report(.not. cv % same_as(t), &
         & "nor are they T - two REALIZATIONS of one structure, and " // &
         & "agreement never made two parties one party", nfail)

    call report(.not. cv % same_as(hv), &
         & "nor V(H_context): THE CONTROL CHAIN IS NOT THE " // &
         & "OPERATION HOST, though both are five-element chains here", &
         & nfail)

    call report(.not. hv % same_as(q) .and. sets % num_members_of(q) .eq. NQ, &
         & "and V(H_context) is still not Q - two members against " // &
         & "five, as at Levels 6 and 7", nfail)

  end subroutine check_the_three_carriers_stay_apart

  !===================================================================!
  ! The step size: a field on E at Level 5, one scalar in the
  ! temporal step. For a UNIFORM h the scalar is an exact specialization,
  ! and this level says so in that direction - the scalar does not
  ! become a general variable-step model by being sufficient here.
  !===================================================================!

  subroutine check_scalar_step_specializes_the_field(nfail)

    integer, intent(inout) :: nfail

    type(temporal_step)   :: one
    type(stored_field)           :: h
    real(dp), allocatable :: hv(:)
    integer               :: i
    logical               :: ok

    one = backward_euler(decay, H_STEP)
    h   = step_sizes(e)
    call h % real_vector(hv)

    ok = .true.
    do i = 1, sets % num_members_of(e)
       ok = ok .and. &
            & (abs(hv(sets % index_in(e, sets % member_of(e, i))) - one % h) &
            &  .lt. TOL)
    end do
    call report(ok, &
         & "h(e) = the step's scalar h at every step: the scalar is an " // &
         & "EXACT SPECIALIZATION of the uniform step field", nfail)

    call report(size(hv) .eq. NE, &
         & "and the field still carries one value per step - four " // &
         & "numbers the scalar happens to agree with, not a variable-" // &
         & "step model the marcher implements", nfail)

  end subroutine check_scalar_step_specializes_the_field

  !===================================================================!
  ! THE explicit march: four forward-euler steps, on Q, with the
  ! five-vertex host carried alongside.
  !
  ! Against the production reviewed at Gate B this FAILED, and the
  ! failure is recorded verbatim in NUCLEUS-OBSERVATIONS.md TI-14:
  ! read_statement built the state on the HOST's vertex set and took
  ! its width as size(q) / num_vertices() = 2/5 = 0.
  !===================================================================!

  subroutine check_forward_march(nfail)

    integer, intent(inout) :: nfail

    real(dp), allocatable :: state(:)
    type(graph)           :: d
    integer               :: n
    logical               :: ok

    ! Every prefix from q0, so the whole trajectory is pinned rather
    ! than only its end.
    ok = .true.
    do n = 1, NSTEPS
       call march(MARCH_FORWARD, n, state, d)
       ok = ok .and. (maxval(abs(state - FE_TRAJECTORY(:, n))) .lt. TOL)
    end do
    call report(ok, &
         & "forward euler marches [2,0] -> [1,1] -> [1/2,1] -> " // &
         & "[1/4,3/4] -> [1/8,1/2], every prefix pinned", nfail)

    call march(MARCH_FORWARD, NSTEPS, state, d)
    call report(d % same_as(q) .and. &
         &      maxval(abs(state - FE_TRAJECTORY(:, NSTEPS))) .lt. TOL, &
         & "and the terminal state is [1/8, 1/2], on Q - THE DRIVEN " // &
         & "MARCH'S STATE DOMAIN IS INDEPENDENT OF ITS HOST", nfail)

  end subroutine check_forward_march

  !===================================================================!
  ! THE implicit march, through the full production composition:
  !
  !      driver -> newton -> difference linearization -> gmres
  !
  ! Newton, not bare GMRES: each step drives the whole residual to
  ! zero, and a bare GMRES matvec has already subtracted the affine
  ! part. Newton reaches the linearization by itself - this level
  ! never names that module, and the import gate refuses it.
  !===================================================================!

  subroutine check_backward_march(nfail)

    integer, intent(inout) :: nfail

    real(dp), allocatable :: state(:)
    type(graph)   :: d
    real(dp)      :: worst
    integer       :: n
    logical       :: ok

    ok = .true.
    worst = 0.0_dp
    do n = 1, NSTEPS
       call march(MARCH_BACKWARD, n, state, d)
       ok = ok .and. d % same_as(q)
       worst = max(worst, maxval(abs(state - BE_TRAJECTORY(:, n))))
       ok = ok .and. &
            & (maxval(abs(state - BE_TRAJECTORY(:, n))) .lt. TOL_MARCH)
    end do

    call report(ok, &
         & "backward euler marches [2,0] -> [4/3,4/9] -> [8/9,16/27] " // &
         & "-> [16/27,16/27] -> [32/81,128/243]", nfail)

    write(*,'(1x,a,es12.5)') "       worst backward-euler error : ", worst

  end subroutine check_backward_march

  !===================================================================!
  ! THE whole road, in one call: structural reach (L2), scheme
  ! coefficients (L6), implicit governance (L7), repeated along the
  ! causal chain - with the first step a backward one, as bdf-2 must
  ! and as production intends.
  !===================================================================!

  subroutine check_bdf2_march(nfail)

    integer, intent(inout) :: nfail

    real(dp), allocatable :: state(:)
    type(graph)   :: d
    real(dp)      :: worst
    integer       :: n
    logical       :: ok

    ok = .true.
    worst = 0.0_dp
    do n = 1, NSTEPS
       call march(MARCH_BDF2, n, state, d)
       ok = ok .and. d % same_as(q)
       worst = max(worst, maxval(abs(state - BDF2_TRAJECTORY(:, n))))
       ok = ok .and. &
            & (maxval(abs(state - BDF2_TRAJECTORY(:, n))) .lt. TOL_MARCH)
    end do

    call report(ok, &
         & "bdf-2 marches [2,0] -> [4/3,4/9] -> [5/6,47/72] -> " // &
         & "[1/2,2/3] -> [7/24,83/144], the first step backward as " // &
         & "the scheme requires", nfail)

    write(*,'(1x,a,es12.5)') "       worst bdf-2 error          : ", worst

    call report(maxval(abs(BDF2_TRAJECTORY(:, 1) - &
         &                 BE_TRAJECTORY(:, 1))) .lt. TOL, &
         & "and its first state IS the backward-euler one: a " // &
         & "two-step scheme cannot reach two steps back on its " // &
         & "first step, which is a structural fact before it is a " // &
         & "numerical one", nfail)

  end subroutine check_bdf2_march

  !===================================================================!
  ! n steps of the rule from q0 by the production driver, with
  ! newton over gmres at every implicit step; the state and its
  ! domain are read from the field the last step wrote.
  !===================================================================!

  subroutine march(rule, n, state, d)

    integer              , intent(in)  :: rule
    integer              , intent(in)  :: n
    real(dp), allocatable, intent(out) :: state(:)
    type(graph)          , intent(out) :: d

    class(field), allocatable :: final

    final = marched(rule, decay, H_STEP, n, hcontext, q_initial, &
         & NEWTON_TOLERANCE, LINEAR_TOLERANCE)
    call final % real_vector(state)
    d = final % domain()

  end subroutine march

end program time_level_8
