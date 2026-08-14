!=====================================================================!
! TIME INTEGRATION TOWER . LEVEL 8 . CONSTITUTION
!
! The level answers one question:
!
!      CAN EVERYTHING THE TOWER HAS EARNED BE CONSTITUTED INTO AN
!      ACTUAL MULTI-STEP TIME MARCH?
!
! This is where production class_graph_marcher is finally earned,
! and where the whole road meets itself: structural reach from
! Level 2, scheme coefficients from Level 6, implicit governance
! from Level 7, repeated along a causal chain.
!
!                    THREE FIVE-ELEMENT OBJECTS
!
! The specimen puts three different five-element things in one
! program, and the level's first duty is to keep them apart:
!
!      T                       the tower's instant CARRIER (L0)
!      V(H_context)            the operation host's vertices
!      chain from instants(4)  the marcher's CONTROL CHAIN
!
! plus the two-member state carrier Q, which is none of them. That
! all three have five elements is a coincidence of THIS specimen,
! and the assertions below refuse to rely on it: no two are the
! same carrier, by identity.
!
!      H_context is the compatibility conduit the graph_operation
!      contract requires. It is NOT the time graph.
!
!      clock % instants(4) is the marcher's own control chain. It
!      is NOT T either - it is a second realization of the same
!      one-step structure.
!
! The bridge is therefore EXTENSIONAL, never identity: the chain's
! incidence says what A1 says, step for step. Two parties who agree
! need not be the same party.
!
!                    TWO SPECIALIZATIONS, NOT DEFECTS
!
! Production regenerates a linear chain from nsteps rather than
! consuming G_time, and carries ONE scalar step rather than the
! field h : E -> R. For this specimen both are EXACT
! specializations - the time graph is a simple chain and h is
! uniform - so no defect is established here. They are recorded as
! frontier, for clients that would supply a nonuniform or
! nonlinear time structure. One tower cannot decide that.
!
!                    THE IMPLICIT GOVERNOR IS NEWTON
!
! Not bare GMRES. The marcher drives the FULL attached residual to
! zero via inner % solve(zeros, q, ...), and a bare GMRES matvec
! subtracts the affine constant - which would solve a different
! question. Newton's semantics are the right ones, and Newton
! reaches difference_linearization on its own.
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
  use graph_carrier         , only : counted_set, member_set
  use graph_binary_relation , only : csr_relation
  use class_graph           , only : stored_graph
  use class_graph_field     , only : field
  use class_graph_gmres     , only : gmres
  use class_graph_newton    , only : newton
  use class_graph_marcher   , only : marcher, MARCH_FORWARD, &
       &                             MARCH_BACKWARD, MARCH_BDF2
  use time_carriers_fixture , only : time_carriers
  use time_relations_fixture, only : tail_relation, head_relation
  use time_algebra_fixture  , only : derive_one_step_reach
  use time_fields_fixture   , only : step_sizes
  use triangular_decay_fixture, only : triangular_decay

  implicit none

  type(counted_set)          :: q, t, e
  type(csr_relation), target :: tail, head
  type(csr_relation)         :: a1
  type(stored_graph)         :: hcontext
  type(triangular_decay)     :: decay
  integer                    :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "time integration tower . level 8 . march"
  write(*,'(1x,a)') "============================================="

  call time_carriers(q, t, e)
  tail = tail_relation(e, t)
  head = head_relation(e, t)
  a1   = derive_one_step_reach(tail, head)

  ! The OPERATION HOST - the conduit, not the clock.
  hcontext = stored_graph(NT, tails=[1,2,3,4], heads=[2,3,4,5])

  decay = triangular_decay(q)

  call check_control_chain_realizes_a1(nfail)
  call check_the_three_carriers_stay_apart(nfail)
  call check_scalar_step_specializes_the_field(nfail)
  call check_forward_march(nfail)
  call check_backward_march(nfail)
  call check_bdf2_march(nfail)

  call verdict(nfail, "level 8")

contains

  !===================================================================!
  ! THE first Rosetta bridge: the marcher's control chain says what
  ! A1 says. Extensionally, step for step - and NOT by identity,
  ! which is not required and would be wrong to demand.
  !===================================================================!

  subroutine check_control_chain_realizes_a1(nfail)

    integer, intent(inout) :: nfail

    type(marcher)      :: clock
    type(stored_graph) :: chain
    integer            :: i
    logical            :: ok

    call clock % instants(NSTEPS, chain)

    call report(chain % num_vertices() .eq. NT .and. &
         &      chain % num_edges() .eq. NE, &
         & "the marcher's control chain has five instants and four " // &
         & "steps, as T and E do", nfail)

    ! Each control step joins the instants A1 joins, read through the
    ! TIME carrier's own members rather than through the integers the
    ! chain happens to use.
    ok = .true.
    do i = 1, chain % num_edges()
       ok = ok .and. (chain % edge_tail(i) .eq. i)
       ok = ok .and. (chain % edge_head(i) .eq. i + 1)
       ok = ok .and. a1 % has([t % member(i), t % member(i + 1)])
    end do
    call report(ok, &
         & "and step i joins instant i to instant i+1 - exactly the " // &
         & "pairs A1 holds: PRODUCTION'S CONTROL CHAIN REALIZES THE " // &
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

    type(marcher)                  :: clock
    type(stored_graph)             :: chain
    class(member_set), allocatable :: cv, hv

    call clock % instants(NSTEPS, chain)
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

    call report(.not. hv % same_as(q) .and. q % size() .eq. NQ, &
         & "and V(H_context) is still not Q - two members against " // &
         & "five, as at Levels 6 and 7", nfail)

  end subroutine check_the_three_carriers_stay_apart

  !===================================================================!
  ! The step size: a field on E at Level 5, one scalar in the
  ! marcher. For a UNIFORM h the scalar is an exact specialization,
  ! and this level says so in that direction - the scalar does not
  ! become a general variable-step model by being sufficient here.
  !===================================================================!

  subroutine check_scalar_step_specializes_the_field(nfail)

    integer, intent(inout) :: nfail

    type(marcher)         :: clock
    type(field)           :: h
    real(dp), allocatable :: hv(:)
    integer               :: i
    logical               :: ok

    clock % step = H_STEP
    h = step_sizes(e)
    call h % get_real_vector(hv)

    ok = .true.
    do i = 1, e % size()
       ok = ok .and. &
            & (abs(hv(e % local_index(e % member(i))) - clock % step) &
            &  .lt. TOL)
    end do
    call report(ok, &
         & "h(e) = clock % step at every step: the scalar is an " // &
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

    type(marcher)         :: clock
    real(dp)              :: state(NQ)
    integer               :: n
    logical               :: ok

    clock % rule = MARCH_FORWARD
    clock % step = H_STEP

    ! Every prefix from q0, so the whole trajectory is pinned rather
    ! than only its end.
    ok = .true.
    do n = 1, NSTEPS
       state = Q0
       call clock % march(decay, hcontext, state, n)
       ok = ok .and. (maxval(abs(state - FE_TRAJECTORY(:, n))) .lt. TOL)
    end do
    call report(ok, &
         & "forward euler marches [2,0] -> [1,1] -> [1/2,1] -> " // &
         & "[1/4,3/4] -> [1/8,1/2], every prefix pinned", nfail)

    state = Q0
    call clock % march(decay, hcontext, state, NSTEPS)
    call report(maxval(abs(state - FE_TRAJECTORY(:, NSTEPS))) .lt. TOL, &
         & "and the terminal state is [1/8, 1/2] - THE MARCHER'S " // &
         & "STATE DOMAIN IS INDEPENDENT OF ITS HOST", nfail)

  end subroutine check_forward_march

  !===================================================================!
  ! THE implicit march, through the full production composition:
  !
  !      marcher -> newton -> difference_linearization -> gmres
  !
  ! Newton, not bare GMRES: the marcher drives the whole residual to
  ! zero, and a bare GMRES matvec has already subtracted the affine
  ! part. Newton reaches the linearization by itself - this level
  ! never names that module, and the import gate refuses it.
  !===================================================================!

  subroutine check_backward_march(nfail)

    integer, intent(inout) :: nfail

    type(marcher) :: clock
    real(dp)      :: state(NQ), worst
    integer       :: n
    logical       :: ok

    call govern(clock, MARCH_BACKWARD)

    ok = .true.
    worst = 0.0_dp
    do n = 1, NSTEPS
       state = Q0
       call clock % march(decay, hcontext, state, n)
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

    type(marcher) :: clock
    real(dp)      :: state(NQ), worst
    integer       :: n
    logical       :: ok

    call govern(clock, MARCH_BDF2)

    ok = .true.
    worst = 0.0_dp
    do n = 1, NSTEPS
       state = Q0
       call clock % march(decay, hcontext, state, n)
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
  ! The implicit governor: newton over gmres, as production
  ! composes it elsewhere. Nothing here names the linearization
  ! newton will build.
  !===================================================================!

  subroutine govern(clock, rule)

    type(marcher), intent(out) :: clock
    integer      , intent(in)  :: rule

    clock % rule = rule
    clock % step = H_STEP

    allocate(clock % inner, source=newton())
    select type (nw => clock % inner)
    type is (newton)
       nw % tolerance = 1.0e-13_dp
       allocate(nw % inner, source=gmres())
       nw % inner % tolerance = 1.0e-14_dp
    end select

  end subroutine govern

end program time_level_8
