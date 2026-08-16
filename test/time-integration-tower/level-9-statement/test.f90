!=====================================================================!
! TIME INTEGRATION TOWER . LEVEL 9 . STATEMENT
!
! The level answers one question: WHAT COMPLETE QUESTION DOES THE
! USER ASK. One initial-value problem, stated once, in the terms
! the tower spent nine levels earning:
!
!      Given
!
!          Q      = { x, y }              the state coordinates
!          q(0)   = [2, 0]                a FIELD on Q
!          qdot   = -S(q)
!          S(q)   = [ x, y - x ]
!          t0 = 0,  t4 = 2,  h = 1/2
!          scheme = bdf-2, started with one backward-euler step
!
!      compute
!
!          q(t4).
!
!                    THE BOUNDARY IS A FIELD
!
! The statement's two ends are FIELDS ON Q. The initial state
! arrives as one and its plain vector is fetched once; the answer
! is written back into one and its domain is checked. That is the
! field contract as the nucleus states it - fetch once, work in
! arrays, write back once - and the marcher's raw-array core is
! left exactly as it is. Nothing was refactored to make a public
! argument prettier.
!
!                    WHAT IS CHECKED BEFORE THE ANSWER IS SAID
!
! The endpoint is not asserted because 5 is the last integer. It is
! established through the tower's own objects: time(t4) = 2 read
! off the Level-5 coordinate field, four steps walked, and the
! control chain's terminal instant reached by following its
! incidence from t0 rather than by indexing.
!
! The marker carries the ACTUAL computed field, not the oracle. The
! oracle is what the test compares against; the marker is what the
! program found.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program time_level_9

  use iso_fortran_env       , only : dp => REAL64
  use time_assert           , only : report, verdict
  use time_assert           , only : NQ, NT, NSTEPS, TOL, TOL_MARCH
  use time_assert           , only : T0, T4, C_X, C_Y, H_STEP
  use time_assert           , only : TIME_COORD, BDF2_TRAJECTORY
  use fractal_graph        , only : set_graph => graph
  use graph_set_map        , only : set_map
  use class_graph           , only : directed_stored_graph
  use class_graph_field     , only : field
  use class_graph_gmres     , only : gmres
  use class_graph_newton    , only : newton
  use class_graph_marcher   , only : marcher, MARCH_BDF2
  use time_carriers_fixture , only : time_carriers
  use time_fields_fixture   , only : state_field, instant_coordinates
  use triangular_decay_fixture, only : triangular_decay

  implicit none

  type(set_graph)      :: q, t, e
  type(set_map)      :: sets
  type(directed_stored_graph)     :: hcontext
  type(triangular_decay) :: decay
  type(marcher)          :: clock
  type(field)            :: q_initial, q_final, tcoord
  real(dp), allocatable  :: state(:)
  integer                :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "time integration tower . level 9 . statement"
  write(*,'(1x,a)') "============================================="

  call time_carriers(sets, q, t, e)
  hcontext = directed_stored_graph(NT, tails=[1,2,3,4], heads=[2,3,4,5])
  decay    = triangular_decay(q, NQ)
  tcoord   = instant_coordinates(t)

  ! ---- the statement's near end: a FIELD on Q.
  q_initial = state_field(q)

  ! ---- the constitution, as Level 8 established it.
  clock % rule = MARCH_BDF2
  clock % step = H_STEP
  allocate(clock % inner, source=newton())
  select type (nw => clock % inner)
  type is (newton)
     nw % tolerance = 1.0e-13_dp
     allocate(nw % inner, source=gmres())
     nw % inner % tolerance = 1.0e-14_dp
  end select

  ! ---- fetch once, march, write back once.
  call q_initial % get_real_vector(state)
  call clock % march(decay, hcontext, state, NSTEPS)

  q_final = field('q(t4)', q, NQ, ncomp=1)
  call q_final % set_real_vector(state)

  call check_the_statement_s_two_ends(nfail)
  call check_the_endpoint_is_earned(nfail)
  call check_the_answer(nfail)

  call verdict(nfail, "level 9")

  call say_the_result()

contains

  !===================================================================!
  ! Both ends of the statement live on Q, and the question was asked
  ! and answered in the same terms it was posed in.
  !===================================================================!

  subroutine check_the_statement_s_two_ends(nfail)

    integer, intent(inout) :: nfail

    type(set_graph) :: d

    d = q_initial % domain()
    call report(d % same_as(q), &
         & "the statement's initial state is a FIELD ON Q", nfail)

    d = q_final % domain()
    call report(d % same_as(q), &
         & "and its answer is a FIELD ON Q - asked and answered in " // &
         & "the same terms", nfail)

    call report(q_final % num_entries() .eq. NQ .and. &
         &      q_final % num_components() .eq. 1, &
         & "two entries, one number each, as Q has two members", nfail)

  end subroutine check_the_statement_s_two_ends

  !===================================================================!
  ! The endpoint, established rather than assumed. Not "member 5 is
  ! the last integer" but: four steps were walked, the control chain
  ! reaches its terminal instant in four, and that instant's
  ! coordinate is 2.
  !===================================================================!

  subroutine check_the_endpoint_is_earned(nfail)

    integer, intent(inout) :: nfail

    type(directed_stored_graph)    :: chain
    real(dp), allocatable :: tv(:)
    integer               :: here, i
    logical               :: ok

    call tcoord % get_real_vector(tv)

    call report(abs(tv(sets % index_in(t, T0))) .lt. TOL .and. &
         &      abs(tv(sets % index_in(t, T4)) - 2.0_dp) .lt. TOL, &
         & "time(t0) = 0 and time(t4) = 2 - the span the statement " // &
         & "names, read from the Level-5 coordinate field", nfail)

    ! Walk the production control chain from its first instant,
    ! following incidence, and see where four steps land.
    call clock % instants(NSTEPS, chain)
    here = 1
    ok = .true.
    do i = 1, NSTEPS
       ok = ok .and. (chain % edge_tail(i) .eq. here)
       here = chain % edge_head(i)
    end do
    call report(ok .and. here .eq. chain % num_vertices(), &
         & "and four steps of the control chain reach its terminal " // &
         & "instant, followed through incidence rather than indexed", &
         & nfail)

    call report(here .eq. sets % index_in(t, T4) .and. &
         &      abs(tv(here) - TIME_COORD(NT)) .lt. TOL, &
         & "which corresponds to t4, at coordinate 2: FOUR STEPS OF " // &
         & "1/2 FROM 0", nfail)

  end subroutine check_the_endpoint_is_earned

  !===================================================================!
  ! THE ANSWER.
  !===================================================================!

  subroutine check_the_answer(nfail)

    integer, intent(inout) :: nfail

    real(dp), allocatable :: v(:)

    call q_final % get_real_vector(v)

    call report(maxval(abs(v - BDF2_TRAJECTORY(:, NSTEPS))) .lt. TOL_MARCH, &
         & "q(t4) = [7/24, 83/144] - the complete initial-value " // &
         & "problem, answered on the domain it was posed on", nfail)

    call report(abs(v(sets % index_in(q, C_X)) - 7.0_dp/24.0_dp) .lt. TOL_MARCH &
         & .and. abs(v(sets % index_in(q, C_Y)) - 83.0_dp/144.0_dp) .lt. TOL_MARCH, &
         & "read coordinate by coordinate at Q's local positions: " // &
         & "x = 7/24, y = 83/144", nfail)

  end subroutine check_the_answer

  !===================================================================!
  ! The result contract: one marker, two real tokens, in Q's
  ! declaration order - and they are the COMPUTED values, never the
  ! oracle. Whether they ARE q(t4) is the test above's business.
  !===================================================================!

  subroutine say_the_result()

    real(dp), allocatable :: v(:)

    call q_final % get_real_vector(v)

    write(*,'(1x,a,2(1x,es23.16))') "TIME_INTEGRATION_RESULT =", &
         & v(sets % index_in(q, C_X)), v(sets % index_in(q, C_Y))

  end subroutine say_the_result

end program time_level_9
