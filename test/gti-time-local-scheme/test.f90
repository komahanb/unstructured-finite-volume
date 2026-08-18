!=====================================================================!
! The time-local motif suite: coefficient rows build the state
! tuple of one evaluation point from time samples,
!
!      q^(m) = sum_i w^(m)_i q_i,
!
! and the residual form answers there through the form evaluator.
!
! Two motifs are exercised. The two-sample motif builds q and
! qdot; the three-sample motif builds q, qdot, and qddot; and the
! toy residual R = q + qdot + qddot + xi reproduces both by hand:
!
!      two samples:    R = [5.5, 8.5, 12.5]
!      three samples:  R = [10.5, 14.5, 19.5]
!
! The motif knows rows, never scheme names; nothing here is BDF or
! DIRK, and no time graph exists.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_gti_time_local_scheme

  use iso_fortran_env      , only : dp => REAL64
  use fractal_graph        , only : graph
  use class_graph_field    , only : field
  use gti_value_buffers    , only : gti_value_buffer
  use gti_state_bundles    , only : GTI_STATE_Q, GTI_STATE_QDOT, GTI_STATE_QDDOT
  use gti_design_bundles   , only : gti_design_bundle
  use gti_evaluation_points, only : gti_evaluation_point
  use gti_time_local_schemes, only : gti_time_sample, gti_time_motif, &
       & gti_time_local_residual_evaluator
  use gti_toy_forms        , only : toy_time_residual_form

  implicit none

  type(graph) :: states, designs
  type(field) :: q0_field, q1_field, q2_field, xi_field

  type(gti_time_sample)   :: samples2(2), samples3(3)
  type(gti_time_motif)    :: motif2, motif3
  type(gti_design_bundle) :: design

  type(gti_time_local_residual_evaluator) :: evaluator
  type(gti_evaluation_point)              :: point
  type(toy_time_residual_form)            :: r_form
  type(gti_value_buffer)                  :: out

  real(dp), allocatable :: rv(:)
  integer :: nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "gti time-local scheme motif suite"
  write(*,'(1x,a)') "============================================="

  !-------------------------------------------------------------------!
  ! Samples q0 = [1,2,4], q1 = [3,5,8], q2 = [6,9,13] on one state
  ! set, and the design xi = [0.5].
  !-------------------------------------------------------------------!

  call states  % declare()
  call designs % declare()

  q0_field = field('q0', states, 3)
  call q0_field % set_real_vector([1.0_dp, 2.0_dp, 4.0_dp])
  q1_field = field('q1', states, 3)
  call q1_field % set_real_vector([3.0_dp, 5.0_dp, 8.0_dp])
  q2_field = field('q2', states, 3)
  call q2_field % set_real_vector([6.0_dp, 9.0_dp, 13.0_dp])

  xi_field = field('xi', designs, 1)
  call xi_field % set_real_vector([0.5_dp])

  allocate(samples2(1) % state % component(1))
  samples2(1) % state % component(1) % value = q0_field
  samples2(1) % time = 0.0_dp
  allocate(samples2(2) % state % component(1))
  samples2(2) % state % component(1) % value = q1_field
  samples2(2) % time = 1.0_dp

  allocate(samples3(1) % state % component(1))
  samples3(1) % state % component(1) % value = q0_field
  allocate(samples3(2) % state % component(1))
  samples3(2) % state % component(1) % value = q1_field
  allocate(samples3(3) % state % component(1))
  samples3(3) % state % component(1) % value = q2_field

  allocate(design % component(1))
  design % component(1) % value = xi_field

  r_form % nequations = 3

  !-------------------------------------------------------------------!
  ! The two motifs: coefficient rows, nothing more.
  !-------------------------------------------------------------------!

  allocate(motif2 % rule(2))
  motif2 % rule(1) % state_component = GTI_STATE_Q
  motif2 % rule(1) % weight = [0.0_dp, 1.0_dp]
  motif2 % rule(2) % state_component = GTI_STATE_QDOT
  motif2 % rule(2) % weight = [-1.0_dp, 1.0_dp]

  allocate(motif3 % rule(3))
  motif3 % rule(1) % state_component = GTI_STATE_Q
  motif3 % rule(1) % weight = [0.0_dp, 0.0_dp, 1.0_dp]
  motif3 % rule(2) % state_component = GTI_STATE_QDOT
  motif3 % rule(2) % weight = [0.0_dp, -1.0_dp, 1.0_dp]
  motif3 % rule(3) % state_component = GTI_STATE_QDDOT
  motif3 % rule(3) % weight = [1.0_dp, -2.0_dp, 1.0_dp]

  call report(motif2 % size() == 2 .and. motif3 % size() == 3, &
       & "a motif counts its coefficient rows", nfail)

  call report(motif2 % has_component(GTI_STATE_Q) .and. &
       &      motif3 % has_component(GTI_STATE_Q), &
       & "both motifs rule the component q", nfail)

  call report(motif2 % has_component(GTI_STATE_QDOT) .and. &
       & .not. motif2 % has_component(GTI_STATE_QDDOT) .and. &
       &      motif3 % has_component(GTI_STATE_QDDOT), &
       & "ruled components answer yes, unruled answer no", nfail)

  !-------------------------------------------------------------------!
  ! The two-sample point: q = q1, qdot = q1 - q0.
  !-------------------------------------------------------------------!

  call evaluator % build_point(motif2, samples2, design, 1.25_dp, point)

  call point % state % component(1 + GTI_STATE_Q) % value % get_real_vector(rv)
  call report(matches(rv, [3.0_dp, 5.0_dp, 8.0_dp]), &
       & "two samples build q = 0 q0 + 1 q1", nfail)

  call point % state % component(1 + GTI_STATE_QDOT) % value % get_real_vector(rv)
  call report(matches(rv, [2.0_dp, 3.0_dp, 4.0_dp]) .and. &
       & point % state % differential_degree() == 1 .and. &
       & .not. point % state % has_component(GTI_STATE_QDDOT), &
       & "two samples build qdot = q1 - q0, and no qddot seat", nfail)

  call evaluator % value(r_form, motif2, samples2, design, 1.25_dp, out)
  call out % get_real(rv)
  call report(matches(rv, [5.5_dp, 8.5_dp, 12.5_dp]), &
       & "two-sample residual: R = q + qdot + xi", nfail)

  !-------------------------------------------------------------------!
  ! The three-sample point: q = q2, qdot = q2 - q1,
  ! qddot = q0 - 2 q1 + q2.
  !-------------------------------------------------------------------!

  call evaluator % build_point(motif3, samples3, design, 2.5_dp, point)

  call point % state % component(1 + GTI_STATE_Q) % value % get_real_vector(rv)
  call report(matches(rv, [6.0_dp, 9.0_dp, 13.0_dp]), &
       & "three samples build q = q2", nfail)

  call point % state % component(1 + GTI_STATE_QDOT) % value % get_real_vector(rv)
  call report(matches(rv, [3.0_dp, 4.0_dp, 5.0_dp]), &
       & "three samples build qdot = q2 - q1", nfail)

  call point % state % component(1 + GTI_STATE_QDDOT) % value % get_real_vector(rv)
  call report(matches(rv, [1.0_dp, 1.0_dp, 1.0_dp]), &
       & "three samples build qddot = q0 - 2 q1 + q2", nfail)

  call evaluator % value(r_form, motif3, samples3, design, 2.5_dp, out)
  call out % get_real(rv)
  call report(matches(rv, [10.5_dp, 14.5_dp, 19.5_dp]), &
       & "three-sample residual: R = q + qdot + qddot + xi", nfail)

  !-------------------------------------------------------------------!
  ! The built point carries the design and the requested time.
  !-------------------------------------------------------------------!

  call point % design % component(1) % value % get_real_vector(rv)
  call report(point % design % has_entries() .and. matches(rv, [0.5_dp]), &
       & "the design rides into the built point", nfail)

  call report(abs(point % time - 2.5_dp) < 1.0e-14_dp .and. &
       & point % window_id == 0 .and. point % step_id == 0 .and. &
       & point % stage_id == 0, &
       & "the point carries the requested time; the graph ids stay zero", nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all time-local motif checks passed"
  else
     error stop
  end if

contains

  pure function matches(values, expected) result(ok)

    real(dp), intent(in) :: values(:), expected(:)
    logical :: ok

    ok = size(values) == size(expected)
    if (ok) ok = all(abs(values - expected) < 1.0e-14_dp)

  end function matches

  subroutine report(ok, label, nfail)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: label
    integer, intent(inout) :: nfail
    if (ok) then
       write(*,'(1x,a,a)') "PASS : ", label
    else
       write(*,'(1x,a,a)') "FAIL : ", label
       nfail = nfail + 1
    end if
  end subroutine report

end program test_gti_time_local_scheme
