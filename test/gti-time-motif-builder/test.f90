!=====================================================================!
! The named motif builder suite: each named scheme formula leaves
! the builder as plain coefficient rows, and the existing
! time-local seat applies them unchanged. Both faces are checked -
! the raw rows against the formula tables, and the applied
! residual R = q + qdot + xi through
!
!      builder -> gti_time_motif -> samples -> U -> R.
!
! The four builds of this phase:
!
!      BDF1   h = 1/2          qdot = [-2, 2]
!      BDF2   h = 1            qdot = [0.5, -2, 1.5]
!      DIRK   h = 2, g = 1/2   qdot = [-1, 1]
!      ABM2   h = 2, b0 = 1/2  qdot = [-1, 1]
!
! DIRK and ABM2 emit the same rows here by arithmetic accident
! (h g = h b0 = 1), which is itself the point: a motif is rows,
! and the name that chose them is already forgotten.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_gti_time_motif_builder

  use iso_fortran_env      , only : dp => REAL64
  use fractal_graph        , only : graph
  use class_graph_field    , only : field
  use gti_value_buffers    , only : gti_value_buffer
  use gti_state_bundles    , only : GTI_STATE_Q, GTI_STATE_QDOT, GTI_STATE_QDDOT
  use gti_design_bundles   , only : gti_design_bundle
  use gti_time_local_schemes, only : gti_time_sample, gti_time_motif, &
       & gti_time_local_residual_evaluator
  use gti_time_motif_builders, only : gti_time_motif_builder
  use gti_toy_forms        , only : toy_time_residual_form

  implicit none

  type(graph) :: states, designs
  type(field) :: q0_field, q1_field, q2_field, base_field, stage_field, xi_field

  type(gti_time_sample)   :: bdf1_samples(2), bdf2_samples(3), pair_samples(2)
  type(gti_design_bundle) :: design

  type(gti_time_motif_builder)            :: builder
  type(gti_time_motif)                    :: bdf1, bdf2, dirk, abm2
  type(gti_time_motif)                    :: vbdf1, vbdf2u, vbdf2, ubdf2
  type(gti_time_motif)                    :: vabm, uabm, ft
  type(gti_time_local_residual_evaluator) :: evaluator
  type(toy_time_residual_form)            :: r_form
  type(gti_value_buffer)                  :: out

  real(dp), allocatable :: rv(:), w(:)
  integer :: nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "gti named time motif builder suite"
  write(*,'(1x,a)') "============================================="

  !-------------------------------------------------------------------!
  ! Fields and samples. BDF walks history [q0, q1(, q2)]; DIRK and
  ! ABM walk [base, newest], the base assembled outside the
  ! builder's sight.
  !-------------------------------------------------------------------!

  call states  % declare()
  call designs % declare()

  q0_field = field('q0', states, 3)
  call q0_field % set_real_vector([1.0_dp, 2.0_dp, 4.0_dp])
  q1_field = field('q1', states, 3)
  call q1_field % set_real_vector([3.0_dp, 5.0_dp, 8.0_dp])
  q2_field = field('q2', states, 3)
  call q2_field % set_real_vector([6.0_dp, 9.0_dp, 13.0_dp])
  base_field = field('base', states, 3)
  call base_field % set_real_vector([2.0_dp, 4.0_dp, 6.0_dp])
  stage_field = field('stage', states, 3)
  call stage_field % set_real_vector([3.0_dp, 5.0_dp, 9.0_dp])

  xi_field = field('xi', designs, 1)
  call xi_field % set_real_vector([0.5_dp])

  allocate(bdf1_samples(1) % state % component(1))
  bdf1_samples(1) % state % component(1) % value = q0_field
  allocate(bdf1_samples(2) % state % component(1))
  bdf1_samples(2) % state % component(1) % value = q1_field

  allocate(bdf2_samples(1) % state % component(1))
  bdf2_samples(1) % state % component(1) % value = q0_field
  allocate(bdf2_samples(2) % state % component(1))
  bdf2_samples(2) % state % component(1) % value = q1_field
  allocate(bdf2_samples(3) % state % component(1))
  bdf2_samples(3) % state % component(1) % value = q2_field

  allocate(pair_samples(1) % state % component(1))
  pair_samples(1) % state % component(1) % value = base_field
  allocate(pair_samples(2) % state % component(1))
  pair_samples(2) % state % component(1) % value = stage_field

  allocate(design % component(1))
  design % component(1) % value = xi_field

  r_form % nequations = 3

  !-------------------------------------------------------------------!
  ! BDF1, h = 1/2: rows, then the applied residual.
  !-------------------------------------------------------------------!

  call builder % bdf_uniform(1, 0.5_dp, bdf1)

  call report(bdf1 % rule(1) % state_component == GTI_STATE_Q .and. &
       & matches(bdf1 % rule(1) % weight, [0.0_dp, 1.0_dp]), &
       & "BDF1 q row is [0, 1]", nfail)

  call report(bdf1 % rule(2) % state_component == GTI_STATE_QDOT .and. &
       & matches(bdf1 % rule(2) % weight, [-2.0_dp, 2.0_dp]), &
       & "BDF1 qdot row is [-1/h, 1/h] = [-2, 2]", nfail)

  call evaluator % value(r_form, bdf1, bdf1_samples, design, 0.5_dp, out)
  call out % get_real(rv)
  call report(matches(rv, [7.5_dp, 11.5_dp, 16.5_dp]), &
       & "BDF1 applied: R = q + qdot + xi = [7.5, 11.5, 16.5]", nfail)

  !-------------------------------------------------------------------!
  ! BDF2, h = 1: rows, then the applied residual.
  !-------------------------------------------------------------------!

  call builder % bdf_uniform(2, 1.0_dp, bdf2)

  call report(bdf2 % rule(1) % state_component == GTI_STATE_Q .and. &
       & matches(bdf2 % rule(1) % weight, [0.0_dp, 0.0_dp, 1.0_dp]), &
       & "BDF2 q row is [0, 0, 1]", nfail)

  call report(bdf2 % rule(2) % state_component == GTI_STATE_QDOT .and. &
       & matches(bdf2 % rule(2) % weight, [0.5_dp, -2.0_dp, 1.5_dp]), &
       & "BDF2 qdot row is [1/2h, -2/h, 3/2h] = [0.5, -2, 1.5]", nfail)

  call evaluator % value(r_form, bdf2, bdf2_samples, design, 1.0_dp, out)
  call out % get_real(rv)
  call report(matches(rv, [10.0_dp, 14.0_dp, 19.0_dp]), &
       & "BDF2 applied: R = [10, 14, 19]", nfail)

  !-------------------------------------------------------------------!
  ! DIRK diagonal stage, h = 2, gamma = 1/2.
  !-------------------------------------------------------------------!

  call builder % dirk_stage(0.5_dp, 2.0_dp, dirk)

  call report(dirk % rule(1) % state_component == GTI_STATE_Q .and. &
       & matches(dirk % rule(1) % weight, [0.0_dp, 1.0_dp]), &
       & "DIRK q row is [0, 1]", nfail)

  call report(dirk % rule(2) % state_component == GTI_STATE_QDOT .and. &
       & matches(dirk % rule(2) % weight, [-1.0_dp, 1.0_dp]), &
       & "DIRK qdot row is [-1/hg, 1/hg] = [-1, 1]", nfail)

  call evaluator % value(r_form, dirk, pair_samples, design, 2.0_dp, out)
  call out % get_real(rv)
  call report(matches(rv, [4.5_dp, 6.5_dp, 12.5_dp]), &
       & "DIRK applied: R = [4.5, 6.5, 12.5]", nfail)

  !-------------------------------------------------------------------!
  ! ABM implicit corrector, order 2, h = 2, beta0 = 1/2.
  !-------------------------------------------------------------------!

  call builder % abm_corrector(2, 2.0_dp, abm2)

  call report(abm2 % rule(1) % state_component == GTI_STATE_Q .and. &
       & matches(abm2 % rule(1) % weight, [0.0_dp, 1.0_dp]), &
       & "ABM2 q row is [0, 1]", nfail)

  call report(abm2 % rule(2) % state_component == GTI_STATE_QDOT .and. &
       & matches(abm2 % rule(2) % weight, [-1.0_dp, 1.0_dp]), &
       & "ABM2 qdot row is [-2/h, 2/h] = [-1, 1]", nfail)

  call evaluator % value(r_form, abm2, pair_samples, design, 2.0_dp, out)
  call out % get_real(rv)
  call report(matches(rv, [4.5_dp, 6.5_dp, 12.5_dp]), &
       & "ABM2 applied: R = [4.5, 6.5, 12.5]", nfail)

  !-------------------------------------------------------------------!
  ! Every builder output of this phase is the same shape of motif:
  ! two rows, q and qdot ruled, qddot unruled.
  !-------------------------------------------------------------------!

  call report(bdf1 % has_component(GTI_STATE_Q) .and. &
       & bdf2 % has_component(GTI_STATE_Q) .and. &
       & dirk % has_component(GTI_STATE_Q) .and. &
       & abm2 % has_component(GTI_STATE_Q), &
       & "every builder output rules q", nfail)

  call report(bdf1 % has_component(GTI_STATE_QDOT) .and. &
       & bdf2 % has_component(GTI_STATE_QDOT) .and. &
       & dirk % has_component(GTI_STATE_QDOT) .and. &
       & abm2 % has_component(GTI_STATE_QDOT) .and. &
       & .not. bdf1 % has_component(GTI_STATE_QDDOT) .and. &
       & .not. abm2 % has_component(GTI_STATE_QDDOT), &
       & "every builder output rules qdot, and none rules qddot", nfail)

  call report(bdf1 % size() == 2 .and. bdf2 % size() == 2 .and. &
       & dirk % size() == 2 .and. abm2 % size() == 2, &
       & "every builder output carries exactly two rows", nfail)

  !-------------------------------------------------------------------!
  ! The variable-step family. Uniform steps must collapse onto the
  ! uniform builders exactly, and the nonuniform BDF2 row must be
  ! the derivative of the interpolating quadratic at t_n:
  !
  !      h0 = h_n = 2,  h1 = h_{n-1} = 3
  !      qdot = [ 2/15, -5/6, 7/10 ]
  !
  ! proven three ways below - the closed row, and polynomial
  ! exactness on constants, linears, and quadratics at the sample
  ! times t = [-5, -2, 0].
  !-------------------------------------------------------------------!

  call builder % bdf_variable(1, [0.5_dp], vbdf1)
  call report(matches(vbdf1 % rule(1) % weight, bdf1 % rule(1) % weight) .and. &
       & matches(vbdf1 % rule(2) % weight, bdf1 % rule(2) % weight), &
       & "bdf_variable(1, [h]) equals bdf_uniform(1, h)", nfail)

  call builder % bdf_uniform(2, 2.0_dp, ubdf2)
  call builder % bdf_variable(2, [2.0_dp, 2.0_dp], vbdf2u)
  call report(matches(vbdf2u % rule(1) % weight, ubdf2 % rule(1) % weight) .and. &
       & matches(vbdf2u % rule(2) % weight, ubdf2 % rule(2) % weight), &
       & "bdf_variable(2, [h, h]) collapses onto bdf_uniform(2, h)", nfail)

  call builder % bdf_variable(2, [2.0_dp, 3.0_dp], vbdf2)
  call report(matches(vbdf2 % rule(1) % weight, [0.0_dp, 0.0_dp, 1.0_dp]), &
       & "nonuniform BDF2 q row is [0, 0, 1]", nfail)

  w = vbdf2 % rule(2) % weight
  call report(matches(w, [2.0_dp / 15.0_dp, -5.0_dp / 6.0_dp, 0.7_dp]), &
       & "nonuniform BDF2 qdot row is [2/15, -5/6, 7/10]", nfail)

  call report(matches([sum(w)], [0.0_dp]), &
       & "nonuniform BDF2 differentiates constants to zero", nfail)

  call report(matches([w(1) * (-5.0_dp) + w(2) * (-2.0_dp) + w(3) * 0.0_dp], &
       & [1.0_dp]), &
       & "nonuniform BDF2 differentiates q = t to exactly one", nfail)

  call report(matches([w(1) * 25.0_dp + w(2) * 4.0_dp + w(3) * 0.0_dp], &
       & [0.0_dp]), &
       & "nonuniform BDF2 differentiates q = t^2 to zero at t_n", nfail)

  !-------------------------------------------------------------------!
  ! The ABM current-step wrapper: the SAME rows as the uniform
  ! corrector, at both orders - variable-step READY because only
  ! h_n enters; the Adams history lives in the externally
  ! assembled base, and no history-ratio coefficients exist here.
  !-------------------------------------------------------------------!

  call builder % abm_corrector(1, 2.0_dp, uabm)
  call builder % abm_corrector_variable(1, 2.0_dp, vabm)
  call report(matches(vabm % rule(1) % weight, uabm % rule(1) % weight) .and. &
       & matches(vabm % rule(2) % weight, uabm % rule(2) % weight), &
       & "abm_corrector_variable(1, h) equals abm_corrector(1, h)", nfail)

  call builder % abm_corrector(2, 2.0_dp, uabm)
  call builder % abm_corrector_variable(2, 2.0_dp, vabm)
  call report(matches(vabm % rule(1) % weight, uabm % rule(1) % weight) .and. &
       & matches(vabm % rule(2) % weight, uabm % rule(2) % weight), &
       & "abm_corrector_variable(2, h) equals abm_corrector(2, h)", nfail)

  !-------------------------------------------------------------------!
  ! Times in, rows out: bdf_from_times mints the same rows its
  ! steps would.
  !-------------------------------------------------------------------!

  call builder % bdf_from_times([0.0_dp, 2.0_dp], ft)
  call builder % bdf_variable(1, [2.0_dp], vbdf1)
  call report(matches(ft % rule(2) % weight, vbdf1 % rule(2) % weight), &
       & "bdf_from_times([0, 2]) equals bdf_variable(1, [2])", nfail)

  call builder % bdf_from_times([0.0_dp, 3.0_dp, 5.0_dp], ft)
  call report(matches(ft % rule(2) % weight, vbdf2 % rule(2) % weight), &
       & "bdf_from_times([0, 3, 5]) equals bdf_variable(2, [2, 3])", nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all named motif builder checks passed"
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

end program test_gti_time_motif_builder
