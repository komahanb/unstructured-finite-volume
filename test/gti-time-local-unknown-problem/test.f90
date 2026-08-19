!=====================================================================!
! The time-local unknown problem suite: one sample is declared
! unknown, a trial q is injected, and residual(q_trial) answers
! through the existing time-local seat -
!
!      trial -> injected samples -> motif rows -> U -> R.
!
! The injection laws are proven directly: the unknown q seat is
! replaced, every other sample, every sample time, and every
! occupied non-q seat is preserved, and the caller's samples are
! never mutated. The residual laws are proven through BDF1 with
! the newest sample unknown, with a changed trial, with a HISTORY
! sample unknown, and through a DIRK stage - plus the design and
! the clock witnessing that xi and t reach the form.
!
! Nothing is solved anywhere in this suite.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_gti_time_local_unknown_problem

  use iso_fortran_env      , only : dp => REAL64
  use fractal_graph        , only : graph
  use class_graph_field    , only : field
  use gti_value_buffers    , only : gti_value_buffer
  use gti_state_bundles    , only : GTI_STATE_Q, GTI_STATE_QDOT
  use gti_design_bundles   , only : gti_design_bundle
  use gti_time_local_schemes, only : gti_time_sample, gti_time_motif
  use gti_time_motif_builders, only : gti_time_motif_builder
  use gti_time_local_unknown_problems, only : gti_time_local_unknown_problem
  use gti_toy_forms        , only : toy_time_residual_form, toy_clock_form

  implicit none

  type(graph) :: states, designs
  type(field) :: q0_field, q1_field, base_field, zero_field
  type(field) :: qd9_field, qd7_field, xi_field, xi2_field

  type(gti_time_sample)   :: samples_a(2), samples_c(2), samples_d(2), samples_e(2)
  type(gti_time_sample), allocatable :: work(:)
  type(gti_design_bundle) :: design, design2

  type(gti_time_motif_builder)         :: builder
  type(gti_time_motif)                 :: bdf1, dirk
  type(gti_time_local_unknown_problem) :: problem
  type(toy_time_residual_form)         :: r_form
  type(toy_clock_form)                 :: clock
  type(gti_value_buffer)               :: trial_a, trial_b, trial_c, trial_d, out

  real(dp), allocatable :: rv(:), rw(:)
  integer :: nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "gti time-local unknown problem suite"
  write(*,'(1x,a)') "============================================="

  !-------------------------------------------------------------------!
  ! Fields, samples, motifs, trials. Placeholders sit where the
  ! unknown will be injected.
  !-------------------------------------------------------------------!

  call states  % declare()
  call designs % declare()

  q0_field = field('q0', states, 3)
  call q0_field % set_real_vector([1.0_dp, 2.0_dp, 4.0_dp])
  q1_field = field('q1', states, 3)
  call q1_field % set_real_vector([3.0_dp, 5.0_dp, 8.0_dp])
  base_field = field('base', states, 3)
  call base_field % set_real_vector([2.0_dp, 4.0_dp, 6.0_dp])
  zero_field = field('placeholder', states, 3)
  call zero_field % set_real_vector([0.0_dp, 0.0_dp, 0.0_dp])
  qd9_field = field('qdot9', states, 3)
  call qd9_field % set_real_vector([9.0_dp, 9.0_dp, 9.0_dp])
  qd7_field = field('qdot7', states, 3)
  call qd7_field % set_real_vector([7.0_dp, 7.0_dp, 7.0_dp])

  xi_field = field('xi', designs, 1)
  call xi_field % set_real_vector([0.5_dp])
  xi2_field = field('xi2', designs, 1)
  call xi2_field % set_real_vector([2.0_dp])

  allocate(samples_a(1) % state % component(1))
  samples_a(1) % state % component(1) % value = q0_field
  samples_a(1) % time = 0.0_dp
  allocate(samples_a(2) % state % component(1))
  samples_a(2) % state % component(1) % value = zero_field
  samples_a(2) % time = 0.5_dp

  allocate(samples_c(1) % state % component(1))
  samples_c(1) % state % component(1) % value = zero_field
  allocate(samples_c(2) % state % component(1))
  samples_c(2) % state % component(1) % value = q1_field

  allocate(samples_d(1) % state % component(1))
  samples_d(1) % state % component(1) % value = base_field
  allocate(samples_d(2) % state % component(1))
  samples_d(2) % state % component(1) % value = zero_field

  allocate(samples_e(1) % state % component(2))
  samples_e(1) % state % component(1) % value = q0_field
  samples_e(1) % state % component(2) % value = qd9_field
  allocate(samples_e(2) % state % component(2))
  samples_e(2) % state % component(1) % value = zero_field
  samples_e(2) % state % component(2) % value = qd7_field

  allocate(design % component(1))
  design % component(1) % value = xi_field
  allocate(design2 % component(1))
  design2 % component(1) % value = xi2_field

  call builder % bdf_uniform(1, 0.5_dp, bdf1)
  call builder % dirk_stage(0.5_dp, 2.0_dp, dirk)

  call trial_a % set_real([3.0_dp, 5.0_dp, 8.0_dp])
  call trial_b % set_real([4.0_dp, 6.0_dp, 9.0_dp])
  call trial_c % set_real([1.0_dp, 2.0_dp, 4.0_dp])
  call trial_d % set_real([3.0_dp, 5.0_dp, 9.0_dp])

  r_form % nequations = 3

  !-------------------------------------------------------------------!
  ! The injection laws.
  !-------------------------------------------------------------------!

  call problem % inject_trial_q(samples_a, 2, trial_a, work)

  call work(2) % state % component(1 + GTI_STATE_Q) % value % get_real_vector(rv)
  call report(matches(rv, [3.0_dp, 5.0_dp, 8.0_dp]), &
       & "the trial replaces the unknown q seat", nfail)

  call work(1) % state % component(1 + GTI_STATE_Q) % value % get_real_vector(rv)
  call report(matches(rv, [1.0_dp, 2.0_dp, 4.0_dp]), &
       & "the other samples' q is preserved", nfail)

  call report(abs(work(1) % time - 0.0_dp) < 1.0e-14_dp .and. &
       &      abs(work(2) % time - 0.5_dp) < 1.0e-14_dp, &
       & "the sample times are preserved", nfail)

  call problem % inject_trial_q(samples_e, 2, trial_a, work)
  call work(2) % state % component(1 + GTI_STATE_Q) % value % get_real_vector(rv)
  call work(2) % state % component(1 + GTI_STATE_QDOT) % value % get_real_vector(rw)
  call report(matches(rv, [3.0_dp, 5.0_dp, 8.0_dp]) .and. &
       & matches(rw, [7.0_dp, 7.0_dp, 7.0_dp]) .and. &
       & work(1) % state % has_component(GTI_STATE_QDOT), &
       & "occupied non-q seats survive the injection, unknown included", nfail)

  !-------------------------------------------------------------------!
  ! residual(q_trial) through BDF1 and DIRK.
  !-------------------------------------------------------------------!

  call problem % value(r_form, bdf1, samples_a, 2, trial_a, design, 0.5_dp, out)
  call out % get_real(rv)
  call report(matches(rv, [7.5_dp, 11.5_dp, 16.5_dp]), &
       & "BDF1, newest unknown: R(q_trial) = [7.5, 11.5, 16.5]", nfail)

  call report(out % nentries == 3 .and. out % ncomp == 1, &
       & "the output shape is held by the existing evaluator", nfail)

  call problem % value(r_form, bdf1, samples_a, 2, trial_b, design, 0.5_dp, out)
  call out % get_real(rv)
  call report(matches(rv, [10.5_dp, 14.5_dp, 19.5_dp]), &
       & "a changed trial changes the residual: [10.5, 14.5, 19.5]", nfail)

  call problem % value(r_form, bdf1, samples_c, 1, trial_c, design, 0.5_dp, out)
  call out % get_real(rv)
  call report(matches(rv, [7.5_dp, 11.5_dp, 16.5_dp]), &
       & "BDF1, history unknown: the same residual by symmetry of the data", nfail)

  call problem % value(r_form, dirk, samples_d, 2, trial_d, design, 2.0_dp, out)
  call out % get_real(rv)
  call report(matches(rv, [4.5_dp, 6.5_dp, 12.5_dp]), &
       & "DIRK stage unknown: R(q_trial) = [4.5, 6.5, 12.5]", nfail)

  !-------------------------------------------------------------------!
  ! The design and the time truly reach the form.
  !-------------------------------------------------------------------!

  call problem % value(r_form, bdf1, samples_a, 2, trial_a, design2, 0.5_dp, out)
  call out % get_real(rv)
  call report(matches(rv, [9.0_dp, 13.0_dp, 18.0_dp]), &
       & "a different design shifts the residual: xi = 2 adds 1.5", nfail)

  call problem % value(clock, bdf1, samples_a, 2, trial_a, design, 2.5_dp, out)
  call out % get_real(rv)
  call report(matches(rv, [2.5_dp]), &
       & "the clock form reads the requested time: R = [t] = [2.5]", nfail)

  !-------------------------------------------------------------------!
  ! The caller's samples were never mutated by any of the above.
  !-------------------------------------------------------------------!

  call samples_a(2) % state % component(1 + GTI_STATE_Q) % value % get_real_vector(rv)
  call samples_a(1) % state % component(1 + GTI_STATE_Q) % value % get_real_vector(rw)
  call report(matches(rv, [0.0_dp, 0.0_dp, 0.0_dp]) .and. &
       & matches(rw, [1.0_dp, 2.0_dp, 4.0_dp]), &
       & "the input samples are never mutated: the placeholder still holds zero", nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all unknown problem checks passed"
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

end program test_gti_time_local_unknown_problem
