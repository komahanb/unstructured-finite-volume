!=====================================================================!
! The time-local Newton driver suite: for one unknown sample,
!
!      residual(q_trial)
!        -> analytic local Jacobian action
!        -> dense local Newton step
!        -> converged q*
!
! proven on three linear motifs with known closed forms -
!
!      BDF1  R = 3 q - 2 q0 + xi        J = 3 I
!      DIRK  R = 2 q - q_base + xi      J = 2 I
!      BDF2  R = 2.5 q + 0.5 q0 - 2 q1 + xi,  J = 2.5 I
!
! - and one nonlinear scalar, R = q^2 - 4, where Newton genuinely
! iterates from q = 3 down to q = 2. Failure to converge is an
! answer, not an error, and the caller's samples are never touched.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_gti_time_local_newton_driver

  use iso_fortran_env      , only : dp => REAL64
  use fractal_graph        , only : graph
  use class_graph_field    , only : field
  use gti_value_buffers    , only : gti_value_buffer
  use gti_state_bundles    , only : GTI_STATE_Q
  use gti_design_bundles   , only : gti_design_bundle
  use gti_time_local_schemes, only : gti_time_sample, gti_time_motif
  use gti_time_motif_builders, only : gti_time_motif_builder
  use gti_time_local_newton_drivers, only : gti_time_local_newton_driver, &
       & gti_time_local_newton_options, gti_time_local_newton_result
  use gti_toy_forms        , only : toy_time_residual_form, toy_square_form

  implicit none

  type(graph) :: states, scalars, designs
  type(field) :: q0_field, q1_field, base_field, zero_field, zero1_field, xi_field

  type(gti_time_sample)   :: bdf1_samples(2), dirk_samples(2), bdf2_samples(3)
  type(gti_time_sample)   :: square_samples(1)
  type(gti_design_bundle) :: design, no_design

  type(gti_time_motif_builder) :: builder
  type(gti_time_motif)         :: bdf1, dirk, bdf2, identity_motif

  type(gti_time_local_newton_driver)  :: driver
  type(gti_time_local_newton_options) :: options, starved
  type(gti_time_local_newton_result)  :: result

  type(toy_time_residual_form) :: r_form
  type(toy_square_form)        :: square

  type(gti_value_buffer) :: trial, dq, q_start, scalar_start, out

  real(dp), allocatable :: rv(:), jac(:,:)
  real(dp) :: root_bdf1(3), root_dirk(3), root_bdf2(3)
  integer :: nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "gti time-local newton driver suite"
  write(*,'(1x,a)') "============================================="

  !-------------------------------------------------------------------!
  ! Fields, samples, motifs, and the closed-form roots.
  !-------------------------------------------------------------------!

  call states  % declare()
  call scalars % declare()
  call designs % declare()

  q0_field = field('q0', states, 3)
  call q0_field % set_real_vector([1.0_dp, 2.0_dp, 4.0_dp])
  q1_field = field('q1', states, 3)
  call q1_field % set_real_vector([3.0_dp, 5.0_dp, 8.0_dp])
  base_field = field('base', states, 3)
  call base_field % set_real_vector([2.0_dp, 4.0_dp, 6.0_dp])
  zero_field = field('placeholder', states, 3)
  call zero_field % set_real_vector([0.0_dp, 0.0_dp, 0.0_dp])
  zero1_field = field('scalar placeholder', scalars, 1)
  call zero1_field % set_real_vector([0.0_dp])

  xi_field = field('xi', designs, 1)
  call xi_field % set_real_vector([0.5_dp])

  allocate(bdf1_samples(1) % state % component(1))
  bdf1_samples(1) % state % component(1) % value = q0_field
  allocate(bdf1_samples(2) % state % component(1))
  bdf1_samples(2) % state % component(1) % value = zero_field

  allocate(dirk_samples(1) % state % component(1))
  dirk_samples(1) % state % component(1) % value = base_field
  allocate(dirk_samples(2) % state % component(1))
  dirk_samples(2) % state % component(1) % value = zero_field

  allocate(bdf2_samples(1) % state % component(1))
  bdf2_samples(1) % state % component(1) % value = q0_field
  allocate(bdf2_samples(2) % state % component(1))
  bdf2_samples(2) % state % component(1) % value = q1_field
  allocate(bdf2_samples(3) % state % component(1))
  bdf2_samples(3) % state % component(1) % value = zero_field

  allocate(square_samples(1) % state % component(1))
  square_samples(1) % state % component(1) % value = zero1_field

  allocate(design % component(1))
  design % component(1) % value = xi_field

  call builder % bdf_uniform(1, 0.5_dp, bdf1)
  call builder % dirk_stage(0.5_dp, 2.0_dp, dirk)
  call builder % bdf_uniform(2, 1.0_dp, bdf2)

  allocate(identity_motif % rule(1))
  identity_motif % rule(1) % state_component = GTI_STATE_Q
  identity_motif % rule(1) % weight = [1.0_dp]

  root_bdf1 = (2.0_dp * [1.0_dp, 2.0_dp, 4.0_dp] - 0.5_dp) / 3.0_dp
  root_dirk = ([2.0_dp, 4.0_dp, 6.0_dp] - 0.5_dp) / 2.0_dp
  root_bdf2 = (2.0_dp * [3.0_dp, 5.0_dp, 8.0_dp] &
       &       - 0.5_dp * [1.0_dp, 2.0_dp, 4.0_dp] - 0.5_dp) / 2.5_dp

  r_form % nequations = 3

  !-------------------------------------------------------------------!
  ! The residual and the exact Jacobian action, BDF1.
  !-------------------------------------------------------------------!

  call trial % set_real([3.0_dp, 5.0_dp, 8.0_dp])
  call driver % residual(r_form, bdf1, bdf1_samples, 2, trial, design, 0.5_dp, out)
  call out % get_real(rv)
  call report(matches(rv, [7.5_dp, 11.5_dp, 16.5_dp], 1.0e-14_dp), &
       & "BDF1 residual(q_trial) = [7.5, 11.5, 16.5]", nfail)

  call dq % set_real([1.0_dp, 0.0_dp, 2.0_dp])
  call driver % jacobian_action(r_form, bdf1, bdf1_samples, 2, trial, dq, &
       & design, 0.5_dp, out)
  call out % get_real(rv)
  call report(matches(rv, [3.0_dp, 0.0_dp, 6.0_dp], 1.0e-14_dp), &
       & "BDF1 exact action: J dq = 3 dq = [3, 0, 6]", nfail)

  !-------------------------------------------------------------------!
  ! The one assembly verb: dense_jacobian probes the action whole,
  ! one basis direction per column - here exactly 3 I.
  !-------------------------------------------------------------------!

  call driver % dense_jacobian(r_form, bdf1, bdf1_samples, 2, trial, &
       & design, 0.5_dp, jac)
  call report(size(jac, 1) == 3 .and. size(jac, 2) == 3 .and. &
       & matches([jac(1,1), jac(2,2), jac(3,3)], &
       &         [3.0_dp, 3.0_dp, 3.0_dp], 1.0e-14_dp) .and. &
       & matches([jac(1,2), jac(1,3), jac(2,1), jac(2,3), jac(3,1), jac(3,2)], &
       &         [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], 1.0e-14_dp), &
       & "dense_jacobian assembles the whole J = 3 I, exactly", nfail)

  !-------------------------------------------------------------------!
  ! Newton on BDF1: one exact step lands on the closed-form root.
  !-------------------------------------------------------------------!

  call q_start % set_real([0.0_dp, 0.0_dp, 0.0_dp])
  call driver % solve(r_form, bdf1, bdf1_samples, 2, q_start, design, &
       & 0.5_dp, options, result)

  call result % q % get_real(rv)
  call report(matches(rv, root_bdf1, 1.0e-9_dp), &
       & "BDF1 Newton lands on q* = (2 q0 - xi)/3", nfail)

  call report(result % residual_norm <= options % residual_tolerance, &
       & "BDF1 residual at the root is below tolerance", nfail)

  call report(result % converged, &
       & "BDF1 reports converged", nfail)

  call report(result % iterations >= 1 .and. &
       & result % iterations <= options % max_iterations, &
       & "BDF1 spends a positive, bounded iteration count", nfail)

  !-------------------------------------------------------------------!
  ! DIRK and BDF2: the same three verbs, different rows.
  !-------------------------------------------------------------------!

  call driver % jacobian_action(r_form, dirk, dirk_samples, 2, trial, dq, &
       & design, 2.0_dp, out)
  call out % get_real(rv)
  call report(matches(rv, [2.0_dp, 0.0_dp, 4.0_dp], 1.0e-14_dp), &
       & "DIRK exact action: J dq = 2 dq = [2, 0, 4]", nfail)

  call driver % solve(r_form, dirk, dirk_samples, 2, q_start, design, &
       & 2.0_dp, options, result)
  call result % q % get_real(rv)
  call report(result % converged .and. matches(rv, root_dirk, 1.0e-9_dp), &
       & "DIRK Newton lands on q* = (q_base - xi)/2", nfail)

  call driver % jacobian_action(r_form, bdf2, bdf2_samples, 3, trial, dq, &
       & design, 1.0_dp, out)
  call out % get_real(rv)
  call report(matches(rv, [2.5_dp, 0.0_dp, 5.0_dp], 1.0e-14_dp), &
       & "BDF2 exact action: J dq = 2.5 dq = [2.5, 0, 5]", nfail)

  call driver % solve(r_form, bdf2, bdf2_samples, 3, q_start, design, &
       & 1.0_dp, options, result)
  call result % q % get_real(rv)
  call report(result % converged .and. matches(rv, root_bdf2, 1.0e-9_dp), &
       & "BDF2 Newton lands on its exact linear root", nfail)

  !-------------------------------------------------------------------!
  ! The solve never touches the caller's samples, and the result
  ! carries lawful shapes.
  !-------------------------------------------------------------------!

  call bdf1_samples(2) % state % component(1 + GTI_STATE_Q) % value % get_real_vector(rv)
  call report(matches(rv, [0.0_dp, 0.0_dp, 0.0_dp], 1.0e-14_dp), &
       & "the input samples are never mutated by a solve", nfail)

  call report(result % q % nentries == q_start % nentries .and. &
       & result % q % ncomp == q_start % ncomp, &
       & "the result q keeps the shape of q_initial", nfail)

  call report(result % residual % ncomp == 1 .and. &
       & result % residual % nentries == 3, &
       & "the result residual is a vector", nfail)

  !-------------------------------------------------------------------!
  ! The nonlinear scalar: genuine iteration, and lawful failure.
  !-------------------------------------------------------------------!

  call scalar_start % set_real([3.0_dp])

  starved % max_iterations = 1
  starved % residual_tolerance = 1.0e-30_dp
  call driver % solve(square, identity_motif, square_samples, 1, scalar_start, &
       & no_design, 0.0_dp, starved, result)
  call report((.not. result % converged) .and. result % iterations == 1, &
       & "a starved solve answers converged=false, never an error stop", nfail)

  call driver % solve(square, identity_motif, square_samples, 1, scalar_start, &
       & no_design, 0.0_dp, options, result)
  call result % q % get_real(rv)
  call report(result % converged .and. matches(rv, [2.0_dp], 1.0e-8_dp) .and. &
       & result % iterations > 1, &
       & "Newton walks q^2 = 4 from 3 down to 2 in several steps", nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all newton driver checks passed"
  else
     error stop
  end if

contains

  pure function matches(values, expected, tolerance) result(ok)

    real(dp), intent(in) :: values(:), expected(:)
    real(dp), intent(in) :: tolerance
    logical :: ok

    ok = size(values) == size(expected)
    if (ok) ok = all(abs(values - expected) < tolerance)

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

end program test_gti_time_local_newton_driver
