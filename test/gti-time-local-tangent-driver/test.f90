!=====================================================================!
! The time-local tangent driver suite: at a converged local state
! R(q*, xi, t) = 0, the tangent along a design direction solves
!
!      J_q q_eta = - R_xi [eta],
!
! with J_q reused verbatim from the Newton driver's exact
! motif-row action:
!
!      BDF1  J = 3 I     eta = 1 -> q_eta = -1/3
!      DIRK  J = 2 I     eta = 2 -> q_eta = -1
!      BDF2  J = 2.5 I   eta = [1,2,3] -> q_eta = [-0.4,-0.8,-1.2]
!
! plus the nonlinear scalar S = q^2 + xi - 4 at q* = 2, xi = 0:
! J = 4, eta = 1, q_eta = -1/4. Every linear residual is held
! below 1e-12, and no input sample is ever touched.
!
! (The task sheet's BDF2 root reads [2.2, 3.4, 5.4]; the exact
! root at xi = 0.5 is [2, 3.4, 5.4], used here - for this linear
! toy J_q and R_xi are constant, so the tangent laws are
! unchanged either way.)
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_gti_time_local_tangent_driver

  use iso_fortran_env      , only : dp => REAL64
  use fractal_graph        , only : graph
  use class_graph_field    , only : field
  use gti_value_buffers    , only : gti_value_buffer
  use gti_state_bundles    , only : GTI_STATE_Q
  use gti_design_bundles   , only : gti_design_bundle
  use gti_time_local_schemes, only : gti_time_sample, gti_time_motif
  use gti_time_motif_builders, only : gti_time_motif_builder
  use gti_time_local_tangent_drivers, only : gti_time_local_tangent_driver, &
       & gti_time_local_tangent_result
  use gti_toy_forms        , only : toy_time_residual_form, toy_square_design_form

  implicit none

  type(graph) :: states, scalars, designs
  type(field) :: q0_field, q1_field, base_field, zero_field, zero1_field
  type(field) :: xi_field, xi0_field

  type(gti_time_sample)   :: bdf1_samples(2), dirk_samples(2), bdf2_samples(3)
  type(gti_time_sample)   :: square_samples(1)
  type(gti_design_bundle) :: design, design0

  type(gti_time_motif_builder) :: builder
  type(gti_time_motif)         :: bdf1, dirk, bdf2, identity_motif

  type(gti_time_local_tangent_driver) :: tangent
  type(gti_time_local_tangent_result) :: result

  type(toy_time_residual_form)  :: r_form
  type(toy_square_design_form)  :: square

  type(gti_value_buffer) :: qstar_bdf1, qstar_dirk, qstar_bdf2, qstar_square
  type(gti_value_buffer) :: eta1, eta2, eta3, eta_scalar, rhs

  real(dp), allocatable :: rv(:)
  real(dp), parameter :: tolpiv = 1.0e-14_dp
  integer :: nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "gti time-local tangent driver suite"
  write(*,'(1x,a)') "============================================="

  !-------------------------------------------------------------------!
  ! Fields, samples, motifs, converged states, and directions.
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
  xi0_field = field('xi0', designs, 1)
  call xi0_field % set_real_vector([0.0_dp])

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
  allocate(design0 % component(1))
  design0 % component(1) % value = xi0_field

  call builder % bdf_uniform(1, 0.5_dp, bdf1)
  call builder % dirk_stage(0.5_dp, 2.0_dp, dirk)
  call builder % bdf_uniform(2, 1.0_dp, bdf2)

  allocate(identity_motif % rule(1))
  identity_motif % rule(1) % state_component = GTI_STATE_Q
  identity_motif % rule(1) % weight = [1.0_dp]

  call qstar_bdf1 % set_real((2.0_dp * [1.0_dp, 2.0_dp, 4.0_dp] - 0.5_dp) / 3.0_dp)
  call qstar_dirk % set_real(([2.0_dp, 4.0_dp, 6.0_dp] - 0.5_dp) / 2.0_dp)
  call qstar_bdf2 % set_real((2.0_dp * [3.0_dp, 5.0_dp, 8.0_dp] &
       & - 0.5_dp * [1.0_dp, 2.0_dp, 4.0_dp] - 0.5_dp) / 2.5_dp)
  call qstar_square % set_real([2.0_dp])

  call eta1 % set_real([1.0_dp, 1.0_dp, 1.0_dp])
  call eta2 % set_real([2.0_dp, 2.0_dp, 2.0_dp])
  call eta3 % set_real([1.0_dp, 2.0_dp, 3.0_dp])
  call eta_scalar % set_real([1.0_dp])

  r_form % nequations = 3

  !-------------------------------------------------------------------!
  ! BDF1: the explicit design RHS and the tangent.
  !-------------------------------------------------------------------!

  call tangent % design_rhs(r_form, bdf1, bdf1_samples, 2, qstar_bdf1, &
       & design, eta1, 0.5_dp, rhs)
  call rhs % get_real(rv)
  call report(matches(rv, [-1.0_dp, -1.0_dp, -1.0_dp], 1.0e-14_dp), &
       & "BDF1 design rhs = -R_xi eta = [-1, -1, -1]", nfail)

  call bdf1_samples(2) % state % component(1 + GTI_STATE_Q) % value % get_real_vector(rv)
  call report(matches(rv, [0.0_dp, 0.0_dp, 0.0_dp], 1.0e-14_dp), &
       & "design_rhs never touches the caller's samples", nfail)

  call tangent % solve_design_tangent(r_form, bdf1, bdf1_samples, 2, qstar_bdf1, &
       & design, eta1, 0.5_dp, tolpiv, result)
  call result % q_tangent % get_real(rv)
  call report(matches(rv, [-1.0_dp, -1.0_dp, -1.0_dp] / 3.0_dp, 1.0e-12_dp), &
       & "BDF1 tangent: J = 3 I gives q_eta = -eta/3", nfail)

  call report(result % linear_residual_norm <= 1.0e-12_dp, &
       & "BDF1 linear residual norm is below 1e-12", nfail)

  call report(result % solved, &
       & "the result reports solved after a successful solve", nfail)

  call report(result % q_tangent % nentries == qstar_bdf1 % nentries .and. &
       & result % q_tangent % ncomp == qstar_bdf1 % ncomp, &
       & "the tangent keeps the shape of q star", nfail)

  call result % rhs % get_real(rv)
  call report(matches(rv, [-1.0_dp, -1.0_dp, -1.0_dp], 1.0e-14_dp), &
       & "the result carries the rhs it answered", nfail)

  !-------------------------------------------------------------------!
  ! DIRK and BDF2: same verbs, different rows.
  !-------------------------------------------------------------------!

  call tangent % design_rhs(r_form, dirk, dirk_samples, 2, qstar_dirk, &
       & design, eta2, 2.0_dp, rhs)
  call rhs % get_real(rv)
  call report(matches(rv, [-2.0_dp, -2.0_dp, -2.0_dp], 1.0e-14_dp), &
       & "DIRK design rhs = [-2, -2, -2]", nfail)

  call tangent % solve_design_tangent(r_form, dirk, dirk_samples, 2, qstar_dirk, &
       & design, eta2, 2.0_dp, tolpiv, result)
  call result % q_tangent % get_real(rv)
  call report(result % solved .and. &
       & matches(rv, [-1.0_dp, -1.0_dp, -1.0_dp], 1.0e-12_dp) .and. &
       & result % linear_residual_norm <= 1.0e-12_dp, &
       & "DIRK tangent: J = 2 I gives q_eta = [-1, -1, -1], residual tiny", nfail)

  call tangent % design_rhs(r_form, bdf2, bdf2_samples, 3, qstar_bdf2, &
       & design, eta3, 1.0_dp, rhs)
  call rhs % get_real(rv)
  call report(matches(rv, [-1.0_dp, -2.0_dp, -3.0_dp], 1.0e-14_dp), &
       & "BDF2 design rhs = [-1, -2, -3]", nfail)

  call tangent % solve_design_tangent(r_form, bdf2, bdf2_samples, 3, qstar_bdf2, &
       & design, eta3, 1.0_dp, tolpiv, result)
  call result % q_tangent % get_real(rv)
  call report(result % solved .and. &
       & matches(rv, [-0.4_dp, -0.8_dp, -1.2_dp], 1.0e-12_dp) .and. &
       & result % linear_residual_norm <= 1.0e-12_dp, &
       & "BDF2 tangent: J = 2.5 I gives q_eta = [-0.4, -0.8, -1.2]", nfail)

  call bdf2_samples(3) % state % component(1 + GTI_STATE_Q) % value % get_real_vector(rv)
  call report(matches(rv, [0.0_dp, 0.0_dp, 0.0_dp], 1.0e-14_dp), &
       & "solve_design_tangent never touches the caller's samples", nfail)

  !-------------------------------------------------------------------!
  ! The nonlinear scalar: S = q^2 + xi - 4 at q* = 2, xi = 0.
  !-------------------------------------------------------------------!

  call tangent % solve_design_tangent(square, identity_motif, square_samples, 1, &
       & qstar_square, design0, eta_scalar, 0.0_dp, tolpiv, result)
  call result % q_tangent % get_real(rv)
  call report(result % solved .and. matches(rv, [-0.25_dp], 1.0e-12_dp), &
       & "the nonlinear scalar tangent: J = 4 gives q_eta = -1/4", nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all tangent driver checks passed"
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

end program test_gti_time_local_tangent_driver
