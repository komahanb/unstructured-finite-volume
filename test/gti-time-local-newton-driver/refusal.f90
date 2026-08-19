!=====================================================================!
! The refusals that must die at the Newton seat, one per
! invocation:
!
!      negiter    zero iterations allowed
!      negrtol    a negative residual tolerance
!      negstep    a negative step tolerance
!      negsing    a zero singular tolerance
!      novalues   an initial q with no values
!      dqshape    a direction whose shape contradicts the trial
!      wide       a residual of two components - not a vector
!      short      two equations against three unknowns
!      shortjac   the same short form probed through the one
!                 assembly verb, dense_jacobian
!      singular   a constant residual, whose Jacobian is zero
!
! Every case must error stop; a case that returns is a failure of
! the suite. Convergence failure is NOT among these - that is an
! answer, and the law suite proves it.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program refusal

  use iso_fortran_env      , only : dp => REAL64
  use fractal_graph        , only : graph
  use class_graph_field    , only : field
  use gti_value_buffers    , only : gti_value_buffer
  use gti_design_bundles   , only : gti_design_bundle
  use gti_time_local_schemes, only : gti_time_sample, gti_time_motif
  use gti_time_motif_builders, only : gti_time_motif_builder
  use gti_time_local_newton_drivers, only : gti_time_local_newton_driver, &
       & gti_time_local_newton_options, gti_time_local_newton_result
  use gti_toy_forms        , only : toy_time_residual_form, toy_constant_form, &
       & toy_wide_form, toy_short_form

  implicit none

  type(graph) :: states, designs
  type(field) :: q0_field, zero_field, xi_field

  type(gti_time_sample)   :: samples(2)
  type(gti_design_bundle) :: design

  type(gti_time_motif_builder)        :: builder
  type(gti_time_motif)                :: bdf1
  type(gti_time_local_newton_driver)  :: driver
  type(gti_time_local_newton_options) :: options
  type(gti_time_local_newton_result)  :: result

  type(toy_time_residual_form) :: r_form
  type(toy_constant_form)      :: constant
  type(toy_wide_form)          :: wide
  type(toy_short_form)         :: short

  type(gti_value_buffer) :: trial, dq, out
  real(dp), allocatable :: jac(:,:)
  character(len=32) :: which

  call get_command_argument(1, which)

  call states  % declare()
  call designs % declare()

  q0_field = field('q0', states, 3)
  call q0_field % set_real_vector([1.0_dp, 2.0_dp, 4.0_dp])
  zero_field = field('placeholder', states, 3)
  call zero_field % set_real_vector([0.0_dp, 0.0_dp, 0.0_dp])
  xi_field = field('xi', designs, 1)
  call xi_field % set_real_vector([0.5_dp])

  allocate(samples(1) % state % component(1))
  samples(1) % state % component(1) % value = q0_field
  allocate(samples(2) % state % component(1))
  samples(2) % state % component(1) % value = zero_field

  allocate(design % component(1))
  design % component(1) % value = xi_field

  call builder % bdf_uniform(1, 0.5_dp, bdf1)

  call trial % set_real([3.0_dp, 5.0_dp, 8.0_dp])

  r_form % nequations = 3

  select case (trim(which))

  case ('negiter')

     options % max_iterations = 0
     call options % validate()

  case ('negrtol')

     options % residual_tolerance = -1.0_dp
     call options % validate()

  case ('negstep')

     options % step_tolerance = -1.0_dp
     call options % validate()

  case ('negsing')

     options % singular_tolerance = 0.0_dp
     call options % validate()

  case ('novalues')

     call trial % clear()
     call driver % solve(r_form, bdf1, samples, 2, trial, design, &
          & 0.5_dp, options, result)

  case ('dqshape')

     call dq % set_real([1.0_dp, 0.0_dp])
     call driver % jacobian_action(r_form, bdf1, samples, 2, trial, dq, &
          & design, 0.5_dp, out)

  case ('wide')

     call driver % residual(wide, bdf1, samples, 2, trial, design, 0.5_dp, out)

  case ('short')

     call driver % solve(short, bdf1, samples, 2, trial, design, &
          & 0.5_dp, options, result)

  case ('shortjac')

     ! the one assembly verb holds each probed column to the
     ! unknown's size - the same law the residual answers to
     call driver % dense_jacobian(short, bdf1, samples, 2, trial, &
          & design, 0.5_dp, jac)

  case ('singular')

     call driver % solve(constant, bdf1, samples, 2, trial, design, &
          & 0.5_dp, options, result)

  case default

     error stop 'refusal: unknown case'

  end select

  write(*,*) 'refusal case survived: ', trim(which)

end program refusal
