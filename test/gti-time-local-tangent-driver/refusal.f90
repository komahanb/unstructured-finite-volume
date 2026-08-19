!=====================================================================!
! The refusals that must die at the tangent seat, one per
! invocation:
!
!      negsing    a zero singular tolerance
!      noqstar    a q star with no values
!      noeta      a design direction with no values
!      wide       a two-component design partial - not a vector
!      rhssize    two equations against three unknowns
!      singular   a design-only residual: R_xi lives, J_q = 0
!
! Every case must error stop; a case that returns is a failure of
! the suite.
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
  use gti_time_local_tangent_drivers, only : gti_time_local_tangent_driver, &
       & gti_time_local_tangent_result
  use gti_toy_forms        , only : toy_time_residual_form, &
       & toy_design_only_form, toy_wide_form, toy_short_form

  implicit none

  type(graph) :: states, designs
  type(field) :: q0_field, zero_field, xi_field

  type(gti_time_sample)   :: samples(2)
  type(gti_design_bundle) :: design

  type(gti_time_motif_builder)        :: builder
  type(gti_time_motif)                :: bdf1
  type(gti_time_local_tangent_driver) :: tangent
  type(gti_time_local_tangent_result) :: result

  type(toy_time_residual_form) :: r_form
  type(toy_design_only_form)   :: design_only
  type(toy_wide_form)          :: wide
  type(toy_short_form)         :: short

  type(gti_value_buffer) :: qstar, eta, rhs
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

  call qstar % set_real([0.5_dp, 1.0_dp, 2.5_dp])
  call eta   % set_real([1.0_dp, 1.0_dp, 1.0_dp])

  r_form % nequations = 3

  select case (trim(which))

  case ('negsing')

     call tangent % solve_design_tangent(r_form, bdf1, samples, 2, qstar, &
          & design, eta, 0.5_dp, 0.0_dp, result)

  case ('noqstar')

     call qstar % clear()
     call tangent % solve_design_tangent(r_form, bdf1, samples, 2, qstar, &
          & design, eta, 0.5_dp, 1.0e-14_dp, result)

  case ('noeta')

     call eta % clear()
     call tangent % design_rhs(r_form, bdf1, samples, 2, qstar, &
          & design, eta, 0.5_dp, rhs)

  case ('wide')

     call tangent % design_rhs(wide, bdf1, samples, 2, qstar, &
          & design, eta, 0.5_dp, rhs)

  case ('rhssize')

     call tangent % solve_design_tangent(short, bdf1, samples, 2, qstar, &
          & design, eta, 0.5_dp, 1.0e-14_dp, result)

  case ('singular')

     call tangent % solve_design_tangent(design_only, bdf1, samples, 2, qstar, &
          & design, eta, 0.5_dp, 1.0e-14_dp, result)

  case default

     error stop 'refusal: unknown case'

  end select

  write(*,*) 'refusal case survived: ', trim(which)

end program refusal
