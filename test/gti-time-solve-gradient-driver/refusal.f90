!=====================================================================!
! The refusals that must die through the composition, one per
! invocation - each PROPAGATED from the seat that owns the law:
!
!      fwdopt    a zero Newton iteration budget, refused by the
!                local Newton options
!      revopt    a zero singular tolerance, refused by the
!                reverse options
!      badgraph  a relation without an unknown, refused by the
!                graph representation
!      zeroterm  a functional without terms, refused by the
!                functional representation - AFTER a converged
!                march
!      noeta     a valueless design direction, refused by the
!                seed driver - after a converged march
!
! The composition duplicates no lower-layer law: what dies here
! dies with the owner's message. Its own two scalar guards on the
! design actions protect against future lower-driver changes; no
! public path can trip them today, and run.sh proves their
! presence in source instead.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program refusal

  use iso_fortran_env      , only : dp => REAL64
  use fractal_graph        , only : graph
  use class_graph_field    , only : field
  use gti_value_buffers    , only : gti_value_buffer
  use gti_design_bundles   , only : gti_design_bundle
  use gti_time_motif_builders, only : gti_time_motif_builder
  use gti_time_graphs      , only : gti_time_graph
  use gti_time_functional_seed_drivers, only : gti_time_functional
  use gti_time_solve_gradient_drivers, only : gti_time_solve_gradient_driver, &
       & gti_time_solve_gradient_options, gti_time_solve_gradient_result
  use gti_toy_forms        , only : toy_time_residual_form, &
       & toy_sum_time_functional

  implicit none

  type(graph) :: states, designs
  type(field) :: q1_field, zero_field, xi_field

  type(gti_time_motif_builder) :: builder
  type(gti_time_graph)         :: time_graph
  type(gti_design_bundle)      :: design

  type(gti_time_functional)             :: functional
  type(gti_time_solve_gradient_driver)  :: driver
  type(gti_time_solve_gradient_options) :: options
  type(gti_time_solve_gradient_result)  :: result

  type(toy_time_residual_form)  :: r_form
  type(toy_sum_time_functional) :: f_form

  type(gti_value_buffer) :: eta
  character(len=32) :: which

  call get_command_argument(1, which)

  call states  % declare()
  call designs % declare()

  q1_field = field('q1', states, 3)
  call q1_field % set_real_vector([1.0_dp, 2.0_dp, 4.0_dp])
  zero_field = field('placeholder', states, 3)
  call zero_field % set_real_vector([0.0_dp, 0.0_dp, 0.0_dp])
  xi_field = field('xi', designs, 3)
  call xi_field % set_real_vector([0.5_dp, 0.5_dp, 0.5_dp])

  allocate(design % component(1))
  design % component(1) % value = xi_field

  !-------------------------------------------------------------------!
  ! The lawful two-vertex, one-relation graph and a one-term
  ! functional, corrupted per case.
  !-------------------------------------------------------------------!

  allocate(time_graph % vertex(2))
  allocate(time_graph % vertex(1) % sample % state % component(1))
  time_graph % vertex(1) % sample % state % component(1) % value = q1_field
  time_graph % vertex(1) % has_solution = .true.
  allocate(time_graph % vertex(2) % sample % state % component(1))
  time_graph % vertex(2) % sample % state % component(1) % value = zero_field
  time_graph % vertex(2) % sample % time = 0.5_dp

  allocate(time_graph % relation(1))
  call builder % bdf_uniform(1, 0.5_dp, time_graph % relation(1) % motif)
  time_graph % relation(1) % sample_vertex   = [1, 2]
  time_graph % relation(1) % unknown_sample  = 2
  time_graph % relation(1) % evaluation_time = 0.5_dp

  allocate(functional % term(1))
  functional % term(1) % vertex_index    = 2
  functional % term(1) % evaluation_time = 0.5_dp

  call eta % set_real([1.0_dp, 1.0_dp, 1.0_dp])

  r_form % nequations = 3

  select case (trim(which))

  case ('fwdopt')

     options % forward % newton % max_iterations = 0
     call driver % solve(r_form, f_form, time_graph, functional, design, eta, &
          & options, result)

  case ('revopt')

     options % reverse % singular_tolerance = 0.0_dp
     call driver % solve(r_form, f_form, time_graph, functional, design, eta, &
          & options, result)

  case ('badgraph')

     time_graph % relation(1) % unknown_sample = 0
     call driver % solve(r_form, f_form, time_graph, functional, design, eta, &
          & options, result)

  case ('zeroterm')

     deallocate(functional % term)
     call driver % solve(r_form, f_form, time_graph, functional, design, eta, &
          & options, result)

  case ('noeta')

     call eta % clear()
     call driver % solve(r_form, f_form, time_graph, functional, design, eta, &
          & options, result)

  case default

     error stop 'refusal: unknown case'

  end select

  write(*,*) 'refusal case survived: ', trim(which)

end program refusal
