!=====================================================================!
! The refusals that must die at the march seat, one per
! invocation:
!
!      norelations      a graph with vertices and no relations
!      idx0             solving relation 0
!      idxhigh          solving past the relations
!      unsolvedhistory  a history vertex nobody has solved
!      noq              an unknown whose q field holds no reals
!      ncomp2           an unknown whose q field is not a flat
!                       vector
!
! Every case must error stop; a case that returns is a failure of
! the suite. Non-convergence is NOT among these - that is an
! answer, and the law suite proves it.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program refusal

  use iso_fortran_env      , only : dp => REAL64
  use fractal_graph        , only : graph
  use class_graph_field    , only : field
  use gti_design_bundles   , only : gti_design_bundle
  use gti_time_motif_builders, only : gti_time_motif_builder
  use gti_time_graphs      , only : gti_time_graph
  use gti_time_forward_drivers, only : gti_time_forward_driver, &
       & gti_time_forward_options, gti_time_forward_step_result, &
       & gti_time_forward_result
  use gti_toy_forms        , only : toy_time_residual_form

  implicit none

  type(graph) :: states, designs
  type(field) :: q1_field, zero_field, int_field, block_field, xi_field

  type(gti_time_motif_builder) :: builder
  type(gti_time_graph)         :: time_graph
  type(gti_design_bundle)      :: design

  type(gti_time_forward_driver)      :: driver
  type(gti_time_forward_options)     :: options
  type(gti_time_forward_step_result) :: step
  type(gti_time_forward_result)      :: march

  type(toy_time_residual_form) :: r_form
  character(len=32) :: which

  call get_command_argument(1, which)

  call states  % declare()
  call designs % declare()

  q1_field = field('q1', states, 3)
  call q1_field % set_real_vector([1.0_dp, 2.0_dp, 4.0_dp])
  zero_field = field('placeholder', states, 3)
  call zero_field % set_real_vector([0.0_dp, 0.0_dp, 0.0_dp])
  int_field = field('integers', states, 3)
  call int_field % set_integer_vector([0, 0, 0])
  block_field = field('block', states, 1, ncomp=3)
  call block_field % set_real_vector([0.0_dp, 0.0_dp, 0.0_dp])
  xi_field = field('xi', designs, 3)
  call xi_field % set_real_vector([0.5_dp, 0.5_dp, 0.5_dp])

  allocate(design % component(1))
  design % component(1) % value = xi_field

  !-------------------------------------------------------------------!
  ! The lawful two-vertex, one-relation graph every case corrupts.
  !-------------------------------------------------------------------!

  allocate(time_graph % vertex(2))
  allocate(time_graph % vertex(1) % sample % state % component(1))
  time_graph % vertex(1) % sample % state % component(1) % value = q1_field
  time_graph % vertex(1) % has_solution = .true.
  allocate(time_graph % vertex(2) % sample % state % component(1))
  time_graph % vertex(2) % sample % state % component(1) % value = zero_field

  allocate(time_graph % relation(1))
  call builder % bdf_uniform(1, 0.5_dp, time_graph % relation(1) % motif)
  time_graph % relation(1) % sample_vertex   = [1, 2]
  time_graph % relation(1) % unknown_sample  = 2
  time_graph % relation(1) % evaluation_time = 0.5_dp

  r_form % nequations = 3

  select case (trim(which))

  case ('norelations')

     deallocate(time_graph % relation)
     call driver % solve_all(r_form, time_graph, design, options, march)

  case ('idx0')

     call driver % solve_relation(r_form, time_graph, 0, design, options, step)

  case ('idxhigh')

     call driver % solve_relation(r_form, time_graph, 2, design, options, step)

  case ('unsolvedhistory')

     time_graph % vertex(1) % has_solution = .false.
     call driver % solve_relation(r_form, time_graph, 1, design, options, step)

  case ('noq')

     time_graph % vertex(2) % sample % state % component(1) % value = int_field
     call driver % solve_relation(r_form, time_graph, 1, design, options, step)

  case ('ncomp2')

     time_graph % vertex(2) % sample % state % component(1) % value = block_field
     call driver % solve_relation(r_form, time_graph, 1, design, options, step)

  case default

     error stop 'refusal: unknown case'

  end select

  write(*,*) 'refusal case survived: ', trim(which)

end program refusal
