!=====================================================================!
! The refusals that must die at the reverse seat, one per
! invocation:
!
!      negsing        a zero singular tolerance
!      norelations    reversing a graph with no relations
!      idx0           adjoining relation 0
!      idxhigh        adjoining past the relations
!      seedsize       a seed array shorter than the vertices
!      noeta          a design direction with no values
!      unsolved       an unknown vertex nobody has solved
!      noseed         an unknown vertex with an empty seed
!      seedshape      an unknown seed of the wrong size
!      nolambda       propagating a valueless lambda
!      propshape      propagating into an occupied seed of the
!                     wrong shape
!      rdsize         a residual design action of two entries
!                     against three unknowns
!      singular       a design-only residual: R_xi lives, J_u = 0
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
  use gti_time_motif_builders, only : gti_time_motif_builder
  use gti_time_graphs      , only : gti_time_graph
  use gti_time_reverse_drivers, only : gti_time_reverse_driver, &
       & gti_time_reverse_options, gti_time_reverse_step_result, &
       & gti_time_reverse_result
  use gti_toy_forms        , only : toy_time_residual_form, &
       & toy_design_only_form, toy_short_form

  implicit none

  type(graph) :: states, designs
  type(field) :: q1_field, q2_field, xi_field

  type(gti_time_motif_builder) :: builder
  type(gti_time_graph)         :: time_graph
  type(gti_design_bundle)      :: design

  type(gti_time_reverse_driver)      :: reverse
  type(gti_time_reverse_options)     :: options
  type(gti_time_reverse_step_result) :: step
  type(gti_time_reverse_result)      :: sweep

  type(toy_time_residual_form) :: r_form
  type(toy_design_only_form)   :: design_only
  type(toy_short_form)         :: short

  type(gti_value_buffer) :: eta, lambda
  type(gti_value_buffer) :: seeds(2), short_seeds(1)
  character(len=32) :: which

  call get_command_argument(1, which)

  call states  % declare()
  call designs % declare()

  q1_field = field('q1', states, 3)
  call q1_field % set_real_vector([1.0_dp, 2.0_dp, 4.0_dp])
  q2_field = field('q2', states, 3)
  call q2_field % set_real_vector([0.5_dp, 1.0_dp, 2.5_dp])
  xi_field = field('xi', designs, 3)
  call xi_field % set_real_vector([0.5_dp, 0.5_dp, 0.5_dp])

  allocate(design % component(1))
  design % component(1) % value = xi_field

  !-------------------------------------------------------------------!
  ! The lawful solved two-vertex, one-relation graph every case
  ! corrupts. The unknown vertex's seed starts occupied.
  !-------------------------------------------------------------------!

  allocate(time_graph % vertex(2))
  allocate(time_graph % vertex(1) % sample % state % component(1))
  time_graph % vertex(1) % sample % state % component(1) % value = q1_field
  time_graph % vertex(1) % has_solution = .true.
  allocate(time_graph % vertex(2) % sample % state % component(1))
  time_graph % vertex(2) % sample % state % component(1) % value = q2_field
  time_graph % vertex(2) % has_solution = .true.

  allocate(time_graph % relation(1))
  call builder % bdf_uniform(1, 0.5_dp, time_graph % relation(1) % motif)
  time_graph % relation(1) % sample_vertex   = [1, 2]
  time_graph % relation(1) % unknown_sample  = 2
  time_graph % relation(1) % evaluation_time = 0.5_dp

  call eta    % set_real([1.0_dp, 1.0_dp, 1.0_dp])
  call lambda % set_real([1.0_dp, 1.0_dp, 1.0_dp])
  call seeds(2) % set_real([1.0_dp, 1.0_dp, 1.0_dp])

  r_form % nequations = 3

  select case (trim(which))

  case ('negsing')

     options % singular_tolerance = 0.0_dp
     call reverse % solve_relation_adjoint(r_form, time_graph, 1, design, eta, &
          & seeds, options, step)

  case ('norelations')

     deallocate(time_graph % relation)
     call reverse % reverse_all(r_form, time_graph, design, eta, seeds, &
          & options, sweep)

  case ('idx0')

     call reverse % solve_relation_adjoint(r_form, time_graph, 0, design, eta, &
          & seeds, options, step)

  case ('idxhigh')

     call reverse % solve_relation_adjoint(r_form, time_graph, 2, design, eta, &
          & seeds, options, step)

  case ('seedsize')

     call short_seeds(1) % set_real([1.0_dp, 1.0_dp, 1.0_dp])
     call reverse % solve_relation_adjoint(r_form, time_graph, 1, design, eta, &
          & short_seeds, options, step)

  case ('noeta')

     call eta % clear()
     call reverse % solve_relation_adjoint(r_form, time_graph, 1, design, eta, &
          & seeds, options, step)

  case ('unsolved')

     time_graph % vertex(2) % has_solution = .false.
     call reverse % solve_relation_adjoint(r_form, time_graph, 1, design, eta, &
          & seeds, options, step)

  case ('noseed')

     call seeds(2) % clear()
     call reverse % solve_relation_adjoint(r_form, time_graph, 1, design, eta, &
          & seeds, options, step)

  case ('seedshape')

     call seeds(2) % set_real([1.0_dp, 1.0_dp])
     call reverse % solve_relation_adjoint(r_form, time_graph, 1, design, eta, &
          & seeds, options, step)

  case ('nolambda')

     call lambda % clear()
     call reverse % propagate_relation(r_form, time_graph, 1, design, lambda, seeds)

  case ('propshape')

     call seeds(1) % set_real([1.0_dp, 2.0_dp])
     call reverse % propagate_relation(r_form, time_graph, 1, design, lambda, seeds)

  case ('rdsize')

     call reverse % solve_relation_adjoint(short, time_graph, 1, design, eta, &
          & seeds, options, step)

  case ('singular')

     call reverse % solve_relation_adjoint(design_only, time_graph, 1, design, &
          & eta, seeds, options, step)

  case default

     error stop 'refusal: unknown case'

  end select

  write(*,*) 'refusal case survived: ', trim(which)

end program refusal
