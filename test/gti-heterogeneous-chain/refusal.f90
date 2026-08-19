!=====================================================================!
! The propagation refusals of the heterogeneous chain, one per
! invocation - no new law, no new message: what dies here dies
! with the message of the seat that has always owned it:
!
!      reorder       DIRK stored before BDF1 reads unsolved v2 as
!                    history - the forward driver's history law
!      badweight     an ABM2 row weighted thrice over two
!                    vertices - the graph representation's law
!      unsolvedterm  a functional term on a vertex no relation
!                    solves - the functional representation's law,
!                    firing AFTER a fully converged march
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
  integer :: v
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

  call eta % set_real([1.0_dp, 1.0_dp, 1.0_dp])

  r_form % nequations = 3

  !-------------------------------------------------------------------!
  ! The lawful heterogeneous graph plus one unreferenced vertex 5,
  ! corrupted per case.
  !-------------------------------------------------------------------!

  allocate(time_graph % vertex(5))
  allocate(time_graph % vertex(1) % sample % state % component(1))
  time_graph % vertex(1) % sample % state % component(1) % value = q1_field
  time_graph % vertex(1) % has_solution = .true.
  do v = 2, 5
     allocate(time_graph % vertex(v) % sample % state % component(1))
     time_graph % vertex(v) % sample % state % component(1) % value = zero_field
  end do

  allocate(time_graph % relation(3))
  call builder % bdf_uniform(1, 0.5_dp, time_graph % relation(1) % motif)
  time_graph % relation(1) % sample_vertex   = [1, 2]
  time_graph % relation(1) % unknown_sample  = 2
  time_graph % relation(1) % evaluation_time = 0.5_dp
  call builder % dirk_stage(0.5_dp, 2.0_dp, time_graph % relation(2) % motif)
  time_graph % relation(2) % sample_vertex   = [2, 3]
  time_graph % relation(2) % unknown_sample  = 2
  time_graph % relation(2) % evaluation_time = 1.0_dp
  call builder % abm_corrector(2, 1.0_dp, time_graph % relation(3) % motif)
  time_graph % relation(3) % sample_vertex   = [3, 4]
  time_graph % relation(3) % unknown_sample  = 2
  time_graph % relation(3) % evaluation_time = 2.0_dp

  allocate(functional % term(1))
  functional % term(1) % vertex_index    = 4
  functional % term(1) % evaluation_time = 2.0_dp

  select case (trim(which))

  case ('reorder')

     ! DIRK first: it reads v2 as history before BDF1 solves it
     time_graph % relation(1) % sample_vertex = [2, 3]
     call builder % dirk_stage(0.5_dp, 2.0_dp, time_graph % relation(1) % motif)
     time_graph % relation(2) % sample_vertex = [1, 2]
     call builder % bdf_uniform(1, 0.5_dp, time_graph % relation(2) % motif)
     call driver % solve(r_form, f_form, time_graph, functional, design, eta, &
          & options, result)

  case ('badweight')

     time_graph % relation(3) % motif % rule(1) % weight = &
          & [1.0_dp, 0.0_dp, 0.0_dp]
     call driver % solve(r_form, f_form, time_graph, functional, design, eta, &
          & options, result)

  case ('unsolvedterm')

     ! a term on the vertex no relation solves: the march converges
     ! everything it owns, and the functional law still refuses
     functional % term(1) % vertex_index = 5
     call driver % solve(r_form, f_form, time_graph, functional, design, eta, &
          & options, result)

  case default

     error stop 'refusal: unknown case'

  end select

  write(*,*) 'refusal case survived: ', trim(which)

end program refusal
