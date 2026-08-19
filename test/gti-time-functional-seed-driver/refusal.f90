!=====================================================================!
! The refusals that must die at the functional seed seat, one per
! invocation:
!
!      zeroterm    a functional with no terms at all
!      vertex0     a term addressing vertex 0
!      vertexhigh  a term addressing past the graph
!      unsolved    a term on a vertex nobody has solved
!      novertexq   a term on a vertex whose sample holds no q
!      seedsize    a seed array shorter than the vertices
!      noeta       a design direction with no values
!      qnovalues   a vertex q holding no real values
!      valshort    a "functional" of two entries at value
!      valwide     a "functional" of two components at value
!      gradshort   the same two entries at the q-gradient
!      gradwide    the same two components at the q-gradient
!      dashort     the same two entries at the design action
!      dawide      the same two components at the design action
!      propshape   seeding into an occupied seed of the wrong
!                  shape
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
  use gti_time_functional_seed_drivers, only : gti_time_functional_term, &
       & gti_time_functional, gti_time_functional_seed_result, &
       & gti_time_functional_seed_driver
  use gti_toy_forms        , only : toy_sum_time_functional, toy_short_form, &
       & toy_scalar_wide_form

  implicit none

  type(graph) :: states, designs
  type(field) :: q1_field, q2_field, int_field, xi_field

  type(gti_time_motif_builder) :: builder
  type(gti_time_graph)         :: time_graph
  type(gti_design_bundle)      :: design

  type(gti_time_functional)             :: functional
  type(gti_time_functional_seed_driver) :: seeder
  type(gti_time_functional_seed_result) :: seeding

  type(toy_sum_time_functional) :: f_form
  type(toy_short_form)          :: short
  type(toy_scalar_wide_form)    :: scalar_wide

  type(gti_value_buffer) :: eta, out
  type(gti_value_buffer) :: seeds(2), short_seeds(1)
  character(len=32) :: which

  call get_command_argument(1, which)

  call states  % declare()
  call designs % declare()

  q1_field = field('q1', states, 3)
  call q1_field % set_real_vector([1.0_dp, 2.0_dp, 4.0_dp])
  q2_field = field('q2', states, 3)
  call q2_field % set_real_vector([0.5_dp, 1.0_dp, 2.5_dp])
  int_field = field('integers', states, 3)
  call int_field % set_integer_vector([0, 0, 0])
  xi_field = field('xi', designs, 3)
  call xi_field % set_real_vector([0.5_dp, 0.5_dp, 0.5_dp])

  allocate(design % component(1))
  design % component(1) % value = xi_field

  !-------------------------------------------------------------------!
  ! The lawful solved two-vertex, one-relation graph, and a
  ! one-term functional on vertex 2, corrupted per case.
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

  allocate(functional % term(1))
  functional % term(1) % vertex_index    = 2
  functional % term(1) % evaluation_time = 0.5_dp

  call eta % set_real([1.0_dp, 1.0_dp, 1.0_dp])

  select case (trim(which))

  case ('zeroterm')

     deallocate(functional % term)
     call seeder % seed_all(f_form, time_graph, functional, design, eta, &
          & seeds, seeding)

  case ('vertex0')

     functional % term(1) % vertex_index = 0
     call functional % validate(time_graph)

  case ('vertexhigh')

     functional % term(1) % vertex_index = 9
     call functional % validate(time_graph)

  case ('unsolved')

     time_graph % vertex(2) % has_solution = .false.
     call functional % validate(time_graph)

  case ('novertexq')

     deallocate(time_graph % vertex(2) % sample % state % component)
     call functional % validate(time_graph)

  case ('seedsize')

     call seeder % seed_all(f_form, time_graph, functional, design, eta, &
          & short_seeds, seeding)

  case ('noeta')

     call eta % clear()
     call seeder % design_action(f_form, time_graph, functional % term(1), &
          & design, eta, out)

  case ('qnovalues')

     time_graph % vertex(2) % sample % state % component(1) % value = int_field
     call seeder % q_gradient(f_form, time_graph, functional % term(1), &
          & design, out)

  case ('valshort')

     call seeder % value(short, time_graph, functional % term(1), design, out)

  case ('valwide')

     call seeder % value(scalar_wide, time_graph, functional % term(1), &
          & design, out)

  case ('gradshort')

     call seeder % q_gradient(short, time_graph, functional % term(1), &
          & design, out)

  case ('gradwide')

     call seeder % q_gradient(scalar_wide, time_graph, functional % term(1), &
          & design, out)

  case ('dashort')

     call seeder % design_action(short, time_graph, functional % term(1), &
          & design, eta, out)

  case ('dawide')

     call seeder % design_action(scalar_wide, time_graph, functional % term(1), &
          & design, eta, out)

  case ('propshape')

     call seeds(2) % set_real([1.0_dp, 2.0_dp])
     call seeder % seed_all(f_form, time_graph, functional, design, eta, &
          & seeds, seeding)

  case default

     error stop 'refusal: unknown case'

  end select

  write(*,*) 'refusal case survived: ', trim(which)

end program refusal
