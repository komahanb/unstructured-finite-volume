!=====================================================================!
! The refusals that must die at the degree-r seat, one per
! invocation:
!
!      deg0             a degree below one - degree has no upper
!                       cap here; the form contract bounds it
!      negsing          a zero singular tolerance
!      norelations      traversing a graph with no relations
!      idx0 / idxhigh   a relation index outside the graph
!      degdim           a derivative array whose degree dimension
!                       disagrees with the options
!      vertdim          a derivative array whose vertex dimension
!                       disagrees with the graph
!      noeta            a design direction with no values
!      unsolved         an unknown vertex nobody has solved
!      unsolvedhistory  a history vertex nobody has solved
!      noq              an unknown q holding no real values
!      badseed          a preoccupied derivative seat of the
!                       wrong size
!      shortres         two assembled equations against one
!                       unknown
!      singular         a design-only residual: J_u = 0
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
  use gti_time_degree_r_directional_drivers, only : &
       & gti_time_degree_r_directional_driver, &
       & gti_time_degree_r_directional_options, &
       & gti_time_degree_r_relation_result, &
       & gti_time_degree_r_directional_result
  use gti_toy_forms        , only : toy_qdot_square_form, &
       & toy_design_only_scalar_form, toy_short_form

  implicit none

  type(graph) :: states, designs
  type(field) :: one_field, int_field, xi_field

  type(gti_time_motif_builder) :: builder
  type(gti_time_graph)         :: time_graph
  type(gti_design_bundle)      :: design

  type(gti_time_degree_r_directional_driver)  :: driver
  type(gti_time_degree_r_directional_options) :: options
  type(gti_time_degree_r_relation_result)     :: step
  type(gti_time_degree_r_directional_result)  :: result

  type(toy_qdot_square_form)        :: r_form
  type(toy_design_only_scalar_form) :: design_only
  type(toy_short_form)              :: short

  type(gti_value_buffer) :: eta
  type(gti_value_buffer) :: vd(1, 2), vd_deep(2, 2), vd_narrow(1, 1)
  integer :: v
  character(len=32) :: which

  call get_command_argument(1, which)

  call states  % declare()
  call designs % declare()

  one_field = field('one', states, 1)
  call one_field % set_real_vector([1.0_dp])
  int_field = field('integers', states, 1)
  call int_field % set_integer_vector([1])
  xi_field = field('xi', designs, 1)
  call xi_field % set_real_vector([1.0_dp])

  allocate(design % component(1))
  design % component(1) % value = xi_field

  call eta % set_real([1.0_dp])

  !-------------------------------------------------------------------!
  ! The lawful solved two-vertex scalar chain, corrupted per case.
  !-------------------------------------------------------------------!

  allocate(time_graph % vertex(2))
  do v = 1, 2
     allocate(time_graph % vertex(v) % sample % state % component(1))
     time_graph % vertex(v) % sample % state % component(1) % value = one_field
     time_graph % vertex(v) % has_solution = .true.
  end do

  allocate(time_graph % relation(1))
  call builder % bdf_uniform(1, 1.0_dp, time_graph % relation(1) % motif)
  time_graph % relation(1) % sample_vertex   = [1, 2]
  time_graph % relation(1) % unknown_sample  = 2
  time_graph % relation(1) % evaluation_time = 1.0_dp

  select case (trim(which))

  case ('deg0')

     options % max_degree = 0
     call driver % solve_all(r_form, time_graph, design, eta, options, result)

  case ('negsing')

     options % singular_tolerance = 0.0_dp
     call driver % solve_relation(r_form, time_graph, 1, design, eta, &
          & vd, options, step)

  case ('norelations')

     deallocate(time_graph % relation)
     call driver % solve_all(r_form, time_graph, design, eta, options, result)

  case ('idx0')

     call driver % solve_relation(r_form, time_graph, 0, design, eta, &
          & vd, options, step)

  case ('idxhigh')

     call driver % solve_relation(r_form, time_graph, 2, design, eta, &
          & vd, options, step)

  case ('degdim')

     call driver % solve_relation(r_form, time_graph, 1, design, eta, &
          & vd_deep, options, step)

  case ('vertdim')

     call driver % solve_relation(r_form, time_graph, 1, design, eta, &
          & vd_narrow, options, step)

  case ('noeta')

     call eta % clear()
     call driver % solve_relation(r_form, time_graph, 1, design, eta, &
          & vd, options, step)

  case ('unsolved')

     time_graph % vertex(2) % has_solution = .false.
     call driver % solve_relation(r_form, time_graph, 1, design, eta, &
          & vd, options, step)

  case ('unsolvedhistory')

     time_graph % vertex(1) % has_solution = .false.
     call driver % solve_relation(r_form, time_graph, 1, design, eta, &
          & vd, options, step)

  case ('noq')

     time_graph % vertex(2) % sample % state % component(1) % value = int_field
     call driver % solve_relation(r_form, time_graph, 1, design, eta, &
          & vd, options, step)

  case ('badseed')

     call vd(1, 1) % set_real([1.0_dp, 2.0_dp])
     call driver % solve_relation(r_form, time_graph, 1, design, eta, &
          & vd, options, step)

  case ('shortres')

     call driver % solve_relation(short, time_graph, 1, design, eta, &
          & vd, options, step)

  case ('singular')

     call driver % solve_relation(design_only, time_graph, 1, design, eta, &
          & vd, options, step)

  case default

     error stop 'refusal: unknown case'

  end select

  write(*,*) 'refusal case survived: ', trim(which)

end program refusal
