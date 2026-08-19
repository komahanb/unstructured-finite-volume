!=====================================================================!
! The refusals that must die at the growth seat, one per
! invocation:
!
!      novertex        growing an empty graph
!      wrongvertex     discarding with a stale vertex index
!      wrongrelation   discarding with a stale relation index
!      nosentinel      a candidate relation that never names its
!                      own vertex
!      wrongunknown    a candidate whose unknown is not the
!                      appended vertex
!      badmotif        a candidate row weighted thrice over two
!                      vertices - the graph's law, propagated
!      unsolvedhistory a candidate reading unsolved history - the
!                      forward driver's law, propagated
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
  use gti_design_bundles   , only : gti_design_bundle
  use gti_time_motif_builders, only : gti_time_motif_builder
  use gti_time_graphs      , only : gti_time_graph
  use gti_time_forward_drivers, only : gti_time_forward_options
  use gti_time_adaptive_growth_drivers, only : gti_time_growth_candidate, &
       & gti_time_growth_step_result, gti_time_adaptive_growth_driver
  use gti_toy_forms        , only : toy_time_residual_form

  implicit none

  type(graph) :: states, designs
  type(field) :: q1_field, zero_field, xi_field

  type(gti_time_motif_builder) :: builder
  type(gti_time_graph)         :: time_graph, empty_graph
  type(gti_design_bundle)      :: design

  type(gti_time_adaptive_growth_driver) :: grower
  type(gti_time_growth_candidate)       :: candidate
  type(gti_time_growth_step_result)     :: step
  type(gti_time_forward_options)        :: options

  type(toy_time_residual_form) :: r_form
  integer :: av, ar
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
  ! The lawful base: one solved vertex, and one lawful BDF1
  ! candidate, corrupted per case.
  !-------------------------------------------------------------------!

  allocate(time_graph % vertex(1))
  allocate(time_graph % vertex(1) % sample % state % component(1))
  time_graph % vertex(1) % sample % state % component(1) % value = q1_field
  time_graph % vertex(1) % has_solution = .true.

  candidate % vertex % sample % time = 0.5_dp
  allocate(candidate % vertex % sample % state % component(1))
  candidate % vertex % sample % state % component(1) % value = zero_field
  call builder % bdf_uniform(1, 0.5_dp, candidate % relation % motif)
  candidate % relation % sample_vertex   = [1, -1]
  candidate % relation % unknown_sample  = 2
  candidate % relation % evaluation_time = 0.5_dp

  r_form % nequations = 3

  select case (trim(which))

  case ('novertex')

     call grower % append_candidate(empty_graph, candidate, av, ar)

  case ('wrongvertex')

     call grower % append_candidate(time_graph, candidate, av, ar)
     call grower % discard_last_candidate(time_graph, 1, ar)

  case ('wrongrelation')

     call grower % append_candidate(time_graph, candidate, av, ar)
     call grower % discard_last_candidate(time_graph, av, 0)

  case ('nosentinel')

     candidate % relation % sample_vertex = [1, 1]
     call grower % append_candidate(time_graph, candidate, av, ar)

  case ('wrongunknown')

     candidate % relation % unknown_sample = 1
     call grower % append_candidate(time_graph, candidate, av, ar)

  case ('badmotif')

     candidate % relation % motif % rule(1) % weight = &
          & [1.0_dp, 0.0_dp, 0.0_dp]
     call grower % append_candidate(time_graph, candidate, av, ar)

  case ('unsolvedhistory')

     time_graph % vertex(1) % has_solution = .false.
     call grower % try_candidate(r_form, time_graph, candidate, design, &
          & options, .true., step)

  case default

     error stop 'refusal: unknown case'

  end select

  write(*,*) 'refusal case survived: ', trim(which)

end program refusal
